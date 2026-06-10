# Plan: Speed up 2D Matérn smoothing for non-identity links

## Context

`eb_smoother(..., family = "matern", link = "softplus" | "log")` is **50–200× slower** than `link = "identity"` at medium scale (n ≈ 500–2000, n_spde ≈ 500–2000). Theoretical Laplace-vs-closed-form gap should be ~10×, so there is at least one O(n_spde²)/O(n_spde³) operation hidden in the path. The user has not profiled yet, and TMB+SPDE are already in place — the slowness is happening despite the existing C++ infrastructure.

After tracing the code, the dominant suspect is **a silent TMB→R fallback** that happens whenever any element of the TMB-computed posterior moments comes out non-finite. The R reference path is dramatically slower (R-side Newton with R callbacks for objective/gradient, INLA precision rebuilt every outer step, triple-nested loops for empirical-Bayes β). A second large suspect is `.compute_diag_A_Qinv_At` densifying the full inverse of an n_spde×n_spde sparse SPD matrix.

The plan is **diagnose first, then fix the highest-ROI items**, then escalate only if needed.

---

## Phase 1 — Diagnostics (do this BEFORE any optimization)

### 1.1 Surface the silent fallback

[R/02_Matern.R:3621-3628](R/02_Matern.R) and [R/02_Matern.R:3721-3728](R/02_Matern.R) currently fall back from TMB to the R reference path with no warning. Make this visible:

```r
# inside .fit_matern_laplace_dispatch_known_noise and ..._unknown_noise
if (!is.null(tmb_invalid_reason) && !identical(backend_use, "laplace_tmb")) {
  if (getOption("EBSmoothr.warn_tmb_fallback", TRUE)) {
    warning("Matern TMB Laplace fit fell back to R reference path: ",
            tmb_invalid_reason, call. = FALSE)
  }
}
return(fit_r_reference(tmb_fallback_error = tmb_invalid_reason))
```

This single change tells us, for any real run, whether (A) is happening. The fallback condition `.matern_laplace_tmb_invalid_reason` ([R/02_Matern.R:2683-2706](R/02_Matern.R)) trips on **any** non-finite entry in `posterior$mean / var / second_moment`. For `link = "log"` the moment formula `exp(2*eta_mean + eta_var) * (exp(eta_var) - 1)` overflows easily when `eta_var` is moderate; one bad element kills the whole TMB result.

### 1.2 Add a profiling script

Create `bench/profile_matern.R` (new file):

```r
library(EBSmoothr); set.seed(1)
loc <- as.matrix(expand.grid(seq(0,1,length.out=32), seq(0,1,length.out=32)))
n <- nrow(loc); truth <- 2 + sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
x <- truth + rnorm(n, sd=0.3); s <- rep(0.3, n)

Rprof(tmp <- tempfile(), interval = 0.01, line.profiling = TRUE)
fit <- eb_smoother(x, s, family="matern", locations=loc, link="softplus")
Rprof(NULL)
print(summaryRprof(tmp, lines="show")$by.self[1:30,])

# Diagnostic flags
print(fit$laplace_diagnostics$tmb_fallback_error)
print(fit$laplace_diagnostics$tmb_mode_source)
```

Run once on a representative problem. The `by.self` table tells us whether time is in `inla.spde2.precision`, `Matrix::solve`, `nlminb` callbacks, or TMB's `MakeADFun`. Don't optimize without this.

---

## Phase 2 — Quick wins (highest ROI, do in order)

### 2.1 Fix the fallback root cause (target: stop falling back)

If 1.1 confirms the fallback is firing, the fix is to make `.matern_response_moments_from_eta` ([R/02_Matern.R:2230-2248](R/02_Matern.R)) and `.matern_laplace_fit_has_finite_posterior` ([R/02_Matern.R:2666-2681](R/02_Matern.R)) robust to large `eta_var` for the **log** link:

- Compute the log-normal moments in log-space and clip / NA-fill the bad indices instead of returning `Inf`.
- Treat any element with `eta_var > THRESHOLD` (say 50) as effectively unbounded and replace with NA in posterior — but do NOT mark the whole fit as invalid for that. Adjust `.matern_laplace_fit_has_finite_posterior` to allow a small number of NA entries (return TRUE if ≥99% of entries are finite, accompanied by a warning).
- For softplus the GH quadrature ([R/00_setup.R:59-85](R/00_setup.R)) is already stable, but double-check `.softplus_stable` doesn't overflow at the GH nodes.

Expected impact: **gets the TMB path actually used**. This alone could close most of the gap.

### 2.2 Add `src/Makevars` and `src/Makevars.win` (free 1.3-2× across the board)

Create [src/Makevars](src/Makevars) and [src/Makevars.win](src/Makevars.win):

```
CXX_STD = CXX17
PKG_CXXFLAGS = -O3 -DNDEBUG -fno-math-errno
```

(Skip `-march=native` for CRAN compatibility; user can override locally via `~/.R/Makevars`.)

Then rebuild: `rm src/EBSmoothr.so src/EBSmoothr.o; devtools::load_all()`.

### 2.3 Replace `.compute_diag_A_Qinv_At` with Takahashi-aware version

[R/01_LGP.R:32-37](R/01_LGP.R) currently solves `Q V = A^T` densely. When `A = Diagonal(n_spde)`, `V` is the **full dense inverse** — O(n_spde²) memory and O(n_spde · nnz(L)) time. The proper fix uses `INLA::inla.qinv` (INLA is already in `Imports`):

```r
.compute_diag_A_Qinv_At <- function(A, Q) {
  factor <- .factorize_spd(Q)
  Q_mat  <- factor$matrix
  use_qinv <- inherits(A, "ddiMatrix") || is(A, "diagonalMatrix") ||
              .pattern_subset_of_Q(A, Q_mat)
  if (isTRUE(use_qinv)) {
    qinv <- tryCatch(INLA::inla.qinv(Q_mat), error = function(e) NULL)
    if (!is.null(qinv)) {
      if (inherits(A, "ddiMatrix") || is(A, "diagonalMatrix")) {
        return(as.numeric(Matrix::diag(qinv)))
      }
      AT <- Matrix::t(A)
      V  <- qinv %*% AT
      return(as.numeric(Matrix::rowSums(A * Matrix::t(V))))
    }
  }
  # fallback: original dense path with a warning
  A <- Matrix::Matrix(A, sparse = TRUE)
  V <- .solve_spd_factor(factor, Matrix::t(A))
  as.numeric(Matrix::rowSums(A * Matrix::t(V)))
}

.pattern_subset_of_Q <- function(A, Q) TRUE  # SPDE projection: 3 nz per row,
# row-pair pattern is contained in mesh adjacency = pattern(Q). Exact check is
# possible but unnecessary — the fallback handles bad cases.
```

Affected callers: [R/02_Matern.R:601](R/02_Matern.R), [R/02_Matern.R:677](R/02_Matern.R), [R/02_Matern.R:682](R/02_Matern.R), [R/02_Matern.R:2391](R/02_Matern.R), [R/02_Matern.R:2399](R/02_Matern.R), [R/02_Matern.R:2967](R/02_Matern.R), [R/02_Matern.R:2975](R/02_Matern.R), [R/01_LGP.R:995](R/01_LGP.R), [R/01_LGP.R:1030](R/01_LGP.R), [R/03_eb_smoother.R:877](R/03_eb_smoother.R), [R/03_eb_smoother.R:920](R/03_eb_smoother.R) — all keep the same signature, no caller changes needed.

Correctness: for any `A` whose row-pair pattern lies in pattern(Q) (the SPDE projection always satisfies this — neighbors in the same triangle are adjacent in Q), `inla.qinv(Q)` returns all entries needed. For `Diagonal`, only `diag(Qinv)` is needed and Takahashi gives that. Tests in [tests/testthat/test-matern.R](tests/testthat/test-matern.R) should still pass. Add a regression test asserting the new and old paths agree numerically on a small mesh.

Expected impact: **5–20× on 2D meshes with n_spde ≥ 1000**. Helps identity, log, AND softplus.

### 2.4 Cache `MakeADFun` across the 5 log-link starts

[R/02_Matern.R:3171-3213](R/02_Matern.R) hard-codes 5 starts for `link = "log"` + EB beta + `!fix_g`, rebuilding `MakeADFun` each time. The TMB tape and symbolic Cholesky pattern don't depend on the starting values — build once, reset `obj$par[]` between starts:

```r
objA_shared <- make_objA(par0)   # build once
stepA_results <- lapply(seq_along(stepA_starts), function(i) {
  par_i <- stepA_starts[[i]]
  # write par_i into objA_shared$par by name
  for (nm in names(par_i)) {
    idx <- which(names(objA_shared$par) == nm)
    if (length(idx)) objA_shared$par[idx] <- as.numeric(par_i[[nm]])[seq_along(idx)]
  }
  bounds <- make_bounds(objA_shared)
  opt <- tryCatch(.matern_laplace_tmb_optimize(
      objA_shared, "TMB Matern Laplace Step A",
      lower = bounds$lower, upper = bounds$upper),
    error = function(e) e)
  if (inherits(opt, "error")) return(NULL)
  list(opt = opt, start_index = i,
       last_par_best = objA_shared$env$last.par.best)
})
# After picking winner, restore its last.par.best for mode extraction
```

Expected impact: ~0.3–1 s saved per start × 4 redundant builds = **1–4 s for log link EB**.

### 2.5 Tighten Step B refinement trigger

[R/02_Matern.R:3266](R/02_Matern.R) and [R/02_Matern.R:3278](R/02_Matern.R) treat any `optA$convergence != 0` as "refine", and run `nlminb(eval.max = 20000, iter.max = 20000)`. nlminb code 1 ("relative convergence") is fine. Tighten:

```r
refine_mode <- !had_stepA_mode || (
  !(as.integer(optA$convergence) %in% c(0L, 1L)) &&
  max(abs(objA$gr(optA$par))) > 1e-3
)
# and lower the budgets
optB <- stats::nlminb(start = objB$par, objective = objB$fn, gradient = objB$gr,
                      control = list(eval.max = 2000, iter.max = 2000))
```

Expected impact: avoids spurious 20000-iter runs.

---

## Phase 3 — Larger refactors (only if Phase 2 isn't enough)

After Phase 2, re-run the benchmark. If softplus/log are still > 10× identity:

- **3.1** Cache the symbolic Cholesky of the SPDE precision in the R reference path. [R/02_Matern.R:192-211](R/02_Matern.R) (`.matern_precision_from_log_params`) re-factorizes from scratch every outer step. Use `Matrix::update(cholQ, Q)` to keep the AMD permutation and skip symbolic factorization. Touch `.factorize_spd` in [R/01_LGP.R:2-25](R/01_LGP.R).
- **3.2** Add `compute_latent_variance` flag to the post-fit summary so users who only need `posterior` (not `posterior_spatial_field`) skip the second `.compute_diag_A_Qinv_At` call ([R/02_Matern.R:2975](R/02_Matern.R)).
- **3.3** Report `eta_mean`, `eta_var` directly via TMB `REPORT()` so the post-fit step doesn't recompute the H factorization in R.

---

## Verification

Create `bench/bench_links.R` (run before and after Phase 2):

```r
library(EBSmoothr); library(microbenchmark); set.seed(1)
sizes <- list(small = 14, medium = 32, large = 56)  # n_spde ~ 200, 1000, 3100
for (nm in names(sizes)) {
  k <- sizes[[nm]]
  loc <- as.matrix(expand.grid(seq(0,1,length.out=k), seq(0,1,length.out=k)))
  n <- nrow(loc); truth <- 2 + sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
  x <- truth + rnorm(n, sd=0.3); s <- rep(0.3, n)
  cat("\n=== ", nm, " (n = ", n, ") ===\n", sep="")
  bm <- microbenchmark(
    identity = eb_smoother(x, s, family="matern", locations=loc, link="identity"),
    softplus = eb_smoother(x, s, family="matern", locations=loc, link="softplus"),
    log      = eb_smoother(x, s, family="matern", locations=loc, link="log"),
    times = if (nm == "large") 1L else 3L
  )
  print(bm, unit = "s")
  # confirm no fallback
  fit_sp <- eb_smoother(x, s, family="matern", locations=loc, link="softplus")
  stopifnot(is.null(fit_sp$laplace_diagnostics$tmb_fallback_error))
}
```

Then run existing tests:
```bash
cd /Users/ziangzhang/Desktop/EBSmooth/EBSmoothr
Rscript -e 'devtools::test()'
```

Numerical regression: in [tests/testthat/test-matern.R](tests/testthat/test-matern.R) add a test that the new `.compute_diag_A_Qinv_At` agrees with the dense-inverse reference within 1e-10 on a small mesh.

**Targets**:
- After Phase 2: softplus and log within 5–10× of identity at n_spde = 3000
- No silent fallback fires under typical inputs
- All existing tests pass

## Critical files

- [R/02_Matern.R](R/02_Matern.R) — backend dispatch, TMB Laplace path, R reference path, post-fit summary
- [R/01_LGP.R](R/01_LGP.R) — `.compute_diag_A_Qinv_At`, `.factorize_spd`, `.solve_spd_factor`
- [R/00_setup.R](R/00_setup.R) — `.softplus_gaussian_moments` (already vectorized, fine)
- [src/EBSmoothr.cpp](src/EBSmoothr.cpp) — TMB objective (already correct; only need build flags)
- [DESCRIPTION](DESCRIPTION) — confirms INLA in Imports (so `inla.qinv` is available)
- New: `src/Makevars`, `src/Makevars.win`, `bench/profile_matern.R`, `bench/bench_links.R`

## Reusable helpers already in the codebase
- `.factorize_spd` ([R/01_LGP.R:2](R/01_LGP.R)) — sparse Cholesky via Matrix::Cholesky, returns logdet
- `.solve_spd_factor` ([R/01_LGP.R:27](R/01_LGP.R)) — sparse triangular solve
- `.matern_observation_terms` ([R/02_Matern.R:2187](R/02_Matern.R)) — vectorized link gradient/Hessian
- `INLA::inla.qinv` — sparse Takahashi recursion, already importable
- `Matrix::update(cholQ, newQ)` — reuse symbolic factorization (Phase 3.1)
