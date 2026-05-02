# Make Exact Matern Fast and Rigorous in `d = 2`

## Summary
- Keep `exact` as the default Matern backend when `pc.penalty = NULL`.
- Treat the `d = 2` slowdown as a sparse linear-algebra implementation bug.
- Replace the current `perm = FALSE` sparse Cholesky path with permutation-aware factorization and reuse the resulting factors across the exact Matern objective and posterior reconstruction.
- Validate the refactor against a test-local legacy reference that reproduces the old slow formulas on small deterministic 2D examples.

## Planned Changes
- Refactor `EBSmoothr/R/02_Matern.R` so `.exact_matern_sufficient_stats()` factorizes `Q` and `Q_post` once via shared permutation-aware helpers and reuses those factors for quadratic forms and log-determinants.
- Refactor `.exact_matern_posterior_from_stats()` so it reuses the stored `Q_post` factorization for posterior means and variances instead of refactorizing.
- Extend the shared sparse helpers in `EBSmoothr/R/01_LGP.R` so they accept either a matrix or a precomputed factorization.
- Add 2D rigor regression tests in `EBSmoothr/tests/testthat/test-matern.R` that compare the new exact path against a legacy `perm = FALSE` reference for:
  - exact log-marginal likelihood
  - profiled beta
  - posterior means and variances
  - latent-field posterior summaries
  - learned-noise exact objective
- Add a dedicated 2D runtime benchmark under `internal/simulations/` that reports:
  - exact known-noise runtime
  - exact learned-noise runtime
  - INLA-PC learned-noise runtime
  - single exact objective evaluation runtime

## Validation
- Run `R CMD INSTALL EBSmoothr`.
- Run the full `testthat` suite.
- Run the new 2D runtime study and record the benchmark outputs under `internal/simulations/results/`.

## Assumptions
- Switching from `perm = FALSE` to `perm = TRUE` is mathematically equivalent provided that solves and derived quantities correctly respect the factorization permutation.
- Remaining output differences after the refactor should be at floating-point noise level only.
