# EBSmoothr 0.2.5 Performance Review

## Scope

Reviewed the uncommitted EBSmoothr 0.2.5 performance changes relative to commit
`1960a80`, focusing on the sparse SPD factorization cache, direct stationary
Matern SPDE precision assembly, exact Matern data-statistics hoisting, and
Fisher-PQL early stopping diagnostics.

## Review outcome

- Approved the stationary alpha-2 SPDE direct precision formula. For a
  stationary template, INLA's diagonal expression
  `D0 * (D1^2 * M0 + D1 * D2 * (M1 + t(M1)) + M2) * D0` reduces to the scalar
  combination used by `.matern_spde_precision_direct()`.
- Approved learned-noise data-statistics scaling. The cached unit-noise
  statistics scale by `1 / noise_sd^2`, and `logdet_D` gains
  `2 * n_obs * log(noise_sd)`.
- Approved CHOLMOD numeric refactorization reuse after exact sparsity-pattern
  matching on `i` and `p`. A local Matrix probe confirmed
  `Matrix::.updateCHMfactor()` returns a separate updated factor without
  mutating the original factor object.
- Adjusted the SPD factorization cache eviction to remove the oldest cached
  entry when the configured maximum is reached, matching the documented
  behavior instead of clearing the entire cache.
- Approved Fisher-PQL early stopping and diagnostic truncation. Callers consume
  `eta_change` and `stepA_log_marginals` by their actual lengths or by
  `inner_iterations`, and `max_inner_iter` preserves the configured cap.

## Validation

- `Rscript -e 'devtools::test()'`
  - Result: 755 passing tests, 0 failures, 0 warnings, 0 skips.
- `git diff --check`
  - Result: clean.
- `Rscript -e 'devtools::check(args = c("--no-manual", "--no-build-vignettes"), error_on = "never")'`
  - Result: 0 errors, 0 warnings, 4 notes.
  - Notes were for current-time verification, existing non-standard top-level
    files, existing `unlockBinding()` calls, and native routine registration.

## Integration notes

- The benchmark scripts in `EBSmoothr/bench/` are suitable to keep as developer
  harnesses, while generated `results_*.rds` files should remain untracked.
- Unrelated EBFPCA workspace changes were intentionally left unstaged.
