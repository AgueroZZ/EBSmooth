# Exact Matern 2D Permutation-Aware Refactor

## What Changed
- Refactored the shared sparse SPD helpers in `EBSmoothr/R/01_LGP.R` so they now:
  - use permutation-aware sparse Cholesky by default
  - accept either a matrix or a precomputed factorization
  - expose a reusable solve path and log-determinant through the same factorization object
- Refactored the exact Matern helpers in `EBSmoothr/R/02_Matern.R` so:
  - `.exact_matern_sufficient_stats()` factorizes `Q` and `Q_post` once and stores both factors in the returned stats object
  - quadratic forms and log-determinants reuse those stored factors
  - `.exact_matern_posterior_from_stats()` reuses the stored `Q_post` factorization for posterior means and posterior variances instead of refactorizing
- Added exact-rigor regression coverage in `EBSmoothr/tests/testthat/test-matern.R` with a test-local legacy `perm = FALSE` reference for small deterministic 2D cases.
- Added `internal/simulations/study_runtime_2d_exact_vs_inla.R` and the corresponding benchmark outputs under `internal/simulations/results/`.

## Accuracy Validation
- Installed the package successfully with:
  - `R CMD INSTALL -l /tmp/ebsmoothr-lib EBSmoothr`
- Ran a targeted 2D exact equivalence check against the test-local legacy `perm = FALSE` reference and required all checks to pass at `1e-8` tolerance with `stopifnot(...)`.
- Observed numerical differences were at floating-point noise level:
  - `max_stats_diff = 1.137e-13`
  - `max_state_mean_diff = 5.551e-16`
  - `max_state_var_diff = 1.735e-17`
  - `max_latent_mean_diff = 5.551e-16`
  - `max_latent_var_diff = 9.437e-16`
  - `unknown_noise_loglik_diff = 1.812e-13`
- The explicit targeted validation command completed with:
  - `2D exact equivalence checks passed`

## Efficiency Validation
- Ran the dedicated 2D runtime study:
  - `Rscript internal/simulations/study_runtime_2d_exact_vs_inla.R`
- Benchmark setup:
  - `n = 3000`
  - `d = 2`
  - regular `60 x 50` grid
  - `max.edge = c(0.1, 0.2)`
  - `n_spde = 3352`
- Recorded runtimes in `internal/simulations/results/runtime_2d_exact_vs_inla.csv`:
  - `matern_exact_known_2d = 7.482 s`
  - `matern_exact_learned_2d = 9.198 s`
  - `matern_pc_learned_2d = 9.147 s`
  - `matern_exact_objective_eval_2d = 0.040 s`
- The single exact objective evaluation is no longer dominated by the old `perm = FALSE` sparse Cholesky bottleneck.

## Notes
- A full installed-package `testthat::test_dir("EBSmoothr/tests/testthat")` run was started, but the suite remained long-running because of existing heavy INLA/Matern tests in the repo. The hard validation recorded for this refactor is therefore:
  - successful package install
  - explicit 2D exact-vs-legacy equivalence pass
  - dedicated 2D runtime benchmark
