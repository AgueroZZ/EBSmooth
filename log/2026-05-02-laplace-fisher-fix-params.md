# Log: Laplace Fisher Backend and Partial Parameter Fixing

Date: 2026-05-02

## Summary

Implemented `backend = "laplace_fisher"` for log-link Matern and L-GP smoothers and added `fix_params` for partial EB parameter fixing.

## Changes

- Matern `backend = "auto"` now resolves to `laplace_fisher` for `link = "log"` across known-noise, learned-noise, and two-dimensional fixed-beta cases.
- L-GP `backend = "auto"` now resolves to `laplace_fisher` for `link = "log"` and keeps TMB for identity-link default fits.
- Fisher Laplace uses the original log-posterior conditional mode objective, but replaces the log-link normal observation curvature used in the Laplace covariance/logdet from `exp(eta) * (2 * exp(eta) - x) / s^2` to `exp(eta)^2 / s^2`.
- Fisher fits return `backend = "laplace_fisher"`, `laplace_curvature = "fisher"`, and `log_likelihood_semantics = "laplace_fisher_<beta_mode>"`.
- Explicit `backend = "laplace"`, `laplace_r`, `laplace_tmb`, and `inla` keep observed-Hessian semantics.
- Added `fix_params` to `eb_smoother()`, `ebnm_Matern_generator()` returned functions, and `ebnm_LGP_generator()` returned functions.
- Matern `fix_params` supports `range`, `sigma`, and `beta`; L-GP supports `scale` and `beta`.
- `fix_g = TRUE` remains supported as a shortcut for Matern range/sigma or L-GP scale.
- Spatial score tables now carry Laplace curvature fields and log-link Matern spatial scoring reports Fisher score semantics by default.
- Permutation scoring with `refit = FALSE` now merges user `fix_params` with the fixed observed Matern range/sigma behavior.
- Added tests for Fisher backends, beta modes, partial fixing, Fisher curvature stability, spatial score semantics, and permutation p-values.

## Validation

- Targeted tests were updated for `test-matern.R`, `test-lgp.R`, `test-eb-smoother.R`, and `test-spatial-scores.R`.
- Roxygen documentation and generated Rd files should be regenerated after this change.

## Follow-up: Matern Fisher TMB Dispatch

- Corrected Matern `backend = "laplace_fisher"` dispatch so it uses the TMB
  implementation whenever the observed-Hessian `backend = "laplace"` TMB path is
  supported.
- The TMB Step A optimizer is unchanged. It still optimizes the original
  observed-Hessian Laplace objective and extracts the conditional mode from TMB.
- The returned `laplace_fisher` fit now rebuilds the Step B posterior precision,
  covariance summaries, and Fisher Laplace score at the TMB mode using
  Fisher/Gauss-Newton observation curvature.
- `laplace_r` remains available as an explicit reference/testing backend; the
  public default API remains `backend = "auto"`, `"laplace"`, or
  `"laplace_fisher"`.
- Regenerated `ebnm_Matern_generator.Rd`.
- Re-ran targeted `test-matern.R`, `test-eb-smoother.R`, and
  `test-spatial-scores.R`; all completed successfully. INLA emitted its known
  retry warnings in the Matern test, but testthat finished with success.

## Documentation and Check Pass

- Regenerated roxygen documentation and NAMESPACE for the package.
- Added `flashier` to `Suggests` because spatial scoring optionally extracts
  loadings from flashier objects via `requireNamespace("flashier")`.
- Removed generated `src/*.o` and `src/*.so` files and local R session files
  before source-package checking.
- `R CMD build EBSmoothr --no-build-vignettes` completed successfully.
- `R CMD check --no-manual --no-build-vignettes --no-tests` on the built
  tarball completed with one existing NOTE about `unlockBinding()` in the
  INLA quiet-default helper and no warnings or errors.
- Targeted tests passed for `test-matern.R`, `test-lgp.R`,
  `test-eb-smoother.R`, and `test-spatial-scores.R`. A full testthat check was
  not used as the release gate because the INLA subprocess emitted runaway
  Hessian diagnostics under `R CMD check`; the targeted files completed
  normally.
