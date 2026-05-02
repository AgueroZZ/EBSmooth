# Plan: Laplace Fisher Backend and Partial Parameter Fixing

Date: 2026-05-02

## Goals

- Add public `backend = "laplace_fisher"` for log-link Matern and L-GP smoothers.
- Make `backend = "auto"` resolve to Fisher Laplace for log-link Matern and L-GP fits.
- Preserve explicit observed-Hessian behavior for `backend = "laplace"`, `laplace_r`, `laplace_tmb`, and INLA comparisons.
- Add `fix_params` for partial EB parameter fixing while preserving `fix_g = TRUE`.

## Implementation Plan

1. Add shared `fix_params` validation helpers and fixed-beta resolution.
2. Factor log-link observation curvature so the conditional mode objective remains unchanged, while Fisher Laplace uses `exp(eta)^2 / s^2` for covariance/logdet curvature.
3. Route log-link Matern `auto` fits to `laplace_fisher`, including known noise, learned noise, and two-dimensional fixed-beta cases.
4. Implement Matern range/sigma partial fixing in exact and Laplace optimizers; reject range/sigma fixing for explicit INLA.
5. Add an R L-GP Laplace path for log-link Fisher fits and for fixed L-GP scale.
6. Propagate backend, curvature, and `log_likelihood_semantics` through `eb_smoother()` and `spatial_*` score tables.
7. Update permutation `refit = FALSE` to merge user `fix_params` with fixed observed Matern range/sigma.
8. Update tests, roxygen documentation, README notes, simulation helper backend checks, and this plan/log record.

## Compatibility Notes

- `backend = "laplace"` remains observed-Hessian Laplace and remains the route for reproducing old behavior.
- `backend = "laplace_fisher"` is valid only with `link = "log"`.
- `fix_g = TRUE` maps to `fix_params = c("range", "sigma")` for Matern and `fix_params = "scale"` for L-GP.
- Fixed learned-noise parameters are out of scope; fixed observation noise remains represented by supplying `s`.

## Dispatch Correction

- Matern `backend = "laplace_fisher"` should use the TMB implementation whenever
  the corresponding observed-Hessian `backend = "laplace"` TMB path is
  supported.
- The Fisher option should not introduce a separate R-only optimization path.
  It reuses the TMB Step A mode/hyperparameter fit and changes the returned
  Step B posterior approximation curvature to Fisher/Gauss-Newton curvature.
- `laplace_r` remains an explicit reference/testing backend, not the default
  public path for Fisher fits.
