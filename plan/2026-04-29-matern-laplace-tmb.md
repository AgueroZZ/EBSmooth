# TMB Acceleration for Matern Laplace Backend

## Summary
- Add a TMB/C++ implementation for the package-owned Matern Laplace backend.
- Keep the existing R Laplace implementation as a correctness reference and fallback.
- Route `backend = "laplace"` to TMB when supported and to the R reference path otherwise.

## Implementation Plan
- Extend the package TMB objective with an internal `model_id` selector:
  - `model_id = 0` keeps the existing L-GP objective.
  - `model_id = 1` adds the Matern objective for known-noise `alpha = 2` SPDE models.
- Build the Matern precision in C++ from INLA's `M0`, `M1`, and `M2` basis:
  - `Q = tau^2 * (kappa^4 M0 + 2 * kappa^2 M1 + M2)`.
  - Convert from `log_range` and `log_sigma` using the existing SPDE parameterization.
- Use `density::GMRF(Q)(w)` for the latent Gaussian prior and TMB's Laplace integration.
- Preserve beta semantics:
  - Fixed beta is mapped/fixed in TMB.
  - Empirical-Bayes beta is optimized as a fixed-effect parameter.
  - Proper-prior beta is integrated as a random effect with its Gaussian prior.
  - Flat-prior beta is integrated as a random effect without a beta prior term.
- Add explicit `backend = "laplace_tmb"` and `backend = "laplace_r"` options.
- Add diagnostics through `laplace_implementation`, with values `"tmb"` or `"r"`.

## Validation Plan
- Compare `laplace_tmb` with `laplace_r` for log-link EB, fixed-beta, and proper-prior beta fits.
- Check that default `backend = "laplace"` uses TMB for fixed, EB, flat-prior, and proper-prior beta modes when otherwise supported.
- Check identity-link `laplace_tmb` against the exact Gaussian backend.
- Confirm explicit `backend = "laplace_tmb"` supports flat-prior beta and matches the R reference objective.
- Run focused Matern, eb_smoother, and L-GP tests to catch routing or TMB objective regressions.
