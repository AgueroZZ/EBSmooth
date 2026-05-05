# Softplus Posterior Moment Audit and Fix

## Goal

Replace first-order delta-method softplus posterior moments with deterministic
marginal Gaussian-transform moments under the Laplace posterior approximation.
For each fitted linear predictor marginal
`eta_i ~ N(eta_mean_i, eta_var_i)`, report:

- `E[softplus(eta_i)]`
- `Var[softplus(eta_i)]`

using fixed Gauss-Hermite quadrature.

## Implementation Plan

- Add shared internal helpers for stable softplus, stable sigmoid, cached
  Gauss-Hermite quadrature, and softplus Gaussian moments.
- Route LGP, Matern, and constant-reference softplus response-moment code
  through the shared deterministic moment helper.
- Keep identity moments and log-link lognormal formulas unchanged.
- Verify TMB uses an AD-safe softplus implementation and that public LGP
  softplus entry points are accepted.
- Update roxygen documentation and generated Rd files to describe softplus
  response-scale moment semantics.

## Validation Plan

- Add direct helper tests against high-accuracy numerical integration over
  `eta_var` values `0`, `0.01`, `0.1`, `0.5`, and `1`.
- Keep posterior sampler consistency tests for LGP and Matern softplus fits,
  using Monte Carlo-aware tolerances.
- Run targeted test files for LGP, Matern, and `eb_smoother`.
- Run source package build and a no-tests package check gate after
  documentation regeneration.
