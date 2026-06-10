# Non-Identity Initialization Scale Fix

## Goal

Fix log-link and softplus-link initialization so Matern and L-GP starts use
compatible response, latent, and raw-noise scales without adding public API.

## Implementation Plan

- Add shared internal helpers for stable inverse-softplus, positive response
  floors, response-mean beta initialization, response means from eta, and raw
  residual noise scales.
- Update Matern initialization so `beta_init` matches the response-scale mean,
  `sigma_data` is computed on the latent eta scale with delta-method known-noise
  weights, and learned `noise_sd_init` is computed from raw response residuals.
- Update L-GP initialization so the known-noise TMB path, R Laplace/Fisher path,
  and learned-noise path share one resolver. When beta is not supplied, use the
  response-scale mean for the intercept and zero for remaining global trend
  coefficients.
- Preserve explicit `beta_fixed`, explicit `g_init$beta`, fixed-scale settings,
  PC-prior overrides, backend routing, likelihood semantics, and
  `pql_inner_iter` behavior.

## Validation

- Add Matern tests for log/softplus response-mean beta starts, delta-method
  latent scale weights, and raw learned-noise starts.
- Add L-GP tests for known-noise default beta starts, learned-noise beta starts,
  explicit beta overrides, and scale-only `g_init` behavior.
- Run targeted Matern and L-GP test files, then the full testthat directory if
  runtime permits.
