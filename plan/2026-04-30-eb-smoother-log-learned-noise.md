# Add Log-Link Learned-Noise Support to `eb_smoother()`

## Scope

- Add `eb_smoother(..., family = "matern", link = "log", s = NULL)`.
- Add `eb_smoother(..., family = "nonspatial", link = "log")`.
- Keep `beta_prec = 0` as a flat prior on the linear-predictor intercept `beta`
  for every link.
- Keep the nonspatial family as the constant-baseline model.

## Design

- Matern log-link learned-noise fits route through the existing Laplace
  machinery with `log_noise_sd` as an optimized hyperparameter.
- The TMB Matern objective supports `learn_noise` and includes the common
  learned noise SD in the Gaussian observation likelihood.
- Matern `backend = "laplace"` uses TMB when supported and falls back to the R
  Laplace reference; explicit `laplace_r`, `laplace_tmb`, `inla`, and
  `inla_pc` remain available.
- INLA log-link learned-noise fits use external local profiling over the common
  noise SD and evaluate the package Laplace objective at the selected INLA mode
  for the primary comparable `log_likelihood`.
- Nonspatial log-link fits model `x_i ~ N(exp(beta), s_i^2)` for known noise and
  `x_i ~ N(exp(beta), sigma^2)` for learned noise.

## Accuracy Notes

- Flat and proper log-link beta modes use the same flat-on-beta and
  Gaussian-on-beta semantics as identity-link fits.
- For log-link flat beta, the implementation follows the requested
  flat-on-beta local Laplace semantics. No Jacobian term for `exp(beta)` is
  added.
- TMB learned-noise empirical-Bayes fits use a small deterministic multi-start
  set to avoid a local mode where the smoother collapses to a constant signal
  with large observation noise.

## Validation

- Compare nonspatial log known-noise fixed beta against the direct Gaussian
  likelihood at `exp(beta)`.
- Compare nonspatial log empirical-Bayes beta against the positive weighted
  mean MLE on the beta scale.
- Compare Matern log learned-noise `laplace_tmb` and `laplace_r` for
  empirical-Bayes, fixed, flat-prior, and proper-prior beta modes.
- Check that Matern log learned-noise PC-noise priors are included in both R
  Laplace and TMB objectives.
- Smoke-check explicit INLA learned-noise log-link fixed-beta fits against the
  package Laplace comparable objective.
