# Matern `clinear` Beta-EB Spike

This internal experiment checks whether an INLA `pcmatern` fit with
`clinear` used for the intercept can reproduce the exact profiled objective
used by the current exact Gaussian Matern implementation.

## Setup

- `range` and `sigma` use PC priors via `inla.spde2.pcmatern()`.
- The intercept is represented through `f(beta_cov, model = "clinear")`
  with `beta_cov = 1` for every observation.
- The `clinear` hyperparameter uses a flat prior, so this spike targets the
  `betaprec < 0` interpretation in which `beta0` is optimized rather than
  integrated out.

## Results

case | stepA_penalized | exact_profile_plus_prior_at_mode | exact_profile_plus_prior_optimum | gap_mode | gap_optimum | beta_clinear | beta_exact | range_clinear | range_exact | sigma_clinear | sigma_exact | stepA_mlik_integration | stepA_mlik_gaussian
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
1d | -2.845455 | -2.845453 | -2.845429 | -1e-06 | -2.6e-05 | 0.343762 | 0.346508 | 6.477811 | 6.491694 | 0.646478 | 0.647123 | -5.806490 | -2.761992
2d | -28.429858 | -28.429849 | -28.429832 | -1e-05 | -2.6e-05 | -0.000494 | -0.000974 | 0.505571 | 0.504768 | 0.485147 | 0.484773 | -33.492478 | -30.447981

## Takeaway

With the correct `clinear` specification, the INLA Step A penalized
objective is numerically almost identical to the exact profiled objective
plus the PC prior contribution. The fitted `(beta0, range, sigma)` values
are also nearly identical to the exact optimum in both the 1D and 2D
examples.

This supports the idea that `betaprec < 0` may be implementable within the
fast INLA framework by representing the intercept as a `clinear`
hyperparameter instead of a fixed effect.
