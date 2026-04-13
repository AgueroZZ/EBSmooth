# Matern `clinear` vs Exact Beta-EB Benchmark

This benchmark compares two implementations of the `betaprec < 0`
interpretation for the Matern smoother:

- `exact_profile`: exact Gaussian profile objective plus the PC prior;
- `clinear_inla`: INLA Step A with `clinear` used to represent the
  intercept as a hyperparameter.

The timings below are fit-phase timings after the spatial setup has already
been built. This mirrors the package pattern in which the generator is built
once for a fixed set of locations and then reused for fitting.

## Summary

case | n | reps | exact_mean_seconds | exact_sd_seconds | clinear_mean_seconds | clinear_sd_seconds | speedup_exact_over_clinear | exact_eval_count | objective_gap | beta_gap | range_gap | sigma_gap
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
1d_n120 | 120.000000 | 3.000000 | 0.463000 | 0.019157 | 2.063000 | 0.026851 | 0.224430 | 106.000000 | -1.6e-05 | -0.002408 | -0.001908 | 0.000081
2d_n252 | 252.000000 | 2.000000 | 4.801500 | 0.634275 | 2.221000 | 0.005657 | 2.161864 | 116.000000 | -5.0e-06 | -0.000834 | 0.000223 | 0.000126
2d_n432 | 432.000000 | 1.000000 | 10.927000 | NA | 2.351000 | NA | 4.647809 | 121.000000 | -1.0e-05 | -0.001806 | 0.000215 | 0.000163

## Interpretation

- `speedup_exact_over_clinear > 1` means the exact/profile method is slower
  than the `clinear`-based INLA fit by that factor.
- The parameter and objective gaps are included to confirm that both methods
  still converge to essentially the same fitted solution.
