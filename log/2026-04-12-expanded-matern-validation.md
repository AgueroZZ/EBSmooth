## Update

Completed the expanded exact-Matern validation pass following the first
simulation summary.

## Completed Additions

- Added a larger 2D recovery experiment with grid sizes `18x14` and `24x18`.
- Increased the saved 2D surface-recovery example to `24x18` observations.
- Added a profile marginal-likelihood surface over `(rho, sigma)` for the large
  2D example, profiling out `beta0` at each grid point.
- Saved the grid values for the profile surface to
  `internal/simulations/results/matern_profile_loglik_surface_2d.csv`.
- Updated the internal markdown and HTML validation summaries.

## Main Findings

- In the larger 2D recovery experiment, the fitted range and sigma are close to
  the truth on average, and the posterior surface correlation is about `0.98`.
- In the large 2D example (`n = 432`), the posterior surface tracks the true
  latent surface well.
- The profile marginal-likelihood surface is plausible: the optimized point
  lies near the true point, and the plotted high-likelihood region contains
  both.

## Notes

- The profile surface corresponds to the actual empirical-Bayes target because
  `beta0` is profiled out at each `(rho, sigma)` pair.
