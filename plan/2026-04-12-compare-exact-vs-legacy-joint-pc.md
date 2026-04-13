# Compare Exact Matern to Legacy INLA with Joint PC Priors

## Goal

Run a lightweight internal sanity check that compares:

- the current exact Gaussian EB Matern implementation, and
- an INLA-based Matern fit with joint PC priors on both range and sigma.

## Planned Changes

- Extend the internal comparison script so the legacy INLA route fits both hyperparameters instead of only a range-like quantity.
- Use a slightly smoother truth by increasing the true range.
- Use noisier observations so the smoothing behavior is easier to compare visually.
- Save a short markdown/HTML note, a comparison figure, and small CSV summaries for hyperparameters, runtime, posterior similarity, and objective values.

## Deliverables

- Updated `internal/simulations/compare_matern_exact_vs_legacy_inla.R`
- Refreshed comparison outputs in `internal/simulations/results/`
