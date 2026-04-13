# Compare Exact Matern to Legacy INLA with Joint PC Priors

## What changed

- Replaced the earlier legacy comparison, which only optimized a range-like hyperparameter, with a joint-PC-prior INLA comparison that fits both range and sigma.
- Increased the simulated truth's range so the latent surface is a bit smoother.
- Increased the observation noise so the denoising comparison is more visible.

## Why

Comparing the exact implementation to a legacy INLA route is more meaningful when both methods are allowed to move both main smoothing hyperparameters. This keeps the sanity check lightweight while making the comparison more interpretable.

## Key findings

- On a single `20 x 16` 2D example with truth `(range, sigma, beta0) = (0.30, 0.35, 0.20)` and `noise sd = 0.22`, the exact fit gave `(0.2914, 0.3527, 0.1959)`.
- The legacy INLA fit with joint PC priors gave `(0.3035, 0.3630, 0.1958)`.
- The two posterior means were nearly identical: correlation `0.999996` and RMSE `0.001070`.
- The exact objective gap between the exact optimum and the legacy fitted `(range, sigma)` pair was only `0.027654`.
- The legacy raw `mlik` matched the exact objective evaluated at the legacy fitted hyperparameters up to numerical rounding.
- Runtime remained in the same general order of magnitude, but the exact implementation was slower on this example (`14.389s` versus `4.857s`).
