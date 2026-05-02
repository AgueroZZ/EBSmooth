# Matern Laplace TMB Acceleration Log

## Changes
- Added a `model_id` selector to the package TMB objective.
- Preserved the existing L-GP TMB path under `model_id = 0`.
- Added a Matern TMB path under `model_id = 1` for known-noise `alpha = 2` fits.
- Added `backend = "laplace_tmb"` and `backend = "laplace_r"` for explicit implementation selection.
- Updated `backend = "laplace"` to use TMB when supported and fall back to the R reference implementation otherwise.
- Added `laplace_implementation` diagnostics to Matern Laplace fits.
- Updated flat-prior beta semantics so the TMB path integrates beta as a random effect without adding a beta prior term, matching the L-GP convention.

## Validation Notes
- Quick local sanity checks showed `laplace_tmb` and `laplace_r` matching for log-link empirical-Bayes, fixed-beta, and proper-prior beta fits.
- Identity-link `laplace_tmb` matched the exact Gaussian backend to numerical tolerance in a small smoke test.
- Explicit `backend = "laplace_tmb"` supports flat-prior beta, and default `backend = "laplace"` uses TMB for flat-prior beta when the rest of the Matern TMB path is supported.
- Focused tests passed for `test-matern.R`, `test-eb-smoother.R`, and `test-lgp.R`.
- A local `n = 60` log-link empirical-Bayes benchmark took about 1.17 seconds with `laplace_tmb` versus 13.98 seconds with `laplace_r`, with log-likelihood difference below `1e-6`.
- Additional TMB-only smoke benchmarks took about 1.51 seconds at `n = 120` and 2.39 seconds at `n = 300`.
