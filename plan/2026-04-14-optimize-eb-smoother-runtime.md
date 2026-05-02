# Optimize `eb_smoother()` Runtime for Matern and Repeated Smoothing

## Summary
- Refactor exact Matern fitting so optimization evaluates only scalar objectives and reconstructs posterior summaries once at the final optimum.
- Add reusable public `Matern_setup()` and thread it through the Matern APIs.
- Reuse cached Matern setup objects in the EBMF post-smoothing helpers.
- Extend the internal runtime study with public-API and repeated-setup benchmarks.

## Key Changes
- In `EBSmoothr/R/02_Matern.R`:
  - add objective-only helpers for profiled exact Matern likelihoods with known and learned noise
  - add a reusable `Matern_setup()` object and setup validation helpers
  - update `ebnm_Matern_generator()` to accept either raw `locations` or a prebuilt `setup`
  - keep exact fixed-parameter evaluation behavior unchanged
  - optimize known-noise and learned-noise exact fits over hyperparameters only, profiling out `beta`
- In `EBSmoothr/R/03_eb_smoother.R`:
  - allow `family = "matern"` to use either `locations` or `setup`
  - reuse cached Matern objects on known-noise and learned-noise paths
- In `attempts_in_EBMF/`:
  - build one `Matern_setup()` per location grid and reuse it inside the component loops
- In `internal/simulations/study_runtime_scaling_1d.R`:
  - add a public-API benchmark at `n = 3000`
  - add a repeated Matern setup-reuse benchmark at `n = 3000`

## Test Plan
- Verify exact objective-only helpers match the full-state exact formulas at fixed hyperparameters.
- Verify `Matern_setup()`-based fits match raw-location fits for exact and INLA-PC Matern paths.
- Keep existing Matern objective-breakdown regression tests passing.
- Run package-level install and `testthat` after regenerating documentation.

## Assumptions
- Preserve the current statistical target for known-noise exact Matern fits.
- Keep `Matern_setup()` as the only new public API added in this runtime pass.
- Keep the existing INLA PC-prior fitting logic unchanged apart from setup reuse.
