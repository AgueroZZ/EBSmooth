# Exact `beta` Semantics and `nonspatial` Baseline

## What Changed
- Added `family = "nonspatial"` to `eb_smoother()` as the exact constant-mean Gaussian baseline with no spatial structure.
- Added public `beta_fixed` handling to `eb_smoother()` for both `matern` and `nonspatial`.
- Extended the Matern implementation so fixed `beta` is supported while still optimizing hyperparameters in all main routes:
  - exact known-noise;
  - exact learned-noise;
  - INLA Step A with PC prior;
  - learned-noise INLA Step A with PC prior.
- Kept the exact Gaussian Matern objectives mathematically aligned with “optimize `beta` unless fixed,” while using analytic profiling internally when that is exact and faster.
- Implemented exact `nonspatial` fitting for:
  - known `s`;
  - `s = NULL` with one learned common scalar noise SD.
- Updated the public `eb_smoother_fit` presentation so the wrapper reports the user-facing `beta` semantics (`optimized` versus `fixed`) instead of leaking exact-Matern implementation details such as profiled-beta handling.
- Added targeted tests covering:
  - fixed and optimized `beta` for `nonspatial`;
  - fixed and optimized `beta` for learned-noise `matern`;
  - fixed-beta support for the Matern PC-prior path;
  - wrapper validation and printed summaries.
- Added internal studies for:
  - null-truth `matern` versus `nonspatial` selection under `eb_smoother`;
  - spatial matrix factorization with direct `matern`-versus-`nonspatial` screening and selective smoothing.
- Updated README and simulation index documentation to describe the new baseline family and studies.

## Validation
- Regenerated the relevant `.Rd` files from roxygen comments.
- Ran targeted package tests:
  - `EBSmoothr/tests/testthat/test-matern.R`
  - `EBSmoothr/tests/testthat/test-eb-smoother.R`
- Smoke-tested the new internal study scripts:
  - `internal/simulations/study_nonspatial_vs_matern_null.R --smoke`
  - `internal/simulations/study_mf_selective_smoothing.R --smoke`
