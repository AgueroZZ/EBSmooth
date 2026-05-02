# Simplify EBSmoothr Matern Backend API

## Notes
- Added a plan to simplify Matern backend names while preserving compatibility
  aliases for existing code.
- `pc.penalty` is the single source of truth for whether INLA uses PC-prior
  penalization.

## Implementation
- Added backend argument matching/canonicalization so `laplace_tmb` maps to
  `laplace`, and `inla_pc` maps to `inla` when `pc.penalty` is supplied.
- Updated public backend defaults and documentation to list the simplified
  backend set.
- Updated INLA-backed fits to report `fit$backend = "inla"` regardless of
  whether `pc.penalty` is supplied.
- Added targeted tests for Laplace implementation reporting and compatibility
  aliases.

## Validation
- Parsed `R/02_Matern.R`, `R/03_eb_smoother.R`,
  `tests/testthat/test-eb-smoother.R`, and `tests/testthat/test-matern.R`.
- Ran `tests/testthat/test-eb-smoother.R`.
- Ran `tests/testthat/test-matern.R`; INLA emitted external binary retry and
  segmentation-fault diagnostics, but the test file completed successfully.
- Parsed SpatialEBMF `code/simulation_functions.R`.
