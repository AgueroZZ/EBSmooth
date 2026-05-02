# EBSmoothr Point-Exponential Screening Update

## Notes
- Added a plan to support a sparse point-exponential nonspatial reference in
  `eb_smoother()` for constrained nonnegative-loading screening workflows.
- The implementation should avoid unrelated dirty files in the EBSmoothr
  working tree.

## Implementation
- Added `family = "point_exponential"` to `eb_smoother()` with fixed-noise and
  profiled scalar-noise support through `ebnm::ebnm_point_exponential()`.
- Added `matern_n_starts`, defaulting to one start on the TMB Matern
  learned-noise path; `matern_n_starts = 5` preserves the earlier multistart
  behavior.
- Updated the automatic two-dimensional Matern mesh default to use a coarser
  outer mesh while preserving observed locations as mesh vertices.
- Regenerated roxygen documentation.

## Validation
- Parsed `R/02_Matern.R`, `R/03_eb_smoother.R`, and
  `tests/testthat/test-eb-smoother.R`.
- Ran `tests/testthat/test-eb-smoother.R`.
- Ran `tests/testthat/test-matern.R`; INLA emitted retry/segfault diagnostics
  from its external binary, but the test file completed successfully.
- Ran a point-exponential profiled-noise smoke test through `devtools::load_all()`.
