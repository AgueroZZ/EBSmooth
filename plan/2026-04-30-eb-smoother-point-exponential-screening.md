# EBSmoothr Point-Exponential Screening Update

## Summary
Add a sparse nonspatial point-exponential family to `eb_smoother()` and make
the default Matern screening path faster for learned-noise log-link fits.

## Implementation
- Add `family = "point_exponential"` to `eb_smoother()` with fixed-noise and
  profiled-noise support through `ebnm::ebnm_point_exponential()`.
- Add `matern_n_starts`, defaulting to one start, while preserving the previous
  five-start behavior when explicitly requested.
- Use a coarser automatic outer mesh for two-dimensional Matern setup when
  `max.edge = NULL`, while preserving observed locations as mesh vertices.
- Update public documentation and tests for the new screening reference.

## Validation
- Parse package R sources.
- Run targeted `testthat` coverage for point-exponential fits, Matern start
  counts, and default mesh vertex inclusion.
