## Update

Separated the expensive exact-objective validation from the main Matern
PC-prior fitting path and added a post-fit objective inspection helper.

## Completed Changes

- Added `matern_objective_breakdown()` to recompute exact Matern objective
  components after fitting.
- Added internal helpers for exact fixed-beta, profiled-beta, and
  integrated-flat objective calculations.
- Added `compute_exact_diagnostic` to `ebnm_Matern_generator()`.
- Changed the default PC-prior Step A path so it no longer automatically
  computes `log_likelihood_exact_at_stepA_mode`.
- Stored Matern objective context and beta-mode metadata in fitted objects.
- Updated tests and internal documentation to use the post-fit inspection
  workflow.

## Validation

- `roxygen2::roxygenize("EBSmoothr")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-matern.R")`

## Notes

- This change is primarily about runtime clarity: fitting and validation are now
  separate operations.
- The earlier runtime benchmark for `Matern + PC prior`, `betaprec = 0` should
  now be interpreted as including an extra exact-objective cross-check that is
  no longer part of the default fit path.
