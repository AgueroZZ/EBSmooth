## Update

Completed an internal simulation-based validation pass for the exact Matern
backend.

## Completed Validation

- Added a reproducible internal simulation script:
  `internal/simulations/run_matern_validation.R`
- Saved internal validation outputs under:
  `internal/simulations/results/`
- Ran a one-dimensional hyperparameter-recovery experiment across increasing
  sample sizes.
- Ran one-dimensional and two-dimensional surface-recovery examples and saved
  figures.
- Checked the exact marginal likelihood against a dense Gaussian calculation on
  a small problem.
- Saved a generated markdown summary and HTML rendering of the validation
  results.

## Key Findings

- In a well-specified one-dimensional simulation, the fitted Matern `range` and
  `sigma` move toward the truth as `n` increases, and fitted-surface RMSE
  decreases while correlation increases.
- In a representative two-dimensional example, the fitted posterior mean tracks
  the true simulated surface well.
- The exact marginal likelihood agrees with a dense Gaussian evaluation to
  machine precision on the checked parameter grid.

## Bug Found and Fixed

- The validation uncovered a log-determinant bug in `.compute_logdet_spd()`.
- The helper was fixed before finalizing the validation results.
- After the fix, the dense Gaussian marginal-likelihood check and the package
  test suite both passed.

## Notes

- This validation answers whether the current implementation is internally
  coherent in the exact-model setting; it is not a formal asymptotic proof.
