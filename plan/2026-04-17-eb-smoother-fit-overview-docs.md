# Add Documentation and Overview Methods for `eb_smoother_fit`

## Summary
- Add compact S3 overview methods for `eb_smoother_fit`.
- Document how to inspect fitted objects through `print()` and `summary()`.
- Keep the change presentation-only: no statistical or structural changes to stored fit contents.

## Planned Changes
- Add `print.eb_smoother_fit()`, `summary.eb_smoother_fit()`, and `print.summary.eb_smoother_fit()` in the `eb_smoother` implementation file.
- Keep `print()` compact and one-screen, showing:
  - family, backend, and noise mode
  - number of observations
  - fitted prior parameters
  - fitted beta
  - fitted noise SD when available
  - log-likelihood
- Make `summary()` return a structured summary object that includes:
  - the same high-level metadata
  - posterior mean/SD quantiles derived from the stored posterior
  - backend-specific diagnostics only when already stored in the fit
- Update public docs:
  - `README.md`
  - `eb_smoother()` help
  - help for the new `eb_smoother_fit` overview methods
- Add targeted tests for dispatch, printed content, and optional diagnostics across the main fit modes.

## Validation
- Regenerate namespace/man entries from roxygen comments.
- Install the package to a temporary library.
- Run the targeted `test-eb-smoother.R` test file.
