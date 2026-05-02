# `eb_smoother_fit` Overview and Documentation Update

## What Changed
- Added S3 overview methods for `eb_smoother_fit`:
  - `print.eb_smoother_fit()`
  - `summary.eb_smoother_fit()`
  - `print.summary.eb_smoother_fit()`
- Kept the new methods presentation-only. They reorganize stored fit contents and do not recompute expensive model quantities.
- Added compact posterior summaries based on stored posterior means and variances.
- Added conditional diagnostics reporting for stored backend-specific fields only.
- Updated the public quick-start documentation in `README.md` so users can see the intended `fit` / `summary(fit)` workflow.
- Updated tests to cover the new overview methods across:
  - Matern exact known noise
  - Matern exact learned noise
  - Matern INLA-PC learned noise
  - LGP known noise
  - LGP learned noise

## Validation
- Regenerated package namespace and manual pages from roxygen comments.
- Installed the package successfully to a temporary library.
- Ran the targeted `EBSmoothr/tests/testthat/test-eb-smoother.R` test file to validate:
  - S3 dispatch
  - compact printed overview content
  - summary object class and printed sections
  - optional diagnostics handling
