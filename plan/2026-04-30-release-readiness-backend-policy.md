# Release-Readiness Backend Consistency and Auto Policy

## Summary

- Make backend choice a runtime implementation detail: a successful backend
  must report the same package-level objective semantics as the reference path.
- Keep identity-link Matern INLA fits on the package exact Gaussian objective,
  evaluated at the INLA mode.
- Keep log-link Matern INLA fits on the package Laplace objective, evaluated at
  the INLA mode.
- Validate difficult explicit INLA log-link learned-noise fits against the
  package Laplace reference and error if validation fails.

## Implementation

- Pass learned-noise state into Matern auto backend resolution so
  `backend = "auto"` does not choose INLA for log-link learned-noise fits.
- Keep auto INLA routing only for the verified unpenalized 2D fixed-beta
  known-noise log-link path.
- Add INLA-vs-Laplace validation diagnostics for explicit log-link
  learned-noise Matern INLA fits.
- Update documentation and internal notes to describe comparable objective
  semantics and validation behavior.

## Validation

- Add regression tests for log-link learned-noise INLA validation behavior.
- Add regression tests that `eb_smoother()` auto routing keeps 2D log-link
  learned-noise fixed-beta fits on the Laplace/TMB path.
- Run targeted Matern, `eb_smoother()`, and LGP tests after installation.
- Run a focused release-readiness benchmark matrix covering 1D/2D, known versus
  learned noise, fixed versus empirical-Bayes beta, and representative
  backends.

## Acceptance

- Explicit backends either return comparable fitted objectives and parameters or
  fail loudly with a validation error.
- `backend = "auto"` uses accuracy-validated paths first and runtime evidence
  second.
- Release notes report backend consistency status, runtime comparisons, and the
  resulting auto backend decision table.
