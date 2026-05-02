# Release-Readiness Backend Consistency and Auto Policy

## Implemented

- Added learned-noise awareness to Matern backend auto routing. Log-link
  learned-noise Matern fits now stay on the package Laplace/TMB path under
  `backend = "auto"`.
- Kept the verified auto-INLA path limited to unpenalized 2D log-link
  fixed-beta known-noise Matern fits.
- Removed the unconditional explicit-INLA rejection for log-link learned-noise
  no-PC Matern fits.
- Added validate-or-error behavior for explicit log-link learned-noise INLA
  fits: the returned INLA mode is checked against the package Laplace reference
  for objective, beta, Matern hyperparameters, learned noise SD, and posterior
  means. A failed validation raises an error instead of returning a
  non-comparable fit.
- Updated Matern and `eb_smoother()` documentation to state backend-comparable
  objective semantics: exact objective for identity-link INLA modes and package
  Laplace objective for log-link INLA modes.
- Added a release-readiness benchmark script:
  `internal/simulations/release_readiness_backend_check.R`.

## Validation

- `R CMD INSTALL EBSmoothr`: passed.
- Targeted tests passed:
  - `EBSmoothr/tests/testthat/test-matern.R`
  - `EBSmoothr/tests/testthat/test-eb-smoother.R`
  - `EBSmoothr/tests/testthat/test-lgp.R`
- `R CMD build EBSmoothr`: passed and produced `EBSmoothr_0.1.0.tar.gz`.
- `R CMD check --no-manual --no-build-vignettes --no-tests EBSmoothr_0.1.0.tar.gz`:
  passed with one existing NOTE about the quiet-INLA helper temporarily
  unlocking INLA namespace bindings.
- A direct `R CMD check` on the source directory was stopped during full
  testthat after it ran much longer than the targeted suite. Before tests, it
  showed source-directory-only warnings for compiled objects and hidden project
  files; the tarball check confirmed `.Rbuildignore` removes those.

## Benchmark Results

- Accuracy results were written to
  `internal/simulations/results/backend_release_accuracy.csv`.
- Runtime results were written to
  `internal/simulations/results/backend_release_runtime.csv`.
- Markdown summary was written to
  `internal/simulations/results/backend_release_summary.md`.
- Accuracy matrix:
  - identity known-noise 1D EB: exact, TMB Laplace, and INLA agreed; INLA
    objective delta was about `3e-5`.
  - log known-noise 1D EB: TMB Laplace, R Laplace, and INLA agreed; INLA
    objective delta was about `2e-5`.
  - log learned-noise 1D EB: TMB Laplace and R Laplace agreed; explicit INLA
    no-PC failed validation with objective delta about `0.11`, so it correctly
    errored instead of returning the fit.
  - log learned-noise 1D EB with PC prior: TMB Laplace, R Laplace, and INLA-PC
    were comparable; INLA-PC passed validation with objective delta about
    `0.072` and posterior-mean max delta about `0.020`.
  - log known-noise 2D fixed beta: TMB Laplace, explicit INLA, and auto-INLA
    agreed with objective deltas below `1e-6` in the small accuracy case.
- Runtime matrix:
  - 1D log learned-noise EB: `auto`/TMB took about `0.48s`, `0.89s`, and
    `2.45s` for `n = 100, 300, 1000`; R Laplace took about `30s` at
    `n = 100`.
  - 2D log known-noise fixed beta: INLA was faster than TMB at every tested
    grid size: about `5.4s` vs `13.0s` at `16x16`, `7.0s` vs `13.7s` at
    `24x24`, and `10.6s` vs `27.6s` at `32x32`.

## Final Auto Policy

- Identity-link Matern: `exact`.
- Nonspatial: `exact`.
- Log-link Matern learned-noise: `laplace`/TMB.
- Log-link Matern known-noise 2D fixed beta without PC prior: `inla`.
- R Laplace remains a reference/debug backend and is not selected by auto when
  TMB is supported.
