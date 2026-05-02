# Matern Auto Backend Policy

## Changes

- Added `.matern_auto_backend()` and `.matern_resolve_backend()` so
  `ebnm_Matern_generator()` and `eb_smoother()` share the same backend policy.
- Updated `backend = "auto"` for Matern:
  - identity link: `exact`;
  - log link: package `laplace` by default;
  - log-link, 2D, fixed-beta, unpenalized fits with optimized Matern
    hyperparameters: `inla`.
- Left explicit backend choices unchanged.
- Kept log-link INLA comparable likelihood evaluation warm-started from the INLA
  spatial mode in both the generator and wrapper routes.
- Updated Rd text to describe the conservative auto policy.

## Tests Added

- Generator auto dispatch checks for 2D log-link fixed beta and empirical-Bayes
  beta.
- Wrapper auto dispatch checks matching the generator policy.

## Validation Run

- `R CMD INSTALL EBSmoothr`
- `Rscript -e 'library(EBSmoothr); testthat::test_file("EBSmoothr/tests/testthat/test-matern.R", reporter = "summary")'`
- `Rscript -e 'library(EBSmoothr); testthat::test_file("EBSmoothr/tests/testthat/test-eb-smoother.R", reporter = "summary")'`
- `Rscript -e 'library(EBSmoothr); testthat::test_file("EBSmoothr/tests/testthat/test-lgp.R", reporter = "summary")'`
- `Rscript -e 'tools::checkRd("EBSmoothr/man/ebnm_Matern_generator.Rd"); tools::checkRd("EBSmoothr/man/eb_smoother.Rd")'`
- `git diff --check` on the touched tracked files.

Small 2D log-link smoke comparison, TMB Laplace minus INLA:

- fixed beta log-likelihood difference: `6.6084e-07`
- fixed beta max posterior mean difference: `0.00853387`
- empirical-Bayes beta difference: `-0.00152751`
- empirical-Bayes log-likelihood difference: `0.000585506`
- empirical-Bayes max posterior mean difference: `0.00703349`

## Rationale

The low-risk performance route is to change only `auto` dispatch for the single
2D case where benchmarks show INLA is clearly faster while matching the package
Laplace objective. Other beta modes continue to default to TMB Laplace so the
package-owned implementation remains the default correctness anchor.
