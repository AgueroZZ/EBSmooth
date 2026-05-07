# Inlabru prior and likelihood semantics cleanup

## Summary

Updated the experimental softplus Matern `backend = "inlabru"` policy so it no
longer creates hidden data-driven PC priors:

- known-noise inlabru fits keep `pc.penalty = NULL` as an unpenalized,
  non-PC-SPDE fit;
- learned-noise inlabru fits now require an explicit `pc.penalty` list with
  `range`, `sigma`, and `noise` entries;
- supplied inlabru PC priors must include all required components, avoiding
  implicit range/sigma/noise prior synthesis.

The primary inlabru `log_likelihood` is now the package's observed-Hessian
manual Laplace objective evaluated at the inlabru fitted hyperparameters. Raw
inlabru/INLA marginal likelihood is retained as
`log_likelihood_inlabru_mlik_integration`.

## Validation targets

- `pkgload::load_all("EBSmoothr", export_all = TRUE)`
- `testthat::test_file("EBSmoothr/tests/testthat/test-matern.R")`
- a small known-noise no-PC softplus smoke comparison between `backend =
  "laplace"` and `backend = "inlabru"`

## Validation results

- `pkgload::load_all("EBSmoothr", export_all = TRUE)` passed after applying
  the changes to the current branch.
- `devtools::document("EBSmoothr")` passed and kept the generated Rd docs in
  sync with the roxygen source after applying the changes to the current branch.
- `testthat::test_file("EBSmoothr/tests/testthat/test-matern.R")` passed after
  applying the changes to the current branch. INLA printed optimizer/segfault
  recovery diagnostics during the run, but the testthat process completed
  successfully.
- Known-noise no-PC smoke comparison (`n = 50`) passed:
  - `pc_penalty_is_null = TRUE`
  - `log_likelihood_semantics = "laplace_at_inlabru_params_empirical_bayes"`
  - max posterior mean difference versus manual Laplace: `0.005927582`
  - max posterior variance difference versus manual Laplace: `0.0002824525`
  - primary comparable likelihood minus raw inlabru mlik: `10.41908`
  - elapsed times: manual Laplace `0.31s`, inlabru `10.416s`
- A final policy smoke check passed for known-noise `pc.penalty = NULL`,
  learned-noise missing-`pc.penalty` rejection, and explicit NULL prior-entry
  rejection.
