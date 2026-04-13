## Update

Implemented an optional `pc.penalty` mode for the Matern backend while keeping
the default exact Gaussian empirical-Bayes implementation intact.

## Completed Changes

- Added a generator-level `pc.penalty` argument to
  `ebnm_Matern_generator()`.
- Added parsing and validation helpers for `pc.penalty`, including defaults for
  omitted `range` or `sigma` anchors and default tail probability `0.5`.
- Kept the default exact Matern branch unchanged when `pc.penalty = NULL`.
- Added an INLA Step A backend for `pc.penalty != NULL` and `fix_g = FALSE`:
  - uses `inla.spde2.pcmatern()`;
  - returns posterior summaries directly from Step A;
  - uses `resA$misc$log.posterior.mode` as the main `log_likelihood`.
- Added an INLA Step B backend for `pc.penalty != NULL` and `fix_g = TRUE`:
  - respects fixed `g_init` hyperparameters and fixed `beta_fixed`;
  - returns a penalized pseudo-objective equal to Step B plus the PC prior
    log-density on `(log(range), log(sigma))`.
- Added backend metadata and diagnostics including:
  - `backend`;
  - `prior_family`;
  - `pc_penalty`;
  - Step A penalized and marginal-likelihood diagnostics;
  - Step B conditional marginal likelihood and PC-prior contribution.
- Left `posterior_sampler` available only in exact mode; the INLA PC-prior mode
  now returns an explicit warning and an `NA` matrix as a compatibility stub.
- Added regression tests for:
  - `pc.penalty` parsing;
  - INLA Step A output semantics;
  - fixed-parameter penalized pseudo-objective semantics;
  - Step B prior invariance when `theta` is fixed;
  - agreement between exact and weak-PC fits on a small 2D example.
- Updated internal math notes and internal vignettes to explain the dual
  backend structure and the backend-dependent interpretation of
  `log_likelihood`.

## Validation

- `Rscript -e 'library(testthat); library(pkgload); pkgload::load_all("EBSmoothr", quiet = TRUE, export_all = TRUE); test_file("EBSmoothr/tests/testthat/test-matern.R")'`
- `Rscript -e 'library(testthat); library(pkgload); pkgload::load_all("EBSmoothr", quiet = TRUE, export_all = TRUE); test_dir("EBSmoothr/tests/testthat", reporter = "summary")'`
- `Rscript -e 'roxygen2::roxygenize("EBSmoothr")'`
- `Rscript internal/vignettes/render_internal_vignettes.R`
- `pandoc internal/math/ebnm-smooth-foundations.md -s -o internal/math/ebnm-smooth-foundations.html`
- `Rscript EBSmoothr/sanity_check_matern.R`
- `R CMD check EBSmoothr --no-manual`

## Notes

- In exact mode, `log_likelihood` remains the exact Gaussian EB marginal
  likelihood.
- In `pc.penalty` mode, `log_likelihood` is backend-dependent by design and is
  now explicitly labeled through additional diagnostic fields.
- `R CMD check` still reports two pre-existing NOTES about hidden files and the
  `LICENSE` stub when checking the source directory directly.
