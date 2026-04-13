## Update

Completed the exact empirical-Bayes rewrite of the Matern backend and the
fixed-beta cleanup for both smoothing backends.

## Completed Changes

- Added `beta_fixed` support to the returned fitter from
  `ebnm_LGP_generator()`.
- Changed L-GP `fix_g = TRUE` semantics so that:
  - `theta` is fixed at `g_init$scale`;
  - `beta` is fixed at `beta_fixed`;
  - `beta_fixed` defaults to zero;
  - Step B only infers `U`.
- Replaced the old Matern `INLA::inla()` Step A / Step B workflow with an exact
  Gaussian marginal-likelihood implementation based on:
  - `INLA::inla.spde2.matern()`;
  - `INLA::inla.spde2.precision()`;
  - outer optimization over `(log(range), log(sigma), beta0)`.
- Extended `Matern()` to store both `theta = log(range)` and `sigma`.
- Restricted `ebnm_Matern_generator()` to `link = "identity"` only.
- Changed Matern `log_likelihood` to return the exact empirical-Bayes marginal
  likelihood.
- Changed Matern `inla_result` to `NULL` as a compatibility placeholder.
- Added exact Gaussian `posterior_spatial_field` summaries for the latent mesh
  field.
- Added `testthat` infrastructure and regression tests for:
  - L-GP fixed-beta behavior;
  - L-GP likelihood diagnostics;
  - Matern identity-only behavior;
  - Matern exact fixed and optimized fits in 1D and 2D.
- Updated package documentation, internal math notes, internal vignettes, and
  rendered HTML outputs.
- Updated `sanity_check_matern.R` to use the new exact identity-link workflow.

## Validation

- `Rscript -e 'pkgload::load_all("EBSmoothr", quiet = TRUE)'`
- `Rscript -e 'pkgload::load_all("EBSmoothr", quiet = TRUE); testthat::test_dir("EBSmoothr/tests/testthat", reporter = "summary")'`
- `Rscript internal/vignettes/render_internal_vignettes.R`
- `Rscript EBSmoothr/sanity_check_matern.R`
- `pandoc internal/math/ebnm-smooth-foundations.md -s -o internal/math/ebnm-smooth-foundations.html`

## Notes

- The Matern rewrite is motivated by the goal of matching the empirical-Bayes
  objective more closely than the previous `INLA::inla()` hyperprior workflow.
- The implementation uses SPDE precision matrices available in the current INLA
  installation and does not depend on the local `rSPDE` runtime.
