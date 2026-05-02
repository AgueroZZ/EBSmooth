# Exact `beta` Semantics and `nonspatial` Baseline

## Summary
- Add `family = "nonspatial"` to `eb_smoother()` as the no-spatial-information Gaussian baseline.
- Make the public `beta` semantics parallel for `matern` and `nonspatial`:
  - default: optimize `beta` by marginal likelihood;
  - optional: hold `beta` fixed through `beta_fixed`.
- Keep analytic profiling of `beta` as an implementation detail for exact Gaussian objectives whenever it is mathematically exact and faster.
- Add internal studies for:
  - null-truth `matern` versus `nonspatial` model selection;
  - matrix factorization with `2` spatial and `2` nonspatial components plus selective smoothing.

## Planned Changes
- Extend `eb_smoother()` with:
  - `family = "nonspatial"`;
  - public `beta_fixed = NULL` for `matern` and `nonspatial`.
- Update the exact Matern code paths so fixed `beta` is allowed while Matern hyperparameters are still being optimized.
- Update the learned-noise Matern paths so fixed `beta` works for both:
  - exact Gaussian fitting;
  - INLA Step A with PC prior.
- Implement the `nonspatial` family for both:
  - known `s`;
  - `s = NULL` with one learned common scalar noise SD.
- Keep the public semantics focused on fitted versus fixed `beta`, while allowing exact Gaussian code to profile `beta` internally.
- Add targeted tests for:
  - fixed and optimized `beta` in the Matern paths;
  - fixed and optimized `beta` in the `nonspatial` paths;
  - public wrapper validation and printed summaries.
- Add new internal studies:
  - `study_nonspatial_vs_matern_null.R`;
  - `study_mf_selective_smoothing.R`.
- Update README and package help pages so the public API describes the new `nonspatial` family and the fitted/fixed intercept semantics.

## Validation
- Regenerate package documentation from roxygen comments.
- Run targeted package tests for:
  - `test-matern.R`
  - `test-eb-smoother.R`
- Smoke-run the new internal studies and confirm they emit summary artifacts.
