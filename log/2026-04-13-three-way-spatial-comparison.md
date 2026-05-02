## Update

Extended `tryingSpatial.R` so it now supports a clean three-way comparison
between integrated spatial decomposition, post-hoc spatial smoothing, and
fully non-spatial decomposition.

## Completed Changes

- Added `extract_normalized_components_from_matrices()` so non-`flashier`
  factor/loadings matrices can reuse the same alignment and normalization path.
- Added `refit_loadings_from_factors()` for the two-stage workflow where
  smoothed factors are turned back into a full matrix decomposition.
- Added `post_smooth_flash_factors()` to:
  - smooth `standard_fit$fit$F_pm` using `ebnm_Matern_generator()`;
  - use `standard_fit$fit$F_psd` as the uncertainty input;
  - summarize the resulting `rho`, `sigma`, `beta`, and log-likelihood values.
- Added top-level plotting and summaries for:
  - `Spatial Flash`;
  - `Standard Flash`;
  - `Standard Flash + Spatial Post-smoothing`.

## Notes

- The two-stage method is still methodologically different from integrated
  spatial flash. It should be interpreted as a post-processing baseline, not as
  an equivalent objective.
