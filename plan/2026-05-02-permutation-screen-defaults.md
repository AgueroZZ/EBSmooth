# Default to fixed-parameter permutation screening

## Goal

Make permutation screening the preferred spatial-selection route and align the
permutation p-value with a direct observed-score-versus-permuted-scores
comparison.

## Planned Changes

- Change `spatial_scores_permutation()` default `refit` from `TRUE` to
  `FALSE`.
- Keep refit support available for users who explicitly want to re-estimate
  Matern parameters on each permuted data set.
- Compute `p_value` as the fraction of valid permutation scores greater than or
  equal to the observed spatial score, without adding the observed score to the
  permutation-null sample.
- Propagate observed Matern fits from `spatial_select(method = "permutation")`
  when `keep_fits = TRUE`, so downstream smoothing workflows can reuse the
  selected spatial fits.
- Update tests and documentation for the default and p-value semantics.
