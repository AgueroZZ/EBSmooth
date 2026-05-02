# Add documented spatial scoring APIs

## Goal

Add public wrappers that let users score fitted factor/loadings components under
Matern spatial smoothers and nonspatial reference families without copying code
from the simulation notebook.

## Scope

- Add `spatial_scores()`, `spatial_select()`, and
  `spatial_scores_permutation()`.
- Rename the high-level Gaussian baseline from `family = "nonspatial"` to
  `family = "constant"` and export `Constant()`.
- Keep `"nonspatial"` as an intentional error so users update code explicitly.
- Support Matern spatial specs with `link = "identity"` and `link = "log"`.
- Support reference families `constant`, `point_exponential`, `point_normal`,
  and `point_laplace`.
- Refactor high-level point-family fitting so all point families share
  fixed-noise and profiled-noise support.
- Add high-level Matern `fix_g` support for permutation scoring with fixed
  observed Matern parameters.
- Regenerate roxygen docs and add targeted tests.

## Implementation Plan

1. Refactor `eb_smoother()` point-family code into shared helpers and update
   validation for point-mass reference families.
2. Replace public `Nonspatial()` with `Constant()` and update user-facing
   score semantics and docs from `nonspatial` to `constant`.
3. Thread `fix_g` through high-level Matern `eb_smoother()` calls, including
   exact and Laplace learned-noise paths.
4. Add `R/04_spatial_scores.R` with extraction helpers, spec validation,
   scoring, selection, and permutation p-value logic.
5. Update the nonnegative log-Matern simulation notebook to call
   `spatial_select()` with a log-link Matern model and point-exponential
   reference.
6. Add tests for rename behavior, point-family support, spatial scoring,
   selection, and deterministic permutation scoring.
7. Regenerate documentation and run targeted parse/Rd and testthat checks.
