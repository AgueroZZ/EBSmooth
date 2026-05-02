# Component-Specific Spatial Screening for Matrix Factorization

## Summary
- Replace the old “smooth every fitted loading, then run all-spatial smooth-EBMF” demo with a component-specific routing workflow.
- Use a `2 spatial + 2 nonspatial` simulation as the main matrix-factorization example.
- Screen each fitted loading with `eb_smoother(matern)` versus `eb_smoother(nonspatial)` once after regular EBMF.
- Use that screening result both for selective initialization and for the subsequent mixed smooth-EBMF backfit.

## Planned Changes
- Update `attempts_in_EBMF/simulation_functions.R` with helpers for:
  - screening summaries that include selected route and spatial extent diagnostics;
  - selective smoother initialization that only changes spatial-labeled loadings;
  - factor-specific mixed backfit using `flashier`’s internal factor-level `ebnm.fn` storage.
- Keep nonspatial components on the original loading-side EBNM route during iterative updates.
- Keep the factor-side EBNM unchanged across all components.
- Rewrite `attempts_in_EBMF/simulation.rmd` to show:
  - regular EBMF;
  - aligned screening labels;
  - selective initialization;
  - mixed smooth-EBMF warm start;
  - all-spatial warm-start baseline.
- Add aligned recovery tables and method comparisons by component type.

## Validation
- Smoke-test the new helper workflow on a small simulation.
- Verify the main Rmd setting still yields four fitted components under the chosen seed and signal level.
- Knit the full `simulation.rmd` document once to ensure chunk order and object names are consistent.
