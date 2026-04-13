## Update

Added an internal iid-vs-Matern selection study for the exact Matern backend.

## Completed Changes

- Added `internal/simulations/study_iid_vs_matern_selection.R`.
- Added a staged experiment design with pilot and confirmatory passes.
- Added exact iid-normal comparison via
  `ebnm::ebnm_normal(mode = "estimate", scale = "estimate")`.
- Added per-replicate outputs for:
  - `loglik_matern`, `loglik_iid`, and `delta`;
  - fitted Matern `range`, `sigma`, and `beta`;
  - iid fitted `mean` and `sigma`;
  - the selection label and `range_over_mesh`.
- Added cell-level Wilson intervals, `roughly_half` flags, and
  Matern-selected `range_over_mesh` summaries.
- Added a mesh-sensitivity rerun on a refined mesh for selected near-boundary
  Matern fits.
- Added summary plots for iid selection probability and `delta` versus fitted
  Matern range.

## Validation

- `Rscript internal/simulations/study_iid_vs_matern_selection.R --smoke`
- `Rscript internal/simulations/study_iid_vs_matern_selection.R --pilot_reps=1 --confirmatory_reps=1 --run_confirmatory=false --mesh_sample_size=0 --output_prefix=iid_vs_matern_selection_tinycheck`

## Notes

- The script defaults to the full pilot/confirmatory study but also supports a
  smoke mode via `--smoke`.
- The baseline confirmatory cell is the fixed-domain design with `tau0 = 0.5`
  and `n = 80`.
- The smoke run completed end-to-end and produced the expected CSV, markdown,
  mesh-sensitivity, and figure artifacts under `internal/simulations/results/`.
- The tiny all-cell check confirmed that the full design grid renders both
  regimes and all summary outputs without enabling the expensive confirmatory
  pass.

## Moderate Run

- Ran a moderate staged experiment:
  - `Rscript internal/simulations/study_iid_vs_matern_selection.R --pilot_reps=20 --confirmatory_reps=60 --mesh_sample_size=10 --output_prefix=iid_vs_matern_selection_moderate`
- Main output artifacts:
  - `internal/simulations/results/iid_vs_matern_selection_moderate_draws.csv`
  - `internal/simulations/results/iid_vs_matern_selection_moderate_summary.csv`
  - `internal/simulations/results/iid_vs_matern_selection_moderate_summary.md`
  - `internal/simulations/results/figures/iid_vs_matern_selection_moderate_selection_probability.png`
  - `internal/simulations/results/figures/iid_vs_matern_selection_moderate_delta_vs_range.png`
- Headline findings from the moderate run:
  - No design cell showed strong evidence that Matern wins systematically; no
    cell had `p_iid <= 0.3`.
  - Several small-to-moderate-`n` cells were consistent with `p_iid` near
    one-half.
  - As `n` increased under the weaker-signal setting `tau0 = 0.2`, iid-normal
    was selected much more often, reaching `p_iid = 1.0` at
    increasing-domain `n = 640`.
  - When Matern was selected in many cells, the fitted range was still well
    below the mesh edge, supporting a near-boundary interpretation.
  - A handful of large-`n`, weak-signal Matern fits failed with
    `non-finite finite-difference value`, so the largest-cell results should be
    read as moderate-run evidence rather than a final asymptotic statement.
