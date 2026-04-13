## Update

Added a 1D runtime scaling study across LGP and Matern variants.

## Completed Changes

- Added `internal/simulations/study_runtime_scaling_1d.R`.
- Included LGP, exact Matern, and Matern + PC-prior variants.
- Included both `betaprec = 0` and `betaprec < 0` interpretations.
- Included both the exact/profile and `clinear` implementations for the
  `Matern + PC prior`, `betaprec < 0` case.
- Added method-level execution with incremental output saving.
- Added a single-method CLI mode so large-`n` cases can be benchmarked with an
  external timeout.
- Saved a raw timing CSV, a runtime plot, a markdown summary, and an HTML
  rendering under `internal/simulations/results/`.

## Validation

- `Rscript internal/simulations/study_runtime_scaling_1d.R`
- `Rscript internal/simulations/study_runtime_scaling_1d.R --single-n=100 --single-method=lgp_b0 --timeout-seconds=30`
- `Rscript internal/simulations/study_runtime_scaling_1d.R --single-n=100 --single-method=matern_pc_blt0_clinear --timeout-seconds=30`
- `Rscript internal/simulations/study_runtime_scaling_1d.R --render-from-csv=/Users/ziangzhang/Desktop/EBSmooth/internal/simulations/results/runtime_scaling_1d_methods.csv --timeout-seconds=120`

## Notes

- The study is intended to show rough scaling trends rather than provide a
  high-precision benchmark with many replicates.
- For `n = 10000`, several slow variants were run one method at a time with an
  external `120`-second timeout so the full comparison could be completed
  without losing already-computed results.
- The final outputs are:
  - `internal/simulations/results/runtime_scaling_1d_methods.csv`
  - `internal/simulations/results/runtime_scaling_1d_methods.png`
  - `internal/simulations/results/runtime_scaling_1d_methods_summary.md`
  - `internal/simulations/results/runtime_scaling_1d_methods_summary.html`
