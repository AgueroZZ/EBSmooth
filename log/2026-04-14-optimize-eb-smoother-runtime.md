# 2026-04-14 Optimize `eb_smoother()` Runtime

## Summary
- Refactored the exact Matern implementation so optimization no longer reconstructs posterior summaries at every objective evaluation.
- Added reusable public `Matern_setup()` and threaded it through both `ebnm_Matern_generator()` and `eb_smoother()`.
- Updated EBMF-side post-smoothing helpers to reuse one cached Matern setup per location grid.
- Extended the internal runtime study with public single-fit and repeated-setup benchmarks.

## Implementation Notes
- Added internal profiled objective helpers for exact Matern fits with known noise and learned common noise.
- Added a small optimizer fallback from `BFGS` to `Nelder-Mead` for exact Matern paths when the BFGS endpoint lands in the finite-penalty region.
- Preserved exact fixed-parameter evaluation semantics and Matern objective-breakdown behavior.
- Regenerated roxygen outputs, including:
  - `EBSmoothr/man/Matern_setup.Rd`
  - updated `EBSmoothr/man/ebnm_Matern_generator.Rd`
  - updated `EBSmoothr/man/eb_smoother.Rd`

## Validation
- `R CMD INSTALL EBSmoothr`
- `Rscript -e 'library(EBSmoothr); testthat::test_dir("EBSmoothr/tests/testthat", reporter = "summary")'`
- `Rscript internal/simulations/study_runtime_scaling_1d.R --timeout-seconds=60`
- `Rscript internal/simulations/study_runtime_scaling_1d.R --render-from-csv=/Users/ziangzhang/Desktop/EBSmooth/internal/simulations/results/runtime_scaling_1d_methods.csv --timeout-seconds=60`

## Benchmark Highlights
- Public API benchmark at `n = 3000`:
  - `matern_exact_known`: `2.553s`
  - `matern_exact_learned`: `2.485s`
  - `matern_pc_known`: `3.430s`
  - `matern_pc_learned`: `3.864s`
- Repeated learned-noise Matern PC smoothing over 4 components at `n = 3000`:
  - rebuild from `locations` each call: `12.225s`
  - reuse one `Matern_setup()`: `11.273s`

## Output Files
- `internal/simulations/results/runtime_scaling_1d_methods.csv`
- `internal/simulations/results/runtime_scaling_1d_methods_summary.md`
- `internal/simulations/results/runtime_public_api_n3000.csv`
- `internal/simulations/results/runtime_repeated_matern_n3000.csv`
