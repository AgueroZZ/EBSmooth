# 2026-04-13 Implement `eb_smoother()`

## Summary
- Added a new public `eb_smoother()` API for both known observation standard errors and learned common noise SDs.
- Preserved `ebnm_LGP_generator()` and `ebnm_Matern_generator()` as known-SE, `ebnm`-compatible low-level interfaces.
- Added learned-noise fitting paths for exact Matern, INLA-PC Matern, and TMB-backed L-GP.

## Code Changes
- Added `EBSmoothr/R/03_eb_smoother.R` with:
  - argument normalization
  - fit-object wrapping
  - learned-noise L-GP fitting
  - exported `eb_smoother()` dispatch
- Extended `EBSmoothr/R/02_Matern.R` with:
  - learned-noise exact Matern state evaluation
  - learned-noise INLA-PC Step A fitting
  - optional `pc.penalty$noise` handling for the new API
- Extended `EBSmoothr/src/EBSmoothr.cpp` and `EBSmoothr/R/01_LGP.R` so the TMB objective can learn one common observation noise SD while keeping the old known-`s` path fixed.
- Exported `eb_smoother()` in `EBSmoothr/NAMESPACE` and added `EBSmoothr/man/eb_smoother.Rd`.
- Updated the EBMF experiment helpers in `attempts_in_EBMF/` to use `eb_smoother()` for posterior-mean post-smoothing.

## Validation
- Reinstalled the package with `R CMD INSTALL EBSmoothr`.
- Ran the package tests with:
  - `Rscript -e 'library(EBSmoothr); testthat::test_dir("EBSmoothr/tests/testthat", reporter = "summary")'`
- Added regression tests covering:
  - known-noise Matern wrapping
  - known-noise L-GP wrapping
  - learned-noise exact Matern
  - learned-noise INLA-PC Matern with `pc.penalty$noise`
  - learned-noise L-GP
  - API validation for unsupported combinations
