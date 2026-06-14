# EB-FPCA Weather Exploration Script

## Summary

Added a first R script for exploring FPCA baselines on the Lin et al. iFPCA
supplement weather dataset.

## Files Added or Updated

- `EBFPCA/scripts/explore_fpca_weather.R`: reads the supplement zip directly,
  fits several FPCA baselines, and writes plots plus metrics.
- `EBFPCA/results/fpca_weather_exploration/`: generated first-run output files.
- `EBFPCA/knowledge.md`: added local supplement-data status and first-run
  exploration notes.

## Data Source

The script reads:

- zip: `EBFPCA/data/biom12457-sup-0002-suppdatacode.zip`
- member: `iFPCAcodes/fdaM/examples/weather/temperature.csv`

The data are 365 daily temperature grid points for 35 Canadian weather
locations.

## Methods in the Script

Always attempted:

- `raw_grid_prcomp`: base R PCA on centered grid values.
- `fda_pca_fd_bspline`: `fda::pca.fd` after B-spline smoothing.
- `fda_pca_fd_fourier`: `fda::pca.fd` after Fourier smoothing.

Optional, skipped if packages are unavailable:

- `fdapace_FPCA`: `fdapace::FPCA`.
- `refund_fpca_sc`: `refund::fpca.sc`.

## First Run

Command:

```sh
Rscript --vanilla EBFPCA/scripts/explore_fpca_weather.R
```

Observed on this machine:

- `fda` was installed and ran successfully.
- `fdapace` was not installed and was skipped.
- `refund` was not installed and was skipped.
- The successful methods gave very similar first three components.
- PC1 explained about 88 percent of variation, PC2 about 8.5 percent, and PC3
  about 2 percent.
- The B-spline and Fourier eigenfunctions had absolute correlations above 0.96
  with the corresponding raw-grid PCA eigenfunctions for the first three PCs.

## Outputs

Generated under `EBFPCA/results/fpca_weather_exploration/`:

- `metrics.csv`
- `fits.rds`
- `run_summary.md`
- `weather_curves.png`
- `variance_explained.png`
- `harmonics.png`
- `scores_pc1_pc2.png`
- `session_info.txt`

## Follow-Up

- Install `fdapace` and rerun the same script to add the PACE baseline.
- Install `refund` and rerun if the smoothed covariance baseline is needed.
- Add a separate MATLAB or Octave wrapper only if exact reproduction of the
  original iFPCA supplement method is required.

## 2026-06-10 Update: PACE, Refund, and iFPCA Replication

Installed legitimate CRAN packages:

- `fdapace`
- `refund`

Reran:

```sh
Rscript --vanilla EBFPCA/scripts/explore_fpca_weather.R
```

The updated exploration now includes successful `fdapace::FPCA` and
`refund::fpca.sc` baselines. PACE is nearly identical to raw grid PCA on this
dense common-grid weather dataset; its first three variance proportions are
about 89.3 percent, 8.6 percent, and 2.1 percent.

Added:

- `EBFPCA/scripts/replicate_lin_ifpca_weather.R`

This script is an R translation of the Lin et al. supplement weather workflow.
It implements:

- B-spline basis coefficients for centered weather curves;
- regularized FPCA through a roughness-penalized generalized eigenproblem;
- greedy L0 support selection for interpretable FPCA;
- 10-fold CV kappa selection;
- PACE overlay for comparison.

Main command:

```sh
Rscript --vanilla EBFPCA/scripts/replicate_lin_ifpca_weather.R
```

Main output directory:

- `EBFPCA/results/lin_ifpca_weather_replication/`

Generated files include:

- `lin_ifpca_weather_eigenfunctions.png`: iFPCA, regularized FPCA, and PACE
  overlay.
- `lin_ifpca_weather_eigenfunctions_paper_style.png`: iFPCA vs regularized
  FPCA plot closer to the supplement weather figure.
- `lin_ifpca_weather_cumulative_variance.png`
- `lin_ifpca_weather_metrics.csv`
- `lin_ifpca_weather_support.csv`
- `lin_ifpca_weather_fits.rds`

The current recommended run uses `nbasis = 181`, `gamma = 3000`, and 10 CV
folds. This is an approximate R replication rather than a bit-for-bit MATLAB
rerun. The qualitative result matches the paper's purpose: iFPCA creates
localized eigenfunctions with explicit inactive intervals, while PACE and
regularized FPCA produce smooth global eigenfunctions.

Selected iFPCA active basis counts:

- PC1: 103
- PC2: 60
- PC3: 13

Detected iFPCA active support intervals:

- PC1: day 1.5-128.5 and day 271.5-363.5
- PC2: day 123.5-181.5 and day 195.5-275.5
- PC3: day 95.5-126.5
