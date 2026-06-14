# FPCA Weather Exploration

## Input

- Supplement zip: `/Users/ziangzhang/Desktop/EBSmooth/EBFPCA/data/biom12457-sup-0002-suppdatacode.zip`
- Data member: `iFPCAcodes/fdaM/examples/weather/temperature.csv`
- Grid points: 365
- Curves: 35

## Methods

- `raw_grid_prcomp`: ok; Base R prcomp on centered daily grid values.
- `fda_pca_fd_bspline`: ok; fda::pca.fd with bspline basis; nbasis=35, lambda=0.01.
- `fda_pca_fd_fourier`: ok; fda::pca.fd with fourier basis; nbasis=35, lambda=0.01.
- `fdapace_FPCA`: ok; fdapace::FPCA on dense list-form trajectories.
- `refund_fpca_sc`: ok; refund::fpca.sc on dense matrix trajectories.

## Output Files

- `metrics.csv`: summary metrics for fitted and skipped methods.
- `fits.rds`: fitted objects and reconstructed curves.
- `weather_curves.png`: raw weather curves plus mean curve.
- `variance_explained.png`: variance proportions by method.
- `harmonics.png`: estimated eigenfunctions.
- `scores_pc1_pc2.png`: PC1/PC2 score scatter, where available.
- `session_info.txt`: R session metadata.
