# Lin iFPCA Weather Replication

## Notes

This is an R translation of the MATLAB supplement workflow. It uses the
same weather data, greedy L0 support-selection algorithm, CV rule, and
roughness penalty idea, but the R `fda` basis setup is not bit-for-bit
identical to the original MATLAB FDA toolbox.

## Settings

- B-spline nbasis: 181
- Number of components: 3
- Gamma: 3000
- CV folds: 10
- Kappa grid: 1 to 181 by 1

## Selected Kappas

- PC1: selected kappa 103 (minimum-CV kappa 181)
- PC2: selected kappa 60 (minimum-CV kappa 158)
- PC3: selected kappa 13 (minimum-CV kappa 162)

## Variance Proportions

 component      ifpca regularized_fpca       pace selected_kappa min_cv_kappa
       PC1 0.77980935       0.88636854 0.89323640            103          181
       PC2 0.12431001       0.08515179 0.08588688             60          158
       PC3 0.01554297       0.02038850 0.02087673             13          162

## iFPCA Active Support Intervals

 component start_day end_day
       PC1       1.5   128.5
       PC1     271.5   363.5
       PC2     123.5   181.5
       PC2     195.5   275.5
       PC3      95.5   126.5

## Output Files

- `lin_ifpca_weather_raw_curves.png`
- `lin_ifpca_weather_eigenfunctions.png`
- `lin_ifpca_weather_eigenfunctions_paper_style.png`
- `lin_ifpca_weather_cumulative_variance.png`
- `lin_ifpca_weather_metrics.csv`
- `lin_ifpca_weather_support.csv`
- `lin_ifpca_weather_fits.rds`
