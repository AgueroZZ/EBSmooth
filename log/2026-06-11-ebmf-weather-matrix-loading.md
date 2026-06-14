# EBM/F Weather Matrix Loading

## Summary

Added data loading and matrix formatting to `EBFPCA/scripts/fit_ebmf.R` for
the Canadian weather temperature data from the Lin et al. iFPCA supplement.

## Data Source

- zip: `EBFPCA/data/biom12457-sup-0002-suppdatacode.zip`
- member: `iFPCAcodes/fdaM/examples/weather/temperature.csv`

## Matrix Objects

The script creates these top-level R objects:

- `Y`: 35 locations x 365 days, suitable for matrix factorization.
- `Y_centered`: same dimensions as `Y`, after subtracting daily means across
  locations.
- `day_grid`: numeric grid from 0.5 to 364.5.
- `location_names`: 35 weather location labels from the supplement.
- `weather_temperature_day_by_location`: original 365 days x 35 locations
  matrix.
- `weather_long`: long table with location, day index, day, and temperature.
- `weather_inputs`: list containing all matrix inputs and metadata.

## Outputs

Running:

```sh
Rscript --vanilla EBFPCA/scripts/fit_ebmf.R
```

writes:

- `EBFPCA/results/ebmf_weather_matrix/weather_matrix_inputs.rds`
- `EBFPCA/results/ebmf_weather_matrix/weather_matrix_location_by_day.csv`
- `EBFPCA/results/ebmf_weather_matrix/weather_matrix_centered_by_day.csv`
- `EBFPCA/results/ebmf_weather_matrix/weather_observations_long.csv`

## Validation

Confirmed:

- `Y` is `35 x 365`.
- `Y_centered` is `35 x 365`.
- `weather_long` has `12775` rows.
- Column means of `Y_centered` are zero up to numerical precision.
