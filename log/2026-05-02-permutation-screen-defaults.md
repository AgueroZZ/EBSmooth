# Permutation screening defaults updated

## Summary

- Changed `spatial_scores_permutation()` to default to `refit = FALSE`.
- Updated the one-sided permutation p-value to use only valid permutation
  scores: `mean(permutation_score >= observed_score)`.
- The observed score is no longer included through a `+1` pseudo-count in the
  numerator or denominator.
- Updated `spatial_select(method = "permutation")` so it propagates
  `matern_fits` and `fits` from the permutation result when available.
- Updated targeted tests to assert the default fixed-parameter behavior and the
  direct p-value calculation.

## Validation

- `Rscript -e 'pkgload::load_all("EBSmoothr", quiet=TRUE, export_all=TRUE); testthat::test_file("EBSmoothr/tests/testthat/test-spatial-scores.R")'`
  - Passed: 46 assertions.
- `Rscript -e 'roxygen2::roxygenize("EBSmoothr")'`
  - Regenerated `spatial_scores_permutation.Rd`.
- `Rscript -e 'r_files <- list.files("EBSmoothr/R", full.names = TRUE); r_files <- r_files[grepl("[.]R$", r_files)]; invisible(lapply(r_files, parse)); rd_files <- list.files("EBSmoothr/man", full.names = TRUE); rd_files <- rd_files[grepl("[.]Rd$", rd_files)]; invisible(lapply(rd_files, tools::checkRd)); cat("package parse/Rd ok\n")'`
  - Passed.
