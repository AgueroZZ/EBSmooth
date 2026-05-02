# Spatial scoring APIs implemented

## Summary

- Added public `spatial_scores()`, `spatial_select()`, and
  `spatial_scores_permutation()` wrappers with roxygen documentation and
  generated Rd pages.
- Renamed the public high-level Gaussian baseline to `family = "constant"` and
  exported `Constant()`. The old `family = "nonspatial"` path now errors.
- Added shared high-level support for `point_exponential`, `point_normal`, and
  `point_laplace`, including fixed-noise and profiled-noise fits.
- Added high-level Matern `fix_g` support so permutation scoring can reuse
  observed Matern range/sigma parameters when `refit = FALSE`.
- Updated the nonnegative log-Matern simulation notebook to use
  `EBSmoothr::spatial_select()` with `spatial = list(family = "matern",
  link = "log")` and `reference = list(family = "point_exponential")`.
- Updated current package docs and examples from `nonspatial` to `constant`.

## Validation

- `Rscript -e 'pkgload::load_all("EBSmoothr", quiet=TRUE, export_all=TRUE); testthat::test_file("EBSmoothr/tests/testthat/test-eb-smoother.R")'`
  - Passed: 104 assertions.
- `Rscript -e 'pkgload::load_all("EBSmoothr", quiet=TRUE, export_all=TRUE); testthat::test_file("EBSmoothr/tests/testthat/test-spatial-scores.R")'`
  - Passed: 44 assertions.
- `Rscript -e 'pkgload::load_all("EBSmoothr", quiet=TRUE, export_all=TRUE); testthat::test_file("EBSmoothr/tests/testthat/test-matern.R")'`
  - Passed: 150 assertions.
  - INLA emitted retry/segfault messages in child processes during existing
    tests, but testthat completed with zero failures.
- `Rscript -e 'r_files <- list.files("EBSmoothr/R", full.names = TRUE); r_files <- r_files[grepl("[.]R$", r_files)]; invisible(lapply(r_files, parse)); cat("R parse ok\n"); rd_files <- list.files("EBSmoothr/man", full.names = TRUE); rd_files <- rd_files[grepl("[.]Rd$", rd_files)]; invisible(lapply(rd_files, tools::checkRd)); cat("Rd check completed\n")'`
  - Passed.
