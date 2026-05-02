# Quiet INLA Defaults for Local Runs

- Diagnosed repeated `/bin/kstat` noise as coming from `INLA::inla.getOption.default()`, which calls `parallel::detectCores(all.tests = TRUE, logical = FALSE)` on this machine and emits shell errors while still returning `NA`.
- Diagnosed `Fail to read A:B[:C] from [NA:1]` as the downstream consequence of the same default-thread path, where INLA receives `num.threads = "NA:1"` and falls back internally.
- Added internal helpers in `EBSmoothr/R/02_Matern.R` to locally override INLA's default-option function during package-controlled INLA calls, using a quiet thread default based on `parallel::detectCores(all.tests = FALSE, logical = FALSE)` with fallback to `1:1`.
- Wrapped Matern mesh/SPDE construction, exact SPDE precision evaluation, and all INLA PC-prior fit helpers in the quiet-default context so both package use and tests avoid the noisy autodetection path.
- Added a regression test asserting that the quiet INLA defaults produce an explicit non-`NA` thread specification.
