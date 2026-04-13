## Objective

Create several internal vignettes that illustrate how to use the `EBSmoothr` package without wiring them into the formal package vignette build.

## Plan

1. Review the exported package interfaces and existing sanity-check scripts to identify stable usage patterns.
2. Create an internal vignette directory under `EBSmoothr/inst/internal/vignettes`.
3. Add a short index document that explains the purpose of the internal vignettes and what each file covers.
4. Add a package overview vignette that introduces the empirical-Bayes smoothing workflow and the two prior families.
5. Add an L-GP vignette with identity-link and log-link examples using `LGP_setup()` and `ebnm_LGP_generator()`.
6. Add a Matern vignette with one-dimensional and two-dimensional examples using `ebnm_Matern_generator()`.
7. Record the major update in the project `log` folder so other agents can see what was added.
