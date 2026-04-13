# EBSmooth

`EBSmooth` is the research repository for the `EBSmoothr` R package and its
supporting validation materials.

The package implements empirical-Bayes smoothers for Gaussian means with two
prior families:

- L-GP smoothers with a TMB backend for one-dimensional problems.
- Matern smoothers with an exact Gaussian backend and an optional
  INLA/SPDE PC-prior workflow for one-dimensional and two-dimensional
  problems.

The package source lives in [`EBSmoothr/`](EBSmoothr), while longer internal
notes, validation studies, and collaborator-facing vignette drafts live under
[`internal/`](internal).

## Repository Layout

- `EBSmoothr/`: installable R package source.
- `internal/math/`: mathematical notes that explain the modeling targets.
- `internal/simulations/`: simulation studies and benchmark scripts.
- `internal/vignettes/`: collaborator-focused walkthroughs and rendered
  outputs.
- `plan/`: implementation plans for major updates.
- `log/`: change logs for major updates.

## Installation

The package depends on `INLA`, which is distributed through the R-INLA
repository rather than CRAN. Install dependencies first, then install the
package from source.

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  )
)

install.packages(c("TMB", "Matrix", "numDeriv", "LaplacesDemon", "ebnm"))
install.packages("remotes")

remotes::install_local("EBSmoothr")
```

For interactive development from this repository, use:

```r
devtools::load_all("EBSmoothr")
```

## Quick Start

### L-GP smoothing

```r
library(EBSmoothr)

t <- seq(0, 1, length.out = 100)
s <- rep(0.1, length(t))
x <- sin(2 * pi * t) + rnorm(length(t), sd = s)

lgp_setup <- LGP_setup(t)
lgp_fit <- ebnm_LGP_generator(lgp_setup)(x, s)

head(lgp_fit$posterior)
```

### Matern smoothing

```r
library(EBSmoothr)

loc <- seq(0, 1, length.out = 100)
s <- rep(0.1, length(loc))
x <- 0.25 + sin(2 * pi * loc) + rnorm(length(loc), sd = s)

matern_fit <- ebnm_Matern_generator(locations = loc)(x, s)

head(matern_fit$posterior)
```

## Additional Documentation

- Public API documentation is generated from roxygen comments in
  [`EBSmoothr/R/`](EBSmoothr/R).
- Internal mathematical notes are indexed in
  [`internal/math/README.md`](internal/math/README.md).
- Internal collaborator walkthroughs are indexed in
  [`internal/vignettes/README.md`](internal/vignettes/README.md).
- Simulation and benchmarking materials are indexed in
  [`internal/simulations/README.md`](internal/simulations/README.md).
