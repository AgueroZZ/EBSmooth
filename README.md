# EBSmooth

`EBSmooth` is the research and writing workspace for the
[`EBSmoothr`](https://github.com/AgueroZZ/EBSmoothr) R package and its
supporting validation materials. In addition to continuous smoothers, the
package includes an exact symmetric binary Markov normal-means solver for
transcript-state problems with intron/exon states encoded as `0`/`1`. Its
`flip_prob = 0.5` limit is the exact equal-probability iid baseline.

The package implements empirical-Bayes smoothers for Gaussian means with two
continuous prior families:

- L-GP smoothers with TMB and Fisher-Laplace log-link backends for
  one-dimensional problems.
- Matern smoothers with an exact Gaussian backend and an optional
  sparse-Laplace / INLA-SPDE workflow for one-dimensional and two-dimensional
  problems, including log-link fits for positive mean functions. Log-link
  Matern and L-GP fits use Fisher Laplace by default under `backend = "auto"`;
  explicit `backend = "laplace"` keeps observed-Hessian Laplace semantics.
- Symmetric binary Markov smoothing for ordered binary states, with exact
  forward-backward inference conditional on the transition probability.

The canonical package source is maintained in the standalone
[`AgueroZZ/EBSmoothr`](https://github.com/AgueroZZ/EBSmoothr) repository and is
included here as the [`EBSmoothr/`](EBSmoothr) git submodule. Longer internal
notes, validation studies, and collaborator-facing vignette drafts live under
[`internal/`](internal).

## Repository Layout

- `EBSmoothr/`: submodule pointing to the canonical installable R package.
- `internal/math/`: mathematical notes that explain the modeling targets.
- `internal/simulations/`: simulation studies and benchmark scripts.
- `internal/vignettes/`: collaborator-focused walkthroughs and rendered
  outputs.
- `plan/`: implementation plans for major updates.
- `log/`: change logs for major updates.

Clone this research workspace together with the package submodule using:

```sh
git clone --recurse-submodules https://github.com/AgueroZZ/EBSmooth.git
```

For an existing clone, initialize the package submodule with:

```sh
git submodule update --init --recursive
```

## Installation

The package depends on `INLA`, which is distributed through the R-INLA
repository rather than CRAN. Copy and paste the following block into an R
console:

```r
install.packages("pak")

pak::repo_add(
  INLA = "https://inla.r-inla-download.org/R/stable"
)

pak::pak("AgueroZZ/EBSmoothr")
```

For interactive development from this repository, use:

```r
devtools::load_all("EBSmoothr")
```

## Choosing the Main Entry Point

`EBSmoothr` has three collaborator-facing workflows:

- If the observation standard errors `s` are known, use the
  `ebnm`-compatible generator interfaces:
  `ebnm_LGP_generator()` for one-dimensional L-GP smoothing and
  `ebnm_Matern_generator()` for one-dimensional or two-dimensional Matern
  smoothing.
- If the observation standard errors are unknown, use `eb_smoother()` with
  `s = NULL` to fit an empirical-Bayes smoother that learns one common
  observation noise SD.
- For ordered binary states with known `s`, use `ebnm_binary_markov()` with a
  fixed transition probability or an empirical-Bayes numerical search.

`eb_smoother()` also supports known `s`, but the top-level documentation below
focuses on the learned-noise case. For `ebnm`, `flashier`, or `ebmf`
integration, the generator interfaces remain the primary documented path.

## Quick Start

### Known `s`: `ebnm` workflow

Use the generator interfaces when `s` is known and you want an
`ebnm`-compatible fitter:

```r
library(EBSmoothr)

loc <- seq(0, 1, length.out = 100)
s <- rep(0.1, length(loc))
x <- sin(2 * pi * loc) + rnorm(length(loc), sd = s)

matern_fit <- ebnm_Matern_generator(locations = loc)(x, s)

head(matern_fit$posterior)
matern_fit$fitted_g
matern_fit$fitted_beta
```

For one-dimensional local Gaussian-process smoothing with known `s`, use
`ebnm_LGP_generator(LGP_setup(...))` in the same way.

For positive Matérn mean functions with known `s`, use the log link:

```r
positive_fit <- ebnm_Matern_generator(
  locations = loc,
  link = "log"
)(exp(sin(2 * pi * loc)) + rnorm(length(loc), sd = s), s)

head(positive_fit$posterior)
```

Use `fix_params` to hold selected EB parameters fixed while estimating the
rest. For example, this holds the Matern marginal SD fixed at 1 while still
estimating range and intercept:

```r
sigma_fixed_fit <- eb_smoother(
  x,
  s = 0.1,
  family = "matern",
  locations = loc,
  g_init = Matern(sigma = 1),
  fix_params = "sigma"
)
```

### Known `s`: symmetric binary Markov model

For ordered transcript positions, `ebnm_binary_markov()` replaces the iid prior
with the symmetric transition matrix

```text
       next state
       0       1
0    1-q       q
1      q     1-q
```

where `q` is the probability of changing state between adjacent positions. The
initial intron/exon probabilities are both `0.5`, so every prior marginal also
remains `0.5`. The expected run length is `1 / q`, and `q = 0.5` reduces exactly
to the iid Bernoulli(0.5) baseline because both transition rows then equal
`(0.5, 0.5)`.

By default, the solver numerically searches for `q` using exact marginal-
likelihood evaluations:

```r
x <- c(0.05, 0.2, 0.55, 0.9, 1.1)
s <- rep(0.2, length(x))

markov_fit <- ebnm_binary_markov(x, s)

markov_fit$fitted_g
markov_fit$posterior[, c("prob_zero", "prob_one")]
markov_fit$posterior_transitions
markov_fit$viterbi_path
```

Supply `flip_prob` to hold the transition probability fixed:

```r
fixed_markov_fit <- ebnm_binary_markov(x, s, flip_prob = 0.05)
```

Use the same solver for the equal-probability iid reference:

```r
iid_fit <- ebnm_binary_markov(x, s, flip_prob = 0.5)
```

Conditional on a fixed `q`, inference and marginal-likelihood evaluation use an
exact log-space forward-backward algorithm. The default empirical-Bayes search
uses `stats::optimize()` and separately evaluates both interval endpoints. The
likelihood need not be concave, so this numerical search is not guaranteed to
find the global MLE. The current symmetric model assigns the same expected run
length to introns and exons; an asymmetric transition model is a planned
biological extension.

### Unknown `s`: empirical-Bayes smoothing

Use `eb_smoother()` when `s` is unknown and you want to learn one common noise
SD. `family = "constant"` provides the Gaussian baseline model for
comparison or for cases where no spatial structure is assumed:

```r
library(EBSmoothr)

loc <- seq(0, 1, length.out = 100)
x <- sin(2 * pi * loc) + rnorm(length(loc), sd = 0.1)

fit_learned <- eb_smoother(
  x = x,
  s = NULL,
  family = "matern",
  locations = loc,
  pc.penalty = list(range = 0.2, sigma = 0.3, noise = 0.1)
)

fit_baseline <- eb_smoother(
  x = x,
  s = NULL,
  family = "constant",
  beta_prec = 0
)

fit_learned
summary(fit_learned)
summary(fit_baseline)

head(fit_learned$posterior)
fit_learned$fitted_g
fit_learned$fitted_beta
```

For the Matern and constant families, the intercept is optimized by marginal
likelihood by default. Supply `beta_fixed` to hold it fixed, or
`beta_prec = 0 / >0` to use a flat or proper zero-mean Gaussian prior on the
intercept. Compact `print()` and `summary()` overview methods are available
for the `eb_smoother_fit` objects returned by `eb_smoother()`.

## Additional Documentation

- Public API documentation is generated from roxygen comments in
  [`EBSmoothr/R/`](EBSmoothr/R).
- Internal mathematical notes are indexed in
  [`internal/math/README.md`](internal/math/README.md).
- Internal collaborator walkthroughs are indexed in
  [`internal/vignettes/README.md`](internal/vignettes/README.md).
- Simulation and benchmarking materials are indexed in
  [`internal/simulations/README.md`](internal/simulations/README.md).
