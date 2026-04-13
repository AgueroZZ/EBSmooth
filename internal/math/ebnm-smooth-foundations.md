---
title: "Smooth EBNM Foundations for EBSmoothr"
---

# Smooth EBNM Foundations for `EBSmoothr`

## Purpose

This note records the basic mathematics behind smooth empirical Bayes normal
means (smooth EBNM), and explains how the current `EBSmoothr` package implements
that framework.

The intended audience is future collaborators on this repository who need a
compact reference for:

- the core EBNM problem;
- the smooth-EBNM extension with Gaussian-process-style smoothing priors;
- the specific L-GP and Matern parameterizations used in `EBSmoothr`;
- the places where the current implementation matches the mathematics;
- the places where the current implementation still needs correction.

## Source Context

This note is based on:

- the classical EBNM summary in the Stephens Lab reading list;
- the `ebnm` package paper describing the standard two-step EBNM workflow;
- the smooth-EBNM formula shown in the screenshot supplied with this project on
  2026-04-11;
- the current package code in `EBSmoothr/R/01_LGP.R`,
  `EBSmoothr/R/02_Matern.R`, and `EBSmoothr/src/EBSmoothr.cpp`.

## 1. Classical EBNM

The classical empirical Bayes normal means problem starts from observations
`x = (x_1, ..., x_n)` with known standard errors `s = (s_1, ..., s_n)`:

```text
x_i | theta_i, s_i ~ N(theta_i, s_i^2),   i = 1, ..., n.
```

The empirical Bayes idea is to assume the unknown means are themselves drawn
from a prior distribution `g` in some family `G`:

```text
theta_i ~ g,   g in G.
```

The standard EBNM workflow has two steps:

```text
Step A: estimate g by maximizing the marginal likelihood

    g_hat = argmax_{g in G} p(x | s, g)
          = argmax_{g in G} product_i int p(x_i | theta_i, s_i) g(theta_i) dtheta_i.

Step B: compute posterior distributions and posterior summaries

    p(theta_i | x_i, s_i, g_hat).
```

In other words, EBNM is not just posterior inference under a fixed prior. It is
posterior inference after first estimating the prior from the full dataset.

## 2. Smooth EBNM

Smooth EBNM replaces the iid prior on `theta_i` with a structured prior over a
function defined on locations `s_i`.

The screenshot supplied with the project uses the following template:

```text
y_i | theta(s_i) ~ N(theta(s_i), sigma_i^2),
theta(.) ~ h(GP),   GP in G.
```

This is the right high-level idea. For implementation work, it is often clearer
to separate the latent Gaussian process from the observation-scale mean:

```text
x_i | eta(s_i), s_i ~ N(h(eta(s_i)), s_i^2),
eta(.) ~ GP_lambda.
```

Here:

- `eta(.)` is a latent smooth Gaussian field indexed by location;
- `h` is a link function mapping the latent field to the observation-scale mean;
- `lambda` denotes the smoothing hyperparameters of the prior family.

In this notation, the observation-scale mean is

```text
theta(s_i) = h(eta(s_i)).
```

This is the formulation that best matches the current package implementation.

### 2.1 Smooth EBNM as a two-step procedure

Smooth EBNM still follows the same EB logic as classical EBNM:

```text
Step A: estimate the smoothing hyperparameters lambda by maximizing
        the marginal likelihood

    lambda_hat = argmax_lambda p(x | s, lambda)
               = argmax_lambda int p(x | eta, s) p(eta | lambda) deta.

Step B: compute posterior summaries under lambda_hat

    p(eta | x, s, lambda_hat)
    and therefore
    p(theta | x, s, lambda_hat).
```

So the main difference from classical EBNM is not the EB structure. The
difference is the form of the prior. Instead of an iid prior `g` on each
`theta_i`, we use a joint smoothing prior over the full latent mean vector.

### 2.2 Identity and log links

The current package is built around Gaussian observations. At the package
level, the useful links are:

```text
Identity link: theta_i = eta_i.
Log link:      theta_i = exp(eta_i).
```

If Step B yields a Gaussian approximation

```text
eta_i | x, s, lambda_hat approx N(m_i, v_i),
```

then:

```text
Identity link:
    E(theta_i | x) = m_i,
    Var(theta_i | x) = v_i.

Log link:
    E(theta_i | x) approx exp(m_i + v_i / 2),
    Var(theta_i | x) approx exp(2 m_i + v_i) (exp(v_i) - 1).
```

These are exactly the moment transforms used in the L-GP code for the `log`
link.

For the current exact Matern implementation, only the identity link is
supported. This is because the closed-form marginal likelihood and posterior
updates are implemented for the Gaussian normal-means case
`x_i | eta_i ~ N(eta_i, s_i^2)`.

## 3. The L-GP Model in `EBSmoothr`

### 3.1 Mathematical model

The L-GP implementation is a one-dimensional smooth EBNM model on locations
`t_1, ..., t_n`. The latent field is represented in basis form:

```text
eta = X beta + B U.
```

Here:

- `X` is a global polynomial design matrix;
- `beta` are global trend coefficients;
- `B` is a local basis matrix;
- `U` are local coefficients.

The package assigns a Gaussian prior to `U` with precision matrix

```text
Q_U(theta) = exp(theta) P,
```

so that

```text
U | theta ~ N(0, Q_U(theta)^(-1)).
```

When `betaprec > 0`, the package also uses a proper Gaussian prior on `beta`:

```text
beta ~ N(0, (betaprec I)^(-1)).
```

When `betaprec = 0`, `beta` is treated as diffuse. When `betaprec < 0`, the
package uses an empirical-Bayes mode in which `beta` is optimized in Step A
rather than integrated out.

The observation model is:

```text
x_i | eta_i, s_i ~ N(eta_i, s_i^2)             for the identity link,
x_i | eta_i, s_i ~ N(exp(eta_i), s_i^2)        for the log link.
```

### 3.2 How this is implemented

The core TMB objective is in `EBSmoothr/src/EBSmoothr.cpp`. It evaluates the
joint log density

```text
log p(x | eta, s) + log p(U | theta) + log p(beta),
```

with:

- `eta = X beta + B U`;
- Gaussian likelihood for `x`;
- Gaussian prior on `U`;
- optional Gaussian prior on `beta`.

This is the right joint target for a smooth EBNM model with a latent Gaussian
basis representation.

The R-side workflow is:

1. `LGP_setup()` builds `X`, `B`, `P`, `logPdet`, and the metadata needed by
   TMB.
2. `ebnm_LGP_generator()` performs Step A and Step B using Laplace
   approximation.
3. Posterior means and variances are returned in an `ebnm`-style object.

### 3.3 Step A / Step B logic

For the main branches, the code matches the intended mathematics:

- if `betaprec >= 0`, Step A integrates out both `U` and `beta` and optimizes
  `theta`;
- if `betaprec < 0`, Step A integrates out `U` and optimizes `(theta, beta)`;
- Step B then fixes the Step A hyperparameter estimates and computes posterior
  summaries for the latent mean.

In the main non-`fix_g` branches, this is mathematically aligned with smooth
EBNM.

### 3.4 `fix_g` semantics in the current code

The current L-GP code now uses the following convention:

```text
If fix_g = TRUE:
    theta is fixed at g_init$scale,
    beta is fixed at beta_fixed,
    and beta_fixed defaults to zero if omitted.
```

This is computationally simple and explicit. It also avoids the earlier
ambiguity in which some branches silently re-estimated `beta` while others
implicitly fixed it at zero.

## 4. The Matern Model in `EBSmoothr`

### 4.1 Mathematical model

The Matern implementation is a smooth EBNM model in one or two spatial
dimensions. Let `w` denote the latent Matern field on a mesh, and let `A` be
the projector from mesh nodes to observation locations. Then the latent mean is

```text
eta = beta0 + A w.
```

The prior on `w` is a stationary Matern Gaussian field implemented through the
INLA/SPDE representation. The current observation model is Gaussian with the
identity link:

```text
x_i | eta_i, s_i ~ N(eta_i, s_i^2).
```

The smoothing hyperparameters optimized by the package are:

- `range`, parameterized internally as `theta = log(range)`;
- `sigma`, the latent-field marginal standard deviation;
- `beta0`, the global intercept.

### 4.2 How this is implemented

The current package now supports two Matern backends.

#### Exact backend

If `pc.penalty = NULL`, the current package code:

1. normalizes the locations and builds a mesh;
2. constructs a stationary SPDE template with `INLA::inla.spde2.matern()`;
3. for any fixed `(range, sigma, beta0)`, builds the latent precision matrix
   `Q(range, sigma)` with `INLA::inla.spde2.precision()`;
4. evaluates the exact Gaussian marginal likelihood

```text
ell = -0.5 * [ n log(2 pi) + log|D| - log|Q| + log|Q_post|
               + r' D^{-1} r - b' Q_post^{-1} b ],
```

   where

```text
D = diag(s_1^2, ..., s_n^2),
r = x - beta0,
b = A' D^{-1} r,
Q_post = Q + A' D^{-1} A;
```

5. maximizes this exact marginal likelihood over `(log(range), log(sigma), beta0)`
   when `fix_g = FALSE`;
6. computes the exact Gaussian posterior

```text
w | x, s, lambda_hat ~ N(Q_post^{-1} b, Q_post^{-1}),
```

   and then maps that posterior to the observation scale via `beta0 + A w`.

This exact branch is the clean empirical-Bayes implementation: the primary
`log_likelihood` is the exact Step A marginal likelihood.

#### PC-prior backend

If `pc.penalty` is supplied, the package switches to an INLA backend based on
`INLA::inla.spde2.pcmatern()`.

- With `fix_g = FALSE`, the fitter runs an INLA Step A optimization with free
  Matern hyperparameters and free intercept.
- The main `log_likelihood` is then taken to be the Step A penalized objective
  `resA$misc$log.posterior.mode`.
- The object also records Step A marginal-likelihood diagnostics such as
  `resA$mlik["log marginal-likelihood (integration)"]` and
  `resA$mlik["log marginal-likelihood (Gaussian)"]`.
- By default, the fitter does not automatically recompute the manual exact
  Gaussian objective at the Step A mode, because that extra check can be
  expensive at large `n`.
- Instead, the fitted object can be inspected afterwards with
  `matern_objective_breakdown(fit)`, which reports the exact fixed-beta,
  profiled-beta, and integrated-flat objectives, as well as the PC-prior term
  when present.
- With `fix_g = TRUE`, the fitter instead runs an INLA Step B computation with
  fixed hyperparameters and fixed intercept.
- In this fixed-parameter PC-prior branch, the main `log_likelihood` is a
  penalized pseudo-objective:

```text
Step B conditional marginal likelihood
+ log p(log(range), log(sigma) | PC prior).
```

This split reflects the empirical fact that Step A quantities are
prior-sensitive, while Step B `mlik` becomes prior-invariant once the
hyperparameters are held fixed.

### 4.3 Where the current implementation matches the intended math

The current Matern code now matches the intended smooth-EBNM mathematics in the
main ways:

- the observation model is the Gaussian normal-means model with known `s_i`;
- the prior is a stationary Matern Gaussian field;
- in exact mode, the primary `log_likelihood` is the exact empirical-Bayes
  marginal likelihood after optimizing the smoothing hyperparameters;
- in `pc.penalty` mode with `fix_g = FALSE`, the primary `log_likelihood` is
  the Step A penalized objective used by the INLA fit;
- in `pc.penalty` mode with `fix_g = TRUE`, the primary `log_likelihood` is an
  explicitly labeled penalized fixed-parameter surrogate;
- posterior means and variances come from the exact Gaussian posterior in exact
  mode and from INLA posterior summaries in PC-prior mode.

### 4.4 Current limitations

The Matern implementation is now much cleaner mathematically, but it still has
scope limitations by design:

- only the identity link is supported;
- `alpha` is fixed by the user and is not optimized;
- the meaning of `log_likelihood` is backend-dependent and must be interpreted
  together with the backend-specific diagnostics;
- `posterior_sampler` is implemented only for the exact backend;
- the compatibility field `inla_result` is `NULL` in exact mode and is
  populated only in the optional PC-prior mode.

## 5. A Unified View of the Package

Both `LGP` and `Matern` should be understood as specific prior families inside a
single smooth-EBNM template:

```text
x_i | eta_i, s_i ~ N(h(eta_i), s_i^2),
eta ~ smooth Gaussian prior with hyperparameters lambda.
```

Then:

- Step A estimates `lambda`;
- Step B computes posterior summaries given `lambda_hat`;
- the returned `fitted_g` should summarize the fitted smoothing prior;
- the returned `log_likelihood` should represent the primary Step A objective
  for the active backend, or an explicitly labeled penalized surrogate when the
  hyperparameters are fixed.

This is the cleanest way to think about the package.

## 6. Implementation Summary by File

### `EBSmoothr/src/EBSmoothr.cpp`

- Defines the joint likelihood and prior terms for the L-GP model.
- Encodes the identity-link and log-link observation models.
- Provides the Laplace target used by the TMB-based fitter.

### `EBSmoothr/R/01_LGP.R`

- Builds the L-GP basis and penalty objects.
- Runs Step A and Step B for the TMB-based L-GP smoother.
- Converts latent Gaussian moments to observation-scale moments for the log
  link.

### `EBSmoothr/R/02_Matern.R`

- Builds the mesh and projector matrix for one-dimensional or two-dimensional
  locations.
- Defines either an exact stationary Matern SPDE template or an INLA
  `pcmatern` template, depending on whether `pc.penalty` is supplied.
- Evaluates and optimizes the exact Gaussian marginal likelihood in exact mode.
- Runs INLA Step A or Step B in the optional PC-prior mode.
- Returns backend-specific likelihood diagnostics, posterior summaries, and
  mesh information.

## 7. Practical Takeaways for Future Development

If the package is to behave as a mathematically clean smooth-EBNM solver, the
following principles should be maintained:

1. `log_likelihood` should always be clearly documented as the primary
   objective for the active backend.
2. `fix_g` should have explicit semantics for both smoothing hyperparameters
   and any fixed global coefficients.
3. `g_init` should continue to function as a real initialization or fixed
   hyperparameter carrier.
4. The exposed `link` choices should match both the mathematical model and the
   implemented backend.
5. Backend-specific diagnostic fields should remain explicit whenever exact EB
   and INLA-PC modes expose different objective quantities.

## References

- Stephens Lab empirical Bayes reading list:
  https://stephenslab.github.io/reading_lists/empirical_bayes.html
- Willwerscheid, Carbonetto, and Stephens (2021), `ebnm` package paper:
  https://arxiv.org/abs/2110.00152
- Smooth-EBNM formula excerpt provided in the local screenshot on 2026-04-11:
  `/var/folders/32/d3z6m_356kj8891z4qvp50340000gn/T/TemporaryItems/NSIRD_screencaptureui_vgxX8e/Screenshot 2026-04-11 at 11.21.24 PM.png`
