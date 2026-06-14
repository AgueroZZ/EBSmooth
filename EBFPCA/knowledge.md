# EB-FPCA Knowledge Document

Last updated: 2026-06-11

## Purpose

This document tracks background reading and design implications for using the
current `EBSmoothr` package as the basis for empirical-Bayes functional
principal component analysis (EB-FPCA). It should be updated iteratively as new
papers are added, as implementation choices are made, and as simulations reveal
which simplifications are acceptable.

## Source Inventory

The Zotero screenshot in `related-paper/list-of-paper.png` lists 14 papers. I
copied 12 PDF attachments that were present in local Zotero storage into
`related-paper/`.

Copied PDFs:

- `Boland et al. - 2023 - Central Posterior Envelopes for Bayesian Functional Principal Component Analysis.pdf`
- `Gertheiss et al. - 2023 - Functional Data Analysis An Introduction and Recent Developments.pdf`
- `James et al. - 2000 - Principal component models for sparse functional data.pdf`
- `Jiang et al. - 2020 - BayesTime Bayesian Functional Principal Components for Sparse Longitudinal Data.pdf`
- `Jiang et al. - 2022 - Bayesian multivariate sparse functional principal components analysis with application to longitudin.pdf`
- `Li et al. - 2013 - Selecting the Number of Principal Components in Functional Data.pdf`
- `Lin et al. - 2016 - Interpretable Functional Principal Component Analysis.pdf`
- `Sang et al. - 2025 - Functional principal component analysis with informative observation times.pdf`
- `Sartini et al. - 2026 - Fast Bayesian Functional Principal Components Analysis.pdf`
- `Suarez and Ghosal - 2017 - Bayesian Estimation of Principal Components for Functional Data.pdf`
- `Yao et al. - 2005 - Functional Data Analysis for Sparse Longitudinal Data.pdf`
- `Ye - Functional principal component models for sparse and irregularly spaced data by Bayesian inference.pdf`

Not copied because Zotero did not currently store a PDF attachment:

- van der Linde (2008), `Variational Bayesian functional PCA`: Zotero stores a
  ScienceDirect HTML snapshot and full-text cache, but no PDF attachment.
- Goldsmith, Zipunnikov, and Schrack (2015), `Generalized Multilevel
  Function-on-Scalar Regression and Principal Component Analysis`: Zotero stores
  metadata and abstract, but no attachment.

## Common Modeling Template

Most papers build on the same sparse FPCA observation model:

```text
y_ij = mu(t_ij) + sum_{k = 1}^K xi_ik phi_k(t_ij) + eps_ij,
eps_ij ~ N(0, sigma^2 or s_ij^2),
xi_i = (xi_i1, ..., xi_iK) ~ N(0, Lambda),
```

where `mu(t)` is the population mean curve, `phi_k(t)` are orthonormal
eigenfunctions, `xi_ik` are subject scores, and observations may be sparse,
irregular, noisy, or partially observed. Generalized versions replace the
Gaussian observation model with a link-function likelihood.

The empirical-Bayes angle for `EBSmoothr` is to estimate the smoothing
hyperparameters, score variances, noise level, and possibly `K` from the
marginal likelihood or a predictive criterion, then report posterior summaries
for curves, scores, mean, and eigenfunctions conditional on the fitted prior.

## Software Leads

The papers point to several usable software paths, but only a subset are direct
R packages rather than supplementary code or Matlab code.

### Directly Usable R Packages or R-Centric Code

- `fdapace`: a CRAN/R-universe R implementation of PACE-style FPCA for sparse
  or dense functional data. This is the most immediately useful frequentist
  baseline for Yao, Mueller, and Wang (2005)-style sparse longitudinal FPCA and
  for score prediction by conditional expectation.
- `fda::pcaPACE`: a smaller PACE-related function in the `fda` package that
  estimates functional principal components from a covariance estimate. It is
  useful as a reference, but `fdapace` is the more complete sparse-FPCA software
  path.
- `BayesTime`: an R package from Jiang et al. (2020), available at
  `github.com/biocore/bayestime`, implementing Bayesian SFPCA with Stan,
  PSIS-LOO model selection, Pareto-shape diagnostics, and posterior predictive
  checks.
- Jiang et al. (2022) implement multivariate SFPCA using R and Stan. This is
  directly relevant for later multivariate extensions, but the paper reads more
  like method-specific R/Stan code than a mature general package.
- Goldsmith, Zipunnikov, and Schrack (2015) use Stan with an R interface and
  public code for generalized multilevel function-on-scalar regression plus
  FPCA. This is useful for future generalized or multilevel EB-FPCA, but not a
  drop-in package for our first prototype.

### Not Direct R Packages

- Sartini et al. (2026) provide Stan code and simulation routines for FAST, but
  not a named R package in the paper. The code should be usable from R through
  `rstan` or `cmdstanr` if the supplementary material is obtained.
- Ye's Bayesian FPCA implementation is reported as Matlab code.
- van der Linde (2008), Suarez and Ghosal (2017), Lin et al. (2016), Boland et
  al. (2023), Sang et al. (2025), and Li, Wang, and Carroll (2013) are more
  important here as methodological references than as direct R software sources.
- Gertheiss et al. (2023) is a software-aware FDA review, not a new FPCA
  package paper.

## Open Data and Simulation Resources

This ranking separates open real datasets from open simulation code. For the
first EB-FPCA benchmark, the best targets are sources that either provide
paper-specific analysis code or expose a simple real functional dataset that we
can sparsify under controlled missingness.

### Best Immediate Benchmark Sources

- Jiang et al. (2020), BayesTime: the paper gives the most direct reproducible
  benchmark path. The `biocore/bayestime` R package includes a `data`
  directory and examples using `data("ECAM")`, and the separate
  `knightlab-analyses/BayesTime-analyses` repository contains the manuscript
  `applications` and `simulations` folders. This should be the first sparse
  longitudinal benchmark against `fdapace` and an EB-FPCA prototype.
- Jiang et al. (2022), multivariate SFPCA: the
  `knightlab-analyses/mfpca-analyses` repository contains `data`,
  `application_1_skin`, `application_2_T2D`, `simulation_studies`, and
  `simulations`. This is a strong second benchmark after the univariate
  workflow is stable, especially for later multivariate EB-FPCA.
- Canadian weather data, used across van der Linde (2008), Suarez and Ghosal
  (2017), and Lin et al. (2016): `fda::CanadianWeather` and `fda::daily`
  provide daily temperature and precipitation at 35 Canadian locations averaged
  over 1960-1994. This is dense and small, so it is ideal for a fast
  smoke-test: estimate a dense reference FPCA, sparsify each curve, add known
  noise, and compare reconstruction and eigenfunction recovery.
- Lin et al. (2016), interpretable FPCA: the EEG example can be connected to
  the public UCI EEG Database. It is larger and preprocessing-heavy, but it can
  be used as a dense functional dataset, then sparsified or downsampled for
  robustness tests.
- Boland et al. (2023), central posterior envelopes: the paper and supplement
  point to public R code for simulation data generation and a simulated-data
  tutorial. This is more useful for checking posterior uncertainty summaries
  and envelope-style displays than for a first real-data benchmark.
- Gertheiss et al. (2023), FDA tutorial: the accompanying GitHub repository has
  `Rcode` and `data` folders for reproducing tutorial examples, including a
  publicly available running dataset. This is useful for general FDA sanity
  checks, but it is less targeted to sparse Bayesian FPCA.

### Usable With Caveats

- Sartini et al. (2026), FAST: all relevant Stan code and simulation routines
  are provided as supplementary material, but the motivating DASH4D CGM data is
  not yet available because trial results are forthcoming. Use this for method
  simulation design, not for open real-data comparison.
- Goldsmith, Zipunnikov, and Schrack (2015): Stan/R code is publicly available
  in supplement or from the first author's site, but the motivating
  accelerometry dataset is not clearly available as a ready public dataset.
  This is better saved for generalized or multilevel simulation work.
- Sang et al. (2025): the GitHub repository `spj1125/FPCA` provides R code for
  simulations with informative observation times, but the CD4/viral-load
  application data is not clearly open. This is a good later stress test for
  informative sampling weights.
- Yao, Mueller, and Wang (2005): the PACE method is implemented in modern R
  software, and the paper's yeast example can be conceptually reproduced by
  sparsifying complete time-course data, but the exact paper data are not the
  cleanest immediate open benchmark in this folder.

### Recommended First Simulation Plan

1. Start with `fda::CanadianWeather` because it is small, deterministic, and
   already available through CRAN data.
2. Add the BayesTime ECAM workflow as the first sparse longitudinal
   paper-reproduction benchmark.
3. Use UCI EEG only after the API and metrics are stable, because preprocessing
   choices can dominate the method comparison.
4. Keep FAST, Boland CPE, Goldsmith multilevel, and Sang informative-time
   setups as later stress tests rather than MVP benchmarks.

### Local Supplement Data Status

The current local supplement file is
`EBFPCA/data/biom12457-sup-0002-suppdatacode.zip`. It is the Lin et al. iFPCA
supplement and contains MATLAB code for EEG, weather, two simulations, and an
included copy of the MATLAB FDA toolbox.

The first R-side exploration script is
`EBFPCA/scripts/explore_fpca_weather.R`. It reads the Canadian weather
temperature data directly from the zip member
`iFPCAcodes/fdaM/examples/weather/temperature.csv`, then compares:

- raw grid PCA via base R `prcomp`;
- smoothed B-spline FPCA via `fda::pca.fd`;
- smoothed Fourier FPCA via `fda::pca.fd`;
- optional `fdapace::FPCA`, if `fdapace` is installed;
- optional `refund::fpca.sc`, if `refund` is installed.

The first run wrote outputs to
`EBFPCA/results/fpca_weather_exploration/`. On the current machine, `fda` was
installed and ran successfully. After installing the legitimate CRAN packages
`fdapace` and `refund`, `fdapace::FPCA` and `refund::fpca.sc` also ran
successfully. The fitted methods agreed closely on the first three components:
PC1 explained about 88-90 percent of variation, PC2 about 8.5 percent, and PC3
about 2 percent. PACE was nearly identical to raw grid PCA on this dense
regular dataset; this is expected because the weather data are fully observed
on a common daily grid.

The interpretable-FPCA replication script is
`EBFPCA/scripts/replicate_lin_ifpca_weather.R`. It translates the core MATLAB
supplement algorithm for the weather analysis into R:

- B-spline representation of centered weather curves;
- roughness-penalized generalized eigenproblem for the regularized FPCA
  comparator;
- greedy L0 basis support selection for iFPCA;
- 10-fold CV kappa selection using the supplement's one-standard-error-style
  rule;
- PACE overlay through `fdapace::FPCA`.

The main replication outputs are in
`EBFPCA/results/lin_ifpca_weather_replication/`. The current recommended run
uses `nbasis = 181`, `gamma = 3000`, and 10 CV folds. It is an approximate R
replication, not a bit-for-bit MATLAB rerun, because modern R `fda` and the
MATLAB FDA toolbox differ in high-dimensional basis setup and numerical
solvers. Still, it reproduces the qualitative paper message: regularized FPCA
and PACE yield smooth global eigenfunctions, while iFPCA truncates components
to localized active intervals. In the current run, the selected active basis
counts are 103, 60, and 13 for PC1-PC3; the detected active support intervals
are approximately day 1.5-128.5 and 271.5-363.5 for PC1, day 123.5-181.5 and
195.5-275.5 for PC2, and day 95.5-126.5 for PC3.

The first EBM/F data-loading scaffold is `EBFPCA/scripts/fit_ebmf.R`. It reads
the same supplement weather temperature data and creates:

- `Y`: a 35 locations x 365 days matrix for matrix factorization;
- `Y_centered`: `Y` after subtracting the daily mean across locations;
- `weather_temperature_day_by_location`: the original 365 days x 35 locations
  matrix;
- `weather_long`: a long observation table with location, day, and
  temperature;
- `observed_mask`: a logical matrix for missingness-aware matrix methods.

Running the script writes reusable matrix inputs to
`EBFPCA/results/ebmf_weather_matrix/`.

## Paper-Level Takeaways

### Foundational Sparse FPCA

James, Hastie, and Sugar (2000) formulate sparse FPCA as a reduced-rank mixed
effects model. Smooth basis functions represent the mean and principal
component curves. The key lesson is that directly estimating a full spline
coefficient covariance matrix can overfit badly when curves are sparse; fitting
a rank-constrained model for the leading components is more stable.

Yao, Mueller, and Wang (2005) introduce PACE for sparse longitudinal data. They
estimate the mean and covariance surface nonparametrically, eigendecompose the
smoothed covariance, and estimate subject scores by conditional expectation
(BLUP under Gaussian assumptions). This is the most natural baseline for an
`EBSmoothr` prototype because it separates mean/covariance smoothing from score
prediction.

Li, Wang, and Carroll (2013) study how to choose the number of functional
principal components. Their marginal BIC is designed for sparse and dense
functional data; their conditional AIC is useful when the effective number of
components may grow. For EB-FPCA, this argues for treating `K` selection as a
first-class part of the API, not an afterthought.

Sang, Kong, and Yang (2025) show that observation times can be informative. If
visit times depend on the latent outcome trajectory, ordinary FPCA estimates can
be biased. They model the visit process through an intensity function and use
inverse-intensity weighting in penalized spline estimation. This is probably not
MVP scope, but the data object should allow optional observation weights so this
extension is not blocked later.

### Bayesian FPCA and EB-FPCA

van der Linde (2008) is conceptually close to what we want: a generative
Bayesian FPCA model for noisy and sparse curves, smooth eigenfunctions through a
Demmler-Reinsch-type basis, variational posterior approximation, model choice
over the number and smoothness of eigenfunctions, and posterior uncertainty for
reconstructed curves. Important warning: smoothness and orthogonality constraints
interact awkwardly, so rotation/sign/order handling must be explicit.

Ye develops Bayesian FPCA for sparse and irregular continuous or binary
responses. The mean and eigenfunctions are represented with penalized splines,
eigenfunctions live on a Stiefel manifold, and MCMC includes a
Langevin-Bingham step plus RJ-MCMC for `K`. This is too heavy for the first
`EBSmoothr` implementation, but it flags the right generalized-response path:
identity-link Gaussian first, then binary/count links later.

Suarez and Ghosal (2017) model the covariance function through an approximate
spectral decomposition and put priors on finite-dimensional covariance
structures. They emphasize model selection over both the number of principal
components and basis dimension, and warn that inverse-Wishart hyperparameters
must be chosen carefully to avoid unrealistic covariance smoothness. For
`EBSmoothr`, this supports a covariance-first EB prototype as a legitimate
baseline, but it is less aligned with the current latent-field smoother code
than a reduced-rank latent factor model.

Jiang et al. (2020) implement Bayesian SFPCA in Stan as BayesTime. The most
useful parts for our workflow are not the sampler itself but the diagnostics:
PSIS-LOO for selecting `K` and basis complexity, Pareto shape diagnostics for
outlying influential curves, and posterior predictive checks for model fit.

Jiang et al. (2022) extend SFPCA to multivariate longitudinal trajectories.
Scores are independent within an outcome for identifiability but can be
correlated across outcomes; a Cholesky parameterization keeps the score
covariance positive definite. This is a later extension for multi-omics or
multi-feature longitudinal settings, not a first implementation target.

Goldsmith, Zipunnikov, and Schrack (2015) combine generalized
function-on-scalar regression with multilevel FPCA in a Bayesian Stan framework.
The practical lesson is that FPCA can be integrated with covariate-adjusted
means and nested random effects. This is useful for future designs with
subject-level and replicate-level variation, but the first EB-FPCA should avoid
multilevel/covariate complexity.

Sartini et al. (2026) propose FAST, a fully Bayesian FPCA that projects
eigenfunctions onto an orthonormal spline basis, samples the orthonormal
coefficient matrix efficiently using polar-decomposition parameter expansion,
and orders eigenvalues during sampling. The key implementation lesson is that an
orthonormal basis can make the constraints much easier. If we implement a joint
model, an orthonormalized spline/grid basis should be considered before trying
to constrain arbitrary smoother outputs.

Boland et al. (2023) focus on uncertainty summaries for Bayesian FPCA. Their
central posterior envelopes use functional depth to summarize posterior samples
of mean functions and eigenfunctions, avoiding only pointwise symmetric
intervals. This is an output-layer idea: once `EBSmoothr` can sample or
approximate posterior functions, we should summarize uncertainty at both
pointwise and functional-envelope levels.

### Interpretability and General FDA Context

Lin, Wang, and Cao (2016) introduce interpretable FPCA with eigenfunctions that
are exactly zero outside important intervals. This can be viewed as support
selection for eigenfunctions. It is not necessary for MVP, but it suggests a
future "sparse eigenfunction" option, possibly through local shrinkage or
support penalties.

Gertheiss et al. (2023) provide the broad FDA context: smoothing, registration,
FPCA, regression, inference, clustering, classification, and software. For this
project, the important reminder is that FPCA is both an exploratory and a
model-building tool; outputs should include reconstructed curves, score plots,
variance explained, and diagnostic plots, not only fitted hyperparameters.

## Implications for `EBSmoothr`

The current package already solves a related problem:

```text
x_i | theta_i, s_i ~ N(theta_i, s_i^2)
theta(.) has a smooth Gaussian prior indexed by 1D or 2D locations
smooth hyperparameters are estimated by marginal likelihood
posterior summaries are returned under the fitted prior
```

EB-FPCA needs a second latent layer:

```text
y_ij | f_i(t_ij), s_ij ~ N(f_i(t_ij), s_ij^2)
f_i(t) = mu(t) + sum_k xi_ik phi_k(t)
mu and phi_k are smooth functions
xi_ik are latent subject scores with variances lambda_k
```

So the reusable pieces are:

- marginal-likelihood EB philosophy;
- Gaussian latent-field algebra;
- L-GP and Matern smoothness priors;
- exact Gaussian posterior updates where the model can be reduced to linear
  Gaussian form;
- learned-noise and known-noise interfaces;
- posterior summary conventions.

The missing pieces are:

- grouped observations by subject/curve;
- basis or grid representation for `mu` and `phi_k`;
- score integration or optimization;
- orthogonality, sign, and ordering constraints for eigenfunctions;
- selection of `K`;
- curve-level diagnostics and visualization.

## Recommended Development Path

### MVP 0: Covariance-First EB-FPCA Baseline

This is the fastest route to something useful:

1. Accept sparse observations as `(subject, time, y, s)` with optional weights.
2. Estimate `mu(t)` using existing one-dimensional `eb_smoother()` or
   `ebnm_*_generator()` on pooled observations, accounting for repeated times
   when possible.
3. Estimate a smoothed covariance surface from residual cross-products.
4. Eigendecompose the covariance on a common grid.
5. Estimate subject scores by conditional Gaussian prediction.
6. Return reconstructed curves, scores, eigenfunctions, eigenvalues, variance
   explained, and model-selection diagnostics over `K`.

This is PACE-like and not yet a fully joint EB factor model, but it will give a
baseline for simulation and API design.

### MVP 1: Joint Gaussian EB-FPCA

For a more principled `EBSmoothr` model, use:

```text
y_i = A_i mu + A_i Phi xi_i + eps_i
xi_i ~ N(0, Lambda)
mu ~ smooth prior
columns(Phi) ~ smooth priors plus orthonormality constraints
```

Start with fixed `K`, known `s`, one-dimensional time, and a common evaluation
grid. Use an orthonormal spline basis or post-update orthonormalization to make
constraints tractable. Estimate smoothing hyperparameters and eigenvalues by
marginal likelihood. Add unknown common noise after the known-noise path is
stable.

### Later Extensions

- Generalized responses: binary/count/log-link models using Fisher-PQL or
  Laplace ideas from the existing Matern non-identity work.
- Multilevel FPCA: separate subject-level and visit-level components.
- Multivariate FPCA: cross-outcome score covariance with a Cholesky
  parameterization.
- Informative observation times: optional inverse-intensity weights.
- Interpretable FPCA: support penalties or shrinkage for localized
  eigenfunctions.
- Bayesian uncertainty displays: posterior samples and central posterior
  envelopes.

## Design Decisions to Revisit

- Should the first public API be a separate `eb_fpca()` function, or an
  `eb_smoother()` extension? Current evidence favors a separate `eb_fpca()`
  because grouped curves, scores, and `K` selection are different from EBNM.
- Should the smoother prior for eigenfunctions be L-GP, Matern, or spline-first?
  For MVP, spline-first may simplify orthonormality; Matern/L-GP can be used for
  smoothing mean and covariance baselines.
- Should `K` be selected by marginal likelihood/BIC, PSIS-LOO, or held-out
  curves? For an EB package, start with marginal likelihood/BIC and add
  predictive checks later.
- Should scores be integrated out exactly or optimized as random effects? For
  Gaussian known-noise MVP, exact integration is preferred if computationally
  feasible.

## Iterative Update Protocol

When this document changes, update:

1. `Last updated`.
2. `Source Inventory` if papers or PDFs are added.
3. `Paper-Level Takeaways` for reading notes.
4. `Implications for EBSmoothr` and `Recommended Development Path` when
   implementation evidence changes.
5. `Update History` below with a one-line note.

## Update History

- 2026-06-11: Added `fit_ebmf.R` weather-data loading and matrix-format
  scaffold for EBM/F experiments.
- 2026-06-10: Added local supplement-data status and the first R weather-data
  FPCA exploration script/results.
- 2026-06-10: Installed CRAN `fdapace` and `refund`, reran PACE/refund
  baselines, and added an R translation of the Lin et al. weather iFPCA
  replication workflow.
- 2026-06-08: Added open-data and simulation-resource ranking for benchmark
  planning, separating immediate datasets from code-only or caveated resources.
- 2026-06-08: Added software inventory distinguishing directly usable R paths
  from supplementary Stan/Matlab code and methodological references.
- 2026-06-08: Created first knowledge document from Zotero screenshot, local
  PDF copies, Zotero metadata/cache, and a quick inspection of the current
  `EBSmoothr` package structure.
