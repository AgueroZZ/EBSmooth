# Add Log-Link Learned-Noise Support to `eb_smoother()`

## Implemented

- Added nonspatial log-link support for known and learned common-noise fits.
- Added Matern log-link learned-noise support through R Laplace and TMB Laplace.
- Extended the Matern TMB objective with `learn_noise`, `log_noise`, and an
  optional PC prior on the learned noise SD.
- Kept `beta_prec = 0` as flat on the linear-predictor beta scale for all links.
- Updated Matern auto routing so log-link learned-noise fits default to the
  package Laplace/TMB path, with explicit `inla` and `inla_pc` still available.
- Added deterministic multi-starts for TMB learned-noise empirical-Bayes fits to
  avoid a local constant-signal/large-noise mode.
- Added package Laplace comparable objective evaluation for INLA log-link
  learned-noise fits using the INLA spatial mode as the warm start.
- Updated tests and documentation for log-link learned-noise Matern and
  nonspatial workflows.

## Validation Run

- `R CMD INSTALL EBSmoothr` completed successfully after the C++ and R changes.
- Smoke checks showed nonspatial log known/learned-noise fits return finite
  positive response summaries.
- Smoke checks showed Matern log learned-noise `laplace_tmb` and `laplace_r`
  agree for empirical-Bayes beta after the TMB multi-start fix.
- Smoke checks showed Matern log learned-noise TMB and R Laplace agree when a
  PC prior includes `noise`.
- Explicit INLA learned-noise log-link fixed-beta smoke fits returned finite
  package Laplace comparable log-likelihoods. INLA can still print external
  binary retry messages for failed candidate noise/profile evaluations, but the
  wrapper discards failed candidates and returns the best finite profile fit.

## Caveat

- Strict integration of a flat prior on beta under a log link is improper in
  some simple constant-mean cases because the likelihood can approach a
  non-zero constant as `beta -> -Inf`. The implementation follows the requested
  backend-independent flat-on-beta semantics and uses local Laplace behavior
  where needed rather than changing to a flat prior on `exp(beta)`.

## Follow-up: Backend Comparability Fixes

- Updated identity-link INLA beta profiling to use the package exact objective
  rather than INLA Step A diagnostics.
- Updated identity-link INLA primary `log_likelihood` for known-noise and
  learned-noise Matern fits to be the package exact objective evaluated at the
  INLA mode. The INLA Step A quantities remain available as diagnostics.
- Marked explicit `backend = "inla"` for `family = "matern", link = "log",
  s = NULL` without a PC prior as unsupported, because INLA can fail at low
  learned-noise candidates and otherwise return a non-comparable local profile
  point. The default `laplace`/`laplace_tmb` backend and `inla_pc` with a PC
  prior remain available.
- Added regression tests for identity-link INLA exact-objective semantics and
  for the unsupported no-PC log-link learned-noise INLA case.

## Superseded Follow-up

- The unconditional no-PC log-link learned-noise INLA rejection above was
  replaced on 2026-04-30 by validate-or-error behavior. Explicit INLA now runs
  and returns only when it validates against the package Laplace reference; see
  `log/2026-04-30-release-readiness-backend-policy.md`.
