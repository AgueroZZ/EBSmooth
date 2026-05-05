# 2026-05-05 Matérn efficiency and equivalence update

## Summary

- Public Matérn `backend = "laplace"` and `backend = "laplace_fisher"` now require a successful TMB fit. They no longer silently fall back to the internal R reference path.
- `backend = "laplace_r"` remains accepted for internal validation, but it is hidden from the primary backend default lists and documented as internal validation only.
- `.compute_diag_A_Qinv_At()` now has a selected-inverse path using `INLA::inla.qinv()` for large sparse systems when the required pattern is covered. Small systems and uncovered patterns keep the previous sparse Cholesky solve path.
- Known-noise log-link `matern_n_starts` now strictly controls the number of
  Step A starts. `matern_n_starts = 1` means exactly one start; explicit
  multistart fits require `matern_n_starts > 1`.

## Equivalence evidence

Baseline fitted outputs were saved to `/tmp/ebsmoothr_matern_efficiency_baseline.rds` before edits and compared after edits on:

- known-noise softplus 2D
- known-noise log 1D
- learned-noise log 1D
- explicit `laplace_r` learned-noise log 1D
- known-noise log 2D

All compared quantities except the intentional known-noise log start-count
semantic fix had zero observed difference after the update:

- `posterior$mean`
- `posterior$var`
- `posterior$second_moment`
- `fitted_g$theta`
- `fitted_g$sigma`
- `fitted_beta`
- `log_likelihood`
- backend diagnostics including `laplace_implementation`, `tmb_mode_source`, and `stepA_n_starts`

The old known-noise log 2D baseline had `stepA_n_starts = 2` and
`stepA_best_start = 2`, even though the user-facing default was
`matern_n_starts = 1`. This was treated as a semantic bug rather than a
speed-only optimization target. After the correction:

- default known-noise log 2D: `stepA_n_starts = 1`, elapsed `22.033s`
- explicit `matern_n_starts = 1`: `stepA_n_starts = 1`, elapsed `23.514s`
- explicit `matern_n_starts = 2`: `stepA_n_starts = 2`, elapsed `44.433s`
- explicit `matern_n_starts = 5`: `stepA_n_starts = 5`, elapsed `118.126s`

Default and explicit one-start outputs were identical for posterior mean,
posterior variance, and log likelihood. Explicit two/five-start fits differed
from one-start by about `2.9e-10` in posterior means on the representative
fixture, while running the requested additional starts.

## Timing evidence

Before/after timing used the same machine and local source loading.

- Small known-noise softplus 2D fit: baseline median `16.383s`, after median `17.169s`. This path does not use the selected-inverse threshold and the small difference is optimizer noise.
- Learned-noise log 1D fit: baseline median `0.227s`, after median `0.280s`. The absolute difference is about `0.05s` and this path does not use the selected-inverse threshold.
- Small diagonal helper `n = 800`: baseline median `0.011s`, after median `0.012s`; threshold avoids the slower qinv path on small systems.
- Large 2D SPDE diagonal helper `n_spde = 3483`: old solve probe `1.213s`, selected-inverse helper after update median `0.602s`.
- Correcting known-noise log default starts removes the unintended second Step
  A optimization from default fits. On the representative 2D log fixture, the
  default elapsed time moved from the old two-start baseline `38.434s` to
  `22.033s`.

## Validation

- `pkgload::load_all("EBSmoothr", export_all = TRUE, quiet = TRUE)` succeeded.
- Custom targeted assertions passed for:
  - public `laplace` / `laplace_fisher` erroring instead of using R when TMB is unavailable,
  - explicit `laplace_r` still working for internal validation,
  - known-noise and learned-noise log-link fits respecting requested
    `matern_n_starts`,
  - diagonal selected-inverse variance matching dense reference,
  - pattern-covered sparse selected-inverse variance matching dense reference,
  - pattern-uncovered sparse designs falling back internally and matching dense reference.
- `testthat::test_file("EBSmoothr/tests/testthat/test-lgp.R")` passed with 49 assertions.
- A full `test-matern.R` run was not used as a final gate because INLA emitted extremely large external-program diagnostics on this machine; targeted checks were used for this change instead.
