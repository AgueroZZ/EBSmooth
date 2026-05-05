# 2026-05-05 Matérn efficiency and equivalence plan

## Goal

Speed up Matérn non-identity link computations where possible, tighten backend
semantics so public Laplace backends never silently route to the internal R
reference implementation, and make `matern_n_starts` mean exactly the requested
number of Step A starts.

## Implementation scope

- Enforce TMB-or-error semantics for public `backend = "laplace"` and `backend = "laplace_fisher"`.
- Keep `backend = "laplace_r"` accepted for internal validation only.
- Add a selected-inverse posterior variance path using `INLA::inla.qinv()` only when:
  - the system is large enough to avoid small-case slowdown,
  - the design is diagonal, or
  - a sparse pattern check proves that all required entries are present.
- Preserve the old sparse Cholesky solve path for small systems and uncovered patterns.
- Make known-noise and learned-noise log-link Step A multistart construction
  strictly respect `matern_n_starts`.
- Defer compiler flags, posterior overflow tolerance, `MakeADFun` multistart caching, and Step B trigger tuning.

## Acceptance criteria

- Baseline fitted outputs and diagnostics match after the update for stable
  fixtures except for the intentional known-noise log-link start-count semantic
  fix.
- Default `matern_n_starts = 1` and explicit `matern_n_starts = 1` both report
  `stepA_n_starts == 1`; explicit `matern_n_starts = 5` reports
  `stepA_n_starts == 5`.
- No posterior `NA` or `Inf` values are introduced.
- Public `laplace` / `laplace_fisher` error when TMB is unavailable or invalid.
- Explicit `laplace_r` remains available for internal validation.
- Large sparse diagonal posterior variance computation is faster than the old solve path.
- Small or uncovered variance cases are not materially slowed by selected-inverse overhead.
