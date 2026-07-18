# Standalone EBSmoothr repository

## Outcome

- Created the public canonical package repository at
  <https://github.com/AgueroZZ/EBSmoothr>.
- Preserved the package-only history with `git subtree split --prefix=EBSmoothr`.
- Published standalone `main` at commit
  `9bec472fc645779e8d97a36d48d76c2d52986710`.
- Published the annotated tag and release
  [`v0.2.6`](https://github.com/AgueroZZ/EBSmoothr/releases/tag/v0.2.6).
- Replaced the package directory in the research repository with a git
  submodule pointing to the standalone package commit.

## Public installation

```r
install.packages("pak")

pak::repo_add(
  INLA = "https://inla.r-inla-download.org/R/stable"
)

pak::pak("AgueroZZ/EBSmoothr")
```

The release-specific form is
`pak::pak("AgueroZZ/EBSmoothr@v0.2.6")`.

## Verification

- `devtools::test()` passed all 827 tests.
- `R CMD build` produced `EBSmoothr_0.2.6.tar.gz`.
- `R CMD check --as-cran --no-manual` completed with 0 errors, 0 warnings,
  and 5 non-blocking notes from existing implementation or local
  network/toolchain checks.
- A clean temporary-library installation from
  `AgueroZZ/EBSmoothr@v0.2.6` succeeded, exported
  `ebnm_binary_markov()`, and returned the expected Viterbi path in a smoke
  test.
- The locally built release tarball has SHA-256
  `286fede005f7cbc147f2cd682a6e9c9ec030456ab8b4574d5bb9556cc6ed1a70`.

## Safety

The pre-existing dirty linked worktree was moved before the repository split
and remains intact at
`.claude/worktrees/wizardly-thompson-263977`. Its seven modified files and
three untracked files were not staged, overwritten, or published as part of
this migration.
