# Known- vs Unknown-Noise Documentation Refresh

## What Changed
- Repositioned the top-level package documentation around two main workflows:
  - known standard errors `s` via `ebnm_LGP_generator()` and
    `ebnm_Matern_generator()`
  - unknown standard errors via `eb_smoother(s = NULL)`
- Simplified the `README.md` quick start so it emphasizes the two collaborator-
  facing entry points rather than showcasing every supported combination.
- Updated the package description and key function help pages to keep the
  workflow split consistent without changing any implementation behavior.
- Revised the overview vignette so it no longer implies all workflows require
  known standard errors, and added short scope notes to the L-GP and Matern
  workflow vignettes.

## Validation Plan
- Regenerate roxygen documentation.
- Re-render internal vignette HTML outputs.
- Run the targeted `test-eb-smoother.R` test file.
- Smoke-run the documented known-`s` generator example and the unknown-`s`
  `eb_smoother()` example.
