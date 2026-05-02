# Refresh Documentation Positioning for Known- vs Unknown-Noise Workflows

## Summary
- Reframe the collaborator-facing documentation around the package's two main
  entry points:
  - known standard errors `s`: `ebnm_LGP_generator()` /
    `ebnm_Matern_generator()`
  - unknown standard errors: `eb_smoother(s = NULL)`
- Keep the implementation and exported APIs unchanged.
- Update both package reference docs and collaborator-facing vignettes so the
  top-level narrative is consistent.

## Planned Changes
- Update `README.md` to add a clear entry-point split and simplify the quick
  start to:
  - one known-`s` generator example
  - one unknown-`s` `eb_smoother()` example
- Update `EBSmoothr/DESCRIPTION` so the package description mentions both the
  known-noise `ebnm` workflow and the learned-noise `eb_smoother()` workflow.
- Update roxygen for:
  - `eb_smoother()`
  - `ebnm_LGP_generator()`
  - `ebnm_Matern_generator()`
  - `Nonspatial()`
  - `LGP_setup()`
- Update `internal/vignettes/00-overview.Rmd` so it no longer implies all
  workflows require known standard errors.
- Add a short scope note to `01-lgp-workflow.Rmd` and `02-matern-workflow.Rmd`
  clarifying that they focus on known-`s` generator workflows.

## Validation
- Regenerate `man/` and `NAMESPACE` from roxygen comments.
- Re-render the internal vignette HTML outputs.
- Run the targeted `test-eb-smoother.R` tests.
- Smoke-run the two core examples documented in the refreshed narrative:
  - known-`s` generator
  - unknown-`s` `eb_smoother()`
