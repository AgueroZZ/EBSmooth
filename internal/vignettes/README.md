# Internal Vignettes

This folder contains internal vignette drafts and rendered outputs for
`EBSmoothr`.

These documents are meant for collaborators who want longer, tutorial-style
examples without registering them as formal package vignettes or shipping them
inside the installed package payload.

## Files

- `00-overview.Rmd`: package overview and model-selection guidance.
- `01-lgp-workflow.Rmd`: one-dimensional smoothing with the L-GP prior.
- `02-matern-workflow.Rmd`: one-dimensional and two-dimensional Matern
  smoothing, including both the exact EB backend and the optional INLA
  PC-prior backend.
- `03-package-sanity-check.Rmd`: a light end-to-end validation document that
  exercises the main public package paths and checks a few key invariants.
- `render_internal_vignettes.R`: render all vignette sources to HTML.
- `rendered/`: saved HTML outputs produced from the `.Rmd` sources.
