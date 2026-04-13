## Objective

Move the internal vignette materials out of the package `inst/` directory and
save rendered outputs alongside the source documents.

## Plan

1. Check whether `inst/` is ignored by Git and whether leaving the documents there would affect package contents.
2. Move the internal vignette sources to a top-level repository folder.
3. Update the vignette setup so the documents can render from source without requiring a prior package installation.
4. Add a small render script that rebuilds all internal vignette outputs.
5. Render the HTML outputs and save them in a dedicated subfolder.
6. Record the change in the project log for other collaborators.
