## Objective

Add a lightweight internal vignette that exercises the main public `EBSmoothr`
interfaces and records a simple end-to-end sanity check for collaborators.

## Plan

1. Add a numbered internal vignette so it participates in the existing
   rendering script.
2. Cover the main public fitting paths:
   - L-GP identity;
   - L-GP log link;
   - L-GP fixed-beta mode;
   - exact Matern in 1D and 2D;
   - Matern with `pc.penalty`, including fixed-parameter mode.
3. Include simple assertions in the vignette so rendering fails if key object
   structure or objective relationships break.
4. Render the vignette to HTML and update the internal vignette index.
5. Record the change in `log/` for other collaborators.
