## Objective

Make the exploratory script compare three distinct workflows:

1. Integrated spatial flash decomposition.
2. Standard flash decomposition followed by spatial post-smoothing of factors.
3. Standard flash decomposition without spatial smoothing.

## Plan

1. Add helpers for post-smoothing flash factors with `ebnm_Matern_generator()`.
2. Optionally refit loadings after factor smoothing so the two-stage method is
   compared as a full decomposition rather than factor denoising alone.
3. Reuse the existing normalization and plotting helpers to visualize all three
   workflows with the same alignment logic.
