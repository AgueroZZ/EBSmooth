# EB-FPCA Background Learning

## Summary

Started the `EBFPCA/` workspace by copying locally stored Zotero PDF
attachments for the background papers and creating the first version of an
iterative EB-FPCA knowledge document.

## Source Handling

- Located the Zotero screenshot at `EBFPCA/related-paper/list-of-paper.png`.
- Queried the local Zotero database in read-only immutable mode because the
  live database was locked by Zotero.
- Copied 12 available PDF attachments from `/Users/ziangzhang/Zotero/storage/`
  into `EBFPCA/related-paper/`.
- Confirmed two screenshot-listed papers did not currently have Zotero PDF
  attachments:
  - van der Linde (2008), `Variational Bayesian functional PCA`, had a local
    ScienceDirect HTML snapshot and Zotero full-text cache.
  - Goldsmith, Zipunnikov, and Schrack (2015), `Generalized Multilevel
    Function-on-Scalar Regression and Principal Component Analysis`, had
    Zotero metadata and abstract but no attachment.

## Knowledge Document

Created `EBFPCA/knowledge.md` with:

- source inventory and copied-PDF status;
- common sparse/Bayesian FPCA model template;
- directly usable software leads;
- open-data and simulation-resource ranking for benchmark planning;
- paper-level takeaways;
- implications for the current `EBSmoothr` L-GP and Matern smoother framework;
- recommended MVP path for covariance-first and joint Gaussian EB-FPCA;
- future extensions and update protocol.

## 2026-06-08 Update: Open Data and Simulation Resources

Added a benchmark-oriented section to `EBFPCA/knowledge.md`.

Immediate benchmark candidates:

- BayesTime / Jiang et al. (2020): open R package plus manuscript analysis
  repository with application and simulation folders.
- mSFPCA / Jiang et al. (2022): open analysis repository with data,
  applications, and simulation folders; better for later multivariate work.
- Canadian weather data: open `fda` package data used by multiple FPCA papers;
  good for a small deterministic dense-to-sparse benchmark.
- UCI EEG Database / Lin et al. (2016): open dense EEG time series; useful
  after preprocessing and benchmark metrics are stable.

Caveated resources:

- FAST / Sartini et al. (2026): open Stan and simulation routines, but DASH4D
  CGM application data are not yet public.
- Boland et al. (2023): open R code and simulated-data tutorial for posterior
  envelope workflows, not a ready real-data benchmark.
- Goldsmith et al. (2015): public Stan/R code, but application accelerometry
  data are not clearly ready for direct reuse.
- Sang et al. (2025): open R simulation code for informative observation times,
  but the CD4/viral-load application data are not clearly public.

## Validation

- Confirmed `EBFPCA/related-paper/` contains 12 copied PDF files plus the
  screenshot.
- Used `pypdf` from the bundled Codex Python runtime to extract text from the
  copied PDFs.
- Inspected the current package README and core R files to avoid writing a
  generic literature summary disconnected from the existing codebase.
- Did not modify package source, tests, documentation, or existing dirty
  worktree files.
