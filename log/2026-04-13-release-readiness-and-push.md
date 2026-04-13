# Update Log: Release Readiness and Publish Pass

Date: 2026-04-13

## Summary

Performed a final release-readiness sweep for `EBSmoothr` before publishing the
current package state.

## Changes

- Added a repository-level `README.md` with package overview, repository
  structure, installation notes, and quick-start examples.
- Updated `EBSmoothr/DESCRIPTION` to remove the unused `LazyData` setting,
  enable roxygen markdown explicitly, and record the `INLA` additional
  repository.
- Tightened `.gitignore` and `EBSmoothr/.Rbuildignore` so hidden files,
  compiled shared objects, object files, and R check artifacts do not pollute
  the repository or source package.
- Converted `EBSmoothr/LICENSE` into a valid R-package MIT license stub and
  added `EBSmoothr/LICENSE.md` with the full license text.
- Removed hidden local files and compiled artifacts that were causing package
  check warnings and notes.

## Verification

- `R CMD check --no-manual EBSmoothr_0.1.0.tar.gz`

## Notes

- The repository currently tracks major research documentation in top-level
  `internal/`, `plan/`, and `log/` directories rather than shipping those
  materials inside the installable package.
