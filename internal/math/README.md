# Internal Math Notes

This folder contains internal mathematical notes for `EBSmoothr`.

The goal is to document the core smooth-EBNM modeling framework, explain how
the package instantiates that framework with the L-GP and Matern priors, and
record how the current code implements the main empirical-Bayes targets across
the exact and optional INLA PC-prior backends.

## Files

- `ebnm-smooth-foundations.md`: a compact math note connecting classical EBNM,
  smooth EBNM, and the current package implementation.
