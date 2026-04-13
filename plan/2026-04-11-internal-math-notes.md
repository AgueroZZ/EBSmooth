## Objective

Create an internal math note that explains the smooth-EBNM framework behind
`EBSmoothr` and maps the mathematics to the current L-GP and Matern
implementations.

## Plan

1. Review the official EBNM references and the supplied smooth-EBNM formula.
2. Restate the smooth-EBNM model in notation that matches the package code.
3. Explain how the L-GP implementation fits the smooth-EBNM framework.
4. Explain how the Matern implementation fits the smooth-EBNM framework.
5. Record the main mismatches between the ideal math and the current code so
   future code changes can target them directly.
6. Store the note in a dedicated `internal/math` folder for future collaborators.
