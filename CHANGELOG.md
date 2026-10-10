# Changelog

## Unified GitHub manuscript — 2026-10-10

- The complete version-12 theory and the active-stage extension are presented in one English `docs/manuscript.md`, with a linked table of contents, proofs, equations, numerical evidence and references.
- A concise README, a scientific overview and vector SVG figures follow the presentation of the companion spatial-discretization and wavelet-clustering repositories.
- Numerical evidence is grouped in `results/`, with CSV tables, SVG figures and a source/result hash manifest; reproduction commands write to these paths.
- The original LaTeX source is retained under `archive/version-12/`. The PDF and superseded split documentation are removed from the current tree and remain available in Git history.
- CI validates the executable experiments without requiring a PDF build.

## Active-stage PH extension — 2026-10-10

- A single research document extends the gauge framework to active temporal efforts, nonlinear Hamiltonians and dissipative/open PH systems.
- An exact elimination result recovers the version-12 mixed cap family in the quadratic conservative case.
- A variance obstruction explains the loss of high order under a physical one-effort dissipative compression.
- A reusable projected-gradient stage solver supports alternating stage gauges and reports separate power, quadrature and nonlinear-solver certificates.
- Reproducible checks cover the cap/stage equivalence, coordinate covariance, passive port elimination, temporal order, nonlinear wave dynamics and forced dissipation.
- The original version-12 manuscript and verification program are unchanged.

## Version 12 — 2026-08-18

Initial public repository release.

- Final 55-page manuscript and self-contained LaTeX source.
- Transmission-line verification program with 28 guarded assertion groups.
- Archived numerical report and reproducibility documentation.
- Continuous integration for the numerical checks and manuscript build.


