# Space-Time Stokes-Dirac Discretization

[![Validation](https://github.com/redhamoulla/space-time-stokes-dirac/actions/workflows/validate.yml/badge.svg)](https://github.com/redhamoulla/space-time-stokes-dirac/actions/workflows/validate.yml)

This repository contains the manuscript and reproducibility package for:

> **Space-Time Stokes-Dirac Discretization: Affine Pairing Gauges, Constitutive Cayley Coordinates, and Active Interior Reduction**  
> Redha Moulla, AXIA, France (2026)

The paper develops a finite-dimensional theory of the freedom left by a discrete Stokes identity and follows that freedom through space-time reduction, constitutive closure, symplecticity, order conditions, and interior realization.

## Main results

- Compatible pairings for a possibly nonsurjective differential form an affine family over an alternating quotient space.
- Temporal, spatial, and genuinely mixed gauges are separated on the discrete cylinder.
- Cap gauges become shifted Cayley coordinates for quadratic Hamiltonian systems, with an explicit symplectic locus and defect identity.
- Active polynomial time cells generate Gauss-Legendre/Padé cap maps by interior elimination.
- A rectangular Stokes theorem reconstructs an open Dirac representation for the 2012 transmission-line discretization.
- Mixed cylinder gauges yield reversible, exactly power-preserving deformations and a recursive higher-order hierarchy.
- Every rank-`2q` alternating cap deformation has a minimal `2q`-dimensional auxiliary Schur generator.
- Gauss coordinates admit finite positive Padé-Stieltjes/Jacobi realizations; the corrected hierarchy has an explicit positivity obstruction.

The distinction between algebraic generation, passive realization, and selection by independently prescribed field dynamics is maintained throughout.

## Repository layout

```text
.
├── manuscript/
│   ├── moulla_2026_space_time_stokes_dirac.tex
│   ├── moulla_2026_space_time_stokes_dirac.pdf
│   └── README.md
├── experiments/
│   ├── verify_transmission_line_v12.py
│   ├── README.md
│   └── results/
│       └── verification_report_v12.txt
├── docs/
│   ├── reproducibility.md
│   └── version-12-summary.md
├── .github/workflows/validate.yml
├── CITATION.cff
├── CHANGELOG.md
├── Makefile
└── requirements.txt
```

## Read the paper

- [Compiled manuscript](manuscript/moulla_2026_space_time_stokes_dirac.pdf)
- [LaTeX source](manuscript/moulla_2026_space_time_stokes_dirac.tex)

## Reproduce the numerical checks

Python 3.10 or later is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python experiments/verify_transmission_line_v12.py
```

A successful run ends with:

```text
ALL 28 ASSERTION GROUPS PASSED
```

The archived reference output is available in [`experiments/results`](experiments/results/verification_report_v12.txt). See [the reproducibility notes](docs/reproducibility.md) for the scope of each check.

## Build the manuscript

A standard TeX Live installation with `latexmk` is sufficient:

```bash
make paper
```

The generated PDF is written to `build/`.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). GitHub can export them directly in BibTeX or plain-text form through the repository's **Cite this repository** menu.

## License

Copyright © 2026 Redha Moulla. All rights reserved. No reuse license has been granted yet; see [`LICENSE`](LICENSE).

