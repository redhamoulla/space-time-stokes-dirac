# Space-Time Stokes–Dirac Discretization

[![Validation](https://github.com/redhamoulla/space-time-stokes-dirac/actions/workflows/validate.yml/badge.svg)](https://github.com/redhamoulla/space-time-stokes-dirac/actions/workflows/validate.yml)

**Geometric discretization of port-Hamiltonian systems in space and time.**  
Redha Moulla · 2026

This research repository connects the 2012 geometric spatial discretization to a temporal construction based on Stokes–Dirac pairing gauges and active polynomial stages. The construction covers conservative quadratic systems and extends to nonlinear systems with boundary inputs and dissipation.

**[Read the manuscript →](docs/manuscript.md)** — one continuous article with the theory, proofs, simulations, limitations and references, rendered directly on GitHub.

![Spatial Stokes discretization, active temporal stages and port-Hamiltonian closure produce a discrete balance of stored energy, boundary work and dissipation.](docs/assets/overview.svg)

## What the construction preserves

The active temporal efforts satisfy the discrete power balance

```math
H(x_{n+1})-H(x_n)=h\sum_{i=1}^{s}b_i\left(u_i^{\mathsf T}y_i-e_i^{\mathsf T}R_i e_i\right).
```

This identity is exact under the stated discrete chain rule and constitutive assumptions; quadrature and nonlinear-solver errors are measured explicitly in the implementation.

- **Spatial structure.** A rectangular Stokes identity reconstructs the open Dirac representation, including its boundary ports.
- **Temporal freedom.** Alternating pairing gauges generate energy-compatible schemes; a common-stage gauge recovers the original mixed cap construction exactly in the quadratic conservative case.
- **Nonlinear dynamics and dissipation.** Retaining active temporal efforts preserves the dissipation variance that a single-effort compression loses. Under the stated scaling and regularity assumptions, the gauge retains the base method’s temporal order.

The zero-gauge projected-gradient integrator builds on existing methods. The manuscript distinguishes those foundations from the gauge construction and states the scope of each result.

## Repository guide

| Location | Contents |
| --- | --- |
| [docs/manuscript.md](docs/manuscript.md) | Complete article, including the original theory and its active-stage extension |
| [experiments/](experiments/) | Reusable temporal solver, algebra checks and transmission-line simulations |
| [results/](results/) | Recorded validation reports, JSON results, CSV tables and SVG figures |
| [archive/version-12/manuscript.tex](archive/version-12/manuscript.tex) | Original version-12 LaTeX source, retained for provenance |
| [CITATION.cff](CITATION.cff) | Citation metadata |

## Reproduce the results

Use Python 3.10 or later. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
make verify
```

The full run executes the **28 original assertion groups**, **14 temporal-engine checks**, **18 independent geometry checks**, and the nonlinear transmission-line simulation assertions. To run only the active-stage extension:

```bash
make verify-extension
```

| Script | Validation |
| --- | --- |
| [verify_transmission_line_v12.py](experiments/verify_transmission_line_v12.py) | Original Stokes, Dirac, cap-map and order identities |
| [ph_gauges.py](experiments/ph_gauges.py) | Temporal engine, coordinate covariance, reversibility and exact cap bridge |
| [verify_stage_geometry.py](experiments/verify_stage_geometry.py) | Port elimination, dissipation variance, order obstruction and gauge scaling |
| [validate_ph_extension.py](experiments/validate_ph_extension.py) | Nonlinear closed, forced and dissipative transmission lines; convergence and counterexamples |

The scripts regenerate their evidence in `results/`, with tables in `results/tables/` and figures in `results/figures/`. The recorded experiments support fourth-order temporal convergence on the tested forced dissipative problem. They also exhibit why energy conservation alone does not establish trajectory accuracy. A general simultaneous space–time PDE convergence theorem remains open.

## Citation and license

Use [CITATION.cff](CITATION.cff) or GitHub’s **Cite this repository** menu. Copyright © 2026 Redha Moulla. All rights reserved; see [LICENSE](LICENSE).
