# Numerical verification

`verify_transmission_line_v12.py` is the executable verification companion to the manuscript. It contains 28 guarded assertion groups and stops with `AssertionError` if a structural identity, tolerance, or convergence law fails.

The checks cover:

- the rectangular spatial Stokes identity and reconstructed open Dirac chart;
- trace-generated and mixed volume-to-cap reductions;
- exact power balance, dissipation, and interface cancellation;
- Gauss-stage condensation and Cayley-coordinate reconstruction;
- the one-term and recursive higher-order gauge hierarchies;
- noncommutation with the scalar functional calculus;
- stable evaluation near the identity and pole classification;
- minimal alternating Schur generation;
- the Gauss-3 Stieltjes realization and the corrected-hierarchy positivity obstruction.

Run from the repository root:

```bash
python -m pip install -r requirements.txt
python experiments/verify_transmission_line_v12.py
```

The reference console output is archived in `results/verification_report_v12.txt`.

