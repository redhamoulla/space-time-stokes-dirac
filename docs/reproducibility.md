# Reproducibility notes

## Environment

The reference run used Python 3.12, NumPy 2.3, and SciPy 1.17. The code uses only public NumPy and SciPy APIs and is expected to work with Python 3.10+, NumPy 1.24+, and SciPy 1.10+.

## One-command checks

```bash
make verify
```

The command must terminate normally and print:

```text
ALL 28 ASSERTION GROUPS PASSED
```

The tests use explicit tolerances. They do not infer theorem statements from plots or from visually fitted slopes.

## What is verified

The executable checks the finite-dimensional algebra used in the paper: Stokes identities, skew and energy defects, exact Schur reductions, convergence branches, pole distinctions, multiprecision commutator asymptotics, dispersion selection, and the two-cell power balance.

The newly emphasized interior-generation results are tested directly:

- a rank-four alternating cap change is recovered from four auxiliary alternating modes;
- the augmented skew defect is zero to roundoff;
- the Schur error is approximately `8.02e-16`;
- the eliminated-response error is approximately `1.83e-15`;
- the Gauss-3 coordinate is reproduced by `w = 1/12` and `lambda = 1/60`;
- the corrected `(s,r) = (2,1)` quotient is exactly `-1/720` at `u = 61`.

## What is not claimed

These checks validate the stated finite-dimensional identities and examples. They are not a proof of simultaneous space-time convergence to an infinite-dimensional PDE, nor do they establish selection of a gauge by field dynamics prescribed independently of the target gauge.

## Manuscript build

```bash
make paper
```

The build requires `latexmk` and a reasonably complete TeX Live distribution. The compiled reference PDF is committed under `manuscript/`.

