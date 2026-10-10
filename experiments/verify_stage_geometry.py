#!/usr/bin/env python3
"""Independent checks of the active-stage PH extension's algebra and limits.

The exact flow uses scipy.linalg.expm; the compression calculation and
linear-fractional structure reduction do not use the Newton implementation.
The negative examples are expected scientific outcomes, not accepted failures.
"""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import block_diag, expm

from ph_gauges import PHModel, common_gauge, scheme


ROOT = Path(__file__).resolve().parent


def main():
    rng = np.random.default_rng(20261010)
    checks = {}
    result = {"parameters": {"seed": 20261010}, "checks": checks}

    # Full state+port reduction: preserve the output and induced feedthrough.
    s, n, m = 2, 3, 2
    b = np.array([0.5, 0.5])
    raw = rng.standard_normal((n, n))
    J = raw - raw.T
    R = np.diag([0.0, 0.2, 0.7])
    B = rng.standard_normal((n, m))
    N = np.array([[0.0, 0.3], [-0.3, 0.0]])
    Ls = block_diag(*([J - R] * s))
    Rs = block_diag(*([R] * s))
    Bs = block_diag(*([B] * s))
    Ns = block_diag(*([N] * s))
    raw = rng.standard_normal((s*n, s*n))
    K = 0.08 * (raw - raw.T)
    wroot = np.sqrt(np.repeat(b, n))
    Kh = K / wroot[:, None] / wroot[None, :]
    Z = np.block([[Ls, Bs], [-Bs.T, -Ns]])
    D = block_diag(Kh, np.zeros((s*m, s*m)))
    T = np.eye(s*(n+m)) + Z @ D
    Zeff = np.linalg.solve(T, Z)
    Rfull = block_diag(Rs, np.zeros((s*m, s*m)))
    expected_sym = -np.linalg.solve(T, Rfull) @ np.linalg.inv(T).T
    congruence_error = np.linalg.norm((Zeff + Zeff.T)/2 - expected_sym)
    g = rng.standard_normal(s*(n+m))
    f = Zeff @ g
    e = g - D @ f
    E, F = e[:s*n], f[:s*n]
    U, Y = g[s*n:], -f[s*n:]
    original_error = np.linalg.norm(f - Z @ e)
    power_error = abs(g @ f + E @ Rs @ E)
    # Retaining the original output after eliminating the gauge is incorrect.
    wrong_Y = Bs.T @ g[:s*n] + Ns @ U
    wrong_power = abs(g[:s*n] @ F - U @ wrong_Y + E @ Rs @ E)
    checks["full_port_congruence"] = congruence_error < 5e-13
    checks["full_port_reconstruction"] = original_error < 5e-13
    checks["full_port_power"] = power_error < 5e-13
    checks["retaining_old_output_is_wrong"] = wrong_power > 1e-4
    checks["effective_dissipation_psd"] = np.linalg.eigvalsh(-0.5*(Zeff+Zeff.T))[0] > -1e-12
    result["port_elimination"] = {
        "congruence_error": float(congruence_error),
        "reconstruction_error": float(original_error),
        "power_error": float(power_error),
        "old_output_power_error": float(wrong_power),
        "condition_number": float(np.linalg.cond(T)),
    }

    # An exact singular conservative gauge, and removal by strict resistance.
    J2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
    singular = np.eye(2) + J2 @ J2
    strict = np.eye(2) + (J2 - 0.2*np.eye(2)) @ J2
    checks["conservative_singularity_counterexample"] = np.linalg.norm(singular) == 0
    checks["strict_resistance_regularizes_example"] = abs(np.linalg.det(strict)) > 0.01

    # Dissipation variance and the minimum number of linear square channels.
    _, bw = leggauss(3)
    bw = bw / 2
    E3 = rng.standard_normal((3, n))
    mean = bw @ E3
    full_diss = np.einsum("i,ij,jk,ik->", bw, E3, R, E3)
    variance = np.einsum("i,ij,jk,ik->", bw, E3-mean, R, E3-mean)
    V = np.kron(np.diag(bw)-np.outer(bw, bw), R)
    variance_error = abs(full_diss - mean @ R @ mean - variance)
    form_error = abs(E3.ravel() @ V @ E3.ravel() - variance)
    rank = int(np.linalg.matrix_rank(V, tol=1e-12))
    checks["variance_identity"] = variance_error < 1e-13
    checks["variance_quadratic_form"] = form_error < 1e-13
    checks["variance_rank"] = rank == 4
    result["dissipation_variance"] = {
        "variance": float(variance), "identity_error": float(variance_error),
        "quadratic_form_error": float(form_error), "rank": rank,
        "expected_rank": 4,
    }

    # Same physical oscillator, one cap effort vs two active efforts.
    hs = [0.1, 0.05, 0.025, 0.0125]
    x0 = np.array([1.0, 0.0])
    oscillator = []
    for gamma in (0.0, 0.2):
        L = J2 - gamma*np.eye(2)
        exact = expm(L) @ x0
        for method in ("single_cap", "two_active"):
            rows = []
            model = PHModel(
                H=lambda x: float(x@x/2), grad=lambda x: x,
                hess=lambda x: np.eye(2), J=J2, R=gamma*np.eye(2))
            # The nonzero h^4 common gauge is retained in the active method.
            itg = scheme(s=2, q=2, gauge=common_gauge(0.3*J2),
                         atol=1e-15, rtol=1e-14)
            for h in hs:
                x = x0.copy()
                maximum_balance = 0.0
                count = int(round(1/h))
                S = h*J2/12
                M = np.linalg.solve(np.eye(2)+h*L@S-h*L/2,
                                    np.eye(2)+h*L@S+h*L/2)
                for k in range(count):
                    if method == "single_cap":
                        new = M @ x
                        effort = (x+new)/2-S@(new-x)
                        balance = (new@new-x@x)/2+h*gamma*(effort@effort)
                    else:
                        new, diag = itg.step(model, x, k*h, h)
                        balance = diag["balance"]
                    maximum_balance = max(maximum_balance, abs(balance))
                    x = new
                error = float(np.linalg.norm(x-exact))
                row = {"h": h, "error": error, "max_power_defect": float(maximum_balance)}
                if rows:
                    row["order"] = float(np.log2(rows[-1]["error"]/error))
                rows.append(row)
            target = 2 if gamma and method == "single_cap" else 4
            checks["oscillator_{}_gamma_{}".format(method, gamma)] = (
                abs(rows[-1]["order"]-target) < 0.02
                and max(row["max_power_defect"] for row in rows) < 1e-12)
            oscillator.append({"gamma": gamma, "method": method, "expected_order": target, "rows": rows})
    result["oscillator"] = oscillator

    # Evaluate the continuous variance on exact endpoints, independently.
    L = J2 - 0.2*np.eye(2)
    nodes, weights = leggauss(12)
    variance_rows = []
    coefficient = 0.2*np.linalg.norm(L@x0)**2/12
    for h in (0.1, 0.05, 0.025, 0.0125):
        times, weights_h = h*(nodes+1)/2, h*weights/2
        values = np.array([expm(t*L)@x0 for t in times])
        mean = (weights_h@values)/h
        variance = float(0.2*np.einsum("i,ij,ij->", weights_h, values-mean, values-mean))
        endpoint = expm(h*L)@x0
        cap_defect = float((endpoint@endpoint-x0@x0)/2 + h*0.2*(mean@mean))
        variance_rows.append({"h": h, "variance_over_h3": variance/h**3,
                              "expected_limit": coefficient,
                              "exact_endpoint_identity_error": abs(cap_defect+variance)})
    checks["continuous_variance_cubic_coefficient"] = abs(variance_rows[-1]["variance_over_h3"]/coefficient-1) < 0.01
    checks["continuous_variance_exact_identity"] = max(r["exact_endpoint_identity_error"] for r in variance_rows) < 1e-13
    result["continuous_variance"] = variance_rows

    # A gauge with correct skewness and reversal can still reduce order.
    scaling_rows = []
    for exponent in (2, 4):
        def gauge(ctx, exponent=exponent):
            return ctx["h"]**exponent*np.kron(np.diag(ctx["b"]), J2)
        itg = scheme(s=2, q=2, gauge=gauge, atol=1e-15, rtol=1e-14)
        model = PHModel(lambda x: float(x@x/2), lambda x: x,
                        lambda x: np.eye(2), J2)
        rows = []
        for h in hs:
            x = x0.copy()
            max_energy = 0.0
            for k in range(int(round(1/h))):
                x, diag = itg.step(model, x, k*h, h)
                max_energy = max(max_energy, abs(model.H(x)-0.5))
            error = float(np.linalg.norm(x-expm(J2)@x0))
            row = {"h": h, "error": error, "max_energy_drift": max_energy}
            if rows:
                row["order"] = float(np.log2(rows[-1]["error"]/error))
            rows.append(row)
        checks["gauge_scaling_h{}".format(exponent)] = (
            abs(rows[-1]["order"]-exponent) < 0.02 and
            max(row["max_energy_drift"] for row in rows) < 1e-12)
        scaling_rows.append({"exponent": exponent, "rows": rows})
    result["insufficient_scaling_counterexample"] = scaling_rows

    # Ordinary bools make the stored evidence strict JSON.
    result["checks"] = {key: bool(value) for key, value in checks.items()}
    result["passed"] = sum(result["checks"].values())
    result["total"] = len(checks)
    result["all_passed"] = all(result["checks"].values())
    destination = ROOT/"results"/"stage_geometry_validation.json"
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    assert result["all_passed"], result["checks"]
    print("ALL {} STAGE GEOMETRY CHECKS PASSED".format(len(checks)))
    print("Single-effort dissipative order:", oscillator[2]["rows"][-1]["order"])
    print("Active-stage dissipative order:", oscillator[3]["rows"][-1]["order"])
    print("Port congruence error:", congruence_error)


if __name__ == "__main__":
    main()
