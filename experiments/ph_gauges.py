"""Active-stage projected-gradient PH integration with weighted-skew gauges.

Model convention::

    xdot = (J(x) - R(x)) grad H(x) + B(x) u(t)
    y    = B(x).T grad H(x) + S(x) u(t)

The numerical output uses the gauged stage effort, not grad H at a
collocation node. This distinction is essential for the power identity.

For a Dirac feedthrough S.T=-S, u.T S u=0.  The implementation also admits
resistive feedthrough: writing S=N+S_a, N=(S+S.T)/2, the reported dissipation
is e.T R e + u.T N u. Passivity then requires sym(R)>=0 and N>=0, in addition
to skew J. No positivity is inferred from the name of a model attribute.

For s Gauss collocation nodes and q Gauss integration nodes on [0, 1],
X(tau) = x0 + h sum_j A_j(tau) F_j. Set
G_i = (1/b_i) integral ell_i(tau) grad H(X(tau)) d tau, W=diag(b) tensor I,
and E=G-W^{-1} K_h F with K_h.T=-K_h. The gauge makes no contribution
to contracted power: F.T W E=F.T W G.  The integral is approximated
with q-point Gaussian quadrature. If H has polynomial
degree nu, q >= ceil(nu*s/2) makes the chain rule exact in exact arithmetic.
For nonpolynomial H, ``chain_defect`` explicitly measures its quadrature error.

This implementation needs only NumPy and is compatible with Python >=3.8.
It deliberately does not silently accept an unconverged nonlinear solve.

The ungauged construction is the established projected-gradient variant
of collocation (Hairer, 2010; Celledoni and Hoiseth, 2017). We do not claim
novelty for that temporal integrator. The additional gauge below retains
the same active temporal cell and exposes its weighted-skew algebra.
"""

from dataclasses import dataclass
from typing import Callable, Optional, Any

import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import leggauss


@dataclass
class PHModel:
    """Convenience model; a duck-typed object with these attributes also works.

    J, R, B and S are arrays or functions of x.  Missing R, B and S are
    zero.  A vector B denotes one port.  u(t) may return a scalar for one
    port. hess=None selects numerical differentiation of the stage residual.
    S may be skew (Dirac feedthrough) or include a positive symmetric
    resistive feedthrough; the latter contributes to ``diss`` explicitly.
    """

    H: Callable
    grad: Callable
    hess: Optional[Callable]
    J: Any
    R: Any = None
    B: Any = None
    S: Any = None
    u: Optional[Callable] = None


class NewtonConvergenceError(RuntimeError):
    """Raised when the stage equations cannot be solved to the stated tolerance."""

    def __init__(self, message, diagnostics):
        super().__init__(message)
        self.diagnostics = diagnostics


@dataclass
class CommonStageGauge:
    """Constant-state gauge K_h=h**exponent * (b b.T) tensor Theta.

    Theta is a d-by-d alternating matrix in state coordinates. If exponent
    is None, use 2*s. The resulting correction is the same at every stage:
    E_i-G_i=-h**exponent Theta sum_j b_j F_j.
    An even exponent is reversal-equivariant; exponent >=2*s retains the
    base order under smoothness and a locally regular stage solve. A lower
    or odd exponent is allowed deliberately, for verification experiments.
    """

    theta: Any
    exponent: Optional[int] = None

    def __post_init__(self):
        self.theta = np.array(self.theta, dtype=float, copy=True)
        if self.theta.ndim != 2 or self.theta.shape[0] != self.theta.shape[1]:
            raise ValueError("Theta must be a square state-space matrix")
        if not np.all(np.isfinite(self.theta)):
            raise ValueError("Theta contains nonfinite values")
        if self.exponent is not None:
            if int(self.exponent) != self.exponent or self.exponent < 0:
                raise ValueError("Gauge exponent must be a nonnegative integer")
            self.exponent = int(self.exponent)

    def __call__(self, context):
        d = len(context["x0"])
        if self.theta.shape != (d, d):
            raise ValueError("Theta dimension must equal the state dimension")
        b, h = context["b"], context["h"]
        exponent = 2*len(b) if self.exponent is None else self.exponent
        return h**exponent * np.kron(np.outer(b, b), self.theta)


def common_gauge(theta, exponent=None):
    """Construct a constant-state common-effort gauge; see CommonStageGauge."""
    return CommonStageGauge(theta, exponent)


def _matrix(value, x, shape, name):
    if value is None:
        return np.zeros(shape, dtype=float)
    a = np.asarray(value(x) if callable(value) else value, dtype=float)
    if name == "B" and a.ndim == 1:
        a = a[:, None]
    if name == "S" and a.ndim == 0 and shape == (1, 1):
        a = a.reshape(1, 1)
    if a.shape != shape:
        raise ValueError("{} has shape {}; expected {}".format(name, a.shape, shape))
    if not np.all(np.isfinite(a)):
        raise ValueError("{} contains nonfinite values".format(name))
    return a


class GaugedProjectedGradient:
    """Gauss projected-gradient scheme, with independent chain quadrature.

    s=1,q=1 is implicit midpoint (not energy-exact for general H).
    s=1,q=2 is AVF for Hamiltonians of degree <=4.
    s=2,q=4 is a fourth-order projected-gradient scheme, exact in energy
    for degree <=4 H and skew J, when the stage solve is exact. These order
    statements apply to gauge=None and gauges O(h**(2*s)), with bounded
    smooth coefficients and a locally regular stage root. Preservation of
    energy alone does not establish order, solvability, or symplecticity.

    gauge is None, a stage matrix of shape (s*d,s*d), or a callback
    gauge(context). The context contains F,X,Xq,x0,x1,t,h,c,b,W,model.
    Stage-major stacking is used. A callback can depend on the entire
    reconstructed path; it is differentiated numerically. Constant arrays
    and CommonStageGauge have an analytic derivative in the stage flows.
    Skewness is checked without repairing the supplied matrix. Set
    check_skew=False only to measure deliberately incompatible gauges.

    The convergence threshold is atol+rtol*max(1, ||F||_inf, ||rhs||_inf).
    Both the actual residual and this threshold are returned.  The power
    certificate uses the signed residual; ``residual_power_bound`` is its
    absolute-term bound and is useful even when cancellation occurs.
    """

    def __init__(self, s=2, q=4, gauge=None, atol=2e-14, rtol=2e-13,
                 maxiter=24, check_skew=True, skew_tolerance=2e-12):
        if int(s) != s or int(q) != q or s < 1 or q < s:
            raise ValueError("s and q must be positive integers with q >= s")
        self.s, self.q = int(s), int(q)
        self.gauge = gauge
        self.check_skew = bool(check_skew)
        self.skew_tolerance = float(skew_tolerance)
        self.atol, self.rtol, self.maxiter = float(atol), float(rtol), int(maxiter)
        if self.atol < 0 or self.rtol < 0 or self.maxiter < 1:
            raise ValueError("Invalid Newton tolerance or iteration limit")
        c, b = leggauss(self.s)
        tau, w = leggauss(self.q)
        self.c, self.b = (c + 1.0) / 2.0, b / 2.0
        self.tau, self.w = (tau + 1.0) / 2.0, w / 2.0
        polynomials = []
        for i in range(self.s):
            p = Polynomial([1.0])
            for j in range(self.s):
                if i != j:
                    p *= Polynomial([-self.c[j], 1.0]) / (self.c[i] - self.c[j])
            polynomials.append(p)
        integrals = [p.integ() for p in polynomials]
        self.ell = np.array([p(self.tau) for p in polynomials]).T
        self.A = np.array([p(self.c) - p(0.0) for p in integrals]).T
        self.Aq = np.array([p(self.tau) - p(0.0) for p in integrals]).T
        self.W = self.ell.T * self.w[None, :] / self.b[:, None]

    def step(self, model, x, t, h, F0=None):
        """Return ``(xnew, diagnostics)`` for one step, or raise on failure.

        Optional F0 has shape (s, dim(x)) and provides a stage predictor.
        Key signed identity (up to reported floating-point closure error)::

            balance = deltaH - work + diss
                    = chain_defect + gauge_power + structure_defect
                      + residual_power_certificate

        ``diss`` equals ``state_diss`` plus ``feedthrough_diss``. The latter
        is h*sum b_i u_i.T S_i u_i, zero for skew S in exact arithmetic.
        ``structure_defect`` is h*sum b_i e_i.T J_i e_i and should vanish
        for skew J. ``chain_defect`` includes quadrature and roundoff error.
        ``gauge_power`` is h*F.T K_h F, zero for an alternating gauge in
        exact arithmetic. ``G`` is the projected constitutive effort and
        ``e`` is the actual gauged effort used in the dynamics and ports.
        No positivity, skewness, or quadrature exactness is silently assumed.
        """
        x = np.asarray(x, dtype=float)
        if x.ndim != 1 or not np.all(np.isfinite(x)):
            raise ValueError("x must be a finite one-dimensional state")
        if not np.isfinite(h) or h == 0 or not np.isfinite(t):
            raise ValueError("t and h must be finite, and h must be nonzero")
        d = x.size
        s = self.s
        Jv = getattr(model, "J")
        Rv = getattr(model, "R", None)
        Bv = getattr(model, "B", None)
        Sv = getattr(model, "S", None)
        if Bv is None:
            p = 0
        else:
            b0 = np.asarray(Bv(x) if callable(Bv) else Bv, dtype=float)
            if b0.ndim == 1:
                b0 = b0[:, None]
            if b0.ndim != 2 or b0.shape[0] != d:
                raise ValueError("B must have shape (state_dimension, port_dimension)")
            p = b0.shape[1]
        uv = getattr(model, "u", None)
        if uv is None:
            us = np.zeros((s, p))
        else:
            us = np.array([np.atleast_1d(uv(t + h * c)) for c in self.c], dtype=float)
            if us.shape != (s, p):
                raise ValueError("u has stage shape {}; expected {}".format(us.shape, (s, p)))
            if not np.all(np.isfinite(us)):
                raise ValueError("u contains nonfinite values")
        constant = not any(callable(v) for v in (Jv, Rv, Bv))
        fixed_gauge = (self.gauge is None or not callable(self.gauge)
                       or isinstance(self.gauge, CommonStageGauge))
        analytic = constant and fixed_gauge and callable(getattr(model, "hess", None))
        wb = np.repeat(self.b, d)
        Wstage = np.diag(wb)
        if constant:
            J0 = _matrix(Jv, x, (d, d), "J")
            R0 = _matrix(Rv, x, (d, d), "R")
            B0 = _matrix(Bv, x, (d, p), "B")
            Js = np.broadcast_to(J0, (s, d, d))
            Rs = np.broadcast_to(R0, (s, d, d))
            Bs = np.broadcast_to(B0, (s, d, p))
            Ks = Js - Rs

        def evaluate(F):
            Xq = x[None, :] + h * self.Aq.dot(F)
            Xi = x[None, :] + h * self.A.dot(F)
            gradients = np.array([model.grad(z) for z in Xq], dtype=float)
            if gradients.shape != (self.q, d):
                raise ValueError("grad H must have state dimension")
            G = self.W.dot(gradients)
            if self.gauge is None:
                gauge_matrix = np.zeros((s*d, s*d))
            elif callable(self.gauge):
                context = {"F": F, "X": Xi, "Xq": Xq, "x0": x,
                           "x1": x+h*self.b.dot(F), "t": t, "h": h,
                           "c": self.c, "b": self.b, "W": Wstage,
                           "model": model}
                gauge_matrix = np.asarray(self.gauge(context), dtype=float)
            else:
                gauge_matrix = np.asarray(self.gauge, dtype=float)
            if gauge_matrix.shape != (s*d, s*d):
                raise ValueError("Gauge matrix must have shape (s*d,s*d)")
            if not np.all(np.isfinite(gauge_matrix)):
                raise ValueError("Gauge matrix contains nonfinite values")
            skew_error = float(np.max(np.abs(gauge_matrix+gauge_matrix.T)))
            gauge_scale = float(np.max(np.abs(gauge_matrix)))
            if self.check_skew and skew_error > self.skew_tolerance*gauge_scale:
                raise ValueError("Gauge is not alternating: relative defect {:.3e}".format(
                    skew_error/max(gauge_scale, np.finfo(float).tiny)))
            e = G - (gauge_matrix.dot(F.reshape(-1))/wb).reshape(s, d)
            if constant:
                jmat, rmat, bmat, kmat = Js, Rs, Bs, Ks
            else:
                jmat = np.array([_matrix(Jv, z, (d, d), "J") for z in Xi])
                rmat = np.array([_matrix(Rv, z, (d, d), "R") for z in Xi])
                bmat = np.array([_matrix(Bv, z, (d, p), "B") for z in Xi])
                kmat = jmat - rmat
            rhs = np.einsum("iab,ib->ia", kmat, e)
            rhs += np.einsum("iap,ip->ia", bmat, us)
            return F - rhs, (Xi, Xq, e, rhs, jmat, rmat, bmat, kmat, G, gauge_matrix)

        if F0 is None:
            j0 = _matrix(Jv, x, (d, d), "J")
            r0 = _matrix(Rv, x, (d, d), "R")
            b0 = _matrix(Bv, x, (d, p), "B")
            fbase = (j0 - r0).dot(np.asarray(model.grad(x), dtype=float))
            F = np.array([fbase + b0.dot(u) for u in us])
        else:
            F = np.asarray(F0, dtype=float).copy()
            if F.shape != (s, d):
                raise ValueError("F0 must have shape (s, state_dimension)")

        nfev, njev, line_search_steps = 0, 0, 0
        norm_history = []
        residual, data = evaluate(F)
        nfev += 1
        for iteration in range(self.maxiter + 1):
            norm = float(np.max(np.abs(residual)))
            scale = max(1.0, float(np.max(np.abs(F))), float(np.max(np.abs(data[3]))))
            threshold = self.atol + self.rtol * scale
            norm_history.append(norm)
            if np.isfinite(norm) and norm <= threshold:
                break
            if iteration == self.maxiter or not np.isfinite(norm):
                raise NewtonConvergenceError(
                    "PH stages failed: residual {:.3e}, tolerance {:.3e}, iterations {}".format(
                        norm, threshold, iteration),
                    {"stage_residual": norm, "threshold": threshold,
                     "iterations": iteration, "residual_history": norm_history})
            if analytic:
                hessians = np.array([model.hess(z) for z in data[1]], dtype=float)
                if hessians.shape != (self.q, d, d):
                    raise ValueError("hess H must be square with state dimension")
                de = h * np.einsum("ik,kj,kab->ijab", self.W, self.Aq, hessians)
                gauge_derivative = (data[9]/wb[:, None]).reshape(s, d, s, d).transpose(0, 2, 1, 3)
                de = de - gauge_derivative
                jac = np.eye(s * d) - np.einsum("iab,ijbc->iajc", data[7], de).reshape(s*d, s*d)
            else:
                flat = F.reshape(-1)
                jac = np.empty((s*d, s*d))
                for k in range(s*d):
                    eps = np.cbrt(np.finfo(float).eps) * max(1.0, abs(flat[k]))
                    plus, minus = flat.copy(), flat.copy()
                    plus[k] += eps
                    minus[k] -= eps
                    rp, _ = evaluate(plus.reshape(s, d))
                    rm, _ = evaluate(minus.reshape(s, d))
                    nfev += 2
                    jac[:, k] = ((rp-rm)/(2*eps)).reshape(-1)
            njev += 1
            try:
                direction = np.linalg.solve(jac, -residual.reshape(-1)).reshape(s, d)
            except np.linalg.LinAlgError as exc:
                raise NewtonConvergenceError("Singular Newton stage Jacobian", {
                    "stage_residual": norm, "iterations": iteration}) from exc
            alpha = 1.0
            accepted = False
            best = None
            for ls in range(16):
                candidate = F + alpha * direction
                rr, dd = evaluate(candidate)
                nfev += 1
                candidate_norm = float(np.max(np.abs(rr)))
                if np.isfinite(candidate_norm):
                    if best is None or candidate_norm < best[0]:
                        best = (candidate_norm, candidate, rr, dd)
                    if candidate_norm <= (1 - 1e-4*alpha) * norm or candidate_norm <= threshold:
                        F, residual, data = candidate, rr, dd
                        accepted = True
                        break
                alpha *= 0.5
                line_search_steps += 1
            if not accepted:
                if best is not None and best[0] < norm:
                    _, F, residual, data = best
                else:
                    raise NewtonConvergenceError("Newton line search could not reduce stage residual", {
                        "stage_residual": norm, "threshold": threshold,
                        "iterations": iteration, "residual_history": norm_history})

        Xi, Xq, es, rhs, jmat, rmat, bmat, kmat, Gs, gauge_matrix = data
        xnew = x + h * self.b.dot(F)
        smat = np.array([_matrix(Sv, z, (p, p), "S") for z in Xi])
        ys = np.einsum("iap,ia->ip", bmat, es) + np.einsum("ipq,iq->ip", smat, us)
        power = np.einsum("ip,ip->i", us, ys)
        state_dissipation = np.einsum("ia,iab,ib->i", es, rmat, es)
        feedthrough_dissipation = np.einsum("ip,ipq,iq->i", us, smat, us)
        dissipation = state_dissipation + feedthrough_dissipation
        work = float(h * self.b.dot(power))
        state_diss = float(h * self.b.dot(state_dissipation))
        feedthrough_diss = float(h * self.b.dot(feedthrough_dissipation))
        diss = float(h * self.b.dot(dissipation))
        H0, H1 = float(model.H(x)), float(model.H(xnew))
        deltaH = H1 - H0
        discrete_chain = float(h * np.einsum("i,ia,ia->", self.b, Gs, F))
        chain_defect = deltaH - discrete_chain
        gauge_power = float(h * F.reshape(-1).dot(gauge_matrix.dot(F.reshape(-1))))
        structure_defect = float(h * np.einsum("i,ia,iab,ib->", self.b, es, jmat, es))
        residual_power = float(h * np.einsum("i,ia,ia->", self.b, es, residual))
        residual_power_bound = float(abs(h) * np.einsum("i,ia,ia->", self.b, np.abs(es), np.abs(residual)))
        balance = deltaH - work + diss
        closure = balance - chain_defect - gauge_power - structure_defect - residual_power
        diagnostics = {
            "H0": H0, "H1": H1, "deltaH": deltaH, "work": work, "diss": diss,
            "state_diss": state_diss, "feedthrough_diss": feedthrough_diss,
            "balance": balance, "chain_defect": chain_defect,
            "gauge_power": gauge_power,
            "gauge_matrix_defect": float(np.max(np.abs(gauge_matrix+gauge_matrix.T))),
            "gauge_matrix_norm": float(np.linalg.norm(gauge_matrix, np.inf)),
            "effort_correction_norm": float(np.max(np.abs(es-Gs))),
            "structure_defect": structure_defect,
            "structure_matrix_defect": float(np.max(np.abs(jmat + jmat.transpose(0, 2, 1)))),
            "residual_power_certificate": residual_power,
            "residual_power_bound": residual_power_bound,
            "certificate_closure": float(closure),
            "stage_residual": float(np.max(np.abs(residual))),
            "residual": residual.copy(), "residual_threshold": float(threshold),
            "iterations": iteration, "nfev": nfev, "njev": njev,
            "line_search_steps": line_search_steps,
            "analytic_jacobian": analytic, "residual_history": norm_history,
            "F": F.copy(), "G": Gs.copy(), "e": es.copy(), "y": ys.copy(), "u": us.copy(),
            "gauge_matrix": gauge_matrix.copy(),
            "X": Xi.copy(), "Xq": Xq.copy(), "c": self.c.copy(), "b": self.b.copy(),
            "stage_power": power, "stage_dissipation": dissipation,
            "stage_state_dissipation": state_dissipation,
            "stage_feedthrough_dissipation": feedthrough_dissipation,
        }
        return xnew, diagnostics


def scheme(s=2, q=4, **kwargs):
    """Construct a GaugedProjectedGradient; kwargs include gauge and tolerances."""
    return GaugedProjectedGradient(s=s, q=q, **kwargs)


def selfcheck():
    """Run deterministic engine checks, including the original v12 cap bridge.

    The reusable integrator above depends only on NumPy. This repository
    selfcheck additionally imports the unchanged v12 verification module.
    Every recorded check raises on failure; no result is merely displayed.
    """
    import json
    from pathlib import Path
    from types import SimpleNamespace
    from verify_transmission_line_v12 import gauss_cap_gauge

    records = []

    def check(name, value, tolerance):
        value = float(value)
        passed = bool(np.isfinite(value) and value <= tolerance)
        records.append({"check": name, "value": value,
                        "tolerance": float(tolerance), "passed": passed})
        assert passed, "{}: {:.6e} > {:.6e}".format(name, value, tolerance)

    j2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
    duffing = PHModel(
        H=lambda x: .5*np.dot(x, x) + .25*x[0]**4,
        grad=lambda x: np.array([x[0] + x[0]**3, x[1]]),
        hess=lambda x: np.diag([1.0 + 3.0*x[0]**2, 1.0]),
        J=j2)
    x0, h = np.array([.7, -.2]), .17
    integrator = scheme(s=2, q=4, gauge=common_gauge(.4*j2))
    x1, diagnostic = integrator.step(duffing, x0, .3, h)
    xback, _ = integrator.step(duffing, x1, .3+h, -h)
    check("nonlinear_common_gauge_reversal", np.linalg.norm(xback-x0), 2e-12)
    check("quartic_energy_exactness", abs(duffing.H(x1)-duffing.H(x0)), 2e-12)
    check("power_certificate_closure", abs(diagnostic["certificate_closure"]), 2e-13)
    numerical_model = PHModel(duffing.H, duffing.grad, None, j2)
    numerical_x, numerical_diagnostic = integrator.step(numerical_model, x0, .3, h)
    assert diagnostic["analytic_jacobian"]
    assert not numerical_diagnostic["analytic_jacobian"]
    check("analytic_and_numerical_newton_endpoints", np.linalg.norm(x1-numerical_x), 2e-12)
    check("analytic_and_numerical_newton_stages",
          np.linalg.norm(diagnostic["F"]-numerical_diagnostic["F"]), 2e-11)

    rejected = False
    try:
        scheme(s=2, q=4, gauge=1e-8*np.eye(4)).step(duffing, x0, 0., h)
    except ValueError as error:
        rejected = "not alternating" in str(error)
    check("nonalternating_gauge_rejected", 0.0 if rejected else 1.0, 0.0)

    # A full stage gauge and state-dependent open PH model exercise covariance
    # beyond constant coefficients and the common-effort subclass.
    rng = np.random.default_rng(2701)
    raw = rng.normal(size=(4, 4))
    k_stage = .3*h**4*(raw-raw.T)
    original = PHModel(
        H=duffing.H, grad=duffing.grad, hess=duffing.hess,
        J=lambda x: (1.0+.2*x[0]**2)*j2,
        R=lambda x: np.diag([.02+.01*x[0]**2, .03+.02*x[1]**2]),
        B=lambda x: np.array([[1.0, .2*x[0]], [.1*x[1], .7]]),
        S=.13*j2,
        u=lambda t: np.array([.2*np.sin(t), -.1*np.cos(2*t)]))
    transform = np.array([[1.7, .3], [-.2, .8]])
    inverse = np.linalg.inv(transform)
    transformed = PHModel(
        H=lambda z: original.H(inverse@z),
        grad=lambda z: inverse.T@original.grad(inverse@z),
        hess=lambda z: inverse.T@original.hess(inverse@z)@inverse,
        J=lambda z: transform@original.J(inverse@z)@transform.T,
        R=lambda z: transform@original.R(inverse@z)@transform.T,
        B=lambda z: transform@original.B(inverse@z),
        S=original.S, u=original.u)
    pullback = np.kron(np.eye(2), inverse)
    new_gauge = pullback.T@k_stage@pullback
    xx, dd = scheme(s=2, q=4, gauge=k_stage).step(original, x0, .3, h)
    zz, dz = scheme(s=2, q=4, gauge=new_gauge).step(transformed, transform@x0, .3, h)
    check("nonorthogonal_covariance_endpoint", np.linalg.norm(zz-transform@xx), 3e-12)
    check("nonorthogonal_covariance_stage_flow", np.linalg.norm(dz["F"]-dd["F"]@transform.T), 3e-11)
    check("nonorthogonal_covariance_stage_effort", np.linalg.norm(dz["e"]-dd["e"]@inverse), 3e-11)
    check("nonorthogonal_covariance_port_output", np.linalg.norm(dz["y"]-dd["y"]), 3e-11)
    check("nonorthogonal_covariance_power",
          max(abs(dz[key]-dd[key]) for key in ("deltaH", "work", "diss")), 3e-12)

    # Compare to the actual unchanged v12 cap gauge, not a second stage solve.
    j4 = np.kron(np.eye(2), j2)
    raw = rng.normal(size=(4, 4))
    q4 = np.eye(4) + .2*raw.T@raw
    raw = rng.normal(size=(4, 4))
    theta = .2*(raw-raw.T)
    generator = j4@q4
    assert np.linalg.norm(theta@generator-generator@theta) > .1
    line = SimpleNamespace(interconnection=j4, hessian=q4)
    quadratic = PHModel(lambda x: .5*x@q4@x,
                        lambda x: q4@x, lambda x: q4, j4)
    cap_results = []
    for stages in (1, 2, 3):
        step = .3
        active = scheme(s=stages, q=stages, gauge=common_gauge(theta))
        stage_map = np.column_stack([
            active.step(quadratic, np.eye(4)[:, column], 0.0, step)[0]
            for column in range(4)])
        s_gauss = gauss_cap_gauge(line, stages, step)
        mixed = s_gauss + step**(2*stages-1)*theta
        common = np.eye(4) + step*j4@mixed
        cap_map = np.linalg.solve(common-.5*step*generator,
                                  common+.5*step*generator)
        error = np.linalg.norm(stage_map-cap_map, ord=np.inf)
        check("v12_mixed_cap_bridge_s{}".format(stages), error, 4e-12)
        cap_results.append({"stages": stages, "step": step, "map_error": float(error)})

    report = {"all_passed": True, "check_count": len(records), "checks": records,
              "cap_bridge": cap_results,
              "noncommuting_theta_generator_norm": float(np.linalg.norm(
                  theta@generator-generator@theta))}
    destination = Path(__file__).resolve().parent / "results" / "ph_gauges_selfcheck.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("PH gauge engine: {} checks passed; results/{}".format(len(records), destination.name))
    for entry in cap_results:
        print("  v12 cap bridge s={}: {:.3e}".format(entry["stages"], entry["map_error"]))
    return report


if __name__ == "__main__":
    selfcheck()

