#!/usr/bin/env python3
"""Executable checks for the v12 space--time transmission-line results.

The spatial data use the modal power chart of the polynomial reduction

    P_N --d/dz--> P_{N-1}

from Moulla--Lefevre--Maschke (2012).  Both spaces are represented in the
shifted, L2-orthonormal Legendre basis.  Thus the rectangular mixed mass
matrix is [I; 0], while endpoint traces remain explicit.

The script verifies twenty-six independent claims and raises AssertionError
when a stated tolerance or convergence law fails.

1. The rectangular spatial Stokes identity is exact to roundoff.
2. The chosen boundary causality produces a finite-dimensional open Dirac
   realization, including its skew feedthrough term.
3. The spatial-derivative cap trace has a numerical right inverse and the
   volume-to-cap pullback recovers an arbitrary alternating cap gauge.
4. The raw trace-generated cylinder system containing
   Omega_st=T.T S_h T realizes exactly the prescribed cap propagator.
5. A cylinder-lifted mixed cap gauge preserves the exact cap--mantle power
   balance, dissipation, and cancellation at a glued interface.
6. The Gauss-s hierarchy with the perturbation h**(2*s-1) Theta has global
   order 2*s for s=1,2,3 while remaining exactly energy preserving.
7. The deformed propagators generically fail to commute with the spatial
   generator, so they are outside its scalar functional calculus.
8. A one-parameter gauge selected at a small reference step reduces modal
   phase error at unseen coarser steps without a material conditioning cost.
9. The temporal gauge is a Cayley chart on the energy-orthogonal group.
10. A P_s -> P_{s-1} Gauss cell with active stage variables condenses to
   the Pade map and induces the announced two-point cap gauge.
11. The rectangular Stokes identity alone reconstructs J, B, B.T, and N.
12. The norm used to scale the gauge is invariant under changes of state
   coordinates compatible with the energy metric.
13. The two polynomial endpoint reconstructions used in the rectangular
   theorem are nonsingular.
14. Cancellation-free formulas evaluate exact-flow and Gauss gauges stably
   near the identity, where the algebraic inverse Cayley chart is sensitive.
15. The explicit coefficient Theta_s^star raises the Gauss-s deformation
   to exact order 2*s+2 and is unique within the fixed one-term ansatz.
16. The coefficient d_s is recovered directly from a high-precision local
   asymptotic test, without estimating it from a global convergence slope.
17. A second, noncommuting term retains order 2*s+2 while moving the map
   outside the scalar functional calculus of the line generator.
18. The second central coefficient Xi_s^star raises the order by two more.
19. Decimal multiprecision resolves the s=3 noncommutation law at step
    sizes for which binary64 roundoff masks the commutator.
20. Every numerical block is protected by an executable assertion.
21. Coordinate poles of C_3 at z=+/-i*sqrt(60) are removable in the
    polynomial propagator Pi(-z)^(-1) Pi(z), whose continued value is one.
22. The corrected (s,r)=(2,1) hierarchy has a left-half-plane pole and is
    therefore not A-stable, despite sixth-order Hamiltonian accuracy.
23. Exact rational arithmetic supports, but does not promote to a theorem,
    the coefficient sign conjecture for 1<=s<=12 and 0<=j<=12.
24. Every tested rank-four alternating cap change is recovered by a
    four-dimensional skew interior Schur generator.
25. The Gauss-3 coordinate has the announced positive one-mode Stieltjes
    realization.
26. The corrected (s,r)=(2,1) coordinate violates positive Stieltjes
    realizability because its shifted quotient is negative for u>60.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, getcontext, localcontext
from fractions import Fraction
from math import factorial

import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import legder, leggauss, legval
from scipy.linalg import block_diag, expm, null_space
from scipy.optimize import minimize_scalar


@dataclass(frozen=True)
class LineData:
    derivative: np.ndarray
    mixed_mass: np.ndarray
    boundary_form: np.ndarray
    left_trace: np.ndarray
    right_trace: np.ndarray
    interconnection: np.ndarray
    input_map: np.ndarray
    output_map: np.ndarray
    feedthrough: np.ndarray
    hessian: np.ndarray


def shifted_legendre_line(
    degree: int = 8,
    length: float = 2.0,
    capacitance: float = 3.0,
    inductance: float = 2.0,
) -> LineData:
    """Build the exact rectangular polynomial transmission-line data."""
    nodes, weights = leggauss(2 * degree + 5)
    physical_weights = 0.5 * length * weights

    values = np.empty((degree + 1, nodes.size))
    derivatives = np.empty_like(values)
    for order in range(degree + 1):
        coefficient = np.zeros(order + 1)
        coefficient[-1] = 1.0
        normalization = np.sqrt((2.0 * order + 1.0) / length)
        values[order] = normalization * legval(nodes, coefficient)
        derivatives[order] = (
            normalization
            * (2.0 / length)
            * legval(nodes, legder(coefficient))
        )

    # D_{kn}=<phi_k, d phi_n/dz>, k=0,...,N-1, n=0,...,N.
    derivative = np.einsum(
        "ki,i,ni->kn",
        values[:degree],
        physical_weights,
        derivatives,
    )
    mixed_mass = np.vstack((np.eye(degree), np.zeros((1, degree))))

    left_trace = np.array(
        [
            np.sqrt((2.0 * order + 1.0) / length) * (-1.0) ** order
            for order in range(degree + 1)
        ]
    )
    right_trace = np.array(
        [
            np.sqrt((2.0 * order + 1.0) / length)
            for order in range(degree + 1)
        ]
    )
    boundary_form = np.outer(right_trace, right_trace) - np.outer(
        left_trace, left_trace
    )

    # Reconstruct an effort polynomial from its N active power coordinates
    # and one prescribed boundary trace.  The q effort is prescribed at the
    # left endpoint and the p effort at the right endpoint, matching the
    # complementary closed conditions used in the 2012 line experiment.
    canonical_lift = mixed_mass.copy()
    radical = np.concatenate((np.zeros(degree), [1.0]))

    def boundary_lift(trace: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        reconstruction = canonical_lift - np.outer(
            radical, trace[:degree] / trace[degree]
        )
        injection = radical / trace[degree]
        return reconstruction, injection

    reconstruct_left, inject_left = boundary_lift(left_trace)
    reconstruct_right, inject_right = boundary_lift(right_trace)

    zero = np.zeros((degree, degree))
    interconnection = np.block(
        [
            [zero, -derivative @ reconstruct_right],
            [-derivative @ reconstruct_left, zero],
        ]
    )

    # Inputs u=(e_q(0), e_p(length)); outputs
    # y=(e_p(0), -e_q(length)).
    input_map = np.zeros((2 * degree, 2))
    input_map[degree:, 0] = -derivative @ inject_left
    input_map[:degree, 1] = -derivative @ inject_right

    output_map = np.zeros((2, 2 * degree))
    output_map[0, degree:] = left_trace @ reconstruct_right
    output_map[1, :degree] = -right_trace @ reconstruct_left
    feedthrough = np.array(
        [
            [0.0, left_trace @ inject_right],
            [-right_trace @ inject_left, 0.0],
        ]
    )

    hessian = np.diag(
        np.concatenate(
            (
                np.full(degree, 1.0 / capacitance),
                np.full(degree, 1.0 / inductance),
            )
        )
    )
    return LineData(
        derivative=derivative,
        mixed_mass=mixed_mass,
        boundary_form=boundary_form,
        left_trace=left_trace,
        right_trace=right_trace,
        interconnection=interconnection,
        input_map=input_map,
        output_map=output_map,
        feedthrough=feedthrough,
        hessian=hessian,
    )


def mixed_generator(data: LineData) -> tuple[np.ndarray, np.ndarray]:
    """Return the Gauss-2 coefficient and an intrinsically scaled two-form.

    If y=Q**(1/2)x are energy-orthonormal coordinates, the matrix of a
    covariant two-form is Q**(-1/2) Theta Q**(-1/2).  We normalize in that
    chart, rather than with a coordinate-dependent Frobenius norm.
    """
    degree = data.derivative.shape[0]
    size = 2 * degree
    raw = np.zeros((size, size))
    couplings = (
        (1, degree + degree - 2, 1.0),
        (2, degree + degree - 3, -0.7),
        (degree - 2, degree + 1, 0.4),
    )
    for first, second, value in couplings:
        raw[first, second] = value
        raw[second, first] = -value

    gauss_coefficient = data.hessian @ data.interconnection @ data.hessian / 12.0
    metric_sqrt = np.diag(np.sqrt(np.diag(data.hessian)))
    metric_inv_sqrt = np.diag(1.0 / np.sqrt(np.diag(data.hessian)))
    gauss_orthonormal = (
        metric_inv_sqrt @ gauss_coefficient @ metric_inv_sqrt
    )
    raw *= np.linalg.norm(gauss_orthonormal, ord="fro") / np.linalg.norm(
        raw, ord="fro"
    )
    mixed = metric_sqrt @ raw @ metric_sqrt
    return gauss_coefficient, mixed


def two_form_norm(form: np.ndarray, metric: np.ndarray) -> float:
    """Metric-induced norm of a covariant two-form."""
    metric_inverse = np.linalg.inv(metric)
    squared = np.trace(
        metric_inverse @ form @ metric_inverse @ form.T
    )
    return float(np.sqrt(max(squared, 0.0)))


def spatial_cap_trace(data: LineData) -> np.ndarray:
    """Return T(z0,z1)=(I_2 tensor D_N)(z1-z0)."""
    spatial_derivative = np.block(
        [
            [data.derivative, np.zeros_like(data.derivative)],
            [np.zeros_like(data.derivative), data.derivative],
        ]
    )
    return np.hstack((-spatial_derivative, spatial_derivative))


def cap_map(
    data: LineData,
    step: float,
    gauge: np.ndarray,
) -> np.ndarray:
    """Closed cap-to-cap map for Delta=h J(Q xbar-S Delta)."""
    identity = np.eye(data.hessian.shape[0])
    generator = data.interconnection @ data.hessian
    common = identity + step * data.interconnection @ gauge
    return np.linalg.solve(
        common - 0.5 * step * generator,
        common + 0.5 * step * generator,
    )


def structural_residuals(data: LineData) -> dict[str, float]:
    stokes = (
        data.mixed_mass @ data.derivative
        + data.derivative.T @ data.mixed_mass.T
        - data.boundary_form
    )
    return {
        "stokes": np.linalg.norm(stokes, ord="fro"),
        "interconnection_skew": np.linalg.norm(
            data.interconnection + data.interconnection.T, ord="fro"
        ),
        "output_input_duality": np.linalg.norm(
            data.output_map - data.input_map.T, ord="fro"
        ),
        "feedthrough_skew": np.linalg.norm(
            data.feedthrough + data.feedthrough.T, ord="fro"
        ),
    }


def volume_cap_reduction_checks(data: LineData) -> dict[str, float]:
    """Check the right-inverse and pullback formulas used in the v8 reduction.

    The raw slab potential has two temporal caps, each carrying q and p
    polynomials in P_N.  The cap trace

        T(z_0,z_1) = (I_2 tensor D_N)(z_1-z_0)

    is onto because D_N:P_N->P_{N-1} is onto.  A Moore--Penrose right
    inverse is used only as a convenient numerical section; the theorem is
    independent of that choice.
    """
    degree = data.derivative.shape[0]
    trace = spatial_cap_trace(data)
    single_cap_derivative = trace[:, trace.shape[1] // 2 :]
    section = single_cap_derivative.T @ np.linalg.inv(
        single_cap_derivative @ single_cap_derivative.T
    )
    lifting = np.vstack((np.zeros_like(section), section))

    _, theta = mixed_generator(data)
    volume_form = trace.T @ theta @ trace
    reduced_form = lifting.T @ volume_form @ lifting

    cap_dimension = 2 * degree
    return {
        "trace_right_inverse": np.linalg.norm(
            trace @ lifting - np.eye(cap_dimension), ord="fro"
        ),
        "volume_form_skew": np.linalg.norm(
            volume_form + volume_form.T, ord="fro"
        ),
        "pullback_recovery": np.linalg.norm(
            reduced_form - theta, ord="fro"
        )
        / np.linalg.norm(theta, ord="fro"),
        "trace_rank": float(np.linalg.matrix_rank(trace)),
        "cap_dimension": float(cap_dimension),
        "boundary_gauge_dimension": float(
            cap_dimension * (cap_dimension - 1) // 2
        ),
    }


def dynamic_condensation_checks(
    data: LineData,
    step: float = 0.1,
) -> dict[str, float]:
    """Assemble the raw cylinder system and compare its Schur map.

    The raw slab variable z has a nonunique component in ker(T).  A basis N
    of that kernel and a multiplier lambda impose N.T z=0 without changing
    the physical cap map.  The effort remains in the rectangular space
    (P_N)^2, with its two boundary constraints imposed explicitly.  Thus
    neither the square interconnection nor the cap gauge is inserted in the
    raw system: its blocks contain D_N, M_N, endpoint traces, and
    Omega_st=T.T S_h T.  Their elimination recovers both objects.
    """
    gauss_coefficient, theta = mixed_generator(data)
    cap_gauge = step * gauss_coefficient + step**3 * theta
    trace = spatial_cap_trace(data)
    volume_form = trace.T @ cap_gauge @ trace
    kernel_basis = null_space(trace)

    cap_dimension, raw_dimension = trace.shape
    kernel_dimension = kernel_basis.shape[1]
    polynomial_dimension = data.derivative.shape[1]
    effort_dimension = 2 * polynomial_dimension
    identity = np.eye(cap_dimension)
    zero_cap_raw = np.zeros((cap_dimension, raw_dimension))
    zero_cap_effort = np.zeros((cap_dimension, effort_dimension))
    zero_cap_kernel = np.zeros((cap_dimension, kernel_dimension))
    zero_kernel_cap = np.zeros((kernel_dimension, cap_dimension))
    zero_kernel_effort = np.zeros((kernel_dimension, effort_dimension))

    effort_projection = block_diag(
        data.mixed_mass.T,
        data.mixed_mass.T,
    )
    zero_derivative = np.zeros_like(data.derivative)
    rectangular_kinematics = np.block(
        [
            [zero_derivative, -data.derivative],
            [-data.derivative, zero_derivative],
        ]
    )
    boundary_constraints = np.zeros((2, effort_dimension))
    boundary_constraints[0, :polynomial_dimension] = data.left_trace
    boundary_constraints[1, polynomial_dimension:] = data.right_trace
    zero_boundary_raw = np.zeros((2, raw_dimension))
    zero_boundary_cap = np.zeros((2, cap_dimension))
    zero_boundary_kernel = np.zeros((2, kernel_dimension))

    # Unknown ordering: raw slab z, full P_N effort, final cap x1, multiplier.
    def assemble(gauge_selector: np.ndarray) -> np.ndarray:
        """Use G.T z=0 with G.T ker(T) invertible to select a representative."""
        return np.block(
            [
                [trace, zero_cap_effort, -identity, zero_cap_kernel],
                [
                    zero_cap_raw,
                    -step * rectangular_kinematics,
                    identity,
                    zero_cap_kernel,
                ],
                [
                    volume_form,
                    trace.T @ effort_projection,
                    -0.5 * trace.T @ data.hessian,
                    kernel_basis,
                ],
                [
                    gauge_selector.T,
                    zero_kernel_effort,
                    zero_kernel_cap,
                    np.zeros((kernel_dimension, kernel_dimension)),
                ],
                [
                    zero_boundary_raw,
                    boundary_constraints,
                    zero_boundary_cap,
                    zero_boundary_kernel,
                ],
            ]
        )

    system = assemble(kernel_basis)
    right_hand_side = np.vstack(
        (
            -identity,
            identity,
            0.5 * trace.T @ data.hessian,
            np.zeros((kernel_dimension, cap_dimension)),
            np.zeros((2, cap_dimension)),
        )
    )
    solution = np.linalg.solve(system, right_hand_side)
    final_slice = slice(
        raw_dimension + effort_dimension,
        raw_dimension + effort_dimension + cap_dimension,
    )
    condensed_map = solution[final_slice]
    cap_reference = cap_map(data, step, cap_gauge)

    # A genuinely different complement of ker(T), not merely a different
    # right inverse used after the solve.  Because T N=0, the selector below
    # still satisfies G.T N=I and therefore fixes one representative.
    rng = np.random.default_rng(20260818)
    alternative_selector = kernel_basis + trace.T @ (
        0.05 * rng.normal(size=(cap_dimension, kernel_dimension))
    )
    alternative_system = assemble(alternative_selector)
    alternative_solution = np.linalg.solve(
        alternative_system, right_hand_side
    )
    alternative_condensed_map = alternative_solution[final_slice]

    canonical_section = trace.T @ np.linalg.inv(trace @ trace.T)
    alternative_section = canonical_section + (
        np.eye(raw_dimension) - canonical_section @ trace
    ) @ (0.1 * rng.normal(size=(raw_dimension, cap_dimension)))
    canonical_pullback = canonical_section.T @ volume_form @ canonical_section
    alternative_pullback = (
        alternative_section.T @ volume_form @ alternative_section
    )
    effort_reconstruction = np.linalg.solve(
        np.vstack((effort_projection, boundary_constraints)),
        np.vstack((identity, np.zeros((2, cap_dimension)))),
    )
    recovered_interconnection = rectangular_kinematics @ effort_reconstruction

    return {
        "raw_system_residual": np.linalg.norm(
            system @ solution - right_hand_side, ord="fro"
        ),
        "condensed_map_relative_error": np.linalg.norm(
            condensed_map - cap_reference, ord="fro"
        )
        / np.linalg.norm(cap_reference, ord="fro"),
        "section_independence": np.linalg.norm(
            canonical_pullback - alternative_pullback, ord="fro"
        )
        / np.linalg.norm(cap_gauge, ord="fro"),
        "gauge_fixing_map_change": np.linalg.norm(
            condensed_map - alternative_condensed_map, ord="fro"
        )
        / np.linalg.norm(cap_reference, ord="fro"),
        "alternative_system_residual": np.linalg.norm(
            alternative_system @ alternative_solution - right_hand_side,
            ord="fro",
        ),
        "gauge_multiplier_norm": np.linalg.norm(
            solution[-kernel_dimension:], ord="fro"
        ),
        "rectangular_generator_recovery": np.linalg.norm(
            recovered_interconnection - data.interconnection, ord="fro"
        ),
        "raw_system_condition": np.linalg.cond(system),
    }


def rectangular_reconstruction_checks(data: LineData) -> dict[str, float]:
    """Recover the open square Dirac blocks from rectangular Stokes data."""
    degree = data.derivative.shape[0]
    polynomial_dimension = degree + 1
    active_dimension = 2 * degree
    effort_dimension = 2 * polynomial_dimension

    effort_projection = block_diag(data.mixed_mass.T, data.mixed_mass.T)
    zero_derivative = np.zeros_like(data.derivative)
    kinematics = np.block(
        [
            [zero_derivative, -data.derivative],
            [-data.derivative, zero_derivative],
        ]
    )
    input_trace = np.zeros((2, effort_dimension))
    input_trace[0, :polynomial_dimension] = data.left_trace
    input_trace[1, polynomial_dimension:] = data.right_trace
    output_trace = np.zeros((2, effort_dimension))
    output_trace[0, polynomial_dimension:] = data.left_trace
    output_trace[1, :polynomial_dimension] = -data.right_trace

    reconstruction = np.vstack((effort_projection, input_trace))
    reconstruct_left = np.vstack((data.mixed_mass.T, data.left_trace))
    reconstruct_right = np.vstack((data.mixed_mass.T, data.right_trace))
    inverse = np.linalg.inv(reconstruction)
    active_injection = inverse[:, :active_dimension]
    input_injection = inverse[:, active_dimension:]
    recovered_j = kinematics @ active_injection
    recovered_b = kinematics @ input_injection
    recovered_c = output_trace @ active_injection
    recovered_n = output_trace @ input_injection
    stokes = (
        effort_projection.T @ kinematics
        + kinematics.T @ effort_projection
        - input_trace.T @ output_trace
        - output_trace.T @ input_trace
    )
    return {
        "rectangular_stokes_block": np.linalg.norm(stokes, ord="fro"),
        "recovered_interconnection": np.linalg.norm(
            recovered_j - data.interconnection, ord="fro"
        ),
        "recovered_input_map": np.linalg.norm(
            recovered_b - data.input_map, ord="fro"
        ),
        "recovered_output_map": np.linalg.norm(
            recovered_c - data.output_map, ord="fro"
        ),
        "recovered_feedthrough": np.linalg.norm(
            recovered_n - data.feedthrough, ord="fro"
        ),
        "recovered_j_skew": np.linalg.norm(
            recovered_j + recovered_j.T, ord="fro"
        ),
        "recovered_port_duality": np.linalg.norm(
            recovered_c - recovered_b.T, ord="fro"
        ),
        "recovered_n_skew": np.linalg.norm(
            recovered_n + recovered_n.T, ord="fro"
        ),
        "left_endpoint_min_singular": np.linalg.svd(
            reconstruct_left, compute_uv=False
        )[-1],
        "right_endpoint_min_singular": np.linalg.svd(
            reconstruct_right, compute_uv=False
        )[-1],
        "reconstruction_condition": np.linalg.cond(reconstruction),
    }


def pade_coefficients(order: int) -> list[float]:
    """Coefficients of the numerator P_s of the diagonal Pade pair."""
    return [
        factorial(2 * order - degree)
        * factorial(order)
        / (
            factorial(2 * order)
            * factorial(degree)
            * factorial(order - degree)
        )
        for degree in range(order + 1)
    ]


def _bernoulli_number(index: int) -> Fraction:
    """Exact Bernoulli number; the B_1 convention is irrelevant here."""
    values = [Fraction(0) for _ in range(index + 1)]
    for degree in range(index + 1):
        values[degree] = Fraction(1, degree + 1)
        for position in range(degree, 0, -1):
            values[position - 1] = position * (
                values[position - 1] - values[position]
            )
    return values[0]


def exact_central_coefficient(stages: int, level: int) -> Fraction:
    """Return alpha_{s,j} exactly by formal power-series division in u=z^2."""
    target_degree = stages + level
    coefficients = [
        Fraction(
            factorial(2 * stages - degree) * factorial(stages),
            factorial(2 * stages)
            * factorial(degree)
            * factorial(stages - degree),
        )
        for degree in range(stages + 1)
    ]
    even = [coefficients[2 * degree] for degree in range(stages // 2 + 1)]
    odd_hat = [
        coefficients[2 * degree + 1]
        for degree in range((stages - 1) // 2 + 1)
    ]
    quotient: list[Fraction] = []
    for degree in range(target_degree + 1):
        residual = even[degree] if degree < len(even) else Fraction(0)
        residual -= sum(
            (
                odd_hat[index] * quotient[degree - index]
                for index in range(1, min(degree, len(odd_hat) - 1) + 1)
            ),
            Fraction(0),
        )
        quotient.append(residual / odd_hat[0])
    gauss_coordinate = Fraction(1, 2) * quotient[target_degree]
    exact_coordinate = _bernoulli_number(2 * target_degree) / factorial(
        2 * target_degree
    )
    return exact_coordinate - gauss_coordinate


def coefficient_sign_checks(
    max_stages: int = 12,
    max_level: int = 12,
) -> dict[str, float]:
    """Exact evidence for, not a proof of, the alpha sign conjecture."""
    coefficients = [
        (stages, level, exact_central_coefficient(stages, level))
        for stages in range(1, max_stages + 1)
        for level in range(max_level + 1)
    ]
    sign_failures = sum(
        coefficient == 0
        or ((-1) ** (stages + level + 1)) * coefficient <= 0
        for stages, level, coefficient in coefficients
    )
    return {
        "tested_coefficients": float(len(coefficients)),
        "sign_failures": float(sign_failures),
        "alpha_3_2_relative_error": abs(
            float(exact_central_coefficient(3, 2) - Fraction(173, 11975040000))
        ),
        "alpha_4_3_relative_error": abs(
            float(exact_central_coefficient(4, 3) - Fraction(1483, 538238769868800))
        ),
    }


def rational_reduction_checks() -> dict[str, float]:
    """Check a removable chart pole and a genuine left-half-plane pole."""
    coordinate_resonance = 1j * np.sqrt(60.0)
    p3 = lambda z: 1.0 + z / 2.0 + z**2 / 10.0 + z**3 / 120.0
    odd_hat3 = lambda z: 0.5 + z**2 / 120.0
    delta31 = lambda z: z**6 / 100800.0
    pi31 = lambda z: p3(z) + 2.0 * odd_hat3(z) * delta31(z)
    continued_value = pi31(coordinate_resonance) / pi31(-coordinate_resonance)

    # Pi_{2,1}(z)=1+z/2+z^2/12-z^4/720.
    roots = np.roots([-1.0 / 720.0, 0.0, 1.0 / 12.0, 0.5, 1.0])
    positive_root = float(
        min(
            root.real
            for root in roots
            if abs(root.imag) < 1.0e-10 and root.real > 0.0
        )
    )
    true_left_pole = -positive_root
    pi21 = lambda z: 1.0 + z / 2.0 + z**2 / 12.0 - z**4 / 720.0
    imaginary_axis = np.linspace(-100.0, 100.0, 4001)
    values = np.array(
        [pi21(1j * xi) / pi21(-1j * xi) for xi in imaginary_axis]
    )
    return {
        "c3_coordinate_denominator": abs(odd_hat3(coordinate_resonance)),
        "pi31_at_coordinate_pole_error": abs(
            pi31(coordinate_resonance) + 5.0
        ),
        "continued_map_at_coordinate_pole_error": abs(continued_value - 1.0),
        "left_half_plane_pole": true_left_pole,
        "left_pole_denominator_residual": abs(pi21(-true_left_pole)),
        "left_pole_numerator_magnitude": abs(pi21(true_left_pole)),
        "imaginary_axis_unit_modulus_defect": float(
            np.max(np.abs(np.abs(values) - 1.0))
        ),
    }


def interior_schur_generation_checks() -> dict[str, float]:
    """Verify a rank-minimal alternating auxiliary Schur generator."""
    rng = np.random.default_rng(20260818)
    cap_dimension = 6
    interior_dimension = 4
    coupling = rng.standard_normal((cap_dimension, interior_dimension))
    while np.linalg.matrix_rank(coupling) < interior_dimension:
        coupling = rng.standard_normal((cap_dimension, interior_dimension))

    j2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
    canonical = block_diag(j2, j2)
    interior_form = -canonical
    delta = coupling @ canonical @ coupling.T

    seed = rng.standard_normal((cap_dimension, cap_dimension))
    theta_zero = seed - seed.T
    theta_one = theta_zero + delta
    full_generator = np.block(
        [
            [theta_zero, coupling],
            [-coupling.T, interior_form],
        ]
    )
    effective = theta_zero + coupling @ np.linalg.solve(
        interior_form, coupling.T
    )

    cap_state = rng.standard_normal(cap_dimension)
    interior_state = np.linalg.solve(
        interior_form, coupling.T @ cap_state
    )
    generated_response = theta_zero @ cap_state + coupling @ interior_state
    return {
        "full_generator_skew_defect": float(
            np.linalg.norm(full_generator + full_generator.T)
        ),
        "effective_gauge_error": float(np.linalg.norm(effective - theta_one)),
        "eliminated_response_error": float(
            np.linalg.norm(generated_response - theta_one @ cap_state)
        ),
        "target_change_rank": float(np.linalg.matrix_rank(delta)),
        "interior_dimension": float(interior_dimension),
    }


def stieltjes_interior_generation_checks() -> dict[str, float]:
    """Check the positive Gauss-3 realization and the (2,1) obstruction."""
    grid = np.linspace(0.0, 120.0, 481)
    gauss_three_f = 5.0 / (grid + 60.0)
    weight = 1.0 / 12.0
    spectral_point = 1.0 / 60.0
    interior_realization = weight / (1.0 + grid * spectral_point)
    gauss_three_c = 6.0 * (grid + 10.0) / (grid + 60.0)
    generated_coordinate = 1.0 + grid * interior_realization

    obstruction_point = 61.0
    corrected_shifted_value = (
        1.0 / 12.0 - obstruction_point / 720.0
    )
    return {
        "gauss3_stieltjes_realization_error": float(
            np.max(np.abs(gauss_three_f - interior_realization))
        ),
        "gauss3_coordinate_generation_error": float(
            np.max(np.abs(gauss_three_c - generated_coordinate))
        ),
        "gauss3_weight": weight,
        "gauss3_spectral_point": spectral_point,
        "corrected_21_shifted_value_at_u61": corrected_shifted_value,
        "corrected_21_zero_at_u60": 1.0 / 12.0 - 60.0 / 720.0,
    }


def pade_polynomial(order: int, argument: np.ndarray) -> np.ndarray:
    """Numerator P_s of the diagonal [s/s] Pade approximant of exp."""
    identity = np.eye(argument.shape[0])
    result = np.zeros_like(argument)
    power = identity.copy()
    for degree, coefficient in enumerate(pade_coefficients(order)):
        result += coefficient * power
        power = power @ argument
    return result


def gauss_pade_map(data: LineData, order: int, step: float) -> np.ndarray:
    """The diagonal Pade [s/s] map, i.e. the s-stage Gauss map."""
    generator = data.interconnection @ data.hessian
    return np.linalg.solve(
        pade_polynomial(order, -step * generator),
        pade_polynomial(order, step * generator),
    )


def gauss_cap_gauge(data: LineData, order: int, step: float) -> np.ndarray:
    """Cancellation-free two-point gauge whose map is exactly Gauss-s.

    If P_s=E_s+O_s and O_s(z)=z O_hat_s(z), then the common Cayley block is
    (1/2) E_s O_hat_s^{-1}.  We form

        common-I = ((1/2) E_s-O_hat_s) O_hat_s^{-1}

    directly.  Its constant coefficient is exactly zero, avoiding both an
    identity subtraction and inversion of R_s-I.
    """
    generator = data.interconnection @ data.hessian
    symplectic_form = np.linalg.inv(data.interconnection)
    identity = np.eye(generator.shape[0])
    squared_argument = (step * generator) @ (step * generator)
    coefficients = pade_coefficients(order)
    odd_reduction = np.zeros_like(generator)
    delta = np.zeros_like(generator)
    even_power = identity.copy()
    for index in range((order + 1) // 2):
        odd_degree = 2 * index + 1
        odd_reduction += coefficients[odd_degree] * even_power
        even_coefficient = 0.5 * coefficients[2 * index]
        delta_coefficient = even_coefficient - coefficients[odd_degree]
        if index:
            delta += delta_coefficient * even_power
        even_power = even_power @ squared_argument
    if order % 2 == 0:
        delta += 0.5 * coefficients[order] * even_power
    common_minus_identity = np.linalg.solve(
        odd_reduction.T, delta.T
    ).T
    gauge = symplectic_form @ common_minus_identity / step
    return 0.5 * (gauge - gauge.T)


def gauge_from_energy_map(
    data: LineData,
    step: float,
    energy_map: np.ndarray,
) -> np.ndarray:
    """Algebraic inverse chart from a Q-orthogonal map to a cap two-form.

    This exact formula is intentionally retained for moderate-step roundtrip
    tests.  It is ill-conditioned near M=I and is not used for small-step
    evaluation of analytically known exact-flow or Gauss gauges.
    """
    generator = data.interconnection @ data.hessian
    symplectic_form = np.linalg.inv(data.interconnection)
    identity = np.eye(generator.shape[0])
    cayley_ratio = np.linalg.solve(
        (energy_map - identity).T,
        (energy_map + identity).T,
    ).T
    common = 0.5 * step * generator @ cayley_ratio
    return symplectic_form @ (common - identity) / step


def exact_flow_gauge_series(data: LineData, step: float) -> np.ndarray:
    """Stable local Bernoulli-series evaluation of the exact-flow gauge.

    The scalar series has radius 2*pi, set by the nearest nonzero pole of
    (z/2)coth(z/2).  This routine is used only in the small-step experiment;
    a general-step implementation should use a Schur-based matrix function
    and detect exact-flow resonances exp(h*lambda)=1.
    """
    generator = data.interconnection @ data.hessian
    symplectic_form = np.linalg.inv(data.interconnection)
    argument_squared = (step * generator) @ (step * generator)
    # B_{2k}/(2k)! through degree twelve for psi(z)=(z/2)coth(z/2).
    coefficients = (
        1.0 / 12.0,
        -1.0 / 720.0,
        1.0 / 30240.0,
        -1.0 / 1209600.0,
        1.0 / 47900160.0,
        -691.0 / 1307674368000.0,
    )
    bracket = np.zeros_like(generator)
    power = argument_squared.copy()
    for coefficient in coefficients:
        bracket += coefficient * power
        power = power @ argument_squared
    gauge = symplectic_form @ bracket / step
    return 0.5 * (gauge - gauge.T)


def cayley_chart_checks(
    data: LineData,
    step: float = 0.1,
) -> dict[str, float]:
    """Check the gauge <-> O(Q) inverse formulas, including the exact flow."""
    _, theta = mixed_generator(data)
    gauge = 0.07 * theta
    energy_map = cap_map(data, step, gauge)
    recovered_gauge = gauge_from_energy_map(data, step, energy_map)
    recovered_map = cap_map(data, step, recovered_gauge)
    generator = data.interconnection @ data.hessian
    exact_map = expm(step * generator)
    exact_gauge = gauge_from_energy_map(data, step, exact_map)
    exact_reconstruction = cap_map(data, step, exact_gauge)
    return {
        "gauge_roundtrip_relative_error": np.linalg.norm(
            recovered_gauge - gauge, ord="fro"
        )
        / np.linalg.norm(gauge, ord="fro"),
        "map_roundtrip_relative_error": np.linalg.norm(
            recovered_map - energy_map, ord="fro"
        )
        / np.linalg.norm(energy_map, ord="fro"),
        "recovered_gauge_skew_defect": np.linalg.norm(
            recovered_gauge + recovered_gauge.T, ord="fro"
        ),
        "energy_orthogonality_defect": np.linalg.norm(
            energy_map.T @ data.hessian @ energy_map - data.hessian,
            ord="fro",
        )
        / np.linalg.norm(data.hessian, ord="fro"),
        "exact_flow_reconstruction": np.linalg.norm(
            exact_reconstruction - exact_map, ord="fro"
        )
        / np.linalg.norm(exact_map, ord="fro"),
        "exact_gauge_skew_defect": np.linalg.norm(
            exact_gauge + exact_gauge.T, ord="fro"
        ),
    }


def cayley_conditioning_checks(
    data: LineData,
    steps: tuple[float, ...] = (1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8),
) -> list[dict[str, float]]:
    """Separate stable analytic evaluation from inverse-chart sensitivity."""
    generator = data.interconnection @ data.hessian
    rows: list[dict[str, float]] = []
    for step in steps:
        exact_map = expm(step * generator)
        stable_gauge = exact_flow_gauge_series(data, step)
        naive_gauge = gauge_from_energy_map(data, step, exact_map)
        stable_reconstruction = cap_map(data, step, stable_gauge)
        naive_relative_error = np.linalg.norm(
            naive_gauge - stable_gauge, ord="fro"
        ) / np.linalg.norm(stable_gauge, ord="fro")
        rows.append(
            {
                "step": step,
                "naive_gauge_relative_error": naive_relative_error,
                "epsilon_scaled_inverse_error": (
                    naive_relative_error
                    * step**3
                    / np.finfo(float).eps
                ),
                "stable_map_relative_error": np.linalg.norm(
                    stable_reconstruction - exact_map, ord="fro"
                )
                / np.linalg.norm(exact_map, ord="fro"),
                "stable_gauge_skew_defect": np.linalg.norm(
                    stable_gauge + stable_gauge.T, ord="fro"
                ),
            }
        )
    return rows


def order_raising_coefficient(
    data: LineData,
    stages: int,
) -> np.ndarray:
    """Theta_s^star=c_s Omega A^(2s), the unique leading-order selector."""
    generator = data.interconnection @ data.hessian
    symplectic_form = np.linalg.inv(data.interconnection)
    coefficient = (
        (-1.0) ** (stages + 1)
        * factorial(stages) ** 2
        / (factorial(2 * stages) * factorial(2 * stages + 1))
    )
    theta = coefficient * symplectic_form @ np.linalg.matrix_power(
        generator, 2 * stages
    )
    return 0.5 * (theta - theta.T)


def second_cayley_coefficient(stages: int) -> float:
    """The exact second coefficient d_s of the Cayley-gauge difference."""
    if stages == 1:
        return -1.0 / 720.0
    leading = (
        (-1.0) ** (stages + 1)
        * factorial(stages) ** 2
        / (factorial(2 * stages) * factorial(2 * stages + 1))
    )
    return -leading * (
        (stages + 1) * (2 * stages - 3)
        / (6 * (2 * stages - 1) * (2 * stages + 3))
    )


def second_order_raising_coefficient(
    data: LineData,
    stages: int,
) -> np.ndarray:
    """Xi_s^star=d_s Omega A^(2s+2), the second central corrector."""
    generator = data.interconnection @ data.hessian
    symplectic_form = np.linalg.inv(data.interconnection)
    xi = second_cayley_coefficient(stages) * symplectic_form @ (
        np.linalg.matrix_power(generator, 2 * stages + 2)
    )
    return 0.5 * (xi - xi.T)


def _decimal_sine_cosine(argument: Decimal) -> tuple[Decimal, Decimal]:
    """Sine and cosine by convergent Decimal Taylor series."""
    sine = argument
    cosine = Decimal(1)
    sine_term = argument
    cosine_term = Decimal(1)
    squared = argument * argument
    tolerance = Decimal(10) ** (-(getcontext().prec - 10))
    for index in range(1, 200):
        cosine_term *= -squared / Decimal((2 * index - 1) * (2 * index))
        cosine += cosine_term
        sine_term *= -squared / Decimal((2 * index) * (2 * index + 1))
        sine += sine_term
        if abs(sine_term) < tolerance and abs(cosine_term) < tolerance:
            break
    return sine, cosine


def asymptotic_coefficient_checks(
    stage_orders: tuple[int, ...] = (1, 2, 3),
    steps: tuple[str, ...] = ("0.1", "0.05", "0.025"),
) -> list[dict[str, float]]:
    """Recover d_s directly on a high-precision harmonic oscillator.

    For A=J_2, every even common block is scalar.  Decimal arithmetic lets
    us divide the local map error by h**(2*s+3) without roundoff masking the
    coefficient for s=2 and s=3.
    """
    rows: list[dict[str, float]] = []
    with localcontext() as context:
        context.prec = 90
        for stages in stage_orders:
            coefficient_c = Decimal((-1) ** (stages + 1)) * Decimal(
                factorial(stages) ** 2
            ) / Decimal(factorial(2 * stages) * factorial(2 * stages + 1))
            if stages == 1:
                coefficient_d = -Decimal(1) / Decimal(720)
            else:
                coefficient_d = -coefficient_c * Decimal(
                    (stages + 1) * (2 * stages - 3)
                ) / Decimal(6 * (2 * stages - 1) * (2 * stages + 3))
            predicted_imaginary = coefficient_d * Decimal(
                (-1) ** (stages + 1)
            )
            coefficients = [
                Decimal(
                    factorial(2 * stages - degree) * factorial(stages)
                )
                / Decimal(
                    factorial(2 * stages)
                    * factorial(degree)
                    * factorial(stages - degree)
                )
                for degree in range(stages + 1)
            ]
            for step_text in steps:
                step = Decimal(step_text)
                even = Decimal(0)
                odd_reduction = Decimal(0)
                for degree, coefficient in enumerate(coefficients):
                    if degree % 2 == 0:
                        even += coefficient * Decimal(
                            (-1) ** (degree // 2)
                        ) * step**degree
                    else:
                        odd_reduction += coefficient * Decimal(
                            (-1) ** ((degree - 1) // 2)
                        ) * step ** (degree - 1)
                common = even / (Decimal(2) * odd_reduction)
                common += coefficient_c * Decimal((-1) ** stages) * (
                    step ** (2 * stages)
                )
                denominator = common * common + step * step / Decimal(4)
                map_real = (
                    common * common - step * step / Decimal(4)
                ) / denominator
                map_imaginary = common * step / denominator
                exact_sine, exact_cosine = _decimal_sine_cosine(step)
                scale = step ** (2 * stages + 3)
                estimated_real = (map_real - exact_cosine) / scale
                estimated_imaginary = (map_imaginary - exact_sine) / scale
                coefficient_error = (
                    (estimated_real * estimated_real
                     + (estimated_imaginary - predicted_imaginary) ** 2)
                    .sqrt()
                    / abs(predicted_imaginary)
                )
                component_error = abs(
                    estimated_imaginary - predicted_imaginary
                ) / abs(predicted_imaginary)
                rows.append(
                    {
                        "stages": float(stages),
                        "step": float(step),
                        "predicted_coefficient": float(predicted_imaginary),
                        "estimated_coefficient": float(estimated_imaginary),
                        "relative_coefficient_error": float(coefficient_error),
                        "relative_component_error": float(component_error),
                    }
                )
    return rows


DecimalMatrix = list[list[Decimal]]


def _dm_identity(size: int) -> DecimalMatrix:
    return [
        [Decimal(int(row == column)) for column in range(size)]
        for row in range(size)
    ]


def _dm_add(left: DecimalMatrix, right: DecimalMatrix) -> DecimalMatrix:
    return [
        [left[i][j] + right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def _dm_sub(left: DecimalMatrix, right: DecimalMatrix) -> DecimalMatrix:
    return [
        [left[i][j] - right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def _dm_scale(scale: Decimal, matrix: DecimalMatrix) -> DecimalMatrix:
    return [[scale * entry for entry in row] for row in matrix]


def _dm_multiply(left: DecimalMatrix, right: DecimalMatrix) -> DecimalMatrix:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                Decimal(0),
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _dm_power(matrix: DecimalMatrix, exponent: int) -> DecimalMatrix:
    result = _dm_identity(len(matrix))
    factor = matrix
    power = exponent
    while power:
        if power % 2:
            result = _dm_multiply(result, factor)
        factor = _dm_multiply(factor, factor)
        power //= 2
    return result


def _dm_solve(left: DecimalMatrix, right: DecimalMatrix) -> DecimalMatrix:
    """Dense Gaussian elimination with Decimal partial pivoting."""
    size = len(left)
    columns = len(right[0])
    augmented = [left[i][:] + right[i][:] for i in range(size)]
    for pivot_column in range(size):
        pivot_row = max(
            range(pivot_column, size),
            key=lambda row: abs(augmented[row][pivot_column]),
        )
        if augmented[pivot_row][pivot_column] == 0:
            raise ArithmeticError("singular Decimal matrix")
        augmented[pivot_column], augmented[pivot_row] = (
            augmented[pivot_row],
            augmented[pivot_column],
        )
        pivot = augmented[pivot_column][pivot_column]
        augmented[pivot_column] = [entry / pivot for entry in augmented[pivot_column]]
        for row in range(size):
            if row == pivot_column:
                continue
            multiplier = augmented[row][pivot_column]
            if multiplier:
                augmented[row] = [
                    augmented[row][column]
                    - multiplier * augmented[pivot_column][column]
                    for column in range(size + columns)
                ]
    return [row[size:] for row in augmented]


def _dm_frobenius(matrix: DecimalMatrix) -> Decimal:
    return sum(
        (entry * entry for row in matrix for entry in row), Decimal(0)
    ).sqrt()


def decimal_noncommutation_checks(
    steps: tuple[str, ...] = ("0.05", "0.025", "0.0125", "0.00625"),
) -> list[dict[str, float]]:
    """Resolve the s=3 two-term commutator beyond binary64 precision.

    The corrected common block is C_3(z)+c_3 z^6+h^8 J Xi.  The theorem
    predicts [M_h,A]=-h^9[J Xi,A]A+O(h^10).  Decimal arithmetic checks the
    exponent and the full leading matrix, not merely a nonzero printout.
    """
    rows: list[dict[str, float]] = []
    with localcontext() as context:
        context.prec = 90
        zero, one = Decimal(0), Decimal(1)
        q = [one, Decimal(2), Decimal(3), Decimal(5)]
        symplectic = [
            [zero, one, zero, zero],
            [-one, zero, zero, zero],
            [zero, zero, zero, one],
            [zero, zero, -one, zero],
        ]
        hessian = [
            [q[i] if i == j else zero for j in range(4)]
            for i in range(4)
        ]
        generator = _dm_multiply(symplectic, hessian)
        xi = [[zero for _ in range(4)] for _ in range(4)]
        xi[0][2], xi[2][0] = one, -one
        k_matrix = _dm_multiply(symplectic, xi)
        k_commutator = _dm_sub(
            _dm_multiply(k_matrix, generator),
            _dm_multiply(generator, k_matrix),
        )
        predicted = _dm_scale(
            -one, _dm_multiply(k_commutator, generator)
        )
        predicted_norm = _dm_frobenius(predicted)
        previous_norm: Decimal | None = None
        for step_text in steps:
            step = Decimal(step_text)
            z = _dm_scale(step, generator)
            z2 = _dm_multiply(z, z)
            identity = _dm_identity(4)
            even = _dm_add(identity, _dm_scale(Decimal(1) / Decimal(10), z2))
            odd_hat = _dm_add(
                _dm_scale(Decimal(1) / Decimal(2), identity),
                _dm_scale(Decimal(1) / Decimal(120), z2),
            )
            common = _dm_scale(Decimal(1) / Decimal(2), _dm_solve(odd_hat, even))
            common = _dm_add(
                common,
                _dm_scale(Decimal(1) / Decimal(100800), _dm_power(z, 6)),
            )
            common = _dm_add(common, _dm_scale(step**8, k_matrix))
            half_z = _dm_scale(Decimal(1) / Decimal(2), z)
            one_step = _dm_solve(
                _dm_sub(common, half_z), _dm_add(common, half_z)
            )
            commutator = _dm_sub(
                _dm_multiply(one_step, generator),
                _dm_multiply(generator, one_step),
            )
            commutator_norm = _dm_frobenius(commutator)
            scaled = _dm_scale(one / step**9, commutator)
            leading_relative_error = _dm_frobenius(
                _dm_sub(scaled, predicted)
            ) / predicted_norm
            observed_exponent = (
                float("nan")
                if previous_norm is None
                else float(
                    (previous_norm / commutator_norm).ln()
                    / Decimal(2).ln()
                )
            )
            rows.append(
                {
                    "step": float(step),
                    "commutator_norm": float(commutator_norm),
                    "observed_exponent": observed_exponent,
                    "leading_relative_error": float(leading_relative_error),
                }
            )
            previous_norm = commutator_norm
    return rows


def order_raising_checks(
    data: LineData,
    final_time: float = 0.5,
    stage_orders: tuple[int, ...] = (1, 2, 3),
    step_counts: tuple[int, ...] = (4, 8, 16, 32, 64),
) -> tuple[list[dict[str, float]], float]:
    """Verify orders 2s+2 or higher for the unique explicit selector."""
    generator = data.interconnection @ data.hessian
    exact_map = expm(final_time * generator)
    rows: list[dict[str, float]] = []
    for stages in stage_orders:
        theta = order_raising_coefficient(data, stages)
        previous_error: float | None = None
        for step_count in step_counts:
            step = final_time / step_count
            gauge = (
                gauss_cap_gauge(data, stages, step)
                + step ** (2 * stages - 1) * theta
            )
            one_step = cap_map(data, step, gauge)
            numerical_map = np.linalg.matrix_power(one_step, step_count)
            error = np.linalg.norm(
                numerical_map - exact_map, ord="fro"
            ) / np.linalg.norm(exact_map, ord="fro")
            rows.append(
                {
                    "stages": float(stages),
                    "steps": float(step_count),
                    "global_error": error,
                    "observed_global_order": (
                        np.nan
                        if previous_error is None
                        else np.log(previous_error / error) / np.log(2.0)
                    ),
                    "energy_defect": np.linalg.norm(
                        one_step.T @ data.hessian @ one_step - data.hessian,
                        ord="fro",
                    )
                    / np.linalg.norm(data.hessian, ord="fro"),
                }
            )
            previous_error = error

    identity_step = 0.1
    raised_one_stage = cap_map(
        data,
        identity_step,
        identity_step * order_raising_coefficient(data, 1),
    )
    gauss_two = gauss_pade_map(data, 2, identity_step)
    gauss_two_identity = np.linalg.norm(
        raised_one_stage - gauss_two, ord="fro"
    ) / np.linalg.norm(gauss_two, ord="fro")
    return rows, gauss_two_identity


def two_term_hierarchy_checks(
    data: LineData,
    final_time: float = 0.5,
    stage_orders: tuple[int, ...] = (1, 2, 3),
    step_counts: tuple[int, ...] = (4, 8, 16, 32, 64),
) -> list[dict[str, float]]:
    """Verify the high-order, noncommuting two-term line hierarchy."""
    _, xi = mixed_generator(data)
    generator = data.interconnection @ data.hessian
    exact_map = expm(final_time * generator)
    rows: list[dict[str, float]] = []
    for stages in stage_orders:
        theta_star = order_raising_coefficient(data, stages)
        previous_error: float | None = None
        for step_count in step_counts:
            step = final_time / step_count
            gauge = (
                gauss_cap_gauge(data, stages, step)
                + step ** (2 * stages - 1) * theta_star
                + step ** (2 * stages + 1) * xi
            )
            one_step = cap_map(data, step, gauge)
            numerical_map = np.linalg.matrix_power(one_step, step_count)
            error = np.linalg.norm(
                numerical_map - exact_map, ord="fro"
            ) / np.linalg.norm(exact_map, ord="fro")
            rows.append(
                {
                    "stages": float(stages),
                    "steps": float(step_count),
                    "global_error": error,
                    "observed_global_order": (
                        np.nan
                        if previous_error is None
                        else np.log(previous_error / error) / np.log(2.0)
                    ),
                    "energy_defect": np.linalg.norm(
                        one_step.T @ data.hessian @ one_step - data.hessian,
                        ord="fro",
                    ) / np.linalg.norm(data.hessian, ord="fro"),
                    "commutator": np.linalg.norm(
                        one_step @ generator - generator @ one_step,
                        ord="fro",
                    ) / (
                        np.linalg.norm(one_step, ord="fro")
                        * np.linalg.norm(generator, ord="fro")
                    ),
                }
            )
            previous_error = error
    return rows


def second_central_hierarchy_checks(
    data: LineData,
    final_time: float = 0.5,
    stage_orders: tuple[int, ...] = (1, 2),
    step_counts: tuple[int, ...] = (4, 8, 16, 32, 64),
) -> list[dict[str, float]]:
    """Verify the additional two-order gain from Xi_s^star."""
    generator = data.interconnection @ data.hessian
    exact_map = expm(final_time * generator)
    rows: list[dict[str, float]] = []
    for stages in stage_orders:
        theta_star = order_raising_coefficient(data, stages)
        xi_star = second_order_raising_coefficient(data, stages)
        previous_error: float | None = None
        for step_count in step_counts:
            step = final_time / step_count
            gauge = (
                gauss_cap_gauge(data, stages, step)
                + step ** (2 * stages - 1) * theta_star
                + step ** (2 * stages + 1) * xi_star
            )
            one_step = cap_map(data, step, gauge)
            numerical_map = np.linalg.matrix_power(one_step, step_count)
            error = np.linalg.norm(
                numerical_map - exact_map, ord="fro"
            ) / np.linalg.norm(exact_map, ord="fro")
            rows.append(
                {
                    "stages": float(stages),
                    "steps": float(step_count),
                    "global_error": error,
                    "observed_global_order": (
                        np.nan
                        if previous_error is None
                        else np.log(previous_error / error) / np.log(2.0)
                    ),
                    "energy_defect": np.linalg.norm(
                        one_step.T @ data.hessian @ one_step - data.hessian,
                        ord="fro",
                    ) / np.linalg.norm(data.hessian, ord="fro"),
                }
            )
            previous_error = error
    return rows


def gauss_collocation_coefficients(
    stages: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gauss nodes, integrated Lagrange matrix, and weights on [0,1]."""
    raw_nodes, raw_weights = leggauss(stages)
    nodes = 0.5 * (raw_nodes + 1.0)
    weights = 0.5 * raw_weights
    coefficients = np.empty((stages, stages))
    integrated_weights = np.empty(stages)
    for column in range(stages):
        basis = Polynomial([1.0])
        for other in range(stages):
            if other == column:
                continue
            basis *= Polynomial([-nodes[other], 1.0]) / (
                nodes[column] - nodes[other]
            )
        primitive = basis.integ()
        coefficients[:, column] = primitive(nodes) - primitive(0.0)
        integrated_weights[column] = primitive(1.0) - primitive(0.0)
    if not np.allclose(integrated_weights, weights, atol=5.0e-14):
        raise RuntimeError("integrated Lagrange weights disagree with Gauss weights")
    return nodes, coefficients, integrated_weights


def stage_condensation_checks(
    data: LineData,
    step: float = 0.1,
    stage_orders: tuple[int, ...] = (1, 2, 3, 4),
) -> list[dict[str, float]]:
    """Eliminate active P_s -> P_{s-1} stages and recover the cap gauge."""
    generator = data.interconnection @ data.hessian
    state_dimension = generator.shape[0]
    identity = np.eye(state_dimension)
    rows: list[dict[str, float]] = []
    for stages in stage_orders:
        _, coefficients, weights = gauss_collocation_coefficients(stages)
        stage_matrix = (
            np.eye(stages * state_dimension)
            - step * np.kron(coefficients, generator)
        )
        injection = np.kron(np.ones((stages, 1)), identity)
        stage_solution = np.linalg.solve(stage_matrix, injection)
        cap_output = np.kron(weights.reshape(1, stages), generator)
        condensed_map = identity + step * cap_output @ stage_solution
        pade_map = gauss_pade_map(data, stages, step)
        induced_gauge = gauge_from_energy_map(data, step, condensed_map)
        cap_gauge = gauss_cap_gauge(data, stages, step)
        rows.append(
            {
                "stages": float(stages),
                "interior_dimension": float(stages * state_dimension),
                "stage_system_residual": np.linalg.norm(
                    stage_matrix @ stage_solution - injection, ord="fro"
                ),
                "pade_map_relative_error": np.linalg.norm(
                    condensed_map - pade_map, ord="fro"
                )
                / np.linalg.norm(pade_map, ord="fro"),
                "induced_gauge_relative_error": np.linalg.norm(
                    induced_gauge - cap_gauge, ord="fro"
                )
                / max(np.linalg.norm(cap_gauge, ord="fro"), 1.0),
                "induced_gauge_skew_defect": np.linalg.norm(
                    induced_gauge + induced_gauge.T, ord="fro"
                ),
                "energy_defect": np.linalg.norm(
                    condensed_map.T @ data.hessian @ condensed_map
                    - data.hessian,
                    ord="fro",
                )
                / np.linalg.norm(data.hessian, ord="fro"),
                "stage_system_condition": np.linalg.cond(stage_matrix),
            }
        )
    return rows


def hierarchy_checks(
    data: LineData,
    final_time: float = 0.5,
    stage_orders: tuple[int, ...] = (1, 2, 3),
    step_counts: tuple[int, ...] = (4, 8, 16, 32),
) -> list[dict[str, float]]:
    """Verify global orders 2, 4, and 6 for mixed Gauss-s gauges."""
    _, theta = mixed_generator(data)
    generator = data.interconnection @ data.hessian
    exact_map = expm(final_time * generator)
    rows: list[dict[str, float]] = []

    for order in stage_orders:
        previous_error: float | None = None
        for step_count in step_counts:
            step = final_time / step_count
            gauge = (
                gauss_cap_gauge(data, order, step)
                + step ** (2 * order - 1) * theta
            )
            one_step = cap_map(data, step, gauge)
            numerical_map = np.linalg.matrix_power(one_step, step_count)
            error = np.linalg.norm(
                numerical_map - exact_map, ord="fro"
            ) / np.linalg.norm(exact_map, ord="fro")
            rows.append(
                {
                    "stages": float(order),
                    "steps": float(step_count),
                    "step": step,
                    "global_error": error,
                    "observed_global_order": (
                        np.nan
                        if previous_error is None
                        else np.log(previous_error / error) / np.log(2.0)
                    ),
                    "energy_defect": np.linalg.norm(
                        one_step.T @ data.hessian @ one_step - data.hessian,
                        ord="fro",
                    )
                    / np.linalg.norm(data.hessian, ord="fro"),
                    "commutator": np.linalg.norm(
                        one_step @ generator - generator @ one_step,
                        ord="fro",
                    )
                    / (
                        np.linalg.norm(one_step, ord="fro")
                        * np.linalg.norm(generator, ord="fro")
                    ),
                }
            )
            previous_error = error
    return rows


def modal_frequencies(generator_or_map: np.ndarray, step: float | None = None) -> np.ndarray:
    """Positive generator frequencies or positive numerical phase rates."""
    eigenvalues = np.linalg.eigvals(generator_or_map)
    values = np.imag(eigenvalues) if step is None else np.angle(eigenvalues) / step
    return np.sort(values[values > 1.0e-9])


def dispersion_design_checks(
    data: LineData,
    reference_step: float = 0.0125,
    test_steps: tuple[float, ...] = (0.05, 0.1, 0.15),
) -> tuple[float, list[dict[str, float]]]:
    """Select one gauge coefficient once, then test unseen coarser steps."""
    gauss_coefficient, theta = mixed_generator(data)
    generator = data.interconnection @ data.hessian
    exact_frequencies = modal_frequencies(generator)
    identity = np.eye(generator.shape[0])

    def phase_error(alpha: float, step: float) -> float:
        gauge = step * gauss_coefficient + step**3 * alpha * theta
        frequencies = modal_frequencies(cap_map(data, step, gauge), step)
        return float(
            np.sqrt(
                np.mean(((frequencies - exact_frequencies) / exact_frequencies) ** 2)
            )
        )

    optimum = minimize_scalar(
        lambda alpha: phase_error(alpha, reference_step),
        bounds=(-100.0, 100.0),
        method="bounded",
        options={"xatol": 1.0e-10},
    )
    alpha = float(optimum.x)
    rows: list[dict[str, float]] = []
    for step in test_steps:
        base_gauge = step * gauss_coefficient
        selected_gauge = base_gauge + step**3 * alpha * theta
        base_matrix = (
            identity
            + step * data.interconnection @ base_gauge
            - 0.5 * step * generator
        )
        selected_matrix = (
            identity
            + step * data.interconnection @ selected_gauge
            - 0.5 * step * generator
        )
        base_error = phase_error(0.0, step)
        selected_error = phase_error(alpha, step)
        rows.append(
            {
                "step": step,
                "gauss_phase_error": base_error,
                "selected_phase_error": selected_error,
                "relative_reduction": 1.0 - selected_error / base_error,
                "gauss_condition": np.linalg.cond(base_matrix),
                "selected_condition": np.linalg.cond(selected_matrix),
            }
        )
    return alpha, rows


def coordinate_invariance_check(data: LineData) -> dict[str, float]:
    """Check the intrinsic two-form norm after a nonorthogonal chart change."""
    _, theta = mixed_generator(data)
    rng = np.random.default_rng(20260819)
    chart = np.eye(theta.shape[0]) + 0.08 * rng.normal(size=theta.shape)
    chart_inverse = np.linalg.inv(chart)
    transformed_metric = chart_inverse.T @ data.hessian @ chart_inverse
    transformed_theta = chart_inverse.T @ theta @ chart_inverse
    original_norm = two_form_norm(theta, data.hessian)
    transformed_norm = two_form_norm(transformed_theta, transformed_metric)
    return {
        "original_intrinsic_norm": original_norm,
        "transformed_intrinsic_norm": transformed_norm,
        "relative_norm_difference": abs(original_norm - transformed_norm)
        / original_norm,
    }


def open_cell_checks(
    data: LineData,
    step: float = 0.04,
) -> dict[str, dict[str, float]]:
    """Check the open power law for three cap gauges."""
    gauss_coefficient, mixed = mixed_generator(data)
    gauges = {
        "midpoint": np.zeros_like(mixed),
        "mixed_order_2": step * mixed,
        "gauss2_mixed_order_4": (
            step * gauss_coefficient + step**3 * mixed
        ),
    }
    rng = np.random.default_rng(20260816)
    initial = rng.normal(size=data.hessian.shape[0])
    boundary_input = np.array([0.2, -0.15])
    results: dict[str, dict[str, float]] = {}

    for name, gauge in gauges.items():
        left = np.eye(initial.size) - step * data.interconnection @ (
            0.5 * data.hessian - gauge
        )
        increment = np.linalg.solve(
            left,
            step
            * (
                data.interconnection @ data.hessian @ initial
                + data.input_map @ boundary_input
            ),
        )
        effort = (
            data.hessian @ (initial + 0.5 * increment) - gauge @ increment
        )
        output = (
            data.output_map @ effort + data.feedthrough @ boundary_input
        )
        final = initial + increment
        energy_increment = 0.5 * (
            final @ data.hessian @ final
            - initial @ data.hessian @ initial
        )
        supplied_energy = step * output @ boundary_input
        results[name] = {
            "energy_increment": energy_increment,
            "supplied_energy": supplied_energy,
            "balance_residual": abs(energy_increment - supplied_energy),
        }
    return results


def gauss_deformation_checks(
    data: LineData,
    steps: tuple[float, ...] = (0.1, 0.05, 0.025, 0.0125),
) -> list[dict[str, float]]:
    """Check order, energy, symplectic defect, and noncommutation."""
    gauss_coefficient, mixed = mixed_generator(data)
    generator = data.interconnection @ data.hessian
    symplectic_form = np.linalg.inv(data.interconnection)
    rows: list[dict[str, float]] = []
    previous_error: float | None = None

    for step in steps:
        base_gauge = step * gauss_coefficient
        mixed_gauge = base_gauge + step**3 * mixed
        base_map = cap_map(data, step, base_gauge)
        mixed_map = cap_map(data, step, mixed_gauge)
        exact_map = expm(step * generator)
        error = np.linalg.norm(mixed_map - exact_map, ord="fro") / np.linalg.norm(
            exact_map, ord="fro"
        )
        rows.append(
            {
                "step": step,
                "energy_defect": np.linalg.norm(
                    mixed_map.T @ data.hessian @ mixed_map - data.hessian,
                    ord="fro",
                )
                / np.linalg.norm(data.hessian, ord="fro"),
                "commutator": np.linalg.norm(
                    mixed_map @ generator - generator @ mixed_map,
                    ord="fro",
                )
                / (
                    np.linalg.norm(mixed_map, ord="fro")
                    * np.linalg.norm(generator, ord="fro")
                ),
                "gauss_commutator": np.linalg.norm(
                    base_map @ generator - generator @ base_map,
                    ord="fro",
                )
                / (
                    np.linalg.norm(base_map, ord="fro")
                    * np.linalg.norm(generator, ord="fro")
                ),
                "symplectic_defect": np.linalg.norm(
                    mixed_map.T @ symplectic_form @ mixed_map
                    - symplectic_form,
                    ord="fro",
                )
                / np.linalg.norm(symplectic_form, ord="fro"),
                "relative_flow_error": error,
                "observed_local_order": (
                    np.nan
                    if previous_error is None
                    else np.log(previous_error / error) / np.log(2.0)
                ),
            }
        )
        previous_error = error
    return rows


def two_cell_check(
    data: LineData,
    step: float = 0.04,
    dissipation: float = 0.02,
) -> dict[str, float]:
    """Solve two glued cells and check interface cancellation and passivity."""
    gauss_coefficient, mixed = mixed_generator(data)
    gauge = step * gauss_coefficient + step**3 * mixed
    resistance = dissipation * np.eye(data.hessian.shape[0])
    size = data.hessian.shape[0]
    rng = np.random.default_rng(20260817)
    initial_left = rng.normal(size=size)
    initial_right = rng.normal(size=size)
    external_left = 0.12
    external_right = -0.08
    total_size = 2 * size + 2

    def residual(unknown: np.ndarray) -> np.ndarray:
        increment_left = unknown[:size]
        increment_right = unknown[size : 2 * size]
        internal_right_left = unknown[-2]
        internal_left_right = unknown[-1]
        effort_left = (
            data.hessian @ (initial_left + 0.5 * increment_left)
            - gauge @ increment_left
        )
        effort_right = (
            data.hessian @ (initial_right + 0.5 * increment_right)
            - gauge @ increment_right
        )
        input_left = np.array([external_left, internal_right_left])
        input_right = np.array([internal_left_right, external_right])
        output_left = (
            data.output_map @ effort_left + data.feedthrough @ input_left
        )
        output_right = (
            data.output_map @ effort_right + data.feedthrough @ input_right
        )
        return np.concatenate(
            (
                increment_left
                - step
                * (
                    (data.interconnection - resistance) @ effort_left
                    + data.input_map @ input_left
                ),
                increment_right
                - step
                * (
                    (data.interconnection - resistance) @ effort_right
                    + data.input_map @ input_right
                ),
                [internal_right_left - output_right[0]],
                [internal_left_right + output_left[1]],
            )
        )

    origin_residual = residual(np.zeros(total_size))
    system = np.column_stack(
        [
            residual(np.eye(total_size)[:, column]) - origin_residual
            for column in range(total_size)
        ]
    )
    unknown = np.linalg.solve(system, -origin_residual)

    increment_left = unknown[:size]
    increment_right = unknown[size : 2 * size]
    effort_left = (
        data.hessian @ (initial_left + 0.5 * increment_left)
        - gauge @ increment_left
    )
    effort_right = (
        data.hessian @ (initial_right + 0.5 * increment_right)
        - gauge @ increment_right
    )
    input_left = np.array([external_left, unknown[-2]])
    input_right = np.array([unknown[-1], external_right])
    output_left = data.output_map @ effort_left + data.feedthrough @ input_left
    output_right = (
        data.output_map @ effort_right + data.feedthrough @ input_right
    )

    final_left = initial_left + increment_left
    final_right = initial_right + increment_right
    energy_increment = 0.5 * (
        final_left @ data.hessian @ final_left
        - initial_left @ data.hessian @ initial_left
        + final_right @ data.hessian @ final_right
        - initial_right @ data.hessian @ initial_right
    )
    external_power = (
        output_left[0] * input_left[0]
        + output_right[1] * input_right[1]
    )
    interface_power = (
        output_left[1] * input_left[1]
        + output_right[0] * input_right[0]
    )
    dissipated_power = (
        effort_left @ resistance @ effort_left
        + effort_right @ resistance @ effort_right
    )
    predicted_increment = step * (external_power - dissipated_power)
    return {
        "linear_solve_residual": np.linalg.norm(residual(unknown)),
        "interface_power": interface_power,
        "energy_increment": energy_increment,
        "external_energy": step * external_power,
        "dissipated_energy": step * dissipated_power,
        "balance_residual": abs(energy_increment - predicted_increment),
    }


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _last_rows_by_stage(rows: list[dict[str, float]]) -> dict[int, dict[str, float]]:
    result: dict[int, dict[str, float]] = {}
    for row in rows:
        result[int(row["stages"])] = row
    return result


def run_assertions() -> int:
    """Execute tolerance and convergence assertions for every report block."""
    data = shifted_legendre_line()
    count = 0

    def check(condition: bool, message: str) -> None:
        nonlocal count
        _require(condition, message)
        count += 1

    check(max(structural_residuals(data).values()) < 1.0e-10,
          "rectangular spatial identities")

    reconstruction = rectangular_reconstruction_checks(data)
    residual_keys = [key for key in reconstruction if "condition" not in key and "singular" not in key]
    check(max(reconstruction[key] for key in residual_keys) < 1.0e-10,
          "rectangular-to-square reconstruction")
    check(reconstruction["left_endpoint_min_singular"] > 0.1 and
          reconstruction["right_endpoint_min_singular"] > 0.1,
          "endpoint reconstruction nonsingularity")

    reduction = volume_cap_reduction_checks(data)
    check(max(reduction[key] for key in (
        "trace_right_inverse", "volume_form_skew", "pullback_recovery"
    )) < 1.0e-10, "volume-to-cap pullback")

    condensation = dynamic_condensation_checks(data)
    check(max(condensation[key] for key in condensation if "condition" not in key) < 1.0e-10,
          "raw trace-generated condensation")

    chart = cayley_chart_checks(data)
    check(max(chart.values()) < 1.0e-10, "Cayley chart roundtrip")
    conditioning = cayley_conditioning_checks(data)
    check(max(row["stable_map_relative_error"] for row in conditioning) < 1.0e-12,
          "stable Cayley evaluation")

    stage_rows = stage_condensation_checks(data)
    check(max(row[key] for row in stage_rows for key in (
        "stage_system_residual", "pade_map_relative_error",
        "induced_gauge_relative_error", "induced_gauge_skew_defect",
        "energy_defect"
    )) < 1.0e-10, "active-stage condensation")

    open_checks = open_cell_checks(data)
    check(max(row["balance_residual"] for row in open_checks.values()) < 1.0e-12,
          "open cap-mantle power")

    deformation = gauss_deformation_checks(data)
    check(deformation[-1]["observed_local_order"] > 4.8 and
          deformation[-1]["commutator"] > 1.0e-10,
          "fourth-order noncommuting deformation")

    hierarchy = _last_rows_by_stage(hierarchy_checks(data))
    check(all(hierarchy[s]["observed_global_order"] > 2 * s - 0.1 for s in hierarchy),
          "mixed Gauss hierarchy orders")

    raised_rows, gauss_two_identity = order_raising_checks(data)
    raised = _last_rows_by_stage(raised_rows)
    check(all(raised[s]["observed_global_order"] > 2 * s + 2 - 0.1 for s in raised),
          "unique order-raising hierarchy")
    check(gauss_two_identity < 1.0e-12, "exceptional Gauss-1/Gauss-2 identity")

    coefficients = asymptotic_coefficient_checks()
    last_coefficients = _last_rows_by_stage(coefficients)
    check(all(row["relative_component_error"] < 5.0e-4 for row in last_coefficients.values()),
          "direct d_s coefficients")

    two_term = _last_rows_by_stage(two_term_hierarchy_checks(data))
    check(all(two_term[s]["observed_global_order"] > 2 * s + 2 - 0.1 for s in two_term),
          "two-term noncommuting hierarchy orders")

    second_central = _last_rows_by_stage(second_central_hierarchy_checks(data))
    check(all(second_central[s]["observed_global_order"] > 2 * s + 4 - 0.1
              for s in second_central), "second central correction orders")

    decimal_rows = decimal_noncommutation_checks()
    check(decimal_rows[-1]["observed_exponent"] > 8.9,
          "multiprecision s=3 commutator exponent")
    check(decimal_rows[-1]["leading_relative_error"] < 0.02,
          "multiprecision s=3 leading commutator matrix")

    rational = rational_reduction_checks()
    check(rational["c3_coordinate_denominator"] < 1.0e-14 and
          rational["pi31_at_coordinate_pole_error"] < 1.0e-12 and
          rational["continued_map_at_coordinate_pole_error"] < 1.0e-12,
          "removable C3 coordinate resonance")
    check(rational["left_half_plane_pole"] < -10.0 and
          rational["left_pole_denominator_residual"] < 1.0e-10 and
          rational["left_pole_numerator_magnitude"] > 1.0,
          "genuine left-half-plane pole of R_2,1")
    check(rational["imaginary_axis_unit_modulus_defect"] < 1.0e-12,
          "imaginary-axis unit amplification after rational reduction")

    sign_evidence = coefficient_sign_checks()
    check(sign_evidence["sign_failures"] == 0.0 and
          sign_evidence["alpha_3_2_relative_error"] == 0.0 and
          sign_evidence["alpha_4_3_relative_error"] == 0.0,
          "exact evidence for the alpha sign conjecture")

    schur_generation = interior_schur_generation_checks()
    check(max(schur_generation[key] for key in (
              "full_generator_skew_defect", "effective_gauge_error",
              "eliminated_response_error")) < 1.0e-12 and
          schur_generation["target_change_rank"] ==
          schur_generation["interior_dimension"],
          "rank-minimal alternating interior Schur generation")

    stieltjes_generation = stieltjes_interior_generation_checks()
    check(stieltjes_generation["gauss3_stieltjes_realization_error"] < 1.0e-14 and
          stieltjes_generation["gauss3_coordinate_generation_error"] < 1.0e-14 and
          stieltjes_generation["gauss3_weight"] > 0.0 and
          stieltjes_generation["gauss3_spectral_point"] > 0.0,
          "positive Gauss-3 Stieltjes interior generation")
    check(stieltjes_generation["corrected_21_shifted_value_at_u61"] < 0.0 and
          abs(stieltjes_generation["corrected_21_zero_at_u60"]) < 1.0e-15,
          "positive-realization obstruction for corrected (2,1)")

    _, dispersion = dispersion_design_checks(data)
    check(all(row["relative_reduction"] > 0 for row in dispersion),
          "out-of-sample dispersion reduction")
    check(coordinate_invariance_check(data)["relative_norm_difference"] < 1.0e-12,
          "intrinsic gauge normalization")

    glued = two_cell_check(data)
    check(glued["linear_solve_residual"] < 1.0e-12 and
          abs(glued["interface_power"]) < 1.0e-12 and
          glued["balance_residual"] < 1.0e-12,
          "dissipative two-cell interface balance")
    return count


def print_report() -> None:
    data = shifted_legendre_line()
    print("RECTANGULAR 2012 SPATIAL DATA")
    for name, value in structural_residuals(data).items():
        print(f"{name:24s} = {value:.6e}")

    print("\nRECTANGULAR-TO-SQUARE DIRAC RECONSTRUCTION")
    for name, value in rectangular_reconstruction_checks(data).items():
        print(f"{name:32s} = {value:.6e}")

    print("\nVOLUME-TO-CAP REDUCTION")
    for name, value in volume_cap_reduction_checks(data).items():
        print(f"{name:24s} = {value:.6e}")

    print("\nRAW TRACE-GENERATED CELL REALIZATION")
    for name, value in dynamic_condensation_checks(data).items():
        print(f"{name:32s} = {value:.6e}")

    print("\nCAYLEY CHART ON THE ENERGY-ORTHOGONAL GROUP")
    for name, value in cayley_chart_checks(data).items():
        print(f"{name:36s} = {value:.6e}")

    print("\nCAYLEY EVALUATION NEAR THE EXCLUDED IDENTITY")
    for row in cayley_conditioning_checks(data):
        print(
            f"h={row['step']:.1e}: "
            f"naive_gauge_err={row['naive_gauge_relative_error']:.3e}, "
            f"scaled={row['epsilon_scaled_inverse_error']:.3f}, "
            f"stable_map_err={row['stable_map_relative_error']:.3e}, "
            f"stable_skew={row['stable_gauge_skew_defect']:.3e}"
        )

    print("\nACTIVE P_s -> P_{s-1} GAUSS CELL CONDENSATION")
    for row in stage_condensation_checks(data):
        print(
            f"s={int(row['stages'])}, interior={int(row['interior_dimension']):2d}: "
            f"stage_res={row['stage_system_residual']:.3e}, "
            f"Pade={row['pade_map_relative_error']:.3e}, "
            f"gauge={row['induced_gauge_relative_error']:.3e}, "
            f"skew={row['induced_gauge_skew_defect']:.3e}, "
            f"energy={row['energy_defect']:.3e}, "
            f"cond={row['stage_system_condition']:.3f}"
        )

    print("\nOPEN CELL POWER")
    for name, values in open_cell_checks(data).items():
        print(
            f"{name:22s}: dH={values['energy_increment']:+.12e}, "
            f"h y.u={values['supplied_energy']:+.12e}, "
            f"res={values['balance_residual']:.6e}"
        )

    print("\nFOURTH-ORDER MIXED DEFORMATION OF GAUSS-2")
    for row in gauss_deformation_checks(data):
        print(
            f"h={row['step']:.5f}: energy={row['energy_defect']:.3e}, "
            f"comm={row['commutator']:.3e}, "
            f"comm_G2={row['gauss_commutator']:.3e}, "
            f"sympl={row['symplectic_defect']:.3e}, "
            f"flow_err={row['relative_flow_error']:.3e}, "
            f"local_order={row['observed_local_order']:.3f}"
        )

    print("\nGAUSS-s MIXED HIERARCHY: GLOBAL CONVERGENCE")
    for row in hierarchy_checks(data):
        print(
            f"s={int(row['stages'])}, n={int(row['steps']):2d}, "
            f"h={row['step']:.6f}: error={row['global_error']:.3e}, "
            f"global_order={row['observed_global_order']:.3f}, "
            f"energy={row['energy_defect']:.3e}, "
            f"comm={row['commutator']:.3e}"
        )

    raised_rows, gauss_two_identity = order_raising_checks(data)
    print("\nUNIQUE ORDER-RAISING GAUGE")
    for row in raised_rows:
        print(
            f"s={int(row['stages'])}, n={int(row['steps']):2d}: "
            f"error={row['global_error']:.3e}, "
            f"global_order={row['observed_global_order']:.3f}, "
            f"energy={row['energy_defect']:.3e}"
        )
    print(f"s=1 raised map vs Gauss-2 = {gauss_two_identity:.3e}")

    print("\nDIRECT LOCAL COEFFICIENT d_s")
    for row in asymptotic_coefficient_checks():
        print(
            f"s={int(row['stages'])}, h={row['step']:.3f}: "
            f"pred={row['predicted_coefficient']:+.12e}, "
            f"est={row['estimated_coefficient']:+.12e}, "
            f"component_rel={row['relative_component_error']:.3e}, "
            f"full_rel={row['relative_coefficient_error']:.3e}"
        )

    print("\nDECIMAL s=3 NONCOMMUTATION AT FINE STEPS")
    for row in decimal_noncommutation_checks():
        print(
            f"h={row['step']:.6f}: comm={row['commutator_norm']:.3e}, "
            f"exponent={row['observed_exponent']:.5f}, "
            f"leading_rel={row['leading_relative_error']:.3e}"
        )

    print("\nREDUCED RATIONAL HIERARCHY")
    for name, value in rational_reduction_checks().items():
        print(f"{name:42s} = {value:.12e}")

    print("\nEXACT COEFFICIENT-SIGN EVIDENCE")
    for name, value in coefficient_sign_checks().items():
        print(f"{name:42s} = {value:.12e}")

    print("\nMINIMAL ALTERNATING INTERIOR SCHUR GENERATION")
    for name, value in interior_schur_generation_checks().items():
        print(f"{name:42s} = {value:.12e}")

    print("\nPOSITIVE STIELTJES INTERIOR GENERATION")
    for name, value in stieltjes_interior_generation_checks().items():
        print(f"{name:42s} = {value:.12e}")

    print("\nTWO-TERM HIGH-ORDER NONCOMMUTING HIERARCHY")
    for row in two_term_hierarchy_checks(data):
        print(
            f"s={int(row['stages'])}, n={int(row['steps']):2d}: "
            f"error={row['global_error']:.3e}, "
            f"global_order={row['observed_global_order']:.3f}, "
            f"energy={row['energy_defect']:.3e}, "
            f"comm={row['commutator']:.3e}"
        )

    print("\nSECOND CENTRAL CORRECTOR")
    for row in second_central_hierarchy_checks(data):
        print(
            f"s={int(row['stages'])}, n={int(row['steps']):2d}: "
            f"error={row['global_error']:.3e}, "
            f"global_order={row['observed_global_order']:.3f}, "
            f"energy={row['energy_defect']:.3e}"
        )

    alpha, dispersion_rows = dispersion_design_checks(data)
    print("\nONE-PARAMETER DISPERSION DESIGN")
    print(f"selected_alpha_at_h=1/80 = {alpha:.12e}")
    for row in dispersion_rows:
        print(
            f"h={row['step']:.5f}: Gauss={row['gauss_phase_error']:.3e}, "
            f"selected={row['selected_phase_error']:.3e}, "
            f"reduction={100.0 * row['relative_reduction']:.2f}%, "
            f"cond_G={row['gauss_condition']:.3f}, "
            f"cond_sel={row['selected_condition']:.3f}"
        )

    print("\nINTRINSIC GAUGE NORMALIZATION")
    for name, value in coordinate_invariance_check(data).items():
        print(f"{name:32s} = {value:.12e}")

    print("\nTWO GLUED CELLS WITH DISSIPATION")
    for name, value in two_cell_check(data).items():
        print(f"{name:24s} = {value:+.12e}")


if __name__ == "__main__":
    assertion_count = run_assertions()
    print_report()
    print(f"\nALL {assertion_count} ASSERTION GROUPS PASSED")
