# Space–Time Stokes–Dirac Discretization

**Affine pairing gauges, active temporal stages and port-Hamiltonian simulation**  
Redha Moulla · AXIA, France · 2026

This single manuscript contains the complete version-12 geometric theory and its
extension to active temporal efforts, nonlinear Hamiltonians, physical dissipation
and boundary ports. Proofs, simulation protocols, results and references are
available directly on GitHub. Installation and reproduction commands are in the
[repository README](../README.md).

Sections 1–12 retain the original manuscript's numbering and hypotheses.
Sections 13–23 develop the active-stage extension; their statements and citations
use the prefix **T** to keep the two numbering systems unambiguous.
The [original LaTeX source](../archive/version-12/manuscript.tex) is retained for
provenance. The geometric construction begins with the spatial method of
[Moulla, Lefèvre and Maschke (2012)](https://doi.org/10.1016/j.jcp.2011.10.008).

**Reading routes:** start with [the original classification](#v12-section-2) for
the geometric theory, [the active-stage formulation](#extension-section-14) for
the nonlinear solver, or [the numerical evidence](#extension-section-22) for
reproducible simulations. Energy identities, temporal accuracy and simultaneous
space-time PDE convergence are distinguished throughout.

![From spatial Stokes geometry to active temporal stages and a discrete PH power balance](assets/overview.svg)

## Contents

**Geometric foundations and conservative cap construction**

1. [Introduction](#v12-section-1)
2. [The affine Stokes-pairing principle](#v12-section-2)
3. [Stokes' theorem on the space–time cylinder](#v12-section-3)
4. [Discrete cylinders and the full space–time gauge](#v12-section-4)
5. [The complete two-point temporal gauge](#v12-section-5)
6. [The exact symplectic gauge locus](#v12-section-6)
7. [Nonlinear Hamiltonians and the universal segment rule](#v12-section-7)
8. [Multiple stages: classification for invertible stage matrices](#v12-section-8)
9. [Space–time Stokes–Dirac closure and local rigidity](#v12-section-9)
10. [Consequences and discriminating computational validation](#v12-section-10)
11. [Discussion and limitations](#v12-section-11)
12. [Conclusion](#v12-section-12)

**Active-stage extension: nonlinear systems, dissipation and validation**

13. [Purpose and precise scope](#extension-section-13)
14. [Active time cell and constitutive effort](#extension-section-14)
15. [Closure, ports and exact energy balance](#extension-section-15)
16. [Exact recovery of the version-12 mixed cap family](#extension-section-16)
17. [The dissipative obstruction to one-effort compression](#extension-section-17)
18. [Elimination preserves passivity, including the ports](#extension-section-18)
19. [Order, reversal and coordinate covariance](#extension-section-19)
20. [The spatial reduction and a geometrically specified test gauge](#extension-section-20)
21. [Interconnection](#extension-section-21)
22. [Numerical evidence and reproduction](#extension-section-22)
23. [What remains open](#extension-section-23)

[References](#references) · [Reproducibility and provenance](#reproducibility-and-provenance)

---

## Abstract

We classify compatible pairings on a finite-dimensional discrete
space–time cylinder.  For a possibly nonsurjective differential
$`D:E\to F`$, compatibility exists precisely when the boundary form
annihilates $`\ker D`$, and the active freedom is
$`\Lambda^2(E/\ker D)^*`$.  On a temporal–spatial complex this separates
temporal, spatial, and intrinsic mixed gauges.  A regular interior solve
then yields an exact sequence for their visible and invisible cap
reductions.

For nonsingular quadratic Hamiltonian systems, cap gauges are shifted
Cayley coordinates on the no-$`+1`$ part of the energy-orthogonal group;
we determine their exact symplectic locus and defect.  An active polynomial
cell $`\mathcal P_s\to\mathcal P_{s-1}`$ generates the Gauss–$`s`$ cap
map by elimination.  On the 2012 polynomial transmission line, mixed
cylinder lifts preserve cap–mantle power, dissipation, and interface
cancellation.  Every alternating cap gauge also admits a minimal
auxiliary-mode Schur generation whose interior dimension equals its rank;
this algebraic realization is distinguished from selection by an
independently prescribed field cell.  For each fixed semidiscrete generator, a unique central
coefficient raises operator convergence from $`2s`$ to $`2s+2`$; a
recursive rational hierarchy gives order at least $`2s+2r`$, exact when
the next coefficient is nonzero, while
noncentral additions retain the attained order and leave the scalar
functional calculus.  Its reduced polynomial form separates removable
Cayley-chart resonances from true poles and shows unit amplification on
oscillatory modes but not general $`A`$-stability.  The exact coordinate is
a Stieltjes function in $`z^2`$, and $`C_s`$ is its diagonal or
near-diagonal Padé approximant.  These Gauss coordinates have positive
Jacobi-chain realizations, whereas the corrected hierarchy need not: the
first nontrivial member already violates Stieltjes positivity.  Executable tests verify the reductions, power
laws, asymptotic coefficients, high-order branches, and a multiprecision
noncommutation law.

**Keywords.**
Port-Hamiltonian systems; Stokes–Dirac structures; space–time complexes;
discrete gradients; Gauss collocation; symplecticity; power balance;
geometric numerical integration.

**MSC 2020.**
37M15; 65P10; 65M60; 93C20.

<a name="v12-section-1"></a>

## 1. Introduction

Distributed port-Hamiltonian systems expose boundary power as part of the
geometry rather than as an afterthought
[[33]](#v12-ref-vanderschaft2002), [[17]](#v12-ref-kotyczka2018).  On a space–time cylinder, the
initial and final energies and the lateral power ports are different
components of one oriented boundary.  A discretization that starts from
this cylinder should therefore not be assembled by merely attaching a
time index to a spatial Stokes matrix.  It should possess a total discrete
differential, a global boundary form, and a compatible pairing before a
constitutive law is selected.

Existing structure-preserving spatial discretizations construct a
finite-dimensional Stokes identity whose boundary term becomes a power
port [[20]](#v12-ref-moulla2012), [[29]](#v12-ref-seslija2012), [[4]](#v12-ref-brugnoli2022).  In particular,
discrete-exterior-calculus constructions produce simplicial Dirac
structures on space while retaining continuous time; treating space and
time on equal footing is a distinct extension already identified in that
literature [[29]](#v12-ref-seslija2012).  Passivity-consistent time discretizations
and energy-consistent
Petrov–Galerkin formulations provide complementary temporal closures
[[5]](#v12-ref-celledoni2017), [[12]](#v12-ref-giesselmann2025).  Symplectic Runge–Kutta theory and
collocation identify Gauss–Legendre methods as the canonical collocation
schemes preserving quadratic invariants
[[8]](#v12-ref-cooper1987), [[28]](#v12-ref-sanzserna1988), [[27]](#v12-ref-reich2000), [[14]](#v12-ref-hairer2006).  Discrete-gradient
methods instead preserve a prescribed Hamiltonian exactly, including for
nonlinear Hamiltonians
[[21]](#v12-ref-mclachlan1999), [[13]](#v12-ref-gonzalez1996), [[26]](#v12-ref-quispel2008), [[10]](#v12-ref-eidnes2022).  High-order
segment and line-integral constructions connect this viewpoint to
collocation and Gaussian quadrature
[[15]](#v12-ref-hairer2010), [[3]](#v12-ref-brugnano2010), [[25]](#v12-ref-miyatake2016), [[16]](#v12-ref-kemmochi2023).  Locally exact
discrete-gradient modifications likewise use step- and
linearization-dependent skew operators to reproduce a target linearized
flow while preserving energy [[9]](#v12-ref-cieslinski2013).  The overlap is the
use of a constitutively adapted skew operator.  The present hierarchy is
not claimed as the first such modification: its additional content is the
affine Stokes classification before closure, the exact volume-to-cap
sequence, the central/noncentral split, and the realization of either
direction by a mixed cylinder two-form.

There are established genuinely space–time routes with different primary
structures.  Covariant variational integrators derive a discrete
multisymplectic form formula from a discrete action
[[23]](#v12-ref-marsden1998); space–time finite elements provide arbitrary-order
energy and momentum laws for multisymplectic PDEs
[[6]](#v12-ref-celledoni2021); and HDG/FEEC formulations characterize local
multisymplecticity on numerical traces
[[22]](#v12-ref-mclachlan2020), [[31]](#v12-ref-stern2026), [[32]](#v12-ref-sternzampa2026).  The present construction is
not a claim of priority over those field discretizations and does not
assume a discrete action.  Its different question is to classify all
pairings compatible with a prescribed discrete Stokes formula, transport
their affine freedom through a cell reduction, and identify which mixed
gauges generate boundary maps outside the scalar method-of-lines class.

The common algebraic fact is simple but consequential: a Stokes identity
or a discrete chain rule fixes only the symmetric part of a pairing.  The
remaining freedom is alternating on the quotient by directions invisible
to the differential.  The central question is whether this observation
can be made directly on the space–time cylinder, so that the spatial
freedom found in [[20]](#v12-ref-moulla2012), the temporal discrete-gradient gauge,
and mixed nonseparable pairings become sectors of one object.

This paper answers that question at eight levels.

1. For an arbitrary, possibly nonsurjective, differential
     $`D:E\to F`$, the active compatible pairings form an affine space over
     $`\Lambda^2(E/\ker D)^*`$.  Applied to the genuine total complex through
     a degree-respecting symmetric double, this gives a bilinear gauge between
     complementary active quotients.  Applied instead to the stacked
     one-bidegree cylinder map, it gives temporal, spatial, and mixed sectors.

1. Once a Dirac relation, constitutive closure, boundary causality,
     and regular interior solve have fixed a reconstruction
     $`L_h:Q_\partial\to Q_\mathrm{st}`$, the exact volume-to-cap reduction is
     $`\omega\mapsto L_h^*\omega`$.  We classify its image and kernel and
     derive an intrinsic rank criterion for the visibility of mixed gauges.
     For trace-generated forms, a raw weak cell containing the same volume
     matrix is assembled and exactly realizes the prescribed cap gauge.  We
     state explicitly that this is a realization theorem, not a selection of
     the gauge by independent interior dynamics.

1. For a single time step and a quadratic Hamiltonian, the family
     $`\Lambda^2V^*`$ gives constitutively shifted Cayley coordinates on the
     no-$`+1`$ cell of $`O(Q)`$.  Exact power holds on the whole chart; its symplectic
     intersection and factorized defect are computed explicitly, and the
     exact-flow gauge is obtained in closed form.

1. For an $`s`$-stage Runge–Kutta reconstruction with invertible
     stage matrix, all compatible bilinear stage pairings are classified.
     Local diagonal pairing reduces compatibility to the Cooper condition,
     and a short spanning argument proves directly that symplectic
     collocation is precisely Gauss collocation.  Conversely, starting from
     the polynomial temporal cell
     $`\mathcal P_s\to\mathcal P_{s-1}`$, we retain the Gauss stages as
     active interior variables and prove that their elimination produces the
     diagonal Padé cap map and its induced two-point gauge.

1. For a space–time Stokes–Dirac discretization, Dirac maximality is
     kept separate from the constitutive chain rule.  This yields exact cap
     energy, lateral boundary power, passivity with dissipation, and exact
     cancellation on glued faces.  Under scalar locality, the stacked tensor
     equations force factorization, and maximal exactness then selects
     Gauss–Legendre on each axis.

1. A general rectangular Stokes theorem shows that an invertible
     effort-and-input reconstruction produces a square open representation
     with skew interconnection, dual input/output blocks, and skew
     feedthrough.  This turns the spatial reduction used in the 2012 line
     from a numerical observation into a direct algebraic consequence.

1. On that rectangular transmission-line reduction, a two-form on
     spatial derivatives of temporal increments gives an explicit mixed
     cylinder gauge.  For every $`s`$, its cap reduction gives a
     reversible, energy-preserving deformation of Gauss–$`s`$.  In the
     fixed one-term ansatz, one central coefficient and only one raises the
     order from $`2s`$ to exactly $`2s+2`$; every other fixed coefficient has
     exactly order $`2s`$.  Adding the next odd-power gauge term preserves
     that higher order while permitting noncommutation, and a recursive
     central hierarchy yields successive two-order lower bounds, exact under
     an explicit non-annulation condition.  A commutator
     criterion proves when the resulting method cannot belong to the scalar
     functional calculus of the semidiscrete line.  A dispersion experiment
     turns the gauge into a design variable.

1. Every finite-rank alternating cap gauge has a minimal auxiliary
     Schur generator with interior dimension equal to its rank.  For the
     central scalar coordinates, Stieltjes–Padé theory gives a stronger
     positive Jacobi-chain realization of every Gauss member.  The corrected
     $`(s,r)=(2,1)`$ coordinate fails this positivity test, separating
     unrestricted algebraic generation from passive interior generation.

The last level exposes the following structural chain without identifying
its distinct notions of positivity:

<a name="v12-eq-gauss-stieltjes-jacobi-chain"></a>
**Equation (1).**

```math

 \boxed{\begin{gathered}
 \text{Gauss}
 \ \Longleftrightarrow\ 
 \text{Padé--Stieltjes}
 \ \Longleftrightarrow\ 
 \text{finite positive measure}\\
 \Longleftrightarrow\ 
 \text{positive Jacobi resolvent}
 \ \Longleftrightarrow\ 
 \text{Stieltjes continued fraction}.
 \end{gathered}}
```

**Scope and novelty.**

The affine Stokes principle, classical Gauss quadrature, symplectic
Runge–Kutta theory, discrete gradients, and the Dirac construction are
established ingredients.  The new contribution is their common
quotient-space formulation on a discrete cylinder: the nonsurjective
classification, the explicit distinction between the total and stacked
differentials, the exact visibility sequence for volume-to-cap reduction,
the intrinsic mixed reduction, the constitutively shifted Cayley
coordinates of cap gauges, the
rectangular reconstruction theorem, the active-stage derivation of the
Gauss cap gauge, the sharp auxiliary Schur existence-and-minimality
construction, the
positive Jacobi realization and its obstruction, the symplectic-locus and defect formulas, the exact
one-term order dichotomy, its multi-term correction hierarchy, tensor
rigidity derived from the global Stokes equations, and the mixed
Gauss–$`s`$ family with its scalar-functional-calculus obstruction.  The
trace-generated theorem is deliberately described as exact realization;
only the higher-degree stage cell derives its cap gauge from an
independent interior solve.  We do not call an
arbitrary compatible bilinear form a multi-Dirac structure.  Multi-Dirac
field theory includes additional graded geometric and integrability
requirements [[34]](#v12-ref-vankerschaver2012).

<a name="v12-fig-paper-architecture"></a>

```mermaid
flowchart TD
    V["Compatible volume gauge"] --> C["Regular pullback to cap gauge"]
    C --> M["Constitutive closure: shifted Cayley map"]
    M --> A["Central series: symplectic, higher order"]
    M --> B["Noncentral directions: noncommuting, same power"]
```

The volume gauge is $`\omega\in\Lambda^2Q_{\mathrm{st}}^*`$ and the cap pullback is $`S_h=L_h^*\omega`$.

**Figure 1.** Architecture of the construction.  Compatibility is classified
in the volume, transported to the caps, and closed through a Cayley chart.
Only after the constitutive generator $`A_N=\mathcal J_NQ_N`$ is fixed
does the cap space split operationally into a central accuracy hierarchy
and noncentral design directions; both retain the same contracted power
identity.

**Order convention.**

Every order statement below concerns matrix/operator convergence to
$`e^{hA_N}`$ as $`h\to0`$ for a *fixed finite-dimensional linear
semidiscrete generator* $`A_N`$, with constants allowed to depend on
$`N`$.  It is neither a universal Runge–Kutta order statement for
nonlinear vector fields nor a simultaneous space–time convergence
theorem.  The nonlinear section establishes the chain rule and energy
identity, not the linear hierarchy's order.

**Definition 1.1. Local matrix defect and global operator order.**

<a name="v12-def-matrix-order"></a>

For fixed $`A_N`$, a one-step family $`M_h`$ has *local matrix defect
of order $`p+1`$* if

<a name="v12-eq-local-matrix-order"></a>
**Equation (2).**

```math

 \left\lVert M_h-e^{hA_N}\right\rVert=O(h^{p+1})\qquad(h\to0)
```

in any matrix norm.  It has *global operator order $`p`$* on bounded
time intervals if, for each fixed $`T>0`$ and $`n=T/h\in\mathbb N`$,

<a name="v12-eq-global-operator-order"></a>
**Equation (3).**

```math

 \left\lVert M_h^n-e^{TA_N}\right\rVert=O(h^p).
```

For a consistent finite-dimensional family satisfying
[Equation (2)](#v12-eq-local-matrix-order), the telescoping identity

```math
M_h^n-e^{TA_N}
 =\sum_{k=0}^{n-1}M_h^{n-1-k}
  (M_h-e^{hA_N})e^{khA_N}
```

gives [Equation (3)](#v12-eq-global-operator-order).  Thus the one-power shift used
below is an explicit fixed-$`A_N`$ matrix estimate.

**Terminology for ``exact.”**

The adjective is always qualified by its object.  *Exact power* or
*exact energy* means an algebraic identity at the discrete level;
*exact elimination* means equality of two finite-dimensional Schur
reductions; *exact order* means that the displayed leading error
coefficient is nonzero; and *exact-flow gauge* refers specifically to
the Cayley coordinate of $`e^{hA_N}`$.  None of these phrases asserts a
simultaneous continuum limit unless such a limit is stated separately.

**Notation.**

All spaces are finite-dimensional and real. A vector space $`V`$ and its
dual $`V^*`$ are distinguished until coordinates are chosen. In
coordinates, covectors are represented by columns and paired by transpose.
For an invertible skew matrix $`J`$, we write $`\Omega=J^{-1}`$.
The state dimension is $`n`$; in symplectic statements $`n=2m`$.

<a name="v12-section-2"></a>

## 2. The affine Stokes-pairing principle

We first isolate the algebra that is shared by spatial integration by parts
and the temporal discrete chain rule.

**Definition 2.1. Discrete Stokes data.**

Let $`E,F`$ be vector spaces, let $`D:E\to F`$ be linear, and let
$`b:E\times E\to\mathbb{R}`$ be a symmetric bilinear boundary form. A bilinear
map $`M:E\times F\to\mathbb{R}`$ is *Stokes compatible* if

<a name="v12-eq-stokes-compatibility"></a>
**Equation (4).**

```math

 M(z,Dw)+M(w,Dz)=b(z,w)
 \qquad (z,w\in E).
```

**Theorem 2.2. Active affine Stokes family.**

<a name="v12-thm-affine-stokes"></a>

Let $`K=\ker D`$, $`I=\operatorname{im}D`$, and
$`r=\operatorname{rank}D`$. A compatible pairing exists if and only if

<a name="v12-eq-kernel-obstruction"></a>
**Equation (5).**

```math

 b(k,\ell)=0\qquad(k,\ell\in K).
```

On the active flow space $`I`$, if $`M_0`$ is one compatible pairing,
all compatible restrictions are

<a name="v12-eq-abstract-gauge"></a>
**Equation (6).**

```math

 M_\omega(z,Dw)=M_0(z,Dw)+\omega([z],[w]),
 \qquad \omega\in\Lambda^2(E/K)^*.
```

Consequently the active solution set is an affine space of dimension

<a name="v12-eq-active-dimension"></a>
**Equation (7).**

```math

 \frac{r(r-1)}2.
```

If pairings are required on all of $`E\times F`$, their values transverse
to $`I`$ are unconstrained.  If $`\mathcal H`$ denotes the full
homogeneous direction, it fits into the short exact sequence

<a name="v12-eq-full-affine-direction"></a>
**Equation (8).**

```math

 0\longrightarrow \operatorname{Bil}(E,F/I;\mathbb{R})
 \longrightarrow\mathcal H
 \longrightarrow\Lambda^2(E/K)^*
 \longrightarrow0.
```

Every choice of a complement of $`I`$ in $`F`$ splits this sequence, but
the splitting is not canonical.  The full direction has dimension
$`r(r-1)/2+\dim(E)(\dim(F)-r)`$.

*Proof.* Necessity follows from [Equation (4)](#v12-eq-stokes-compatibility) with $`z,w\in K`$.
Assume [Equation (5)](#v12-eq-kernel-obstruction). Choose a section
$`R:I\to E`$ of $`D:E\to I`$, and write $`z=k_z+RDz`$. Define, for
$`g\in I`$,

```math
M_0(z,g)=\tfrac12 b(RDz,Rg)+b(k_z,Rg).
```

Substitution and [Equation (5)](#v12-eq-kernel-obstruction) verify compatibility on
$`I`$.  Choosing a complement of $`I`$ in $`F`$ and extending arbitrarily
produces a pairing on all of $`E\times F`$.

If $`M`$ and $`M_0`$ are compatible, their difference $`N=M-M_0`$ obeys

```math
N(z,Dw)+N(w,Dz)=0.
```

Taking $`w\in K`$ gives $`N(w,Dz)=0`$ for every $`z`$, hence $`N`$ factors
through $`E/K`$ in its first argument.  The induced bijection
$`\bar D:E/K\to I`$ identifies the two active arguments, so the induced
bilinear form is alternating. Conversely,
every alternating form in [Equation (6)](#v12-eq-abstract-gauge) cancels from the
symmetric identity.  This proves the active classification.  Finally,
compatibility evaluates the second argument only on $`I`$.  Restriction to
$`I`$ maps $`\mathcal H`$ onto $`\Lambda^2(E/K)^*`$, and its kernel
consists precisely of bilinear forms vanishing on $`I`$, equivalently
forms on $`E\times(F/I)`$.  This proves the exact sequence and the
dimension formula.

∎

**Remark 2.3.**

The Stokes-visible active gauge is intrinsically an alternating two-form on
the quotient $`E/\ker D`$, not a skew matrix until a basis and an
identification with its dual are selected.  The kernel term in
[Equation (8)](#v12-eq-full-affine-direction) is invisible to Stokes compatibility and
should not be counted as an active gauge.  This distinction becomes
essential for a space–time differential, whose image is constrained by
integrability.

**Corollary 2.4. Recovery of the polynomial spatial freedom.**

<a name="v12-cor-2012-dimension"></a>

For $`E=\mathcal{P}_N`$, $`F=\mathcal{P}_{N-1}`$, and
$`D=\,\mathrm{d}/\,\mathrm{d} z`$, the derivative is onto, $`\dim F=N`$, and the active
compatible family has dimension

```math
\frac{N(N-1)}2.
```

Thus the freedom in the polynomial Stokes construction of
[[20]](#v12-ref-moulla2012) is the spatial instance of
[Theorem 2.2](#v12-thm-affine-stokes), rather than only an analogy based on Gaussian
orthogonality.

<a name="v12-section-3"></a>

## 3. Stokes' theorem on the space–time cylinder

Let $`\Omega`$ be an oriented $`d`$-dimensional manifold with boundary
and let

```math
Z=[0,T]\times\Omega
```

carry the product orientation.  Its oriented boundary is

<a name="v12-eq-cylinder-boundary"></a>
**Equation (9).**

```math

 \partial Z=
 (\{T\}\times\Omega)-(\{0\}\times\Omega)
 \;\cup\;([0,T]\times\partial\Omega).
```

The first two components will be called the temporal caps and the last one
the lateral mantle.

**Proposition 3.1. One identity on the cylinder.**

<a name="v12-prop-continuous-cylinder"></a>

For $`\alpha\in\Omega^k(Z)`$ and
$`\beta\in\Omega^{d-k}(Z)`$,

<a name="v12-eq-continuous-stokes-cylinder"></a>
**Equation (10).**

```math

 \int_Z \,\mathrm{d}\alpha\wedge\beta
 +(-1)^k\int_Z\alpha\wedge\,\mathrm{d}\beta
 =\int_{\partial Z}\operatorname{tr}\alpha\wedge
                   \operatorname{tr}\beta.
```

The right-hand side splits according to
[Equation (9)](#v12-eq-cylinder-boundary) into a final-minus-initial cap term and a
lateral boundary term.

*Proof.* The graded Leibniz rule gives
$`\,\mathrm{d}(\alpha\wedge\beta) =\,\mathrm{d}\alpha\wedge\beta+(-1)^k\alpha\wedge\,\mathrm{d}\beta`$.
Apply Stokes' theorem and use the product orientation of
[Equation (9)](#v12-eq-cylinder-boundary).

∎

**Corollary 3.2. Port-Hamiltonian balance as a cylinder law.**

<a name="v12-cor-continuous-ph-cylinder"></a>

Suppose that a distributed port-Hamiltonian system satisfies

```math
\partial_t x=(\mathcal J-\mathcal R)e,
 \qquad e=\frac{\delta H}{\delta x},
 \qquad \mathcal R\succeq0,
```

and that the Green identity for $`\mathcal J`$ is

```math
\left\langle e,\mathcal J e\right\rangle_{\Omega}
 =\left\langle y,u\right\rangle_{\partial\Omega}.
```

Then

<a name="v12-eq-continuous-power-cylinder"></a>
**Equation (11).**

```math

 H(x(T))-H(x(0))
 =\int_0^T\left\langle y,u\right\rangle_{\partial\Omega}\,\,\mathrm{d} t
  -\int_0^T\left\langle e,\mathcal R e\right\rangle_{\Omega}\,\,\mathrm{d} t.
```

Thus stored energy on the temporal caps and supplied power through the
lateral mantle are components of one oriented boundary balance.

*Proof.* The variational chain rule and the evolution equation give

```math
\frac{\,\mathrm{d}}{\,\mathrm{d} t}H(x(t))
 =\left\langle e,(\mathcal J-\mathcal R)e\right\rangle_{\Omega}
 =\left\langle e,\mathcal J e\right\rangle_{\Omega}
  -\left\langle e,\mathcal R e\right\rangle_{\Omega}.
```

Use the Green identity and integrate in time.  This is the
port-Hamiltonian specialization of the cap–mantle split in
[Proposition 3.1](#v12-prop-continuous-cylinder).

∎

<a name="v12-section-4"></a>

## 4. Discrete cylinders and the full space–time gauge

### 4.1. The temporal–spatial total complex

**Definition 4.1. Finite graded Stokes complex.**

<a name="v12-def-graded-stokes"></a>

A finite cochain complex $`(C^\bullet,d)`$ of formal dimension $`m`$ is
equipped with a graded Stokes pairing if it has bilinear maps

```math
I_k:C^k\times C^{m-k}\to\mathbb{R}
```

and boundary maps $`B_k`$ such that, for
$`a\in C^k`$ and $`b\in C^{m-k-1}`$,

<a name="v12-eq-graded-discrete-stokes"></a>
**Equation (12).**

```math

 I_{k+1}(da,b)+(-1)^kI_k(a,db)=B_k(a,b).
```

The boundary maps include the appropriate discrete traces.  This is the
cochain form of integration by parts; no metric or constitutive operator
is included in the definition.

Let $`(C_t^\bullet,d_t)`$ and $`(C_x^\bullet,d_x)`$ be finite cochain
complexes.  Their total complex is

<a name="v12-eq-total-complex"></a>
**Equation (13).**

```math

 C_\mathrm{st}^q=\bigoplus_{r+s=q}C_t^r\otimes C_x^s,
```

with, on a homogeneous temporal cochain $`a\in C_t^r`$,

<a name="v12-eq-total-differential"></a>
**Equation (14).**

```math

 d_\mathrm{st}(a\otimes b)
 =d_ta\otimes b+(-1)^r a\otimes d_xb.
```

In matrices this is

<a name="v12-eq-total-differential-matrix"></a>
**Equation (15).**

```math

 D_\mathrm{st}=D_t\otimes I+\Gamma_t\otimes D_x,
```

where $`\Gamma_t`$ acts by $`(-1)^r`$ on temporal degree $`r`$.

**Proposition 4.2. Complex property and signed product Stokes formula.**

<a name="v12-prop-total-complex"></a>

The total differential satisfies $`d_\mathrm{st}^2=0`$.
Assume that the factor complexes have formal dimensions $`m_t,m_x`$ and
carry the pairings of [Definition 4.1](#v12-def-graded-stokes).  For homogeneous simple
tensors of complementary bidegrees, define

<a name="v12-eq-signed-tensor-pairing"></a>
**Equation (16).**

```math

 I_\mathrm{st}(a_t\otimes a_x,c_t\otimes c_x)
 =(-1)^{|a_x||c_t|}
 I_t(a_t,c_t)I_x(a_x,c_x).
```

Let $`|a_t|=r`$, $`|a_x|=s`$.  A homogeneous element of total degree
$`m_t+m_x-r-s-1`$ has the two complementary components

```math
b^{(t)}=c_t\otimes c_x,\qquad
 (|c_t|,|c_x|)=(m_t-r-1,m_x-s),
```

and

```math
b^{(x)}=e_t\otimes e_x,\qquad
 (|e_t|,|e_x|)=(m_t-r,m_x-s-1).
```

Then the total boundary form is

<a name="v12-eq-signed-product-boundary"></a>
**Equation (17).**

```math
\begin{aligned}
B_\mathrm{st}(a,b)
 &=(-1)^{s(m_t-r-1)}
   B_t(a_t,c_t)I_x(a_x,c_x)\\
 &\quad+(-1)^{m_t+s(m_t-r)}
   I_t(a_t,e_t)B_x(a_x,e_x).

\end{aligned}
```

The first term is the temporal-cap contribution and the second is the
spatial-mantle contribution, including the product-orientation sign.

*Proof.* For a homogeneous simple tensor,

```math
\begin{aligned}
d_\mathrm{st}^2(a\otimes b)
 &=d_t^2a\otimes b
   +(-1)^{r+1}d_ta\otimes d_xb
   +(-1)^r d_ta\otimes d_xb
   +a\otimes d_x^2b=0.
\end{aligned}
```

For $`b^{(t)}`$, only the two temporal derivatives can pair factorwise.
Using [Equation (16)](#v12-eq-signed-tensor-pairing), their sum is

```math
(-1)^{s(m_t-r-1)}
 \bigl[I_t(d_ta_t,c_t)+(-1)^rI_t(a_t,d_tc_t)\bigr]
 I_x(a_x,c_x),
```

which is the first term of [Equation (17)](#v12-eq-signed-product-boundary).  For
$`b^{(x)}`$, the sign $`(-1)^r`$ in the spatial part of
[Equation (14)](#v12-eq-total-differential) combines with the interchange sign in
[Equation (16)](#v12-eq-signed-tensor-pairing).  The result is

```math
(-1)^{m_t+s(m_t-r)}I_t(a_t,e_t)
 \bigl[I_x(d_xa_x,e_x)+(-1)^sI_x(a_x,d_xe_x)\bigr],
```

which gives the mantle term.  Linearity handles sums of bidegrees.  This
also displays explicitly why the two boundary pieces add while the mixed
interior derivatives cancel.

∎

The affine freedom of the genuine total complex can now be stated without
replacing it by a directional gradient.  Let $`m`$ be its formal dimension
and fix $`q\in\{0,\ldots,m-1\}`$.  Set

<a name="v12-eq-graded-double-spaces"></a>
**Equation (18).**

```math

 \mathbb E_q=C_\mathrm{st}^q\oplus C_\mathrm{st}^{m-q-1},\qquad
 \mathbb F_q=C_\mathrm{st}^{q+1}\oplus C_\mathrm{st}^{m-q},
```

and define

<a name="v12-eq-graded-double-differential"></a>
**Equation (19).**

```math

 \mathbb D_q(a,b)=(d_\mathrm{st}a,d_\mathrm{st}b).
```

For $`z=(a,b)`$ and $`\widehat z=(\widehat a,\widehat b)`$, introduce the
symmetric doubled boundary form

<a name="v12-eq-graded-double-boundary"></a>
**Equation (20).**

```math

 \mathbf b_q(z,\widehat z)
 =B_q(a,\widehat b)+B_q(\widehat a,b).
```

**Theorem 4.3. Affine classification for the genuine total complex.**

<a name="v12-thm-graded-affine"></a>

Put

```math
U_q=C_\mathrm{st}^q/\ker d_\mathrm{st}^q,
 \qquad
 V_q=C_\mathrm{st}^{m-q-1}/\ker d_\mathrm{st}^{m-q-1}.
```

The pairing

<a name="v12-eq-graded-double-pairing"></a>
**Equation (21).**

```math

 \mathbb M_q\bigl((a,b),(f,g)\bigr)
 =I_{q+1}(f,b)+(-1)^qI_q(a,g)
```

is compatible with $`\mathbb D_q`$ and $`\mathbf b_q`$.  Without a degree
restriction, the active compatible family is affine over

<a name="v12-eq-graded-unrestricted-gauge"></a>
**Equation (22).**

```math

 \Lambda^2(U_q\oplus V_q)^*.
```

If the pairing is required to preserve complementary degree—so that it
pairs $`C_\mathrm{st}^q`$ only with $`dC_\mathrm{st}^{m-q-1}`$, and
$`C_\mathrm{st}^{m-q-1}`$ only with $`dC_\mathrm{st}^{q}`$—its homogeneous
freedom is instead exactly

<a name="v12-eq-graded-active-gauge"></a>
**Equation (23).**

```math

 \boxed{\quad \operatorname{Bil}(U_q,V_q;\mathbb R)
       \simeq U_q^*\otimes V_q^*.\quad}
```

Thus cancellations between different bidegrees are incorporated through
the actual kernels of the total differential.

*Proof.* For $`z=(a,b)`$ and $`\widehat z=(\widehat a,\widehat b)`$, the left-hand
side of the symmetric compatibility identity is

```math
\begin{aligned}
&\mathbb M_q(z,\mathbb D_q\widehat z)
  +\mathbb M_q(\widehat z,\mathbb D_qz)\\
 &=I_{q+1}(d\widehat a,b)+(-1)^qI_q(\widehat a,db)
  +I_{q+1}(da,\widehat b)+(-1)^qI_q(a,d\widehat b)\\
 &=B_q(\widehat a,b)+B_q(a,\widehat b),
\end{aligned}
```

which is [Equation (20)](#v12-eq-graded-double-boundary).  Since
$`\ker\mathbb D_q=\ker d_\mathrm{st}^q\oplus \ker d_\mathrm{st}^{m-q-1}`$, [Theorem 2.2](#v12-thm-affine-stokes) gives
[Equation (22)](#v12-eq-graded-unrestricted-gauge).

For the degree-preserving subfamily, let a homogeneous correction have
off-diagonal components $`N_1(a,d\widehat b)`$ and
$`N_2(b,d\widehat a)`$.  After quotienting the two kernels, compatibility
with zero boundary requires

```math
N_1(u,\widehat v)+N_2(v,\widehat u)
 +N_1(\widehat u,v)+N_2(\widehat v,u)=0.
```

Taking $`\widehat u=0`$ and $`v=0`$ gives
$`N_2(\widehat v,u)=-N_1(u,\widehat v)`$.  Hence one arbitrary bilinear
form $`\eta\in\operatorname{Bil}(U_q,V_q;\mathbb R)`$ determines both
components through

<a name="v12-eq-graded-homogeneous-correction"></a>
**Equation (24).**

```math

 N_\eta\bigl((a,b),\mathbb D_q(\widehat a,\widehat b)\bigr)
 =\eta([a],[\widehat b])-\eta([\widehat a],[b]),
```

and every degree-preserving correction is obtained in this way.

∎

**Remark 4.4. The total and stacked kernels are different.**

<a name="v12-rem-total-versus-stacked"></a>

Take two copies of the complex
$`0\to\mathbb R\xrightarrow{1}\mathbb R\to0`$.  At total degree one,
the total differential is, up to orientation,
$`d_\mathrm{st}(a,b)=b-a`$; its kernel is the diagonal even though both
factor differentials have zero kernel.  Therefore the formula
$`\ker D_\mathrm{st}=K_t\otimes K_x`$ below is not a formula for an arbitrary
total-complex differential.  It belongs to the stacked one-bidegree map,
where the two directional derivatives are retained as separate output
components.

The preceding theorem treats the genuine graded total complex.  The
finite tensor cells used later employ a related but distinct stacked map.
Let
$`(E_a,F_a,D_a,b_a,M_a)`$, $`a\in\{t,x\}`$, be compatible Stokes data,
and let $`G_a:E_a\times E_a\to\mathbb{R}`$ be symmetric interior pairings.  Put

<a name="v12-eq-est"></a>

<a name="v12-eq-fst"></a>

<a name="v12-eq-dst-gradient"></a>
**Equations (25)–(27).**

```math
\begin{aligned}
E_\mathrm{st}&=E_t\otimes E_x,\\
 F_\mathrm{st}&=(F_t\otimes E_x)\oplus(E_t\otimes F_x),\\
 D_\mathrm{st}z&=((D_t\otimes I)z,(I\otimes D_x)z),
\end{aligned}
```

where the sign $`\Gamma_t`$ is inserted in the second component when the
chosen bidegree requires it.  Unlike [Equation (14)](#v12-eq-total-differential), the two
directional images in [Equation (27)](#v12-eq-dst-gradient) are not added.

**Theorem 4.5. Product representative and cylinder boundary.**

<a name="v12-thm-product-representative"></a>

There is a compatible product representative $`M_\mathrm{st}^0`$ defined on
simple tensors by

<a name="v12-eq-product-representative"></a>
**Equation (28).**

```math
\begin{aligned}
&M_\mathrm{st}^0
 \bigl(z_t\otimes z_x,
       (f_t\otimes e_x,e_t\otimes f_x)\bigr)\\
 &\qquad=
 M_t(z_t,f_t)G_x(z_x,e_x)
 +G_t(z_t,e_t)M_x(z_x,f_x),

\end{aligned}
```

and extended bilinearly.  Its boundary form is

<a name="v12-eq-product-boundary"></a>
**Equation (29).**

```math

 b_\mathrm{st}=b_t\otimes G_x+G_t\otimes b_x.
```

The first summand is supported on the temporal caps and the second on the
spatial mantle.

*Proof.* It suffices to test pure tensors $`z=z_t\otimes z_x`$ and
$`w=w_t\otimes w_x`$.  Substitution of [Equation (27)](#v12-eq-dst-gradient) into
[Equation (28)](#v12-eq-product-representative) gives

```math
\begin{aligned}
&\bigl[M_t(z_t,D_tw_t)+M_t(w_t,D_tz_t)\bigr]
   G_x(z_x,w_x)\\
 &\quad+G_t(z_t,w_t)
   \bigl[M_x(z_x,D_xw_x)+M_x(w_x,D_xz_x)\bigr].
\end{aligned}
```

Factor compatibility turns this expression into
$`b_t(z_t,w_t)G_x(z_x,w_x)+G_t(z_t,w_t)b_x(z_x,w_x)`$.
Pure tensors span $`E_\mathrm{st}`$, proving the result.  In the graded
version the Koszul sign in [Equation (14)](#v12-eq-total-differential) produces the
oriented sign in the mantle term.

∎

### 4.2. Classification and sector decomposition

Write $`K_a=\ker D_a`$, $`n_a=\dim E_a`$, and
$`k_a=\dim K_a`$ for $`a\in\{t,x\}`$.

**Theorem 4.6. Full active gauge of the stacked cylinder cell.**

<a name="v12-thm-space-time-gauge"></a>

For [Equation (27)](#v12-eq-dst-gradient),

<a name="v12-eq-space-time-kernel"></a>
**Equation (30).**

```math

 K_\mathrm{st}:=\ker D_\mathrm{st}=K_t\otimes K_x.
```

Consequently

<a name="v12-eq-space-time-rank"></a>
**Equation (31).**

```math

 r_\mathrm{st}:=\operatorname{rank}D_\mathrm{st}
 =n_tn_x-k_tk_x.
```

Every active compatible space–time pairing is

<a name="v12-eq-space-time-affine-family"></a>
**Equation (32).**

```math

 M_\mathrm{st}^{\omega}(z,D_\mathrm{st}w)
 =M_\mathrm{st}^0(z,D_\mathrm{st}w)+\omega([z],[w]),
 \qquad
 \omega\in\mathcal G_\mathrm{st},
```

where

<a name="v12-eq-intrinsic-space-time-gauge"></a>
**Equation (33).**

```math

 \boxed{\quad
 \mathcal G_\mathrm{st}
 =\Lambda^2(E_\mathrm{st}/K_\mathrm{st})^*.
 \quad}
```

Its dimension is

<a name="v12-eq-space-time-gauge-dimension"></a>
**Equation (34).**

```math

 \dim\mathcal G_\mathrm{st}
 =\frac{r_\mathrm{st}(r_\mathrm{st}-1)}2.
```

*Proof.* The kernel of $`D_t\otimes I`$ is $`K_t\otimes E_x`$, while the kernel
of $`I\otimes D_x`$ is $`E_t\otimes K_x`$.  Their intersection is
$`K_t\otimes K_x`$, proving [Equation (30)](#v12-eq-space-time-kernel); rank–nullity
gives [Equation (31)](#v12-eq-space-time-rank).  The product representative in
[Theorem 4.5](#v12-thm-product-representative) proves existence, and
[Theorem 2.2](#v12-thm-affine-stokes) applied to $`D_\mathrm{st}`$ gives
[Equation (32)](#v12-eq-space-time-affine-family)–[Equation (34)](#v12-eq-space-time-gauge-dimension).

∎

**Proposition 4.7. Canonical pure–mixed filtration.**

<a name="v12-prop-canonical-sector-sequence"></a>

Let

```math
Q_\mathrm{st}=(E_t\otimes E_x)/(K_t\otimes K_x),\qquad
 \overline E_t=E_t/K_t,\quad \overline E_x=E_x/K_x.
```

There is a canonical short exact sequence

<a name="v12-eq-canonical-sector-sequence"></a>
**Equation (35).**

```math

 0\longrightarrow
 (\overline E_t\otimes K_x)\oplus(K_t\otimes\overline E_x)
 \xrightarrow{\;\jmath\;}Q_\mathrm{st}
 \xrightarrow{\;\pi_{tx}\;}
 \overline E_t\otimes\overline E_x
 \longrightarrow0.
```

The subspace on the left contains the pure temporal and pure spatial
directions.  The quotient on the right is the canonical mixed sector.
A choice of complements of $`K_t`$ and $`K_x`$ supplies a natural
factorized splitting of this filtration.  General linear sections of
$`\pi_{tx}`$ need not arise from such a pair of complements.

*Proof.* The tensor product of the quotient maps
$`E_t\to\overline E_t`$ and $`E_x\to\overline E_x`$ vanishes on
$`K_t\otimes K_x`$, hence induces $`\pi_{tx}`$.  Its kernel is the image
of $`E_t\otimes K_x+K_t\otimes E_x`$ modulo
$`K_t\otimes K_x`$.  Quotienting the intersection
$`(E_t\otimes K_x)\cap(K_t\otimes E_x)=K_t\otimes K_x`$
identifies this kernel with the direct sum on the left of
[Equation (35)](#v12-eq-canonical-sector-sequence).  Surjectivity follows from that of
the two factor quotient maps.

∎

**Corollary 4.8. Temporal, spatial, and mixed sectors.**

<a name="v12-cor-gauge-sectors"></a>

Choose complements $`E_t=K_t\oplus U_t`$ and
$`E_x=K_x\oplus U_x`$.  Then

<a name="v12-eq-quotient-sector-split"></a>
**Equation (36).**

```math

 E_\mathrm{st}/K_\mathrm{st}
 \simeq Q_t\oplus Q_x\oplus Q_{tx},
```

where

```math
Q_t=U_t\otimes K_x,
 \qquad Q_x=K_t\otimes U_x,
 \qquad Q_{tx}=U_t\otimes U_x
 \simeq\overline E_t\otimes\overline E_x.
```

Accordingly,

<a name="v12-eq-wedge-sector-split"></a>
**Equation (37).**

```math

 \mathcal G_\mathrm{st}
 \simeq
 \bigoplus_{a\in\{t,x,tx\}}\Lambda^2Q_a^*
 \;\oplus\!
 \bigoplus_{a<b}Q_a^*\otimes Q_b^*.
```

The summands contain pure temporal, pure spatial, mixed, and cross-sector
gauges.  The total space $`\mathcal G_\mathrm{st}`$, the pure subspace, and
the mixed quotient in [Equation (35)](#v12-eq-canonical-sector-sequence) are intrinsic;
only the displayed direct-sum splitting depends on the chosen complements.

*Proof.* Expanding
$`(K_t\oplus U_t)\otimes(K_x\oplus U_x)`$ and quotienting by
$`K_t\otimes K_x`$ gives [Equation (36)](#v12-eq-quotient-sector-split).  For a direct
sum $`Q=Q_1\oplus Q_2\oplus Q_3`$, the standard exterior-square identity
is
$`\Lambda^2Q^*=\bigoplus_i\Lambda^2Q_i^* \oplus\bigoplus_{i<j}Q_i^*\otimes Q_j^*`$, which proves [Equation (37)](#v12-eq-wedge-sector-split).

∎

**Proposition 4.9. Pure-axis embeddings and boundary reduction.**

<a name="v12-prop-boundary-reductions"></a>

For nonzero $`k_t\in K_t`$ and $`k_x\in K_x`$, the maps

```math
\begin{aligned}
\iota_t:E_t/K_t&\longrightarrow E_\mathrm{st}/K_\mathrm{st},
 &[z_t]&\longmapsto[z_t\otimes k_x],\\
 \iota_x:E_x/K_x&\longrightarrow E_\mathrm{st}/K_\mathrm{st},
 &[z_x]&\longmapsto[k_t\otimes z_x]
\end{aligned}
```

are injective.  Pullback by $`\iota_t`$ and $`\iota_x`$ sends a
space–time gauge to its pure temporal and pure spatial restrictions.  In
particular, the time-constant pullback contains the polynomial spatial
gauges of [[20]](#v12-ref-moulla2012).  A cap-to-cap temporal method on the full
spatial state space is obtained only after eliminating the interior
spatial flows; it generally receives contributions from
$`Q_t\oplus Q_{tx}`$.  The next sections classify that reduced boundary
gauge.  A gauge supported entirely on $`Q_{tx}`$ can be nonzero while both
pure-axis restrictions vanish.

*Proof.* If $`z_t\otimes k_x\in K_t\otimes K_x`$ and $`k_x\ne0`$, then
$`z_t\in K_t`$; hence $`\iota_t`$ is injective.  The spatial statement is
identical.  Alternating forms pull back functorially.  Finally a form
supported on $`Q_{tx}`$ vanishes on the images of both pure embeddings.

∎

### 4.3. Regular volume-to-cap reduction

The pullbacks in [Proposition 4.9](#v12-prop-boundary-reductions) restrict a gauge along
pure axes.  A numerical time step requires a different operation: the
interior equations must first be solved and the resulting active volume
field restricted to its cap data.  This operation depends on the Dirac
relation, the constitutive closure, the boundary causality, and the step
size; it is not determined by the Stokes pairing alone.

Let

```math
Q_\mathrm{st}=E_\mathrm{st}/K_\mathrm{st}
```

be the active volume space and let $`Q_\partial`$ be the chosen active cap
space.  After fixing the four ingredients above, call the reduction
*regular* if the interior equations determine a linear reconstruction

<a name="v12-eq-regular-lifting"></a>
**Equation (38).**

```math

 L_h:Q_\partial\longrightarrow Q_\mathrm{st}.
```

The notation permits $`L_h`$ to depend on the selected constitutive law
and on $`h`$, but it is held fixed while the following affine gauge
reduction is classified.

**Theorem 4.10. Exact visibility sequence for volume-to-cap gauges.**

<a name="v12-thm-volume-cap-reduction"></a>

For a regular reconstruction $`L_h`$, the homogeneous volume-to-cap reduction is
the pullback map

<a name="v12-eq-volume-cap-pullback"></a>
**Equation (39).**

```math

 \mathcal R_{L_h}:\Lambda^2Q_\mathrm{st}^*
 \longrightarrow\Lambda^2Q_\partial^*,
 \qquad
 \mathcal R_{L_h}(\omega)=L_h^*\omega,
```

where

```math
(L_h^*\omega)(\xi,\zeta)
 =\omega(L_h\xi,L_h\zeta).
```

Put

```math
W_h=\operatorname{im}L_h,\qquad
 K_{\partial,h}=\ker L_h,\qquad r_h=\operatorname{rank}L_h.
```

Then

<a name="v12-eq-volume-cap-image"></a>
**Equation (40).**

```math

 \operatorname{im}\mathcal R_{L_h}
 =
 \left\{\eta\in\Lambda^2Q_\partial^*:
 \iota_k\eta=0\ \text{for every }k\in K_{\partial,h}\right\}
 \simeq\Lambda^2(Q_\partial/K_{\partial,h})^*,
```

and

<a name="v12-eq-volume-cap-kernel"></a>
**Equation (41).**

```math

 \ker\mathcal R_{L_h}
 =
 \left\{\omega\in\Lambda^2Q_\mathrm{st}^*:
 \omega|_{W_h\times W_h}=0\right\}.
```

Equivalently, restriction to $`W_h`$ gives the short exact sequence

<a name="v12-eq-volume-cap-exact-sequence"></a>
**Equation (42).**

```math

 0\longrightarrow\ker\mathcal R_{L_h}
 \longrightarrow\Lambda^2Q_\mathrm{st}^*
 \xrightarrow{\;\mathcal R_{L_h}\;}
 \Lambda^2(Q_\partial/K_{\partial,h})^*
 \longrightarrow0.
```

Consequently,

<a name="v12-eq-volume-cap-image-dimension"></a>

<a name="v12-eq-volume-cap-kernel-dimension"></a>
**Equations (43)–(44).**

```math
\begin{aligned}
\dim\operatorname{im}\mathcal R_{L_h}
 &=\frac{r_h(r_h-1)}2,\\
 \dim\ker\mathcal R_{L_h}
 &=\binom{\dim Q_\mathrm{st}}2-\binom{r_h}2.

\end{aligned}
```

If $`L_h`$ is injective, every cap gauge has a volume extension.  If
$`L_h`$ is surjective, no nonzero volume gauge is invisible on all cap
data.  In particular, an isomorphism $`L_h`$ induces an isomorphism of
gauge spaces.

*Proof.* If $`k\in K_{\partial,h}`$, then
$`(L_h^*\omega)(k,\xi)=\omega(0,L_h\xi)=0`$, proving the inclusion from
left to right in [Equation (40)](#v12-eq-volume-cap-image).  Conversely, every
$`\eta`$ satisfying the annihilation condition descends uniquely to an
alternating form on $`Q_\partial/K_{\partial,h}`$.  The induced map

```math
\overline L_h:Q_\partial/K_{\partial,h}\longrightarrow W_h
```

is an isomorphism.  Transport $`\eta`$ to $`W_h`$ by this isomorphism and
extend the resulting alternating form from $`W_h`$ to $`Q_\mathrm{st}`$
using any linear complement of $`W_h`$.  Its pullback is $`\eta`$, so
[Equation (40)](#v12-eq-volume-cap-image) and surjectivity in
[Equation (42)](#v12-eq-volume-cap-exact-sequence) follow.

The pullback vanishes exactly when $`\omega`$ vanishes on every pair of
vectors in $`W_h`$, which is [Equation (41)](#v12-eq-volume-cap-kernel).  Since
$`\dim W_h=r_h`$, the two dimension formulas follow.  Injectivity of
$`L_h`$ means $`K_{\partial,h}=0`$, while surjectivity means
$`W_h=Q_\mathrm{st}`$; the final claims are immediate.

∎

**Corollary 4.11. Affine reduction and gauge-invariant power.**

<a name="v12-cor-affine-volume-cap"></a>

Assume that a reference compatible volume pairing reduces, under the same
fixed reconstruction, to a cap pairing $`m_{\partial,0}`$.  Then the full reduced
family is

<a name="v12-eq-reduced-affine-family"></a>
**Equation (45).**

```math

 m_{\partial,\omega}
 =m_{\partial,0}+L_h^*\omega.
```

For every cap increment $`\xi`$,

```math
(L_h^*\omega)(\xi,\xi)=0.
```

Thus the visible gauge changes the cap method but not its contracted power
identity.

**Corollary 4.12. Intrinsic reduction of the mixed quotient.**

<a name="v12-cor-mixed-visibility"></a>

Let
$`Q_{tx}=\overline E_t\otimes\overline E_x`$ and let
$`\pi_{tx}:Q_\mathrm{st}\to Q_{tx}`$ be the canonical quotient map in
[Equation (35)](#v12-eq-canonical-sector-sequence).  Pullback embeds the canonical mixed
gauges as $`\pi_{tx}^*\Lambda^2Q_{tx}^*\subset\Lambda^2Q_\mathrm{st}^*`$, and

<a name="v12-eq-intrinsic-mixed-pullback"></a>
**Equation (46).**

```math

 \boxed{\quad
 \mathcal R_{L_h}\circ\pi_{tx}^*
 =(\pi_{tx}L_h)^*.
 \quad}
```

If $`r_{tx,h}=\operatorname{rank}(\pi_{tx}L_h)`$, then the visible
canonical mixed family has dimension

<a name="v12-eq-mixed-visibility-dimension"></a>
**Equation (47).**

```math

 \dim\operatorname{im}(\mathcal R_{L_h}\circ\pi_{tx}^*)
 =\binom{r_{tx,h}}2.
```

Its kernel consists exactly of the forms on $`Q_{tx}`$ that vanish on
$`\operatorname{im}(\pi_{tx}L_h)\times \operatorname{im}(\pi_{tx}L_h)`$, and has dimension

```math
\binom{\dim Q_{tx}}2-\binom{r_{tx,h}}2.
```

In particular, there exists a nonzero canonical mixed gauge whose cap
pullback is nonzero if and only if

<a name="v12-eq-mixed-visibility-rank"></a>
**Equation (48).**

```math

 \operatorname{rank}(\pi_{tx}L_h)\ge2.
```

A prescribed nonzero mixed form may nevertheless belong to the kernel.

*Proof.* Functoriality gives
$`L_h^*\pi_{tx}^*=(\pi_{tx}L_h)^*`$.  Apply
[Theorem 4.10](#v12-thm-volume-cap-reduction) to the single map $`\pi_{tx}L_h`$.  The
image, kernel, and dimension statements follow verbatim and require no
choice of tensor complements.

∎

**Proposition 4.13. Schur realization of the reconstruction.**

<a name="v12-prop-schur-lifting"></a>

Suppose a linear Dirac–constitutive cell has been partitioned into
interior variables $`z_i`$ and cap variables $`\xi`$, and its interior
equations are

<a name="v12-eq-interior-block-system"></a>
**Equation (49).**

```math

 K_{ii}(h)z_i+K_{i\partial}(h)\xi=0.
```

If $`K_{ii}(h)`$ is invertible, the regular reconstruction is

<a name="v12-eq-schur-lifting"></a>
**Equation (50).**

```math

 L_h\xi=
 \begin{bmatrix}
  -K_{ii}(h)^{-1}K_{i\partial}(h)\xi\\[1mm]\xi
 \end{bmatrix}.
```

In coordinates, a volume gauge matrix $`\Omega_\mathrm{st}^\top =-\Omega_\mathrm{st}`$ reduces to

<a name="v12-eq-schur-gauge"></a>
**Equation (51).**

```math

 \Omega_\partial=L_h^\top\Omega_\mathrm{st}L_h.
```

*Proof.* [Equation (50)](#v12-eq-schur-lifting) is the unique solution of
[Equation (49)](#v12-eq-interior-block-system); [Equation (51)](#v12-eq-schur-gauge) is the coordinate
form of [Equation (39)](#v12-eq-volume-cap-pullback).

∎

**Theorem 4.14. Trace-generated cell and exact raw realization.**

<a name="v12-thm-trace-generated-condensation"></a>

Let $`X`$ be the active cap state space and $`Z_h`$ a raw cell space.  Let
$`T_h:Z_h\to X`$ be onto, let $`Q:X\to X^*`$ be symmetric, and let
$`J:X^*\to X`$ be skew.  For an alternating cap map
$`S_h:X\to X^*`$, put

<a name="v12-eq-trace-generated-volume-form"></a>
**Equation (52).**

```math

 \Omega_h=T_h^*S_hT_h:Z_h\longrightarrow Z_h^*.
```

It descends to $`\overline\Omega_h`$ on
$`\overline Z_h=Z_h/\ker T_h`$.  Given $`x_0`$, consider the weak cell
equations for $`([z],e,x_1)\in\overline Z_h\times X^*\times X`$:

<a name="v12-eq-raw-cell-trace"></a>

<a name="v12-eq-raw-cell-dirac"></a>

<a name="v12-eq-raw-cell-constitutive"></a>
**Equations (53)–(55).**

```math
\begin{aligned}
\overline T_h[z]&=x_1-x_0,\\
 x_1-x_0&=hJe,\\
 \overline T_h^*\bigl(e-Q\overline x\bigr)
 +\overline\Omega_h[z]&=0,
 \qquad \overline x=\tfrac12(x_0+x_1).

\end{aligned}
```

Then the cap closure encoded by $`\Omega_h`$ is recovered exactly by
elimination:

<a name="v12-eq-condensed-cap-closure"></a>
**Equation (56).**

```math

 e=Q\overline x-S_h(x_1-x_0)
```

and its exact Schur complement is

<a name="v12-eq-trace-generated-cap-map"></a>
**Equation (57).**

```math

 \boxed{\quad
 \bigl(I+hJS_h-\tfrac h2JQ\bigr)x_1
 =\bigl(I+hJS_h+\tfrac h2JQ\bigr)x_0.
 \quad}
```

[Equation (57)](#v12-eq-trace-generated-cap-map) is always a linear relation.
If

<a name="v12-eq-trace-generated-invertibility"></a>
**Equation (58).**

```math

 I+hJS_h-\tfrac h2JQ
 \quad\text{is invertible},
```

it defines a unique map from $`x_0`$ to $`x_1`$.  The physical variables
$`(e,x_1)`$ are independent of the representative
of $`[z]`$ and of the right inverse used to coordinatize the quotient.

In Euclidean coordinates, let the columns of $`N_h`$ form an orthonormal
basis of $`\ker T_h`$.  A square raw assembly, with a multiplier fixing
the null component, is

<a name="v12-eq-raw-cell-block-system"></a>
**Equation (59).**

```math

 \begin{bmatrix}
 T_h&0&-I&0\\
 0&-hJ&I&0\\
 \Omega_h&T_h^*&-\tfrac12T_h^*Q&N_h\\
 N_h^*&0&0&0
 \end{bmatrix}
 \begin{bmatrix}z\\e\\x_1\\\lambda\end{bmatrix}
 =
 \begin{bmatrix}-I\\I\\\tfrac12T_h^*Q\\0\end{bmatrix}x_0.
```

Under [Equation (58)](#v12-eq-trace-generated-invertibility), this square system is
invertible and its Schur map from $`x_0`$ to $`x_1`$ is exactly
[Equation (57)](#v12-eq-trace-generated-cap-map).  Thus the prescribed cap gauge is not
attached as an additional reduced equation: it is encoded by the
trace-generated volume form and recovered exactly from the raw cell.

*Proof.* Because $`T_h`$ is onto, the induced map
$`\overline T_h:\overline Z_h\to X`$ is an isomorphism and
$`\overline T_h^*`$ is injective.  Using
[Equation (52)](#v12-eq-trace-generated-volume-form), [Equation (53)](#v12-eq-raw-cell-trace), the third cell
equation becomes

```math
\overline T_h^*
 \bigl(e-Q\overline x+S_h(x_1-x_0)\bigr)=0.
```

Injectivity gives [Equation (56)](#v12-eq-condensed-cap-closure).  Substitution into
[Equation (54)](#v12-eq-raw-cell-dirac) gives

```math
x_1-x_0=hJ\bigl(Q\overline x-S_h(x_1-x_0)\bigr),
```

which rearranges to [Equation (57)](#v12-eq-trace-generated-cap-map).  If
$`R_h:X\to Z_h`$ is any right inverse, then

```math
R_h^*\Omega_hR_h
 =(T_hR_h)^*S_h(T_hR_h)=S_h.
```

Two right inverses differ by a map into $`\ker T_h`$, on which
$`\Omega_h`$ vanishes.  Finally, the last row of
[Equation (59)](#v12-eq-raw-cell-block-system) selects the orthogonal representative;
the third row is the unquotiented weak equation.  Eliminating
$`(z,e,\lambda)`$ reproduces the preceding calculation.  Under
[Equation (58)](#v12-eq-trace-generated-invertibility), a homogeneous right-hand side
first gives $`x_1=e=0`$ from the condensed equations, then $`z=0`$ from
$`T_hz=0`$ and $`N_h^*z=0`$, and finally $`\lambda=0`$.  The square raw
system is therefore invertible.

∎

<a name="v12-fig-dynamic-condensation"></a>

```mermaid
flowchart TD
    R["Raw cylinder: caps and mantle"] --> T["Trace assembly: local cell"]
    T --> C["Exact Schur condensation: cap propagator"]
    R -->|"Quotient and right inverse"| C
```

The raw cylinder $`Z_h`$ is considered modulo $`\ker T_h`$. Trace assembly $`T_h`$ defines the local cell $`(z,e,x_1,\lambda)`$ with $`\Omega_h=T_h^*S_hT_h`$, whose condensed cap propagator is $`M_{S,h}`$. The lower route expresses $`R_h^*\Omega_hR_h=S_h`$.

**Figure 2.** Exact raw-cell realization of a prescribed cap gauge.  The lower
identity is independent of the chosen section; the upper path shows that
the gauge encoded in $`\Omega_h=T_h^*S_hT_h`$ is recovered exactly by
elimination.  No claim of selecting $`S_h`$ from independent interior
dynamics is made here.

**Remark 4.15. Singular interior problems.**

If $`K_{ii}(h)`$ is singular, the interior solution may be nonunique or
subject to compatibility constraints.  The reduction is then a linear
relation rather than a map $`L_h`$, and composition of Dirac relations or
presymplectic reduction is required.  Likewise, if changing the volume
gauge changes the interior solve, one obtains $`L_h^\omega`$ and the
global map $`\omega\mapsto(L_h^\omega)^*\omega`$ need not be affine.
[Theorem 4.10](#v12-thm-volume-cap-reduction) classifies the exact gauge freedom within
each fixed regular reduction chart; it does not hide these hypotheses.

**Proposition 4.16. Tangent reduction when the reconstruction depends on the gauge.**

<a name="v12-prop-gauge-dependent-tangent"></a>

Suppose $`\omega\mapsto L_h^\omega`$ is differentiable and define

```math
\mathcal F_h(\omega)=(L_h^\omega)^*\omega.
```

For $`\xi,\zeta\in Q_\partial`$,

<a name="v12-eq-gauge-dependent-differential"></a>
**Equation (60).**

```math
\begin{aligned}
D\mathcal F_h(\omega)[\delta\omega](\xi,\zeta)
 &=\delta\omega(L_h^\omega\xi,L_h^\omega\zeta)\\
 &\quad+\omega(DL_h^\omega[\delta\omega]\xi,L_h^\omega\zeta)\\
 &\quad+\omega(L_h^\omega\xi,DL_h^\omega[\delta\omega]\zeta).

\end{aligned}
```

At the zero gauge,

<a name="v12-eq-tangent-pullback"></a>
**Equation (61).**

```math

 D\mathcal F_h(0)[\delta\omega]=L_h^{0*}\delta\omega.
```

Thus the fixed-reconstruction pullback theorem is also the exact tangent
model of a gauge-dependent nonlinear reduction.  For the trace-generated
family in [Theorem 4.14](#v12-thm-trace-generated-condensation), the stronger finite
identity holds because the quotient trace itself is gauge independent.

*Proof.* Differentiate
$`\omega(L_h^\omega\xi,L_h^\omega\zeta)`$ by the product rule.  At
$`\omega=0`$, the two terms containing the undifferentiated form vanish.

∎

<a name="v12-section-5"></a>

## 5. The complete two-point temporal gauge

We now pass to the cap-to-cap problem obtained after eliminating interior
spatial flows, while keeping the spatial state space arbitrary.  In the
language of [Theorem 4.10](#v12-thm-volume-cap-reduction), its gauge belongs to the image
of $`\mathcal R_{L_h}`$ and can therefore contain the trace of mixed
space–time sectors precisely under [Equation (48)](#v12-eq-mixed-visibility-rank).  This
is also the smallest setting in which the distinction between
gauge-invariant power and gauge-dependent symplecticity is visible.

### 5.1. Temporal Stokes data

Let $`V=\mathbb{R}^n`$, let $`Q=Q^\top`$, and define

```math
H(x)=\tfrac12x^\top Qx,\qquad
 z=(x_0,x_1)\in E_t=V\oplus V,\qquad
 D_tz=x_1-x_0=:\Delta x.
```

The polarized temporal boundary form is

<a name="v12-eq-temporal-boundary"></a>
**Equation (62).**

```math

 b_Q(z,w)=x_1^\top Qy_1-x_0^\top Qy_0,
 \qquad w=(y_0,y_1).
```

Its restriction to $`K_t=\ker D_t=\{(x,x):x\in V\}`$ vanishes. Hence
[Theorem 2.2](#v12-thm-affine-stokes) applies.

A linear two-point effort is a linear map
$`\bar e:E_t\to V^*`$. It defines

```math
M_{\bar e}(z,v)=\bar e(z)[v].
```

Taking $`w=z`$ in [Equation (4)](#v12-eq-stokes-compatibility) shows that compatibility is
exactly the polarized form of

<a name="v12-eq-discrete-chain-rule"></a>
**Equation (63).**

```math

 \bar e(x_0,x_1)[\Delta x]=H(x_1)-H(x_0).
```

**Theorem 5.1. Complete linear temporal gauge.**

<a name="v12-thm-two-point-gauge"></a>

The linear efforts satisfying [Equation (63)](#v12-eq-discrete-chain-rule) for all
$`x_0,x_1`$ are precisely

<a name="v12-eq-effort-intrinsic"></a>
**Equation (64).**

```math

 \bar e_{\mathsf S}(x_0,x_1)
 =Q\bar x+\iota_{\Delta x}\mathsf S,
 \qquad
 \bar x=\tfrac12(x_0+x_1),\quad
 \mathsf S\in\Lambda^2V^*.
```

In coordinates, if $`S^\top=-S`$ represents $`\mathsf S`$, then

<a name="v12-eq-effort-matrix"></a>
**Equation (65).**

```math

 \boxed{\ \bar e_S=Q\bar x-S\Delta x\ }.
```

The affine family has dimension $`n(n-1)/2`$.

*Proof.* The midpoint effort $`Q\bar x`$ is compatible because

```math
(Q\bar x)^\top\Delta x
 =\tfrac12(x_1+x_0)^\top Q(x_1-x_0)
 =H(x_1)-H(x_0).
```

Moreover, $`E_t/K_t\simeq V`$ through $`D_t`$. Applying
[Theorem 2.2](#v12-thm-affine-stokes), the difference from the midpoint pairing is an
arbitrary $`\mathsf S\in\Lambda^2V^*`$ evaluated on
$`(\Delta x,v)`$. Since
$`\mathsf S(\Delta x,v)=(-S\Delta x)^\top v`$, this is
[Equation (65)](#v12-eq-effort-matrix). Alternation makes its contraction with
$`\Delta x`$ vanish.

∎

### 5.2. Gauge-invariant power balance

Consider the port-Hamiltonian step

<a name="v12-eq-ph-step"></a>
**Equation (66).**

```math

 \frac{\Delta x}{h}
 =\bigl(J-R\bigr)\bar e_S+G\bar u,
 \qquad
 \bar y=G^\top\bar e_S,
```

where $`J^\top=-J`$ and $`R^\top=R\succeq0`$.

**Corollary 5.2. Exact gauge-invariant power law.**

<a name="v12-cor-power"></a>

Every $`S^\top=-S`$ gives the exact identity

<a name="v12-eq-power-law"></a>
**Equation (67).**

```math

 H(x_1)-H(x_0)
 =h\,\bar y^\top\bar u-h\,\bar e_S^\top R\bar e_S.
```

For $`R=0`$, the discrete boundary supply is exact. For $`R\succeq0`$, the
step is passive.

*Proof.* Pair [Equation (66)](#v12-eq-ph-step) with $`\bar e_S`$, use
[Equation (63)](#v12-eq-discrete-chain-rule), skewness of $`J`$, and the definition of
$`\bar y`$.

∎

### 5.3. Constitutively shifted Cayley coordinates

The classical Cayley transform gives rational coordinates between a Lie
algebra and an open cell of its matrix group
[[19]](#v12-ref-kostant2001), [[24]](#v12-ref-maslovskaya2024).  The point here is not to rename that
transform, but to identify the *constitutively shifted* variable
$`I+hJS`$ selected by the discrete chain rule.  At fixed $`h`$, the cap
gauges thereby parametrize the energy-orthogonal maps having no
eigenvalue $`+1`$.  This convention differs from the common Cayley cell
that excludes $`-1`$: the group identity is not in the present cell, and
the zero gauge represents the midpoint map rather than the identity.

**Theorem 5.3. Gauge coordinate on the no-$`+1`$ Cayley cell of $`O(Q)`$.**

<a name="v12-thm-gauge-cayley-chart"></a>

Assume that $`J^\top=-J`$ and $`Q=Q^\top`$ are nonsingular, put
$`A=JQ`$, $`\Omega=J^{-1}`$, and fix $`h\ne0`$.  For an alternating
matrix $`S`$, define

<a name="v12-eq-cap-cayley-blocks"></a>
**Equation (68).**

```math

 B_\pm(S)=I+hJS\mathbin{\pm}\frac h2A,
 \qquad M_S=B_-(S)^{-1}B_+(S),
```

whenever $`B_-(S)`$ is invertible.  Then

<a name="v12-eq-gauge-cayley-bijection"></a>
**Equation (69).**

```math

 \left\{S:S^\top=-S,\ B_-(S)\text{ invertible}\right\}
 \xrightarrow[\hspace{12mm}]{\ S\mapsto M_S\ }
 \left\{M:M^\top QM=Q,\ M-I\text{ invertible}\right\}
```

is a bijection, denoted by $`\Phi_h`$.  Its inverse is

<a name="v12-eq-energy-cayley-ratio"></a>

<a name="v12-eq-inverse-gauge-cayley"></a>
**Equations (70)–(71).**

```math
\begin{aligned}
R_M&=(M+I)(M-I)^{-1},\\
 C_M&=\frac h2A R_M,\qquad
 S_M=\frac1h\Omega(C_M-I).
\end{aligned}
```

*Proof.* The step equations are
$`\Delta x=hJ(Q\bar x-S\Delta x)`$, so the discrete chain rule gives
$`x_1^\top Qx_1=x_0^\top Qx_0`$.  Hence every well-defined $`M_S`$ is
$`Q`$-orthogonal.  Moreover,

```math
M_S-I=B_-(S)^{-1}hA,
```

which is invertible because $`A`$ is.

Conversely, let $`M^\top QM=Q`$ and suppose $`M-I`$ is invertible.
The Cayley ratio in [Equation (70)](#v12-eq-energy-cayley-ratio) is $`Q`$-skew:

<a name="v12-eq-q-skew-cayley-ratio"></a>
**Equation (72).**

```math

 R_M^\top Q+QR_M=0.
```

Indeed, multiplying the left-hand side by $`(M-I)^\top`$ and $`M-I`$
reduces it to $`2(M^\top QM-Q)`$.  Since
$`A^\top\Omega=-Q`$ and $`\Omega A=Q`$,

```math
C_M^\top\Omega=-\frac h2R_M^\top Q
 =\frac h2QR_M=\Omega C_M.
```

It follows that
$`S_M^\top=-S_M`$.  Finally,

```math
C_M-\frac h2A=hA(M-I)^{-1},\qquad
 C_M+\frac h2A=hAM(M-I)^{-1},
```

so $`B_-(S_M)`$ is invertible and
$`B_-(S_M)^{-1}B_+(S_M)=M`$.  Solving
$`B_-(S)M_S=B_+(S)`$ for $`I+hJS`$ gives
$`C_{M_S}`$, proving that the two constructions are inverse.

∎

**Corollary 5.4. The exact-flow gauge.**

<a name="v12-cor-exact-flow-gauge"></a>

Where $`e^{hA}-I`$ is invertible, the exact Hamiltonian flow is represented
by

<a name="v12-eq-exact-flow-gauge"></a>
**Equation (73).**

```math

 S_h^\mathrm{exact}
 =\frac1h\Omega
 \left[
  \frac{hA}{2}\coth\!\left(\frac{hA}{2}\right)-I
 \right].
```

The matrix function in brackets is regular at the origin and

<a name="v12-eq-exact-flow-gauge-series"></a>
**Equation (74).**

```math

 S_h^\mathrm{exact}
 =\frac h{12}\Omega A^2
  -\frac{h^3}{720}\Omega A^4+O(h^5).
```

*Proof.* For $`M=e^{hA}`$, functional calculus gives
$`(M+I)(M-I)^{-1}=\coth(hA/2)`$.  Apply
[Equation (71)](#v12-eq-inverse-gauge-cayley); the series follows from the Bernoulli
expansion of $`(z/2)\coth(z/2)`$.

∎

**Remark 5.5. Stable evaluation and chart conditioning.**

<a name="v12-rem-cayley-conditioning"></a>

The inverse formula [Equation (70)](#v12-eq-energy-cayley-ratio) is an exact algebraic
identity, but evaluating it from a rounded matrix $`M`$ is ill-conditioned
as $`M\to I`$: both the subtraction $`M-I`$ and its inverse amplify
perturbations.  This is intrinsic to recovering a finite gauge from the
chosen group coordinate near the excluded identity.  It is not necessary
when the map is known analytically.  For the exact flow one evaluates
$`(z/2)\coth(z/2)-1`$ by its Bernoulli series or a cancellation-free
matrix-function algorithm.  For Gauss–$`s`$, one evaluates
[Equation (99)](#v12-eq-stage-induced-gauss-gauge) directly from the even and odd Padé
polynomials, without forming either $`R_s-I`$ or the inverse Cayley ratio.
For the exact flow, the chart has a genuine resonance whenever
$`e^{hA}`$ has eigenvalue $`+1`$, equivalently when some spectral value
satisfies $`h\lambda\in2\pi i\mathbb Z`$.  The local Bernoulli series has
radius $`2\pi`$, set by the nearest nonzero poles of
$`(z/2)\coth(z/2)`$; away from the local disk one should use a
Schur–Parlett or rational matrix-function evaluation and test the
distance to these resonances.  Since $`S_h^\mathrm{exact}=O(h)`$, a rounded
evaluation of $`M-I`$ generically gives a *relative* inverse-chart
error of size $`O(\epsilon_\mathrm{mach}/h^3)`$.  The scaled values in
[Table 8](#v12-tab-cayley-conditioning) confirm this law.
The numerical tests in [Table 8](#v12-tab-cayley-conditioning) separate this
stable forward evaluation from the condition-sensitive inverse chart.

### 5.4. Canonical selectors

**Proposition 5.6. Symmetry, localization, and order.**

<a name="v12-prop-selectors"></a>

Within [Equation (65)](#v12-eq-effort-matrix):

1. endpoint symmetry
     $`\bar e(x_1,x_0)=\bar e(x_0,x_1)`$ holds if and only if $`S=0`$;

1. dependence on the midpoint alone holds if and only if $`S=0`$;

1. if $`Q`$ is nonsingular, localization as
     $`\bar e=Q((1-\theta)x_0+\theta x_1)`$, with a fixed
     $`\theta`$, holds if and only if $`\theta=\tfrac12`$ and $`S=0`$;

1. for the closed equation $`\dot x=JQx`$, if both $`J`$ and $`Q`$
     are nonsingular, the associated one-step map is second order if and only
     if $`S=0`$. Every nonzero fixed $`S`$ is first order.

*Proof.* Interchanging the endpoints changes $`\Delta x`$ to $`-\Delta x`$, proving
the first two claims. Point localization would require
$`-S=(\theta-\tfrac12)Q`$. Its left-hand side is skew and its right-hand
side symmetric, so both vanish.

For the order statement put $`A=JQ`$. The closed method can be written

```math
\Delta x=hF_S\bar x,\qquad
 F_S=(I+hJS)^{-1}A.
```

Its map is the Cayley transform in
[Equation (68)](#v12-eq-cap-cayley-blocks), equivalently
$`M_S=(I-hF_S/2)^{-1}(I+hF_S/2)`$.
Expansion gives

<a name="v12-eq-order-expansion"></a>
**Equation (75).**

```math

 M_S=I+hA+h^2\left(\tfrac12A^2-JSA\right)+O(h^3).
```

The exact flow has coefficient $`A^2/2`$. Under the stated
invertibility assumptions, $`JSA=0`$ implies $`S=0`$.

∎

**Remark 5.7. Necessary hypotheses.**

If $`Q`$ is singular, a nonzero $`S`$ may satisfy $`SA=0`$, so the
unqualified claim that every nonzero gauge is first order is false.
The hypotheses in [Proposition 5.6](#v12-prop-selectors) are therefore substantive.

<a name="v12-section-6"></a>

## 6. The exact symplectic gauge locus

Assume in this section that $`J`$ is invertible, put
$`\Omega=J^{-1}`$, and keep $`Q=Q^\top`$. The continuous matrix
$`A=JQ`$ is Hamiltonian:

```math
A^\top\Omega+\Omega A=0.
```

**Theorem 6.1. Symplectic-locus criterion.**

<a name="v12-thm-symplectic-locus"></a>

Assume that $`I+hJS`$ and $`B_-(S)`$ are invertible.  Then the map $`M_S`$
is $`\Omega`$-symplectic if and only if

<a name="v12-eq-locus"></a>
**Equation (76).**

```math

 \boxed{\ SJQ=QJS\ }
 \qquad\Longleftrightarrow\qquad
 A^\top S+SA=0.
```

Thus the symplectic gauges form the linear space

```math
\mathcal{L}=\{S\in\Lambda^2V^*:A^\top S+SA=0\}
```

of alternating two-forms invariant under the continuous flow.

*Proof.* A Cayley transform of $`F_S`$ is $`\Omega`$-symplectic exactly when
$`F_S^\top\Omega+\Omega F_S=0`$. Since

```math
\Omega F_S=(I+hSJ)^{-1}Q,
```

this condition is the symmetry of $`(I+hSJ)^{-1}Q`$. Equivalently,

```math
(I+hSJ)^{-1}Q=Q(I+hJS)^{-1}.
```

Multiplication on the left and right yields $`SJQ=QJS`$. Finally,
$`A^\top=-QJ`$, so this is equivalent to $`A^\top S+SA=0`$.

∎

**Corollary 6.2. The symplectic intersection of the Cayley chart.**

<a name="v12-cor-cayley-symplectic-intersection"></a>

Under the hypotheses of [Theorem 5.3](#v12-thm-gauge-cayley-chart), the bijection
[Equation (69)](#v12-eq-gauge-cayley-bijection) restricts to

```math
\mathcal{L}\cap\operatorname{dom}\Phi_h
 \longleftrightarrow
 \left\{M\in O(Q)\cap Sp(\Omega):M-I\text{ invertible}\right\}.
```

Thus the centralizer equation $`A^\top S+SA=0`$ is precisely the
gauge-coordinate description of the intersection between energy
orthogonality and symplecticity.

**Theorem 6.3. Exact symplectic defect.**

<a name="v12-thm-defect"></a>

Assume $`I+hJS`$ is invertible, put
$`F_S=(I+hJS)^{-1}A`$, and let $`D_S=I-\tfrac h2F_S`$. Then

<a name="v12-eq-defect-cayley"></a>

<a name="v12-eq-defect-commutator"></a>
**Equations (77)–(78).**

```math
\begin{aligned}
D_S^\top(M_S^\top\Omega M_S-\Omega)D_S
 &=h(F_S^\top\Omega+\Omega F_S),\\
 F_S^\top\Omega+\Omega F_S
 &=-h(I+hSJ)^{-1}(SJQ-QJS)(I+hJS)^{-1}.

\end{aligned}
```

Hence an off-locus gauge has a one-step symplectic defect of order
$`h^2`$, with leading coefficient controlled by $`SJQ-QJS`$.

*Proof.* Writing $`N_S=I+\tfrac h2F_S`$, commutativity of $`N_S`$ and $`D_S`$
gives $`M_SD_S=N_S`$. Therefore

```math
D_S^\top(M_S^\top\Omega M_S-\Omega)D_S
 =N_S^\top\Omega N_S-D_S^\top\Omega D_S,
```

which expands to [Equation (77)](#v12-eq-defect-cayley). Direct calculation gives

```math
F_S^\top\Omega+\Omega F_S
 =(I+hSJ)^{-1}Q-Q(I+hJS)^{-1}.
```

Multiplying the middle difference by its two denominators proves
[Equation (78)](#v12-eq-defect-commutator).

∎

**Theorem 6.4. Centralizer representation and generic dimension.**

<a name="v12-thm-generic-dimension"></a>

Let $`n=2m`$, $`C=JS`$, and define the $`\Omega`$-adjoint by
$`C^\dagger=\Omega^{-1}C^\top\Omega`$. Then

<a name="v12-eq-centralizer"></a>
**Equation (79).**

```math

 S\in\mathcal{L}
 \quad\Longleftrightarrow\quad
 [C,A]=0,\qquad C^\dagger=C,\qquad S=\Omega C.
```

If $`A`$ is nonderogatory, i.e. its minimal polynomial has degree $`2m`$,
then

<a name="v12-eq-generic-basis"></a>
**Equation (80).**

```math

 \mathcal{L}=\operatorname{span}\{\Omega,\Omega A^2,\ldots,\Omega A^{2m-2}\},
 \qquad \dim\mathcal{L}=m.
```

*Proof.* Substituting $`S=\Omega C`$ in [Equation (76)](#v12-eq-locus), and using
$`Q=\Omega A`$ and $`J\Omega=I`$, gives $`CA=AC`$. Moreover,
$`(\Omega C)^\top=-\Omega C`$ is equivalent to
$`C^\top\Omega=\Omega C`$, namely $`C^\dagger=C`$.

If $`A`$ is nonderogatory, its centralizer consists of
$`p(A)`$, $`\deg p<2m`$. Since $`A^\dagger=-A`$,

```math
p(A)^\dagger=p(-A).
```

The self-adjoint condition implies $`p(-A)=p(A)`$. The polynomial
$`p(-z)-p(z)`$ has degree below the minimal-polynomial degree and
annihilates $`A`$, hence it is zero. Thus $`p`$ is even. The $`m`$ even
powers in [Equation (80)](#v12-eq-generic-basis) are independent by the same
minimal-polynomial argument.

∎

**Corollary 6.5. Two dimensions and time reparametrization.**

<a name="v12-cor-two-dimensional"></a>

The form $`\Omega`$ always belongs to $`\mathcal{L}`$. If $`S=\sigma\Omega`$,
then the method is the midpoint rule with effective step

```math
h_{\mathrm{eff}}=\frac{h}{1+h\sigma}.
```

For $`n=2`$, every skew $`S`$ is a multiple of $`\Omega`$; consequently
the full gauge family is symplectic, although only $`S=0`$ is second order
for a fixed nondegenerate system.

*Proof.* For $`S=\sigma\Omega`$, $`JS=\sigma I`$, and the step equation becomes
$`(1+h\sigma)\Delta x=hA\bar x`$. In two dimensions
$`\Lambda^2V^*`$ is one-dimensional.

∎

**Example 6.6. Two uncoupled modes.**

<a name="v12-ex-two-modes"></a>

Let $`J=\operatorname{diag}(J_2,J_2)`$,
$`J_2=\left[\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right]`$, and
$`Q=\operatorname{diag}(\omega_1I_2,\omega_2I_2)`$, with
$`\omega_1\ne\omega_2>0`$. Then

```math
\mathcal{L}=\left\{\operatorname{diag}(\sigma_1J_2,\sigma_2J_2):
 \sigma_1,\sigma_2\in\mathbb{R}\right\}.
```

Every cross-mode skew block lies outside $`\mathcal{L}`$ unless it satisfies a
frequency-dependent degeneracy. This four-dimensional example is the
smallest setting in which exact energy and symplecticity are genuinely
separated by the gauge.

**Corollary 6.7. Semidiscrete wave equation.**

<a name="v12-cor-wave-locus"></a>

Let a wave semidiscretization with $`N`$ configuration degrees of freedom
have

```math
J=\begin{bmatrix}0&I_N\\-I_N&0\end{bmatrix},
 \qquad
 Q=\begin{bmatrix}K&0\\0&I_N\end{bmatrix},
 \qquad K=K^\top\succ0.
```

If $`K`$ has $`N`$ distinct eigenvalues, then the state dimension is
$`2N`$ and

```math
\dim\mathcal{L}=N,\qquad
 \mathcal{L}=\operatorname{span}\{\Omega,\Omega A^2,\ldots,\Omega A^{2N-2}\}.
```

If eigenvalues are repeated, the locus can be larger and its dimension is
not determined by $`N`$ alone.

*Proof.* The Hamiltonian matrix
$`A=\left[\begin{smallmatrix}0&I_N\\-K&0\end{smallmatrix}\right]`$
has eigenvalues $`\pm i\sqrt{\lambda_j(K)}`$. Under the simple-spectrum
hypothesis they are all distinct, so $`A`$ is nonderogatory.
[Theorem 6.4](#v12-thm-generic-dimension) applies with $`m=N`$. Repeated eigenvalues
enlarge the centralizer in general.

∎

<a name="v12-section-7"></a>

## 7. Nonlinear Hamiltonians and the universal segment rule

Let $`H\in C^1(V,\mathbb{R})`$. A two-point map
$`\bar\nabla H:V\times V\to V^*`$ is a discrete gradient if

<a name="v12-eq-dg"></a>
**Equation (81).**

```math

 \bar\nabla H(x_0,x_1)[x_1-x_0]=H(x_1)-H(x_0),
```

with the consistency condition
$`\bar\nabla H(x,x)=\,\mathrm{d} H(x)`$.

**Proposition 7.1. Pointwise nonlinear gauge.**

<a name="v12-prop-nonlinear-gauge"></a>

For a fixed pair $`x_0\ne x_1`$, if $`g_0`$ is one compatible covector,
then every compatible covector is

```math
g_0+r,\qquad r\in\operatorname{Ann}(\Delta x)
 =\{r\in V^*:r[\Delta x]=0\}.
```

The affine direction has dimension $`n-1`$. Moreover every such $`r`$
can be written $`r=\iota_{\Delta x}\mathsf S`$ for some
$`\mathsf S\in\Lambda^2V^*`$, although this representation is nonunique
for $`n>2`$.

*Proof.* The first assertion is immediate from subtracting two copies of
[Equation (81)](#v12-eq-dg). Choose $`\alpha\in V^*`$ with
$`\alpha[\Delta x]=1`$. Since $`r[\Delta x]=0`$,

```math
\iota_{\Delta x}(\alpha\wedge r)=r.
```

∎

The average vector field (AVF) discrete gradient is

<a name="v12-eq-avf"></a>
**Equation (82).**

```math

 \bar\nabla_{\operatorname{AVF}}H(x_0,x_1)
 =\int_0^1\,\mathrm{d} H(x_0+\tau\Delta x)\,\,\mathrm{d}\tau.
```

**Theorem 7.2. Uniqueness of the universal segment measure.**

<a name="v12-thm-avf-unique"></a>

Let $`\mu`$ be a finite signed Borel measure on $`[0,1]`$, independent of
$`H,x_0,x_1`$, and set

```math
\bar\nabla_\mu H(x_0,x_1)
 =\int_0^1\,\mathrm{d} H(x_0+\tau\Delta x)\,\,\mathrm{d}\mu(\tau).
```

This formula satisfies the discrete chain rule for every $`C^1`$
Hamiltonian and every segment if and only if $`\mu`$ is Lebesgue measure.
Thus AVF is the unique universal discrete gradient in this
segment-measure class.

*Proof.* Sufficiency is the fundamental theorem of calculus. For necessity,
restrict to a one-dimensional affine line and write
$`g(\tau)=H(x_0+\tau\Delta x)`$. Compatibility says

```math
\int_0^1 g'(\tau)\,\,\mathrm{d}\mu(\tau)=g(1)-g(0)
 =\int_0^1g'(\tau)\,\,\mathrm{d}\tau
```

for every $`g\in C^1[0,1]`$. Every continuous function is the derivative
of such a $`g`$; hence the two measures have equal integrals against all
continuous functions. The uniqueness part of the Riesz representation
theorem gives $`\mu=\,\mathrm{d}\tau`$.

∎

**Corollary 7.3.**

For quadratic $`H`$, AVF equals $`Q\bar x`$, the zero-gauge effort of
[Theorem 5.1](#v12-thm-two-point-gauge).

**Remark 7.4. What is and is not selected.**

Symmetry alone does not select AVF: the Gonzalez discrete gradient is also
symmetric [[13]](#v12-ref-gonzalez1996). AVF is affine-covariant because it
integrates a covector field along an affine segment; the standard Gonzalez
formula uses a chosen inner product through
$`\left\lVert \Delta x\right\rVert^{-2}`$. The uniqueness in
[Theorem 7.2](#v12-thm-avf-unique) is therefore precise but deliberately restricted to
universal segment-measure formulas. It is consistent with the established
order theory and structural characterization of AVF
[[26]](#v12-ref-quispel2008), [[10]](#v12-ref-eidnes2022).

**Remark 7.5. Symplecticity versus nonlinear energy preservation.**

No universal claim that the nonlinear symplectic gauge locus is empty is
made here. Classical no-go results require hypotheses concerning
symplecticity, exact energy preservation, consistency, and the absence of
additional first integrals. They imply, roughly, that a method satisfying
both properties for a generic nonintegrable Hamiltonian must coincide with
an exact-flow reparametrization [[11]](#v12-ref-ge1988), [[14]](#v12-ref-hairer2006). This is a
restriction on method classes, not a proof that every particular nonlinear
problem has an empty locus.

<a name="v12-section-8"></a>

## 8. Multiple stages: classification for invertible stage matrices

### 8.1. Bilinear chain-rule pairings

Let $`\mathcal A\in\mathbb{R}^{s\times s}`$ be a Runge–Kutta matrix,
let $`b\in\mathbb{R}^s`$, and let

```math
X=\boldsymbol{1}\otimes x_0+h(\mathcal A\otimes I)F,
 \qquad
 x_1=x_0+h(b^\top\otimes I)F.
```

For a quadratic Hamiltonian, let $`E=(I_s\otimes Q)X`$. A matrix
$`P\in\mathbb{R}^{s\times s}`$ defines the stage pairing

```math
\langle\!\langle E,F\rangle\!\rangle_P
 =\sum_{i,j=1}^sP_{ij}e_i^\top f_j.
```

**Proposition 8.1. Stage chain-rule conditions.**

<a name="v12-prop-p-conditions"></a>

The identity

```math
H(x_1)-H(x_0)
 =h\langle\!\langle E,F\rangle\!\rangle_P
```

holds for all $`x_0,F`$ and all symmetric $`Q`$ if and only if

<a name="v12-eq-p-conditions"></a>
**Equation (83).**

```math

 \boxed{\ \boldsymbol{1}^\top P=b^\top,\qquad
 \mathcal A^\top P+P^\top\mathcal A=bb^\top.\ }
```

*Proof.* Substitution gives

```math
h\langle\!\langle E,F\rangle\!\rangle_P
 =h\,x_0^\top Q(\boldsymbol{1}^\top P\otimes I)F
 +h^2F^\top(\mathcal A^\top P\otimes Q)F.
```

Only the symmetric part of $`\mathcal A^\top P`$ contributes. Comparing
this expression with

```math
\tfrac12(x_0+h(b^\top\otimes I)F)^\top Q
 (x_0+h(b^\top\otimes I)F)-\tfrac12x_0^\top Qx_0
```

for arbitrary data yields [Equation (83)](#v12-eq-p-conditions).

∎

**Theorem 8.2. All compatible stage pairings for invertible $`\mathcal A`$.**

<a name="v12-thm-p-family"></a>

Assume that $`\mathcal A`$ is invertible.  Let

```math
v=\mathcal A^{-1}\boldsymbol{1},\qquad \beta=b^\top v.
```

A compatible $`P`$ exists if and only if

<a name="v12-eq-beta"></a>
**Equation (84).**

```math

 \beta\in\{0,2\}.
```

When it exists, every solution is

<a name="v12-eq-p-general"></a>
**Equation (85).**

```math

 P=\mathcal A^{-\top}\left(\tfrac12bb^\top+K\right),
 \qquad K^\top=-K,
```

where

<a name="v12-eq-k-constraint"></a>
**Equation (86).**

```math

 v^\top K=(1-\tfrac12\beta)b^\top.
```

The solution set is affine of dimension

<a name="v12-eq-p-dimension"></a>
**Equation (87).**

```math

 \frac{(s-1)(s-2)}2.
```

Equivalently, if $`R(z)=1+z\,b^\top(I-z\mathcal A)^{-1}\boldsymbol{1}`$ is the
stability function, existence holds exactly when
$`|R(\infty)|=1`$.

*Proof.* Put $`X=\mathcal A^\top P`$. The second equation in
[Equation (83)](#v12-eq-p-conditions) says

```math
X=\tfrac12bb^\top+K,\qquad K^\top=-K,
```

which gives [Equation (85)](#v12-eq-p-general). The first equation becomes
[Equation (86)](#v12-eq-k-constraint). For every skew $`K`$,
$`v^\top Kv=0`$. Thus a necessary condition is

```math
(1-\tfrac12\beta)b^\top v
 =(1-\tfrac12\beta)\beta=0,
```

which is [Equation (84)](#v12-eq-beta). It is sufficient: every row $`w^\top`$ with
$`w^\top v=0`$ can be represented as $`v^\top K=w^\top`$, for example

```math
K=\frac{vw^\top-wv^\top}{v^\top v}.
```

The linear map $`K\mapsto v^\top K`$ from skew matrices has image
$`v^\perp`$, hence rank $`s-1`$. Its kernel, and therefore the affine
solution family, has dimension

```math
\frac{s(s-1)}2-(s-1)=\frac{(s-1)(s-2)}2.
```

Finally $`R(\infty)=1-b^\top\mathcal A^{-1}\boldsymbol{1}=1-\beta`$.

∎

**Remark 8.3. Singular stage matrices.**

<a name="v12-rem-singular-stage"></a>

[Proposition 8.1](#v12-prop-p-conditions) remains the complete solvability system when
$`\mathcal A`$ is singular, but [Theorem 8.2](#v12-thm-p-family) does not classify that
case.  It is nonempty: for the trapezoidal rule

```math
\mathcal A=
 \begin{bmatrix}0&0\\[1mm]\tfrac12&\tfrac12\end{bmatrix},
 \qquad b=\begin{bmatrix}\tfrac12\\[1mm]\tfrac12\end{bmatrix},
 \qquad
 P=\frac14\begin{bmatrix}1&1\\1&1\end{bmatrix}
```

satisfies [Equation (83)](#v12-eq-p-conditions).  A full singular classification is the
finite-dimensional range problem for

```math
\Phi(P)=\bigl(P^\top\boldsymbol{1},\,
 \mathcal A^\top P+P^\top\mathcal A\bigr);
```

it is deliberately not claimed here.

**Remark 8.4. Nondegeneracy and power skewness.**

Existence of $`P`$ does not imply $`\det P\ne0`$. Nor does the chain rule
alone make the pointwise stage operator $`I_s\otimes J`$ skew for the
$`P`$-pairing. For pointwise Hamiltonian stages $`f_i=Je_i`$, the
contraction vanishes automatically when $`P=P^\top`$. These are separate
requirements. The local Gauss pairing below has
$`P=\operatorname{diag}(b)\succ0`$, so both issues disappear.

### 8.2. Locality, the Cooper condition, and Gauss rigidity

Take distinct real collocation nodes $`c_1,\ldots,c_s`$, with Lagrange
polynomials $`\ell_j`$, and set

<a name="v12-eq-collocation"></a>
**Equation (88).**

```math

 a_{ij}=\int_0^{c_i}\ell_j(\tau)\,\,\mathrm{d}\tau,\qquad
 b_j=\int_0^1\ell_j(\tau)\,\,\mathrm{d}\tau,\qquad
 B=\operatorname{diag}(b).
```

Local diagonal stage pairing is $`P=B`$. The compatibility equations then
reduce to

<a name="v12-eq-cooper"></a>
**Equation (89).**

```math

 B\mathcal A+\mathcal A^\top B=bb^\top,
```

the Cooper condition.

**Theorem 8.5. Symplectic collocation is exactly Gauss collocation.**

<a name="v12-thm-gauss-rigidity"></a>

For a collocation method with $`s`$ distinct real nodes, the following
are equivalent:

1. the local diagonal pairing $`P=B`$ satisfies the chain rule;

1. the Cooper condition [Equation (89)](#v12-eq-cooper) holds;

1. the interpolatory quadrature is exact on $`\mathcal{P}_{2s-1}`$;

1. the nodal polynomial
     $`\pi(\tau)=\prod_{i=1}^s(\tau-c_i)`$ is orthogonal to
     $`\mathcal{P}_{s-1}`$ in $`L^2(0,1)`$;

1. the nodes and weights are the shifted Gauss–Legendre rule.

*Proof.* The first two statements are the specialization of
[Proposition 8.1](#v12-prop-p-conditions). Put

```math
L_j(t)=\int_0^t\ell_j(\tau)\,\,\mathrm{d}\tau.
```

The quadrature functional $`\mathcal Q(p)=\sum_kb_kp(c_k)`$ satisfies

```math
\mathcal Q\bigl((L_iL_j)'\bigr)
 =b_i a_{ij}+b_j a_{ji}.
```

Under the Cooper condition this equals
$`b_ib_j=L_i(1)L_j(1)=\int_0^1(L_iL_j)'`$.
The $`L_i`$ span

```math
\{p\in\mathcal{P}_s:p(0)=0\}=t\mathcal{P}_{s-1}.
```

Their pairwise products span $`t^2\mathcal{P}_{2s-2}`$, and their derivatives
span $`t\mathcal{P}_{2s-2}`$, i.e. all monomials
$`t,t^2,\ldots,t^{2s-1}`$. Constants are already integrated exactly.
Thus Cooper implies exactness on $`\mathcal{P}_{2s-1}`$. The converse follows by
applying exactness to each $`(L_iL_j)'`$.

An $`s`$-node interpolatory rule is exact on $`\mathcal{P}_{2s-1}`$ exactly when,
for every $`q\in\mathcal{P}_{s-1}`$,

```math
\int_0^1\pi(\tau)q(\tau)\,\,\mathrm{d}\tau=0.
```

The unique monic degree-$`s`$ polynomial with this orthogonality property is
the shifted Legendre polynomial. Its roots are automatically distinct and
interior to $`(0,1)`$; they are the Gauss nodes, and the
interpolatory weights are then the Gauss weights.

∎

**Corollary 8.6. Polynomial energy proof.**

<a name="v12-cor-polynomial-energy"></a>

Let $`u\in\mathcal{P}_s([0,1];V)`$ be a collocation polynomial for
$`u'=hJQu`$. Then Gauss nodes imply exact conservation of
$`\tfrac12u^\top Qu`$.

*Proof.* The polynomial $`p(\tau)=(Qu(\tau))^\top u'(\tau)`$ has degree at most
$`2s-1`$ and vanishes at every collocation node. Hence
$`p=\pi q`$, with $`q\in\mathcal{P}_{s-1}`$. Orthogonality in
[Theorem 8.5](#v12-thm-gauss-rigidity) gives

```math
H(u(1))-H(u(0))=\int_0^1p(\tau)\,\,\mathrm{d}\tau=0.
```

∎

**Proposition 8.7. From the local Stokes pairing to active stages.**

<a name="v12-prop-local-stokes-to-stages"></a>

Let the nodes $`c_i`$ be distinct and real, let
$`Q=Q^\top\ne0`$, and let
$`x_h\in\mathcal P_s([0,1];X)`$.  Write its scaled flow as

<a name="v12-eq-polynomial-flow-interpolation"></a>
**Equation (90).**

```math

 h^{-1}x_h'(\tau)=k_h(\tau)
 =\sum_{j=1}^s\ell_j(\tau)k_j,
 \qquad X_i=x_h(c_i).
```

Then polynomial integration gives, without invoking a Runge–Kutta
ansatz,

<a name="v12-eq-stokes-generated-stages"></a>
**Equation (91).**

```math

 X_i=x_0+h\sum_{j=1}^sa_{ij}k_j,\qquad
 x_1=x_0+h\sum_{j=1}^sb_jk_j.
```

For the quadratic effort $`e_h=Qx_h`$, the local diagonal pairing

<a name="v12-eq-local-gauss-stokes-pairing"></a>
**Equation (92).**

```math

 h\sum_{i=1}^sb_i(QX_i)^\top k_i
```

is the exact cap difference
$`H(x_1)-H(x_0)`$ for every $`x_h\in\mathcal P_s`$ if and only if the
Cooper condition holds, equivalently if and only if the nodes and weights
are the shifted Gauss rule.  If additionally $`J^\top=-J`$, closing the
active variables by $`k_i=JQX_i`$ annihilates
[Equation (92)](#v12-eq-local-gauss-stokes-pairing) and turns
[Equation (91)](#v12-eq-stokes-generated-stages) into the Gauss stage equations.

*Proof.* Integrating [Equation (90)](#v12-eq-polynomial-flow-interpolation) from $`0`$ to
$`c_i`$ and to $`1`$ gives [Equation (91)](#v12-eq-stokes-generated-stages).  Moreover,

```math
H(x_1)-H(x_0)
 =h\int_0^1(Qx_h)^\top k_h\,\,\mathrm{d}\tau .
```

The integrand has degree at most $`2s-1`$, so Gauss exactness proves
sufficiency.  Conversely, suppose the pairing is exact for every
$`x_h`$.  Since $`Q`$ is real symmetric and nonzero, some $`v\in X`$
satisfies $`v^\top Qv\ne0`$.  Taking $`x_h(\tau)=p(\tau)v`$ proves
quadrature exactness for $`pp'=\tfrac12(p^2)'`$, for every
$`p\in\mathcal P_s`$.  Polarization with $`p+q`$ then proves exactness
for

```math
(pq)'=p'q+pq'\qquad(p,q\in\mathcal P_s).
```

The products $`pq`$ span $`\mathcal P_{2s}`$, hence their derivatives
span $`\mathcal P_{2s-1}`$.  Thus the quadrature is exact on the whole
space $`\mathcal P_{2s-1}`$, and [Theorem 8.5](#v12-thm-gauss-rigidity) gives the
Cooper condition and the shifted Gauss rule.  Under the constitutive law
each summand is
$`(QX_i)^\top J(QX_i)=0`$.  Thus the active stage equations follow from
the local polynomial Stokes pairing plus constitutive closure, and their
later Schur reduction is not an arbitrary rewriting of a preselected
Runge–Kutta method.

∎

**Theorem 8.8. Active Gauss cell and induced cap gauge.**

<a name="v12-thm-active-gauss-cell"></a>

Let $`X=\mathbb{R}^n`$, let $`A=JQ`$, and take the shifted $`s`$-point
Gauss–Legendre nodes, Lagrange polynomials, matrix $`\mathcal A_s`$, and
weights $`b`$ from [Equation (88)](#v12-eq-collocation).  The polynomial temporal cell

```math
\mathcal P_s([0,1];X)
 \xrightarrow{\ \partial_\tau\ }
 \mathcal P_{s-1}([0,1];X)
```

is represented by active stage states $`X_i`$ and stage flows $`k_i`$:

<a name="v12-eq-active-gauss-stages"></a>

<a name="v12-eq-active-gauss-output"></a>
**Equations (93)–(94).**

```math
\begin{aligned}
X_i&=x_0+h\sum_{j=1}^s a_{ij}k_j,
 &k_i&=AX_i,\\
 x_1&=x_0+h\sum_{j=1}^s b_jk_j.

\end{aligned}
```

Put

<a name="v12-eq-active-stage-block"></a>
**Equation (95).**

```math

 \mathcal K_s(h,A)=I_s\otimes I-h\mathcal A_s\otimes A.
```

Whenever this block is invertible, elimination of the $`sn`$ stage
variables gives

<a name="v12-eq-active-stage-schur-map"></a>
**Equation (96).**

```math

 M_s^\mathrm{cell}(hA)
 =I+h(b^\top\otimes A)\mathcal K_s(h,A)^{-1}(\boldsymbol{1}\otimes I)
 =R_s(hA),
```

where

<a name="v12-eq-pade-polynomial"></a>

<a name="v12-eq-diagonal-pade-map"></a>
**Equations (97)–(98).**

```math
\begin{aligned}
P_s(z)&=\sum_{k=0}^s
 \frac{(2s-k)!\,s!}{(2s)!\,k!\,(s-k)!}\,z^k,
 \\
 R_s(z)&=P_s(-z)^{-1}P_s(z).

\end{aligned}
```

Assume in addition that $`J,Q`$, $`R_s(hA)-I`$, and the displayed
matrix functions below are invertible.  Write

```math
E_s(z)=\tfrac12\bigl(P_s(z)+P_s(-z)\bigr),\qquad
 O_s(z)=\tfrac12\bigl(P_s(z)-P_s(-z)\bigr)
       =z\widehat O_s(z).
```

The cap gauge induced by this interior Schur reduction is

<a name="v12-eq-stage-induced-gauss-gauge"></a>
**Equation (99).**

```math

 \boxed{\quad
 S_{s,h}^\mathrm{cell}
 =\frac1h\Omega
 \left[
  \frac12E_s(hA)\widehat O_s(hA)^{-1}-I
 \right],
 \qquad \Omega=J^{-1}.
 \quad}
```

It is alternating, and its two-point map is exactly the $`s`$-stage
Gauss method.

*Proof.* The polynomial reconstruction and the two stage relations follow from
[Proposition 8.7](#v12-prop-local-stokes-to-stages); conversely they reconstruct a unique
$`x_h\in\mathcal P_s`$.  Stacking the
stage equations gives

```math
\mathcal K_s(h,A)X=(\boldsymbol{1}\otimes I)x_0,
```

and substitution into the endpoint equation proves the first equality in
[Equation (96)](#v12-eq-active-stage-schur-map).

The scalar stability function of Gauss collocation has numerator and
denominator of degree at most $`s`$ and agrees with $`e^z`$ through degree
$`2s`$; it is therefore the diagonal Padé approximant
$`P_s(-z)^{-1}P_s(z)`$.  The same rational identity holds for matrices on
the common resolvent set, proving the second equality.

It remains to identify the cap two-form generated by the stages.  Since
all polynomials in $`A`$ commute,

```math
(R_s+I)(R_s-I)^{-1}
 =E_s(hA)O_s(hA)^{-1}
 =(hA)^{-1}E_s(hA)\widehat O_s(hA)^{-1}.
```

Substitution into the inverse Cayley formula
[Equation (71)](#v12-eq-inverse-gauge-cayley) gives
[Equation (99)](#v12-eq-stage-induced-gauss-gauge).  Alternation follows either from
[Theorem 5.3](#v12-thm-gauge-cayley-chart) and $`R_s^\top QR_s=Q`$, or directly from
$`A^\top\Omega+\Omega A=0`$.

∎

**Remark 8.9. Representation versus interior generation.**

<a name="v12-rem-representation-generation"></a>

The trace-generated theorem [Theorem 4.14](#v12-thm-trace-generated-condensation)
starts from $`S_h`$ and realizes it in a redundant raw cell.  By contrast,
[Theorem 8.8](#v12-thm-active-gauss-cell) starts from the standard local Gauss pairing
on $`\mathcal P_s\to\mathcal P_{s-1}`$ and the constitutive stage law
$`k_i=AX_i`$; the cap gauge is computed only after the active stage block
has been eliminated.  This is a genuine interior generation result.

**Remark 8.10.**

The same monic-polynomial orthogonality is the scalar core of the spatial
$`\mathcal{P}_N/\mathcal{P}_{N-1}`$ Stokes-pairing construction
[[20]](#v12-ref-moulla2012). The common mechanism is not merely that both methods
use Gauss nodes: in both variables, scalar locality turns compatibility
into maximal-degree quadrature exactness, and maximal-degree exactness is
equivalent to orthogonality of the nodal polynomial.

<a name="v12-section-9"></a>

## 9. Space–time Stokes–Dirac closure and local rigidity

### 9.1. Maximality, constitutive closure, and glued faces

Let active state flows and efforts lie in dual spaces of dimension $`d`$,
and let a boundary port have dimension $`r`$.  Coordinates may be chosen
after restricting to the active quotient in
[Equation (33)](#v12-eq-intrinsic-space-time-gauge).  Consider

<a name="v12-eq-open-dirac"></a>
**Equation (100).**

```math

 \mathcal{D}_h=\left\{(f,u;e,y):
 f=-\mathcal J_he-Bu,\quad y=B^\top e+Nu\right\},
 \qquad \mathcal J_h^\top=-\mathcal J_h,\quad N^\top=-N.
```

The symmetric split-power pairing is

<a name="v12-eq-plus-pairing"></a>
**Equation (101).**

```math

 \left\langle\!\left\langle
 (f,u;e,y),(f',u';e',y')
 \right\rangle\!\right\rangle_+
 =e^\top f'+e'^\top f+y^\top u'+y'^\top u.
```

**Proposition 9.1. Finite-dimensional open Dirac structure with feedthrough.**

<a name="v12-prop-dirac-maximal"></a>

The subspace $`\mathcal{D}_h`$ in [Equation (100)](#v12-eq-open-dirac) is maximally isotropic for
[Equation (101)](#v12-eq-plus-pairing).

*Proof.* For any element of $`\mathcal{D}_h`$,

```math
e^\top f+y^\top u
 =-e^\top\mathcal J_he-e^\top Bu+e^\top Bu+u^\top Nu=0.
```

Polarization proves isotropy. The free variables $`(e,u)`$ give
$`\dim\mathcal{D}_h=d+r`$, half the dimension of the ambient split space.
Therefore the isotropic subspace is maximal.

∎

**Theorem 9.2. Rectangular Stokes-to-Dirac reduction.**

<a name="v12-thm-rectangular-stokes-dirac"></a>

Let $`E`$ be a raw effort space, $`X`$ an active state space, and $`U`$
a port-input space.  Consider linear maps

```math
P:E\to X^*,\qquad K:E\to X,\qquad
 \Gamma:E\to U,\qquad Y:E\to U^*.
```

Assume the rectangular Stokes identity

<a name="v12-eq-abstract-rectangular-stokes"></a>
**Equation (102).**

```math

 \left\langle Pe,Ke\right\rangle=\left\langle Ye,\Gamma e\right\rangle
 \qquad(e\in E)
```

and suppose that

<a name="v12-eq-effort-input-reconstruction"></a>
**Equation (103).**

```math

 \mathcal R=\begin{bmatrix}P\\\Gamma\end{bmatrix}:
 E\longrightarrow X^*\oplus U
```

is an isomorphism.  Define the square blocks by

<a name="v12-eq-rectangular-square-blocks"></a>
**Equation (104).**

```math

 \begin{bmatrix}\mathcal J&B\end{bmatrix}=K\mathcal R^{-1},
 \qquad
 \begin{bmatrix}C&N\end{bmatrix}=Y\mathcal R^{-1}.
```

Then

<a name="v12-eq-rectangular-dirac-identities"></a>
**Equation (105).**

```math

 \mathcal J^*=-\mathcal J,\qquad C=B^*,\qquad N^*=-N.
```

Consequently the raw equations reduce exactly to

<a name="v12-eq-rectangular-open-reduction"></a>
**Equation (106).**

```math

 f=\mathcal J\epsilon+Bu,\qquad
 y=B^*\epsilon+Nu,\qquad
 \epsilon=Pe,\quad u=\Gamma e,\quad f=Ke,\quad y=Ye,
```

and the oriented tuple $`(-f,u;\epsilon,y)`$ belongs to the maximally
isotropic relation [Equation (100)](#v12-eq-open-dirac) with
$`\mathcal J_h=\mathcal J`$.

*Proof.* For arbitrary $`(\epsilon,u)`$, write
$`e=\mathcal R^{-1}(\epsilon,u)`$.  Then
[Equation (102)](#v12-eq-abstract-rectangular-stokes) becomes

```math
\left\langle \epsilon,\mathcal J\epsilon+Bu\right\rangle
 =\left\langle C\epsilon+Nu,u\right\rangle.
```

The pure quadratic terms imply $`\mathcal J^*=-\mathcal J`$ and
$`N^*=-N`$; polarization of the mixed term gives $`C=B^*`$.
[Equation (106)](#v12-eq-rectangular-open-reduction) follows from the definitions,
and [Proposition 9.1](#v12-prop-dirac-maximal) gives maximal isotropy.

∎

**Theorem 9.3. Cap-to-mantle power balance.**

<a name="v12-thm-space-time-step"></a>

Let $`H_h`$ be a discrete Hamiltonian on a temporal cap and let
$`\bar e(x_0,x_1)`$ be any compatible temporal boundary reduction of a
space–time pairing, so that

```math
\bar e^\top\Delta x=H_h(x_1)-H_h(x_0).
```

Impose

<a name="v12-eq-dirac-step"></a>
**Equation (107).**

```math

 \left(-\frac{\Delta x}{h},\bar u;\bar e,\bar y\right)\in\mathcal{D}_h.
```

Then

<a name="v12-eq-dirac-power"></a>
**Equation (108).**

```math

 H_h(x_1)-H_h(x_0)=h\,\bar y^\top\bar u.
```

The result holds for every active space–time gauge whose pullback gives
the chosen temporal effort.  The relation $`\mathcal{D}_h`$ supplies maximal
isotropy, whereas the compatible pairing supplies constitutive closure;
neither statement implies the other.

*Proof.* Self-isotropy of [Equation (107)](#v12-eq-dirac-step) gives
$`-\bar e^\top\Delta x/h+\bar y^\top\bar u=0`$.
The chain rule converts the cap term into the Hamiltonian increment.
Gauge independence follows because an alternating two-form vanishes when
evaluated twice on the same active increment.

∎

**Corollary 9.4. Dissipation and gluing.**

<a name="v12-cor-composition"></a>

If a resistive relation adds $`R_h=R_h^\top\succeq0`$, then

<a name="v12-eq-dissipative-power"></a>
**Equation (109).**

```math

 H_h(x_1)-H_h(x_0)
 =h\,\bar y^\top\bar u-h\,\bar e^\top R_h\bar e.
```

If two space–time cells share a face and are interconnected with opposite
orientations,
$`\bar u_1=-\bar y_2`$ and $`\bar u_2=\bar y_1`$, their internal supplies
cancel.  Summing [Equation (109)](#v12-eq-dissipative-power) over a mesh leaves only the
external mantle and the resistive volume terms, exactly as internal faces
cancel in a cellular boundary.

### 9.2. Local diagonal rigidity derived from the total equations

Let $`c_i^t`$, $`1\le i\le N_t`$, and
$`c_j^x`$, $`1\le j\le N_x`$, be distinct nodes with Lagrange bases
$`\ell_i^t`$ and $`\ell_j^x`$.  A pointwise diagonal pairing on the tensor
grid is determined by a weight array $`W=(W_{ji})`$ through

<a name="v12-eq-local-tensor-pairing"></a>
**Equation (110).**

```math

 \left\langle p,q\right\rangle_W
 =\sum_{i=1}^{N_t}\sum_{j=1}^{N_x}
 W_{ji}\,p(c_i^t,c_j^x)q(c_i^t,c_j^x).
```

Suppose the one-dimensional local compatible weight spaces are the lines
$`\operatorname{span}\{w^t\}`$ and $`\operatorname{span}\{w^x\}`$.  This conclusion follows from
[Theorem 8.5](#v12-thm-gauss-rigidity) for maximal polynomial collocation and from its
spatial Stokes analogue.

**Theorem 9.5. Tensor rigidity from global Stokes compatibility.**

<a name="v12-thm-tensor-rigidity"></a>

Assume that the full tensor trial space contains all products
$`\ell_i^t p_x`$ and $`p_t\ell_j^x`$, and that the diagonal pairing
[Equation (110)](#v12-eq-local-tensor-pairing) satisfies the two components of the global
Stokes identity associated with $`D_\mathrm{st}`$.  Then

<a name="v12-eq-tensor-weight-scaled"></a>
**Equation (111).**

```math

 W=\alpha\,w^x(w^t)^\top
```

for a scalar $`\alpha`$.  If the one-dimensional weights and $`W`$ are
normalized, $`\alpha=1`$.  If the factor rules are polynomial rules of
maximal degree, the unique normalized local pairing is

<a name="v12-eq-tensor-weight"></a>
**Equation (112).**

```math

 \boxed{\quad W=w^x(w^t)^\top
 =\mathrm{GL}_x\otimes\mathrm{GL}_t.\quad}
```

*Proof.* Fix a temporal node $`c_i^t`$ and test the spatial component of the global
identity with functions $`\ell_i^t p_x`$ and $`\ell_i^t q_x`$.  Since
$`\ell_i^t(c_k^t)=\delta_{ik}`$, every temporal slice except $`i`$
vanishes in [Equation (110)](#v12-eq-local-tensor-pairing).  The column $`W_{:i}`$ must
therefore be a one-dimensional spatial compatible weight vector.  Hence

```math
W_{:i}=\lambda_iw^x,
 \qquad W=w^x\lambda^\top.
```

Now fix a spatial node and use $`p_t\ell_j^x`$, $`q_t\ell_j^x`$ in the
temporal component.  Every row of $`W`$ belongs to
$`\operatorname{span}\{(w^t)^\top\}`$.  Choosing $`j`$ with $`w_j^x\ne0`$ gives
$`\lambda=\alpha w^t`$, proving
[Equation (111)](#v12-eq-tensor-weight-scaled).  Normalization fixes $`\alpha=1`$, and
maximal polynomial exactness selects the Gauss factors by
the nodal-polynomial argument in [Theorem 8.5](#v12-thm-gauss-rigidity), applied on
each axis.

∎

**Remark 9.6. What locality removes.**

The conclusion is not a uniqueness statement in the full affine family
[Equation (32)](#v12-eq-space-time-affine-family).  The local diagonal ansatz removes the
nonlocal sectors in [Equation (37)](#v12-eq-wedge-sector-split); the global Stokes
equations then force separability.  A nonzero mixed gauge is therefore a
precise algebraic obstruction to representing the pairing by pointwise
product weights, not a violation of Stokes compatibility.

### 9.3. A wave system with boundary ports

Consider the first-order wave system on a spatial nodal grid.  Let
$`W_x=W_x^\top\succ0`$ and $`D_x`$ satisfy the summation-by-parts identity

<a name="v12-eq-sbp-wave"></a>
**Equation (113).**

```math

 W_xD_x+D_x^\top W_x=B_x,
 \qquad B_x=\operatorname{diag}(-1,0,\ldots,0,1).
```

For $`x=(q,p)`$, set

<a name="v12-eq-wave-qj"></a>
**Equation (114).**

```math

 Q_h=\operatorname{diag}(W_x,W_x),\qquad
 \mathcal J_h=
 \begin{bmatrix}
 0&D_xW_x^{-1}\\ D_xW_x^{-1}&0
 \end{bmatrix},
 \qquad
 \mathcal R_h=\rho Q_h^{-1},\quad\rho\ge0.
```

Here $`\mathcal J_h`$ is the unclosed Stokes operator, not a skew matrix:
its symmetric defect is exactly the boundary term encoded by
[Equation (113)](#v12-eq-sbp-wave).  Introducing the two endpoint ports converts this
defect into an open Dirac relation of the form [Equation (100)](#v12-eq-open-dirac).
For any skew gauge $`S`$, define

<a name="v12-eq-wave-gauge-effort"></a>
**Equation (115).**

```math

 g=Q_h\bar x-S\Delta x,
 \qquad v=Q_h^{-1}g=(v_q,v_p),
```

and close the cell by

<a name="v12-eq-wave-step"></a>
**Equation (116).**

```math

 \frac{\Delta x}{h}=(\mathcal J_h-\mathcal R_h)g.
```

**Proposition 9.7. Gauge-independent wave power.**

<a name="v12-prop-wave-power"></a>

Every well-defined solution of [Equation (116)](#v12-eq-wave-step) satisfies

<a name="v12-eq-wave-power"></a>
**Equation (117).**

```math
\begin{aligned}
H_h(x_1)-H_h(x_0)
 &=h\bigl(v_{q,N}v_{p,N}-v_{q,0}v_{p,0}\bigr)\\
 &\quad-h\rho
 \bigl(v_q^\top W_xv_q+v_p^\top W_xv_p\bigr),

\end{aligned}
```

where $`H_h(x)=\tfrac12x^\top Q_hx`$.  The first line is the right-face
power minus the left-face power.  On two adjacent cells, continuity of
the two boundary efforts makes the shared-face contributions cancel with
opposite orientations.

*Proof.* Skewness of $`S`$ gives

```math
H_h(x_1)-H_h(x_0)=g^\top\Delta x.
```

Moreover, [Equation (113)](#v12-eq-sbp-wave) yields

```math
g^\top\mathcal J_hg
 =v_q^\top W_xD_xv_p+v_p^\top W_xD_xv_q
 =v_q^\top B_xv_p,
```

which is the oriented boundary term in [Equation (117)](#v12-eq-wave-power); and
$`g^\top\mathcal R_hg =\rho(v_q^\top W_xv_q+v_p^\top W_xv_p)`$.
Substitute [Equation (116)](#v12-eq-wave-step).  Shared faces occur twice with opposite
entries of $`B_x`$, proving cancellation.

∎

### 9.4. A direct rectangular transmission-line cylinder

We now apply the cylinder gauge to the polynomial reduction of the ideal
transmission line in [[20]](#v12-ref-moulla2012), rather than first selecting a
spatial state-space model and then applying an unrelated time integrator.
Let

```math
E_N=\mathcal{P}_N([0,\ell]),\qquad F_N=\mathcal{P}_{N-1}([0,\ell]),
 \qquad D_N=\partial_z:E_N\longrightarrow F_N,
```

and let the full-rank rectangular pairing $`M_N:E_N\times F_N\to\mathbb{R}`$
satisfy

<a name="v12-eq-rectangular-stokes-line"></a>
**Equation (118).**

```math

 M_ND_N+D_N^\top M_N^\top
 =r_\ell r_\ell^\top-r_0r_0^\top,
```

where $`r_0,r_\ell\in E_N^*`$ are endpoint evaluation.  This is the
coordinate form of the exact polynomial Stokes identity used in the 2012
construction.  No square same-grid differentiation matrix is assumed.
The temporal factor is the minimal exact complex
$`\mathcal{P}_1([0,h])\xrightarrow{\partial_t}\mathcal{P}_0([0,h])`$: its traces are the
two caps and its single flow coordinate is the cap difference divided by
$`h`$.  Consequently the following equations are the simultaneous
closure of one tensor cylinder cell, even though they may subsequently be
eliminated to a cap map.

On a temporal slab, take cap states $`x_a=(q_a,p_a)\in F_N\oplus F_N`$,
$`a\in\{0,1\}`$, and spatial effort polynomials
$`e_q,e_p\in E_N`$.  Put

<a name="v12-eq-active-line-efforts"></a>
**Equation (119).**

```math

 \epsilon=
 \begin{bmatrix}M_N^\top e_q\\M_N^\top e_p\end{bmatrix},
 \qquad
 H_N(x)=\tfrac12x^\top Q_Nx,
 \qquad Q_N=Q_N^\top\succ0.
```

For a skew cap gauge $`S_h^\top=-S_h`$, impose simultaneously

<a name="v12-eq-rectangular-line-flow"></a>

<a name="v12-eq-rectangular-line-closure"></a>
**Equations (120)–(121).**

```math
\begin{aligned}
\frac{q_1-q_0}{h}&=-D_Ne_p,
 &\frac{p_1-p_0}{h}&=-D_Ne_q,
 \\
 \epsilon&=Q_N\frac{x_0+x_1}{2}-S_h(x_1-x_0).

\end{aligned}
```

Choose the boundary causality

<a name="v12-eq-line-port-causality"></a>
**Equation (122).**

```math

 u=\begin{bmatrix}e_q(0)\\e_p(\ell)\end{bmatrix},
 \qquad
 y=\begin{bmatrix}e_p(0)\\-e_q(\ell)\end{bmatrix}.
```

We denote these oriented components by
$`u=(u_L,u_R)^\top`$ and $`y=(y_L,y_R)^\top`$, so that the left supply is
$`y_Lu_L`$ and the right supply is $`y_Ru_R`$.  This convention is used
unchanged when neighboring cells are interconnected.
The homogeneous conditions $`u=0`$ are the complementary endpoint
conditions used for the closed ideal line in [[20]](#v12-ref-moulla2012).

**Proposition 9.8. Direct cap–mantle balance for the rectangular line.**

<a name="v12-prop-rectangular-line-power"></a>

Every well-defined solution of
[Equation (120)](#v12-eq-rectangular-line-flow), [Equation (121)](#v12-eq-rectangular-line-closure) satisfies

<a name="v12-eq-rectangular-line-power"></a>
**Equation (123).**

```math

 H_N(x_1)-H_N(x_0)=h\,y^\top u.
```

If two such cells are glued by
$`u_{1,R}=y_{2,L}`$ and $`u_{2,L}=-y_{1,R}`$, their interface
supplies cancel.  Adding a symmetric positive semidefinite active
resistance $`R_N`$ to the flow relation changes the right-hand side to

```math
h\,y^\top u-h\,\epsilon^\top R_N\epsilon.
```

*Proof.* Skewness of $`S_h`$ and [Equation (121)](#v12-eq-rectangular-line-closure) give

```math
\begin{aligned}
H_N(x_1)-H_N(x_0)
 &=\epsilon_q^\top(q_1-q_0)+\epsilon_p^\top(p_1-p_0)\\
 &=-h\bigl(e_q^\top M_ND_Ne_p
            +e_p^\top M_ND_Ne_q\bigr)\\
 &=-h\,e_q^\top
   \bigl(M_ND_N+D_N^\top M_N^\top\bigr)e_p.
\end{aligned}
```

Using [Equation (118)](#v12-eq-rectangular-stokes-line), the last expression is

```math
h\bigl(e_q(0)e_p(0)-e_q(\ell)e_p(\ell)\bigr)=h\,y^\top u.
```

The interconnection contribution is
$`y_{1,R}y_{2,L}-y_{2,L}y_{1,R}=0`$.  The resistive statement follows
by contracting the additional flow $`-R_N\epsilon`$ with $`\epsilon`$.

∎

The cap gauge in [Equation (121)](#v12-eq-rectangular-line-closure) can now be placed
inside the exact visibility sequence of
[Theorem 4.10](#v12-thm-volume-cap-reduction).  Write
$`\mathcal D_x=I_2\otimes D_N`$ and define the cap trace

<a name="v12-eq-line-cap-trace"></a>
**Equation (124).**

```math

 \mathcal T_N(z_0,z_1)=\mathcal D_x(z_1-z_0),
 \qquad
 \mathcal T_N:(E_N\oplus E_N)^2\to F_N\oplus F_N.
```

For an alternating form
$`\Theta\in\Lambda^2(F_N\oplus F_N)^*`$, define

<a name="v12-eq-mixed-cylinder-lift"></a>
**Equation (125).**

```math

 \omega_{tx}^{\Theta}
 =\mathcal T_N^*\Theta,\qquad
 \omega_{tx}^{\Theta}(z,w)
 =\Theta(\mathcal T_Nz,\mathcal T_Nw).
```

**Proposition 9.9. Surjective mixed realization of all cap gauges.**

<a name="v12-prop-mixed-cylinder-lift"></a>

The form $`\omega_{tx}^{\Theta}`$ is well defined on the active
space–time quotient.  It vanishes on time-constant fields and on fields
whose spatial part lies in $`\ker\mathcal D_x`$; hence it is supported on
the mixed sector $`Q_{tx}`$ of [Equation (36)](#v12-eq-quotient-sector-split).
Equivalently, after passage to the active quotient, there is a unique onto
map $`\widetilde{\mathcal T}_N`$ such that

<a name="v12-eq-line-trace-mixed-factorization"></a>
**Equation (126).**

```math

 \boxed{\quad
 \mathcal T_N=\widetilde{\mathcal T}_N\pi_{tx}.
 \quad}
```

Moreover, $`\mathcal T_N`$ is onto.  If
$`\sigma_x:F_N\oplus F_N\to E_N\oplus E_N`$ is any right inverse of
$`\mathcal D_x`$, then

<a name="v12-eq-line-geometric-lifting"></a>
**Equation (127).**

```math

 \mathcal L_N\xi=(0,\sigma_x\xi)
```

is an explicit geometric right inverse with
$`\mathcal T_N\mathcal L_N=I`$.  More generally, every regular interior
reconstruction $`L_h`$ consistent with the cap trace satisfies
$`\mathcal T_NL_h=I`$, and

<a name="v12-eq-mixed-cap-reduction"></a>
**Equation (128).**

```math

 \boxed{\quad
 \mathcal R_{L_h}(\omega_{tx}^{\Theta})
 =L_h^*\mathcal T_N^*\Theta=\Theta.
 \quad}
```

Consequently every alternating cap gauge on the active line state has an
explicit purely mixed cylinder realization.

*Proof.* Changing either cap potential by an element of
$`\ker\mathcal D_x`$ leaves [Equation (125)](#v12-eq-mixed-cylinder-lift) unchanged.
If $`z_1=z_0`$, the temporal increment vanishes; if both spatial
increments lie in the spatial kernel, their images under
$`\mathcal D_x`$ vanish.  Alternation follows from
$`\Theta^\top=-\Theta`$.  These are precisely the quotient and support
properties defining the $`Q_{tx}`$ sector.  The universal property of the
quotient gives [Equation (126)](#v12-eq-line-trace-mixed-factorization).  Since
$`D_N:\mathcal P_N\to\mathcal P_{N-1}`$ is onto,
$`\mathcal D_x`$ and $`\mathcal T_N`$ are onto and a right inverse
$`\sigma_x`$ exists.  Finally,

```math
L_h^*\mathcal T_N^*\Theta
 =(\mathcal T_NL_h)^*\Theta=\Theta.
```

This is the injective-reconstruction case of
[Theorem 4.10](#v12-thm-volume-cap-reduction).

∎

### 9.5. Interior Schur generation versus encoded realization

The pullback in [Proposition 9.9](#v12-prop-mixed-cylinder-lift) places a prescribed
two-form in the cylinder.  A different question is whether eliminating
genuine auxiliary variables can produce the same cap form.  The following
finite-dimensional result gives a sharp algebraic answer.

**Theorem 9.10. Minimal auxiliary Schur generation.**

<a name="v12-thm-minimal-schur-gauge-generation"></a>

Let $`V`$ be an $`n`$-dimensional cap space and let
$`\Theta_0,\Theta_1\in\Lambda^2V^*`$.  Put
$`\Delta=\Theta_1-\Theta_0`$ and
$`\operatorname{rank}\Delta=2q`$.  Then there exist an interior space
$`W`$ of dimension $`2q`$, a nondegenerate alternating matrix
$`K:W\to W^*`$, and a coupling $`B:W\to V^*`$ such that

<a name="v12-eq-skew-interior-generator"></a>
**Equation (129).**

```math

 \mathbb G=
 \begin{bmatrix}
  \Theta_0&B\\-B^\top&K
 \end{bmatrix},
 \qquad \mathbb G^\top=-\mathbb G,
```

and elimination of the interior variable from

<a name="v12-eq-skew-interior-elimination"></a>
**Equation (130).**

```math

 p=\Theta_0x+B\eta,
 \qquad 0=-B^\top x+K\eta
```

gives exactly

<a name="v12-eq-schur-generated-cap-gauge"></a>
**Equation (131).**

```math

 p=\Theta_1x.
```

Moreover, $`2q`$ is the smallest possible interior dimension among all
representations of the form [Equation (129)](#v12-eq-skew-interior-generator) with
invertible $`K`$.

*Proof.* The real congruence normal form for an alternating matrix gives a full
column-rank matrix $`U\in\mathbb R^{n\times2q}`$ such that

```math
\Delta=U\mathbb J_qU^\top,
 \qquad
 \mathbb J_q=\operatorname{diag}\!\left(
 \begin{bmatrix}0&1\\-1&0\end{bmatrix},\ldots,
 \begin{bmatrix}0&1\\-1&0\end{bmatrix}\right).
```

Take $`B=U`$ and $`K=-\mathbb J_q`$, so that
$`K^{-1}=\mathbb J_q`$.  The second equation in
[Equation (130)](#v12-eq-skew-interior-elimination) gives
$`\eta=K^{-1}B^\top x`$, and substitution in the first gives

```math
p=(\Theta_0+BK^{-1}B^\top)x
   =(\Theta_0+\Delta)x=\Theta_1x.
```

The displayed block matrix is alternating by construction.  Conversely,
every such elimination satisfies

```math
\Delta=BK^{-1}B^\top,
 \qquad
 \operatorname{rank}\Delta\le\dim W.
```

Thus no smaller interior space can generate a rank-$`2q`$ change.

∎

**Remark 9.11. What the theorem does not select.**

<a name="v12-rem-auxiliary-generation-not-selection"></a>

[Theorem 9.10](#v12-thm-minimal-schur-gauge-generation) is stronger than a volume
pullback: the target gauge appears only after an interior variable is
eliminated, and the number of auxiliary modes is minimal.  It is still an
*encoded auxiliary generation*, because $`B`$ and $`K`$ are built
from the desired $`\Delta`$.  It does not prove that a spatial–temporal
field complex, specified independently of $`\Theta_1`$, selects that gauge.
That field-selection problem remains the genuinely mixed open step.
The theorem's role is therefore to close the unrestricted algebraic
existence question and identify its sharp dimension.  Its block
factorization is an application of the standard real normal form for an
alternating matrix, not a claim of a new general realization theory.

The rectangular maps entering [Theorem 9.2](#v12-thm-rectangular-stokes-dirac) are

<a name="v12-eq-line-rectangular-maps"></a>
**Equation (132).**

```math

 \mathcal M_N^*=\operatorname{diag}(M_N^\top,M_N^\top),\qquad
 \mathcal K_N=
 \begin{bmatrix}0&-D_N\\-D_N&0\end{bmatrix},
```

while $`\mathcal C_Ne=u`$ and $`\mathcal Y_Ne=y`$ are the input and
output traces in [Equation (122)](#v12-eq-line-port-causality).

**Lemma 9.12. Polynomial endpoint reconstruction.**

<a name="v12-lem-polynomial-endpoint-reconstruction"></a>

For either $`a=0`$ or $`a=\ell`$, the single-field map

<a name="v12-eq-single-field-endpoint-reconstruction"></a>
**Equation (133).**

```math

 \begin{bmatrix}M_N^\top\\ r_a^\top\end{bmatrix}:
 E_N\longrightarrow F_N^*\oplus\mathbb R
```

is an isomorphism.  Consequently, with the complementary causality in
[Equation (122)](#v12-eq-line-port-causality), the two-field effort/input reconstruction
$`\mathcal R_N`$ below is an isomorphism.

*Proof.* Suppose first that $`M_N^\top e=0`$ and $`r_0^\top e=0`$.  Contracting
[Equation (118)](#v12-eq-rectangular-stokes-line) on both sides with $`e`$ gives
$`(r_\ell^\top e)^2-(r_0^\top e)^2=0`$, hence
$`r_\ell^\top e=0`$.  Left multiplication of the same identity by
$`e^\top`$ now gives

```math
(D_Ne)^\top M_N^\top=0,
```

because both trace terms and $`e^\top M_N`$ vanish.  Thus
$`M_ND_Ne=0`$.  The pairing matrix $`M_N`$ has full column rank, so
$`D_Ne=0`$.  Therefore $`e`$ is constant, and its zero endpoint trace
forces $`e=0`$.  The proof with $`r_\ell^\top e=0`$ is identical, with
the endpoints exchanged.  Each map in
[Equation (133)](#v12-eq-single-field-endpoint-reconstruction) is injective between
spaces of the same dimension $`N+1`$, hence bijective.  Taking the block
sum of the $`r_0`$ reconstruction for $`e_q`$ and the $`r_\ell`$
reconstruction for $`e_p`$ proves the last claim.

∎

The map

<a name="v12-eq-line-effort-input-reconstruction"></a>
**Equation (134).**

```math

 \mathcal R_N=
 \begin{bmatrix}\mathcal M_N^*\\\mathcal C_N\end{bmatrix}
```

is therefore invertible: each full effort polynomial is determined by
its $`N`$ active power coordinates and its prescribed complementary
endpoint trace.  The rectangular identity [Equation (118)](#v12-eq-rectangular-stokes-line) is
exactly [Equation (102)](#v12-eq-abstract-rectangular-stokes).  Therefore
[Theorem 9.2](#v12-thm-rectangular-stokes-dirac) proves, rather than assumes, the
finite-dimensional representation

<a name="v12-eq-reduced-open-line"></a>
**Equation (135).**

```math

 f=\mathcal J_N\epsilon+B_Nu,
 \qquad
 y=C_N\epsilon+N_Nu,
 \qquad
 C_N=B_N^\top,\quad
 \mathcal J_N^\top=-\mathcal J_N,\quad N_N^\top=-N_N.
```

The blocks are explicitly

<a name="v12-eq-line-square-block-recovery"></a>
**Equation (136).**

```math

 \begin{bmatrix}\mathcal J_N&B_N\end{bmatrix}
 =\mathcal K_N\mathcal R_N^{-1},\qquad
 \begin{bmatrix}B_N^\top&N_N\end{bmatrix}
 =\mathcal Y_N\mathcal R_N^{-1}.
```

This square representation is a consequence of the rectangular Stokes
identity; it is not the starting point.  For the closed line put
$`A_N=\mathcal J_NQ_N`$.  If $`f_\mathrm{phys}=\Delta x/h`$, the correct
oriented Dirac statement is

<a name="v12-eq-line-oriented-dirac-identification"></a>
**Equation (137).**

```math

 (-f_\mathrm{phys},u;\epsilon,y)\in\mathcal D_h,
 \qquad \mathcal J_h=\mathcal J_N,\quad B=B_N,\quad N=N_N,
```

in agreement with [Equation (107)](#v12-eq-dirac-step).

**Theorem 9.13. Mixed deformations of the stage-induced Gauss hierarchy.**

<a name="v12-thm-mixed-gauss-deformation"></a>

Assume that $`Q_N`$ and $`\mathcal J_N`$ are nonsingular.  Let
$`\Theta^\top=-\Theta`$ be the cap reduction of
[Equation (125)](#v12-eq-mixed-cylinder-lift).  For $`s\ge1`$, let
$`P_s,E_s,\widehat O_s`$, and $`R_s`$ be those of
[Theorem 8.8](#v12-thm-active-gauss-cell).  Where the displayed inverses exist, define
the stage-induced Gauss–$`s`$ cap gauge by

<a name="v12-eq-gauss-s-cap-gauge"></a>
**Equation (138).**

```math

 S^\mathrm{G}_{s,h}
 =\frac1h\Omega_N
 \left[
  \frac12E_s(hA_N)\widehat O_s(hA_N)^{-1}-I
 \right],
 \qquad \Omega_N=\mathcal J_N^{-1}.
```

Define its mixed deformation by

<a name="v12-eq-mixed-gauss-gauge"></a>
**Equation (139).**

```math

 \boxed{\quad
 S_{s,h}^{\Theta}=S^\mathrm{G}_{s,h}+h^{2s-1}\Theta.
 \quad}
```

For $`u=0`$, the step

<a name="v12-eq-mixed-gauss-step"></a>
**Equation (140).**

```math

 \Delta x
 =h\mathcal J_N
 \left(Q_N\bar x-S_{s,h}^{\Theta}\Delta x\right)
```

has the following properties whenever its matrices are invertible.

1. At $`\Theta=0`$, $`S^\mathrm{G}_{s,h}`$ is generated by elimination
    of the $`s`$ active polynomial stages in
    [Theorem 8.8](#v12-thm-active-gauss-cell).  The added mixed perturbation has the
    volume representative
    $`\mathcal T_N^*(h^{2s-1}\Theta)\mathcal T_N`$, and exact raw-cell
    elimination recovers the sum $`S_{s,h}^{\Theta}`$.

1. It preserves $`H_N`$ exactly, is self-adjoint, and has global order
    at least $`2s`$.  At $`\Theta=0`$, its
    propagator is the diagonal Padé map $`R_s(hA_N)`$, equivalently the
    $`s`$-stage Gauss method.

1. It is symplectic if and only if

   

   <a name="v12-eq-mixed-gauss-symplectic"></a>
   **Equation (141).**

   ```math

     [\mathcal J_N\Theta,A_N]=0.
   ```

1. If $`[\mathcal J_N\Theta,A_N]\ne0`$, then for all sufficiently
    small nonzero $`h`$ its propagator does not commute with $`A_N`$.  It
    therefore lies outside the scalar functional calculus of $`A_N`$.  In
    particular, it cannot be the result of applying any Runge–Kutta or
    collocation method with scalar coefficients to $`\dot x=A_Nx`$.

*Proof.* [Theorem 8.8](#v12-thm-active-gauss-cell), applied with
$`A=A_N`$ and $`\Omega=\Omega_N`$, proves that
$`S^\mathrm{G}_{s,h}`$ is alternating and generated by the interior Gauss
cell, with cap map $`R_s(hA_N)`$.  The mixed lift follows from
[Proposition 9.9](#v12-prop-mixed-cylinder-lift); its exact raw realization is
[Theorem 4.14](#v12-thm-trace-generated-condensation).

Put

```math
B_{s,\pm}=I+h\mathcal J_NS^\mathrm{G}_{s,h}
             \mathbin{\pm}\frac h2A_N,
 \qquad E_{s,h}=h^{2s}\mathcal J_N\Theta.
```

Eliminating the effort from [Equation (140)](#v12-eq-mixed-gauss-step) gives

<a name="v12-eq-mixed-gauss-map-general"></a>
**Equation (142).**

```math

 M_{s,\Theta,h}
 =(B_{s,-}+E_{s,h})^{-1}(B_{s,+}+E_{s,h}),
 \qquad B_{s,-}^{-1}B_{s,+}=R_s(hA_N).
```

The discrete chain rule proves exact energy preservation.  Both
$`S^\mathrm{G}_{s,h}`$ and $`h^{2s-1}\Theta`$ are odd in $`h`$, so reversing
$`h`$ and exchanging the caps inverts the map; hence the method is
self-adjoint.

The exact resolvent identity

```math
M_{s,\Theta,h}-R_s
 =(B_{s,-}+E_{s,h})^{-1}E_{s,h}(I-R_s)
```

is $`O(h^{2s+1})`$.  Since diagonal Padé has local error
$`O(h^{2s+1})`$, the deformation has local error
$`O(h^{2s+1})`$ and therefore global order at least $`2s`$.  The leading
coefficient and the unique exceptional value of $`\Theta`$ are computed
exactly in [Theorem 9.14](#v12-thm-unique-order-raising-gauge) below.

The symplectic criterion [Equation (76)](#v12-eq-locus), applied at each fixed $`h`$,
is automatic for $`S^\mathrm{G}_{s,h}`$, because its map is Gauss and
$`\mathcal J_NS^\mathrm{G}_{s,h}`$ is a rational function of $`A_N`$.
The remaining condition is
$`\Theta A_N=Q_N\mathcal J_N\Theta`$; multiplication by
$`\mathcal J_N`$ gives [Equation (141)](#v12-eq-mixed-gauss-symplectic).  Finally,
$`E_{s,h}=h^{2s}\mathcal J_N\Theta`$ and
$`I-R_s=-hA_N+O(h^2)`$ yield

<a name="v12-eq-mixed-map-commutator"></a>
**Equation (143).**

```math

 [M_{s,\Theta,h},A_N]
 =-h^{2s+1}[\mathcal J_N\Theta,A_N]A_N
  +O(h^{2s+2}).
```

Nonsingularity of $`A_N`$ proves the noncommutation claim.  Every scalar
rational, polynomial, or analytic function of $`A_N`$ commutes with it.

∎

**Theorem 9.14. Exact one-term order dichotomy.**

<a name="v12-thm-unique-order-raising-gauge"></a>

Under the assumptions of [Theorem 9.13](#v12-thm-mixed-gauss-deformation), set

<a name="v12-eq-pade-leading-coefficient"></a>
**Equation (144).**

```math

 c_s=(-1)^{s+1}\frac{(s!)^2}{(2s)!(2s+1)!}.
```

Within the one-term ansatz [Equation (139)](#v12-eq-mixed-gauss-gauge), let $`\Theta`$
be alternating, independent of $`h`$, and fixed as $`h\to0`$.  Then

<a name="v12-eq-mixed-gauss-leading-error"></a>
**Equation (145).**

```math

 M_{s,\Theta,h}-e^{hA_N}
 =h^{2s+1}\bigl(c_sA_N^{2s}-\mathcal J_N\Theta\bigr)A_N
  +O(h^{2s+2}).
```

Consequently the unique $`h`$-independent leading coefficient in this
ansatz that cancels the leading error is

<a name="v12-eq-unique-order-raising-gauge"></a>
**Equation (146).**

```math

 \boxed{\quad
 \Theta_s^\star=c_s\Omega_NA_N^{2s}.
 \quad}
```

It is alternating.  If $`\Theta\ne\Theta_s^\star`$, the method has
exact global order $`2s`$; for $`\Theta=\Theta_s^\star`$, it has global
order exactly $`2s+2`$.  More precisely,

<a name="v12-eq-raised-branch-leading-error"></a>
**Equation (147).**

```math

 M_{s,\Theta_s^\star,h}-e^{hA_N}
 =d_sh^{2s+3}A_N^{2s+3}+O(h^{2s+4}),
```

where

<a name="v12-eq-second-cayley-coefficient"></a>
**Equation (148).**

```math

 d_1=-\frac1{720},\qquad
 d_s=-c_s\frac{(s+1)(2s-3)}{6(2s-1)(2s+3)}\quad(s\ge2).
```

In particular, every noncommuting coefficient in the one-term ansatz has
exactly order $`2s`$.  The cap coefficient
$`\Theta_s^\star`$ depends on the constitutive generator, but by
[Proposition 9.9](#v12-prop-mixed-cylinder-lift) it possesses a pure-mixed cylinder
realization.  For $`s=1`$, the order-raised method is identically the
two-stage Gauss method, not merely a fourth-order approximation to it.

*Proof.* Write $`P_s(z)=\sum_{k=0}^sa_kz^k`$, with $`a_k`$ given in
[Equation (97)](#v12-eq-pade-polynomial).  For $`r\ge0`$, coefficient extraction and
the beta integral give

<a name="v12-eq-pade-residual-coefficients"></a>
**Equation (149).**

```math
\begin{aligned}
g_{s,r}
 &:=[z^{2s+1+r}]\bigl(P_s(-z)e^z-P_s(z)\bigr)\\
 &=\frac{s!}{(2s)!r!}
   \int_0^1(1-t)^r
   \sum_{k=0}^s\frac{(-1)^kt^{2s-k}}{k!(s-k)!}\,\,\mathrm{d} t\\
 &=(-1)^s\frac{s!(s+r)!}
 {(2s)!r!(2s+r+1)!}.

\end{aligned}
```

Indeed,

```math
\sum_{k=0}^s\frac{(-1)^kt^{2s-k}}{k!(s-k)!}
 =\frac{(-1)^s}{s!}t^s(1-t)^s,
```

which proves [Equation (149)](#v12-eq-pade-residual-coefficients).  Since

```math
a_1=\frac12,\qquad
 a_2=\frac{s-1}{4(2s-1)},
```

division by
$`P_s(-z)=1-a_1z+a_2z^2+O(z^3)`$, using $`r=0,1,2`$, yields

<a name="v12-eq-pade-three-term-error"></a>
**Equation (150).**

```math

 R_s(z)-e^z
 =c_sz^{2s+1}\bigl(1+z+\rho_sz^2+O(z^3)\bigr),
 \qquad
 \rho_s=\frac{8s^2+10s-5}{4(2s-1)(2s+3)}.
```

In particular, its first term is

<a name="v12-eq-pade-leading-error"></a>
**Equation (151).**

```math

 R_s(z)-e^z=c_sz^{2s+1}+O(z^{2s+2}).
```

For the mixed perturbation, the exact resolvent identity in the proof of
[Theorem 9.13](#v12-thm-mixed-gauss-deformation), together with
$`E_{s,h}=h^{2s}\mathcal J_N\Theta`$, gives

```math
M_{s,\Theta,h}-R_s(hA_N)
 =-h^{2s+1}\mathcal J_N\Theta A_N+O(h^{2s+2}).
```

Adding [Equation (151)](#v12-eq-pade-leading-error) proves
[Equation (145)](#v12-eq-mixed-gauss-leading-error).  Since $`A_N`$ and
$`\mathcal J_N`$ are invertible, its coefficient vanishes if and only if
$`\Theta=c_s\Omega_NA_N^{2s}`$.  Moreover,

```math
(\Omega_NA_N^{2s})^\top=-\Omega_NA_N^{2s},
```

so this unique coefficient is an admissible two-form.  Every other choice
has nonzero local error at degree $`2s+1`$, hence exact global order
$`2s`$.  Cancellation gives global order at least $`2s+1`$; self-adjoint
one-step methods have even order, so the order is at least $`2s+2`$.
If the gauge is noncommuting, it cannot equal
$`\Theta_s^\star`$, because
$`\mathcal J_N\Theta_s^\star=c_sA_N^{2s}`$ commutes with $`A_N`$.
Finally, $`S^\mathrm{G}_{1,h}=0`$, while
$`h\Theta_1^\star=(h/12)\Omega_NA_N^2=S^\mathrm{G}_{2,h}`$; the two maps
are therefore identical.

It remains to determine the exact order on the raised branch.  Introduce
the even common blocks

```math
C_s(z)=\tfrac12E_s(z)\widehat O_s(z)^{-1},\qquad
 C_\mathrm{ex}(z)=\frac z2\coth\!\left(\frac z2\right).
```

The scalar identity

<a name="v12-eq-cayley-common-block-difference"></a>
**Equation (152).**

```math

 C_\mathrm{ex}(z)-C_s(z)
 =\frac{z\bigl(R_s(z)-e^z\bigr)}
 {(e^z-1)(R_s(z)-1)}
```

transfers the Padé residual to the gauge coordinate.  For $`s\ge2`$,
[Equation (150)](#v12-eq-pade-three-term-error) and

```math
\frac{z^2}{(e^z-1)(R_s(z)-1)}
 =\left(\frac z{e^z-1}\right)^2+O(z^{2s})
 =1-z+\frac5{12}z^2+O(z^3)
```

give

<a name="v12-eq-two-cayley-coefficients"></a>
**Equation (153).**

```math

 C_\mathrm{ex}(z)-C_s(z)
 =c_sz^{2s}+d_sz^{2s+2}+O(z^{2s+4}),
```

where

```math
d_s=c_s\left(\rho_s-\frac7{12}\right)
 =-c_s\frac{(s+1)(2s-3)}{6(2s-1)(2s+3)}.
```

For $`s=1`$, $`C_1=1`$, so the Bernoulli expansion of
$`C_\mathrm{ex}`$ gives $`d_1=-1/720`$ directly.  Thus the corrected common
block $`C_s+c_sz^{2s}`$ differs from $`C_\mathrm{ex}`$ by
$`-d_sz^{2s+2}+O(z^{2s+4})`$.  The same resolvent identity used above,
now with $`z=hA_N`$, proves
[Equation (147)](#v12-eq-raised-branch-leading-error).  Every $`d_s`$ is nonzero, and
$`A_N`$ is invertible, so the local order is exactly $`2s+3`$ and the
global order exactly $`2s+2`$.

∎

**Corollary 9.15. Order, centrality, and nonseparability at one term.**

<a name="v12-cor-one-term-order-centrality"></a>

Within [Equation (139)](#v12-eq-mixed-gauss-gauge), the following are equivalent:

1. the global order is greater than $`2s`$;

1. $`\Theta=\Theta_s^\star`$;

1. the global order is exactly $`2s+2`$.

This unique higher-order branch lies in the symplectic centralizer and is
a scalar rational function of $`A_N`$.  Conversely, every coefficient
with $`[\mathcal J_N\Theta,A_N]\ne0`$ has exactly order $`2s`$, is
nonsymplectic for small nonzero $`h`$, and lies outside the scalar
functional calculus of $`A_N`$.

*Proof.* The order equivalences follow from
[Theorem 9.14](#v12-thm-unique-order-raising-gauge).  Moreover,
$`\mathcal J_N\Theta_s^\star=c_sA_N^{2s}`$, so the exceptional
coefficient commutes with $`A_N`$ and satisfies
[Equation (141)](#v12-eq-mixed-gauss-symplectic).  The remaining statements follow from
[Equation (143)](#v12-eq-mixed-map-commutator) and
[Theorem 9.13](#v12-thm-mixed-gauss-deformation).

∎

**Theorem 9.16. Two-term high-order noncommuting family.**

<a name="v12-thm-two-term-gauge-hierarchy"></a>

Let $`\Xi^\top=-\Xi`$ be independent of $`h`$, and set

<a name="v12-eq-two-term-gauge"></a>
**Equation (154).**

```math

 S_{s,h}^{\star,\Xi}
 =S^\mathrm{G}_{s,h}
 +h^{2s-1}\Theta_s^\star+h^{2s+1}\Xi.
```

The associated method is energy preserving, self-adjoint, and has global
order at least $`2s+2`$.  Its leading error is

<a name="v12-eq-two-term-leading-error"></a>
**Equation (155).**

```math

 M_{s,\star,\Xi,h}-e^{hA_N}
 =h^{2s+3}\bigl(d_sA_N^{2s+2}-\mathcal J_N\Xi\bigr)A_N
 +O(h^{2s+4}).
```

Hence it has exact order $`2s+2`$ unless

<a name="v12-eq-second-central-corrector"></a>
**Equation (156).**

```math

 \Xi=\Xi_s^\star:=d_s\Omega_NA_N^{2s+2},
```

in which case its order is at least $`2s+4`$.  If
$`[\mathcal J_N\Xi,A_N]\ne0`$, the method has exact order $`2s+2`$
while being nonsymplectic and outside the scalar functional calculus of
$`A_N`$ for all sufficiently small nonzero $`h`$.

*Proof.* The first correction gives [Equation (147)](#v12-eq-raised-branch-leading-error).  The
additional cap term changes the common block by
$`h^{2s+2}\mathcal J_N\Xi`$.  The resolvent identity therefore contributes
$`-h^{2s+3}\mathcal J_N\Xi A_N+O(h^{2s+4})`$, proving
[Equation (155)](#v12-eq-two-term-leading-error).  Invertibility of $`A_N`$ and
$`\mathcal J_N`$ gives the unique second corrector.  Cancellation and
self-adjointness imply order at least $`2s+4`$.  Finally,

```math
[M_{s,\star,\Xi,h},A_N]
 =-h^{2s+3}[\mathcal J_N\Xi,A_N]A_N+O(h^{2s+4}),
```

and the symplectic criterion reduces to the same commutator.

∎

**Corollary 9.17. Recursive central correction hierarchy.**

<a name="v12-cor-recursive-gauge-hierarchy"></a>

Near the origin, define the even Taylor coefficients by

<a name="v12-eq-recursive-cayley-coefficients"></a>
**Equation (157).**

```math

 C_\mathrm{ex}(z)-C_s(z)
 =\sum_{j\ge0}\alpha_{s,j}z^{2s+2j},
 \qquad \alpha_{s,0}=c_s,\quad\alpha_{s,1}=d_s.
```

For $`r\ge1`$, the truncated cap gauge

<a name="v12-eq-recursive-central-gauge"></a>
**Equation (158).**

```math

 S_{s,h}^{[r]}
 =S^\mathrm{G}_{s,h}
 +\sum_{j=0}^{r-1}\alpha_{s,j}
 h^{2s-1+2j}\Omega_NA_N^{2s+2j}
```

is energy preserving, self-adjoint, symplectic, and has global order at
least $`2s+2r`$.  If $`\alpha_{s,r}\ne0`$, that order is exact.  Within
the common disk of analyticity, the full convergent series is the
exact-flow gauge.  At level $`r`$, an arbitrary alternating addition
$`h^{2s-1+2r}\Xi`$ preserves order at least $`2s+2r`$; it introduces
noncommutation whenever $`[\mathcal J_N\Xi,A_N]\ne0`$.

*Proof.* Multiplication by $`h\mathcal J_N`$ turns the $`j`$-th cap correction in
[Equation (158)](#v12-eq-recursive-central-gauge) into

```math
\alpha_{s,j}(hA_N)^{2s+2j}.
```

Thus the common block agrees with
$`C_\mathrm{ex}(hA_N)`$ through degree $`2s+2r-2`$, and the Cayley map has
local error $`O(h^{2s+2r+1})`$.  Every central correction is alternating
and commutes with $`A_N`$, proving the structural claims.  The last
statement follows from the same one-order shift between a cap-gauge term
and its map perturbation used in the preceding theorem.

∎

**Map of the central notation.**
The construction passes from the diagonal Padé polynomial to a Cayley
coordinate, adds the central Taylor correction, clears the coordinate
denominator, and only then forms the intrinsic propagator:

```math
P_s
 \ \longrightarrow\ (E_s,\widehat O_s)
 \ \longrightarrow\ C_s=\frac{E_s}{2\widehat O_s}
 \ \xrightarrow{\ +\delta_{s,r}\ }\ C_{s,r}
 \ \longrightarrow\ \Pi_{s,r}
 \ \longrightarrow\
 R_{s,r}=\frac{\Pi_{s,r}(z)}{\Pi_{s,r}(-z)}.
```

Here $`P_s(z)/P_s(-z)`$ is the Gauss–$`s`$ map,
$`\delta_{s,r}`$ is the truncated central correction,
$`C_{s,r}`$ is a chart coordinate and may diverge, whereas
$`\Pi_{s,r}`$ is the denominator-cleared polynomial whose coprime quotient
determines the true poles of $`R_{s,r}`$.

**Theorem 9.18. Central rational hierarchy.**

<a name="v12-thm-central-rational-hierarchy"></a>

Let the even functions $`C_s`$, $`C_\mathrm{ex}`$, and the coefficients
$`\alpha_{s,j}`$ be defined by
[Equation (157)](#v12-eq-recursive-cayley-coefficients).  For $`r\ge1`$, set

<a name="v12-eq-central-delta-pi"></a>
**Equation (159).**

```math

 \delta_{s,r}(z)=\sum_{j=0}^{r-1}\alpha_{s,j}z^{2s+2j},
 \qquad
 \Pi_{s,r}(z)=P_s(z)+2\widehat O_s(z)\delta_{s,r}(z).
```

Then the central propagator has the two equivalent representations

<a name="v12-eq-central-rational-family"></a>
**Equation (160).**

```math

 C_{s,r}(z)=C_s(z)+\delta_{s,r}(z),
 \qquad
 R_{s,r}(z)=
 \bigl(C_{s,r}(z)-z/2\bigr)^{-1}
 \bigl(C_{s,r}(z)+z/2\bigr)
 =\Pi_{s,r}(-z)^{-1}\Pi_{s,r}(z).
```

The last expression is understood after cancellation of common scalar
polynomial factors.  It has the following properties.

1. $`R_{s,r}(-z)=R_{s,r}(z)^{-1}`$.  Hence the induced one-step
    map is self-adjoint.

1. If $`A_N^\top Q_N+Q_NA_N=0`$, then
    $`R_{s,r}(hA_N)^\top Q_NR_{s,r}(hA_N)=Q_N`$.  Thus the map preserves
    the quadratic Hamiltonian exactly.

1. For fixed $`A_N`$,

   

   <a name="v12-eq-central-rational-order"></a>
   **Equation (161).**

   ```math

     R_{s,r}(hA_N)-e^{hA_N}=O(h^{2s+2r+1}).
   ```

   
    In the terminology of [Definition 1.1](#v12-def-matrix-order), the local matrix defect
    is of order at least $`2s+2r+1`$, so the global operator order is at
    least $`2s+2r`$.  It is
    exactly $`2s+2r`$ when $`\alpha_{s,r}\ne0`$ and
    $`A_N^{2s+2r+1}\ne0`$.

1. Let $`N_{s,r}/D_{s,r}`$ be the coprime scalar reduction of
    $`\Pi_{s,r}(z)/\Pi_{s,r}(-z)`$.  The true finite-step poles of the
    matrix propagator are exactly the failures of invertibility of
    $`D_{s,r}(hA_N)`$.  Zeros of $`\widehat O_s`$ are singularities of the
    coordinate $`C_s`$, not additional poles of $`R_{s,r}`$: the polynomial
    representation extends the propagator through them, with $`R_{s,r}=1`$.

1. The reduced scalar propagator has no pole on the imaginary axis.
    For every real $`\xi`$, its continuous rational extension satisfies

   

   <a name="v12-eq-central-rational-imaginary-axis"></a>
   **Equation (162).**

   ```math

     |R_{s,r}(i\xi)|=1,
     \qquad
     \arg R_{s,r}(i\xi)
     =2\operatorname{atan2}\!\left(
       \Im\Pi_{s,r}(i\xi),\Re\Pi_{s,r}(i\xi)
      \right)\pmod{2\pi}.
   ```

   
    In a finite $`C_{s,r}`$-chart this becomes

   

   <a name="v12-eq-central-rational-phase"></a>
   **Equation (163).**

   ```math

     \arg R_{s,r}(i\xi)
     =2\operatorname{atan2}\!\left(\xi,2C_{s,r}(i\xi)\right)
     \quad\pmod{2\pi}.
   ```

The cap gauge generating $`R_{s,r}(hA_N)`$ is exactly
[Equation (158)](#v12-eq-recursive-central-gauge).

*Proof.* Since

```math
C_{s,r}(z)\mathbin{\pm}\frac z2
 =\frac{P_s(\pm z)+2\widehat O_s(z)\delta_{s,r}(z)}
        {2\widehat O_s(z)},
```

the coordinate denominator cancels and gives the polynomial quotient in
[Equation (160)](#v12-eq-central-rational-family).  This proves in particular that poles
of $`C_s`$ are not automatically poles of the propagator.  At a zero of
$`\widehat O_s`$, the coprimality of the diagonal Padé pair gives
$`P_s(z)=P_s(-z)=E_s(z)\ne0`$; hence the continued value is $`1`$.

The function $`C_{s,r}`$ is even, which gives the reversibility identity.
From $`A_N^\top Q_N=-Q_NA_N`$, every even rational function satisfies
$`C_{s,r}(hA_N)^\top Q_N=Q_NC_{s,r}(hA_N)`$.  Expanding the two quadratic
products then gives

```math
(C_{s,r}+hA_N/2)^\top Q_N(C_{s,r}+hA_N/2)
 =(C_{s,r}-hA_N/2)^\top Q_N(C_{s,r}-hA_N/2),
```

which proves exact $`Q_N`$-orthogonality.  By construction,
$`C_\mathrm{ex}-C_{s,r}=O(z^{2s+2r})`$.  The Cayley resolvent identity adds
one power of $`z`$, proving [Equation (161)](#v12-eq-central-rational-order); its first
coefficient is nonzero under the stated hypotheses.  Cancelling the
greatest common polynomial factor gives the asserted intrinsic pole
condition.

All coefficients of $`\Pi_{s,r}`$ are real.  Therefore
$`\Pi_{s,r}(-i\xi)=\overline{\Pi_{s,r}(i\xi)}`$.  If either factor
vanishes on the imaginary axis, so does the other, with the same
multiplicity; that factor disappears in the coprime reduction.  The
reduced denominator consequently has no imaginary-axis zero, and the
quotient has unit modulus.  Its argument is twice the argument of
$`\Pi_{s,r}(i\xi)`$, proving
[Equation (162)](#v12-eq-central-rational-imaginary-axis).  Where the $`C_{s,r}`$-chart
is finite, the same calculation applied to
$`C_{s,r}(i\xi)+i\xi/2`$ gives
[Equation (163)](#v12-eq-central-rational-phase).  Finally, multiplication of
[Equation (158)](#v12-eq-recursive-central-gauge) by $`h\mathcal J_N`$ produces exactly
the added even powers in $`C_{s,r}(hA_N)`$.

∎

**Remark 9.19. The Gauss–1 coincidence is exceptional.**

<a name="v12-rem-gauss-coincidence-exception"></a>

For $`s=1,r=1`$, $`C_{1,1}(z)=1+z^2/12=C_2(z)`$, so the corrected
one-stage map is exactly Gauss–2.  There is no general identity
``corrected Gauss–$`s`$ equals Gauss–$`(s+1)`$.”  Already

```math
C_{2,1}(z)=1+\frac{z^2}{12}-\frac{z^4}{720},
 \qquad
 C_3(z)=\frac{1+z^2/10}{1+z^2/60},
```

which agree only through the required Taylor degree and define distinct
rational maps.

**Proposition 9.20. Coordinate resonances and loss of $`A`$-stability.**

<a name="v12-prop-rational-resonance-stability"></a>

The distinction between the gauge chart and the reduced propagator is
already visible in the first nontrivial members.

1. For $`s=3`$,

   

   
   ```math
   C_3(z)=6\frac{z^2+10}{z^2+60}.
   ```

   
    At $`z=\pm i\sqrt{60}`$, the coordinate diverges, whereas for every
    $`r\ge1`$,
    $`\Pi_{3,r}(z)=\Pi_{3,r}(-z)=-5`$ and the continued propagator is
    $`R_{3,r}(z)=1`$.

1. For $`s=2,r=1`$,

   

   <a name="v12-eq-pi-21"></a>
   **Equation (164).**

   ```math

     \Pi_{2,1}(z)=1+\frac z2+\frac{z^2}{12}-\frac{z^4}{720}.
   ```

   
    It has a positive real zero
    $`\rho_\star=10.1275149065\ldots`$.  Hence the reduced denominator
    $`\Pi_{2,1}(-z)`$ has a genuine pole at
    $`z=-\rho_\star`$ in the open left half-plane.  The sixth-order map
    $`R_{2,1}`$ is therefore not $`A`$-stable.

*Proof.* For $`s=3`$,
$`\widehat O_3(z)=1/2+z^2/120`$; at its two zeros the correction in
[Equation (159)](#v12-eq-central-delta-pi) vanishes and direct substitution in $`P_3`$
gives $`-5`$ for both signs.  For $`s=2`$,
$`2\widehat O_2=1`$, $`\alpha_{2,0}=-1/720`$, and
[Equation (164)](#v12-eq-pi-21) follows.  Direct polynomial root isolation gives
$`\Pi_{2,1}(10.1275)>0>\Pi_{2,1}(10.1276)`$, and Newton refinement gives
$`\rho_\star`$.  Since
$`\Pi_{2,1}(z)-\Pi_{2,1}(-z)=z`$, the numerator is nonzero at
$`-\rho_\star`$, so the pole does not cancel.

∎

**Remark 9.21. Scope of the stability statement.**

The absence of imaginary-axis poles in
[Theorem 9.18](#v12-thm-central-rational-hierarchy) is the property relevant to the
positive-energy Hamiltonian line.  It does not imply left-half-plane
stability for dissipative, parabolic, or split extensions.  Root
localization of $`\Pi_{s,r}`$ is therefore a separate design constraint,
not a consequence of high order or exact quadratic conservation.

**Remark 9.22. Relation with the exact-flow gauge.**

<a name="v12-rem-order-gauge-bernoulli"></a>

Equivalently,

```math
S_h^\mathrm{exact}-S^\mathrm{G}_{s,h}
 =h^{2s-1}\Theta_s^\star+O(h^{2s+1}).
```

Thus $`c_s`$ is the leading coefficient of the *difference*
between the exact-flow gauge and the Gauss–$`s`$ gauge.  It is not, in
general, the raw Bernoulli coefficient of
$`(z/2)\coth(z/2)-1`$.  For example, at $`s=3`$ the exact $`z^6`$
coefficient is $`1/30240`$, the Gauss–3 gauge already contributes
$`1/43200`$, and their difference is
$`c_3=1/100800`$.

**Proposition 9.23. Stieltjes representation and Padé identification.**

<a name="v12-prop-stieltjes-pade-coordinate"></a>

With $`u=z^2`$, define

```math
F(u)=\frac{C_\mathrm{ex}(\sqrt u)-1}{u}.
```

Then

<a name="v12-eq-stieltjes-coth"></a>
**Equation (165).**

```math

 F(u)=2\sum_{k\ge1}\frac1{u+4\pi^2k^2}
 =\int_{[0,1/(4\pi^2)]}\frac{\,\mathrm{d}\nu(\lambda)}{1+u\lambda},
 \qquad
 \nu=\sum_{k\ge1}\frac{2}{4\pi^2k^2}
       \delta_{1/(4\pi^2k^2)}.
```

Thus $`F`$ is a Stieltjes function.  Moreover,
$`C_s(\sqrt u)`$ is the Padé approximant of type

```math
\left[\left\lfloor\frac s2\right\rfloor
 \Big/\left\lfloor\frac{s-1}{2}\right\rfloor\right]
```

to $`C_\mathrm{ex}(\sqrt u)`$ at $`u=0`$.  Equivalently,
$`(C_s(\sqrt u)-1)/u`$ is of type $`[m-1/m-1]`$ when $`s=2m`$, and
of type $`[m-1/m]`$ when $`s=2m+1`$.

*Proof.* The Mittag–Leffler expansion of $`\coth`$ gives the first equality in
[Equation (165)](#v12-eq-stieltjes-coth); rewriting each summand with
$`\lambda_k=(4\pi^2k^2)^{-1}`$ gives the positive discrete measure.
The numerator and denominator degrees of
$`C_s=E_s/(2\widehat O_s)`$, expressed in $`u`$, are respectively
$`\lfloor s/2\rfloor`$ and $`\lfloor(s-1)/2\rfloor`$.  Finally,
[Equation (152)](#v12-eq-cayley-common-block-difference) gives
$`C_\mathrm{ex}(\sqrt u)-C_s(\sqrt u)=O(u^s)`$, which is exactly the Padé
matching condition for those degrees.  Division of the constant-free
difference by $`u`$ gives the two stated types.  This places the coordinate
hierarchy within classical Stieltjes–Padé theory
[[36]](#v12-ref-wynn1968).

∎

**Theorem 9.24. Positive interior generation of the Gauss coordinates.**

<a name="v12-thm-stieltjes-interior-generation"></a>

Let

```math
F_s(u)=\frac{C_s(\sqrt u)-1}{u},
 \qquad q_s=\left\lfloor\frac{s-1}{2}\right\rfloor.
```

There exist $`a_s\ge0`$, positive numbers
$`w_{s,k},\lambda_{s,k}>0`$, $`1\le k\le q_s`$, and hence

```math
L_s=\operatorname{diag}(\lambda_{s,1},\ldots,\lambda_{s,q_s})\succ0,
 \qquad
 b_s=(\sqrt{w_{s,1}},\ldots,\sqrt{w_{s,q_s}})^\top,
```

such that

<a name="v12-eq-finite-stieltjes-realization"></a>
**Equation (166).**

```math

 F_s(u)=a_s+\sum_{k=1}^{q_s}
 \frac{w_{s,k}}{1+u\lambda_{s,k}}
 =a_s+b_s^\top(I+uL_s)^{-1}b_s.
```

The sum is absent when $`q_s=0`$.  Consequently $`C_s`$ is generated by
the positive interior system

<a name="v12-eq-positive-coordinate-interior-system"></a>
**Equation (167).**

```math

 (I+uL_s)\eta=b_sx,
 \qquad
 p=a_sx+b_s^\top\eta,
 \qquad
 C_s(\sqrt u)x=x+up.
```

After Lanczos tridiagonalization with starting vector $`b_s`$, the same
realization has a positive Jacobi matrix and repeated scalar Schur
elimination gives a finite Stieltjes continued fraction.

Here *positive Jacobi-chain realization* means precisely that
$`L_s\succ0`$, the residues $`w_{s,k}`$ are positive, the associated
rational Stieltjes function is nonnegative on $`u\ge0`$, and the Jacobi
off-diagonal coefficients can be chosen positive.  It does not assert
that the complete augmented matrix is positive in the energy metric of
the original port-Hamiltonian system.

This positive mechanism does not generate the entire corrected hierarchy.
Already

<a name="v12-eq-corrected-stieltjes-obstruction"></a>
**Equation (168).**

```math

 \frac{C_{2,1}(\sqrt u)-1}{u}
 =\frac1{12}-\frac{u}{720}<0
 \qquad (u>60),
```

so $`C_{2,1}`$ admits no realization of the positive form
[Equation (166)](#v12-eq-finite-stieltjes-realization).

*Proof.* By [Proposition 9.23](#v12-prop-stieltjes-pade-coordinate), $`F_s`$ is a diagonal or
subdiagonal Padé approximant of the Stieltjes function $`F`$.  The
classical Stieltjes–Padé theorem gives real negative simple poles and
positive residues; a diagonal member may additionally have a nonnegative
constant term at infinity [[36]](#v12-ref-wynn1968).  Rewriting every pole as
$`-1/\lambda_{s,k}`$ gives [Equation (166)](#v12-eq-finite-stieltjes-realization).
This is the finite positive-measure realization of a rational Stieltjes
function.

Eliminating $`\eta`$ in
[Equation (167)](#v12-eq-positive-coordinate-interior-system) gives $`p=F_s(u)x`$, hence
the last identity.  Orthogonal Lanczos reduction of $`L_s`$ with starting
vector $`b_s`$ preserves the scalar resolvent and produces a symmetric
tridiagonal Jacobi matrix; successive block elimination of its last
coordinate yields the Stieltjes continued fraction [[35]](#v12-ref-wall1948).
The corresponding positive-network interpretation belongs to classical
network synthesis [[1]](#v12-ref-basubose1983).  Finally, the right-hand side
of [Equation (166)](#v12-eq-finite-stieltjes-realization) is nonnegative for every
$`u\ge0`$, whereas direct substitution of $`C_{2,1}`$ gives
[Equation (168)](#v12-eq-corrected-stieltjes-obstruction).  The claimed obstruction
follows.

∎

**Remark 9.25. Two distinct notions of interior generation.**

<a name="v12-rem-two-interior-generations"></a>

[Theorem 9.10](#v12-thm-minimal-schur-gauge-generation) uses a nondegenerate alternating
interior block and can generate every finite-rank cap two-form, but it may
be indefinite and is designed from the target gauge.  By contrast,
[Theorem 9.24](#v12-thm-stieltjes-interior-generation) uses a positive constitutive
interior solve and gives the canonical Jacobi chain behind Gauss–Padé.
The failure in [Equation (168)](#v12-eq-corrected-stieltjes-obstruction) shows that the
Taylor-corrected central hierarchy cannot simply be reinterpreted as the
same passive chain with extra stages.  A field-generated realization of
those corrections must either use a different rational hierarchy or admit
indefinite/nonpassive internal structure.

**Remark 9.26. A coefficientwise sign conjecture.**

<a name="v12-rem-alpha-sign-conjecture"></a>

Exact arithmetic for $`1\le s\le12`$ and $`0\le j\le12`$ gives

<a name="v12-eq-alpha-sign-conjecture"></a>
**Equation (169).**

```math

 (-1)^{s+j+1}\alpha_{s,j}>0.
```

If proved for all $`s,j`$, this would remove the non-annulation hypothesis
from [Corollary 9.17](#v12-cor-recursive-gauge-hierarchy), [Theorem 9.18](#v12-thm-central-rational-hierarchy).
The Stieltjes representation makes [Equation (169)](#v12-eq-alpha-sign-conjecture)
natural, but standard real-axis Padé error inequalities do not by
themselves establish this coefficientwise statement.  We therefore retain
it as a conjecture, and every exact-order claim continues to state the
hypothesis $`\alpha_{s,r}\ne0`$.

**Corollary 9.27. The Gauss–2 member.**

<a name="v12-cor-mixed-gauss-two"></a>

For $`s=2`$,

<a name="v12-eq-mixed-gauss-two-gauge"></a>
**Equation (170).**

```math

 S_{2,h}^{\Theta}
 =\frac h{12}Q_N\mathcal J_NQ_N+h^3\Theta.
```

Writing

```math
D_2=I-\frac h2A_N+\frac{h^2}{12}A_N^2,
 \qquad
 N_2=I+\frac h2A_N+\frac{h^2}{12}A_N^2,
 \qquad E_h=h^4\mathcal J_N\Theta,
```

the propagator is

<a name="v12-eq-mixed-gauss-map"></a>
**Equation (171).**

```math

 M_{2,\Theta,h}=(D_2+E_h)^{-1}(N_2+E_h).
```

It is the fourth-order example used below.

*Proof.* For $`P_2(z)=I+z/2+z^2/12`$,
$`E_2(z)=I+z^2/12`$ and $`\widehat O_2(z)=I/2`$.  Substitution in
[Equation (138)](#v12-eq-gauss-s-cap-gauge), [Equation (142)](#v12-eq-mixed-gauss-map-general) gives the stated
formulas.

∎

**Remark 9.28. What is and is not new in the line example.**

The case $`\Theta=0`$ is only Gauss–$`s`$ applied to the 2012 spatial
reduction and lies within the architecture of [[18]](#v12-ref-kotyczka2019).
The difference is the lifted $`Q_{tx}`$ term.  It preserves the proved
order and power but, under the commutator hypothesis, destroys both the
scalar-coefficient Runge–Kutta method-of-lines representation and
symplecticity.  Thus the construction does not relabel Gauss as
space–time; it contains the whole Gauss hierarchy as the separable
zero-mixed-gauge family, realizes every deformation through the raw
cylinder system, and gives a falsifiable algebraic criterion for leaving
the broader scalar functional-calculus subclass.  The one-term
classification also exposes a sharp tradeoff: its unique higher-order
coefficient is central.  The two-term family in
[Theorem 9.16](#v12-thm-two-term-gauge-hierarchy) removes that apparent obstruction by
retaining order $`2s+2`$ with an arbitrary noncommuting next coefficient.
High-order
energy-preserving methods exist in other forms
[[25]](#v12-ref-miyatake2016), [[10]](#v12-ref-eidnes2022), [[16]](#v12-ref-kemmochi2023); the specific contribution
here is the mixed cylinder realization and the nonseparability criterion.
Numerical advantage still requires a selector; below, modal dispersion
provides one such test.

**Energy currents and multisymplectic forms.**

Polarization of the quadratic energy law gives a symmetric bilinear
current.  A multisymplectic form formula is instead an alternating
two-form evaluated on two variational solutions
[[2]](#v12-ref-bridges2001), [[27]](#v12-ref-reich2000).  Gauss collocation may preserve both because
the Cooper condition preserves quadratic invariants, but the two tensor
parities must not be identified.  In the present theory, symmetric power
balance holds throughout the compatible affine family, while the
alternating temporal symplectic law holds exactly on $`\mathcal L`$.

<a name="v12-section-10"></a>

## 10. Consequences and discriminating computational validation

### 10.1. Classification tables

<a name="v12-tab-gauge"></a>

**Table 1. Properties of a fixed, $`h`$-independent one-step linear gauge
for nonsingular $`J,Q`$.  Step-dependent odd gauges such as
[Equation (139)](#v12-eq-mixed-gauss-gauge) are not covered by the order column.**

| Gauge | Exact energy | Symplectic | Order |
| --- | --- | --- | --- |
| $`S=0`$ | yes | yes | two |
| $`S\in\mathcal{L}\setminus\{0\}`$ | yes | yes | one |
| $`S\notin\mathcal{L}`$ | yes | no | one |

<a name="v12-tab-nonlinear"></a>

**Table 2. Nonlinear two-point methods: generic structural properties.**

| Method | Exact prescribed energy | Symplectic generically | Affine-covariant construction |
| --- | --- | --- | --- |
| AVF discrete gradient | yes | no | yes |
| Gonzalez discrete gradient | yes | no | no, metric-dependent |
| Implicit midpoint | only for quadratic $`H`$ | yes | yes |

### 10.2. Verification on a stacked space–time polynomial cell

We first test the stacked, nonsurjective differential rather than the
genuine total-complex map or a time-stepping reduction.  Take

```math
E_t=\mathcal{P}_3,\qquad E_x=\mathcal{P}_4,
```

with monomial coefficients on $`[0,1]^2`$, exact $`L^2`$ product
pairings, and

```math
D_\mathrm{st}z=
 \bigl((\partial_t\otimes I)z,(I\otimes\partial_x)z\bigr).
```

Then $`\dim E_\mathrm{st}=20`$, $`\dim F_\mathrm{st}=31`$, and the computed
rank is $`19`$, exactly as predicted by
[Equation (31)](#v12-eq-space-time-rank); hence
$`\dim\mathcal G_\mathrm{st}=19\cdot18/2=171`$.

Starting from the exact product representative, we add unit-Frobenius
gauges supported respectively in the temporal, spatial, mixed, and
temporal–mixed cross sectors of [Equation (37)](#v12-eq-wedge-sector-split).  The matrix
residual

<a name="v12-eq-st-residual"></a>
**Equation (172).**

```math

 \varepsilon_\mathrm{st}
 =\left\lVert M_\mathrm{st}D_\mathrm{st}
 +(M_\mathrm{st}D_\mathrm{st})^\top-B_\mathrm{st}\right\rVert_F
```

and the alternating self-contraction are reported in
[Table 3](#v12-tab-space-time-residuals).

<a name="v12-tab-space-time-residuals"></a>

**Table 3. Compatibility residuals on the $`\mathcal{P}_3\otimes\mathcal{P}_4`$
space–time cell.**

| Gauge sector | $`\varepsilon_\mathrm{st}`$ | $`\vert \omega(z,z)\vert `$ |
| --- | --- | --- |
| zero representative | $`5.93\times10^{-16}`$ | $`0`$ |
| temporal | $`6.71\times10^{-16}`$ | $`5.55\times10^{-17}`$ |
| spatial | $`6.71\times10^{-16}`$ | $`0`$ |
| mixed | $`6.71\times10^{-16}`$ | $`3.99\times10^{-18}`$ |
| temporal–mixed cross | $`5.93\times10^{-16}`$ | $`0`$ |

This check confirms two claims that cannot be seen in a method-of-lines
experiment: the active rank is smaller than $`\dim F_\mathrm{st}`$, and a
genuinely mixed nonseparable gauge can satisfy the same global Stokes
identity as the product representative.

### 10.3. The 2012 ideal line: reconstruction, active stages, and design

We use the physical parameters and complementary boundary causality of the
ideal-line experiment in [[20]](#v12-ref-moulla2012):

```math
\ell=2,\qquad C=3,\qquad L=2,
```

with $`N=8`$.  In the shifted $`L^2(0,\ell)`$-orthonormal Legendre basis,
the exact rectangular power chart is particularly transparent:
$`M_N=[I_N;0]`$, while $`D_N`$ and both endpoint traces are evaluated
exactly.  The residual in [Equation (118)](#v12-eq-rectangular-stokes-line) is
$`3.97\times10^{-13}`$.  Eliminating the degree-$`N`$ radical under the
causality [Equation (122)](#v12-eq-line-port-causality) gives
[Equation (135)](#v12-eq-reduced-open-line), with

```math
\left\lVert \mathcal J_N+\mathcal J_N^\top\right\rVert_F=3.08\times10^{-13},
 \quad
 \left\lVert C_N-B_N^\top\right\rVert_F=7.66\times10^{-14},
 \quad
 \left\lVert N_N+N_N^\top\right\rVert_F=0.
```

The stronger blockwise test of
[Theorem 9.2](#v12-thm-rectangular-stokes-dirac), [Equation (136)](#v12-eq-line-square-block-recovery) is shown
in [Table 4](#v12-tab-rectangular-reconstruction).  Every square open-system block
is reconstructed from $`(M_N,D_N)`$ and the four endpoint traces.

<a name="v12-tab-rectangular-reconstruction"></a>

**Table 4. Rectangular-to-square reconstruction of the open line.**

| Quantity | residual or condition number |
| --- | --- |
| polarized rectangular Stokes block | $`5.61\times10^{-13}`$ |
| recovery of $`\mathcal J_N`$ | $`1.10\times10^{-14}`$ |
| recovery of $`B_N`$ | $`7.35\times10^{-16}`$ |
| recovery of $`B_N^\top`$ from the output trace | $`3.26\times10^{-15}`$ |
| recovery of $`N_N`$ | $`1.11\times10^{-16}`$ |
| minimum singular value of each endpoint reconstruction | $`4.54\times10^{-1}`$ |
| $`\kappa_2(\mathcal R_N)`$ | $`1.42\times10^{1}`$ |

The transmission-line computation is organized as a falsifiable
validation of the rectangular Stokes theorem, the volume-to-cap pullback,
exact trace-generated realization, independence of the raw gauge fixing,
the Cayley bijection and its conditioning, active-stage Gauss
condensation, the open power law, the exact order dichotomy and its unique
order-raising coefficient, the scalar-functional-calculus obstruction,
dispersion-based gauge selection, and cancellation after solving two
coupled cells.  We choose a rank-six alternating $`\Theta`$
coupling distinct $`q`$- and $`p`$-Legendre modes.  Its scale is fixed in
energy-orthonormal coordinates:

<a name="v12-eq-intrinsic-theta-normalization"></a>
**Equation (173).**

```math

 \left\lVert Q_N^{-1/2}\Theta Q_N^{-1/2}\right\rVert_F
 =\left\lVert Q_N^{-1/2}
 (Q_N\mathcal J_NQ_N/12)Q_N^{-1/2}\right\rVert_F.
```

Equivalently, this is the norm on covariant two-forms induced by $`Q_N`$,
and is invariant under a change of state coordinates.
For the reduction test we use the following numerical section:

```math
\sigma_x=\mathcal D_x^\top
 (\mathcal D_x\mathcal D_x^\top)^{-1}
 .
```

It is inserted in [Equation (127)](#v12-eq-line-geometric-lifting) and is not used in the
time-step solve and therefore tests [Proposition 9.9](#v12-prop-mixed-cylinder-lift)
independently; in coordinates the lifted form is
$`\Omega_\mathrm{st}=\mathcal T_N^\top\Theta\mathcal T_N`$.
The cap dimension is $`16`$, so surjectivity predicts all
$`\binom{16}{2}=120`$ cap gauges to be reachable.  The residuals are
reported in [Table 5](#v12-tab-line-pullback).

<a name="v12-tab-line-pullback"></a>

**Table 5. Volume-to-cap pullback on the degree-eight transmission line.**

| Quantity | residual or rank |
| --- | --- |
| $`\Vert \mathcal T_N\mathcal L_N-I\Vert _F`$ | $`1.29\times10^{-14}`$ |
| $`\Vert \Omega_\mathrm{st}+\Omega_\mathrm{st}^\top\Vert _F`$ | $`4.77\times10^{-14}`$ |
| $`\Vert \mathcal L_N^\top\Omega_\mathrm{st}\mathcal L_N-\Theta\Vert _F/ \Vert \Theta\Vert _F`$ | $`5.73\times10^{-15}`$ |
| $`\operatorname{rank}\mathcal T_N`$ | $`16`$ |
| reachable cap-gauge dimension | $`120`$ |

The preceding pullback is followed by an exact raw-realization test.  At
$`h=0.1`$, set $`S_h=S_{2,h}^{\Theta}`$ and assemble the full raw system
[Equation (59)](#v12-eq-raw-cell-block-system) with
$`\Omega_h=\mathcal T_N^\top S_h\mathcal T_N`$.  The raw potential has
dimension $`36`$; the effort remains in the full rectangular space
$`(\mathcal P_8)^2`$ of dimension $`18`$; and the two endpoint conditions
$`e_q(0)=e_p(\ell)=0`$ are imposed explicitly.  Since
$`\ker\mathcal T_N`$ has dimension $`20`$, including the final cap and
the null-space multiplier gives a $`90\times90`$ system.  Its blocks use
$`D_N`$, $`M_N`$, the endpoint traces, and $`\Omega_h`$; neither
$`\mathcal J_N`$ nor a cap propagator is inserted.  After solving the
system for the sixteen coordinate initial states, its $`x_1`$-block is
compared with [Equation (171)](#v12-eq-mixed-gauss-map).

In coordinates, the input trace in
[Equation (134)](#v12-eq-line-effort-input-reconstruction) is

```math
\mathcal C_N=
 \begin{bmatrix}r_0^\top&0\\0&r_\ell^\top\end{bmatrix}.
```

If the columns of $`N_T`$ span $`\ker\mathcal T_N`$, the implemented
rectangular cell is

<a name="v12-eq-raw-rectangular-line-system"></a>
**Equation (174).**

```math

\begin{bmatrix}
 \mathcal T_N&0&-I&0\\
 0&-h\mathcal K_N&I&0\\
 \Omega_h&\mathcal T_N^\top\mathcal M_N^*&
   -\tfrac12\mathcal T_N^\top Q_N&N_T\\
 N_T^\top&0&0&0\\
 0&\mathcal C_N&0&0
\end{bmatrix}
\begin{bmatrix}z\\e\\x_1\\\lambda\end{bmatrix}
=
\begin{bmatrix}
 -I\\I\\\tfrac12\mathcal T_N^\top Q_N\\0\\0
\end{bmatrix}x_0.
```

The square map $`\mathcal J_N`$ is itself recovered only during
elimination:

<a name="v12-eq-rectangular-generator-recovery"></a>
**Equation (175).**

```math

 \mathcal J_N=
 \mathcal K_N
 \begin{bmatrix}\mathcal M_N^*\\\mathcal C_N\end{bmatrix}^{-1}
 \begin{bmatrix}I\\0\end{bmatrix}.
```

<a name="v12-tab-dynamic-condensation"></a>

**Table 6. Exact realization by the raw trace-generated line cell.**

| Quantity | value |
| --- | --- |
| raw block-system residual | $`6.80\times10^{-15}`$ |
| relative error of the condensed cap map | $`1.92\times10^{-16}`$ |
| relative change under a different right inverse | $`1.05\times10^{-14}`$ |
| cap-map change under a different raw gauge fixing | $`1.17\times10^{-16}`$ |
| alternative raw-system residual | $`6.73\times10^{-15}`$ |
| null-space multiplier norm | $`7.75\times10^{-16}`$ |
| recovery of $`\mathcal J_N`$ from rectangular blocks | $`1.10\times10^{-14}`$ |
| condition number of the raw system | $`1.53\times10^{2}`$ |

This test has the precise scope of
[Theorem 4.14](#v12-thm-trace-generated-condensation): the prescribed cap gauge is
encoded in the raw volume matrix and recovered by elimination.  It is not
presented as selection of that gauge by independent interior dynamics.
The additional gauge-fixing test solves the complete raw system with a
nonorthogonal complement of $`\ker\mathcal T_N`$; the physical cap map is
unchanged to roundoff.

The genuinely interior test is instead the polynomial Gauss cell of
[Theorem 8.8](#v12-thm-active-gauss-cell).  We assemble its $`sn`$-dimensional stage
block for $`s=1,\ldots,4`$, eliminate the stages, and only then recover the
cap gauge by [Equation (71)](#v12-eq-inverse-gauge-cayley).  The results are in
[Table 7](#v12-tab-active-stage-condensation); the largest energy defect is
$`2.69\times10^{-15}`$.

<a name="v12-tab-active-stage-condensation"></a>

**Table 7. Active temporal-stage elimination at $`h=0.1`$, $`n=16`$.**

| $`s`$ | $`sn`$ | stage residual | map–Padé | induced gauge | $`\kappa_2`$ |
| --- | --- | --- | --- | --- | --- |
| 1 | 16 | $`4.73\times10^{-16}`$ | $`1.77\times10^{-16}`$ | $`1.88\times10^{-14}`$ | 1.608 |
| 2 | 32 | $`9.87\times10^{-16}`$ | $`2.59\times10^{-16}`$ | $`2.90\times10^{-14}`$ | 2.976 |
| 3 | 48 | $`1.49\times10^{-15}`$ | $`2.45\times10^{-16}`$ | $`2.31\times10^{-14}`$ | 3.660 |
| 4 | 64 | $`2.00\times10^{-15}`$ | $`1.35\times10^{-15}`$ | $`4.75\times10^{-14}`$ | 3.805 |

The Cayley coordinate is first checked independently at the moderate step
$`h=0.1`$.  A prescribed alternating gauge is mapped to $`O(Q_N)`$,
inverted by [Equation (71)](#v12-eq-inverse-gauge-cayley), and mapped forward again.  The
gauge and map round trips have relative errors $`6.07\times10^{-13}`$
and $`1.93\times10^{-16}`$; the energy-orthogonality defect is
$`2.34\times10^{-15}`$.  Applying the same algebraic inverse to
$`\exp(hA_N)`$ and reconstructing the flow gives
$`2.99\times10^{-16}`$ relative error.

That test does not conceal the singular conditioning of the inverse chart
near $`I`$.  In [Table 8](#v12-tab-cayley-conditioning), the reference exact-flow
gauge is evaluated directly by the Bernoulli series in
[Equation (74)](#v12-eq-exact-flow-gauge-series).  The ``naive” column instead forms
$`\exp(hA_N)-I`$ and applies [Equation (71)](#v12-eq-inverse-gauge-cayley).  Its rapidly
growing error is predicted by [Remark 5.5](#v12-rem-cayley-conditioning), while the
map reconstructed from the direct gauge remains at roundoff.

<a name="v12-tab-cayley-conditioning"></a>

**Table 8. Stable exact-flow gauge evaluation near the excluded identity.**

| $`h`$ | relative naive gauge error | error$`\,h^3/\epsilon_\mathrm{mach}`$ | stable map error |
| --- | --- | --- | --- |
| $`10^{-2}`$ | $`2.36\times10^{-10}`$ | 1.062 | $`1.43\times10^{-16}`$ |
| $`10^{-4}`$ | $`2.60\times10^{-4}`$ | 1.169 | $`1.47\times10^{-16}`$ |
| $`10^{-6}`$ | $`2.62\times10^{2}`$ | 1.181 | $`1.55\times10^{-16}`$ |
| $`10^{-8}`$ | $`1.31\times10^{8}`$ | 0.591 | $`1.11\times10^{-16}`$ |

The same form is lifted to the cylinder by
[Equation (125)](#v12-eq-mixed-cylinder-lift).  At $`h=0.04`$, the open
balance is tested for midpoint, the second-order mixed deformation
$`S_h=h\Theta`$, and the fourth-order deformation
[Equation (170)](#v12-eq-mixed-gauss-two-gauge).

<a name="v12-tab-rectangular-line-power"></a>

**Table 9. Direct rectangular transmission-line cap–mantle balance.**

| Gauge | $`\Delta H_N`$ | $`h\,y^\top u`$ | residual |
| --- | --- | --- | --- |
| midpoint | $`2.67994318\times10^{-2}`$ | $`2.67994318\times10^{-2}`$ | $`5.13\times10^{-16}`$ |
| mixed, order two | $`2.67401314\times10^{-2}`$ | $`2.67401314\times10^{-2}`$ | $`6.31\times10^{-16}`$ |
| Gauss–2 plus mixed, order four | $`2.74765368\times10^{-2}`$ | $`2.74765368\times10^{-2}`$ | $`2.71\times10^{-16}`$ |

For the homogeneous line, the zero-mixed member of
[Equation (170)](#v12-eq-mixed-gauss-two-gauge) is Gauss–2.  At $`h=0.05`$, its relative
commutator with $`A_N`$ is $`6.92\times10^{-17}`$.  The mixed deformation
has energy defect $`1.51\times10^{-15}`$, relative commutator
$`2.59\times10^{-6}`$, and relative symplectic defect
$`1.20\times10^{-6}`$.  Refinement gives [Table 10](#v12-tab-mixed-gauss-order);
the one-step error has the expected fifth-power scaling, while the
commutator remains decisively above roundoff.

<a name="v12-tab-mixed-gauss-order"></a>

**Table 10. Fourth-order mixed deformation of Gauss–2 on the closed line.**

| $`h`$ | relative one-step error | one-step error slope | relative $`\left\lVert [M_{\Theta,h},A_N]\right\rVert_F`$ |
| --- | --- | --- | --- |
| 0.1000 | $`8.9501\times10^{-3}`$ | --- |$`7.3650\times10^{-5}`$ |
| 0.0500 | $`3.3624\times10^{-4}`$ | 4.734 | $`2.5947\times10^{-6}`$ |
| 0.0250 | $`1.1167\times10^{-5}`$ | 4.912 | $`8.3450\times10^{-8}`$ |
| 0.0125 | $`3.5483\times10^{-7}`$ | 4.976 | $`2.6264\times10^{-9}`$ |

To distinguish a one-step slope from convergence of the integrator, we
also fix the final time $`T=0.5`$, compose the map, and compare it with
$`\exp(TA_N)`$.  For each $`s`$, the mixed coefficient $`\Theta`$ is held
fixed and the perturbation is exactly $`h^{2s-1}\Theta`$.  The final
refinement slopes in [Table 11](#v12-tab-gauss-hierarchy-global) approach the proved
global orders $`2,4,6`$; the largest one-step energy defect over the table
is $`3.08\times10^{-15}`$.

<a name="v12-tab-gauss-hierarchy-global"></a>

**Table 11. Global convergence of the stage-induced, mixed Gauss–$`s`$ hierarchy
at $`T=0.5`$.  The four error columns correspond to $`4,8,16,32`$ steps.**

| $`s`$ | expected order | $`n=4`$ | $`n=8`$ | $`n=16`$ | $`n=32`$ | final slope |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 2 | $`6.398\!\times10^{-1}`$ | $`2.870\!\times10^{-1}`$ | $`8.393\!\times10^{-2}`$ | $`2.179\!\times10^{-2}`$ | 1.946 |
| 2 | 4 | $`9.816\!\times10^{-2}`$ | $`8.058\!\times10^{-3}`$ | $`5.368\!\times10^{-4}`$ | $`3.407\!\times10^{-5}`$ | 3.978 |
| 3 | 6 | $`4.462\!\times10^{-3}`$ | $`8.240\!\times10^{-5}`$ | $`1.340\!\times10^{-6}`$ | $`2.115\!\times10^{-8}`$ | 5.986 |

The preceding coefficient is noncommuting and therefore cannot equal the
exceptional value in [Equation (146)](#v12-eq-unique-order-raising-gauge).  The same
experiment is now repeated with $`\Theta=\Theta_s^\star`$.  The observed
orders in [Table 12](#v12-tab-order-raising-gauge) are $`4,6,8`$, and the largest
energy defect is $`3.27\times10^{-15}`$.  For $`s=1`$, the relative
difference between the raised map and Gauss–2 at $`h=0.1`$ is
$`2.42\times10^{-16}`$, confirming the exact identity in
[Theorem 9.14](#v12-thm-unique-order-raising-gauge).

<a name="v12-tab-order-raising-gauge"></a>

**Table 12. Global convergence of the unique order-raising coefficient in
the fixed one-term ansatz at $`T=0.5`$.  The error columns correspond to
$`n`$ time steps.**

| $`s`$ | predicted order | $`n=4`$ | $`n=8`$ | $`n=16`$ | $`n=32`$ | $`n=64`$ | final slope |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 4 | $`9.950\!\times10^{-2}`$ | $`8.184\!\times10^{-3}`$ | $`5.454\!\times10^{-4}`$ | $`3.462\!\times10^{-5}`$ | $`2.172\!\times10^{-6}`$ | 3.994 |
| 2 | 6 | $`1.335\!\times10^{-2}`$ | $`2.689\!\times10^{-4}`$ | $`4.466\!\times10^{-6}`$ | $`7.084\!\times10^{-8}`$ | $`1.111\!\times10^{-9}`$ | 5.995 |
| 3 | 8 | $`1.027\!\times10^{-3}`$ | $`4.985\!\times10^{-6}`$ | $`2.053\!\times10^{-8}`$ | $`8.123\!\times10^{-11}`$ | $`3.189\!\times10^{-13}`$ | 7.993 |

The order statement is not inferred from these slopes.  To test the new
coefficient in [Equation (147)](#v12-eq-raised-branch-leading-error) directly, we use the
unit harmonic oscillator $`A=J_2`$ and 90-digit arithmetic.  After
dividing the one-step matrix error by $`h^{2s+3}`$, its component along
$`J_2`$ must converge to $`(-1)^{s+1}d_s`$.  At $`h=0.025`$ the recovered
values in [Table 13](#v12-tab-direct-ds-coefficient) already agree to about
$`3.5\times10^{-4}`$ relatively.  The full normalized remainder decreases
linearly with $`h`$, exactly as the $`O(h^{2s+4})`$ term predicts.

<a name="v12-tab-direct-ds-coefficient"></a>

**Table 13. Direct local asymptotic test of the exact coefficient $`d_s`$.**

| $`s`$ | $`d_s`$ | predicted $`J_2`$ coefficient | recovered coefficient at $`h=0.025`$ |
| --- | --- | --- | --- |
| 1 | $`-1/720`$ | $`-1.388888889\times10^{-3}`$ | $`-1.388403230\times10^{-3}`$ |
| 2 | $`1/30240`$ | $`-3.306878307\times10^{-5}`$ | $`-3.305724437\times10^{-5}`$ |
| 3 | $`-1/2268000`$ | $`-4.409171076\times10^{-7}`$ | $`-4.407653975\times10^{-7}`$ |

We next retain $`\Theta_s^\star`$ and add the same fixed noncommuting line
two-form $`h^{2s+1}\Xi`$ for $`s=1,2,3`$.  The global orders remain
$`4,6,8`$, but the commutators in
[Table 14](#v12-tab-two-term-noncommuting) are nonzero at resolved step sizes.  This
is the experimental distinction absent from the one-term exceptional
branch: high order and departure from scalar functional calculus now
coexist.  The largest one-step energy defect in the experiment is
$`3.11\times10^{-15}`$.

<a name="v12-tab-two-term-noncommuting"></a>

**Table 14. Two-term high-order noncommuting hierarchy on the transmission
line.  Errors are global at $`T=0.5`$; commutators are reported at
$`n=16`$, before the higher-stage values reach roundoff.**

| $`s`$ | exact order | error, $`n=32`$ | final slope | relative commutator, $`n=16`$ |
| --- | --- | --- | --- | --- |
| 1 | 4 | $`3.407\times10^{-5}`$ | 3.994 | $`2.533\times10^{-7}`$ |
| 2 | 6 | $`7.071\times10^{-8}`$ | 5.995 | $`2.474\times10^{-10}`$ |
| 3 | 8 | $`8.120\times10^{-11}`$ | 7.993 | $`2.416\times10^{-13}`$ |

For $`s=3`$, binary64 eventually masks the noncommuting signal.  We
therefore repeat the one-step calculation in 90-digit decimal arithmetic
on the following fully specified four-dimensional Hamiltonian test:

<a name="v12-eq-decimal-test-matrices"></a>
**Equation (176).**

```math

 \mathcal J=
 \begin{bmatrix}
 0&1&0&0\\-1&0&0&0\\0&0&0&1\\0&0&-1&0
 \end{bmatrix},\qquad
 Q=\operatorname{diag}(1,2,3,5),\qquad
 A=\mathcal JQ,
 \qquad
 \Xi=
 \begin{bmatrix}
 0&0&1&0\\0&0&0&0\\-1&0&0&0\\0&0&0&0
 \end{bmatrix}.
```

The common block is evaluated directly as

```math
C_3(hA)+\frac{(hA)^6}{100800}+h^8\mathcal J\Xi,
 \qquad
 C_3(z)=\frac12\left(I+\frac{z^2}{10}\right)
                 \left(\frac I2+\frac{z^2}{120}\right)^{-1},
```

using dense decimal Gaussian elimination with partial pivoting.  Thus the
experiment requires no hidden random matrix or fitted coefficient.
The theorem predicts

```math
[M_h,A]=-h^9[\mathcal J\Xi,A]A+O(h^{10}).
```

[Table 15](#v12-tab-decimal-noncommutation) verifies both the exponent and the full
leading matrix after division by $`h^9`$; this is independent of the
global-order slope above.

<a name="v12-tab-decimal-noncommutation"></a>

**Table 15. Multiprecision resolution of the $`s=3`$ noncommutation law.**

| $`h`$ | $`\Vert [M_h,A]\Vert _F`$ | observed exponent | relative error in the leading matrix |
| --- | --- | --- | --- |
| 0.050000 | $`4.044\times10^{-11}`$ | -- |$`8.436\times10^{-2}`$ |
| 0.025000 | $`7.901\times10^{-14}`$ | 8.99967 | $`4.223\times10^{-2}`$ |
| 0.012500 | $`1.543\times10^{-16}`$ | 8.99992 | $`2.112\times10^{-2}`$ |
| 0.006250 | $`3.014\times10^{-19}`$ | 8.99998 | $`1.056\times10^{-2}`$ |

The reduced-polynomial checks are independent of the small-step order
tests.  At the $`C_3`$ chart resonance the computed
$`|\widehat O_3|`$ is $`1.11\times10^{-16}`$, while
$`|\Pi_{3,1}+5|=4.44\times10^{-16}`$ and the continued map differs from
one by $`1.78\times10^{-16}`$.  Root computation for
[Equation (164)](#v12-eq-pi-21) gives the genuine left-half-plane pole
$`-10.12751490652`$, with denominator residual
$`1.78\times10^{-14}`$ and numerator magnitude $`10.1275`$.  On 4001
points $`\xi\in[-100,100]`$, the largest defect in
$`|R_{2,1}(i\xi)|=1`$ is $`4.44\times10^{-16}`$.  These three tests
separate a removable coordinate resonance, a true propagator pole, and
the protected oscillatory axis.

Finally, exact `Fraction` arithmetic evaluates all 156 coefficients
$`\alpha_{s,j}`$ with $`1\le s\le12`$, $`0\le j\le12`$.  No zero or sign
failure relative to [Equation (169)](#v12-eq-alpha-sign-conjecture) occurs.  This is
reported only as finite evidence for the conjecture, not as its proof.

The interior-generation assertions directly target the two new results.
For a rank-four random alternating cap change in dimension six, the
four-dimensional construction in
[Theorem 9.10](#v12-thm-minimal-schur-gauge-generation) has zero skew defect, Schur
error $`8.02\times10^{-16}`$, and eliminated-response error
$`1.83\times10^{-15}`$.  For Gauss–3, the one-mode data

```math
a_3=0,\qquad w_{3,1}=\frac1{12},\qquad
 \lambda_{3,1}=\frac1{60}
```

reproduce both $`F_3(u)=5/(u+60)`$ and $`C_3`$ to roundoff on
$`u\in[0,120]`$.  At $`u=61`$, the corrected value in
[Equation (168)](#v12-eq-corrected-stieltjes-obstruction) is exactly $`-1/720`$, so the
positive-realization obstruction is not a floating-point inference.

Choosing the second central corrector
$`\Xi_s^\star=d_s\Omega_NA_N^{2s+2}`$ gives observed global orders six
and eight for $`s=1,2`$, respectively; see
[Table 16](#v12-tab-second-central-corrector).  These are the first two nontrivial
levels of the recursive exact-gauge expansion in
[Corollary 9.17](#v12-cor-recursive-gauge-hierarchy).

<a name="v12-tab-second-central-corrector"></a>

**Table 16. Second central correction on the transmission line.**

| $`s`$ | proved lower bound | error, $`n=32`$ | final slope |
| --- | --- | --- | --- |
| 1 | 6 | $`7.084\times10^{-8}`$ | 5.995 |
| 2 | 8 | $`1.522\times10^{-10}`$ | 7.994 |

The affine freedom can also be used as a design space.  We minimize the
relative root-mean-square error of the eight positive modal frequencies
within the one-dimensional family
$`S_{2,h}^{\alpha\Theta}=S^\mathrm{G}_{2,h}+\alpha h^3\Theta`$, using only
the small reference step $`h_\star=1/80`$.  This gives
$`\alpha_\star=57.3601`$, which is then frozen.  The tests in
[Table 17](#v12-tab-dispersion-selector) use only unseen, coarser steps.

<a name="v12-tab-dispersion-selector"></a>

**Table 17. Dispersion selection in a fixed one-parameter mixed family.**

| $`h`$ | Gauss RMS phase error | selected error | reduction | $`\kappa_2`$ Gauss | $`\kappa_2`$ selected |
| --- | --- | --- | --- | --- | --- |
| 0.05 | $`3.605\times10^{-4}`$ | $`1.143\times10^{-4}`$ | 68.28% | 1.202 | 1.202 |
| 0.10 | $`4.864\times10^{-3}`$ | $`1.875\times10^{-3}`$ | 61.45% | 1.448 | 1.443 |
| 0.15 | $`1.845\times10^{-2}`$ | $`9.589\times10^{-3}`$ | 48.03% | 1.861 | 1.941 |

The selected map remains exactly energy preserving; only its phase is
optimized.  The improvement persists away from the fitting step and the
linear-system condition number changes mildly.  Repeating the norm
calculation after a fixed nonorthogonal coordinate change gives a relative
difference below displayed precision, confirming the invariance of
[Equation (173)](#v12-eq-intrinsic-theta-normalization).  This experiment does not claim
universal superiority of the selected gauge, but it does answer a sharper
question: the mixed degree of freedom can correct dispersion without
sacrificing the exact power identity or the proved temporal order.

Finally, two degree-eight cells are solved simultaneously with the
power-preserving interface constraints of
[Proposition 9.8](#v12-prop-rectangular-line-power), the fourth-order mixed gauge, and
$`R_N=0.02I`$.  The algebraic solve residual is
$`5.06\times10^{-16}`$, the computed interface power is
$`3.33\times10^{-16}`$, and

```math
\Delta H_{N,1}+\Delta H_{N,2}
 =2.75921723\times10^{-2},
```

while external supply minus dissipation is

```math
3.29176951\times10^{-2}-5.32552281\times10^{-3}
 =2.75921723\times10^{-2},
```

with residual $`1.78\times10^{-15}`$.  These calculations test exactly
the structural chain asserted in this paper: the rectangular Stokes
identity produces the open Dirac chart; active temporal stages generate
the Gauss cap gauge; the mixed volume form realizes its prescribed
perturbation; contracted power is unchanged; the hierarchy attains orders
exactly $`2s`$ away from the unique one-term coefficient
$`\Theta_s^\star`$, which raises the order to exactly $`2s+2`$; a second
noncommuting term retains that higher order and leaves the scalar
functional-calculus subclass; the central series raises the order further;
a fixed gauge can reduce
modal phase error; and internal supplies cancel in a solved, dissipative
two-cell system.  Additional assertions verify a rank-minimal skew Schur
generator, the positive one-mode realization of $`C_3`$, and the sign
obstruction $`(C_{2,1}(\sqrt u)-1)/u<0`$ beyond $`u=60`$.

All values are generated by the accompanying script
`verify_transmission_line_v12.py`.  Its 28 assertion groups use
explicit tolerances and convergence thresholds; the program terminates
with an `AssertionError` if any reported structural claim fails.

### 10.4. Numerical verification of the linear locus

We verify the identities in double precision on the two-mode system of
[Example 6.6](#v12-ex-two-modes), with

```math
Q=\operatorname{diag}(I_2,2I_2),\qquad
 S_{\mathrm{in}}=\operatorname{diag}(0.3J_2,-0.2J_2),\qquad
 S_{\mathrm{out}}=
 \begin{bmatrix}0&I_2\\-I_2&0\end{bmatrix}.
```

Here $`S_{\mathrm{in}}\in\mathcal{L}`$ and
$`S_{\mathrm{out}}\notin\mathcal{L}`$. At $`h=0.07`$, the Frobenius residuals

```math
\varepsilon_H=\left\lVert M_S^\top QM_S-Q\right\rVert_F,\qquad
 \varepsilon_\Omega=\left\lVert M_S^\top\Omega M_S-\Omega\right\rVert_F
```

are shown in [Table 18](#v12-tab-residuals). Energy is at roundoff for every gauge,
whereas the symplectic residual separates the locus.

<a name="v12-tab-residuals"></a>

**Table 18. One-step structural residuals at $`h=0.07`$.**

| Gauge | $`\varepsilon_H`$ | $`\varepsilon_\Omega`$ |
| --- | --- | --- |
| $`S=0`$ | $`5.14\times10^{-16}`$ | $`2.94\times10^{-16}`$ |
| $`S_{\mathrm{in}}`$ | $`2.90\times10^{-16}`$ | $`1.58\times10^{-16}`$ |
| $`S_{\mathrm{out}}`$ | $`7.19\times10^{-16}`$ | $`9.82\times10^{-3}`$ |

For the off-locus gauge, decreasing $`h`$ from $`0.1`$ to $`0.0125`$
gives

```math
\frac{\varepsilon_\Omega}{h^2}
 =2.0074,\ 2.0019,\ 2.0005,\ 2.0001,
```

confirming [Theorem 6.3](#v12-thm-defect). To test order, we integrate to $`T=1`$ from
$`x_0=(0.7,-0.2,0.4,0.9)^\top`$ and compare with
$`\exp(A)x_0`$. The last two mesh refinements in
[Table 19](#v12-tab-order) show the predicted orders.

<a name="v12-tab-order"></a>

**Table 19. Global error and observed order at $`T=1`$.**

| Gauge | error, $`h=1/80`$ | error, $`h=1/160`$ | observed order |
| --- | --- | --- | --- |
| $`S=0`$ | $`1.0302\times10^{-4}`$ | $`2.5757\times10^{-5}`$ | $`1.9999`$ |
| $`S_{\mathrm{in}}`$ | $`5.7094\times10^{-3}`$ | $`2.8350\times10^{-3}`$ | $`1.0100`$ |
| $`S_{\mathrm{out}}`$ | $`2.5020\times10^{-2}`$ | $`1.2547\times10^{-2}`$ | $`0.9957`$ |

These checks are algebraic rather than qualitative long-time plots: any
contrary result falsifies either the implementation or one of the stated
nondegeneracy assumptions.

### 10.5. Summary of established results and hypotheses

The conjectural statements motivating this work are resolved as follows.

1. For an arbitrary linear map $`D:E\to F`$, the unrestricted active
     gauge is exactly $`\Lambda^2(E/\ker D)^*`$:
     [Theorem 2.2](#v12-thm-affine-stokes).

1. For the genuine total complex, symmetric doubling gives
     $`\Lambda^2(U_q\oplus V_q)^*`$; imposing complementary degree reduces
     this to $`U_q^*\otimes V_q^*`$:
     [Theorem 4.3](#v12-thm-graded-affine).  The tensor pairing and both oriented boundary
     terms, including their Koszul signs, are explicit in
     [Equation (16)](#v12-eq-signed-tensor-pairing), [Equation (17)](#v12-eq-signed-product-boundary).  This includes
     cancellations between bidegrees through the actual total-complex kernels.

1. For the distinct stacked tensor-cylinder map,
     $`\ker D_\mathrm{st}=K_t\otimes K_x`$, the rank is
     $`n_tn_x-k_tk_x`$, and the gauge decomposes into temporal, spatial, and
     mixed sectors through a canonical filtration:
     [Theorem 4.6](#v12-thm-space-time-gauge), [Proposition 4.7](#v12-prop-canonical-sector-sequence), [Corollary 4.8](#v12-cor-gauge-sectors).

1. In every fixed regular reduction chart, the volume-to-cap map is
     $`L_h^*`$, and its image, kernel, and dimensions are exactly
     [Equation (40)](#v12-eq-volume-cap-image)–[Equation (44)](#v12-eq-volume-cap-kernel-dimension).  The
     canonical mixed part reduces through $`(\pi_{tx}L_h)^*`$, with no chosen
     splitting: [Corollary 4.12](#v12-cor-mixed-visibility).

1. For a surjective trace, the gauge-bearing weak cell
     [Equation (59)](#v12-eq-raw-cell-block-system) realizes exactly the prescribed cap
     gauge in [Equation (57)](#v12-eq-trace-generated-cap-map), provided the stated map
     block is invertible.  This is an encoding-and-recovery result, not a
     claim of interior selection.  The fixed-reconstruction theorem is also
     the tangent map of a gauge-dependent reduction:
     [Proposition 4.16](#v12-prop-gauge-dependent-tangent).

1. The complete temporal linear gauge is
     $`\Lambda^2V^*`$; which members arise from a selected volume model is
     determined by $`\operatorname{im}\mathcal R_{L_h}`$:
     [Theorem 4.10](#v12-thm-volume-cap-reduction), [Theorem 5.1](#v12-thm-two-point-gauge).

1. For nonsingular $`J,Q`$, the temporal family gives a bijective,
     constitutively shifted Cayley coordinate on
     $`\{M\in O(Q):M-I\text{ invertible}\}`$.  Its inverse and the exact-flow
     gauge are [Equation (71)](#v12-eq-inverse-gauge-cayley), [Equation (73)](#v12-eq-exact-flow-gauge); the
     inverse coordinate is condition-sensitive near $`I`$, as stated in
     [Remark 5.5](#v12-rem-cayley-conditioning), and the
     symplectic subfamily maps to $`O(Q)\cap Sp(\Omega)`$.

1. The gauge-invariant power balance is exact, including ports and
     dissipation: [Corollary 5.2](#v12-cor-power), [Corollary 9.4](#v12-cor-composition).

1. The midpoint rule is uniquely selected by endpoint symmetry,
     midpoint localization, and—for nonsingular $`J,Q`$—order two:
     [Proposition 5.6](#v12-prop-selectors).

1. Symplecticity does not select the zero gauge. Its exact locus is
     [Equation (76)](#v12-eq-locus); generically it has dimension $`m`$, and in two
     dimensions it is the entire gauge:
     [Theorem 6.4](#v12-thm-generic-dimension), [Corollary 6.5](#v12-cor-two-dimensional).

1. AVF is uniquely selected only in the stated universal
     segment-measure class: [Theorem 7.2](#v12-thm-avf-unique).

1. For invertible stage matrices, the stage-pairing existence
     criterion, dimension, and relation to $`R(\infty)`$ are given by
     [Theorem 8.2](#v12-thm-p-family); the singular case remains the range problem stated
     in [Remark 8.3](#v12-rem-singular-stage).

1. Among collocation methods, local compatible pairing is equivalent
     to Gauss–Legendre collocation: [Theorem 8.5](#v12-thm-gauss-rigidity).

1. Conversely, the active temporal cell
     $`\mathcal P_s\to\mathcal P_{s-1}`$ generates the Gauss–$`s`$ cap map
     by eliminating $`sn`$ interior stage variables; the even–odd Padé
     factorization gives its two-point gauge directly:
     [Proposition 8.7](#v12-prop-local-stokes-to-stages), [Theorem 8.8](#v12-thm-active-gauss-cell).

1. The global tensor Stokes equations force a local diagonal pairing
     to factor, and maximal polynomial exactness then gives the product Gauss
     rule; this is not a uniqueness statement among nonlocal gauges:
     [Theorem 9.5](#v12-thm-tensor-rigidity).

1. A genuine finite-dimensional Dirac claim requires isotropy and
     maximality. These are proved in [Proposition 9.1](#v12-prop-dirac-maximal), including an
     antisymmetric feedthrough; the temporal effort is correctly identified
     as a constitutive closure.

1. An invertible rectangular effort-and-input reconstruction turns
     the polarized Stokes identity into a square open Dirac representation,
     forcing skew interconnection, input/output duality, and skew
     feedthrough: [Theorem 9.2](#v12-thm-rectangular-stokes-dirac).  For the polynomial
     line, the required reconstruction is proved invertible in
     [Lemma 9.12](#v12-lem-polynomial-endpoint-reconstruction).

1. A mixed gauge can preserve the full cylinder identity while having
     zero pure temporal and spatial pullbacks; the polynomial cell check in
     [Table 3](#v12-tab-space-time-residuals) exhibits this possibility directly.

1. On the rectangular 2012 transmission-line complex, the cap trace
     is onto and every cap two-form has the explicit pure-mixed lift
     [Equation (125)](#v12-eq-mixed-cylinder-lift).  The hierarchy
     [Equation (139)](#v12-eq-mixed-gauss-gauge) retains exact energy and reversibility.  It
     has exactly order $`2s`$ for every
     $`\Theta\ne\Theta_s^\star`$, whereas the unique fixed coefficient in
     the one-term ansatz [Equation (146)](#v12-eq-unique-order-raising-gauge) raises the order
     to exactly $`2s+2`$.  The two-term family
     [Equation (154)](#v12-eq-two-term-gauge) retains that order with an arbitrary next
     coefficient, and the central choice [Equation (156)](#v12-eq-second-central-corrector)
     raises it by two more.  Moreover,
     [Equation (143)](#v12-eq-mixed-map-commutator) detects departure
     from the scalar functional calculus of the semidiscrete line.

1. The central truncation has the intrinsic reduced form
     $`R_{s,r}=\Pi_{s,r}(-z)^{-1}\Pi_{s,r}(z)`$:
     [Theorem 9.18](#v12-thm-central-rational-hierarchy).  Coordinate poles of $`C_s`$
     are removable $`+1`$ resonances, the reduced map has no pole on the
     imaginary axis, but it need not be $`A`$-stable; the explicit
     $`(s,r)=(2,1)`$ member has a left-half-plane pole:
     [Proposition 9.20](#v12-prop-rational-resonance-stability).  The exact coordinate is a
     Stieltjes function in $`z^2`$, and $`C_s`$ is its stated Padé
     approximant: [Proposition 9.23](#v12-prop-stieltjes-pade-coordinate).  The strict sign
     law for every $`\alpha_{s,j}`$ remains the explicitly labelled
     conjecture [Remark 9.26](#v12-rem-alpha-sign-conjecture).

1. At the unrestricted algebraic level, every rank-$`2q`$ change of
     cap gauge is generated by eliminating exactly $`2q`$ auxiliary
     alternating modes, and no smaller generator of that block form is
     possible: [Theorem 9.10](#v12-thm-minimal-schur-gauge-generation).  This sharp
     existence-and-dimension statement is distinct from a field-selection
     mechanism.  The stronger positive
     Jacobi-chain mechanism generates every Gauss coordinate, but not the
     corrected $`(2,1)`$ member:
     [Theorem 9.24](#v12-thm-stieltjes-interior-generation).  Thus algebraic generation,
     passive generation, and selection by an independent field cell are
     three distinct statements.
     [Table 4](#v12-tab-rectangular-reconstruction), [Table 6](#v12-tab-dynamic-condensation), [Table 7](#v12-tab-active-stage-condensation), [Table 8](#v12-tab-cayley-conditioning), [Table 11](#v12-tab-gauss-hierarchy-global), [Table 12](#v12-tab-order-raising-gauge), [Table 13](#v12-tab-direct-ds-coefficient), [Table 14](#v12-tab-two-term-noncommuting), [Table 15](#v12-tab-decimal-noncommutation), [Table 16](#v12-tab-second-central-corrector), [Table 17](#v12-tab-dispersion-selector)
     verify the three reductions, stable coordinate evaluation, the direct
     asymptotic coefficient, the one- and two-term order branches, the
     multiprecision noncommutation exponent, and a
     concrete dispersion gain; the solved two-cell test verifies interface
     power.

<a name="v12-section-11"></a>

## 11. Discussion and limitations

Four layers should remain distinct in future extensions.

1. *Stokes compatibility* is an affine constraint on a pairing.
     Its homogeneous direction is alternating on a quotient space.

1. *Dirac maximality* concerns a subspace of flows and efforts
     equipped with a symmetric split pairing. Isotropy alone is
     insufficient; the dimension or orthogonal-complement condition must also
     be checked.

1. *Constitutive closure* inserts the Hamiltonian through a
     continuous or discrete gradient. It supplies energy balance but is not
     itself a proof of maximal isotropy.

1. *Boundary reduction* composes the preceding layers and
     eliminates interior variables.  Only after a regular reconstruction $`L_h`$
     has been established is the reduced gauge $`L_h^*\omega`$; singular
     interior blocks require a reduction of relations.

Compatible gauges should also be distinguished by their origin.

1. A *geometric or premetric gauge* depends only on the
     differential, traces, orientation, and approximation spaces.  The local
     Gauss stage pairing is of this type.

1. A *constitutively adapted gauge* is obtained only after
     $`J,Q`$, and $`h`$ are supplied.  The induced Gauss cap gauge
     [Equation (99)](#v12-eq-stage-induced-gauss-gauge) and the exact-flow gauge
     [Equation (73)](#v12-eq-exact-flow-gauge) belong to this class.

1. A *design gauge* uses the remaining compatible freedom to
     optimize dispersion, conditioning, sparsity, or another numerical
     objective.  The one-parameter mixed experiment is an example, not a
     universal selector.

All three classes satisfy the same affine Stokes constraint, but only the
first can be attributed to the discrete geometry alone.

The central rational hierarchy is tailored to Hamiltonian generators, not
to the whole left half-plane.  Its reduced polynomial form removes the
artificial poles of the $`C_s`$-coordinate and proves unit amplification
on the oscillatory axis, but [Proposition 9.20](#v12-prop-rational-resonance-stability)
shows that high order and exact quadratic conservation do not imply
$`A`$-stability.  Consequently dissipative extensions require an
independent root analysis of $`\Pi_{s,r}`$.  Likewise, the Stieltjes–Padé
identification is proved, whereas the coefficientwise sign and
non-annulation law [Equation (169)](#v12-eq-alpha-sign-conjecture) is not; all hierarchy
orders beyond the computed levels retain the stated conditional wording.

The cylinder and double-complex results are finite-dimensional discrete
statements.  They do not by themselves prove convergence to an
infinite-dimensional field theory.  Such a result requires approximation,
stability, trace, and compactness estimates for the selected trial spaces.
Recent space–time FEEC work combines fourth-order Gauss–Legendre time
integration with a spatial finite-element complex and proves a local
multisymplectic conservation law for nonlinear wave equations
[[7]](#v12-ref-cheng2026).  That PDE convergence and local variational statement is
complementary to, and stronger in a different direction than, the exact
finite-dimensional quadratic cap-power law proved here.
Descriptor systems and state-dependent mass matrices require additional
constraint hypotheses; discrete-gradient pairs provide one relevant
extension [[30]](#v12-ref-schulze2023).  Likewise, the compatible graded complex
constructed here supplies the Stokes and power layers of a discrete field
theory, but full multi-Dirac integrability is a further geometric property
[[34]](#v12-ref-vankerschaver2012).  Finally, the product-Gauss uniqueness theorem is
deliberately local and diagonal.  The mixed sectors in
[Equation (37)](#v12-eq-wedge-sector-split) show explicitly why nonlocal space–time
Petrov–Galerkin pairings remain possible.

The visibility theorem is finite-dimensional and exact, but conditional
on a fixed regular reconstruction.  It classifies every volume two-form once the
Dirac relation, constitutive law, boundary causality, and interior solve
are fixed; it does not assert that these data are gauge independent.
[Proposition 4.16](#v12-prop-gauge-dependent-tangent) identifies its precise first variation,
while [Theorem 4.14](#v12-thm-trace-generated-condensation) proves exact raw realization
for a gauge already encoded as $`T_h^*S_hT_h`$.  The latter is not a
derivation of $`S_h`$ from independent volume dynamics.  The missing
interior mechanism is supplied in one direction by
[Theorem 8.8](#v12-thm-active-gauss-cell): the polynomial stage block is specified
before the cap gauge and its Schur elimination produces that gauge.  A
second direction is now exact at the algebraic level:
[Theorem 9.10](#v12-thm-minimal-schur-gauge-generation) proves that every rank-$`2q`$
cap deformation is the Schur reduction of exactly $`2q`$ auxiliary
alternating modes.  Because that generator is constructed from the target
deformation, it is not yet a selector.  The positive Stieltjes realization
in [Theorem 9.24](#v12-thm-stieltjes-interior-generation) is more canonical and produces
the Gauss sequence, but [Equation (168)](#v12-eq-corrected-stieltjes-obstruction) proves
that it cannot generate the full Taylor-corrected hierarchy.  A fully
coupled spatial–temporal cell in which independently specified mixed
bubble modes select $`\Theta`$, rather than encode it, therefore remains
open in this sharper form.
Extending [Equation (42)](#v12-eq-volume-cap-exact-sequence) to singular descriptor cells
requires composition of linear Dirac relations.  Extending it over a
mesh with nontrivial topology requires separating exact, closed, and
cohomological modes rather than labelling all transverse directions
inactive.

The mixed transmission-line deformation proves structural
nonseparability; [Table 17](#v12-tab-dispersion-selector) additionally demonstrates
a concrete, but problem-specific, numerical gain.  Its generic loss of
symplecticity is quantified by [Equation (141)](#v12-eq-mixed-gauss-symplectic).  A general
selector balancing dispersion, causality, sparsity, and conditioning
remains open.  The global convergence table varies $`h`$ at fixed spatial
degree; it is not a uniform-in-$`N`$ field convergence theorem.  Nor is
the scalar cap power identity by itself a local multisymplectic form
formula on variations; deriving such a law requires an additional tangent
cell analysis.  Finally,
the commutator test excludes the scalar functional calculus only after the
stated spatial chart and boundary causality are fixed.  It does not exclude
matrix-valued coefficients, splittings with additional operators, or
augmented problem-dependent formulations.

<a name="v12-section-12"></a>

## 12. Conclusion

The discrete fundamental theorem of calculus and spatial integration by
parts are boundary components of one affine Stokes problem on the
discrete cylinder.  The genuine total differential and the stacked
directional differential encode this statement in different linear maps.
For the former, a degree-respecting symmetric double has gauge
$`U_q^*\otimes V_q^*`$; for the latter, the intrinsic active freedom is
$`\Lambda^2(E_\mathrm{st}/\ker D_\mathrm{st})^*`$.  The canonical filtration
[Equation (35)](#v12-eq-canonical-sector-sequence) then identifies pure directions and
the mixed quotient without pretending that a chosen direct-sum splitting
is intrinsic.

A fixed regular interior solve transports this volume geometry to the
caps by the exact map $`\omega\mapsto L_h^*\omega`$.
[Equation (42)](#v12-eq-volume-cap-exact-sequence) determines all visible and invisible
volume gauges and proves when every temporal gauge possesses a volume
extension; [Equation (46)](#v12-eq-intrinsic-mixed-pullback) isolates its mixed part
canonically.  A trace-generated volume form can be inserted in the raw
weak cell, and [Theorem 4.14](#v12-thm-trace-generated-condensation) proves exact
recovery of the prescribed cap gauge.  This result is an exact
realization theorem, not a hidden selection principle.

Within the finite-dimensional, nonsingular, quadratic setting,
[Theorem 5.3](#v12-thm-gauge-cayley-chart) identifies alternating gauges bijectively
with the open subset of $`O(Q)`$ having no eigenvalue $`+1`$.  These are
constitutively shifted coordinates built from the classical Cayley
transform; the group identity is excluded and the zero gauge represents
midpoint.  The exact flow, Padé maps, and mixed deformations are points
in this cell; the symplectic centralizer is exactly the preimage of
$`O(Q)\cap Sp(\Omega)`$.  More importantly,
[Theorem 8.8](#v12-thm-active-gauss-cell) reverses the direction of construction:
starting from the geometric Gauss pairing on
$`\mathcal P_s\to\mathcal P_{s-1}`$, it eliminates $`sn`$ active stage
variables and derives both the diagonal Padé map and its constitutively
adapted cap gauge.  Thus Gauss–$`s`$ is no longer merely encoded in a
minimal two-point cell.

In space, [Theorem 9.2](#v12-thm-rectangular-stokes-dirac) proves that the raw
rectangular Stokes identity and an invertible effort/input reconstruction
force the complete square open Dirac representation, including
feedthrough.  Applied to the 2012 transmission line, it combines with the
active temporal cell and the canonical mixed lift
$`\mathcal T_N=\widetilde{\mathcal T}_N\pi_{tx}`$.  Boundary power,
dissipation, and interface cancellation remain exact.  At the
unrestricted algebraic level, every desired rank-$`2q`$ change of cap
gauge can be generated by
a minimal $`2q`$-dimensional alternating auxiliary block.  This proves
finite-dimensional interior generability, but the block still encodes the
target gauge and therefore does not constitute a field selector.  For every $`s`$,
the perturbation $`h^{2s-1}\Theta`$ is self-adjoint, energy preserving,
and has exactly order $`2s`$ unless
$`\Theta=\Theta_s^\star=c_s\Omega_NA_N^{2s}`$.  This unique fixed
coefficient within the one-term ansatz raises the order to exactly
$`2s+2`$; at $`s=1`$ it gives exactly Gauss–2.  It is central, but this
is not a high-order obstruction: a second arbitrary alternating term
retains order $`2s+2`$ and, when noncommuting, separates the map from the
scalar functional calculus of the semidiscrete line.  The central choice
$`\Xi_s^\star=d_s\Omega_NA_N^{2s+2}`$ raises the order by two more, and
successive Taylor coefficients generate the explicit rational hierarchy
[Theorem 9.18](#v12-thm-central-rational-hierarchy), with global order at least
$`2s+2r`$ and exactness conditional on $`\alpha_{s,r}\ne0`$, converging
locally to the exact-flow gauge.  Its reduced form
$`\Pi_{s,r}(-z)^{-1}\Pi_{s,r}(z)`$ separates removable $`+1`$ chart
resonances from true propagator poles: there are no reduced poles on the
imaginary axis, but the family is not generally $`A`$-stable.  The
Stieltjes representation of the exact coordinate proves the Padé status
of $`C_s`$ and gives each Gauss coordinate the sequence
[Equation (1)](#v12-eq-gauss-stieltjes-jacobi-chain), with positivity referring to the
measure, residues, resolvent, and Jacobi coefficients rather than to the
complete augmented port-Hamiltonian energy matrix.  The corrected
$`(2,1)`$ member is negative in the associated
Stieltjes variable for $`u>60`$, proving that this passive mechanism does
not extend to the whole Taylor-corrected hierarchy.  Strict non-annulation of all remaining coefficients is left as
a precisely stated conjecture.  The computations independently verify
rectangular reconstruction, trace-generated realization, active-stage
generation, stable Cayley evaluation, the coefficient $`d_s`$, the
one-term dichotomy, the noncommuting high-order branch—including its
$`s=3`$ multiprecision asymptotics—the reduced-pole distinction, exact
coefficient-sign evidence through $`s,j\le12`$, and dispersion design.
All numerical blocks are guarded by executable assertions.  The remaining
fundamental step is no longer algebraic existence of an interior generator;
it is a fully coupled field cell whose mixed deformation is selected by
spatial–temporal interior modes specified independently of the target
gauge, together with a uniform
space–time convergence and local multisymplectic analysis.

## Acknowledgment

The author thanks the port-Hamiltonian and geometric-integration communities
whose independent spatial, temporal, and field-theoretic constructions
motivate the unified algebraic viewpoint developed here.


---

<a name="extension-section-13"></a>

## 13. Purpose and precise scope

The version-12 manuscript classifies Stokes-compatible pairings and constructs quadratic conservative cap maps by eliminating active temporal stages. Its energy identity also accommodates dissipation, but the high-order conservative cap correction cannot simply be inserted into a dissipative one-effort step without losing order.

This extension retains the active stages. It adds an alternating correction to their projected efforts, proves the resulting port power law, and identifies the common-stage subclass that reduces **exactly** to the mixed cap family of version 12 in its quadratic conservative setting. Thus the new formulation extends that family while retaining the stage information needed for nonlinear energy and physical dissipation.

The zero-gauge time integrator is an established projected-gradient/energy-preserving collocation method [T2–T4]. Skew modifications of discrete-gradient operators are also established [T5–T7]. No priority claim is made for high-order energy-consistent integration itself. The contributions investigated here are the stage-gauge classification, its exact connection to the repository's cap family, and the obstruction caused by compressing dissipative effort variation into one effort.

Families of high-order energy-preserving continuous-stage Runge–Kutta methods [T9] and discontinuous Galerkin methods with interior effort variables [T10] are further relevant antecedents. A different gauge parametrization alone does not prove that the resulting integrator is different after elimination. The exact cap bridge below is a result within this repository's framework, not a claim that retaining interior efforts was previously unknown.

All convergence statements concern a fixed finite-dimensional spatial discretization and the regular local solution branch. The experiments use the repository's 2012 polynomial transmission-line reduction, including its physical boundary ports and skew feedthrough. They do not establish uniform convergence as space and time are refined simultaneously, shock resolution, descriptor consistency, general positivity, or universal symplecticity.

<a name="extension-section-14"></a>

## 14. Active time cell and constitutive effort

Consider the finite-dimensional PH realization

```math
\dot x=(J(x)-R(x))\nabla H(x)+B(x)u(t),
\qquad y=B(x)^\top\nabla H(x)+N(x)u(t),
```

where $`J^\top=-J`$, $`R^\top=R\succeq0`$ and $`N^\top=-N`$. The matrices may depend smoothly on the state. A resistive feedthrough can be treated separately by adding its symmetric part to the boundary dissipation.

On a time interval of length $`h>0`$, let $`c_i,b_i`$ be the $`s`$ Gauss nodes and positive weights on $`[0,1]`$, with $`\sum_i b_i=1`$. Let $`\ell_i`$ be their cardinal polynomials and set

```math
a_i(\tau)=\int_0^\tau\ell_i(\sigma)\,d\sigma,
\qquad
X(\tau)=x_0+h\sum_{i=1}^s a_i(\tau)F_i.
```

Here $`F_i`$ are physical time derivatives, not increments. Define the projected gradient

```math
G_i=\frac1{b_i}\int_0^1\ell_i(\tau)\nabla H(X(\tau))\,d\tau.
```

Write $`F,G,E`$ for the stage-stacked vectors and $`W=\operatorname{diag}(b_1,\ldots,b_s)\otimes I_n`$. The exact chain rule is

```math
H(X(1))-H(x_0)=hF^\top WG.
```

The independent temporal quadrature must integrate this constitutive projection accurately. For a polynomial Hamiltonian of degree $`\nu`$, the integrand defining $`G_i`$ has degree at most $`\nu s-1`$; $`q`$-point Gauss quadrature is exact if $`2q\geq\nu s`$. This requirement is distinct from the number $`s`$ of active stages. In particular, $`s=2,q=4`$ handles quartic energy exactly in exact arithmetic. Finite quadrature generally leaves a nonzero chain defect for nonpolynomial energy.

### Proposition T1 — Complete linear stage-gauge family

For a correction linear in $`F`$ and required to preserve the chain pairing for **every** stage vector, the admissible efforts are exactly

```math
E=G-W^{-1}K_hF,\qquad K_h^\top=-K_h.
```

For fixed cell data the family has dimension $`sn(sn-1)/2`$.

**Proof.** Write $`E-G=CF`$. The required identity is $`F^\top WCF=0`$ for all $`F`$. Polarization gives $`WC+C^\top W=0`$. Thus $`K_h=-WC`$ is alternating, and the converse follows immediately. This classifies linear corrections on the full stage space; it does not assert that every nonlinear correction has one fixed matrix representation. At a nonzero particular $`F`$, nonlinear corrections instead belong to the annihilator of $`F`$ under the weighted pairing. ∎

A state-dependent matrix $`K_h(F,X)`$ still gives exact power pointwise if it is alternating at every evaluation. Its regularity, step-size scaling and reversal properties must be checked separately.

<a name="extension-section-15"></a>

## 15. Closure, ports and exact energy balance

The stage equations are

```math
F_i=(J(X_i)-R(X_i))E_i+B(X_i)U_i,
\qquad
Y_i=B(X_i)^\top E_i+N(X_i)U_i,
```

where $`X_i=X(c_i)`$ and $`U_i=u(t_0+c_i h)`$. The updated state is $`x_1=x_0+h\sum_i b_iF_i`$.

### Theorem T2 — Gauged PH power law

With exact constitutive integration and an exact stage solve, every alternating stage gauge satisfies

```math
H(x_1)-H(x_0)
=h\sum_i b_iU_i^\top Y_i-h\sum_i b_iE_i^\top R(X_i)E_i.
```

**Proof.** Alternation gives $`F^\top WG=F^\top WE+F^\top K_hF=F^\top WE`$. Substitute the stage equations; the terms with $`J`$ and $`N`$ vanish. ∎

Consequently closed conservative systems preserve $`H`$, and unforced dissipative systems have nonincreasing $`H`$. This is an identity for the computed stage efforts. Substituting $`\nabla H(X_i)`$ for $`E_i`$ in the output generally destroys it.

### Computable certificate

Let $`r_i=F_i-[(J_i-R_i)E_i+B_iU_i]`$ and

```math
\delta_{\rm chain}=H(x_1)-H(x_0)-hF^\top WG.
```

For symmetric $`R`$ and skew feedthrough, direct algebra gives

```math
\Delta H-\mathcal W+\mathcal D
=\delta_{\rm chain}+hF^\top K_hF
+h\sum_i b_iE_i^\top J_iE_i
+h\sum_i b_iE_i^\top r_i,
```

with $`\mathcal W=h\sum_i b_iU_i^\top Y_i`$ and $`\mathcal D=h\sum_i b_iE_i^\top R_iE_i`$. The implementation records these contributions separately, including their floating-point closure. It rejects a supplied nonalternating gauge rather than silently replacing it by its skew part. Newton failure raises an exception.

<a name="extension-section-16"></a>

## 16. Exact recovery of the version-12 mixed cap family

The common-stage gauge is

```math
K_h=h^{2s}(bb^\top\otimes\Theta),
\qquad \Theta^\top=-\Theta,
\qquad b=(b_1,\ldots,b_s)^\top.
```

It produces the same correction at each stage:

```math
E_i=G_i-h^{2s}\Theta v,
\qquad v=\sum_jb_jF_j=\frac{x_1-x_0}{h}.
```

### Theorem T3 — Active-stage extension of the mixed Gauss cap map

Suppose $`H(x)=x^\top Qx/2`$, $`J,Q`$ are constant and nonsingular, and $`R=0`$, $`u=0`$. Wherever the cap chart is regular, eliminating the active stages gives exactly the version-12 cap gauge

```math
S_{s,h}=S^{\rm G}_{s,h}+h^{2s-1}\Theta.
```

The equality holds for arbitrary alternating $`\Theta`$, including coefficients that do not commute with the conservative generator.

**Proof.** Gaussian exactness gives $`G_i=QX_i`$: the projection integrand has degree at most $`2s-1`$. Put $`A=JQ`$ and let $`R_s(z)`$ be the diagonal Padé stability function of Gauss–$`s`$. Define

```math
\varphi_s(z)=\frac{R_s(z)-1}{z},
\qquad \varphi_s(0)=1.
```

For the constant stage source $`-h^{2s}J\Theta v`$, elimination of the Gauss stage equations gives

```math
\Delta x=(R_s(hA)-I)x_0-h^{2s}\varphi_s(hA)J\Theta\Delta x.
```

Let $`B_-^{\rm G}=I+hJS^{\rm G}_{s,h}-hA/2`$ and $`B_+^{\rm G}=B_-^{\rm G}+hA`$. The original cap identity implies

```math
B_-^{\rm G}(R_s(hA)-I)=hA,
\qquad B_-^{\rm G}\varphi_s(hA)=I.
```

Multiplication therefore yields

```math
(B_-^{\rm G}+h^{2s}J\Theta)\Delta x=hAx_0,
```

and hence the cap propagator

```math
M_h=(B_-^{\rm G}+h^{2s}J\Theta)^{-1}
     (B_+^{\rm G}+h^{2s}J\Theta),
```

which is precisely the version-12 mixed deformation. No commutation of $`\Theta`$ with $`A`$ was used. ∎

The theorem supplies the exact bridge to the original family. For general nonlinear energy or nonzero dissipation, the active-stage equations remain the definition of the extension. Replacing all stage efforts by a single effective effort would lose information, as the next section proves.

<a name="extension-section-17"></a>

## 17. The dissipative obstruction to one-effort compression

### Theorem T4 — Temporal effort variance obstructs order above two

Consider a closed smooth PH system $`\dot x=L\nabla H(x)`$, with constant invertible $`L=J-R`$, $`J^\top=-J`$ and $`R^\top=R\succeq0`$. A method constrained to satisfy simultaneously

```math
\Delta x=hLE,\qquad
H(x_1)-H(x_0)=-hE^\top RE
```

with one effort $`E`$ cannot have order at least three at a state for which

```math
\dot e(0)^\top R\dot e(0)>0,\qquad e(t)=\nabla H(x(t)).
```

The statement assumes the usual smooth local error expansion and does not cover methods that change the physical resistance, introduce additional dissipative channels, or use more than one effort.

**Proof.** On the exact endpoints, the first equation forces

```math
E=\bar e=\frac1h\int_0^h e(t)\,dt.
```

But the exact energy loss is

```math
H(x(h))-H(x(0))
=-h\bar e^\top R\bar e
-\int_0^h(e(t)-\bar e)^\top R(e(t)-\bar e)\,dt.
```

For $`e\in C^2`$, the last integral equals

```math
\frac{h^3}{12}\dot e(0)^\top R\dot e(0)+O(h^4).
```

An order-$`p`$ method with $`p\geq3`$ changes the exact endpoint by $`O(h^{p+1})`$ and its forced effort by $`O(h^p)`$. Both corresponding changes to the energy balance are $`O(h^{p+1})=O(h^4)`$, which cannot cancel the nonzero cubic variance term. ∎

This is a restriction on an unmodified physical one-effort closure, not an impossibility theorem for high-order dissipative discrete-gradient methods. Such methods exist, including approaches that modify and complete the dissipative operator [T7].

### Quadratic coefficient check

For $`H=\|x\|^2/2`$, $`L=J-R`$, and a cap effort $`E=\bar x-S_h\Delta x`$ with $`S_h=hS_1+O(h^3)`$, expansion gives

```math
M_h-e^{hL}
=h^3\left(\frac{L^3}{12}-LS_1L\right)+O(h^4).
```

If $`L`$ is invertible, cancellation requires $`S_1=L/12`$, whose symmetric part is $`-R/12`$. It cannot be alternating when $`R\ne0`$. This reproduces the obstruction without a trajectory argument.

### Proposition T5 — Minimal information in stage dissipation

For constant $`R`$, the stage mean $`\bar E=\sum_i b_iE_i`$ satisfies

```math
\sum_i b_i E_i^\top RE_i
=\bar E^\top R\bar E
+\sum_i b_i(E_i-\bar E)^\top R(E_i-\bar E).
```

The discarded variance has matrix

```math
V=(\operatorname{diag}(b)-bb^\top)\otimes R,
\qquad \operatorname{rank}V=(s-1)\operatorname{rank}R.
```

**Proof.** Expand the right-hand side and use $`\sum_i b_i=1`$. The temporal matrix is congruent to $`I-\sqrt b\,\sqrt b^\top`$, a rank-$`(s-1)`$ orthogonal projector. The rank formula follows from the Kronecker product. ∎

Any exact representation of this variance by a sum of squares of linear scalar channels needs at least $`\operatorname{rank}V`$ channels, since a factorization $`V=C^\top C`$ has $`\operatorname{rank}C\leq`$ number of rows. An alternating gauge cannot restore this missing symmetric quadratic form.

<a name="extension-section-18"></a>

## 18. Elimination preserves passivity, including the ports

Normalize the weighted variables by $`\widehat F=W^{1/2}F`$, $`\widehat E=W^{1/2}E`$, $`\widehat G=W^{1/2}G`$ and $`\widehat K=W^{-1/2}K_hW^{-1/2}`$. Port variables use the corresponding weights with the port dimension. Let the block-diagonal stage matrices be denoted by $`\mathcal J,\mathcal R,\mathcal B,\mathcal N`$.

Define the full flow and effort vectors and their structure matrices by

```math
f=\begin{bmatrix}\widehat F\\-\widehat Y\end{bmatrix},
\quad
g=\begin{bmatrix}\widehat G\\\widehat U\end{bmatrix},
\quad
Z=\begin{bmatrix}\mathcal J-\mathcal R&\mathcal B\\-\mathcal B^\top&-\mathcal N\end{bmatrix},
\quad
D=\begin{bmatrix}\widehat K&0\\0&0\end{bmatrix}.
```

The gauged closure is $`f=Z(g-Df)`$. Whenever $`T=I+ZD`$ is invertible,

```math
f=Z_{\rm eff}g,\qquad Z_{\rm eff}=T^{-1}Z.
```

### Proposition T6 — Dissipative congruence

```math
\operatorname{sym}Z_{\rm eff}
=T^{-1}\operatorname{sym}Z\,T^{-\top}
=-T^{-1}
\begin{bmatrix}\mathcal R&0\\0&0\end{bmatrix}
T^{-\top}\preceq0.
```

**Proof.** Multiply $`Z_{\rm eff}+Z_{\rm eff}^\top`$ by $`T`$ and $`T^\top`$. The additional terms cancel because $`D^\top=-D`$. ∎

For $`\mathcal R=0`$ the eliminated structure is skew. Its graph is isotropic for the split power pairing and has half the ambient dimension, so it is maximally isotropic. Thus the conservative algebraic reduction gives an open Dirac structure; the nonlinear Hamiltonian enters through a separate constitutive closure.

The elimination generally changes input, output and feedthrough blocks together. Keeping only a modified state matrix and retaining the old output is not the same PH system.

For fixed matrices, strictly positive $`\mathcal R`$ ensures invertibility of the state block for every alternating $`\widehat K`$. Indeed, a nonzero null vector would give $`\widehat E=-\widehat K\widehat F`$, $`\widehat F=(\mathcal J-\mathcal R)\widehat E`$, hence $`0=\widehat E^\top\widehat F=-\widehat E^\top\mathcal R\widehat E`$. Semidefinite resistance does not give this guarantee: with no resistance and $`\mathcal J=\widehat K`$ equal to the canonical $`2\times2`$ skew matrix, $`I+\mathcal J\widehat K=0`$. This algebraic fact is not a global existence or uniqueness theorem for the nonlinear stage solve.

<a name="extension-section-19"></a>

## 19. Order, reversal and coordinate covariance

### Proposition T7 — Order inherited by a sufficiently small stage gauge

Assume the zero-gauge scheme has local error $`O(h^{p+1})`$. Assume the stage residual, coefficients and gauge are sufficiently smooth jointly in $`(h,x_0,F)`$ near $`h=0`$, with uniform derivative bounds on the neighborhood considered, and the zero-gauge stage Jacobian has a uniformly bounded inverse near the consistent branch. If $`K_h=O(h^p)`$ uniformly there, the gauged scheme has local error $`O(h^{p+1})`$ and global order at least $`p`$ on a fixed smooth finite-time problem.

**Proof.** Substituting the zero-gauge stages into the gauged equations leaves residual $`O(h^p)`$. The local implicit-function estimate gives a stage change $`O(h^p)`$. The endpoint multiplies this change by $`h`$. Standard local Lipschitz stability and discrete Gronwall then give the global estimate. For a smooth vanishing gauge, local solvability follows directly at $`h=0`$: the residual Jacobian with respect to the derivative stages is the identity. ∎

For the known Gauss projected-gradient base of order $`2s`$ under its smoothness and quadrature assumptions, the sufficient choice is $`K_h=O(h^{2s})`$. The common-stage gauge uses exactly this scaling. This proposition inherits an order; it does not prove automatic order elevation for dissipative or nonlinear systems from the conservative coefficients of version 12.

Antisymmetry alone is insufficient. For example, Gauss–2 with $`H=\|x\|^2/2`$, $`J^2=-I`$, and $`K_h=h^2(\operatorname{diag}(b)\otimes J)`$ yields the map $`R_2(hJ/(1-h^2))`$. It preserves energy and is reversible, but is only second order.

### Reversal

Let $`P`$ reverse the order of Gauss stages and leave state coordinates unchanged. For an autonomous model, a state-independent gauge gives a self-adjoint method if

```math
K_{-h}=PK_hP^\top.
```

The stage derivatives are physical derivatives; they are permuted, not negated, when the interval is traversed with negative $`h`$. The common-stage gauge satisfies the condition because $`Pb=b`$ and $`h^{2s}`$ is even. For state-dependent gauges, the equality must hold with the full reversed path and exchanged endpoints. Dissipation does not prevent algebraic self-adjointness; the backward step is not a forward dissipative evolution.

For even powers of $`h`$, admissible alternating matrices preserve the two eigenspaces of $`P`$. Their dimension is $`n_+(n_+-1)/2+n_-(n_--1)/2`$, where $`n_+=n\lceil s/2\rceil`$ and $`n_-=n\lfloor s/2\rfloor`$. Odd powers exchange those eigenspaces and have dimension $`n_+n_-`$.

### Linear coordinates

For $`z=Tx`$, $`T`$ nonsingular, use

```math
\widetilde H(z)=H(T^{-1}z),\quad
\widetilde J=TJT^\top,\quad
\widetilde R=TRT^\top,\quad
\widetilde B=TB.
```

The gauge must transform as a covariant two-form:

```math
\widetilde K_h=(I_s\otimes T^{-\top})K_h(I_s\otimes T^{-1}).
```

Substitution shows that every exact stage solution transforms accordingly. The method therefore commutes with the nodal/modal changes of spatial coordinates used by the 2012 discretization. Recomputing a Euclidean normalization after a nonorthogonal coordinate change would instead define a different gauge.

<a name="extension-section-20"></a>

## 20. The spatial reduction and a geometrically specified test gauge

The experiments import the unchanged version-12 rectangular transmission-line construction. Efforts lie in $`\mathcal P_N`$, flow/state coefficients in $`\mathcal P_{N-1}`$, and the reconstruction fixes the same two physical boundary inputs. The rectangular identity implies the reduced skew $`J_N`$, dual input/output matrices and skew feedthrough $`N_N`$. No antisymmetrization is applied to repair these raw spatial blocks.

In the orthonormal Legendre state coordinates, write $`D_{\rm a}=DM`$, the derivative restricted to the active polynomial coefficients. The test gauge is a normalized version of

```math
\Theta_0=\begin{bmatrix}0&D_{\rm a}\\-D_{\rm a}^\top&0\end{bmatrix}.
```

Thus the correction depends on a specified spatial derivative and on temporal stage flows. It is nonzero and is not fitted to a target cap propagator or to the convergence reference. The normalization and amplitude are recorded in the executable experiment.

This is a concrete geometric choice within the admissible family. It is not a proof that the original PDE uniquely selects this gauge, nor that an independent physical bubble dynamics selects it. The version-12 distinction between compatible design, algebraic realization and independent field selection remains in force.

<a name="extension-section-21"></a>

## 21. Interconnection

On a common time grid, impose opposing stage port interconnections. The two terms $`h\sum_i b_iU_i^\top Y_i`$ then cancel exactly when adding the cell balances, regardless of their internal alternating gauges.

For different grids, let $`W_A,W_B`$ be the time-weight matrices including local step lengths, and let $`C`$ contain the common-refinement integrals of the two temporal port bases. The adjoint transfers

```math
U_A=-W_A^{-1}CY_B,\qquad
U_B=W_B^{-1}C^\top Y_A
```

give $`U_A^\top W_AY_A+U_B^\top W_BY_B=0`$. The stage-gauge power law is therefore compatible with this interface construction. This is an algebraic extension; the current validation campaign does not claim a new multirate convergence theorem or a computational speedup.

<a name="extension-section-22"></a>

## 22. Numerical evidence and reproduction

The numerical tables below are generated from the committed result files. The complete protocol, failures raised by the nonlinear solver, reference tolerances and model parameters are retained with the code.

Run from the repository root, using the dependencies in `requirements.txt`:

```bash
make verify-extension
```

The three programs respectively check the implementation, the independent algebraic identities and counterexamples, and the nonlinear transmission-line dynamics. `make verify` additionally runs the unchanged 28 version-12 assertion groups. The reported extension run used Python 3.9.13, NumPy 1.26.4, SciPy 1.13.1 and float64 arithmetic; environment details are included in the result files.

### Nonlinear transmission line

The Hamiltonian is the spatial integral of $`(q^2+p^2)/2+\alpha q^4/4`$ on the unit interval, with $`\alpha=2`$. The constitutive spatial Gauss rule has at least $`2N+2`$ points, sufficient for exact integration of this polynomial energy. All trajectory errors below are against an independently integrated, fixed-$`N`$ DOP853 solution (`rtol=2.5e-13`, `atol=2.5e-15`). A second reference with looser tolerances checks reference sensitivity. These are temporal, not spatial-refinement, measurements.

For the closed nonlinear case, $`N=6`$, $`T=1`$, the initial field is the projection of $`q(z)=0.35\sin(\pi z/2)`$ with $`p=0`$. The last refinement is $`h=0.0125`$ to $`0.00625`$.

| Scheme | Active stages / constitutive quadrature | Last observed order | Final spatial $`L^2`$ error |
|---|---:|---:|---:|
| Projected-gradient base | 1 / 2 | 1.9395 | $`9.1594\times10^{-6}`$ |
| Same base with spatial gauge | 1 / 2 | 1.8759 | $`1.0986\times10^{-5}`$ |
| Projected-gradient base | 2 / 4 | 3.9684 | $`8.8556\times10^{-9}`$ |
| Same base with spatial gauge | 2 / 4 | 3.9685 | $`8.8944\times10^{-9}`$ |

The gauge has unit Frobenius norm in the specified orthonormal spatial chart and the prescribed $`h^{2s}`$ scaling. It is not tuned to lower the trajectory error; indeed, these results do not show an accuracy improvement over the base scheme. Both gauged and ungauged runs retain the energy to floating-point accuracy with sufficiently accurate constitutive quadrature.

![Nonlinear temporal convergence](../results/figures/ph_extension_order.svg)

For the open dissipative case, $`N=8`$, $`R=0.12I`$, $`T=4`$, the initial amplitude is $`0.12`$, and the physical boundary inputs are $`(0.06\sin(1.7t),0.04\sin(2.3t))`$. The tested steps are $`0.04,0.02,0.01,0.005,0.0025`$.

| Scheme | Last observed order | Final spatial $`L^2`$ error at $`h=0.0025`$ |
|---|---:|---:|
| Two-stage projected-gradient base | 3.98584 | $`1.57420\times10^{-7}`$ |
| Same base with spatial gauge | 3.98584 | $`1.57437\times10^{-7}`$ |

The two reference solutions differ by $`7.97\times10^{-12}`$ at the endpoint. Over all forced runs, the maximum cumulative power defect is $`1.29\times10^{-16}`$, the maximum one-step balance defect is $`2.85\times10^{-17}`$, and reconstruction of the physical boundary outputs differs by at most $`2.23\times10^{-16}`$. The integrated density equations have maximum residual $`6.66\times10^{-14}`$. The unforced dissipative run has no positive energy increment and loses $`0.0179689`$ units of energy, equal to the accumulated discrete dissipation to roundoff.

The coarse forced steps initially show rates below two; refinement resolves the expected fourth-order regime. Their small power defects alone would not have revealed this trajectory error.

![Forced energy balance and independent trajectory comparison](../results/figures/ph_extension_open_damped.svg)

### Independent geometry checks and expected counterexamples

All 14 engine self-checks pass. Direct comparison with the unchanged version-12 cap implementation at one, two and three active stages gives map errors below $`4.03\times10^{-16}`$ for a noncommuting gauge. Nonlinear reversal has error $`1.08\times10^{-15}`$; covariance under a nonorthogonal coordinate change, including state-dependent structure and physical outputs, has errors below $`5.4\times10^{-16}`$.

All 18 checks in `verify_stage_geometry.py` pass. The full state/port elimination satisfies the dissipative congruence to $`4.68\times10^{-15}`$; keeping the old output after elimination instead produces a power error of $`1.73`$ in the same test. The variance form has the predicted rank four for three stages and rank-two resistance.

For the damped oscillator with $`J^2=-I`$, $`Q=I`$ and $`R=0.2I`$, the conservative fourth-order cap correction drops to observed order $`1.99999`$ under single-effort physical dissipation. The active two-stage scheme with a nonzero common gauge gives $`3.99995`$. Both still have one-step power defects below $`2\times10^{-16}`$: energy balance does not distinguish their orders.

Two further negative controls prevent overinterpreting preservation:

- Using one instead of two constitutive quadrature points on a quartic one-stage trajectory produces a maximum relative energy drift $`2.98\times10^{-3}`$; sufficient quadrature reduces it to $`2.64\times10^{-15}`$.
- An unresolved linear spatial mode with $`\omega h=2.301`$ has a final phase lag $`1.274`$ radians despite relative energy drift below $`3\times10^{-15}`$.

The deliberately underscaled alternating gauge from Section 19 has observed order $`2.00068`$, while the $`h^4`$ gauge gives $`4.00000`$. The noncommuting spatial gauges preserve the quadratic energy but have nonzero symplectic defects, consistent with the distinction already made in version 12.

### Machine-readable evidence

- [Engine self-checks](../results/ph_gauges_selfcheck.json).
- [Independent geometry and compression checks](../results/stage_geometry_validation.json).
- [Complete nonlinear/open-system results](../results/ph_extension_validation.json) and [tabular export](../results/tables/ph_extension_validation.csv).

All numerical claims here concern the specified examples. Exact power identities and the stated algebraic results follow from their proofs; numerical agreement is an implementation check, not their proof.

<a name="extension-section-23"></a>

## 23. What remains open

The results establish an active-stage extension with exact discrete power and a proved bridge to the conservative quadratic cap family. They leave several substantive questions:

- A useful gauge-selection principle independent of a prescribed target propagator, with demonstrated benefits beyond an individual test.
- Nonlinear order elevation by nonzero gauges, rather than only retention of the base order.
- Uniform approximation and stability estimates in both spatial degree and time step.
- General descriptor systems, state-dependent mass operators, shocks and positivity constraints.
- An additional tangent analysis for local multisymplectic conservation; scalar energy balance is not such a law.
- Efficient sparse or matrix-free implementations of the implicit stage equations.

Energy accuracy and trajectory accuracy are reported separately. Neither a zero power residual nor a passed algebraic test establishes phase accuracy, PDE convergence or bibliographic novelty.


## Reproducibility and provenance

The original 28-group algebraic verification remains in
[verify_transmission_line_v12.py](../experiments/verify_transmission_line_v12.py),
with its [archived output](../results/verification_report_v12.txt).
The executable active-stage extension adds 14 engine checks, 18 independent
geometry checks and the nonlinear simulation assertions described in Section 22.
The [manifest](../results/manifest.json) records the numerical environment and
SHA-256 hashes of the computational sources and recorded evidence.
The original report and the extension report come from distinct runs; their
floating-point residuals need not coincide across numerical environments.

The version-12 editorial distinctions are retained: an auxiliary alternating
Schur realization is not a field-based selection principle, and a positive
Jacobi realization of a Stieltjes coordinate is not positivity of the full
augmented PH structure. The extension does not claim a new projected-gradient
base integrator, a universal accuracy benefit from its chosen gauge, or a
simultaneous space-time convergence theorem.

The repository is presented in Markdown with vector figures. The original PDF
and earlier split documentation remain recoverable through Git history.

## References

### Geometric foundations and original manuscript


<a name="v12-ref-basubose1983"></a>

**[1]** 
S. Basu and N. K. Bose,
 Matrix Stieltjes series and network models,
 *SIAM Journal on Mathematical Analysis* 14 (1983),
209–222.
 [doi:10.1137/0514018](https://doi.org/10.1137/0514018).

<a name="v12-ref-bridges2001"></a>

**[2]** 
T. J. Bridges and S. Reich,
 Multi-symplectic integrators: numerical schemes for Hamiltonian
PDEs that conserve symplecticity,
 *Physics Letters A* 284 (2001), 184–193.
 [doi:10.1016/S0375-9601(01)00294-8](https://doi.org/10.1016/S0375-9601(01)00294-8).

<a name="v12-ref-brugnano2010"></a>

**[3]** 
L. Brugnano, F. Iavernaro, and D. Trigiante,
 Hamiltonian boundary value methods (energy preserving discrete
line integral methods),
 *Journal of Numerical Analysis, Industrial and Applied
Mathematics* 5 (2010), 17–37.
 <https://arxiv.org/abs/0910.3621>.

<a name="v12-ref-brugnoli2022"></a>

**[4]** 
A. Brugnoli, R. Rashad, and S. Stramigioli,
 Dual field structure-preserving discretization of
port-Hamiltonian systems using finite element exterior calculus,
 *Journal of Computational Physics* 471 (2022), 111601.
 [doi:10.1016/j.jcp.2022.111601](https://doi.org/10.1016/j.jcp.2022.111601).

<a name="v12-ref-celledoni2017"></a>

**[5]** 
E. Celledoni and E. H. Høiseth,
 Energy-preserving and passivity-consistent numerical
discretization of port-Hamiltonian systems,
 arXiv:1706.08621, 2017.
 <https://arxiv.org/abs/1706.08621>.

<a name="v12-ref-celledoni2021"></a>

**[6]** 
E. Celledoni and J. Jackaman,
 Discrete conservation laws for finite element discretisations of
multisymplectic PDEs,
 *Journal of Computational Physics* 444 (2021), 110520.
 [doi:10.1016/j.jcp.2021.110520](https://doi.org/10.1016/j.jcp.2021.110520).

<a name="v12-ref-cheng2026"></a>

**[7]** 
P. Cheng,
 A high-order multi-symplectic finite element method for nonlinear
wave equations,
 *Physica D: Nonlinear Phenomena* 496 (2026), 135329.
 [doi:10.1016/j.physd.2026.135329](https://doi.org/10.1016/j.physd.2026.135329).

<a name="v12-ref-cooper1987"></a>

**[8]** 
G. J. Cooper,
 Stability of Runge–Kutta methods for trajectory problems,
 *IMA Journal of Numerical Analysis* 7 (1987), 1–13.
 [doi:10.1093/imanum/7.1.1](https://doi.org/10.1093/imanum/7.1.1).

<a name="v12-ref-cieslinski2013"></a>

**[9]** 
J. L. Cieśliński,
 Locally exact modifications of discrete gradient schemes,
 arXiv:1304.6533, 2013.
 <https://arxiv.org/abs/1304.6533>.

<a name="v12-ref-eidnes2022"></a>

**[10]** 
S. Eidnes,
 Order theory for discrete gradient methods,
 *BIT Numerical Mathematics* 62 (2022), 1207–1255.
 [doi:10.1007/s10543-022-00909-z](https://doi.org/10.1007/s10543-022-00909-z).

<a name="v12-ref-ge1988"></a>

**[11]** 
Z. Ge and J. E. Marsden,
 Lie–Poisson Hamilton–Jacobi theory and Lie–Poisson integrators,
 *Physics Letters A* 133 (1988), 134–139.
 [doi:10.1016/0375-9601(88)90773-6](https://doi.org/10.1016/0375-9601(88)90773-6).

<a name="v12-ref-giesselmann2025"></a>

**[12]** 
J. Giesselmann, A. Karsai, and T. Tscherpel,
 Energy-consistent Petrov–Galerkin time discretization of
port-Hamiltonian systems,
 *The SMAI Journal of Computational Mathematics*
11 (2025), 335–367.
 [doi:10.5802/smai-jcm.127](https://doi.org/10.5802/smai-jcm.127).

<a name="v12-ref-gonzalez1996"></a>

**[13]** 
O. Gonzalez,
 Time integration and discrete Hamiltonian systems,
 *Journal of Nonlinear Science* 6 (1996), 449–467.
 [doi:10.1007/BF02440162](https://doi.org/10.1007/BF02440162).

<a name="v12-ref-hairer2006"></a>

**[14]** 
E. Hairer, C. Lubich, and G. Wanner,
 *Geometric Numerical Integration: Structure-Preserving
Algorithms for Ordinary Differential Equations*, 2nd ed.,
 Springer, Berlin, 2006.
 [doi:10.1007/3-540-30666-8](https://doi.org/10.1007/3-540-30666-8).

<a name="v12-ref-hairer2010"></a>

**[15]** 
E. Hairer,
 Energy-preserving variant of collocation methods,
 *Journal of Numerical Analysis, Industrial and Applied
Mathematics* 5 (2010), 73–84.
 <https://www.unige.ch/ hairer/preprints/crete.html>.

<a name="v12-ref-kemmochi2023"></a>

**[16]** 
T. Kemmochi,
 Higher order discrete gradient method by the discontinuous
Galerkin time-stepping method,
 arXiv:2308.02334, 2023.
 <https://arxiv.org/abs/2308.02334>.

<a name="v12-ref-kotyczka2018"></a>

**[17]** 
P. Kotyczka, B. Maschke, and L. Lefèvre,
 Weak form of Stokes–Dirac structures and geometric
discretization of port-Hamiltonian systems,
 *Journal of Computational Physics* 361 (2018), 442–476.
 [doi:10.1016/j.jcp.2018.02.006](https://doi.org/10.1016/j.jcp.2018.02.006).

<a name="v12-ref-kotyczka2019"></a>

**[18]** 
P. Kotyczka and L. Lefèvre,
 Discrete-time port-Hamiltonian systems: a definition based on
symplectic integration,
 *Systems & Control Letters* 133 (2019), 104530.
 [doi:10.1016/j.sysconle.2019.104530](https://doi.org/10.1016/j.sysconle.2019.104530).

<a name="v12-ref-kostant2001"></a>

**[19]** 
B. Kostant and P. W. Michor,
 The generalized Cayley map from an algebraic group to its Lie
algebra,
 in *The Orbit Method in Geometry and Physics*, Progress in
Mathematics 213, Birkhäuser, 2003, 259–296.
 <https://arxiv.org/abs/math/0109066>.

<a name="v12-ref-moulla2012"></a>

**[20]** 
R. Moulla, L. Lefèvre, and B. Maschke,
 Pseudo-spectral methods for the spatial symplectic reduction of
open systems of conservation laws,
 *Journal of Computational Physics* 231 (2012), 1272–1292.
 [doi:10.1016/j.jcp.2011.10.008](https://doi.org/10.1016/j.jcp.2011.10.008).

<a name="v12-ref-mclachlan1999"></a>

**[21]** 
R. I. McLachlan, G. R. W. Quispel, and N. Robidoux,
 Geometric integration using discrete gradients,
 *Philosophical Transactions of the Royal Society A*
357 (1999), 1021–1045.
 [doi:10.1098/rsta.1999.0363](https://doi.org/10.1098/rsta.1999.0363).

<a name="v12-ref-mclachlan2020"></a>

**[22]** 
R. I. McLachlan and A. Stern,
 Multisymplecticity of hybridizable discontinuous Galerkin
methods,
 *Foundations of Computational Mathematics* 20 (2020),
35–69.
 [doi:10.1007/s10208-019-09415-1](https://doi.org/10.1007/s10208-019-09415-1).

<a name="v12-ref-marsden1998"></a>

**[23]** 
J. E. Marsden, G. W. Patrick, and S. Shkoller,
 Multisymplectic geometry, variational integrators, and nonlinear
PDEs,
 *Communications in Mathematical Physics* 199 (1998),
351–395.
 [doi:10.1007/s002200050505](https://doi.org/10.1007/s002200050505).

<a name="v12-ref-maslovskaya2024"></a>

**[24]** 
B. Wembe, C. Offen, S. Maslovskaya, S. Ober-Blöbaum, and P. Singh,
 Commutator-free Cayley methods,
 *Journal of Computational and Applied Mathematics* 477
(2026), 117184.
 [doi:10.1016/j.cam.2025.117184](https://doi.org/10.1016/j.cam.2025.117184).

<a name="v12-ref-miyatake2016"></a>

**[25]** 
Y. Miyatake and J. C. Butcher,
 A characterization of energy-preserving methods and the
construction of parallel integrators for Hamiltonian systems,
 *SIAM Journal on Numerical Analysis* 54 (2016), 1993–2013.
 [doi:10.1137/15M1020861](https://doi.org/10.1137/15M1020861).

<a name="v12-ref-quispel2008"></a>

**[26]** 
G. R. W. Quispel and D. I. McLaren,
 A new class of energy-preserving numerical integration methods,
 *Journal of Physics A: Mathematical and Theoretical*
41 (2008), 045206.
 [doi:10.1088/1751-8113/41/4/045206](https://doi.org/10.1088/1751-8113/41/4/045206).

<a name="v12-ref-reich2000"></a>

**[27]** 
S. Reich,
 Multi-symplectic Runge–Kutta collocation methods for Hamiltonian
wave equations,
 *Journal of Computational Physics* 157 (2000), 473–499.
 [doi:10.1006/jcph.1999.6372](https://doi.org/10.1006/jcph.1999.6372).

<a name="v12-ref-sanzserna1988"></a>

**[28]** 
J. M. Sanz-Serna,
 Runge–Kutta schemes for Hamiltonian systems,
 *BIT Numerical Mathematics* 28 (1988), 877–883.
 [doi:10.1007/BF01954907](https://doi.org/10.1007/BF01954907).

<a name="v12-ref-seslija2012"></a>

**[29]** 
M. Seslija, A. J. van der Schaft, and J. M. A. Scherpen,
 Discrete exterior geometry approach to structure-preserving
discretization of distributed-parameter port-Hamiltonian systems,
 *Journal of Geometry and Physics* 62 (2012), 1509–1531.
 [doi:10.1016/j.geomphys.2012.02.006](https://doi.org/10.1016/j.geomphys.2012.02.006).

<a name="v12-ref-schulze2023"></a>

**[30]** 
P. Schulze,
 Structure-preserving time discretization of port-Hamiltonian
systems via discrete gradient pairs,
 arXiv:2311.00403, 2023.
 <https://arxiv.org/abs/2311.00403>.

<a name="v12-ref-stern2026"></a>

**[31]** 
A. Stern and E. Zampa,
 Multisymplecticity in finite element exterior calculus,
 *Foundations of Computational Mathematics* 26 (2026),
1809–1858.
 [doi:10.1007/s10208-025-09720-y](https://doi.org/10.1007/s10208-025-09720-y).

<a name="v12-ref-sternzampa2026"></a>

**[32]** 
A. Stern and E. Zampa,
 Finite element exterior calculus for time-dependent Hamiltonian
partial differential equations,
 arXiv:2601.00103, 2026.
 <https://arxiv.org/abs/2601.00103>.

<a name="v12-ref-vanderschaft2002"></a>

**[33]** 
A. J. van der Schaft and B. M. Maschke,
 Hamiltonian formulation of distributed-parameter systems with
boundary energy flow,
 *Journal of Geometry and Physics* 42 (2002), 166–194.
 [doi:10.1016/S0393-0440(01)00083-3](https://doi.org/10.1016/S0393-0440(01)00083-3).

<a name="v12-ref-vankerschaver2012"></a>

**[34]** 
J. Vankerschaver, H. Yoshimura, and M. Leok,
 The Hamilton–Pontryagin principle and multi-Dirac structures for
classical field theories,
 *Journal of Mathematical Physics* 53 (2012), 072903.
 [doi:10.1063/1.4731481](https://doi.org/10.1063/1.4731481).

<a name="v12-ref-wall1948"></a>

**[35]** 
H. S. Wall,
 *Analytic Theory of Continued Fractions*,
 Chelsea, New York, 1948.

<a name="v12-ref-wynn1968"></a>

**[36]** 
P. Wynn,
 Upon the Padé table derived from a Stieltjes series,
 *SIAM Journal on Numerical Analysis* 5 (1968), 805–834.
 [doi:10.1137/0705060](https://doi.org/10.1137/0705060).

### Active-stage extension

Citations marked T refer to the following entries.


<a name="extension-ref-1"></a>

**[T1]** R. Moulla, L. Lefèvre and B. Maschke, *Pseudo-spectral methods for the spatial symplectic reduction of open systems of conservation laws*, Journal of Computational Physics 231 (2012), 1272–1292. [DOI](https://doi.org/10.1016/j.jcp.2011.10.008).
<a name="extension-ref-2"></a>

**[T2]** D. Cohen and E. Hairer, *Linear energy-preserving integrators for Poisson systems*, BIT 51 (2011), 91–101. [Author manuscript](https://www.unige.ch/~hairer/preprints/ch-poisson.pdf).
<a name="extension-ref-3"></a>

**[T3]** E. Celledoni and E. H. Høiseth, *Energy-Preserving and Passivity-Consistent Numerical Discretization of Port-Hamiltonian Systems* (2017). [arXiv:1706.08621](https://arxiv.org/abs/1706.08621).
<a name="extension-ref-4"></a>

**[T4]** J. Giesselmann, A. Karsai and T. Tscherpel, *Energy-consistent Petrov–Galerkin time discretization of port-Hamiltonian systems*, SMAI Journal of Computational Mathematics (2025). [DOI](https://doi.org/10.5802/smai-jcm.127).
<a name="extension-ref-5"></a>

**[T5]** J. L. Cieśliński, *Locally exact modifications of discrete gradient schemes* (2013). [arXiv:1304.6533](https://arxiv.org/abs/1304.6533).
<a name="extension-ref-6"></a>

**[T6]** S. Eidnes, *Order theory for discrete gradient methods*, BIT 62 (2022), 1207–1255. [DOI](https://doi.org/10.1007/s10543-022-00909-z).
<a name="extension-ref-7"></a>

**[T7]** L. Lautwein, N. Marheineke, H. N. Myhr and K. Schäfers, *Fourth-Order Energy-Consistent Average Vector Field Discrete Gradient Methods for Dissipative Hamiltonian Systems* (2026). [arXiv:2609.24692](https://arxiv.org/abs/2609.24692).
<a name="extension-ref-8"></a>

**[T8]** R. Moulla, *Space-Time Stokes-Dirac Discretization: Affine Pairing Gauges, Constitutive Cayley Coordinates, and Active Interior Reduction*, version 12 (2026). [Repository manuscript](../archive/version-12/manuscript.tex).
<a name="extension-ref-9"></a>

**[T9]** Y. Miyatake and J. C. Butcher, *A characterization of energy-preserving methods and the construction of parallel integrators for Hamiltonian systems*, SIAM Journal on Numerical Analysis 54 (2016), 1993–2013. [DOI](https://doi.org/10.1137/15M1020861), [author manuscript](https://arxiv.org/abs/1505.02537).
<a name="extension-ref-10"></a>

**[T10]** T. Kemmochi, *Higher order discrete gradient method by the discontinuous Galerkin time-stepping method* (2023). [arXiv:2308.02334](https://arxiv.org/abs/2308.02334).

For comparison with [T5], the exact-flow cap gauge of version 12 rewrites the closed quadratic step as $`\Delta x/h=J_hQ\bar x`$, where $`J_h=(J^{-1}+hS_h)^{-1}=\operatorname{tanhc}(hJQ/2)J`$. In canonical coordinates this is the specialization of [T5, equation (5.1)]. The present result concerns the placement of such gauges in the active-stage and space-time construction; it does not claim a new exact linear propagator.
