# Applied Mathematics Inventory

**What this is.** A map of the *applied numerical mathematics* used in the
`nmr_shielding` C++ library (`src/`), organised by *Numerical Recipes*
category. For each method: the operation, where the work is actually done
(`file:line`), and **who does the heavy lifting** — a library routine
(Eigen / nanoflann / libtorch), an external binary (APBS / MOPAC / GROMACS /
libdssp), or hand-rolled closed-form code.

**Threshold.** SVD-level and up. Trivial vector arithmetic (dot products,
distances, simple multiply-accumulate without special structure), table
lookups, and pure I/O parsing are excluded — noted only where the heavy math
lives *outside* our code.

**Audience.** Anyone documenting the maths (e.g. building Mathematica
worked examples against these algorithms), reviewing numerical hygiene, or
porting a kernel. Code locations are line-cited but drift; treat `file:line`
as a starting point and confirm against the function name.

---

## Two facts that shape everything below

1. **No eigendecomposition / principal-axis extraction exists anywhere.**
   The rank-2 shielding tensor is *never* diagonalised to (V_zz, η). The
   irreducible T0/T1/T2 split (`Types.cpp`) is a closed-form SO(3) projection
   in a norm-preserving real-spherical-harmonic basis, **not** a Jacobi
   eigensolve. This is deliberate — the full T2 is the thesis object
   (T2-is-sacred). If you reach for `SelfAdjointEigenSolver`, stop: it is not
   used and was not omitted by accident.

2. **No Ewald / PME anywhere in our code.** Every *local* electrostatic
   kernel is a truncated real-space cutoff sum over a k-d tree neighbour
   list. The only reciprocal-space electrostatics in the pipeline is the PME
   term GROMACS computes and writes to the `.edr`, which we read back.

**The five sources of numerical machinery.** Almost everything physics-
specific is hand-rolled closed-form on top of one of these:

| Heavy lifter | Primitive | Used by |
|---|---|---|
| **Eigen 3.4** | `JacobiSVD`, `LLT` (Cholesky) | ring plane-fit, Kabsch ×2, EEQ solve |
| **nanoflann** | k-d tree radius / k-NN | every cutoff sum (SASA, Coulomb, water, AIMNet2) |
| **libtorch** | NN forward + reverse-mode autograd | AIMNet2 charge + charge-response gradient |
| **external binaries** | APBS multigrid PB, MOPAC MOZYME SCF, GROMACS PME, libdssp | computed elsewhere, read back |
| **hand-rolled** | closed-form kernels, tensors, quadrature, interpolation, perception | the rest |

---

## 1. Dense linear-algebra factorizations

All via Eigen. Every use of SVD is on a fixed 3×3 — hence `JacobiSVD`
(one-sided Jacobi, accurate at small size) rather than `BDCSVD`
(divide-and-conquer, only wins for large matrices).

### 1.1 SVD — `Eigen::JacobiSVD`

- **Best-fit plane normal (total least squares).** Centre ring-vertex
  coordinates; the plane normal is the right singular vector of the smallest
  singular value. `Ring.cpp:31` — `JacobiSVD(coords, ComputeFullV)`,
  `matrixV().col(2)`. Sign fixed by right-hand rule. **This single SVD feeds
  the ring geometry consumed by all five ring calculators** (BiotSavart,
  HaighMallion, PiQuadrupole, RingSusceptibility, Dispersion).
- **Kabsch / orthogonal Procrustes rigid superposition.** Optimal rotation
  minimising ‖R(p−p̄) − (q−q̄)‖² from the 3×3 cross-covariance
  `H = P·Qᵀ`: `R = V·diag(1,1,sign det(V·Uᵀ))·Uᵀ`. Two sites:
  - `TripeptidePoseAssembler.cpp:46` — backbone (N, CA, C) alignment of the
    DFT/perceived pose onto the protein frame.
  - `RmsdTrackingTrajectoryResult.cpp:77` — M-point trajectory RMSD tracking.
  Both carry the **reflection guard on `sign(det)`, not raw `det`** (a raw
  det injects O(ε) anisotropic scaling and breaks orthogonality) — aligned
  across both sites in the 2026-05-21 maths review.

### 1.2 Cholesky — `Eigen::LLT`

- **EEQ charge-equilibration solve** (D4 EEQ, Caldeweyher 2019).
  `EeqResult.cpp:166`. Dense N×N symmetric-positive-definite Coulomb matrix
  (Ohno-Klopman off-diagonals `1/√(R² + 1/(η_iη_j))`, hardness + self-Coulomb
  diagonal — the self term is included precisely to make A diagonally
  dominant → SPD). The charge-neutrality constraint is handled by **bordered-
  system / KKT reduction**, not by factoring the indefinite augmented matrix:
  solve `A·u = −χ` and `A·v = 1` (two RHS reusing one factorization), then
  `λ = −(Q + 1ᵀu)/(1ᵀv)`, `q = −(u + λv)`, with a final uniform shift to
  clean the O(cond A) residual. **This is the only first-class O(N³) dense
  factorization in the library** and its dominant cost.

> Note: `JacobiSVD` and `LLT` are the *only* Eigen decomposition routines
> invoked anywhere in `src/`. There is no LU, QR, or eigensolver call.

---

## 2. Spherical-tensor / irreducible decomposition (closed-form)

The load-bearing operation of the whole project, and a closed-form algebraic
projection — **not** a diagonalization.

- **Cartesian → T0/T1/T2 irreducible spherical-tensor decomposition.**
  `Types.cpp:25-94` (`SphericalTensor::Decompose` / `Reconstruct`),
  hand-rolled, no library.
  - T0 = tr σ / 3 (rank-0 isotropic scalar).
  - T1 from the antisymmetric part via the Levi-Civita dual (rank-1
    pseudovector, 3 components).
  - T2 = the 5 real-spherical-harmonic components of the traceless-symmetric
    part, in an **isometric (Frobenius-norm-preserving)** basis: the
    `√2`, `√(3/2)`, `1/√2` coefficients are chosen so Σ|T2_m|² = ΣS_ij².
    This lets the T2 cosine / inner-product (`Types.cpp:103,114`) be a plain
    Frobenius dot on the 5-vectors.
  - Every calculator routes its `Mat3` through this before emission. O(1)
    arithmetic, but it is *the* representation choice.

---

## 3. Rank-2 tensor algebra above threshold

- **Similarity transform R·σ·Rᵀ** (frame rotation of a rank-2 tensor). Two
  rotation provenances feed the identical congruence operation:
  - **Kabsch SVD rotation** — `TripeptidePoseAssembler.cpp:72`
    (`RotateTensor`), applied to DFT shielding tensors at `:170` and `:521`.
  - **Analytic Gram–Schmidt donor frame** — `LarsenHBondGrid.cpp:438`
    (`RotateTensorToProteinLabFrame`), transpose multiply direction because
    grid tensors are *stored* as `σ_canon = R_log·σ_log·R_logᵀ`.
- **EFG outer-product construction + deviatoric (traceless) projection.**
  `V_ab = q(3 r_a r_b / r⁵ − δ_ab / r³)` accumulated over neighbours, then
  `σ −= (tr σ / 3)·I`. Each term is traceless analytically (Gauss/Laplace),
  but float accumulation breaks it, so the projection is applied before
  decomposition. Recurs across Coulomb (`CoulombResult.cpp:164,199`), AIMNet2
  (`AIMNet2Result.cpp:408,425`), APBS (`ApbsFieldResult.cpp:105,116`), and
  the ring kernels. This is the recurring **numerical-hygiene** pattern —
  symmetrize `½(M+Mᵀ)`, de-trace, then (in Coulomb) NaN/Inf-sanitize and
  magnitude-clamp runaway near-contact fields.

---

## 4. Numerical integration / quadrature

- **Haigh–Mallion dipolar surface integral — the only genuine quadrature in
  the physics suite.** `HaighMallionResult.cpp:44-148`. Evaluates
  `∫ (3 ρ_a ρ_b/ρ⁵ − δ_ab/ρ³) dS` over the ring interior, fan-triangulated
  from the centroid, by a **7-point degree-5 symmetric triangle cubature**
  (Dunavant 1985 / Stroud T2:5-1; the `(155±√15)/1200` weights and
  `(6±√15)/21` barycentric orbits computed at runtime from `√15`), with
  **depth-2 adaptive subdivision** — 4-way midpoint split triggered by
  proximity to the 1/ρ³ singularity, giving 7 / 28 / 112 quadrature points
  per fan triangle. Degenerate-triangle and per-point singularity guards.
  This adaptive refinement *is* the regularization of the near-singular
  kernel as the field point approaches the ring surface.
- **SASA — Shrake–Rupley surface sampling.** `SasaResult.cpp:19-98`. Each
  atom's expanded sphere (vdW + 1.4 Å probe) is sampled by a **Fibonacci
  golden-spiral lattice** (`θ = acos(1 − 2(i+½)/n)`, `φ = 2πi/φ_golden`) — a
  low-discrepancy quasi-Monte-Carlo sphere quadrature, lower variance than
  random sampling. SASA = 4πr²·(exposed fraction), exposure by occlusion
  test against k-d-tree neighbours. **Not** LCPO, **not** analytic
  Lee–Richards.

---

## 5. Interpolation

- **Trilinear (3-D) over regular grids.** 8-corner unit-cube blend.
  - **APBS electrostatic potential** — `ApbsFieldResult.cpp:32`
    (`GridCache::Interpolate`); `std::floor` chosen for correct negative-cell
    indexing.
  - **Larsen H-bond tensor grid** — `LarsenHBondGrid.cpp:219`
    (`TrilinearMat3`), per matrix element over axes (r, θ, ρ). The **ρ axis
    is periodic** (wraps at ±180°, `LookupAxis` periodic path,
    `LarsenHBondGrid.cpp:296`); r and θ are bounded with a ±1e-9 FP clamp.
    Grid shape `(Nr, Nθ, Nρ, 3, 3)`; lookup key `(donor class, acceptor
    class)` selects one of 6 archive grids. Cubic-spline smoothness is
    **pre-baked offline in Python** — runtime is purely trilinear.
- **Bilinear (2-D).** **CMAP** backbone φ/ψ cross-term correction —
  `BondedEnergyResult.cpp:86-114`, over a (φ,ψ) grid on [−π,π]².

---

## 6. Numerical differentiation

- **Central finite differences.** `ApbsFieldResult.cpp:68-119`. The electric
  field `E = −∇φ` is a central difference of the interpolated potential
  across one grid spacing (`:75`); the **EFG is the Hessian of φ**, a 3×3
  stencil of central differences of the interpolated field (`:81-119`), then
  explicitly symmetrized `½(EFG+EFGᵀ)` (kills interpolation-noise
  antisymmetric residue) and traceless-projected (removes the self-potential
  Laplacian-delta artifact). This is the only place we differentiate a field
  numerically; everywhere else the field gradient is an analytic kernel.

---

## 7. Automatic differentiation

- **Reverse-mode autograd (backpropagation).**
  `AIMNet2ChargeResponseGradientResult.cpp:162`. The charge-response gradient
  `dL/dr_i` of an L2-of-charges proxy objective `L = Σ_j q_j²` w.r.t. atomic
  positions, in **one backward pass**: coordinate tensor made a leaf
  (`requires_grad_(true)` after `.detach()`, `:122`), forward records the
  tape (`:139`), `loss.backward()` (`:162`), gradient read from
  `coord_gpu.grad()`. Explicitly **analytic autograd, not finite
  differences** — the exact `dq_i/dr_i` diagonal would need N backward
  passes; the L2 proxy gives a non-trivial gradient in one. The prior
  random-perturbation / finite-difference path was removed (design note at
  `AIMNet2Result.cpp:443-462`). Non-finite gradient components reject the
  whole frame to protect downstream Welford accumulators.

---

## 8. Special-function / analytic kernels

Closed-form, but above the trivial line — these carry the physics and have
real numerical structure (singularity behaviour, windowing, high inverse
powers).

- **Analytic finite-wire Biot–Savart.** `BiotSavartResult.cpp:43-63`. The
  Johnson–Bovey double-loop ring current: each straight segment's B-field is
  the closed-form finite-wire result (`cosθ₁ − cosθ₂`, written in cross/dot
  form), summed exactly over polygon edges (`:93-108`, two lobes at
  ±lobe_offset each carrying I/2). **Not** quadrature — the discretization is
  the physical polygon geometry, integrated analytically per edge. Wire-
  endpoint and on-axis singularity guards.
- **Stone rank-4 multipole T-tensor (1/r⁹).** `PiQuadrupoleResult.cpp:54-102`.
  EFG from a point axial quadrupole, `G_ab = T_abcd n_c n_d`, fully
  pre-contracted to closed form — the highest inverse-power kernel in the
  suite. Tracelessness analytic (Laplace).
- **McConnell / ring-susceptibility / H-bond dyadic kernels.** The
  `[9 cosθ d̂⊗b̂ − 3 b̂⊗b̂ − (3 d̂⊗d̂ − I)]/r³` family:
  `McConnellResult.cpp:54` (bond direction b̂; `MopacMcConnellResult`
  scales each bond by the Wiberg bond order), `RingSusceptibilityResult.cpp:54`
  (b̂ → ring normal n̂), `HBondResult.cpp:72` (b̂ → N···O H-bond axis, source
  at the N···O midpoint, DSSP-resolved pairs). Series-summed over sources
  with per-category symmetric-traceless projection.
- **C¹ switching / windowing polynomial.** `DispersionResult.cpp:62-72`. The
  London `(3 d⊗d/r⁸ − I/r⁶)` dispersion kernel multiplied by a CHARMM-form
  smooth cutoff S(r) (Brooks 1983) that is C¹-continuous at both
  R_switch and R_cut. A numerical-stability device: it makes the truncated
  1/r⁶ sum vary smoothly across MD frames rather than jump as sources cross
  the cutoff.
- **EEQ kernels.** `EeqResult.cpp:90-147`. Coordination number via a smooth
  `erfc` pairwise counting sum; Ohno–Klopman γ off-diagonals (see §1.2).

---

## 9. Geometry primitives above threshold (hand-rolled)

- **Signed dihedral via atan2 of cross products.** The ubiquitous torsion
  primitive: `n1 = b1×b2`, `n2 = b2×b3`, `m1 = n1×b̂2`,
  `θ = atan2(m1·n2, n1·n2)`. Sites: `PlanarGeometryResult.cpp:35` (ω, χ₂),
  `BondedEnergyResult.cpp:46` (FF dihedrals), `LarsenHBondGrid.cpp:380`
  (H-bond ρ), tripeptide φ/ψ/χ (`TripeptideBackboneShieldingResult.cpp:28`).
  Angle-wrap via `std::remainder(·, 2π)` (IEEE round-half-to-even).
- **Cremer–Pople 5-ring puckering.** `PlanarGeometryResult.cpp:100-165`. The
  mean-plane normal from sin/cos-weighted vertex displacements
  (`R1 = Σ r_j sin(2πj/5)`, `R2 = Σ r_j cos(2πj/5)`, `n = R1×R2`), then the
  (Q₂, θ₂) amplitude/phase of the m=2 pucker mode — a finite trigonometric-
  basis (discrete Fourier-mode) projection of out-of-plane displacement,
  phase via atan2. The adversarial-reviewed θ-inversion fix lives here (the
  naive edge-cross accumulator inverts the phase by 180°). Q < 1e-6 → NaN θ.
- **Gram–Schmidt orthonormal frame.** `LarsenHBondGrid.cpp:394-435`. The
  Larsen donor reference frame from 3 donor atoms: `z = normalize(H−anchor)`,
  `x = normalize((H−third) − ((H−third)·z) z)` (one G-S projection),
  `y = z×x`. Degeneracy guards at each step return identity + warn.
- **sp2 pyramidalization.** `PlanarGeometryResult.cpp:70-78`. Signed
  out-of-plane displacement of an atom from the plane of its 3 bonded
  neighbours (plane normal by cross product since exactly 3 points define
  it — no SVD needed), right-hand-rule sign, collinear guard.

---

## 10. Discrete / combinatorial algorithms

- **k-d tree nearest-neighbour & radius search.** `SpatialIndexResult.cpp:24`
  — nanoflann `KDTreeSingleIndexAdaptor` (leaf size 10) over atom positions,
  ring centres, bond midpoints. Radius search (`:68`) and k-NN (`:129`). The
  spatial-query workhorse underpinning every O(N·k) cutoff sum (SASA,
  Coulomb, WaterField, AIMNet2 neighbour matrix and Coulomb EFG).
- **Weisfeiler–Lehman K=3 colour-refinement graph isomorphism.**
  `LarsenResidue.cpp:565-712`. Atom-identity perception: match a perceived
  bond-graph against canonical amino-acid chemistry. 3 rounds of signature
  refinement — each node's new colour = hash of `(element, degree, sorted
  multiset of (neighbour element, neighbour prior colour))`, **FNV-1a 64-bit**
  hashing (`:543`), multiset sorted for permutation invariance. Matching
  groups nodes by final colour; **singleton class → strict identity binding**
  (`canonical_ambiguous=false`), **multi-atom class → flagged ambiguous** and
  resolved by **nearest-spatial tiebreak in the pose frame**
  (`AssembleCentralTyped:479`) — because graph-automorphic pairs (PHE
  CD1/CD2, ARG NH1/NH2, methyl Hs) provably cannot be split by any number of
  WL rounds. K=3 chosen so 3-hop neighbourhoods distinguish CG/CD/CE/CZ.
  Perception returns canonical node *indices*; identity is stamped from a
  generated topology table (FATAL on lookup miss).
  - Supporting graph ops: bond graph by distance cutoff (`BuildBondGraph`),
    peptide-amide cut + connected components by DFS (`ConnectedComponents`).
- **Grid-snap quantized lookup (not true k-NN).**
  `TripeptideDftTable::QueryNearest` (`:445`). Retrieve the DFT tripeptide
  row "nearest" to a query (φ, ψ, χ₁..χ₄): each angle is **independently
  snapped to its grid** (per-axis round-to-multiple with ±180 wrap; 2° for
  the AAA baseline, 20° otherwise) → exact relational equality SQL lookup.
  When dropped chi axes leave ties, selection is **deterministic ordinal**
  (`ORDER BY calc_id ASC LIMIT 1`), explicitly *not* a geometric chi-distance
  metric (documented hot-path tradeoff). Chi-fallback (re-query at shallower
  chi depth on perception failure) is driven by the calculators.

---

## 11. External heavy math (read back, not computed here)

Documented for honesty about where the expensive numerics actually happen.
None of this is our code; we marshal input and parse output.

- **APBS / FETK — linearized Poisson–Boltzmann.** `apbs_bridge.c:160`
  (`Vpmg_solve`). A **multigrid** elliptic-PDE solve on a 161³ grid, SDH
  boundary conditions, cubic B-spline (`VCM_BSPL2`) charge-to-grid mapping.
  We extract the potential grid and do the interpolation + finite differences
  of §5 / §6 on top of it.
- **MOPAC — semiempirical SCF.** `MopacResult.cpp:311` (subprocess). PM7 /
  MOZYME linear-scaling localized-orbital SCF, solved entirely by the
  external binary. Our side is regex text-parsing of the `.out` (Mulliken
  charges, Wiberg bond orders, dipole, heat of formation).
- **GROMACS — MD + PME.** `GromacsEnergyResult.cpp`. Per-frame energies read
  from the `.edr`, **including the `coulomb_recip` PME reciprocal-space term**
  — the only Ewald/PME electrostatics in the whole pipeline, and it is
  external.
- **libdssp + cif++ — secondary structure + SASA.** `DsspResult.cpp:112`
  (`dssp dssp_calc(structure, 3, true)`). Kabsch–Sander electrostatic
  H-bond-energy SS model (8-class) and Lee–Richards solvent-accessible
  surface (1.4 Å probe), both computed in the library. We serialise a temp
  PDB and read back φ/ψ/accessibility/donor/acceptor. **No DSSP math is
  implemented here** — it is a marshalling boundary.
- **libtorch — AIMNet2 forward inference.** `AIMNet2Result.cpp:273`
  (`module.forward()` under `NoGradGuard`). Per-atom charges + `aim`
  embedding from a TorchScript model on GPU. (The autograd *gradient* of §7
  is the one place we run `backward`.)

---

## Quick lookup table

| Method | NR category | Heavy lifter | Site |
|---|---|---|---|
| Ring plane normal | SVD / TLS | Eigen JacobiSVD | `Ring.cpp:31` |
| Kabsch superposition (×2) | SVD / orthogonal Procrustes | Eigen JacobiSVD | `TripeptidePoseAssembler.cpp:46`, `RmsdTrackingTrajectoryResult.cpp:77` |
| EEQ charge solve | Cholesky (LLT) + KKT reduction | Eigen LLT | `EeqResult.cpp:166` |
| T0/T1/T2 split | spherical-harmonic tensor decomp | hand-rolled | `Types.cpp:25` |
| R·σ·Rᵀ frame rotation | rank-2 similarity transform | hand-rolled | `TripeptidePoseAssembler.cpp:72`, `LarsenHBondGrid.cpp:438` |
| EFG build + de-trace | outer product + deviatoric projection | hand-rolled Eigen | `CoulombResult.cpp:201`, `AIMNet2Result.cpp:408` |
| Haigh–Mallion integral | triangle Gauss cubature (Dunavant 7-pt) + adaptive subdivision | hand-rolled | `HaighMallionResult.cpp:44` |
| SASA | Shrake–Rupley QMC (Fibonacci lattice) | hand-rolled + nanoflann | `SasaResult.cpp:19` |
| APBS potential sampling | trilinear interp + central-difference ∇/Hessian | hand-rolled | `ApbsFieldResult.cpp:32` |
| Larsen H-bond grid | trilinear interp (3-D, ρ periodic) | hand-rolled | `LarsenHBondGrid.cpp:219` |
| CMAP | bilinear interp (2-D) | hand-rolled | `BondedEnergyResult.cpp:86` |
| AIMNet2 charge response | reverse-mode autograd | libtorch | `AIMNet2ChargeResponseGradientResult.cpp:162` |
| Biot–Savart | analytic finite-wire + edge sum | hand-rolled | `BiotSavartResult.cpp:43` |
| Pi-quadrupole | rank-4 multipole tensor (1/r⁹) | hand-rolled | `PiQuadrupoleResult.cpp:54` |
| McConnell / ring χ / H-bond | dyadic kernel + series sum + de-trace | hand-rolled | `McConnellResult.cpp:54`, `RingSusceptibilityResult.cpp:54`, `HBondResult.cpp:72` |
| Dispersion | London 1/r⁶ + C¹ switching window | hand-rolled | `DispersionResult.cpp:62` |
| Coulomb EFG | truncated lattice sum + de-trace + clamp | hand-rolled | `CoulombResult.cpp:163` |
| Dihedral / torsion | atan2(cross, cross) | hand-rolled | `PlanarGeometryResult.cpp:35` (+ others) |
| Cremer–Pople pucker | discrete Fourier-mode projection | hand-rolled | `PlanarGeometryResult.cpp:100` |
| Gram–Schmidt frame | orthonormalization | hand-rolled | `LarsenHBondGrid.cpp:394` |
| sp2 pyramidalization | plane normal + signed projection | hand-rolled | `PlanarGeometryResult.cpp:70` |
| Spatial queries | k-d tree NN | nanoflann | `SpatialIndexResult.cpp:24` |
| Atom perception | Weisfeiler–Lehman K=3 + FNV-1a | hand-rolled | `LarsenResidue.cpp:565` |
| DFT row lookup | grid-snap quantization + ordinal tiebreak | hand-rolled + Postgres | `TripeptideDftTable.cpp:445` |
| APBS PB | multigrid elliptic PDE | external (APBS/FETK) | `apbs_bridge.c:160` |
| MOPAC | semiempirical SCF | external (MOPAC) | `MopacResult.cpp:311` |
| GROMACS PME | Ewald reciprocal-space | external (GROMACS) | `GromacsEnergyResult.cpp` |
| DSSP / SASA | Kabsch–Sander + Lee–Richards | external (libdssp) | `DsspResult.cpp:112` |
| AIMNet2 charges | NN forward inference | libtorch | `AIMNet2Result.cpp:273` |
