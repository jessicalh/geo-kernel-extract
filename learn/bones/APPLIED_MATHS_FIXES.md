# Applied Maths Fixes: Numerical Stability for Calculators

Running pad. Each calculator gets a section BEFORE implementation.
Review the old project's code, critiques (bones/), Gemini reviews,
and session logs for numerical tricks that solved real problems.

This is not optional. Structurally correct physics with numerically
naive implementation produces garbage at critical atom positions.

---

## General (all calculators)

- MIN_DISTANCE = 0.1 Å singularity guard on all 1/r^n terms (prevents
  numerical 1/0). This is a numerical guard, not a physics filter.
- **KernelEvaluationFilter framework** (KernelEvaluationFilter.h):
  every calculator holds a KernelFilterSet with physics-based filters.
  Four concrete filters:
  - DipolarNearFieldFilter: multipole expansion invalid inside source
    distribution (distance > source_extent/2)
  - SelfSourceFilter: field undefined at source atom itself
  - SequentialExclusionFilter: through-bond coupling not modelled
  - RingBondedExclusionFilter: topological exclusion of ring vertices
    and atoms bonded to ring vertices (through-bond regime). Added
    2026-04-02 during calculator audit. Uses bond graph, not distance.
    Fixed HM fused ring artifact (1.127→1.000) and PQ max|T2|
    reduced 11x (7.39→0.66 A^-5).
  Per-filter rejection counts logged.
- SVD for ring normals (already in Ring::ComputeGeometry — stable
  even for non-planar rings)
- Small PDB differences (0.001 Å) can destroy angular signal:
  cos(theta) computed from dot products on nearly-parallel vectors
  loses precision. Consider atan2 or cross-product magnitude where
  angle stability matters.

## McConnell Bond Anisotropy

- KernelFilterSet: SelfSourceFilter (atom is bond endpoint) +
  DipolarNearFieldFilter (source_extent = bond length ~1.2-1.5A).
  Self-exclusion was previously done via atom_bond_membership lookup;
  now unified into the filter framework.
- Bond length as source_extent: the nearest non-bonded atom to a bond
  midpoint is typically ~1.5A. DipolarNearFieldFilter threshold is
  bond_length/2 ≈ 0.6-0.75A. This rarely triggers beyond what
  SelfSourceFilter already catches, but is honest.
- Category T2 totals: extract symmetric part of accumulated M, apply
  traceless projection, then decompose. Floating-point accumulation
  across many bonds can break tracelessness of the symmetric part.

## Biot-Savart Ring Current

- Wire segment field: near-field divergence when atom is close to
  segment endpoint. Old code checks lenA < 1e-25 and crossSq < 1e-70.
  These thresholds are in SI (metres), so they correspond to
  effectively-zero distances in Angstroms.
- Johnson-Bovey: two loops at ±d. When atom is in the ring plane
  (z ≈ 0), contributions from upper and lower loops nearly cancel
  for the z-component of B. Precision loss. The old code handles
  this by working in SI (metres) to avoid Angstrom-scale underflow.
- **Sign convention**: G_ab = -n_b * B_a * PPM_FACTOR. The minus
  sign comes from sigma_ab = -dB_a^sec / dB_{0,b}. With this
  convention, sigma = I * G gives the correct physical sign using
  literature ring current intensities (I < 0 for diamagnetic).
  Verified analytically: I=-12, atom 3A above PHE -> sigma = +1.40
  ppm (shielded, matching Case 1995). In-plane at 5A -> sigma =
  -0.16 ppm (deshielded). Unit chain: A -> m -> SI Tesla -> x PPM
  -> x literature I -> ppm.
- KernelFilterSet: DipolarNearFieldFilter (source_extent = ring
  diameter). 9 rejections on 1UBQ (4 rings, 602 atoms).

## Haigh-Mallion Surface Integral

- Adaptive subdivision: old code subdivides triangles when atom is
  within 2.0 Å (level 1) or 1.0 Å (level 2) of the ring.
  2 levels max recursion. Each level splits each triangle into 4
  at edge midpoints. Preserved in new code (HM_SUBDIV_THRESHOLD_L1,
  HM_SUBDIV_THRESHOLD_L2).
- 7-point Gaussian quadrature (Stroud T2:5-1, Dunavant degree-5).
  Weights: centroid 9/40, orbit 1 (near vertices) ≈ 0.1259,
  orbit 2 (near edges) ≈ 0.1324.
- **TWO tensors stored**: (1) H_ab = ∫K dS, the raw surface integral
  (symmetric, traceless, A^-1). Tracelessness verified at 4.3e-14
  across 906K pairs. (2) G_ab = -n_b(H·n)_a, the full shielding
  kernel (rank-1, same sign convention as BS). Boyd & Skrynnikov
  JACS 2002 construction.
- **Sign convention**: Same as BS. sigma = I * G. The minus sign
  in G = -n⊗(H·n) parallels G = -n⊗B in BS.
- KernelFilterSet: DipolarNearFieldFilter (source_extent = ring
  diameter). Same rejection profile as BS.

## Coulomb EFG

### Unit chain (explicit, from first principles)

The Coulomb E-field at atom i from all other charges j:

```
E_a(i) = sum_{j!=i} q_j * (r_i - r_j)_a / |r_i - r_j|^3
```

With q in elementary charges (e) and r in Angstroms (A), the raw sum
has units e/A^2. To get the physical electric field in V/A, multiply
by Coulomb's constant in {e, A, eV} units:

```
ke = e / (4 pi epsilon_0) = 14.3996 eV*A / e
   = 14.3996 V*A                              (since eV/e = V)
```

So: E(V/A) = 14.3996 * E_raw(e/A^2).

For the EFG tensor V_ab = dE_a/dx_b:
V(V/A^2) = 14.3996 * V_raw(e/A^3).

CoulombResult stores E in V/A and EFG in V/A^2.

### APBS comparison

APBS returns potential in kT/e. E-field from grid differences is in
kT/(e*A). To convert to V/A:

```
kT/e at 298.15K = k_B * T / e = 0.025693 V
E(V/A) = 0.025693 * E_apbs(kT/(e*A))
```

Both CoulombResult and ApbsFieldResult are in V/A after conversion.
The solvent contribution is: E_solvent = E_apbs - E_coulomb_total.

### Buckingham coefficients

With E in V/A:
- sigma_iso = A * E_z + B * E_z^2    (ppm)
- sigma_T2 = gamma * EFG_T2          (ppm)

A in ppm/(V/A), B in ppm/(V/A)^2, gamma in ppm/(V/A^2).
These have direct literature meaning: A_H from Case 1995 is the
hydrogen shielding polarisability. Learnable, but interpretable.

### Numerical stability

- Tracelessness: sum of traceless terms should be traceless. In
  practice, floating point accumulation breaks this. Apply traceless
  projection after summation: V -= (Tr(V)/3) * I.
- Already done in ApbsFieldResult for APBS EFG. Must do same for
  vacuum Coulomb.
- MIN_DISTANCE = 0.1 A singularity guard on 1/r^3 and 1/r^5 terms.
- Sanitise NaN/Inf from near-zero denominators.
- Clamp E magnitude to APBS_SANITY_LIMIT (100 V/A) for rare pathological geometries.

## Pi-Quadrupole

- **Corrected EFG tensor from Stone T-tensor formalism.**
  The catalogue's original formula had non-zero trace (violated Laplace).
  Correct kernel (EFG from point axial quadrupole at ring center):
  ```
  G_ab = 105 dn² d_a d_b / r⁹
       - 30 dn (n_a d_b + n_b d_a) / r⁷
       - 15 d_a d_b / r⁷
       + 6 n_a n_b / r⁵
       + δ_ab (3/r⁵ - 15 dn²/r⁷)
  ```
  Derivation: V_ab = -(1/3)Θ_cd T_abcd where T_abcd = ∂⁴(1/r)/∂x_a∂x_b∂x_c∂x_d.
  Stone, Theory of Intermolecular Forces, OUP 2013, Ch. 3.
- Properties verified on 723 proteins, 1.03M atom-ring pairs:
  max |Tr| = 5.1e-15 (machine precision traceless),
  max asymmetry = 2.2e-16 (machine precision symmetric).
- Pure T2 (T0 = 0, T1 = 0). The isotropic shift comes from the
  E-field scalar (3cos²θ-1)/r⁴ via Buckingham A-term, stored separately.
- DipolarNearFieldFilter with source_extent = ring diameter. 9 rejections
  on 1UBQ (same as other ring calculators).
- MIN_DISTANCE = 0.1 A singularity guard. Terms go as 1/r⁹ — steeper
  than dipolar (1/r⁵ at leading order). Near-field stability is critical.
- Five-membered rings (TRP5, HIE) produce larger signals than
  six-membered (PHE, TYR, TRP6) because atoms can approach the smaller
  ring center more closely. TRP9 (fused perimeter) produces the smallest
  signal per pair — larger ring averages down the quadrupole EFG.

## Ring Susceptibility

- Same full tensor formula as McConnell with b_hat → ring normal.
  Same numerical properties: T0 = f identity holds at machine precision.
- DipolarNearFieldFilter with source_extent = ring diameter (2 * radius).
  Ring atoms sit at the ring boundary (~1.4A from center for benzene),
  so the nearest non-ring atom is typically >2A from center. The filter
  rarely triggers but is honest about model validity.
- Ring normal from SVD (already in Ring::ComputeGeometry). Stable even
  for non-planar rings (TRP pyrrole slightly puckered).
- Per-ring-type validation: six-membered rings (PHE, TYR, TRP6) should
  produce consistent mean |T0|. Five-membered (TRP5, HIE) should be
  ~30% higher (smaller ring → atoms closer to center). TRP9 (fused
  perimeter) should be ~2x higher. Verified on 720 proteins.

## London Dispersion

- Distance filter R_MIN = 1.5 Å, R_SWITCH = 4.3 Å, R_CUT = 5.0 Å
  per vertex. **CHARMM smooth switching function** (Brooks et al. 1983):
  S(r) = 1 for r <= R_SWITCH, smooth C1-continuous taper for
  R_SWITCH < r < R_CUT, S(r) = 0 for r >= R_CUT.
  Formula: S(r) = (Rc²-r²)²(Rc²+2r²-3Rs²) / (Rc²-Rs²)³.
  Prevents feature discontinuities across MD ensemble frames.
- DipolarNearFieldFilter at ring level (source_extent = ring diameter)
  plus through-bond vertex exclusion (atoms bonded to ring vertices
  excluded, inline check in Compute()).
- Kernel per vertex: K_ab = S(r) * (3 d_a d_b / r⁸ - δ_ab / r⁶)
  (traceless per vertex). Sum over ring vertices. Unit C6=1 —
  learnable per ring type.
- Scalar per vertex: S(r) / r⁶. Sum is the isotropic dispersion signal.
- Small effect: max |T2| = 0.46 A⁻⁶ across 723 proteins. Dominated
  by ring current and bond anisotropy T2 (which go as 1/r³).
- Moderately correlated with BS (0.71) and PQ (0.66) — expected,
  all three come from ring vertex geometry. Near-random vs McConnell
  (0.41) and Coulomb (0.47).

## H-Bond Dipolar

- Same full tensor formula as McConnell with b_hat → H-bond direction
  (donor N → acceptor O).
- **Critical: DipolarNearFieldFilter required.** The H-bond N...O
  distance (~2.8A) means the midpoint is only 1.4A from each endpoint.
  Atoms bonded to the donor N or acceptor O (the donor H at ~0.5A from
  midpoint, the acceptor C at ~0.2A) are inside the source distribution.
  Without filter: max |T2| = 1908 A^-3. With filter: max |T2| = 0.78.
  The filter threshold is source_extent/2 = N...O distance / 2.
- SelfSourceFilter excludes the donor N and acceptor O atoms themselves.
- H-bonds identified from DSSP (Kabsch-Sander energy criterion).
  Backbone only (N-H...O=C). Resolved to atom positions via residue
  backbone index cache.
- Deduplication by (donor_N, acceptor_O) pair. DSSP reports each H-bond
  from both sides (as acceptor of residue i and donor of residue j);
  std::set prevents double-counting.
- Sequential exclusion: sequence_separation < 2 excluded during H-bond
  identification (not in filter set, because it's a property of the
  H-bond, not of the kernel evaluation).
