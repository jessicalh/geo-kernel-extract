# The Dimension Tree

2026-04-13.  From the T2 shielding perturbation down to the
contributing physics, with cosines showing independence at each
branch.  Numbers from full_space_analysis.py on 400 proteins.

---

## Top level

**T2 shielding perturbation (delta):** WT minus ALA.  The rank-2
traceless symmetric tensor measuring how removing an aromatic ring
changes the anisotropic shielding at each atom.  5 components per
atom.  446,006 atoms across 720 proteins.

Branches into 3 primary dimensions + 2 secondary + controls.
The branching is not a taxonomy — it's measured.  Cross-branch
cosines show the empirical angular independence.

---

## Branch 1: Ring current magnetic field

**Physics:** Circulating pi-electrons → magnetic dipole field → T2.
**Equation:** G_ab = -n_b B_a (Pople 1956, Johnson & Bovey 1958).
**Distance:** r^-3.  Measured: -3.04 (BS), -3.06 (HM).
**Tensor:** Rank-1 outer product n ⊗ B.  T0+T1+T2.
**Independence from EFG:** |cos| = 0.44.  Near-random (0.45).

Sub-branches (within-group cosine range [0.36, 1.00]):

### 1a. Model variant: BS vs HM for same ring type
- BS_PHE ↔ HM_PHE: cos ≈ 1.000
- Same physics, same ring, different numerics.  Identical T2.
- Case 1995 found no difference for isotropic; confirmed for T2.

### 1b. Ring type: PHE vs TYR vs TRP vs HIE
- HM_PHE ↔ HM_TYR: cos ≈ 0.65.  Similar 6-membered rings.
- HM_PHE ↔ HM_TRP_pyrrole: cos ≈ 0.36.  6-ring vs 5-ring.
- HM_PHE ↔ HM_HIE: cos ≈ 0.42.  6-ring vs imidazole.
- Ring geometry (center, vertices, radius) determines the
  B-field direction.  Different rings → different T2 directions.
- Case 1995: ring current intensity factors I_PHE=1.46,
  I_TYR=1.24, I_HIE=1.35 (different magnitudes too).

### 1c. Ring susceptibility: RingSusc_total
- Point magnetic dipole at ring center (McConnell-like).
- RingSusc ↔ BS_PHE: cos ≈ 0.55.  Moderate.  The point-dipole
  approximation shares the ring-normal dependence but misses
  the wire-loop near-field structure.

---

## Branch 2: Electric field gradient

**Physics:** Point charges → Coulomb EFG tensor → T2.
**Equation:** V_ab = sum q_j K_ab(d_j) (Buckingham 1960).
**Distance:** r^-3 per charge (sum over all atoms).
**Tensor:** Symmetric, traceless (Gauss's law).  Pure T2.
**Independence from ring current:** |cos| = 0.44.

Sub-branches:

### 2a. Charge source: ff14SB vs MOPAC vs AIMNet2
- ff14SB ↔ MOPAC: cos = 0.58-0.67.  Same geometry, different
  charges.  MOPAC captures charge polarisation.
- ff14SB ↔ AIMNet2: cos = 0.34.  Nearly random.  Neural network
  charges access a different projection of charge polarisation.
- MOPAC R² = 0.39 (C) vs ff14SB = 0.22.  The polarisation
  dimension adds +0.17 for carbon.
- AIMNet2 R² = 0.07 (C).  Wrong projection for mutation delta.

### 2b. Source decomposition: backbone vs aromatic
- EFG_bb ↔ EFG_aro: cos = 0.37-0.85 (element-dependent).
- EFG_aro carries the mutation signal (aromatic charges removed).
- EFG_bb is a control (backbone charges unchanged).

### 2c. Solvation: APBS vs vacuum (DeltaAPBS)
- APBS_EFG ↔ EFG: cos ≈ 0.55.  Correction to the vacuum EFG.
- DeltaAPBS R² ≈ 0.01.  Solvation is negligible for the mutation
  perturbation (solvent boundary barely changes).

---

## Branch 3: Bond magnetic anisotropy

**Physics:** Bond susceptibility Δχ → dipolar field → T2.
**Equation:** M_ab/r^3 (McConnell 1957).
**Distance:** r^-3.
**Tensor:** Asymmetric, non-traceless.  T0+T1+T2.
**Independence from ring current:** |cos| = 0.40.  Independent.
**Independence from EFG:** |cos| = 0.39-0.42.  Independent.

Sub-branches (within-group cosine range [0.27, 0.81]):

### 3a. Bond category
- MC_backbone ↔ MC_sidechain: cos = 0.27.  Nearly orthogonal.
  Different bond populations point in different directions.
- MC_backbone ↔ MC_total: cos = 0.81.  Total dominated by
  backbone (most bonds).
- MC_aromatic carries the mutation signal (removed bonds).

### 3b. Bond-order weighting: MC vs MopacMC
- MC ↔ MopacMC: cos ≈ 0.44-0.47.  Modestly correlated.
- Bond-order weighting changes which bonds dominate the sum.
- MopacMC_total enters N forward selection at rank 4: peptide
  bond partial double-bond character (Wiberg 1.3-1.5) couples
  N electronically to the mutation site.

---

## Branch 4: Pi-quadrupole

**Physics:** Ring pi-cloud quadrupole moment → EFG → T2.
**Equation:** Stone T-tensor, leading r^-5 (Stone 2013 Ch. 3).
**Tensor:** Symmetric, traceless.  Pure T2.
**Independence from ring current:** |cos| = 0.45.  Random.
  Different multipolar orders (dipole vs quadrupole) are
  angularly orthogonal.  Empirical confirmation of Stone 2013.
**Measured distance slope:** -5.05 (theory -5).

Sub-branches:
- PQ_TRP_benzene ↔ PQ_TRP_perimeter: cos = 0.97.  Fused rings.
- PQ_PHE ↔ PQ_TRP_benzene: cos = 0.37.  Different ring sizes.
- PQ_total enters N forward selection at rank 2 (delta +0.04).
  Nitrogen atoms near rings are in the regime where dipole AND
  quadrupole contribute comparably (r^-3 vs r^-5 crossover).

---

## Branch 5: Dispersion

**Physics:** London van der Waals → proximity-gated angular signal.
**Equation:** C6 (3 d_a d_b / r^8 - delta_ab / r^6) (London 1937).
**Tensor:** Dipolar-like but r^-8 anisotropic.
**Key kernel:** DispChi = dispersion_scalar × ring_susceptibility_T2.
  This is NOT independent angular information — it's ring
  susceptibility GATED by dispersion proximity (r^-6 cutoff at 5Å).

**Independence from EFG:** |cos| = 0.53.  Moderate.
**Independence from ring current:** |cos| = 0.50.  Moderate.
  Both depend on ring geometry, producing partial correlation.

Sub-branches:
- Disp_TRP_pyrrole ↔ Disp_TRP_perimeter: cos = 0.97.  Fused.
- Disp_PHE ↔ Disp_TRP_perimeter: cos = 0.30-0.33.  Independent.

**The normalisation revelation:** raw R² = 0.06 (O), normalised
R² = 0.23.  4x increase.  Normalisation strips the per-protein
magnitude scale from the r^-6 falloff, revealing the angular
structure of dispersion for oxygen.  Dispersion becomes the
DOMINANT group for O after normalisation (0.23 vs ring current
0.19 vs EFG 0.16).

---

## Controls (correctly zero)

- **HBond_total:** R² = 0.002.  Backbone H-bonds unchanged in
  mechanical mutants.
- **DeltaAPBS_EFG:** R² = 0.005.  Solvation boundary unchanged.
- **Scalar features:** SASA, SS8, valence, bond order, dipole,
  charges — all < 0.001.  The kernels encode everything.

---

## The tree structure

```
T2 shielding perturbation (delta)
├── Ring current [cos→target: H=0.70, O=0.60]
│   ├── BS model [cos(BS,HM)=1.00 same ring]
│   │   ├── per ring type [cos range 0.36-1.00]
│   │   └── per ring distance rank (ring0..ring5)
│   ├── HM model [identical to BS in T2]
│   └── Ring susceptibility [cos(RS,BS)=0.55]
├── Electric field gradient [cos→target: H=0.93, O=0.77]
│   ├── ff14SB charges (static)
│   ├── MOPAC charges (polarised) [cos(ff14,MOP)=0.60]
│   ├── AIMNet2 charges (learned) [cos(ff14,AIM)=0.34]
│   └── backbone / aromatic decomposition
├── Bond anisotropy [cos→target: 0.40, independent]
│   ├── per bond category [cos range 0.27-0.81]
│   └── MOPAC bond-order weighted [cos(MC,MopMC)=0.45]
├── Pi-quadrupole [cos(PQ,BS)=0.45, orthogonal]
│   └── per ring type [cos range 0.37-0.97]
├── Dispersion [cos→target: H=0.85!, but sparse]
│   ├── per ring type [cos range 0.30-0.97]
│   └── DispChi = disp_scalar × RingSusc T2
└── Controls
    ├── H-bond (0.002, correctly zero)
    ├── Solvation (0.005, correctly zero)
    └── Scalars (< 0.001, encoded by kernels)
```

Each node is a physical interaction with a citation.  Each edge
has a cosine showing angular independence.  The tree is the same
for all elements — what changes is how much each branch contributes
(determined by element electronic structure, Saito 2010).
