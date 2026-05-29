# nmr_extract API

Reader for NMR shielding tensor extraction data.  Loads NPY files
produced by the C++ extractor into typed Python objects with e3nn
Irreps, torch tensors, and numpy arrays.

## Install

```bash
pip install -e /path/to/nmr-shielding/python
```

Requires numpy, torch, e3nn.

## Load a protein

```python
from nmr_extract import load

p = load("path/to/extraction/directory")
```

Returns a `Protein` with all extracted features.  Errors if required
files are missing.  Warns on unrecognised files.

## Protein fields

```
p.protein_id            str
p.n_atoms               int
p.pos                   VectorField (N, 3)
p.element               ndarray (N,) int32 — atomic number
p.residue_type          ndarray (N,) int32
p.residue_index         ndarray (N,) int32
```

### Ring current calculators

```
p.biot_savart           BiotSavartGroup
  .shielding            ShieldingTensor (N, 9)
  .per_type_T0          PerRingTypeT0 (N, 8)
  .per_type_T2          PerRingTypeT2 (N, 40)
  .total_B              VectorField (N, 3)
  .ring_counts          RingCounts (N, 4)

p.haigh_mallion         RingKernelGroup
p.pi_quadrupole         RingKernelGroup
p.dispersion            RingKernelGroup
  .shielding            ShieldingTensor (N, 9)
  .per_type_T0          PerRingTypeT0 (N, 8)
  .per_type_T2          PerRingTypeT2 (N, 40)

p.ring_susceptibility   ShieldingTensor (N, 9)
```

### Per-ring sparse data

```
p.ring_contributions    RingContributions (P, 58)
p.ring_geometry         RingGeometry (R, 10)
```

P = number of (atom, aromatic_ring) pairs evaluated.  R = number of aromatic rings.

### Bond calculators

```
p.mcconnell             McConnellGroup
  .shielding            ShieldingTensor (N, 9)
  .category_T2          PerBondCategoryT2 (N, 25)
  .scalars              McConnellScalars (N, 6)

p.coulomb               CoulombGroup
  .shielding            ShieldingTensor (N, 9)
  .E                    VectorField (N, 3)
  .efg_backbone         EFGTensor (N, 5) — T2-only, symmetric-traceless
  .efg_aromatic         EFGTensor (N, 5) — T2-only
  .scalars              CoulombScalars (N, 4)

p.hbond                 HBondGroup
  .shielding            ShieldingTensor (N, 9)
  .scalars              HBondScalars (N, 4)  — nearest_dist, 1/r^3, count, mcconnell_scalar

p.dssp                  DsspScalars (N, 5)
p.dssp_ss8              ndarray (N, 8) | None — 8-class SS one-hot (H/G/I/E/B/T/S/C)
p.dssp_hbond_energy     ndarray (N, 4) | None — H-bond energies (acc0/acc1/don0/don1, kcal/mol)
p.dssp_chi              ndarray (N, 12) | None — chi1-4 cos/sin/exists (4 x 3 cols)
p.sasa                  ndarray (N,) | None — per-atom Shrake-Rupley SASA (A^2)
p.sasa_normal           VectorField (N, 3) | None — outward surface normal (unit vector, zero for buried)
```

### Optional groups (None if not extracted)

```
p.mopac                 MopacGroup | None
  .core.charges         ndarray (N,)
  .core.scalars         MopacScalars (N, 4)
  .core.bond_orders     BondOrders (B, 3)
  .core.global_         MopacGlobal (4,)    — [hof, dipole_x, dipole_y, dipole_z]
  .coulomb              MopacCoulombGroup — same as CoulombGroup
  .mcconnell            MopacMcConnellGroup — same as McConnellGroup

p.apbs                  APBSGroup | None
  .E                    VectorField (N, 3)
  .efg                  EFGTensor (N, 5) — T2-only, symmetric-traceless

p.orca                  OrcaGroup | None
  .total                ShieldingTensor (N, 9)
  .diamagnetic          ShieldingTensor (N, 9)
  .paramagnetic         ShieldingTensor (N, 9)

p.delta                 DeltaGroup | None
  .shielding            ShieldingTensor (N, 9)
  .scalars              DeltaScalars (N, 6)
  .apbs                 DeltaAPBS (N, 12) | None
  .ring_proximity       DeltaRingProximity (N, 6*n_removed_rings)

p.aimnet2               AIMNet2Group | None
  .charges              AIMNet2Charges (N,) — Hirshfeld charges
  .aim                  AIMNet2AimEmbedding (N, 256) — electronic structure embedding
  .efg                  EFGTensor (N, 5) — Coulomb EFG from AIMNet2 charges (T2 only)
  .efg_aromatic         EFGTensor (N, 5) — aromatic decomposition
  .efg_backbone         EFGTensor (N, 5) — backbone decomposition
  .charge_response_gradient         AIMNet2ChargeResponseGradient (N, 3) | None — d(sum q_j^2)/d(r_i)
  .charge_response_gradient_scalar  ndarray (N,) | None — L2 norm of the gradient vector

p.water_field           WaterFieldGroup | None — trajectory path only
  .efield               VectorField (N, 3) — total water E-field (V/A)
  .efield_first         VectorField (N, 3) — first-shell E-field (<3.5A)
  .efg                  EFGTensor (N, 5) — total water EFG (T2 only)
  .efg_first            EFGTensor (N, 5) — first-shell EFG (T2 only)
  .shell_counts         ndarray (N, 2) — [n_first_shell, n_second_shell]

p.hydration             HydrationGroup | None — trajectory path only
  .data                 ndarray (N, 4) — [asymmetry, dipole_cos, ion_dist, ion_charge]

p.water_polarization    WaterPolarizationGroup | None — trajectory path only
  .data                 ndarray (N, 10) — packed columns (see below)
  .dipole_vector        ndarray (N, 3) — net first-shell water dipole (e·Å)
  .surface_normal       ndarray (N, 3) — SASA-derived outward surface normal (unit vector)
  .asymmetry            ndarray (N,) — half-shell asymmetry using SASA normal
  .dipole_alignment     ndarray (N,) — cos(net dipole, SASA normal)
  .coherence            ndarray (N,) — dipole coherence |sum d_i| / n
  .shell_count          ndarray (N,) — first-shell water count

p.eeq                   EeqGroup | None
  .charges              ndarray (N,) — EEQ partial charges (elementary charges)
  .cn                   ndarray (N,) — coordination number (erfc counting)

p.gromacs_energy        ndarray (1, 43) | None — single-frame GROMACS energy, 43 cols
                        (electrostatic, bonded, VdW, thermo, box, virial,
                        pressure tensor, per-group T); for the per-frame
                        timeline use traj.energy.gromacs (load_trajectory)
p.bonded_energy         ndarray (N, 7) | None — per-atom GROMACS bonded terms
                        (bond, angle, Urey-Bradley, proper, improper, CMAP, total)
                        kJ/mol — trajectory path

p.planar_geometry       PlanarGeometryGroup | None
  .pyramidalization     ndarray (N,) — signed sp2 out-of-plane displacement (A)
  .omega_actual         ndarray (n_residues,) — per-residue omega (Ca-C-N-Ca'), radians
  .omega_deviation      ndarray (n_residues,) — omega - pi, wrapped
  .omega_is_xpro        ndarray (n_residues,) — X->Pro mask (cis/trans is real signal there)
  .aromatic_chi2        ndarray — per-aromatic-ring chi2 (ring-flip observable)
  .pucker_Q             ndarray — Cremer-Pople amplitude (saturated 5-rings, A)
  .pucker_theta         ndarray — Cremer-Pople phase (degrees)
  (omega_*_per_atom(p.residue_index) project residue arrays onto atom rows)

p.topology              TopologyGroup — sidecar tables (always present on load)
  .residues             Residues — per-residue records (names, links, flags)
  .bonds                Bonds — per-bond records (order, category, flags)
  .rings                Rings — per-ring records (kind, type, residue)
  .ring_membership      RingMembership — (ring, vertex-atom) rows
  .manifest             ExtractionManifest — axis sizes, validated at load()

p.category_info         CategoryInfo | None — per-atom AMBER/IUPAC/BMRB names
                        + typed substrate fields (CategoryInfoProjection)

p.tripeptide            TripeptideGroup | None — ProCS15/Larsen DFT (tensorcs15 DSN)
  .bb_shielding         ShieldingTensor (N, 9) — sigma_BB^i, ppm
  .bb_residual_vec      VectorField (N, 3) — central-match residual (Vec3 feature)
  .bb_match_distance    ndarray (N,) — |residual_vec| (A)
  .bb_method_tag        ndarray (N,) — 0=none, 1=OPBE Gaussian, 2=PBE ORCA
  .neighbor_shielding   ShieldingTensor (N, 9) — delta-sigma_BB^{i±1} (Larsen Eq 3)
  .neighbor_residual_vec_prev/_next   VectorField (N, 3) — per-direction cap residuals

p.larsen_hbond          LarsenHBondGroup | None — Larsen 2015 H-bond terms (ppm)
  .shielding            ShieldingTensor (N, 9) — sum of the four Table-2 classes
  .pHB_1 .pHB_2         ShieldingTensor (N, 9) — primary/secondary amide-H (HB)
  .pHaB_1 .pHaB_2       ShieldingTensor (N, 9) — primary/secondary Halpha (HaB)
  .diagnostic_CB        ShieldingTensor (N, 9) — Cbeta reality-check (should be ~0)
  .water_term           ndarray (N,) — 2.07 ppm on amide H with zero geometric H-bond candidates
  .count                ndarray (N,) — contributing H-bond pair count
```

## Tensor types

### SphericalTensor (and ShieldingTensor)

9-component packing: [T0, T1[3], T2[5]].

```python
st = p.biot_savart.shielding

st.data                 # ndarray (N, 9)
st.torch()              # torch.Tensor (N, 9)
st.irreps               # Irreps("1x0e+1x1o+1x2e")

st.T0                   # ndarray (N, 1)
st.T1                   # ndarray (N, 3)
st.T2                   # ndarray (N, 5) — m = -2, -1, 0, +1, +2
st.isotropic            # ndarray (N,) — T0 squeezed
st.T2_magnitude         # ndarray (N,) — L2 norm of T2
```

T2 component ordering matches e3nn.  No permutation needed.

### EFGTensor

A SEPARATE 5-component class (NOT a SphericalTensor) — symmetric-traceless
T2 only, Irreps `1x2e`.  Every EFG in the codebase (Coulomb, MOPAC Coulomb,
water, AIMNet2, APBS) is symmetric-traceless by construction, so T0 (trace)
and T1 (antisymmetric) are structural zeros and are not stored.

```python
efg = p.coulomb.efg_backbone

efg.data                # ndarray (N, 5) — m = -2, -1, 0, +1, +2
efg.torch()             # torch.Tensor (N, 5)
efg.irreps              # Irreps("1x2e")
```

Construction raises on a 9-column array (the pre-2026-05-18 packing);
migrate older NPYs in place with `T2 = old[..., 4:9]`.

### VectorField

```python
v = p.coulomb.E

v.data                  # ndarray (N, 3)
v.torch()               # torch.Tensor (N, 3)
v.irreps                # Irreps("1x1o")
v.x, v.y, v.z           # ndarray (N,)
v.magnitude             # ndarray (N,)
```

### PerRingTypeT2

8 ring types x 5 T2 components = 40.

```python
t2 = p.biot_savart.per_type_T2

t2.irreps               # Irreps("8x2e")
t2.for_type(RingType.PHE)   # ndarray (N, 5)
t2.as_block()           # ndarray (N, 8, 5)
t2.total                # ndarray (N, 5) — sum over types
```

### PerBondCategoryT2

5 bond categories x 5 T2 components = 25.

```python
mc = p.mcconnell.category_T2

mc.irreps               # Irreps("5x2e")
mc.for_category(BondCategory.CO_nearest)  # ndarray (N, 5)
mc.as_block()           # ndarray (N, 5, 5)
```

## RingContributions

Sparse (P, 58) table — one row per (atom, aromatic_ring) pair.

```python
rc = p.ring_contributions

rc.n_pairs              # int
rc.atom_index           # ndarray (P,) int
rc.ring_index           # ndarray (P,) int
rc.ring_type            # ndarray (P,) int — RingType enum values
rc.distance             # ndarray (P,) — Angstroms
rc.rho, rc.z, rc.theta  # cylindrical coordinates in ring frame
rc.cos_phi, rc.sin_phi  # azimuthal angle in ring plane (relative to vertex 0)

# Physics kernels — each is a SphericalTensor (P, 9)
rc.bs                   # Biot-Savart shielding kernel G
rc.hm_H                 # Haigh-Mallion raw integral H (pure T2)
rc.hm                   # Haigh-Mallion shielding kernel G = -n⊗V (rank-1)
rc.pq                   # Pi-quadrupole
rc.chi                  # Ring susceptibility
rc.disp_scalar          # ndarray (P,) — 1/r^6

# T2 from any kernel
rc.bs.T2                # ndarray (P, 5)
rc.bs.torch()           # torch.Tensor (P, 9)

# Filter
rc.for_atom(42)         # RingContributions — rows for atom 42
rc.for_ring_type(RingType.PHE)  # rows for PHE rings only
```

### Column layout (58 columns)

```
[0]     atom_index
[1]     ring_index
[2]     ring_type           0-7
[3]     distance            Angstroms
[4]     rho                 Angstroms
[5]     z                   Angstroms (signed)
[6]     theta               radians
[7]     mcconnell_factor    (3cos^2 theta - 1) / r^3
[8]     exp_decay           exp(-distance / 4.0)
[9:18]  bs_G                SphericalTensor — BS shielding kernel
[18:27] hm_H                SphericalTensor — HM raw integral (pure T2)
[27:36] hm_G                SphericalTensor — HM shielding kernel G = -n⊗V (rank-1)
[36:45] pq_G                SphericalTensor
[45:54] chi_G               SphericalTensor
[54]    disp_scalar         1/r^6
[55]    disp_contacts       vertex contact count
[56]    cos_phi             azimuthal angle cosine (relative to vertex 0)
[57]    sin_phi             azimuthal angle sine (relative to vertex 0)
```

## RingGeometry

Reference table for ring identity and spatial properties.

```python
rg = p.ring_geometry

rg.n_rings              # int
rg.ring_index           # ndarray (R,) int
rg.ring_type            # ndarray (R,) int
rg.residue_index        # ndarray (R,) int
rg.center               # ndarray (R, 3)
rg.normal               # ndarray (R, 3) — unit vectors
rg.radius               # ndarray (R,)
```

## Enums

```python
from nmr_extract import RingType, BondCategory

RingType.PHE            # 0
RingType.TYR            # 1
RingType.TRP_benzene    # 2
RingType.TRP_pyrrole    # 3
RingType.TRP_perimeter  # 4
RingType.HIS            # 5
RingType.HID            # 6
RingType.HIE            # 7

BondCategory.backbone_total   # 0
BondCategory.sidechain_total  # 1
BondCategory.aromatic_total   # 2
BondCategory.CO_nearest       # 3
BondCategory.CN_nearest       # 4
```

## WaterPolarizationGroup

Water polarisation features in the SASA-normal reference frame.
Produced by HydrationGeometryResult from explicit-solvent trajectory.

### Column layout (10 columns)

```
[0:3]   dipole_vector      net first-shell water dipole (e·Å)
[3:6]   surface_normal     SASA-derived outward surface normal (unit vector)
[6]     asymmetry          half-shell asymmetry using SASA normal (0-1)
[7]     dipole_alignment   cos(net dipole, SASA normal) (-1 to +1)
[8]     coherence          dipole coherence |sum d_i| / n_shell (e·Å)
[9]     shell_count        first-shell water count
```

## EeqGroup

Geometry-dependent partial charges from Extended Electronegativity
Equilibration (Caldeweyher et al. 2019, DOI: 10.1063/1.5090222).

```python
eeq = p.eeq
eeq.charges             # ndarray (N,) — partial charges (elementary charges)
eeq.cn                  # ndarray (N,) — coordination number (erfc counting)
```

D4 element parameters (chi, gam, kappa, rcov, rad) are in
`PhysicalConstants.h` (`D4EeqParamsFor()`), from Table S1 of
Caldeweyher et al. (2019).  All in atomic units (Hartree, Bohr).
Fixed reference data, not tuneable.

## Load a trajectory

```python
from nmr_extract import load_trajectory

traj = load_trajectory("output/trajectory.h5")                  # rollups + timelines
traj = load_trajectory("output/trajectory.h5",
                       load_optional_large=True)                # + 256-dim AIMNet2 embedding
```

Returns a `TrajectoryData` for one protein from the analysis H5 master
file (TrajectoryResults that override WriteH5Group write their own
`/trajectory/<group>/`; selections under `/trajectory/selections/<kind>/`).
The reader detects the schema (current per-TR analysis vs legacy ensemble)
and normalises positions to `(T, N, 3)`.

```
traj.protein_id         str
traj.n_atoms / n_frames int
traj.positions          ndarray (T, N, 3) — atom-major per frame
traj.frame_times        ndarray (T,) — ps
```

### Group accessors

Most TrajectoryResult families have an accessor (some H5 groups have no SDK
field yet, e.g. aimnet2_charge_time_series).  A field (or its sub-fields) is
`None` when that TR did not run for the extraction that produced the H5 —
the MOPAC groups in particular are absent unless the run was FullFat
(`--mopac`).

```
traj.welford            WelfordAccess — per-atom mean/variance over frames
                        (bs / hm / mcconnell / eeq / sasa / hbond-count)
traj.energy             EnergyAccess — per-frame GROMACS + bonded energy timeline
traj.water_field        WaterFieldAccess — water E-field/EFG (time-series + Welford)
traj.hydration_geometry HydrationGeometryAccess — SASA-normal water polarisation
traj.hydration_shell    HydrationShellAccess — COM-based water shell features
traj.dssp8              Dssp8Access — 8-state SS time-series + transitions
traj.dihedrals          DihedralTimeSeriesGroup | None — per-residue phi/psi/chi
traj.dihedral_bin_transitions  | None — rotamer / Ramachandran transitions
traj.ring_pucker        RingPuckerTimeSeriesGroup | None — per-ring Cremer-Pople timeline
traj.j_coupling         JCouplingTimeSeriesGroup | None — Karplus 3J observables
traj.ring_neighbourhood_trajectory_stats  | None — per-(atom,ring) geometric residual
traj.rmsd_tracking      RmsdTrackingGroup | None — backbone RMSD vs frame 0
traj.apbs_efg           ApbsEfgTimeSeriesGroup | None — solvated EFG TS (T2-only)

# Dynamics observables — kernel autocorrelation/spectrum + model-free layer
traj.kernel_dynamics            KernelDynamicsGroup | None — per-atom ACF + Parzen power
                                spectrum of the 13 shielding-kernel channels (the instrument)
traj.kernel_coherence           KernelCoherenceGroup | None — zero-lag kernel×kernel Pearson matrix
traj.reorientational_dynamics   ReorientationalDynamicsGroup | None — backbone P2 TCF, Henry-Szabo
                                S2, area-method tau_e, global tau_m (+ tau_m_converged flag), and
                                the 15N R1/R2/NOE + J(omega) relaxation layer on the N-H vectors
traj.ired_order_parameters      IRedOrderParameterGroup | None — reference-free iRED S2 (amide N-H)
traj.dihedral_autocorrelation   DihedralAutocorrelationGroup | None — circular ACF of phi/psi/chi
                                + 1/e decorrelation times

# AIMNet2 fleet trio
traj.aimnet2_embedding                          256-dim per-atom (load_optional_large=True)
traj.aimnet2_charge_response_gradient           Vec3 + scalar per frame
traj.aimnet2_charge_response_gradient_welford   rollup of the above

# MOPAC family — FullFat (--mopac) runs only; absent otherwise
traj.mopac_charge_welford / traj.mopac_bond_order_welford
traj.mopac_coulomb_shielding_time_series        (T2-only)
traj.mopac_mc_shielding_time_series             (T0+T1+T2)
traj.mopac_vs_ff14sb_reconciliation             signed cos(MOPAC_T2, FF14SB_T2)

# Frame selection + legacy
traj.selections         dict[kind -> list[SelectionRecordPy]] (RMSD spikes, chi rotamers, DFT poses)
traj.rollup / traj.bonds  legacy ensemble-schema rollups only (None on analysis H5)
```

The classical per-frame shielding time-series (bs / hm / mcconnell /
pi-quadrupole / ring-susceptibility / dispersion / hbond) and the
tripeptide / Larsen shielding time-series are written to
`/trajectory/*_time_series/` by their C++ TRs; consult `CATALOG` and
`nmr_extract._trajectory` for the current set of typed accessors.

## Catalog

The format contract.  Maps every NPY filename to its metadata.

```python
from nmr_extract import CATALOG

for stem, spec in CATALOG.items():
    print(f"{stem:30s} {spec.group:20s} required={spec.required}")
```

## HM: two representations

Haigh-Mallion provides both forms per ring:

- `rc.hm_H` — raw surface integral H (pure T2, geometry only)
- `rc.hm` — shielding kernel G = -n⊗V (rank-1; V = H·n), T0+T1+T2

`rc.hm` sums match `p.haigh_mallion.per_type_T2`.  `rc.hm_H` gives
the unscaled geometric kernel for analysis.
