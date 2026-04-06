# Extraction Order: NMR Shielding Tensor Prediction System

This document defines the ordering of extraction via the dependency graph
of ConformationResult types. Each result declares its dependencies; the
pipeline resolves the order. The manual phase numbering from the prior
version is replaced by the dependency graph. Every operation from the
v1 categorical analysis (ops 1-470) is accounted for.

Aligned with CONSTITUTION.md (2026-03-29 revision). Key change: manual
Phases 0-6 replaced by dependency-graph-driven extraction via
ConformationResult singletons.

---

## Sign Convention (stated once, used everywhere)

NMR convention: shielding is positive. sigma > 0 = diamagnetic shielding.
Ring current above an aromatic ring produces positive sigma (more shielded).
In the ring plane, sigma is negative (deshielded).

Biot-Savart geometric kernel: G_ab = -n_b * B_a * PPM_FACTOR (negative sign,
matching NMR literature).

Verification test: a proton 3 Angstroms above a PHE ring on the normal axis
has sigma > 0.

---

## Dependency Graph

The pipeline resolves extraction order from the dependency graph:

1. Check which results are attached to the ProteinConformation
2. Find results whose dependencies are all satisfied
3. Compute and attach them
4. Repeat until all desired results are attached

The order emerges from the science (can't run Coulomb before charges)
not from a numbered list. An agent adding a new result type declares
dependencies and the pipeline slots it in.

### Dependency graph (canonical)

```
GeometryResult               requires: nothing
DsspResult                   requires: nothing
ChargeAssignmentResult       requires: nothing
EnrichmentResult             requires: nothing
OrcaShieldingResult          requires: nothing
    |
    v
SpatialIndexResult           requires: GeometryResult
ApbsFieldResult              requires: ChargeAssignmentResult
    |
    v
MolecularGraphResult         requires: SpatialIndexResult
BiotSavartResult             requires: SpatialIndexResult, GeometryResult
HaighMallionResult           requires: SpatialIndexResult, GeometryResult
McConnellResult              requires: SpatialIndexResult, GeometryResult
PiQuadrupoleResult           requires: SpatialIndexResult, GeometryResult
RingSusceptibilityResult     requires: SpatialIndexResult, GeometryResult
DispersionResult             requires: SpatialIndexResult, GeometryResult
CoulombResult                requires: ChargeAssignmentResult, SpatialIndexResult
HBondResult                  requires: DsspResult, SpatialIndexResult
    |
    v
FeatureExtractionResult      requires: all physics results above
    |
    v
PredictionResult             requires: FeatureExtractionResult
```

### Resolved order (one valid topological sort)

The following is one valid execution order. Other topological sorts
are equally valid -- the pipeline only guarantees that dependencies
are satisfied before a result runs.

**Tier 0 (no dependencies -- can run in parallel):**
- GeometryResult
- DsspResult
- ChargeAssignmentResult
- EnrichmentResult
- OrcaShieldingResult

**Tier 1 (depends on Tier 0 results):**
- SpatialIndexResult (requires GeometryResult)
- ApbsFieldResult (requires ChargeAssignmentResult)

**Tier 2 (depends on Tier 0-1 results -- can run in parallel):**
- MolecularGraphResult (requires SpatialIndexResult)
- BiotSavartResult (requires SpatialIndexResult, GeometryResult)
- HaighMallionResult (requires SpatialIndexResult, GeometryResult)
- McConnellResult (requires SpatialIndexResult, GeometryResult)
- PiQuadrupoleResult (requires SpatialIndexResult, GeometryResult)
- RingSusceptibilityResult (requires SpatialIndexResult, GeometryResult)
- DispersionResult (requires SpatialIndexResult, GeometryResult)
- CoulombResult (requires ChargeAssignmentResult, SpatialIndexResult)
- HBondResult (requires DsspResult, SpatialIndexResult)

**Tier 3 (depends on all physics results):**
- FeatureExtractionResult

**Tier 4 (depends on features):**
- PredictionResult

---

## Construction: Load and Build Objects

Before any ConformationResult can be attached, the protein and its
first ProteinConformation must exist.

### What runs
- Parse PDB/mmCIF via cifpp
- Construct Protein (sequence, residue types, canonical atom templates)
- Construct first ProteinConformation via factory method:
  - Crystal structure: `protein.AddCrystalConformation(positions, metadata)`
  - NMR ensemble: `protein.AddNMRConformation(positions, metadata)` per model
  - MD trajectory: `protein.AddMDFrame(positions, metadata)` per frame
  - Prediction: `protein.AddPrediction(positions, metadata)`
- Detect bonds via OpenBabel bond perception
- Detect rings from residue type + atom presence + protonation state
- Assign protonation state from hydrogen atoms present
- Construct ProteinBuildContext from source metadata

### Prerequisites
- PDB/mmCIF file
- Source metadata (PDB ID, deposition date, organism, etc.)
- Experimental metadata for conformation subtype (resolution, R-factor,
  temperature for crystal; ensemble member for NMR; etc.)

### What it produces

**On Protein:**
- Residue list with types, sequence numbers, chain IDs, insertion codes
  (type: AminoAcid enum, per residue)
- Atom list with elements, PDB atom names
  (type: Element enum, std::string for display only)
- Ring definitions (topology only, no geometry): vertex atom indices,
  type class, parent residue, fused partner index
  (type: per ring type class, vector<size_t> atom indices)
- Bond definitions (topology only): atom index pairs, bond order,
  bond category
  (type: BondOrder enum, BondCategory enum)
- Chain addressing: chain ID per residue, multi-chain separation markers
  (type: std::string chain_id, int sequence_number, std::string insertion_code)
- ProteinBuildContext: source, protonation tool, force field, assumptions

**On ProteinConformation (typed subclass):**
- Atom positions: Vec3 per atom, in Angstroms (const after construction)
- Protonation state: per-residue protonation decisions
- Subtype-specific metadata (resolution, R-factor, etc.)

### Categorical analysis ops accounted for
- ops 1-7 (backbone atom identity): will be resolved as typed Role at enrichment,
  but the PDB atom name is stored here as metadata
- op 21 (ring type index): ring type class assigned from residue name + atom
  presence + protonation state
- op 272 (HIS normalization for residue encoding): protonation variant preserved
  through ring type assignment (Bug 1 fix: HID/HIE/HIP NOT normalized to HIS)

### Insertion code handling
Insertion codes (PDB column 27, e.g. "12A", "12B") are stored on each residue.
Sequence separation computations use (chain_id, sequence_number, insertion_code)
triples, not bare sequence numbers. Two residues with the same sequence number
but different insertion codes are distinct positions.

### Multi-chain correctness
Chain ID is an addressing property on each residue, not an entity.
Sequence separation between atoms on different chains is infinity
(represented as INT_MAX). No bond, no BFS path, no sequence-based
feature crosses a chain boundary unless the bond graph explicitly
connects them (e.g., disulfide).

---

## GeometryResult

**Requires:** nothing
**Accessor:** `conformation.GeometryData()`

### What it computes
Ring geometry and bond geometry from const positions. Pre-built
type/category collections. Ring pair properties.

### What it stores

**Ring geometry (per ring):**
- center: Vec3, Angstroms (mean of vertex positions)
- normal: Vec3, normalised (from SVD of vertex positions)
- radius: double, Angstroms (mean vertex distance from center)
- vertices: vector<Vec3>, Angstroms (vertex atom positions)

**Bond geometry (per bond):**
- length: double, Angstroms
- direction: Vec3, normalised
- midpoint: Vec3, Angstroms

**Global geometry:**
- bounding_min, bounding_max: Vec3, Angstroms
- center_of_geometry: Vec3, Angstroms
- radius_of_gyration: double, Angstroms

**Pre-built collections:**
- rings_by_type: map<RingTypeIndex, vector<size_t>>
- bonds_by_category: map<BondCategory, vector<size_t>>
- residues_by_type: map<AminoAcid, vector<size_t>>

**Ring pair properties:**
- ring_pairs: for each pair of rings:
  - center_distance: double, Angstroms
  - normal_dot: double (n1 . n2)
  - normal_cross_mag: double (|n1 x n2|)
  - is_fused: bool
  (type: vector<RingPair>, O(R^2) where R ~ 20)

Replaces ops 148-149 (ring normal dot/cross products).

---

## DsspResult

**Requires:** nothing
**Accessor:** `conformation.Dssp()`

### What it computes
Secondary structure assignment, backbone angles, SASA, H-bond
partners via libdssp.

### Input
ProteinConformation atom positions, Protein residue/atom structure.

### What it stores (per residue, via DsspPerResidue)
- secondary_structure: char (H/G/I/E/B/T/S/C)
- phi: double, radians
- psi: double, radians
- sasa: double, Angstroms^2
- hbond_acceptors[2]: residue index + energy
- hbond_donors[2]: residue index + energy

### Physics query methods
- `SecondaryStructure(residue_index) -> char`
- `Phi(residue_index) -> double`
- `Psi(residue_index) -> double`
- `SASA(residue_index) -> double`
- `HBondPartners(residue_index) -> vector<HBondPartner>`

---

## ChargeAssignmentResult

**Requires:** nothing
**Accessor:** `conformation.ChargeAssignment()`

### What it computes
Partial charges and VdW radii from force field (tleap/prmtop, ff14SB).

### Input
Protonation state, force field from ProteinBuildContext.

### What it stores (per atom)
- partial_charge: double, elementary charge units (e)
- vdw_radius: double, Angstroms

### Physics query methods
- `ChargeAt(atom_index) -> double`
- `VdWRadiusAt(atom_index) -> double`
- `TotalCharge() -> double`

---

## EnrichmentResult

**Requires:** nothing
**Accessor:** `conformation.Enrichment()`

### What it computes

**Hybridisation (via OpenBabel):**
- Input: Bond topology, atom positions
- Output per atom: hybridisation (Hybridisation enum: sp, sp2, sp3, unassigned)

**Atom role assignment (typed, no string parsing):**
- Input: Element, bond connectivity, residue type, backbone atom indices
- Output per atom: role (AtomRole enum)

The role assignment replaces ALL string-based identity checks from v1:
- ops 1-7 (is_backbone_atom): role in {BackboneN, BackboneCA, BackboneC,
  BackboneO, AmideH, AlphaH}
- ops 8-9 (check_amide_H): role == AmideH
- ops 10-12 (check_alpha_H): role == AlphaH
- ops 13-18 (check_methyl): role == MethylH
- ops 30-33 (is_sp2_parent): hybridisation == sp2 (from OpenBabel)
- op 47 (is_hydrogen): element == Element::H
- op 52 (is_aromatic_residue): residue.IsAromatic()
- op 324 (is_hbond_donor): from role + element + bond connectivity
- op 325 (is_hbond_acceptor): from role + element

**Categorical atom properties (all booleans, set once):**
- Input: Role, element, residue type, bond connectivity
- Output per atom:
  - is_backbone: bool (role in backbone set)
  - is_amide_H: bool (role == AmideH)
  - is_alpha_H: bool (role == AlphaH)
  - is_methyl: bool (role == MethylH)
  - is_aromatic_H: bool (role == AromaticH)
  - is_on_aromatic_residue: bool (residue.IsAromatic())
  - is_hbond_donor: bool
  - is_hbond_acceptor: bool
  - parent_is_sp2: bool (for hydrogen atoms: parent hybridisation == sp2)
  - parent_atom_index: size_t (for hydrogen atoms: nearest bonded heavy atom)

Replaces ops 47-52, 60, 326-327.

**Pre-built collection:**
- atoms_by_role: map<AtomRole, vector<size_t>>

### Categorical analysis ops accounted for
- ops 1-7: replaced by Role taxonomy
- ops 8-12: replaced by Role taxonomy
- ops 13-18: replaced by Role taxonomy
- ops 19-20 (ring_has_nitrogen): replaced by ring type class NitrogenCount()
- ops 30-33 (is_sp2_parent): replaced by hybridisation from OpenBabel
- ops 47-52: replaced by element enum + role + categorical bools
- ops 53-60: replaced by bond connectivity + parent_atom_index
- op 364 (backbone_plus_cb set): replaced by Role taxonomy -- "sidechain beyond CB"
  is atoms whose role is SidechainC/N/O/S or AromaticC/N/AromaticH

---

## SpatialIndexResult

**Requires:** GeometryResult
**Accessor:** `conformation.SpatialIndex()`

### What it computes

**Spatial index construction (nanoflann):**
- Input: ProteinConformation atom positions, ring centers, bond midpoints
- Output:
  - atom_kd_tree: nanoflann KD-tree over atom positions
  - ring_center_kd_tree: nanoflann KD-tree over ring centers
  - bond_midpoint_kd_tree: nanoflann KD-tree over bond midpoints

**Neighbour list construction:**
- Input: Spatial index, all atom positions
- Output per atom:
  - spatial_neighbours: for each atom, all atoms within 15A with stored
    distance (double, Angstroms) and direction (Vec3, normalised)
    (type: vector<AtomNeighbour>)

Replaces the O(N) brute-force scans in ops 54-57 (parent search),
ops 61-65 (acceptor search), ops 66-70 (C=O search), ops 71-73 (packing density),
ops 74-77 (ring neighbour minimum vertex distance).

### Categorical analysis ops accounted for
- ops 36-46: replaced by spatial query + ring neighbourhood
- ops 61-65: replaced by spatial query for H-bond acceptors
- ops 66-70: replaced by spatial query for nearest C=O
- ops 71-73: replaced by spatial query for heavy atom count
- ops 74-77: replaced by ring neighbour minimum vertex distance

### Enrichment summary
After GeometryResult, EnrichmentResult, and SpatialIndexResult are attached,
the ProteinConformation has:
- Every atom enriched with role, categoricals, hybridisation, partial charge,
  spatial neighbourhood
- Every ring with computed geometry and pair relationships
- Pre-built type/role/category collections for O(1) filtered queries
- KD-trees for all spatial queries

---

## ApbsFieldResult

**Requires:** ChargeAssignmentResult
**Accessor:** `conformation.ApbsField()`

### What it computes
Solvated electric field and field gradient via APBS (Poisson-Boltzmann).

### Input
Atom positions, partial charges, VdW radii from ChargeAssignmentResult,
experimental conditions from ProteinConformation metadata.

### What it stores (per atom)
- apbs_efield: Vec3, solvated E-field at atom position (V/Angstrom)
- apbs_efg: Mat3, solvated electric field gradient tensor (V/Angstrom^2)
- apbs_efg_spherical: SphericalTensor decomposition of apbs_efg
  (T0: double, T1: array<3>, T2: array<5>)

APBS values are sanitised: NaN/Inf clamped to zero, magnitude > 100 V/A
flagged as suspect (op 314, op 468).

Replaces ops 313-323 (apbs_field_features input).

### Physics query methods
- `EFieldAt(atom_index) -> Vec3`
- `EFGAt(atom_index) -> Mat3`
- `EFGSphericalAt(atom_index) -> SphericalTensor`
- `IsSanitised(atom_index) -> bool`

---

## MolecularGraphResult

**Requires:** SpatialIndexResult
**Accessor:** `conformation.MolecularGraph()`

### What it computes

**Molecular graph construction and BFS:**
- Input: Bond topology, atom positions, ring atom indices
- Output per atom:
  - graph_dist_ring: int, BFS bond count to nearest ring atom
    (-1 if unreachable, sentinel 99.0 at feature time)
  - graph_dist_N: int, BFS bond count to nearest nitrogen
  - graph_dist_O: int, BFS bond count to nearest oxygen
  - eneg_sum_1: double, sum of Pauling electronegativities of 1-bond neighbours
  - eneg_sum_2: double, sum of electronegativities of 2-bond neighbours
  - n_pi_bonds_3: int, number of pi bonds within 3 bonds
  - is_conjugated: bool, part of a conjugated bond chain

**Through-bond BFS distance:**
- Input: Bond topology, ring atom indices
- Output per atom:
  - bfs_to_nearest_ring_atom: int, minimum bond count to any aromatic ring atom
  - bfs_decay: double, e^{-lambda * n} where lambda is a learned parameter
    and n is the BFS distance. Lambda initialized to ln(2)/3 (half at 3 bonds).

This is the through-bond proxy feature. It provides through-bond distance
information that through-space distance alone cannot capture (e.g., cis vs
trans isomers have different through-space distance to the same ring but
identical through-bond distance).

Replaces ops 290-302 (graph_features) with richer, correctly computed BFS.
Fixes Bug 2 (MolecularGraph BFS differences) by using the typed bond graph
from construction rather than re-detecting bonds with different thresholds.

### Physics query methods
- `BondDistToRing(atom_index) -> int`
- `BondDistToElement(atom_index, Element) -> int`
- `BFSDecay(atom_index) -> double`
- `IsConjugated(atom_index) -> bool`

### Categorical analysis ops accounted for
- ops 290-302 (graph_features): replaced by MolecularGraphResult
- ops 434-441 (optional graph/APBS guards): replaced by framework's
  dependency checking -- if a result isn't attached, the accessor throws

---

## BiotSavartResult

**Requires:** SpatialIndexResult, GeometryResult
**Accessor:** `conformation.BiotSavart()`

### What it computes
For each ring in the protein:
  For each atom within range (15A cutoff from ring center, via spatial index):

**Biot-Savart ring current:**
- Input: Ring geometry (center, normal, vertices, JB lobe offset), atom position
- Computation: Johnson-Bovey double-loop B-field model
- Output per atom per ring (RingNeighbourhood entry):
  - G_tensor: Mat3, full geometric kernel G_ab (dimensionless)
  - G_spherical: SphericalTensor of G_tensor
  - B_field: Vec3, magnetic field at atom from this ring (dimensionless)
  - B_cylindrical: Vec3 (B_n, B_rho, B_phi) in ring frame
  - gaussian_density: double, learned per-type spatial envelope

### Ring property post-pass (after all atoms processed for one ring)
- R.total_B_field_at_center: Vec3, sum of B-field from all OTHER rings
  at R's center position
- R.intensity_used: double, the ring type's effective intensity
- R.total_G_T0_diagnostic: double, sum of |G_T0| over all nearby atoms
- Ring-to-ring mutual field: for each ring pair (R, R'), the B-field
  from R at R'.center and from R' at R.center

### Atom post-pass (after all rings processed for one atom)
- ring_neighbourhood sorted by distance to ring center
- n_rings_within_3A, 5A, 8A, 12A: int counts
- mean_ring_distance: double (average distance to all ring centers)
- nearest_ring_atom_distance: double (minimum distance to any vertex of
  nearest ring)
- total_B_field: Vec3, sum of B_field from all rings
- total_G_tensor: Mat3, sum of G_tensor from all rings
- total_G_spherical: SphericalTensor of total_G_tensor
- per_type_G_T0_sum[8]: double per ring type, accumulated G_T0
- per_type_G_T2_sum[8]: array<5> per ring type, accumulated G_T2
- G_iso_exp_sum: double, exp-weighted sum of G_T0 (decay length 4.0A)
- G_T2_exp_sum: array<5>, exp-weighted sum of G_T2
- G_iso_var_8A: double, variance of G_T0 within 8A

### Physics query methods
- `TensorAtAtom(atom_index, ring_index) -> Mat3`
- `BFieldAtAtom(atom_index, ring_index) -> Vec3`
- `SumT0ByRingType(atom_index, RingTypeIndex) -> double`
- `TotalBField(atom_index) -> Vec3`
- `NearestRingDistance(atom_index) -> double`
- `RingNeighbourhoodAt(atom_index) -> const vector<RingNeighbourhood>&`

### Categorical analysis ops accounted for
- ops 78-93: Biot-Savart per-ring accumulation
- ops 94-104: Type-summed isotropic factors, B-field magnitude, ring counts,
  ring distances
- ops 105-132: Nearest ring features (ring-frame scalars, cylindrical B,
  one-hot type)
- ops 133-146: Second nearest ring features
- ops 147-156: Multi-ring geometry (normal dot/cross, variance, ratio)
- ops 157-175: L1 B-fields, normals, displacements
- ops 176-183: L2 type-summed T2, nearest ring T2

Replaces ops 78-93, 105-132, 133-146, 157-163, 164-175, 176-183.

---

## HaighMallionResult

**Requires:** SpatialIndexResult, GeometryResult
**Accessor:** `conformation.HaighMallion()`

### What it computes
For each ring in the protein:
  For each atom within range:

**Haigh-Mallion surface integral:**
- Input: Ring geometry, atom position
- Computation: 7-point Gaussian quadrature on fan triangulation
- Output per atom per ring (RingNeighbourhood entry):
  - hm_tensor: Mat3, full rank-2 surface integral tensor (dimensionless)
  - hm_spherical: SphericalTensor of hm_tensor
  - hm_B_field: Vec3, HM B-field from this ring

### Atom post-pass
- per_type_hm_T0_sum[8]: double per ring type
- per_type_hm_T2_sum[8]: array<5> per ring type

### Physics query methods
- `TensorAtAtom(atom_index, ring_index) -> Mat3`
- `BFieldAtAtom(atom_index, ring_index) -> Vec3`
- `SumT0ByRingType(atom_index, RingTypeIndex) -> double`
- `NearestRingContribution(atom_index) -> double`

### Categorical analysis ops accounted for
- ops 274-289: Haigh-Mallion features (all)

---

## PiQuadrupoleResult

**Requires:** SpatialIndexResult, GeometryResult
**Accessor:** `conformation.PiQuadrupole()`

### What it computes
For each ring in the protein:
  For each atom within range (15.0A from ring center):

**Pi-quadrupole field:**
- Input: Ring geometry (center, normal), atom position
- Computation: Second derivative of axial quadrupole potential
  quad_factor = (3*cos^2(theta) - 1) / r^4
  Full EFG tensor: 15*dn^2*d_a*d_b/r^7 - 3*dn*(n_a*d_b+n_b*d_a)/r^5
                   - 3*d_a*d_b/r^5 + delta_ab*(3*dn^2/r^5 - 1/r^3)
- Output per atom per ring (RingNeighbourhood entry):
  - quad_tensor: Mat3, quadrupole EFG tensor
  - quad_spherical: SphericalTensor of quad_tensor
  - quad_scalar: double, (3*cos^2(theta)-1)/r^4

### Physics query methods
- `TensorAtAtom(atom_index, ring_index) -> Mat3`
- `ScalarSum(atom_index) -> double`

### Categorical analysis ops accounted for
- ops 393-408: quadrupole_field_features (all)

---

## RingSusceptibilityResult

**Requires:** SpatialIndexResult, GeometryResult
**Accessor:** `conformation.RingSusceptibility()`

### What it computes
For each ring in the protein:
  For each atom within range (15.0A from ring center):

**Ring susceptibility anisotropy (McConnell from ring center):**
- Input: Ring geometry (center, normal), atom position
- Computation: Dipolar tensor from ring center:
  T_ab = 3*d_a*d_b/r^5 - delta_ab/r^3
  McConnell scalar: (3*cos^2(theta)-1)/r^3
- Output per atom per ring (RingNeighbourhood entry):
  - chi_tensor: Mat3, dipolar tensor from ring center
  - chi_spherical: SphericalTensor of chi_tensor
  - chi_scalar: double, (3*cos^2(theta)-1)/r^3

### Physics query methods
- `TensorAtAtom(atom_index, ring_index) -> Mat3`
- `ScalarSum(atom_index) -> double`
- `NearestContribution(atom_index) -> double`

### Categorical analysis ops accounted for
- ops 409-424: ring_chi_aniso_features (all)

---

## DispersionResult

**Requires:** SpatialIndexResult, GeometryResult
**Accessor:** `conformation.Dispersion()`

### What it computes
For each ring in the protein:
  For each atom within range:

**London dispersion contacts:**
- Input: Ring geometry (vertices), atom position
- Computation: Per-vertex 1/r^6 sum with anisotropic tensor:
  T_disp_ab += 3*d_a*d_b/r^8 - delta_ab/r^6
  disp_scalar += 1/r^6
- Distance filter: R_MIN=1.5A to R_CUT=5.0A from each vertex
- Early exit: ring center distance > 8.0A (R_CUT + max ring radius)
- Output per atom per ring (RingNeighbourhood entry):
  - disp_tensor: Mat3, accumulated dispersion tensor over ring vertices
  - disp_spherical: SphericalTensor of disp_tensor
  - disp_scalar: double, sum of 1/r^6
  - disp_contacts: int, number of vertex contacts in range

### Physics query methods
- `TensorAtAtom(atom_index, ring_index) -> Mat3`
- `ContactCount(atom_index) -> int`
- `InvR6Sum(atom_index) -> double`

### Categorical analysis ops accounted for
- ops 350-362: dispersion_tensor_features (all)

---

## McConnellResult

**Requires:** SpatialIndexResult, GeometryResult
**Accessor:** `conformation.McConnell()`

### What it computes
For each bond in the protein:
  For each atom within range (via spatial index on bond midpoints):

**McConnell bond anisotropy:**
- Computation: Full dipolar tensor from bond midpoint:
  T_ab = 3*d_a*d_b/r^5 - delta_ab/r^3
  McConnell scalar = (3*cos^2(theta)-1)/r^3 (derived FROM tensor trace)
- Per-bond-category accumulation:
  PeptideCO, PeptideCN, SidechainCO, Aromatic (separate totals)
- Output per atom per bond (BondNeighbourhood entry):
  - bond_dipolar_tensor: Mat3
  - bond_dipolar_spherical: SphericalTensor
  - bond_mcconnell_scalar: double (derived from tensor)
  - bond_ref: index + BondCategory
  - distance_to_midpoint: double

The v1 code had TWO McConnell methods: one computing only scalar factors,
another using the molecular graph for full category-decomposed tensors.
In this rewrite, there is ONE McConnell calculator that produces full
tensors per bond, categorised by bond type. The scalar factor is DERIVED
from the tensor trace, never computed instead of the tensor.

### Atom post-pass
- mcconnell_co_sum: double, sum of McConnell factors for all PeptideCO bonds
- mcconnell_cn_sum: double, sum for PeptideCN bonds
- mcconnell_sidechain_sum: double, sum for SidechainCO bonds
- mcconnell_aromatic_sum: double, sum for Aromatic bonds
- mcconnell_co_nearest: double, McConnell factor of nearest PeptideCO bond
- T2_CO_nearest: SphericalTensor, nearest PeptideCO dipolar tensor decomposed
- T2_CN_nearest: SphericalTensor, nearest PeptideCN dipolar tensor decomposed
- T2_backbone: SphericalTensor, sum of all backbone bond tensors
- T2_sidechain: SphericalTensor, sum of all sidechain bond tensors
- T2_aromatic: SphericalTensor, sum of all aromatic bond tensors
- nearest_CO_midpoint: Vec3
- nearest_CO_dist: double
- nearest_CN_dist: double
- dir_nearest_CO: Vec3, normalised

### Physics query methods
- `TensorAtAtom(atom_index, bond_index) -> Mat3`
- `CategorySum(atom_index, BondCategory) -> double`
- `NearestCOContribution(atom_index) -> double`
- `BondNeighbourhoodAt(atom_index) -> const vector<BondNeighbourhood>&`

### Categorical analysis ops accounted for
- ops 185-210: Bond anisotropy (both scalar and tensor versions unified)
- ops 303-312: Full bond anisotropy (unified with McConnellResult)

---

## CoulombResult

**Requires:** ChargeAssignmentResult, SpatialIndexResult
**Accessor:** `conformation.Coulomb()`

### What it computes
For each atom with |charge| > 1e-10:
  For each atom within 20A (via spatial index):

**Coulomb electric field gradient:**
- Computation:
  E-field: E += q * dr / r^3
  EFG tensor: V_ab += q * (3*dr_a*dr_b/r^5 - delta_ab/r^3)
- Decomposed by source:
  backbone atoms: E_backbone, EFG_backbone
  sidechain atoms: E_sidechain, EFG_sidechain
  aromatic sidechain atoms (beyond CB): E_aromatic, EFG_aromatic
  solvent contribution: E_solvent = apbs_efield - coulomb_E_total

Note: long-range Coulomb (beyond 20A) handled by APBS via ApbsFieldResult.
The Coulomb calculator handles short-range vacuum contributions; APBS
handles solvated long-range. Both are available as separate features.

Note: coulomb_E_ring_proj (E_total . nearest_ring_normal) is NOT computed
here. It is computed at feature extraction time because it depends on
ring geometry data from BiotSavartResult, which CoulombResult does not
depend on. This avoids a circular dependency.

### What it stores (per atom)
- coulomb_E_total: Vec3
- coulomb_E_backbone: Vec3
- coulomb_E_sidechain: Vec3
- coulomb_E_aromatic: Vec3
- coulomb_EFG_total: Mat3 + SphericalTensor
- coulomb_EFG_backbone: Mat3 + SphericalTensor
- coulomb_EFG_sidechain: Mat3
- coulomb_EFG_aromatic: Mat3 + SphericalTensor
- coulomb_E_solvent: Vec3 (derived: APBS - vacuum)
- coulomb_EFG_solvent: Mat3 (derived: APBS - vacuum)
- coulomb_E_magnitude: double
- coulomb_E_bond_proj: double (for H atoms)
- coulomb_E_backbone_frac: double
- aromatic_E_magnitude: double
- aromatic_E_bond_proj: double (for H atoms)
- aromatic_n_sidechain_atoms: int

Fixes Bug 3 (charge loading differences) by using charges from
ChargeAssignmentResult (single authority: force field via prmtop).

### Physics query methods
- `EFieldAt(atom_index) -> Vec3`
- `BackboneField(atom_index) -> Vec3`
- `SidechainField(atom_index) -> Vec3`
- `EFGAt(atom_index) -> Mat3`
- `BackboneFraction(atom_index) -> double`

### Categorical analysis ops accounted for
- ops 211-239: Coulomb field (all)
- ops 363-383: Aromatic field features (unified with Coulomb decomposition)

---

## HBondResult

**Requires:** DsspResult, SpatialIndexResult
**Accessor:** `conformation.HBond()`

### What it computes
For each atom that is a donor or acceptor:
  Find nearest valid H-bond partner (spatial index, filtered by role):
- Partner validity: donor-acceptor pairing, not same residue,
  not sequential backbone N...O (|resnum diff| < 2, same chain)
- Computation:
  Dipolar tensor: T_ab = 3*d_a*d_b/r^5 - delta_ab/r^3
  Geometry: distance, D-H-A angle, donor/acceptor classification

### What it stores (per atom)
- hbond_nearest_dist: double, Angstroms (or 99.0 sentinel if none)
- hbond_nearest_dir: Vec3, normalised direction to partner
- hbond_nearest_tensor: Mat3, dipolar tensor
- hbond_nearest_spherical: SphericalTensor
- hbond_inv_d3: double, 1/d^3 (0 if d > 50A)
- hbond_is_backbone: bool, both atoms in backbone
- hbond_count_within_3_5A: int, count of valid partners within 3.5A
- hbond_is_donor: bool (this atom is a donor)
- hbond_is_acceptor: bool (this atom is an acceptor)

### Physics query methods
- `NearestDistance(atom_index) -> double`
- `NearestTensor(atom_index) -> Mat3`
- `IsDonor(atom_index) -> bool`
- `IsAcceptor(atom_index) -> bool`
- `CountWithin(atom_index, double radius) -> int`

### Categorical analysis ops accounted for
- ops 324-349: hbond_features (all)

---

## OrcaShieldingResult

**Requires:** nothing (loaded from files)
**Accessor:** `conformation.OrcaShielding()`

### What it computes
- Loading wild-type ORCA output files (shielding tensors per atom)
- Loading ALA mutant ORCA output files (shielding tensors per atom)
- Atom matching between WT and ALA structures
- Computing delta tensors: WT shielding minus ALA shielding per atom

This is the ground truth that the ML model learns to predict.

### What it stores (per atom)
- wt_shielding: Mat3 + SphericalTensor (ppm)
- ala_shielding: Mat3 + SphericalTensor (ppm)
- delta_shielding: Mat3 + SphericalTensor (ppm)

### Physics query methods
- `DeltaT0At(atom_index) -> double`
- `DeltaTensorAt(atom_index) -> Mat3`
- `HasTarget(atom_index) -> bool`

---

## FeatureExtractionResult

**Requires:** all physics ConformationResult types
**Accessor:** `conformation.FeatureExtraction()`

### What it computes
Each feature is a subclass of Feature with:
- name: string (for manifest and serialisation)
- irrep: IrrepType enum (L0, L1e, L1o, L2e)
- compute(atom_index, conformation): returns the feature value

Features do NOT compute physics. They READ from ConformationResult objects.

**Ring current features (47 L0, 7 L1, 10 L2 = 64 features)**

L0 scalars (47):
- per_type_G_T0[0..7]: 8 doubles, from conformation.BiotSavart().SumT0ByRingType()
  (PHE, TYR, TRP6, TRP5, TRP9, HIS, HID, HIE)
  Note: v1 had 7 types (no TRP9). Rewrite has 8.
- B_total_mag: double, from conformation.BiotSavart().TotalBField().norm()
- ring_count_3A, ring_count_5A, ring_count_8A, ring_count_12A: 4 ints cast to double
- ring_dist_1st, ring_dist_2nd, ring_dist_3rd, ring_mean_dist: 4 doubles
  (sentinel 99.0 if not enough rings)
- near1_r, near1_rho, near1_z, near1_theta, near1_inv_r3: 5 doubles
- near1_mcconnell: double, (3*cos^2(theta)-1)/r^3
- near1_Bn, near1_Brho, near1_Bphi: 3 doubles (cylindrical B-field)
- near1_G_iso: double (T0 of nearest ring G)
- near1_radius: double
- near1_has_N: double (0 or 1, from ring type NitrogenCount() > 0)
- near1_type_onehot[0..7]: 8 doubles (one-hot of nearest ring type, 8 classes)
- near1_exp_decay: double, exp(-r1/4.0)
- near2_r, near2_z, near2_G_iso, near2_mcconnell: 4 doubles
- ring_n1_dot_n2, ring_n1_cross_n2: 2 doubles
- G_iso_sum_exp: double
- G_iso_var_8A: double
- d2_d1_ratio: double

L1 vectors (7):
- jb_B_total: Vec3 (L1e, pseudovector)
- jb_B_near1: Vec3 (L1e)
- jb_B_near2: Vec3 (L1e)
- ring_normal_1: Vec3 (L1e)
- ring_normal_2: Vec3 (L1e)
- dir_to_ring_1: Vec3 (L1o, polar vector)
- dir_to_ring_2: Vec3 (L1o)

L2 tensors (10):
- per_type_G_T2[0..7]: 8 x array<5> (all L2e)
- jb_G_T2_near1: array<5> (L2e)
- jb_G_T2_exp_sum: array<5> (L2e)

Replaces ops 94-183.

**Bond anisotropy features (5 L0, 1 L1, 5 L2 = 11 features)**

L0 scalars (5):
- mcconnell_co_sum: from conformation.McConnell().CategorySum(atomIdx, PeptideCO)
- mcconnell_co_nearest: from conformation.McConnell().NearestCOContribution()
- mcconnell_cn_sum: from conformation.McConnell().CategorySum(atomIdx, PeptideCN)
- mcconnell_sidechain_sum
- mcconnell_cn_nearest_dist (or sentinel 99.0)

L1 vectors (1):
- mcconnell_dir_nearest_co: Vec3 (L1o)

L2 tensors (5):
- peptide_T2_CO_nearest: array<5> (L2e)
- peptide_T2_CN_nearest: array<5> (L2e)
- bond_T2_backbone: array<5> (L2e)
- bond_T2_sidechain: array<5> (L2e)
- bond_T2_aromatic: array<5> (L2e)

Replaces ops 185-210, 303-312.

**Coulomb field features (4 L0, 2 L1, 2 L2 = 8 features)**

L0 scalars (4):
- coulomb_E_mag: from conformation.Coulomb().EFieldAt().norm()
- coulomb_E_bond_proj: double (0.0 for non-H atoms)
- coulomb_E_ring_proj: E_total . nearest_ring_normal (computed HERE at feature
  time using data from both CoulombResult and BiotSavartResult; 0.0 if no rings)
- coulomb_E_backbone_frac: from conformation.Coulomb().BackboneFraction()

L1 vectors (2):
- coulomb_E_total: Vec3 (L1o)
- coulomb_E_backbone: Vec3 (L1o)

L2 tensors (2):
- coulomb_EFG_T2_total: array<5> (L2e)
- coulomb_EFG_T2_backbone: array<5> (L2e)

Replaces ops 211-239.

**Atom identity features (10 L0 = 10 features)**

- atom_is_H, atom_is_C, atom_is_N, atom_is_O, atom_is_S: 5 doubles (0/1)
- atom_Z: double (atomic number)
- atom_is_backbone, atom_is_amide_H, atom_is_alpha_H, atom_is_methyl: 4 doubles (0/1)

Replaces ops 240-249.

**Structural features (10 L0, 2 L1 = 12 features)**

L0 scalars (10):
- struct_acceptor_dist: double (or sentinel 99.0)
- struct_co_dist: double (or sentinel 99.0)
- struct_co_mcconnell: double (McConnell factor of nearest C=O)
- struct_n_heavy_8A: double (count of heavy atoms within 8A)
- struct_min_ring_atom_dist: double (or sentinel 99.0)
- struct_min_center_ratio: double (min_ring_atom_dist / ring_center_dist)
- struct_B_bond_cos: double (bond_direction . B_total.normalized())
- struct_ch_pi_cos: double (|bond_direction . nearest_ring_normal|)
- struct_is_aromatic_res: double (0/1)
- struct_parent_sp2: double (0/1)

L1 vectors (2):
- struct_dir_acceptor: Vec3 (L1o)
- struct_bond_direction: Vec3 (L1o)

Replaces ops 250-270.

**Residue encoding (20 L0 = 20 features)**

- res_onehot[0..19]: 20 doubles, one-hot for 20 standard amino acids
  Order: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE,
         LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL

IMPORTANT: HID/HIE/HIP map to HIS in the residue one-hot (op 272).
The protonation variant is captured by the ring type one-hot (8 classes),
not by the amino acid encoding.

Replaces ops 271-273.

**Haigh-Mallion features (9 L0, 1 L1, 9 L2 = 19 features)**

L0 scalars (9):
- per_type_hm_T0[0..7]: 8 doubles from conformation.HaighMallion().SumT0ByRingType()
- hm_T0_nearest: double

L1 vectors (1):
- hm_B_total: Vec3 (L1e)

L2 tensors (9):
- per_type_hm_T2[0..7]: 8 x array<5> (L2e)
- hm_T2_nearest: array<5> (L2e)

Replaces ops 274-289.

**Graph features (10 L0 = 10 features)**

- graph_dist_ring: double (sentinel 99.0 if -1)
- graph_dist_N: double (sentinel 99.0 if -1)
- graph_dist_O: double (sentinel 99.0 if -1)
- graph_eneg_sum_1: double
- graph_eneg_sum_2: double
- graph_is_sp, graph_is_sp2, graph_is_sp3: 3 doubles (0/1)
- graph_n_pi_bonds_3: double
- graph_is_conjugated: double (0/1)

Replaces ops 290-302.

**APBS field features (2 L0, 1 L1, 1 L2 = 4 features)**

L0 scalars (2):
- apbs_E_mag: from conformation.ApbsField().EFieldAt().norm()
- apbs_E_bond_proj: double (0.0 for non-H atoms)

L1 vectors (1):
- apbs_E_total: Vec3 (L1o)

L2 tensors (1):
- apbs_EFG_T2: array<5> (L2e)

Replaces ops 313-323.

**H-bond features (6 L0, 1 L1, 1 L2 = 8 features)**

L0 scalars (6):
- hbond_dist, hbond_inv_d3, hbond_is_donor, hbond_is_acceptor,
  hbond_is_backbone, hbond_count_3_5A

L1 vectors (1):
- hbond_dir: Vec3 (L1o)

L2 tensors (1):
- hbond_T2: array<5> (L2e)

Replaces ops 324-349.

**Dispersion features (2 L0, 0 L1, 1 L2 = 3 features)**

- disp_sum_inv_r6, disp_n_contacts: 2 L0
- disp_T2_contact: array<5> (L2e)

Replaces ops 350-362.

**Aromatic sidechain field features (3 L0, 1 L1, 1 L2 = 5 features)**

- arom_E_mag, arom_E_bond_proj, arom_n_sidechain_atoms: 3 L0
- arom_E_total: Vec3 (L1o)
- arom_EFG_T2: array<5> (L2e)

Replaces ops 363-383.

**Per-type nearest ring features (4 L0, 0 L1, 4 L2 = 8 features)**

v1 used 4 collapsed types (PHE, TYR, TRP, HIS-family). Rewrite uses
4 collapsed types for backward compatibility:
- near_dist_PHE, near_dist_TYR, near_dist_TRP, near_dist_HIS: 4 L0
  (sentinel 99.0 if no ring of that type)
- jb_G_T2_near_PHE, _TYR, _TRP, _HIS: 4 L2

The collapse: PHE->PHE, TYR->TYR, TRP5+TRP6+TRP9->TRP,
HIS+HID+HIE->HIS. This matches v1 op 387.

Replaces ops 384-392.

**Quadrupole features (1 L0, 0 L1, 1 L2 = 2 features)**

- quad_scalar_sum: from conformation.PiQuadrupole().ScalarSum()
- quad_EFG_T2: array<5> (L2e)

Replaces ops 393-408.

**Ring susceptibility features (2 L0, 0 L1, 1 L2 = 3 features)**

- chi_aniso_sum, chi_aniso_nearest: 2 L0
- chi_aniso_T2: array<5> (L2e)

Replaces ops 409-424.

**Through-bond features (2 L0 = 2 features)**

- bfs_dist_ring: from conformation.MolecularGraph().BondDistToRing()
- bfs_decay: from conformation.MolecularGraph().BFSDecay()

**Corrected sign convention features**
No additional features needed. The sign convention is enforced in the
calculators, not at feature time. All features inherit the correct convention.

### Feature count summary

| Group | L0 | L1e | L1o | L2e | Total features |
|-------|-----|------|------|------|----------------|
| Ring current | 47 | 5 | 2 | 10 | 64 |
| Bond anisotropy | 5 | 0 | 1 | 5 | 11 |
| Coulomb field | 4 | 0 | 2 | 2 | 8 |
| Atom identity | 10 | 0 | 0 | 0 | 10 |
| Structural | 10 | 0 | 2 | 0 | 12 |
| Residue encoding | 20 | 0 | 0 | 0 | 20 |
| Haigh-Mallion | 9 | 1 | 0 | 9 | 19 |
| Graph | 10 | 0 | 0 | 0 | 10 |
| APBS | 2 | 0 | 1 | 1 | 4 |
| H-bond | 6 | 0 | 1 | 1 | 8 |
| Dispersion | 2 | 0 | 0 | 1 | 3 |
| Aromatic field | 3 | 0 | 1 | 1 | 5 |
| Per-type nearest | 4 | 0 | 0 | 4 | 8 |
| Quadrupole | 1 | 0 | 0 | 1 | 2 |
| Ring chi aniso | 2 | 0 | 0 | 1 | 3 |
| Through-bond | 2 | 0 | 0 | 0 | 2 |
| **Total** | **137** | **6** | **10** | **36** | **189** |

Final e3nn irreps string: `137x0e+6x1e+10x1o+36x2e`

---

## PredictionResult

**Requires:** FeatureExtractionResult
**Accessor:** `conformation.Prediction()`

### What it computes
- Pack features into tensors (L0, L1e, L1o, L2e vectors)
- Run TorchScript model (inference via LibTorch)
- Unpack predictions

### What it stores (per atom)
- predicted_T0: double, isotropic shielding prediction (ppm)
- predicted_T2: array<5>, anisotropic shielding prediction
- confidence: double, heteroscedastic learned uncertainty (sigma)
- tier: HeuristicTier enum (REPORT, PASS, SILENT)

### Tier classification criteria
- REPORT: confidence < threshold_report AND element == H
  AND is_amide_H (the primary training target)
- PASS: confidence < threshold_pass (broader threshold, any atom type
  that has training data)
- SILENT: everything else (atoms with no training data, or predictions
  the model is not confident about)

### Physics query methods
- `PredictedT0At(atom_index) -> double`
- `PredictedT2At(atom_index) -> array<double, 5>`
- `ConfidenceAt(atom_index) -> double`
- `TierAt(atom_index) -> HeuristicTier`

---

## Operations from v1 with NO direct home (flagged)

These v1 operations are accounted for but warrant explicit documentation:

1. **op 466 (spatial matching with 1.0A threshold)**: This was used in
   extract_all to match DFT target atoms to PDB atoms when atom names
   were missing. In the rewrite, atom correspondence is handled at data
   loading (construction) using cifpp's canonical atom naming, not at
   extraction time. If spatial matching is still needed for legacy data,
   it belongs in a data import utility, not in the extraction pipeline.

2. **ops 457-458 (target T0/T2 assignment)**: Target values (DFT shielding
   tensors) are loaded via OrcaShieldingResult, not computed by the
   extraction pipeline.

3. **ops 459-461 (extract_atom_no_target)**: The rewrite does not distinguish
   "with target" and "without target" extraction. Feature extraction is the
   same regardless of whether targets exist. Targets are separate metadata
   in OrcaShieldingResult.

4. **op 462 (distance filter in extract_all)**: min_distance/max_distance
   filtering of atoms by ring proximity is a training data selection criterion,
   not an extraction operation. It belongs in the training data preparation
   pipeline.

5. **op 470 (e3nn_irreps_string)**: Generated automatically from the Feature
   registry. Not a manual string.

---

## Magic numbers consolidated

All magic numbers from v1, with their homes in the rewrite:

| Value | v1 ops | Rewrite location | Handling |
|-------|--------|------------------|----------|
| 99.0 | sentinel | Named constant: NO_DATA_SENTINEL | Unchanged |
| 1.8 A | op 54 | Bond detection (construction, OpenBabel) | Replaced by bond graph |
| 1e-10 | ops 26, 35, 150 | Named constant: NEAR_ZERO_NORM | Unchanged |
| 1e-15 | ops 228, 321 | Named constant: NEAR_ZERO_FIELD | Unchanged |
| 0.1 A | multiple | Named constant: MIN_DISTANCE | Unchanged |
| 3.0 A | op 40 | Named constant: RING_COUNT_SHELL_1 | Unchanged |
| 5.0 A | ops 41, 350 | Named constants: RING_COUNT_SHELL_2, DISP_R_CUT | Unchanged |
| 8.0 A | ops 42, 73 | Named constants: RING_COUNT_SHELL_3, PACKING_RADIUS | Unchanged |
| 12.0 A | op 43 | Named constant: RING_COUNT_SHELL_4 | Unchanged |
| 15.0 A | ops 394, 410 | Named constant: RING_CALC_CUTOFF | Unchanged |
| 4.0 A | op 90 | Named constant: EXP_DECAY_LENGTH | Unchanged |
| 3.5 A | op 335 | Named constant: HBOND_COUNT_RADIUS | Unchanged |
| 1.5 A | op 350 | Named constant: DISP_R_MIN | Unchanged |
| 50.0 A | op 340 | Named constant: HBOND_MAX_DIST | Unchanged |
| 100.0 V/A | op 468 | Named constant: APBS_SANITY_LIMIT | Unchanged |
| 2 | op 332 | Named constant: SEQUENTIAL_EXCLUSION_THRESHOLD | Unchanged |
