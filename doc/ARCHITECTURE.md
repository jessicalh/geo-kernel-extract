# Architecture: NMR Shielding Tensor System

This document describes the system as it exists in code. Every claim
is derived from the headers in src/ and verified against the test
suite (287 tests, 0 failures as of 2026-04-07). Diagrams are in
doc/diagrams/ and were generated from the code structure.

See also: doc/generated/doxygen/ for full class-level UML with members,
collaboration graphs, and include dependency trees.

---

## What this system does

This system computes classical electromagnetic geometric kernels —
ring currents, bond anisotropy, electric field gradients, dispersion,
quadrupole EFGs — at every atom of a protein, as full rank-2
tensors. Each tensor decomposes into T0 (isotropic), T1
(antisymmetric), and T2 (symmetric traceless, 5 components).

The geometric kernels are calibrated against DFT shielding tensors
from ORCA on 723 WT-ALA mutant pairs. The WT-ALA delta isolates
the aromatic contribution. ~80 calculator parameters are tuned to
minimise the T2 residual between the classical kernels and the DFT
deltas. The T2 R² — how well the calibrated classical kernels
reproduce the angular structure of DFT — is the primary measure.

This is not a chemical shift predictor. T0 (the isotropic shift) is
what existing predictors compute. This system's contribution is T2:
the angular pattern that reveals where classical models succeed and
fail, which directions the residual points, and what physics is
missing. The calibrated parameters and the T2 residual map across
hundreds of proteins are the thesis results.

The full tensor output (T0+T1+T2 per calculator per atom) is exported
as NPY arrays. These serve both the calibration pipeline and an
upstream equivariant prediction model. T2 completeness is load-bearing
for both consumers. Do not reduce calculator output to scalars.

---

## 1. Protein Model

**Diagram:** doc/diagrams/01_protein_model (svg, png, mmd source)

A Protein holds what a molecule IS, independent of geometry. It owns:
- Atoms (element, PDB atom name, residue index)
- Residues (amino acid type, protonation variant, backbone index cache)
- Bonds (two atom indices, a BondCategory enum)
- Rings (typed class hierarchy, vertex atom indices)

Protein is non-movable and non-copyable (Protein.h:32-33). This is
the back-pointer safety guarantee: conformations hold a raw
`const Protein*` that is valid for the Protein's lifetime because the
Protein cannot be moved to a new address.

### Conformations

A ProteinConformation is one geometric instance of the protein.
Positions are const after construction — new geometry means a new
conformation, not mutation of an existing one. Protein owns its
conformations via `vector<unique_ptr<ProteinConformation>>`.

Four typed subclasses record provenance:

| Type | Created by | Metadata |
|------|-----------|----------|
| CrystalConformation | BuildFromPdb, BuildFromProtonatedPdb | resolution, R-factor, temperature, PDB ID |
| PredictionConformation | BuildFromOrca | method, confidence |
| MDFrameConformation | BuildFromGromacs | walker, time (ps), Boltzmann weight, RMSD, Rg |
| DerivedConformation | Protein::AddDerived | description string |

Every conformation holds N ConformationAtoms (one per protein atom)
with const positions and public computed fields. Each field has exactly
one writer — the ConformationResult that owns it. This is enforced by
the singleton guarantee: one result per type per conformation, attached
once, never replaced (ProteinConformation.h, results_ map keyed by
type_index).

### Ring type hierarchy

Ring types are classes, not data. Calculator code is ring-type-agnostic:
`ring.Intensity()`, `ring.JBLobeOffset()`, `ring.NitrogenCount()`.
No switch statements on ring type in calculator code.

Eight types in three size categories (Ring.h):
- SixMemberedRing: PheBenzene, TyrPhenol, TrpBenzene
- FiveMemberedRing: TrpPyrrole, HisImidazole, HidImidazole, HieImidazole
- FusedRing: IndolePerimeter (TRP 9-atom)

The three histidine ring types (His, Hid, Hie) encode protonation
state in the type system. Protonation is not a string annotation —
it determines which Ring subclass is constructed during
DetectAromaticRings.

### Topology

CovalentTopology is resolved once during FinalizeConstruction
(Protein.h:146). It uses OpenBabel for bond detection from 3D
coordinates and covalent radii, then classifies bonds into
BondCategory enums (PeptideCO, PeptideCN, BackboneOther,
SidechainCO, Aromatic, Disulfide, SidechainOther, Unknown)
using the backbone index
cache — typed properties, not string comparisons.

---

## 2. Provenance Paths

**Diagram:** doc/diagrams/02_provenance_paths (svg, png, mmd source)

Five kinds of data enter the system. Each has a different protonation
strategy and charge source. All produce the same BuildResult
(BuildResult.h): a Protein, a ChargeSource, and an integer net charge.
After building, the same OperationRunner::Run pipeline processes every
protein identically.

### Path 1: PDB of unknown provenance

**Builder:** `BuildFromPdb(path, pH)` (PdbFileReader.h:24)

The PDB may or may not have hydrogens. reduce (Richardson lab,
linked as a C++ library) strips any existing H atoms and rebuilds
them at the specified pH, optimising flips and rotamers. Charges
come from a ParamFileChargeSource (ff14SB parameter file).

This is the path for bare PDB files where protonation state is
unknown or unreliable.

### Path 2: PDB that is already protonated

**Builder:** `BuildFromProtonatedPdb(path)` (PdbFileReader.h:29)

The PDB already contains all hydrogen atoms (e.g., output from a
previous reduce run, or a PDB from GROMACS that was externally
reprotonated for experimental comparison). reduce is not called.
ProtonationDetectionResult reads existing H atoms to determine
protonation states. Charges from ParamFileChargeSource (ff14SB).

### Path 3: ORCA from tleap (calibration proteins)

**Builder:** `BuildFromOrca(files)` (OrcaRunLoader.h:52)

Calibration data: tleap chose the protonation, built the prmtop,
produced XYZ coordinates for ORCA. The prmtop is the authority for
atom names, residues, and charges (PrmtopChargeSource). All 723
proteins in the consolidated training set use this path.

### Path 4: ORCA from GROMACS pose

Same builder as Path 3 (`BuildFromOrca`), but the origin of the
structure is a GROMACS MD snapshot rather than an AlphaFold
prediction. The ORCA run and prmtop were still prepared via tleap.
From the loader's perspective, the data format is identical.

### Path 5: GROMACS ensemble (trajectory)

**Builder:** `BuildFromGromacs(paths)` (GromacsEnsembleLoader.h:50)

The TPR binary topology (read via libgromacs) provides atom names
(CHARMM naming), residues, and per-atom charges. Charges are
wrapped in a PreloadedChargeSource (ChargeSource.h:174).
ForceField is declared by the caller (FleetPaths.force_field),
not inferred from the data.

Pose PDBs provide positions only — one MDFrameConformation per
pose. A typical ensemble has 10 frames. The full pipeline
(including MOPAC) runs independently on each frame.

### Convergence

All five paths produce the same types: a Protein with typed
conformations, a ChargeSource matched to that Protein's atoms, and
a net charge. OperationRunner::Run does not know or care which
builder created the protein.

The ChargeSource hierarchy (ChargeSource.h) enforces this:
ParamFileChargeSource (ff14SB), PrmtopChargeSource (AMBER prmtop),
GmxTprChargeSource (GROMACS TPR), PreloadedChargeSource (pre-read
values). Each is a distinct type that knows its force field
(ForceField enum) and how to match charges to atoms.

---

## 3. Calculator Pipeline

**Diagram:** doc/diagrams/03_calculator_pipeline (svg, png, mmd source)

OperationRunner::Run (OperationRunner.h:68) sequences all results
in dependency order. It does not hold state, cache results, or
decide what to compute. If a step's prerequisites are on the
conformation, the step runs. If not, it is skipped.

### Tier 0: Foundation and external tools

These provide the substrate that calculators need:

| Result | What it does | Formal dependencies | Runner prerequisite |
|--------|-------------|--------------------|--------------------|
| GeometryResult | Ring geometry (SVD normals), bond geometry, global measures | none | — |
| SpatialIndexResult | 15A neighbour lists via nanoflann k-d tree | GeometryResult | — |
| EnrichmentResult | AtomRole, hybridisation, boolean properties | none | runs after GeometryResult |
| ProtonationDetectionResult | Reads H atoms to assign protonation variants | none | — |
| ChargeAssignmentResult | Per-atom charges from the ChargeSource | none | — |
| DsspResult | Secondary structure, phi/psi, SASA, H-bond partners (libdssp) | none | runs after GeometryResult |
| ApbsFieldResult | Solvated E-field and EFG via Poisson-Boltzmann (APBS) | ChargeAssignmentResult | — |
| MolecularGraphResult | Through-bond distances to rings, N, O via BFS | SpatialIndexResult | — |
| MopacResult | PM7+MOZYME: Mulliken charges, orbital populations, Wiberg bond orders | none | gated on charges |

MopacResult runs on every conformation. On a 167-atom protein
(Q9UR66), it takes ~14 seconds. On 1231 atoms (1UBQ), longer.
This is a semiempirical quantum calculation, not a lookup.

### Tier 1: Classical geometric kernel calculators

Eight calculators, each producing full rank-2 tensor output
(Mat3 + SphericalTensor) at every atom. A SphericalTensor holds
the irreducible decomposition: T0 (isotropic scalar, 1 component),
T1 (antisymmetric, 3 components), T2 (symmetric traceless, 5
components). T0 is the classical chemical shift. T2 is the angular
structure — the thesis's primary analytical result.

Each kernel is a closed-form function of geometry with a literature
derivation. Equations below are exactly what the code computes,
verified against the .cpp files.

**Biot-Savart** (Johnson & Bovey 1958). Induced B-field from ring
current modelled as two circular loops offset above and below the
ring plane. The geometric kernel G relates the B-field to the
shielding tensor via the ring normal.

    G_ab = -n_b · B_a · PPM_FACTOR

PPM_FACTOR = 10⁶ (converts SI field ratios to parts per million).
B is computed by numerical integration of the Biot-Savart law
around both loops. Source geometry: ring vertices (double loop).

**Haigh-Mallion** (Haigh & Mallion 1973). Shielding from ring
current via surface integral over the ring face, using 7-point
Gaussian quadrature on a fan triangulation.

    H_ab = ∫∫ [3 ρ_a ρ_b / |ρ|⁵ - δ_ab / |ρ|³] dS
    G_ab = -n_b · (H · n)_a

H is the raw integral (symmetric, traceless). The effective
B-field V = H·n reduces it to the rank-1 shielding kernel G.

**McConnell** (McConnell 1957). Magnetic anisotropy of a bond
modelled as a point dipole at the bond midpoint. The full
asymmetric tensor, not the commonly cited traceless K_ab.

    M_ab / r³ = [9 cosθ d̂_a b̂_b - 3 b̂_a b̂_b - (3 d̂_a d̂_b - δ_ab)] / r³

where d̂ is the unit vector from midpoint to atom, b̂ is the
bond direction, θ is their angle. Source geometry: bond midpoints.

**Coulomb** (Buckingham 1960). Electric field and field gradient
from point charges, with Coulomb constant k_e = 14.3996 V·A.

    E_a = k_e · Σ_j q_j (r_i - r_j)_a / |r_i - r_j|³
    V_ab = k_e · Σ_j q_j [3 (r_i-r_j)_a (r_i-r_j)_b / |r_i-r_j|⁵ - δ_ab / |r_i-r_j|³]

EFG is traceless (Gauss's law). Decomposed into backbone,
sidechain, aromatic, and solvent (APBS minus vacuum) contributions.
Source geometry: all charged atoms.

**Ring susceptibility** (McConnell 1957). Whole-ring diamagnetic
susceptibility modelled as a point dipole at the ring center.
Same tensor form as bond McConnell with b̂ replaced by n̂.

    M_ab / r³ = [9 cosθ d̂_a n̂_b - 3 n̂_a n̂_b - (3 d̂_a d̂_b - δ_ab)] / r³

Source geometry: ring centers.

**H-bond dipolar** (Barfield & Karplus 1969, Cornilescu & Bax 1999).
Dipolar tensor to H-bond partner from DSSP Kabsch-Sander criterion.
Same tensor form as McConnell with b̂ replaced by ĥ (H-bond direction).

    M_ab / r³ = [9 cosθ d̂_a ĥ_b - 3 ĥ_a ĥ_b - (3 d̂_a d̂_b - δ_ab)] / r³

Source geometry: DSSP H-bond donor-acceptor pairs.

**Pi-quadrupole** (Buckingham 1959, Stone T-tensor). EFG from an
axial quadrupole representing the π-electron density above and
below the ring plane. Fourth derivative of 1/r.

    G_ab = 105 d_n² d_a d_b / r⁹ - 30 d_n (n_a d_b + n_b d_a) / r⁷
         - 15 d_a d_b / r⁷ + 6 n_a n_b / r⁵ + δ_ab (3/r⁵ - 15 d_n²/r⁷)

where d_n = d · n̂ (displacement projected onto ring normal,
in Angstroms; d is unnormalized). Pure T2 — trace is zero.
Source geometry: ring centers.

**Dispersion** (London 1937). Van der Waals 1/r⁶ interaction summed
over ring vertices, modulated by a CHARMM switching function
(Brooks et al. 1983) between R_switch = 4.3 A and R_cut = 5.0 A.

    K_ab = S(r) · (3 d_a d_b / r⁸ - δ_ab / r⁶)

Traceless per vertex. Source geometry: ring vertices.

### MOPAC-derived calculators

Two calculators use the same geometric kernels with quantum-chemical
weights from MopacResult:

**MopacCoulomb.** Same Coulomb kernel as above, but charges are MOPAC
PM7 Mulliken charges instead of ff14SB. Captures polarisation effects
that the fixed-charge force field does not.

**MopacMcConnell.** Same McConnell kernel, but each bond's contribution
is weighted by the MOPAC Wiberg bond order. Bonds with higher electron
density contribute proportionally more.

### DFT reference (when available)

OrcaShieldingResult loads per-atom shielding tensors from an ORCA
NMR output file (dia + para + total, full Mat3 + SphericalTensor).
MutationDeltaResult compares WT and ALA conformations to produce
per-atom delta tensors at all irrep levels (T0, T1, T2).

These are not calculators — they are DFT ground truth. The WT-ALA
delta isolates the aromatic contribution: what changes when you
mutate an aromatic residue to alanine. The calibration pipeline
tunes ~80 calculator parameters to minimize the T2 residual between
classical geometric kernels and these DFT deltas across 723 proteins.

### Dependency enforcement

Each ConformationResult declares its dependencies as type_index
values (ConformationResult.h:31). AttachResult checks them at
attach time and fails with a diagnostic if any are missing. This
is the compile-time-like guarantee that results cannot be attached
in wrong order.

---

## 4. Filters and GeometryChoice

**Diagram:** doc/diagrams/04_geometry_choice_flow (svg, png, mmd source)

### KernelFilterSet

Every calculator that evaluates a geometric kernel holds a
KernelFilterSet (KernelEvaluationFilter.h). Before computing a
kernel at a field point, the calculator builds a
KernelEvaluationContext from already-computed distances and topology,
then passes it to AcceptAll(). All filters must accept.

Five filters, each with a physics reason documented in the header:

| Filter | Criterion | Physics |
|--------|----------|---------|
| MinDistanceFilter | distance >= 0.1 A | 1/r^n kernels diverge numerically |
| DipolarNearFieldFilter | distance > source_extent/2 | Multipole expansion invalid inside source |
| SelfSourceFilter | atom != source atom | Field undefined at source itself |
| SequentialExclusionFilter | sequence separation >= 2 | Through-bond coupling dominates |
| RingBondedExclusionFilter | atom not bonded to ring vertex | Through-bond, not through-space |

Each filter has a Name(), Description(), and RejectReason(ctx) that
states specific values at the point of rejection.

### GeometryChoice recording

During Compute(), calculators record geometric decisions via
GeometryChoiceBuilder (GeometryChoice.h). The builder enforces that
all entity and number additions happen inside a lambda — keeping
recording code visually separate from physics code.

A GeometryChoice records:
- Which entities (atoms, rings, bonds) were involved
- Their roles (Source, Target, Context)
- Their outcomes (Included, Excluded, Triggered, NotTriggered)
- Named numbers (distance, intensity, threshold — with units)
- Which filter rejected, if any
- An optional field sampler (closure that evaluates the physics
  at any 3D point, for UI visualisation)

Choices are appended to the conformation's flat geometry_choices
vector (ProteinConformation.h:143). The UI walks this list,
follows pointers back into the live model objects, and draws.

### OperationLog (UDP)

All operations log structured JSON over UDP (OperationLog.h).
The log is always on. Every result logs its start, completion time,
key metrics, and filter rejection counts. A listener
(ui/udp_listen.py) captures the stream for real-time diagnostics.

---

## 5. Relationship to PATTERNS.md

PATTERNS.md describes what holds this system together and what
breaks it. This document describes what the system IS.

The core correspondences:

| PATTERNS.md pattern | Where it lives in the architecture |
|--------------------|------------------------------------|
| ConformationAtom: private construction, typed fields | ConformationAtom.h:299 — private constructor, friend ProteinConformation |
| Template result access | ProteinConformation.h:66-101 — Result\<T\>(), HasResult\<T\>() |
| Singleton guarantee | ProteinConformation.h:150 — results_ map keyed by type_index |
| Dependency declaration | ConformationResult.h:31 — Dependencies() returns type_index vector |
| Compute() factory | Every result: static Compute() returns unique_ptr |
| Both representations | ConformationAtom.h — every tensor field has Mat3 + SphericalTensor |
| Ring type virtual interface | Ring.h — Intensity(), JBLobeOffset(), TypeName() are virtual |
| Equations in comments, Eigen in code | Every calculator .cpp — equation in comment, Eigen expression below |
| Return codes, not exceptions | BuildResult.h:33 — Ok(), error string. Four catch blocks exist at external library boundaries (cif++, DSSP, stoi); result codes are the pattern for all internal code. |
| Protein is non-movable | Protein.h:32-33 — move constructor deleted |
| Protonation flows through Protein | Five builder paths, Ring type hierarchy encodes protonation |
| Typed charge sources | ChargeSource.h — four concrete implementations, no string dispatch |
| No adapter/wrapper classes | UI reads Protein and ConformationAtom directly (ui/CLAUDE.md) |
| Configuration objects: physics constants at namespace scope | PhysicalConstants.h — constexpr values with comments citing sources |

---

## 6. Files

54 headers in src/. Key entry points:

| File | Purpose |
|------|---------|
| Protein.h | Sequence, topology, conformation ownership |
| ProteinConformation.h | Positions, ConformationAtom, result attachment |
| ConformationAtom.h | All per-atom computed fields (one writer each) |
| ConformationResult.h | ABC for results, WriteFeatures, dependency contract |
| BuildResult.h | Common output of all builders |
| PdbFileReader.h | BuildFromPdb, BuildFromProtonatedPdb |
| OrcaRunLoader.h | BuildFromOrca |
| GromacsEnsembleLoader.h | BuildFromGromacs, RunAllFrames |
| OperationRunner.h | Run, RunMutantComparison, RunEnsemble |
| ChargeSource.h | Four typed charge source implementations |
| KernelEvaluationFilter.h | Five filters, KernelFilterSet |
| GeometryChoice.h | Runtime recording of calculator decisions |
| Ring.h | Ring type hierarchy (8 types, 3 size categories) |
| OperationLog.h | Structured JSON/UDP logging |
| Types.h | Vec3, Mat3, SphericalTensor, all enums |
