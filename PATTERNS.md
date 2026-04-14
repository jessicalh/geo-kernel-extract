# Patterns: What Holds This System Together

This is a physics model, not a software system. If the code reads like
infrastructure, it is wrong. If it reads like equations with data flow,
it is right. Maximise physics depth. Keep the code simple, readable,
and traceable.

You as an AI read code and act upon the name. Like ancient magic the
name is the thing. Humans have similar wiring. If the name of a model
is "HIE" then the thing is a meaningless code. That is why we build
models in programming languages. It is nice they do things quickly but
it is more important they describe the world. Otherwise what they do
quickly has no meaning. HieImidazoleRing with NitrogenCount() == 2
and Aromaticity() == Weak describes the world. "HIE" is a string that
only means something if you already know what it means.

---

**Trajectory streaming pattern:** GromacsProtein (adapter +
accumulators), GromacsFrameHandler (frame lifecycle), free-standing
conformations, two-pass scan/extract. GromacsRunContext holds
trajectory-level state (bonded params from TPR, preloaded EDR
energy frames, cursor position). Owned by GromacsProtein, advanced
by the frame handler per frame. Fully documented in
spec/ENSEMBLE_MODEL.md. All patterns below apply to every path
including trajectory — the streaming classes are infrastructure
around the same ConformationResult / OperationRunner pipeline.

## The System in One Paragraph

This is a one-way analysis system. It loads a protein, creates typed
conformations, and performs successive calculations on them. Each
calculation is a ConformationResult that attaches to a conformation
and writes to typed fields on ConformationAtom. Classical calculators
produce geometric kernels — full rank-2 tensors decomposed into
irreducible representations: T0 (isotropic, 1 component), T1
(antisymmetric, 3 components), T2 (symmetric traceless, 5 components).
The T2 angular structure is the thesis's primary result. It shows
where and in which direction each classical model fails to match DFT.
A calibration pipeline (~80 parameters tuned against DFT WT-ALA
deltas) provides the connection between geometric kernels and quantum
chemistry. Nothing flows backward. Everything accumulates forward.

---

## Core Patterns

### 1. ConformationAtom: private construction, typed fields

ConformationAtom has a private constructor. Only ProteinConformation
can create them. The set of atoms is determined at construction from
the Protein's atom list. Position is const. Computed fields are public
and default-initialised. Identity (element, bonds, residue) goes
through the Protein back-pointer, never duplicated on ConformationAtom.

```cpp
// Identity — from Protein (never changes with geometry)
const auto& identity = conf.Protein().AtomAt(42);
Element elem = identity.element;

// Computed data — from ConformationAtom (geometry-dependent)
const auto& ca = conf.AtomAt(42);
Vec3 B = ca.total_B_field;
SphericalTensor bs = ca.bs_shielding_contribution;
```

### 2. Template result access

`conf.Result<T>()`, `conf.HasResult<T>()`, `conf.AllResults()`.
Adding a new ConformationResult type never modifies ProteinConformation.
The template mechanism is the ONE place templates are used in this system.

```cpp
auto& bs = conf.Result<BiotSavartResult>();
if (conf.HasResult<MopacResult>()) {
    auto& mopac = conf.Result<MopacResult>();
}
```

### 3. Singleton guarantee

One ConformationResult per type per conformation. Computed once,
attached once, never replaced. No other result type writes to its
fields on ConformationAtom. Different parameters means a different
conformation (copy-and-modify), not a second result on the same
conformation. This is why field accumulation on ConformationAtom is
safe — each field has exactly one writer.

### 4. Dependency declaration

Each result declares its dependencies as type_index values. Checked
at AttachResult() time. If a dependency is missing, the attach fails
with a diagnostic message stating what is attached, what is missing.

```cpp
std::vector<std::type_index> Dependencies() const override {
    return { typeid(SpatialIndexResult), typeid(GeometryResult) };
}
```

### 5. Compute() factory

The only way to create a ConformationResult. Static method, takes
the conformation, returns unique_ptr. Caller attaches it.

```cpp
auto bs = BiotSavartResult::Compute(conf);
conf.AttachResult(std::move(bs));
```

### 6. Shielding contribution contract

Every classical calculator stores BOTH geometric output (natural units,
for features and reuse) AND a shielding_contribution SphericalTensor
(ppm, for calibration comparison). Both are exported via WriteFeatures()
as NPY arrays. The calibration pipeline reads these alongside DFT
delta tensors to tune parameters and measure the T2 residual.

### 7. Fields typed by functional analysis

ConformationAtom fields are organised by the result type that writes
them. BiotSavartResult fields, McConnellResult fields, CoulombResult
fields. Not physics-category groups. If a new layer introduces new
results, add new fields to ConformationAtom. Do not subdivide by
physics categories.

### 8. Diagnostic error messages

Error messages state the world at the point of failure. Not advice.
State, values, what was tested.

```
AttachResult: McConnellResult requires SpatialIndexResult.
Attached results: [GeometryResult, DsspResult, ChargeAssignmentResult].
Missing: SpatialIndexResult.
```

```
BiotSavartResult::AtomSummary: atom_index=5042, atom_count=5041.
```

### 9. Return codes, not exceptions

Functions that can fail return a status. Callers check it. No exception
hierarchies. No catch blocks. No unwinding. If APBS fails, the function
returns an error status and the caller decides what to do. If a bond
detection step produces unexpected results, it returns a diagnostic.
Check your return codes like an adult.

Do not create exception class hierarchies. Do not use std::expected.
Do not use try/catch in calculator or result code. Simple return values,
simple checks, simple control flow.

### 10. Equations in comments, Eigen in code

Every calculator has a mathematical equation. The code and the equation
live side by side:

```cpp
// McConnell: T_ab = Delta_chi * (3 * d_a * d_b / r^5 - delta_ab / r^3)
Mat3 T = delta_chi * (3.0 * d * d.transpose() / std::pow(r, 5)
                     - Mat3::Identity() / std::pow(r, 3));
```

A physicist verifies the implementation by comparing the two lines.
The Eigen expression reads like the math. No intermediate variables
hiding the structure. No utility functions obscuring the formula.

### 11. Both representations always present

Every tensor is stored as BOTH Mat3 and SphericalTensor. No conversion
at point of use. No lazy evaluation. No "derive it later." The cost is
14 extra doubles per tensor. At 128 GB, this is noise. Every field that
stores a scalar derived from a tensor MUST have a companion
SphericalTensor. Upstream models will explore this data at every irrep
level. A ConformationAtom with scalar-only fields is INCOMPLETE.

```cpp
atom.total_G_tensor = G;
atom.total_G_spherical = SphericalTensor::Decompose(G);
```

### 12. Ring type virtual interface

Calculator code is ring-type-agnostic. Use `ring.Intensity()`, not a
switch statement. When TOML-calibrated intensities override literature
defaults, they are indexed by `ring.TypeIndexAsInt()`. The calculator
does not need to know which ring type it is processing.

```cpp
double I = ring.Intensity();         // -12.0 for PHE, -5.16 for HIS
double d = ring.JBLobeOffset();    // 0.64 for PHE, 0.50 for HIS
```

---

## Anti-Patterns: Code That Will Be Rejected

### String traversals for atom or structure identity

Traversals of atoms, residues, or rings must ALWAYS be traversals of
objects which contain their properties. NEVER rely on naming to
determine identity, origin, or state. If a name is being used to
distinguish — for example — a protonation state, that comparison must
become an explicit, commented method which is called.

Translations from one naming system to another (e.g., CHARMM to
Standard atom names) must be explicit, done in one place, clearly
labelled with the source and destination environments, and performed
all at once to the shared library standard. The NamingRegistry pattern
in the old code is the model for this. Once translated, the typed
objects carry their properties. No further string work.

Our objective is not to move at the blazing fast speed of string search
but to perform understandable and traceable operations on typed objects.
Code that uses string comparisons to identify atoms, residue types, ring
types, bond categories, or protonation states will be rejected.

### Objects answer questions about themselves

Every piece of information an atom or residue needs must be a method
on the atom or residue object. Not a library call with the atom's
name as a parameter. Not a lookup table indexed by a string extracted
from the atom. A method on the object that returns a typed value.

```cpp
// WRONG: opaque library dive, unreadable, string-based
double r = GetCovalentRadius(GetElementSymbol(compound.atom(name).type()));
bool arom = CheckAromaticity(LookupCompound(residue.code()), atomName);

// RIGHT: the object knows its own properties
double r = atom.CovalentRadius();
bool arom = residue.IsAromatic();
```

The first pattern is efficient optimisation catnip: the data is
technically accessible through a chain of library calls. But a reader
sees five nested function calls and has no idea what the data path is.
Future agents reading this code will not understand what happened and
will not maintain it correctly.

If the information comes from cifpp, from the AminoAcidType table,
or from IupacAtomIdentity — it was resolved at construction time and
stored as a typed property on the object. At runtime, the object
answers questions about itself. Period.

Grep markers: nested function calls where an atom name or residue code
string is passed as a parameter to a lookup function inside a
calculator, result, or feature extractor.

### cifpp and CCD: construction-time only

cifpp is the authority for PDB parsing and Chemical Component Dictionary
data. Its API is string-based — you give it a residue code string and
get back atom name strings and bond pair strings. This is correct and
necessary AT LOADING TIME.

After loading, no code ever calls cifpp again. The typed objects
(AminoAcidType, IupacAtomIdentity, NamingRegistry, Ring type classes,
AtomRole enum) are the runtime authorities. They were built from cifpp
data during PDB loading. They ARE the typed boundary.

An agent that calls `cif::compound_factory` inside a calculator or
feature extractor is diving back through the string boundary. That
code will be rejected. If you need to know whether an atom is aromatic,
call `ring.NitrogenCount()` or check `identity.role`. If you need to
know what protonation variant a HIS is, read the ring type class. The
typed objects already have the answer.

### The naming boundary: cross once, then typed objects forever

Translations between naming systems (PDB, AMBER, CHARMM, ORCA) happen
in one place: the NamingRegistry. It translates at the tool boundary —
when we read a PDB, when we prepare input for tleap, when we parse ORCA
output. Each translation is explicit, labelled with source and
destination environments, and performed all at once to the shared
library standard.

After translation, typed objects carry their properties. Element is an
enum. AtomRole is an enum. BondCategory is an enum. RingTypeIndex is an
enum. AminoAcid is an enum. No further string work.

The AminoAcidType table is the single authority for amino acid
chemistry: atoms, rings, chi angles, protonation variants. Every agent
that needs to know "what atoms does this residue have?" reads
AminoAcidType, not a string table and not cifpp.

Grep markers: `== "H"`, `== "CA"`, `== "PHE"`, `find("HIS")`,
`atom_name`, `res_name` in any calculator, result, or feature code.
`compound_factory` or `cif::` anywhere outside PDB loading code.

### Adapter, wrapper, or bridge classes

In this system there is almost nothing to adapt. The Protein has
identity. The ConformationResult has computed data. Two lines of access,
not a class. Every adapter in this system is a layer between you and
the physics. If you find yourself writing a class whose name contains
Adapter, Wrapper, Proxy, Helper, or Bridge, you are solving a problem
that does not exist here.

### Template metaprogramming beyond Result<T>()

The Result<T>() mechanism is the ONE template pattern. Everything else
is direct typed access: `atom.total_B_field` is a Vec3 member.
`ring.Intensity()` is a virtual method. No tag dispatch. No policy
templates. No CRTP. No SFINAE. If you are writing `template<` anywhere
outside ProteinConformation.h, stop and use a plain function.

### Shadow data structures

Do not create internal vectors, structs, or maps inside a
ConformationResult that mirror fields from ConformationAtom. The data
already has a home. Compute directly into the target fields on
ConformationAtom. A parallel structure means data exists in two places
and the copy step is where bugs live.

### Premature optimisation

Almost all optimisation implies simplifying for speed or memory. Our
goal is the opposite: maximise physics depth. We have 128 GB of RAM
and one protein at a time. Do not use memory pools, SIMD intrinsics,
cache-line alignment, custom allocators, or packed representations.
Use std::vector, direct member variables, and move on.

Grep markers: `alignas`, `pool`, `arena`, `SIMD`, `__attribute__` in
non-external-library code.

### Exception hierarchies

Do not create exception classes. Do not inherit from std::exception.
Do not use std::expected or std::error_code. Return a status. Check it.
If something fails, say what failed and what the values were. Simple
control flow. The system is one-way and processes one protein at a time.
There is nothing to unwind. Four catch blocks exist at external library
boundaries (cif++ parsing, DSSP, stoi, UDP socket) — these are
acceptable at the boundary. Do not add catch blocks in calculator,
result, or pipeline code.

### The utility namespace

Do not create files named Utility.h, Helpers.h, Utils.cpp, Common.h.
Do not create namespaces named util, helpers, common, or misc. If a
function belongs to a type, put it on the type. If it does not belong
to any type, it probably does not need to exist. SphericalTensor has
Decompose(). Ring has intensity(). Vec3 has normalized(). These exist.

### Configuration objects for physics constants

The cutoff radius is 15A because that is where ring current effects
become negligible. The decay length is 4A because that is the
characteristic scale. These are physics constants defined in the
OBJECT_MODEL, not user configuration. Named constants at namespace
scope with a comment citing the source:

```cpp
// Constitution: 15A cutoff for ring current calculations
constexpr double RING_CURRENT_CUTOFF_A = 15.0;
```

The ~80-93 tuneable calculator parameters are the exception — they
are genuinely tuneable, calibrated against T2 R² on DFT WT-ALA
deltas. These are what make the system a calibrated scientific
instrument rather than a collection of literature defaults. The
calibration pipeline (Python e3nn in learn/c_equivariant/) discovers
optimal values; they enter the C++ system as TOML configuration.
The geometric kernels and their T2 structure are the foundation.

### Where citable constants live

Two homes, no exceptions:

- **`calculator_params.toml`** (via `CalculatorConfig::Get()`) — tuneable
  parameters: cutoffs, intensities, guard thresholds.  These are what
  the calibration pipeline optimises.
- **`PhysicalConstants.h`** — reference data from published sources:
  Bondi radii, D4 EEQ element parameters, unit conversions, quadrature
  weights.  Not tuneable.  Every entry has a citation comment with
  author, year, and DOI or ISBN.

If a number comes from a paper and the thesis must cite it, it goes in
one of these two places.  Nowhere else.  No local `constexpr`, no
switch statements with unlabelled values, no "reference TOML" files.

### Over-decomposing into tiny functions

The Johnson-Bovey B-field calculation is ONE physics operation. It
takes a ring and an atom position and produces a B-field, a geometric
kernel, and cylindrical coordinates. When you decompose this into 15
five-line functions, the physics narrative is destroyed. You cannot
read the computation linearly. One function per physics concept. The
Biot-Savart inner loop can be 30-40 lines of Eigen math with the
equation in a comment. That IS the appropriate size.

---

## T2 Completeness — This Is Not Optional

Every calculator must produce full rank-2 tensor output at all irrep
levels (T0, T1, T2). T0-only results are incomplete and will be
rejected. The T2 angular structure is not a nice-to-have — it is
the thesis's primary analytical result.

A geometric kernel evaluated at an atom produces a 3x3 tensor. That
tensor decomposes into T0 (how much), T1 (which way), and T2 (what
angular shape). T0 is the classical chemical shift — every existing
NMR predictor computes this. T2 is what this system adds: the angular
pattern that reveals where and why the classical model breaks down.

If you find yourself reaching for a scalar summary of a calculator's
output, stop. The scalar is T0. You are discarding T2. The five T2
components per calculator per atom are the features the calibration
pipeline tunes against DFT AND the tensor inputs that the upstream
equivariant prediction model consumes. T2 completeness is load-bearing
for both. Without it, the system is a worse version of existing
chemical shift predictors.

### Per-calculator T2 comparison

Each calculator stores its T2 contribution. After all calculators
attach, you can compare: does Biot-Savart or Haigh-Mallion better
explain the T2 pattern near PHE rings? This comparison is possible
only because both store full tensors.

### BS-HM T2 redundancy

Biot-Savart (line integral, rank-1) and Haigh-Mallion (surface
integral, rank-2) make the same T2 prediction (cosine similarity
0.999 across 279K atoms). This is a finding, not a bug — it shows
the two mathematical approximations converge in their angular
structure despite different formulations.

### Per-type T2 decomposition

Per-ring-type T2 arrays (8 ring types × 5 T2 components) are stored
for BiotSavart, HaighMallion, PiQuadrupole, RingSusceptibility, and
Dispersion. These are the features that distinguish which ring type
produces which angular pattern at which atom.

---

## C++ Rules

**Standard:** C++17. No C++20 features.

**Ownership:** unique_ptr for polymorphic objects (ConformationResult,
ProteinConformation). Values for data structs (ConformationAtom fields,
SphericalTensor, RingNeighbourhood). Raw pointers for non-owning
back-references (conformation → protein). No shared_ptr anywhere.

**auto:** Use for iterator types, make_unique results, and range-for
loop variables. Do NOT use for function return types. Spell out
parameter types.

**const:** Positions are const. Element is const. Ring type is const
after protonation. Bond category is const. Const after construction
is a statement about the physics, not defensive programming.

**Naming:**
- Types: PascalCase (BiotSavartResult, ConformationAtom)
- Methods: PascalCase (TotalBField, AtomSummary)
- Private members: snake_case_ with trailing underscore (atom_data_)
- Constants: UPPER_CASE (RING_CURRENT_CUTOFF_A)
- Enums: PascalCase values (RingTypeIndex::PheBenzene)
- Files: PascalCase.h / PascalCase.cpp matching the primary class

**Headers:** One class per header. Forward-declare when a pointer or
reference suffices. No `using namespace std;` in headers.

**CMake:** Single library target (nmr_shielding) with all .cpp files
listed in CMakeLists.txt. Adding a new result = add the .cpp to the
source list in CMakeLists.txt. Always build from `build/`, never from
root (`cd build && cmake .. && make -j$(nproc)`).

---

## What To Read

- **CONSTITUTION.md** — principles, sign conventions, inviolable rules.
- **OBJECT_MODEL.md** — the concrete types. This is what you code against.
- **This file** — what breaks the system and why.
- **The code** — src/ is the ground truth. The calculators demonstrate
  every pattern described here. Read one (e.g. McConnellResult.cpp)
  before writing another.
- **learn/bones/** — historical design docs. Reference only.

---

## Lessons Learned

### 1. Back-pointer safety: Protein is non-movable (Protein.h:32-33)

Protein owns conformations via `vector<unique_ptr<ProteinConformation>>`.
Each conformation holds a raw `const Protein*` back-pointer. If the
Protein moved, the back-pointers would dangle.

The fix: `Protein(Protein&&) = delete`. Proteins live on the heap
via unique_ptr and never move. Since every Protein is constructed by
a builder function that returns it inside a BuildResult, there is no
use case for moving a Protein after construction.

**Pattern:** If children hold raw back-pointers to a parent, either
fix the pointers on move or delete the move constructor. We chose
delete because the construction pattern makes move unnecessary.

### 2. Bond classification at the typed boundary

Bond detection produces bonds. Bond CLASSIFICATION determines what
those bonds ARE (PeptideCO, PeptideCN, etc.). The classification
must use typed properties, not string comparisons.

**Wrong:** `if (atom.pdb_atom_name == "C" && other.pdb_atom_name == "O")`
**Right:** `if (i == res.C && j == res.O)` (backbone index cache)

The backbone index cache (res.N, res.CA, res.C, res.O, etc.) was
populated at the PDB loading boundary from string names. After that,
these are atom INDICES. Using them is using typed properties.

**For sidechain C=O:** Use element pair (C + O) and bond length
(< 1.35A indicates double bond character). This is a physical
criterion, not a name criterion.

### 3. Per-atom backbone lookup via index cache

The old code checked backbone membership by scanning ALL residues:
```cpp
for (const auto& res : residues_)
    if (i == res.N || i == res.CA || ...) a_bb = true;
```

This is O(N * R) per bond pair. The fix: pre-build a `vector<bool>
is_backbone` from the backbone index cache once, then use O(1) lookup.

### 4. PDB loading boundary comments

Any function that uses string comparisons for atom/residue identity
MUST be clearly commented as "PDB LOADING BOUNDARY" with an
explanation of what is being translated from strings to typed objects.
After the boundary, no string work.

Grep for: `PDB LOADING BOUNDARY` to find all boundary code.

---

### 13. Protonation flows through Protein

Protonation determines which atoms exist, which charges apply, which
ring types are assigned. Changing protonation means a new Protein
(copy-and-modify), not mutation of an existing one. ConformationResults
attached to a conformation depend on the protonation that was in effect
when they were computed. There is no "re-protonate" — there is "copy
with new protonation and re-run the pipeline."

### 14. Typed charge sources, not string-dispatched

ChargeAssignmentResult takes a typed ChargeSource, not a string path.
Different force fields are different types (ParamFileChargeSource,
PrmtopChargeSource, GmxTprChargeSource), not the same function with a
different string parameter. The ForceField enum records provenance in
ProteinBuildContext; the code path is determined by type.

### 15. Variant index contract

The position of each ProtonationVariant in AminoAcidType::variants is
load-bearing. Reordering silently breaks every protonation assignment.
ValidateVariantIndices() runs at test startup and aborts on mismatch.
The contract is documented in AminoAcidType.h with explicit indices.

### 16. Element verification at loading boundaries

When loading data from external tools (ORCA NMR output, prmtop),
verify element matching between the loaded data and the Protein's
atoms. If the ordering is wrong, refuse to load — silent misalignment
means every value goes to the wrong atom.

### 17. Each protein is complete on its own

WT and mutant are separate Proteins, each with their own conformations
and results. The library does not know they are a pair.

The mutation delta is a ConformationResult on the WT: it records what
we learned about THIS conformation by comparing it to a mutant. The
mutant is the instrument, not the owner. MutationDeltaResult attaches
to the WT conformation and stores the deltas internally.

```cpp
auto delta = MutationDeltaResult::Compute(wt_conf, mutant_conf);
wt_conf.AttachResult(std::move(delta));
// Now: wt_conf.Result<MutationDeltaResult>().DeltaT0At(atom_idx)
```

---

### 18. FinalizeConstruction: one call, correct order

Every loader must call `protein->FinalizeConstruction(positions)`
after adding all atoms and residues. This performs:
1. CacheResidueBackboneIndices (needs residues)
2. DetectAromaticRings (needs backbone cache for residue types)
3. DetectCovalentBonds (needs rings for aromatic bond classification)

Order matters. The OrcaRunLoader originally missed this entirely,
producing proteins with zero rings and zero bonds. Both PdbFileReader
and OrcaRunLoader now call FinalizeConstruction. Any future loader
must do the same.

### 19. The full McConnell tensor is NOT Δχ × K

The commonly implemented McConnell formula stores only the symmetric
traceless dipolar kernel K_ab = (3 d̂_a d̂_b - δ_ab) / r³. This is
INCOMPLETE — its trace is zero, so it predicts no isotropic shift.

The full McConnell shielding tensor, derived from the magnetic dipole
interaction (see GEOMETRIC_KERNEL_CATALOGUE.md), is:

```
M_ab = [9 cosθ d̂_a b̂_b - 3 b̂_a b̂_b - (3 d̂_a d̂_b - δ_ab)] / r³
```

This tensor is asymmetric (T1 ≠ 0) and non-traceless (T0 = McConnell
scalar f = (3cos²θ - 1)/r³). The same derivation applies to ring
susceptibility (b̂ → ring normal) and H-bond dipolar (b̂ → D-H...A
direction).

Store BOTH: the symmetric K (for features matching the old code) and
the full M (for shielding contribution and residual subtraction).

### 20. Calculators compute geometric kernels, not parameterised output

Classical calculators store the GEOMETRIC KERNEL (the shape of the
effect as a function of position), not the parameterised shielding.
The kernel is what the physics produces. The parameters (Δχ, I, A, γ)
are weights that the model learns. Literature values are starting
points, not gospel — the context of the original measurement rarely
matches what we're computing.

The calculator's job: compute M_ab / r³ accurately with correct
units and numerical stability. The parameter correction model's job:
learn the weights that make it match DFT.

### 21. Calculator validation: the full analytical process

A calculator is not validated by compiling, passing toy tests, or
producing non-zero output. It is validated by demonstrating that it
produces physically meaningful angular structure on real proteins
with DFT ground truth. The process, in order:

**Step 1: Formula verification on a real protein.**
Not a toy protein. Toy proteins cannot participate in the full
pipeline (no DFT, no real charges, no rings, no mutation delta).
Use the ORCA test protein (A0A7C5FAR6) or 1UBQ. Verify mathematical
identities (e.g., T0 = f for McConnell/RingSusceptibility, EFG
tracelessness for Coulomb) at machine precision.

**Step 2: Batch run on the full 723-pair working set.**
Every clean WT+ALA pair. Zero failures. This tests whether the
calculator handles all protein sizes, ring counts, charge states,
and geometries without NaN, Inf, or crash.

**Step 3: Physical magnitude checks.**
The numbers must be interpretable. E-field magnitudes should match
published ranges (Case 1995: 1–10 V/A at backbone amide H). Ring
susceptibility T0 should scale with 1/r³ and ring count. If a
number cannot be explained from the physics, investigate before
proceeding.

**Step 4: Per-source-type breakdown.**
For Coulomb: backbone vs sidechain vs aromatic decomposition.
For ring susceptibility: per-ring-type table (PHE, TYR, TRP6,
TRP5, TRP9, HIE). Six-membered rings should behave similarly.
Five-membered rings should show higher mean signal (atoms can
get closer to smaller ring centers). Fused rings (TRP9) should
show the largest signal. If the pattern is wrong, the calculator
or the ring geometry is wrong.

**Step 5: DFT proximity analysis.**
The calculator's signal should be STRONGER near mutation sites
(where the physical source was removed in the ALA mutant) than
far away. This tests whether the spatial decay matches where DFT
shielding changes actually occur. For ring susceptibility on 720
proteins: 7x near/far ratio at 8A test threshold. If this ratio
is ~1.0, the calculator produces signal that does not correlate
with where the physics happens.

**Step 6: T2 independence between calculators.**
The thesis goal is to derive underlying physics from the T2
angular residual. This requires that different calculators produce
DIFFERENT angular patterns at the same atoms. Measure the 5D
cosine similarity between each pair of calculators' T2 components
across all atoms. Random 5D vectors give mean |cos| ≈ 0.36.
Values near 1.0 mean the two calculators are redundant — the
model cannot distinguish them. Results on 720 proteins, 386K atoms:

  McConnell vs Coulomb EFG:   0.47 (partially correlated — both
    sum over point sources at every atom/bond)
  McConnell vs Ring Chi:      0.41 (near-random — bonds everywhere,
    rings sparse with different orientation)
  Coulomb EFG vs Ring Chi:    0.40 (near-random — different source
    geometries produce different angular patterns)

All three pairs are well below 0.9. Each calculator gives the
model independent angular information. If a new calculator's T2
is parallel to an existing one (|cos| > 0.9), it adds no new
information and should be questioned.

**Step 7: Unit consistency.**
Every field stored on ConformationAtom must have documented units
in the header comment. Coulomb E-field is in V/A (raw sum × ke =
14.3996). APBS E-field is in V/A (native kT/(e·A) × kT/e =
0.025693 V). The Coulomb constant is real physics (the strength
of the electromagnetic interaction), not a model scaling factor.
Do not remove it; do not add arbitrary scaling. If two fields are
compared (e.g., solvent = APBS − vacuum), they must be in the
same units.

**No fallbacks.** If an external tool (APBS, MOPAC, DSSP) fails,
the result is absent, not faked. Substituting a different physical
quantity (vacuum Coulomb for solvated PB) silently corrupts
downstream analysis. Return nullptr. The pipeline checks
HasResult<T>().

### 22. Ring current sign convention: sigma = I * G

The shielding tensor sigma_ab = -dB_a^sec / dB_{0,b} has a minus
sign from the definition. Ring current geometric kernels absorb this:

```
BS:  G_ab = -n_b * B_a * PPM_FACTOR
HM:  G_ab = -n_b * (H . n)_a
```

With this convention, sigma = I * G gives the correct physical sign
using literature ring current intensities (I < 0 for diamagnetic).
Verified analytically: I=-12, atom 3A above PHE -> sigma = +1.40
ppm (shielded). In-plane at 5A -> sigma = -0.16 ppm (deshielded).
Magnitudes match Case (1995).

The sign was initially wrong (G = +n⊗B, giving sigma < 0 above
ring with diamagnetic I). Caught by analytical test, not by
compilation or unit tests. **Always verify sign with a known
physical scenario before declaring a calculator correct.**

### 23. BS-HM T2 redundancy — measure, don't assume

BiotSavart and HaighMallion model the same physics (pi-electron
circulation) with different mathematical approximations. Before
batch validation, it was an open question whether they produce
independent angular information.

Result: T2 cosine similarity = 0.999 across 279K atoms on 465
proteins. They are effectively parallel — the model gains no new
T2 information from having both. This is a FINDING, not a bug.

Both are stored because:
(a) The finding itself is a thesis result.
(b) The raw HM integral (pure T2) is a different geometric object
    from BS G (rank-1), even though the full HM kernel G is parallel.
(c) The TRP fused ring superposition is identical for both: BS ratio
    1.000, HM ratio 1.000. Both are perfectly additive after the
    RingBondedExclusionFilter excludes shared-edge atoms. (The
    previously reported HM 13% excess at ratio 1.127 was an artifact
    from evaluating the surface integral at ring atoms inside the
    source distribution — corrected 2026-04-02.)

### 24. Fused ring representation: document, don't hide

TRP produces three rings (TRP5, TRP6, TRP9). Their intensities
sum exactly: I(TRP9) = I(TRP5) + I(TRP6) = -19.2. For both BS
and HM, the geometry is also additive (ratio 1.000) after
RingBondedExclusionFilter excludes ring vertices and their bonded
neighbours from through-space evaluation.

All three are stored as independent per-type features. The model
sees the component decomposition (TRP5+TRP6) and the whole-system
representation (TRP9).

---

## Numerical Stability

Structurally correct physics with numerically naive implementation
produces garbage at critical atom positions. These rules are
non-negotiable.

### PhysicalConstants.h: the single home for constants and units

All universal constants, unit conversions, and numerical thresholds
live in PhysicalConstants.h. Only universal constants go here — no
model-specific parameters (those live on ring type classes and
calculator parameter structs).

What is there now:
- SI electromagnetic: mu_0, Biot-Savart prefactor
- Unit conversions: A→m, nA→A, PPM_FACTOR, COULOMB_KE (14.3996 V·A),
  KT_OVER_E_298K (0.025693 V for APBS)
- Numerical thresholds: MIN_DISTANCE (0.1 A), NEAR_ZERO_NORM,
  NEAR_ZERO_FIELD, APBS_SANITY_LIMIT (100 V/A)
- Spatial shells: RING_COUNT_SHELL 3/5/8/12 A, RING_CALC_CUTOFF 15 A
- H-bond: HBOND_COUNT_RADIUS (3.5 A), HBOND_MAX_DIST
- Sequence: SEQUENTIAL_EXCLUSION_THRESHOLD (2)

Calculator-specific thresholds (e.g., dispersion R_SWITCH/R_CUT) stay
in the calculator .cpp with their physics documentation. They are not
global because their values depend on the specific switching function.

If you need a new constant, put it here with a comment citing the
source. If you need a new threshold, document why that value.

### Singularity guard

All 1/r^n terms use MIN_DISTANCE = 0.1 A (PhysicalConstants.h).
Below this, the kernel value is numerically meaningless. The
MinDistanceFilter in KernelFilterSet enforces this before any
kernel evaluation begins.

### Filter before computing, not after

Every calculator holds a KernelFilterSet. Call AcceptAll() BEFORE
computing the kernel. Do not compute the kernel and then discard
the result — the computation itself may overflow. Five concrete
filters exist (KernelEvaluationFilter.h). Use them.

### Tracelessness after accumulation

Floating-point accumulation across many source terms breaks the
tracelessness of tensors that are analytically traceless (Coulomb
EFG, Pi-quadrupole EFG). Apply traceless projection after
summation: `V -= (V.trace() / 3.0) * Mat3::Identity()`. Do this
for any tensor where Gauss's law or symmetry guarantees zero trace.

### Unit chains: state them explicitly

Every calculator must document its unit chain in comments. The
chain from raw sum to stored field must be traceable:

    Coulomb: q in e, r in A → raw sum in e/A² → × ke (14.3996 V·A) → V/A
    APBS: native kT/(e·A) → × kT/e (0.025693 V) → V/A
    Biot-Savart: A → m → SI Tesla → × PPM_FACTOR → × I → ppm

If two fields are compared (e.g., solvent = APBS − vacuum), they
must be in the same units. The Coulomb constant k_e = 14.3996 V·A
is real physics. Do not remove it. Do not add arbitrary scaling.

### Sign convention: verify analytically

The ring current sign convention is G_ab = -n_b * B_a * PPM_FACTOR.
The minus sign comes from sigma_ab = -dB_a^sec / dB_{0,b}. With
this convention, sigma = I * G gives the correct physical sign
using literature intensities (I < 0 for diamagnetic rings).

The sign was initially wrong. It was caught by an analytical test
(I=-12, atom 3A above PHE → sigma = +1.40 ppm, shielded), not by
compilation or unit tests. Always verify sign with a known physical
scenario before declaring a calculator correct.

### Near-field stability per kernel

- **Dipolar (McConnell, RingSusceptibility, HBond):** 1/r³ leading
  term. DipolarNearFieldFilter with source_extent = bond length or
  ring diameter or N...O distance.
- **Biot-Savart:** Wire segment divergence near endpoints. Thresholds
  lenA < 1e-25 and crossSq < 1e-70 (SI metres). Work in SI to
  avoid Angstrom-scale underflow. Johnson-Bovey loops at ±d nearly
  cancel in the ring plane (z ≈ 0) — precision loss in B_z.
- **Haigh-Mallion:** Adaptive subdivision at 2.0 A (level 1) and
  1.0 A (level 2) for atoms near the ring face. 7-point Gaussian
  quadrature (Stroud T2:5-1). Two subdivision levels max.
- **Pi-quadrupole:** 1/r⁹ leading term — steepest divergence.
  DipolarNearFieldFilter plus RingBondedExclusionFilter required.
  Without: max |T2| = 7.39 A⁻⁵. With: 0.66 A⁻⁵.
- **Coulomb EFG:** Traceless per term but accumulation breaks it.
  Apply traceless projection after total. Clamp E magnitude to
  100 V/A for rare pathological geometries.
- **Dispersion:** CHARMM switching function (Brooks et al. 1983)
  tapers smoothly between R_switch=4.3 A and R_cut=5.0 A.
  Prevents feature discontinuities across MD ensemble frames.
  Formula: S(r) = (Rc²-r²)²(Rc²+2r²-3Rs²) / (Rc²-Rs²)³.
- **H-bond:** N...O distance ~2.8 A means midpoint is 1.4 A from
  endpoints. Without DipolarNearFieldFilter: max |T2| = 1908 A⁻³.
  With: 0.78 A⁻³. The filter is not optional.

### NaN/Inf: absent, not faked

Sanitise NaN/Inf from near-zero denominators. If an external tool
(APBS, MOPAC, DSSP) fails, the result is absent (nullptr), not
faked. Substituting a different physical quantity (vacuum Coulomb
for solvated PB) silently corrupts downstream analysis. Return
nullptr. The pipeline checks HasResult<T>().

### When implementing a calculator

You must record GeometryChoices via GeometryChoiceBuilder during
Compute(). Every inclusion, exclusion, and triggered event gets a
Record() call inside a lambda. Attach the entities involved (atoms,
rings, bonds) with their roles and outcomes. Add named numbers
(distance, intensity, threshold) with units. If a filter rejects,
call LastRejectorName() and record it. See any existing calculator's
Compute() for the pattern. See GeometryChoice.h for the API.

---

### Data-driven accumulator columns (AllWelfords pattern)

GromacsProteinAtom has ~45 Welford accumulators. WriteCatalog (CSV)
and WriteH5 both need the column list. A switch statement indexed by
column number will go wrong the moment someone adds a column.

The fix: `AllWelfords()` returns a `vector<NamedWelford>` with
`{name, pointer-to-Welford}` pairs. WriteCatalog and WriteH5 iterate
this — no switch, no manual indexing, no parallel names vector.
Adding a new Welford = add the field + one line in AllWelfords().
The CSV header, H5 rollup, and SDK column names all derive from it.

This is not ideal — the accumulation in AccumulateFrame is still
manual, and the pointers-into-self pattern prevents moving
GromacsProteinAtom. But it eliminates the 45-case switch statement
and the duplicate column list, which are the two places where silent
corruption from misindexing would be hardest to detect. The
AccumulateFrame code is write-once (change it when you add a new
calculator); the column serialisation runs on every protein and must
not drift.

**Pattern:** When you have N named fields that need serialisation,
return them as a `{name, pointer}` vector from one method. Iterate
that method in every serialisation path. Never maintain a parallel
list of names and a switch statement indexed by position.

---

## This Document Is Living

When something breaks or a new pattern emerges, it goes here.
