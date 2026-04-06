# Patterns: Implementation Guide for All Agents

This document is given to every implementation agent. It describes how
to write code in this system. Follow it.

This is a physics model, not a software system. If the code reads like
infrastructure, it is wrong. If it reads like equations with data flow,
it is right. You are writing a model that a physicist would recognise,
not an enterprise application. Maximise physics depth. Keep the code
simple, readable, and traceable.

You as an AI read code and act upon the name. Like ancient magic the
name is the thing. Humans have similar wiring. If the name of a model
is "HIE" then the thing is a meaningless code. That is why we build
models in programming languages. It is nice they do things quickly but
it is more important they describe the world. Otherwise what they do
quickly has no meaning. HieImidazoleRing with NitrogenCount() == 2
and Aromaticity() == Weak describes the world. "HIE" is a string that
only means something if you already know what it means.

---

## The System in One Paragraph

This is a one-way analysis system. It loads a protein, creates typed
conformations, and performs successive calculations on them. Each
calculation is a ConformationResult that attaches to a conformation
and writes to typed fields on ConformationAtom. A trained e3nn model
(the ParameterCorrectionResult) provides corrected parameters to the
classical calculators. The L=2 residual after all calculators shows
where classical angular physics breaks down. Nothing flows backward.
Everything accumulates forward.

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
for features and reuse) AND shielding contribution (ppm SphericalTensor,
for residual tracking). At the end of Attach(), if the
ParameterCorrectionResult is present, the calculator calls
SubtractCalculatorContribution with its shielding output. The calculator
works correctly without the correction model — the residual update is
a diagnostic, not a dependency.

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

Calculator code is ring-type-agnostic. Use `ring.intensity()`, not a
switch statement. When the ParameterCorrectionResult provides a
corrected intensity, it is indexed by `ring.TypeIndexAsInt()`. The
calculator does not need to know which ring type it is processing.

```cpp
double I = ring.intensity();         // -12.0 for PHE, -5.16 for HIS
double d = ring.jb_lobe_offset();    // 0.64 for PHE, 0.50 for HIS
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
call `ring.nitrogen_count()` or check `identity.role`. If you need to
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
`ring.intensity()` is a virtual method. No tag dispatch. No policy
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
There is nothing to unwind.

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

The 93 tuneable parameters from the ParameterCorrectionResult are the
exception — they are genuinely tuneable and are passed via typed
CorrectedParameters structs.

### Over-decomposing into tiny functions

The Johnson-Bovey B-field calculation is ONE physics operation. It
takes a ring and an atom position and produces a B-field, a geometric
kernel, and cylindrical coordinates. When you decompose this into 15
five-line functions, the physics narrative is destroyed. You cannot
read the computation linearly. One function per physics concept. The
Biot-Savart inner loop can be 30-40 lines of Eigen math with the
equation in a comment. That IS the appropriate size.

---

## T2 Completeness

Every calculator must produce full tensor output at all irrep levels
(L=0, L=1, L=2). T0-only results are incomplete and will be rejected.
The T2 angular structure is a first-class thesis result.

### Per-calculator T2 explained

Each calculator stores how much T2 residual it reduced. After all
calculators attach, you can compare: does Biot-Savart or Haigh-Mallion
better explain the T2 pattern near PHE rings? Without per-calculator
T2 explained fields, you cannot answer this after the pipeline completes.

### BS-HM T2 disagreement

Biot-Savart (line integral, rank-1) and Haigh-Mallion (surface integral,
rank-2) make opposing T2 predictions at the same geometry. The
disagreement vector at each atom shows where the two models diverge
angularly. This is stored on ConformationAtom and is a feature and
a diagnostic.

### Missing T2 features to include

Per-type PiQuadrupole T2 (8 L2 features), per-type RingSusceptibility
T2 (8 L2 features). MOPAC orbital populations (s_pop, p_pop) provide
per-atom electronic structure but are L=0 features, not L=2 tensors.

### Additional tuneable parameters

McConnell CO midpoint shift split by H-bond state (+1 parameter).
H-bond angular exponent (+1 parameter). Total 93 becomes 95.

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

**CMake:** Modular. Each result type is its own target. Adding a new
result = new subdirectory + one add_subdirectory line.

---

## What To Read For Your Pass

- **Every agent:** CONSTITUTION.md (principles), OBJECT_MODEL.md
  (concrete types), this file (patterns).
- **Your pass specification** from LAYER0_PLAN.md.
- **Previous pass code** in nmr-shielding/src/.
- Do NOT read FEEDBACK.md (historical, mostly resolved).
- Do NOT read DESIGN_REVIEW_PRELIM.md (historical).

---

## Lessons from Pass 0 Correction

### 1. Back-pointer fixup on move (BUG: dangling Protein*)

Protein owns conformations via `vector<unique_ptr<ProteinConformation>>`.
Each conformation holds a raw `const Protein*` back-pointer set at
construction time. When a Protein is moved (e.g. `make_unique<Protein>
(std::move(r.protein))`), the default move constructor transfers the
conformation vector but does NOT update the back-pointers -- they
still point to the old (now-destroyed) Protein location.

The fix: explicit Protein move constructor that calls
`conf->FixProteinBackPointer(this)` on every owned conformation.

**Pattern:** Any class that holds children with raw back-pointers to
the parent MUST have an explicit move constructor that fixes those
pointers. The compiler-generated move is insufficient.

**Detection:** Crashes (segfault, bad_alloc) that appear only when
the owning object is moved (e.g. `make_unique<T>(std::move(value))`
in test fixtures) but not when accessed in-place.

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

## Lessons from Protonation Pipeline (2026-04-01)

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

## Lessons from Calculator Implementation (2026-04-02)

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

## This Document Is Living

After every correction pass, PATTERNS.md is updated with lessons
learned. If a correction agent finds a new anti-pattern or a pattern
that worked well, it goes here. The document gets better with each pass.
