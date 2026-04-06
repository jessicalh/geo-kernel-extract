# Library Interface for the Viewer

What the library exposes. What the viewer consumes. No adapters.

---

## The traversal as it stands today

The object model already supports clean traversal. Here's how it reads:

```cpp
// === Load ===
auto lr = nmr::LoadProtein("/path/to/protein.pdb");
auto& protein = *lr.protein;
auto& conf = protein.CrystalConf();

// === Run classical physics (dependency order) ===
conf.AttachResult(GeometryResult::Compute(conf));
conf.AttachResult(SpatialIndexResult::Compute(conf));
conf.AttachResult(EnrichmentResult::Compute(conf));
conf.AttachResult(DsspResult::Compute(conf));
conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));
// ... all 8 calculators in any order (deps satisfied) ...
conf.AttachResult(BiotSavartResult::Compute(conf));
conf.AttachResult(HaighMallionResult::Compute(conf));
conf.AttachResult(McConnellResult::Compute(conf));
conf.AttachResult(CoulombResult::Compute(conf));
conf.AttachResult(RingSusceptibilityResult::Compute(conf));
conf.AttachResult(PiQuadrupoleResult::Compute(conf));
conf.AttachResult(DispersionResult::Compute(conf));
conf.AttachResult(HBondResult::Compute(conf));

// === Traverse the structure ===

// Protein level: identity (does not change with conformation)
for (size_t i = 0; i < protein.AtomCount(); i++) {
    const auto& id   = protein.AtomAt(i);  // element, name, bonds
    const auto& atom  = conf.AtomAt(i);     // position + all computed fields
    const auto& res   = protein.ResidueAt(id.ResidueIndex());

    // Everything is right here:
    Vec3 pos          = atom.Position();
    Element elem      = id.element;
    AtomRole role     = atom.role;           // set by EnrichmentResult
    bool backbone     = atom.is_backbone;

    // Per-calculator shielding (SphericalTensor, has .T0 and .T2)
    auto bs_sigma     = atom.bs_shielding_contribution;
    auto hm_sigma     = atom.hm_shielding_contribution;
    auto mc_sigma     = atom.mc_shielding_contribution;
    auto coulomb_sig  = atom.coulomb_shielding_contribution;
    auto hbond_sig    = atom.hbond_shielding_contribution;
    auto pq_sigma     = atom.piquad_shielding_contribution;
    auto chi_sigma    = atom.ringchi_shielding_contribution;
    auto disp_sigma   = atom.disp_shielding_contribution;

    // Vector fields
    Vec3 B            = atom.total_B_field;
    Vec3 E            = atom.coulomb_E_total;

    // Full tensors (for glyphs)
    Mat3 G            = atom.total_G_tensor;

    // Ring neighbourhood (per-ring detail)
    for (const auto& rn : atom.ring_neighbours) {
        rn.G_tensor;          // BS kernel from this ring
        rn.hm_tensor;         // HM kernel from this ring
        rn.chi_tensor;        // susceptibility kernel from this ring
        rn.quad_tensor;       // PQ kernel from this ring
        rn.disp_tensor;       // dispersion kernel from this ring
        rn.distance_to_center;
        rn.rho; rn.z; rn.theta;  // cylindrical coords
    }

    // Bond neighbourhood (per-bond detail)
    for (const auto& bn : atom.bond_neighbours) {
        bn.dipolar_tensor;    // McConnell kernel from this bond
        bn.mcconnell_scalar;
        bn.bond_category;
        bn.distance_to_midpoint;
    }
}

// Rings
for (size_t i = 0; i < protein.RingCount(); i++) {
    const auto& ring = protein.RingAt(i);
    const auto& geo  = conf.ring_geometries[i];
    // ring.TypeIndex(), ring.Intensity(), ring.VertexAtomIndices()
    // geo.center, geo.normal, geo.radius
}

// Bonds
for (size_t i = 0; i < protein.BondCount(); i++) {
    const auto& bond = protein.BondAt(i);
    // bond.atom_index_a, bond.atom_index_b, bond.category
    // conf.bond_midpoints[i], conf.bond_directions[i], conf.bond_lengths[i]
}

// === Check before you ask ===
if (conf.HasResult<BiotSavartResult>()) {
    // safe to read bs_shielding_contribution
}
// If you call conf.Result<T>() and T isn't attached:
// FATAL message naming what's missing and what IS attached, then abort.
// Not a silent wrong answer. Not an exception to catch. A programming error.
```

That's the traversal. It already works. The viewer doesn't need an
adapter — it needs this example and the field name table from
VIEWER_MIGRATION.md.

---

## What's missing: grid evaluation for contour plots

The traversal above gives you values **at atom positions**. Contour
plots and isosurfaces need values **on a 3D grid**. The old viewer
called `calculator.ShieldingTensorAtPoint(pt, protein, conf)` in a
loop. The new architecture doesn't have that — results compute once
for all atoms.

### Proposal: SampleAt() on results that support it

Each calculator that can meaningfully evaluate at an arbitrary point
gets a const method:

```cpp
// On BiotSavartResult:
SphericalTensor SampleShieldingAt(Vec3 point) const;
Vec3            SampleBFieldAt(Vec3 point) const;

// On HaighMallionResult:
SphericalTensor SampleShieldingAt(Vec3 point) const;

// On McConnellResult:
SphericalTensor SampleShieldingAt(Vec3 point) const;

// On RingSusceptibilityResult:
SphericalTensor SampleShieldingAt(Vec3 point) const;

// On PiQuadrupoleResult:
SphericalTensor SampleShieldingAt(Vec3 point) const;

// On DispersionResult:
SphericalTensor SampleShieldingAt(Vec3 point) const;

// On CoulombResult:
Vec3            SampleEFieldAt(Vec3 point) const;
Mat3            SampleEfgAt(Vec3 point) const;

// On HBondResult:
SphericalTensor SampleShieldingAt(Vec3 point) const;
```

**Why this works**: the result already holds everything it needs.
BiotSavartResult has access to the protein's rings and the
conformation's ring geometries. It can evaluate the kernel at any
point using the same code that computed atom values. The `const`
means it doesn't mutate state — it's a pure query on an already-
computed result.

**Why not a separate FieldSampler class**: that would just be an
adapter holding a reference to the result. No value added. The result
IS the sampler — it has the state.

**Why not a free function**: it would need the protein, conformation,
and ring geometries as arguments. The result already has those via
its conformation back-pointer. Less to get wrong.

### What the viewer does with SampleAt

```cpp
// Build a scalar field grid for isosurfaces
auto& bs = conf.Result<BiotSavartResult>();

int G = 20;
double extent = 7.0;  // Angstroms
for (iz...) for (iy...) for (ix...) {
    Vec3 pt = origin + Vec3(ix, iy, iz) * spacing;
    double T0 = bs.SampleShieldingAt(pt).T0;
    grid[ix + iy*G + iz*G*G] = T0;
}
// Hand grid to VTK for marching cubes. Done.
```

```cpp
// Build a B-field grid for streamlines
for (...) {
    Vec3 B = bs.SampleBFieldAt(pt);
    // Hand to VTK streamline tracer
}
```

```cpp
// Build a tensor glyph field for McConnell
auto& mc = conf.Result<McConnellResult>();
for (sample points...) {
    SphericalTensor st = mc.SampleShieldingAt(pt);
    // st.T0 ≈ 0 (traceless), st.T2 is the angular pattern
    // Render SH glyph at pt with lobes from T2
}
```

---

## How to show geometric kernels: the visualization modes

### Mode 1: Scalar isosurface (T0)

**What it shows**: surfaces of constant isotropic shielding.
The classic "shielding cone" above a ring.

**Works for**: BS, HM, RingSusceptibility, HBond, Dispersion
(any calculator with nonzero T0).

**Doesn't work for**: McConnell, PiQuadrupole (pure T2, T0 ≈ 0).

**How**: SampleShieldingAt() on a 3D grid → VTK isosurface.
Two surfaces: +threshold (shielded, blue) and -threshold
(deshielded, red). Layered thresholds show field strength falloff.

### Mode 2: Tensor glyph field (T2)

**What it shows**: the angular pattern of anisotropic shielding at
sample points. Lobed shapes from spherical harmonic expansion of T2.

**Works for**: ALL calculators. This is the universal visualization.
McConnell and PQ are ONLY visible this way.

**How**: SampleShieldingAt() at sample points → SphericalTensor →
render f(theta,phi) = sum_m T2_m * Y_2^m as a surface glyph.
Red lobes = deshielded direction, blue lobes = shielded direction.

### Mode 3: Vector field (B, E)

**What it shows**: field direction and magnitude.

**B-field**: BS magnetic field streamlines (butterfly), arrows.
**E-field**: Coulomb electric field arrows.

**How**: SampleBFieldAt() / SampleEFieldAt() on a grid → VTK
streamlines or arrow glyphs.

### Mode 4: Per-atom bubbles

**What it shows**: accumulated calculator contributions at each atom.

**How**: read ConformationAtom fields directly. Color sphere by T0.
Size or opacity by magnitude. Per-calculator toggles (the 8 checkboxes).

No grid sampling needed — this uses atom positions only.

### Mode 5: Per-atom tensor glyphs

**What it shows**: the full tensor at each atom as a glyph.

**How**: read SphericalTensor from ConformationAtom. Render glyph.
This is what the old viewer's TensorGlyph and EllipsoidGlyph do.

### Mode 6: Single-source isolation

**What it shows**: the kernel from ONE ring or ONE bond, in isolation.
User clicks a ring → show that ring's shielding cone. User clicks a
bond → show that bond's McConnell pattern.

**How**: per-ring data is in RingNeighbourhood (already computed per
atom). For grid sampling, SampleAt could take an optional source
filter: "evaluate only ring #3." Or simpler: the viewer can
re-evaluate the kernel formula directly from ring geometry —
it's one formula, documented in GEOMETRIC_KERNEL_CATALOGUE.md.

---

## Flow of control: what the viewer calls

```
1. LoadProtein(pdb_path) → unique_ptr<Protein>
   One call. Protein is fully constructed (atoms, residues, bonds, rings).
   One CrystalConformation with const positions.

2. Run classical physics:
   RunClassicalCalculators(conf, charge_source)
   ← this is the ONE new function we'd add

3. Optionally load DFT:
   auto orca_lr = LoadOrcaRun(files);
   conf.AttachResult(OrcaShieldingResult::Compute(conf, orca_path));

4. Optionally load MD frames (stub):
   // protein.AddMDFrame(positions, metadata);
   // RunClassicalCalculators(md_conf, charge_source);
   // ... repeat per frame

5. Query: traverse atoms, sample grids, render.
```

### RunClassicalCalculators

One function that enforces dependency order. Not an architecture —
a convenience:

```cpp
// In a new header, e.g. Pipeline.h
namespace nmr {

struct PipelineOptions {
    ChargeSource* charge_source = nullptr;  // null = skip Coulomb
    bool skip_dssp = false;                 // true = skip HBond too
};

// Attaches all tier-0 + tier-1 results in correct order.
// Returns list of what was attached (for logging).
// If a dependency fails, stops and reports what failed.
std::vector<std::string> RunClassicalCalculators(
    ProteinConformation& conf,
    const PipelineOptions& opts = {});

}  // namespace nmr
```

Implementation:
```cpp
std::vector<std::string> RunClassicalCalculators(
        ProteinConformation& conf,
        const PipelineOptions& opts) {
    std::vector<std::string> attached;

    // Tier 0
    conf.AttachResult(GeometryResult::Compute(conf));
    attached.push_back("GeometryResult");
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    attached.push_back("SpatialIndexResult");
    conf.AttachResult(EnrichmentResult::Compute(conf));
    attached.push_back("EnrichmentResult");
    if (!opts.skip_dssp) {
        conf.AttachResult(DsspResult::Compute(conf));
        attached.push_back("DsspResult");
    }
    if (opts.charge_source) {
        conf.AttachResult(ChargeAssignmentResult::Compute(conf, *opts.charge_source));
        attached.push_back("ChargeAssignmentResult");
    }

    // Tier 1: all 8 calculators
    conf.AttachResult(BiotSavartResult::Compute(conf));       attached.push_back("BiotSavart");
    conf.AttachResult(HaighMallionResult::Compute(conf));     attached.push_back("HaighMallion");
    conf.AttachResult(McConnellResult::Compute(conf));        attached.push_back("McConnell");
    conf.AttachResult(RingSusceptibilityResult::Compute(conf)); attached.push_back("RingSusceptibility");
    conf.AttachResult(PiQuadrupoleResult::Compute(conf));     attached.push_back("PiQuadrupole");
    conf.AttachResult(DispersionResult::Compute(conf));       attached.push_back("Dispersion");
    if (opts.charge_source) {
        conf.AttachResult(CoulombResult::Compute(conf));      attached.push_back("Coulomb");
    }
    if (!opts.skip_dssp) {
        conf.AttachResult(HBondResult::Compute(conf));        attached.push_back("HBond");
    }

    return attached;
}
```

That's ~30 lines. Not an abstraction — just the dependency order
written down once instead of copied everywhere.

---

## What we DON'T add

- No ViewerProtein adapter wrapping Protein
- No ViewerConformation wrapping ProteinConformation
- No string-based result lookup ("biot_savart")
- No observer/callback pattern for result updates
- No generic "field sampler" factory
- No runtime plugin system for calculators

The viewer includes the library headers it needs, calls LoadProtein,
calls RunClassicalCalculators, traverses ConformationAtom fields,
calls SampleAt for grids. That's the interface.
