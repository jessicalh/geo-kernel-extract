# Azimuthal Angle φ — Specification

## What and why

The ring_contributions sparse table stores cylindrical coordinates
(rho, z, theta) for each (atom, ring) pair.  theta is the polar
angle — the cone angle from the ring normal.  There is no azimuthal
angle φ — the in-plane direction is lost.

For symmetric rings (PHE benzene, 6-fold), φ averages out across
many atoms.  For asymmetric rings (HIE imidazole), the nitrogen
positions break rotational symmetry.  An atom at the same (rho, z)
but different φ relative to the C-N axis sees a different local
field from the asymmetric π-electron distribution.

The BS/HM kernels compute the full tensor correctly (the integration
goes around the actual ring vertices), so the T2 components *encode*
φ implicitly.  But the model has no scalar feature telling it where
the atom sits in the ring plane.  Without that, the MLP cannot
learn that atoms on the nitrogen side of HIE behave differently
from atoms on the carbon side.

Secondary analysis shows HIE self-fit R² = 0.065 (lowest of all
ring types), while EFG dominates the hie_only stratum (R²_efg =
0.969 vs R²_ring = 0.038).  Adding φ as a scalar feature gives
the MLP the context to weight HIE ring kernels differently by
in-plane position — the one geometric degree of freedom that is
currently invisible to the scalar side.

## Definition

For an (atom, ring) pair with ring center C, ring normal n̂, and
ring vertex 0 position V₀:

```
d       = atom_pos - C               (displacement vector)
d_plane = d - (d · n̂) n̂             (projection onto ring plane)
ref     = V₀ - C                     (reference direction: center → vertex 0)
ref_plane = ref - (ref · n̂) n̂       (project onto ring plane — nearly a no-op)

cos φ = (d_plane · ref_plane) / (|d_plane| · |ref_plane|)
sin φ = ((d_plane × ref_plane) · n̂) / (|d_plane| · |ref_plane|)
```

When |d_plane| < ε (atom is on the ring normal axis), cos φ = 1,
sin φ = 0 — the angle is undefined but the value is irrelevant
because rho ≈ 0 means the atom is directly above/below the ring
where all in-plane directions are equivalent.

### Reference direction: vertex 0

Ring::atom_indices[0] is the first vertex in the detection order.
This is consistent across conformations of the same protein (same
topology → same detection → same ordering).  It is NOT consistent
across different ring types — PHE vertex 0 and HIE vertex 0 are
different atoms with different chemical meaning.  This is correct:
φ is a ring-frame quantity, and different ring types have different
frames.  The model learns what φ means for each type through the
ring_type one-hot and per-type kernel decomposition.

### Physical meaning by ring type

- **PHE, TYR** (6-membered): near 6-fold symmetry.  φ carries
  little information — confirmed by PHE having higher ring kernel
  R² without φ than HIE has with everything.
- **TRP_benzene, TRP_pyrrole** (6/5-membered): moderate asymmetry
  from the C-C shared edge.  φ relative to the shared bond encodes
  "near the fused edge" vs "far side."
- **HIS, HID, HIE** (5-membered imidazole): strong asymmetry from
  two nitrogen positions.  φ is most valuable here — it encodes
  the atom's position relative to the nitrogen lone pairs.

## C++ changes

### 1. RingNeighbourhood: add two fields

**File:** `src/ConformationAtom.h`

```cpp
struct RingNeighbourhood {
    // ... existing fields ...
    double cos_phi = 1.0;   // azimuthal: cos of in-plane angle to vertex 0
    double sin_phi = 0.0;   // azimuthal: sin of in-plane angle to vertex 0
};
```

Default (1, 0) corresponds to φ = 0 — atom lies along the vertex 0
direction.  This is the correct default for atoms on the ring axis
(rho = 0) where φ is degenerate.

### 2. BiotSavartResult::Compute: compute φ

**File:** `src/BiotSavartResult.cpp`, in the block that creates new
RingNeighbourhood entries (around line 246-264).

After the existing cylindrical coordinate computation:

```cpp
// Cylindrical coordinates in ring frame
double z = d.dot(geom.normal);
Vec3 d_plane = d - z * geom.normal;
double rho = d_plane.norm();
double theta = std::atan2(d_plane.norm(), std::abs(z));
new_rn.z = z;
new_rn.rho = rho;
new_rn.theta = theta;

// Azimuthal angle in ring plane: angle from center→vertex0 direction.
// Encodes in-plane position relative to ring frame — distinguishes
// nitrogen side from carbon side on asymmetric rings (HIE, TRP).
Vec3 ref = geom.vertices[0] - geom.center;
Vec3 ref_plane = ref - ref.dot(geom.normal) * geom.normal;
double ref_norm = ref_plane.norm();
if (rho > 1e-10 && ref_norm > 1e-10) {
    Vec3 d_hat = d_plane / rho;
    Vec3 ref_hat = ref_plane / ref_norm;
    new_rn.cos_phi = d_hat.dot(ref_hat);
    new_rn.sin_phi = d_hat.cross(ref_hat).dot(geom.normal);
}
// else: defaults (1, 0) — atom on axis or degenerate ring
```

This is 6 lines of Eigen math, no new includes, no new dependencies.

### 3. WriteAllFeatures: add two columns

**File:** `src/ConformationResult.cpp`, in the ring_contributions block.

Change column count from 57 to 59.  Insert cos_phi and sin_phi
after gaussian_density (or anywhere — but end is cleanest for
backward compatibility of column indices 0-56):

```cpp
const size_t C = 59;
// ... existing 57 columns unchanged ...
r[57] = rn.cos_phi;
r[58] = rn.sin_phi;
```

Update the column comment at line 56-58 to document columns 57-58.

### 4. ring_geometry.npy: no change

Ring center, normal, radius, and per-vertex positions are already
written.  No change needed to the ring reference table.

## SDK changes

### 1. _catalog.py: update column count

**File:** `python/nmr_extract/_catalog.py`

The RingContributions ArraySpec expects 57 columns.  Change to 59.
Old 57-column extractions are stale — re-extract.  The SDK reads
59 or fails.  That is the trust boundary doing its job.

### 2. _ring.py: add cos_phi, sin_phi properties

**File:** `python/nmr_extract/_ring.py`

```python
@property
def cos_phi(self) -> np.ndarray:
    """Cosine of azimuthal angle in ring plane (relative to vertex 0)."""
    return self._data[..., 57]

@property
def sin_phi(self) -> np.ndarray:
    """Sine of azimuthal angle in ring plane (relative to vertex 0)."""
    return self._data[..., 58]
```

No fallbacks, no conditional column counts.  59 columns is the
format after this change.

### 3. API.md: document new columns

Add columns 57-58 to the ring_contributions table in `python/API.md`.

## Calibration pipeline changes

### 1. scalars.py: add cos_phi, sin_phi to ring_proximity block

**File:** `learn/src/mutation_set/scalars.py`, in `_ring_proximity()`.

Add two columns per ring to the per-ring scalar block:

```python
COLS = 10  # was 8: + cos_phi + sin_phi
# ...
rp[j, base + 8] = atom_rc.cos_phi[ri]
rp[j, base + 9] = atom_rc.sin_phi[ri]
```

This adds 12 scalars (6 rings × 2) to the 247 total → 259.

### 2. kernels.py: optional φ-modulated kernels

A φ-modulated kernel multiplies the ring T2 by cos(φ) or sin(φ):

    K_φ = cos(φ) × BS_T2

This is a genuinely new angular direction — scalar × L=2 = L=2,
equivariant by construction.  It cannot be obtained from existing
kernels because they don't contain φ.  It lets the model learn
that "BS ring current from the nitrogen side of HIE" differs from
"BS ring current from the carbon side."

Whether this helps is empirical.  Add it as an option gated by a
config flag, not as a default.  The scalar cos_phi/sin_phi features
alone may be sufficient — the MLP can learn the modulation from
scalars × existing kernels.  Test scalars first.

## Test strategy

### 1. Unit: verify φ geometry

Add a test in the BS test suite: place an atom at known (rho, z, φ)
relative to a PHE ring.  Verify cos_phi and sin_phi match.  Check:
- Atom on normal axis (rho=0): cos_phi=1, sin_phi=0
- Atom along vertex 0 direction: cos_phi=1, sin_phi=0
- Atom perpendicular to vertex 0 in ring plane: cos_phi=0, sin_phi=±1
- Atom on opposite side: cos_phi=-1, sin_phi≈0

### 2. Extraction: re-extract one protein

Run nmr_extract on A0A7C5FAR6 (the ORCA test protein).  Verify
ring_contributions.npy has 59 columns.  Load via SDK, check
cos_phi² + sin_phi² ≈ 1 for all pairs with rho > 0.1.

### 3. Batch: full 723 re-extraction

Re-extract all 723 proteins.  φ is computed during BS, no MOPAC
needed — but the full pipeline runs MOPAC anyway.  Use a fresh
extraction run name (not GatedCalibration).  Verify no failures,
no NaN in the new columns.  Point calibration.toml at the new
extraction.

### 4. Training: measure impact

Add cos_phi/sin_phi to scalars.py as described.  Retrain on the
full 723 proteins.  Compare R² overall and per-ring-type stratum.
The hypothesis: HIE-stratum R² improves, PHE/TYR barely change.

### 5. UI: update ring_contributions reader

If the viewer reads ring_contributions, update it to expect 59
columns.  The extra columns are scalars, not tensors — no
visualization change needed unless we add a φ-angle display.

## What this does NOT do

- Does not add new tensor kernels (φ-modulated kernels are optional)
- Does not change any existing column index (0-56 unchanged)
- Does not affect the ring current physics — BS/HM already compute
  correctly over asymmetric rings.  This adds *context* for the
  model to use the existing kernels better.

## What this DOES require

- Full re-extraction of all 723 proteins (fresh run name)
- SDK update: _catalog.py expects 59 columns, _ring.py adds accessors
- Re-extract test fixtures for python/tests/
- scalars.py update: ring_proximity block grows by 2 columns per ring
- Update calibration.toml to point at new extraction

## Impact audit: what touches ring_contributions columns

**SDK (must change — owns the column layout):**

| File | What | Action |
|------|------|--------|
| `_catalog.py:58` | Hardcodes `57` in ArraySpec | Change to `59` |
| `_ring.py:3` | Docstring says `(P, 57)` | Update to `(P, 59)` |
| `_ring.py:146` | `gaussian_density` at column `56` | Unchanged — columns 0-56 intact |
| `_ring.py` | Constructor | Add `cos_phi` at `57`, `sin_phi` at `58` |
| `API.md` | Column reference table | Add columns 57-58 |

**Learning code (does NOT break — uses named SDK accessors):**

| File | Accesses | Breaks? |
|------|----------|---------|
| `scalars.py:135-155` | `.mcconnell_factor`, `.exp_decay`, `.z`, `.rho`, `.theta`, `.disp_scalar`, `.disp_contacts` — all named properties | No. Add `.cos_phi`, `.sin_phi` here when ready. |
| `kernels.py:224-259` | `.for_atom()`, `.ring_type`, `.disp_scalar`, `.bs.T2`, `.chi.T2` — all named properties | No |
| `secondary/loader.py:66-74` | `.atom_index`, `.ring_type`, `.n_pairs` — named properties | No |
| `dataset.py:39` | Checks `ring_contributions.npy` exists (file presence only) | No |

**Tests (must re-extract fixtures):**

| File | What | Action |
|------|------|--------|
| `python/tests/test_load.py` | 79 tests against real extraction data | Re-extract test fixtures with 59-column format |

**UI (check if applicable):**

| File | What | Action |
|------|------|--------|
| Viewer ring_contributions reader (if any) | May expect 57 columns | Update to 59 |

**Summary:** The learning and analysis code never touches raw column
indices on ring_contributions — it goes through the SDK's named
accessors.  The SDK is the only code that knows the column layout.
Zero silent breakage.  The named accessor pattern did its job.
