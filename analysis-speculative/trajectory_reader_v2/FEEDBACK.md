# Feedback: Points of Confusion, Missing Info, Design Suggestions

Written after building the SDK and running it against
`1B1V_4292_analysis.h5` (626 frames, 335 atoms, 23 residues, 3 rings).

---

## 1. What the H5 has vs. what the spec says

### 1.1 Missing groups (spec describes them, H5 does not have them)

The following groups from the spec are absent from the actual H5:

- **`dihedrals/`** -- phi, psi, omega, chi1-4, cos/sin pairs. The spec
  lists this as a `(T, R)` group. Presumably not yet implemented in
  AnalysisWriter.

- **`dssp/`** -- ss8 and hbond_energy `(T, R)`. Same.

- **`predictions/`** -- All ridge/MLP prediction datasets. Depends on
  Workstream A (model export) completing first, which is expected.

- **`projections/`** -- Subspace coords, residual projections,
  Mahalanobis distance. Same dependency.

- **`events/`** -- Intentionally placeholder per spec (Phase 3).

**Impact on SDK:** The loader uses `if "group" in h5` guards, so missing
groups are silently None. Consumers should check before accessing.

### 1.2 Extra datasets in H5 not described in spec

The H5 contains datasets not listed in the spec:

- **`atoms/bfs_decay`** `(N,)` float64
- **`atoms/bfs_to_nearest_ring`** `(N,)` int32
- **`atoms/eneg_sum_1`**, **`atoms/eneg_sum_2`** `(N,)` float64
- **`atoms/graph_dist_N`**, **`atoms/graph_dist_O`** `(N,)` int32
- **`atoms/is_amide_H`**, **`atoms/is_alpha_H`**, **`atoms/is_methyl`**,
  **`atoms/is_aromatic_H`**, **`atoms/is_on_aromatic_residue`** `(N,)` int8
- **`atoms/is_hbond_donor`**, **`atoms/is_hbond_acceptor`** `(N,)` int8
- **`atoms/n_pi_bonds_3`** `(N,)` int32
- **`atoms/parent_is_sp2`** `(N,)` int8
- **`atoms/partial_charge`**, **`atoms/vdw_radius`** `(N,)` float64

These are all from ConformationAtom enrichment fields. The spec's
`atoms/` schema only lists element, residue_index, atom_name, atom_role,
hybridisation, n_bonded, graph_dist_ring, is_backbone, is_conjugated.
The actual writer dumps many more enrichment/topology fields.

**Suggestion:** Either update the spec to document the full set, or
decide which are truly per-protein constants vs. artifacts of writing
frame-0 enrichment data. Currently `partial_charge` and `vdw_radius`
are ff14SB constants (from ChargeAssignmentResult), not per-frame.
The booleans (is_amide_H, is_methyl, etc.) are also frame-independent
enrichment that could be documented as topology.

### 1.3 Extra groups in H5 not in spec

- **`per_ring/`** -- The spec does not describe this group at all. It
  is present in the H5 with `K=6` nearest rings per atom, carrying
  geometry and T2 tensors from 6 calculators. The SDK surfaces it as
  `PerRingGroup`.

- **`ring_geometry/`** -- Not in spec. Contains per-ring, per-frame
  geometry (center, normal, radius) as `(T, n_rings, 7)`. The SDK
  surfaces it as `RingGeometryGroup`.

**Suggestion:** Document per_ring/ and ring_geometry/ in the spec. They
carry significant physics data (per-ring T2 decomposition is exactly
where the angular residual lives).

### 1.4 Naming discrepancy: `efg/coulomb_shielding`

The spec does not list `efg/coulomb_shielding`. The spec lists
`efg/coulomb_total` (which is the raw EFG SphericalTensor). The H5
has BOTH `efg/coulomb_total` and `efg/coulomb_shielding` -- these are
different: `coulomb_total` is the raw Coulomb EFG, `coulomb_shielding`
is the Coulomb shielding contribution (the C++ `coulomb_shielding_contribution`
field). The spec should document both.

---

## 2. Atom element encoding

The spec says `atoms/element` stores "atomic number: 1=H, 6=C, 7=N,
8=O, 16=S". The C++ enum `Element { H, C, N, O, S, Unknown }` uses
ordinals 0-5. The H5 uses atomic numbers, not ordinals.

**Confirmed:** The SDK uses `Element(IntEnum)` with values 1, 6, 7, 8,
16 matching the H5 atomic numbers. This is the correct mapping. The C++
writer calls `AtomicNumberForElement()` before writing.

This is a clear case where the H5 is more self-documenting than the C++
enum ordinals -- good design choice.

---

## 3. Ring atom indices are not sorted

`topology/ring_atom_indices` stores ring atoms in the order they appear
in the ring (bonding order around the ring), not sorted by index. For
PHE ring 0: `[38, 39, 41, 43, 47, 45]`. This means fancy indexing into
h5py datasets with ring atom indices requires sorting first (h5py
requires increasing order for fancy indexing).

**Not a bug** -- the bonding order is semantically meaningful (it
defines the ring polygon). But consumers doing h5py fancy indexing must
sort. The SDK documents this.

---

## 4. Per-ring sentinel value

`per_ring/ring_type` uses -1 (as int8) to indicate "no ring at this
slot" when an atom has fewer than K=6 nearby rings. This is visible in
the data: atom 18 has 2 nearby rings and 4 empty slots.

**Suggestion:** Document the sentinel in the spec or as an H5 attribute.
Currently discoverable only by inspection.

---

## 5. SphericalTensor layout

The 9-element SphericalTensor layout is `[T0, T1_x, T1_y, T1_z,
T2_m-2, T2_m-1, T2_m0, T2_m+1, T2_m+2]`. This is not documented
anywhere in the H5. I inferred it from `SphericalTensor::Decompose()`
in the C++ code and from the struct layout (T0:1, T1:3, T2:5 = 9).

**Suggestion:** Store the layout as an attribute on any dataset that
uses it, or document it once in `meta/`. R users reading the H5 need
to know this. The isometric normalization convention (real spherical
harmonics) should also be documented -- it affects how to interpret T2
magnitudes.

---

## 6. Ring geometry data layout

`ring_geometry/data` has shape `(T, n_rings, 7)` with field names in
`ring_geometry/fields`: `['center_x', 'center_y', 'center_z',
'normal_x', 'normal_y', 'normal_z', 'radius']`.

The fields metadata dataset is a good pattern. The per_ring/ group does
the same with `geometry_fields`.

---

## 7. Absent from H5: per-residue dihedral data

The spec describes an extensive `dihedrals/` group with phi, psi, omega,
chi1-4 and their cos/sin pairs. This is absent. For Phase 1 analysis
(DFT frame selection), dihedrals are important for detecting rotamer
flips. This is presumably a priority for the next C++ implementation
pass.

---

## 8. Design suggestions for the SDK

### 8.1 Lazy vs. eager loading

Currently the loader reads all topology arrays eagerly (into Python
objects) but keeps physics data as h5py dataset references. This is the
right split -- 335 atom objects are cheap, but 626x335x256 float64
arrays are not.

For larger proteins (5000 atoms), the topology construction loop may
become noticeable. If that matters, the enrichment fields could be kept
as raw numpy arrays instead of per-object Python attributes. But 5000
iterations with simple assignments should still be under 50ms.

### 8.2 Frame view object (not implemented)

The task mentions "protein.frame(t) returns a lightweight view" as a
design option. I chose NOT to implement this because:

1. It would encourage materializing per-frame state that the h5py lazy
   loading already provides efficiently.
2. The natural access pattern is `protein.ring_current.bs_T0_per_type[:, atom, :]`
   (time series for one atom) or `protein.ring_current.bs_T0_per_type[t, :, :]`
   (all atoms at one frame). Both are one h5py slice.
3. A frame view object would need to wrap every physics group, creating
   complexity without benefit.

If consumers strongly want `protein.frame(t).atom(i).bs_T0`, it could
be added as a thin wrapper that just passes `t` into the slice calls.

### 8.3 h5py dataset vs numpy array

All physics group attributes hold h5py Dataset objects, not numpy arrays.
This means:
- Slicing triggers disk I/O on each access.
- No in-memory caching.
- Thread-safe for reads (h5py handles the HDF5 lock).

For time series analysis that reads the same dataset repeatedly (e.g.
iterating over atoms), the consumer should do `data = dataset[:]` once
to get a numpy array. The SDK does not do this automatically because
the consumer knows their access pattern and the SDK should not assume
it.

---

## 9. What was not confusing

- The protein/conformation split is clear and the H5 embodies it well.
  Constant groups (atoms/, residues/, topology/) vs. time-indexed groups
  (everything else with T dimension) is exactly right.

- The enum integer values match C++ when checked against the headers.
  AtomRole, BondCategory, BondOrder, RingTypeIndex all use sequential
  ordinals starting from 0.

- The CSR encoding for ring atom indices (ring_offsets/ring_atom_indices)
  is standard and works correctly.

- The per_ring K=6 structure with typed geometry fields and per-calculator
  T2 tensors is rich and well-organized. This is clearly designed for
  the kernel assembly path.

- The ring_geometry separate from per_ring is a good split: ring_geometry
  is per-ring per-frame (3 rings), while per_ring is per-atom per-frame
  (335 atoms x 6 slots). Different access patterns, different shapes.
