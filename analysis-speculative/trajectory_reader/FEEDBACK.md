# Trajectory Reader SDK: Feedback on the H5 and Spec

**Date:** 2026-04-14
**H5 file:** `/shared/2026Thesis/fleet_calibration-working/1B1V_4292/analysis_output/md_analysis.h5`
**Spec:** `spec/ANALYSIS_TRAJECTORY_2026-04-14.md`

---

## 1. Element encoding: spec says atomic number, H5 stores enum index

The spec says:

    atoms/element  (N,)  int32  # 1=H, 6=C, 7=N, 8=O

The actual H5 stores **C++ enum ordinal values**: H=0, C=1, N=2, O=3, S=4.
These are the values of `enum class Element { H, C, N, O, S, Unknown }`,
not atomic numbers.

**Impact:** Anyone reading the spec and writing an R script would get
wrong elements. The spec needs to say "0=H, 1=C, 2=N, 3=O, 4=S"
or the C++ writer needs to convert to atomic numbers before writing.

**Recommendation:** Store atomic numbers. The enum ordinal is an
implementation detail of the C++ code. Atomic number is the physics
identity. An R user who reads `element == 6` knows it's carbon.
An R user who reads `element == 1` thinks it's hydrogen but it's
actually carbon. This will cause silent bugs.

## 2. protein_id is an attribute, not a dataset

The spec says:

    meta/protein_id  scalar  string

The actual H5 stores `protein_id` as an **attribute** on the `meta/` group,
not as a dataset. The attribute value is `"md"` rather than the expected
protein identifier like `"1B1V"`.

**Impact:** The SDK had to check both locations (attribute and dataset).
The value `"md"` is not useful as a protein identifier.

**Recommendation:** Either write it as a dataset per the spec, or update
the spec to say it's an attribute. Either way, use the actual protein
PDB ID (e.g., "1B1V"), not a generic string. The `protein_id` from the
GromacsProtein or from the directory name (1B1V_4292) would be correct.

## 3. Ring topology uses offset+flat-index encoding

The ring atom indices are stored as:

    topology/ring_offsets       (n_rings+1,)  int32
    topology/ring_atom_indices  (total_ring_atoms,)  int32

This is a standard CSR-style ragged array. It works, but the spec doesn't
document this encoding. It just says the topology group exists.

**One oddity:** `ring_offsets` has shape `(n_rings+1,)` which is `(4,)` for
3 rings. This is correct CSR with a trailing sentinel, but the spec should
document it. The last offset is the total number of ring atom entries
(here: 18 = 3 rings x 6 atoms each).

## 4. n_atoms / n_frames / n_residues are attributes, not datasets

These are stored as attributes on the `meta/` group:

    meta.attrs["n_atoms"] = 335
    meta.attrs["n_frames"] = 626
    meta.attrs["n_residues"] = 23

Not mentioned in the spec at all. Useful for pre-allocation without reading
the full arrays, but the spec should document them.

## 5. Hybridisation is all Unassigned

Every atom in the H5 has `hybridisation == 3` (Unassigned). The C++ enrichment
step that assigns sp/sp2/sp3 either does not run in the trajectory path or
does not write to the H5.

**Impact:** Consumers who need hybridisation for chemical context (e.g.,
distinguishing sp2 aromatic carbons from sp3 CH2 groups) cannot rely on
this field. They would need to infer it from AtomRole instead.

**Recommendation:** Either populate this from the enrichment step or
remove it from the H5 to avoid the false impression it's available.
AtomRole already distinguishes AromaticC from SidechainC, so the
information IS there in a different field.

## 6. No chain_id in residues

The C++ `Residue` class has `chain_id` and `insertion_code`, but neither
appears in the H5. For single-chain proteins like 1B1V this doesn't
matter, but multi-chain proteins would need it.

**Recommendation:** Add `residues/chain_id` for future multi-chain support,
or document that the analysis H5 is single-chain only.

## 7. atom_name is stored as bytes, not fixed-length string

`atoms/atom_name` has `dtype=object` (variable-length h5py strings). The
values are `b'CA'`, `b'HB2'`, etc. This is correct but means consumers
must decode bytes to str.

**Minor:** Consider using h5py's fixed-length string type for consistent
behavior across languages. R's `rhdf5` handles both, but fixed-length
is simpler.

## 8. Bond ordering is arbitrary

The bonds are not sorted in any particular order. The first bond in the
H5 is atom 205-206 (a sidechain bond deep in the protein), not the first
peptide bond. This is fine for graph traversal but makes visual inspection
of the bond list confusing.

**Minor:** Sorting bonds by `(atom_index_a, atom_index_b)` would make
the topology easier to inspect manually.

## 9. SphericalTensor convention: 9-component layout is undocumented

Several datasets are `(T, N, 9)` representing a SphericalTensor:
`bs_shielding`, `hm_shielding`, `rs_shielding`, etc.

The 9 components are presumably `[T0, T1[0], T1[1], T1[2], T2[0], T2[1], T2[2], T2[3], T2[4]]`.
This layout matches the C++ `SphericalTensor` struct (T0 scalar, T1 array[3],
T2 array[5]). But this is not documented anywhere in the H5 or spec.

**Impact:** An R user needs to know that index 0 is the isotropic part,
indices 1-3 are the antisymmetric pseudovector, and indices 4-8 are the
traceless symmetric T2 in isometric normalization. Without this, the
9-component vector is opaque.

**Recommendation:** Add an H5 attribute on each 9-component dataset:
`layout = "T0,T1[3],T2[5]"` or similar. Or add a `conventions/` group
with documentation strings.

## 10. Spec groups missing from the H5

The actual H5 only contains: `meta/`, `atoms/`, `residues/`, `positions/`,
`ring_current/`, `topology/`.

Missing from the spec:
- `efg/` (electric field gradients)
- `bond_aniso/` (McConnell)
- `quadrupole/`
- `dispersion/`
- `hbond/`
- `sasa/`
- `water/`
- `charges/`
- `aimnet2_embedding/`
- `predictions/`
- `projections/`
- `dihedrals/`
- `dssp/`
- `events/`

This is expected (the analysis writer is being built incrementally).
The SDK handles missing groups gracefully (returns None). But it means
we could only test the ring_current path against real data.

## 11. topology/ is an undocumented group

The spec's H5 schema does not mention a `topology/` group. The actual H5
has:

    topology/bond_atoms           (B, 2)   int32
    topology/bond_category        (B,)     int32
    topology/bond_order           (B,)     int32
    topology/parent_atom_index    (N,)     int32
    topology/ring_atom_indices    (K,)     int32
    topology/ring_offsets         (n_rings+1,)  int32
    topology/ring_residue         (n_rings,)    int32
    topology/ring_type            (n_rings,)    int32

**Impact:** The spec says the H5 is self-contained, but topology/ is not
in the schema. Consumers need to know it exists to reconstruct bonds and
rings.

**Recommendation:** Add topology/ to the spec schema.

## 12. parent_atom_index sentinel is -1, not SIZE_MAX

The C++ code uses `SIZE_MAX` (`std::numeric_limits<size_t>::max()`) as
the sentinel for "no parent atom" on heavy atoms. The H5 stores int32,
so the sentinel is `-1`.

This is a reasonable conversion (SIZE_MAX as int32 would be -1 anyway),
but it should be documented. The SDK interprets `parent_atom_index < 0`
as "no parent".

## 13. Ring data has no fused_partner_index

The C++ `Ring` class tracks `fused_partner_index` (which ring it's
fused to, e.g., TRP's 5+6 ring pair). The H5 doesn't store this.

**Impact:** Cannot reconstruct TRP's fused ring topology from the H5 alone.
For the current data (1B1V has no TRP), this doesn't matter. For proteins
with TRP, the fused partner relationship would need to be inferred from
shared atoms.

## 14. No is_rotatable on bonds

The C++ `Bond` struct has `is_rotatable`, which is not stored in the H5.
Chi dihedral analysis needs to know which bonds are rotatable.

**Minor for now:** dihedrals are stored per-residue, so the rotatable bond
information is implicit in the dihedral data.

## 15. Frame stride is implicit

The H5 stores `frame_indices` (the original XTC frame numbers) and
`frame_times` (the simulation timestamps). The stride between samples
is implicit (frame_indices goes 0, 2, 4, 6, 8...). It would be useful
to store the stride as a metadata attribute.

---

## Design decisions made in the SDK

1. **Eager loading:** All datasets are read into numpy arrays at load time.
   For 1B1V (335 atoms, 626 frames, ring_current only), this takes <1s
   and ~200 MB. When all groups are populated, the full protein will be
   ~1.2 GB. This is fine for 64 GB machines. For memory-constrained use,
   the loader could be modified to hold h5py Dataset objects (lazy reads)
   instead of numpy arrays.

2. **AtomView on every property access:** `atom.ring_current.bs_T0_per_type`
   creates a new AtomView each time. This is cheap (no data copy -- it
   just stores an index) but creates garbage. An alternative would be
   to cache AtomViews, but the overhead is negligible for analysis scripts.

3. **FrameView vs atom-centric access:** Both are provided. Time series
   analysis (rolling statistics, correlation) uses atom-centric:
   `atom.ring_current.bs_T0_per_type` returns (T, 8). Per-frame cross-atom
   analysis uses FrameView: `protein.frame(t).ring_current("bs_T0_per_type")`
   returns (N, 8). The FrameView also supports `fv[atom_index]` for
   single-atom single-frame access.

4. **No materialized frame objects:** As specified, there are no 626
   conformation objects. A "frame" is an integer index.

5. **Enum values match C++ ordinals:** Since the H5 stores C++ enum
   ordinals (verified against the data), the Python enums use the same
   integer values. This means `Element.H == 0`, not `Element.H == 1`.
   If the H5 switches to atomic numbers, the Element enum would need
   a mapping layer.
