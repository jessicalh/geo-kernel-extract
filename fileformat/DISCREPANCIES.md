# H5 vs NPY Discrepancies — What Stage 2 Code Must Know

**Date:** 2026-04-15
**Context:** 10 reference proteins with analysis H5 files (626 frames
each) and NPY snapshot directories at ~1ns intervals (~26 per protein).
260 ORCA DFTs running on PDB snapshots from those intervals.

Stage 2 analysis will compare kernel time series (from H5) against
DFT shielding tensors (from ORCA, matched to NPY snapshots).  The H5
and NPY were written by different code paths in the same C++ run.
They are *almost* byte-identical — but not quite.  This document
exists so that code working across both paths does not mistake a
known discrepancy for a signal or a bug.

---

## The two data paths

```
                          ┌─ WriteAllFeatures ─→ NPY (at ~1ns snapshots)
Calculator::Compute() ──┤
                          └─ AnalysisWriter ───→ H5  (every sampled frame)
```

Both paths read the same ConformationAtom fields from the same
computed frame.  The calculator output is identical.  The
discrepancies below are all in the *harvest* step — how
AnalysisWriter serialises the data differs from how
WriteAllFeatures serialises the same data.

The kernel T2 tensors, shielding contributions, E-fields, EFGs,
charges, SASA, water features, and bonded energies are
**byte-identical** between the two paths at every snapshot frame.
The cross-validation (29,250 checks, 0 failures) proves this.

The following 7 items are the *complete* list of differences.
Nothing else differs.

---

## 1. Ring geometry normals: different algorithm, same ring

| | NPY | H5 |
|---|---|---|
| **Writer** | `Ring::ComputeGeometry` | `AnalysisWriter` inline |
| **Algorithm** | Cross-product of two ring diagonals | Cross-product of edges [0]→[1] and [0]→[2] |
| **Center** | byte-identical | byte-identical |
| **Radius** | byte-identical | byte-identical |
| **Normal** | up to ~18° difference on thermally distorted rings | |

**Why it doesn't confuse you:**  Neither of these normals enters the
kernel computation.  The calculators (Biot-Savart, Haigh-Mallion,
etc.) use `GeometryResult`'s SVD-fitted normal, which is a *third*
algorithm — and that normal is baked into the kernel T2 tensors that
*are* byte-identical between paths.  The `ring_geometry/` datasets
are diagnostic metadata: "where was the ring, roughly which way did
it point."  They do not participate in any regression or comparison
against DFT.

**Rule:**  Do not use `ring_geometry/normal` from the H5 as a
reference frame for physics comparisons.  If you need the ring
normal that the calculators actually used, it is implicit in the
kernel output.  If you need a consistent normal for visualisation,
pick one algorithm and stick with it — do not mix NPY and H5
normals in the same plot.

---

## 2. Chi dihedral cos/sin: 1 ULP difference (~1e-16)

| | NPY | H5 |
|---|---|---|
| **Writer** | `DsspResult` (libdssp chi computation) | `AnalysisWriter::Dihedral()` (raw atom positions) |
| **Max difference** | 1e-16 (1 ULP of float64) | |

**Why it doesn't confuse you:**  This is two independent `cos()`
calls on the same angle.  1 ULP is the smallest possible difference
for float64.  No analysis code should be sensitive to this.  Any
comparison that uses `==` on a chi cos/sin value has a different
problem.

**Rule:**  Use toleranced comparison (1e-14 or wider) when matching
chi values between paths.  Or just use one path consistently.

---

## 3. Four GROMACS energy terms missing from H5

Missing in H5 `energy/` group (present in NPY `gromacs_energy.npy`):

| NPY column | Name | Why it's absent from H5 |
|---|---|---|
| 2 | `coulomb_14` | 1-4 intramolecular Coulomb — small, nearly constant for a given topology |
| 10 | `lj_14` | 1-4 LJ — same reason |
| 11 | `disper_corr` | Long-range dispersion correction — single scalar, nearly constant |
| 14 | `total_energy` | `potential + kinetic` — exactly reconstructable |

**Why it doesn't confuse you:**  These terms are system-level
thermodynamic bookkeeping from the GROMACS EDR.  They are not
per-atom, not per-residue, and not features in any kernel model.
If you need them for a frame, they are in the NPY at the snapshot
intervals, or reconstructable from the 39 terms that *are* in the
H5.

**Rule:**  Do not assume the H5 `energy/` group is a complete copy
of the EDR.  It has the 39 terms listed in the fileformat object
model.  If you need `total_energy`, compute `potential + kinetic`.
If you need 1-4 terms or dispersion correction, read the NPY.

---

## 4. `aromatic_E_magnitude` not in H5

The scalar magnitude of the aromatic-only Coulomb E-field
(`coulomb_scalars.npy` col 3) is not harvested into the H5.

**Why it doesn't confuse you:**  The aromatic E-field *vector*
(3 components) IS in the H5 as `efg/E_aromatic (T, N, 3)`.  The
magnitude is `np.linalg.norm(h5["efg/E_aromatic"][t], axis=-1)`.
It is a derived scalar from a stored vector.

**Rule:**  Compute it from `efg/E_aromatic` if you need it.  Do not
create a separate dataset for it — that would be a shadow copy of
information already present.

---

## 5. Ion distance sentinel: +inf (NPY) vs -1.0 (H5)

When no ion is within the cutoff radius:

| | NPY (`hydration_shell.npy` col 2) | H5 (`water/nearest_ion_dist`) |
|---|---|---|
| **Sentinel** | `+inf` | `-1.0` |

The AnalysisWriter transforms at harvest time
(`AnalysisWriter.cpp:217-218`).

**Why it doesn't confuse you — if you know about it:**  Both
sentinels are unambiguous (distances are non-negative, so -1.0 is
impossible; +inf is impossible physically).  But code that reads
both paths and does not account for the sentinel difference will
compute nonsensical deltas.

**Rule:**  When comparing `nearest_ion_dist` across paths,
normalise sentinels first.  The cross-validation script does this
correctly (lines 298-312 of `npy_cross_validate.py`).  Pick a
convention and convert on load.  The H5 convention (-1.0) is
recommended for new code because it avoids inf arithmetic.

---

## 6. DSSP data: per-residue (H5) vs per-atom broadcast (NPY)

| | NPY | H5 |
|---|---|---|
| **Shape** | `(N,)` — residue values broadcast to every atom in the residue | `(T, R)` — stored at natural per-residue granularity |
| **Affected datasets** | `dssp_backbone.npy` (phi, psi, helix, sheet), `dssp_chi.npy`, `dssp_ss8.npy` | `dihedrals/`, `dssp/` |

**Why it doesn't confuse you:**  The values are identical —
`h5_value[residue_index[atom]]` == `npy_value[atom]` at every
position.  The H5 is the cleaner representation (no redundancy).
Broadcast via `atoms/residue_index` to get per-atom if needed.

**Rule:**  When joining H5 DSSP/dihedral data against per-atom
kernel data, broadcast explicitly using `atoms/residue_index`.
Do not assume H5 arrays are always `(T, N)` — check the shape.
The fileformat object model documents the shape of every dataset.

---

## 7. Coulomb NPY files absent (--no-coulomb fleet runs)

The fleet runs use `--no-coulomb` because APBS replaces vacuum
Coulomb.  No `coulomb_shielding.npy`, `coulomb_E.npy`, etc. in the
snapshot directories.

The H5 `efg/coulomb_*` datasets ARE populated (CoulombResult runs
regardless; the `--no-coulomb` flag only suppresses NPY output).

**Why it doesn't confuse you:**  The Coulomb data exists in the H5.
It does not exist in the NPY.  There is nothing to compare.  The
calibration uses APBS EFG (`efg/apbs_efg`), not vacuum Coulomb.

**Rule:**  Do not expect `coulomb_*.npy` files in snapshot
directories from `--no-coulomb` runs.  The H5 is the sole source
for frame-level Coulomb data in these runs.  If you need Coulomb
at a snapshot point, read it from the H5 at the matching frame
index.

---

## Summary table for quick reference

| # | Discrepancy | Affects kernel T2? | Affects DFT comparison? | Action needed |
|---|---|---|---|---|
| 1 | Ring normal algorithm | No | No | Don't mix normals across paths |
| 2 | Chi cos/sin 1 ULP | No | No | Use toleranced comparison |
| 3 | 4 missing energy terms | No | No | Reconstruct if needed |
| 4 | aromatic_E_magnitude absent | No | No | Compute from E_aromatic vector |
| 5 | Ion sentinel +inf vs -1.0 | No | No | Normalise on load |
| 6 | DSSP per-residue vs per-atom | No | No | Broadcast via residue_index |
| 7 | Coulomb NPY absent | No | No | Read Coulomb from H5 only |

**None of these discrepancies affect the kernel T2 tensors or the
DFT shielding comparison.**  The physics data that enters any
regression — kernel SphericalTensors, E-fields, EFGs, charges,
SASA, water features, bonded energies — is byte-identical between
NPY and H5 at every validated snapshot frame.

Code that reads kernel time series from H5 and compares against
ORCA DFT at snapshot frames is comparing the same underlying
calculator output through two serialisation paths that agree on
every load-bearing dataset.
