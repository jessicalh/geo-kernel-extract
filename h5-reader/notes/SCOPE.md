# h5-reader — scope

## What this is

A standalone Qt6/VTK reference implementation that opens an analysis H5 file
produced by `nmr_extract --trajectory --analysis` and presents the protein
and its full per-frame physics — positions, kernel contributions, vector
fields, solvent environment, DSSP, dihedrals, bonded energies — as an
animated, inspectable 3D scene.

This program is distributed to Jessica's professors and advisers to run
locally on Windows, macOS, or Linux. Build targets and code quality hold
to the same bar as the main nmr-shielding library: defensible to a
reviewer reading it cold, operable through the UDP log when something
misbehaves, no platform shortcuts.

## What this is NOT

- A writer. It never emits H5 files. It never triggers re-extraction.
- A duplicate protein model. It does not link the `nmr_shielding` library.
  The typed classes inside (`QtProtein`, `QtConformation`, `QtFrame`,
  `QtRing` hierarchy, `QtNamingRegistry`, `QtSphericalTensor`) are its own,
  matched to the H5's serialised shape, not a mirror of the library's types.
- A calculator. The frozen analysis H5 already contains every per-atom
  kernel value at every sampled frame. The reader reads; it does not
  re-derive per-atom results.

The one place the reader evaluates physics is the volumetric Biot-Savart
and Haigh-Mallion field grids around each aromatic ring — the iconic
"butterfly" picture chemists recognise. These grids are not in the H5;
the reader evaluates the closed-form kernels at open-space grid points
using the H5's per-frame ring geometry and the same TOML-calibrated
intensities the extractor used. See `VTK_ANIMATION.md` (future) for the
per-frame update strategy.

## Relationship to the rest of the tree

```
nmr-shielding/
├── src/                 — frozen library (C++). h5-reader does NOT link.
├── ui/                  — the original single-conformation viewer. Frozen
│                          as the library-object-model witness. h5-reader
│                          copies rendering patterns; does not share code.
├── fileformat/          — frozen H5 serialiser/deserialiser (header +
│                          analysis_file.cpp). h5-reader compiles this
│                          directly into its own target, same pattern as
│                          the existing viewer.
├── extern/HighFive/     — vendored header-only HDF5 C++ wrapper.
│                          h5-reader points at this via
│                          HIGHFIVE_INCLUDE_DIR.
└── h5-reader/           — this directory. Standalone Qt/VTK app.
```

The reader is independent of everything except `fileformat/analysis_file.cpp`
and the HighFive headers. When Jessica hands a zip of `h5-reader/` plus
those two dependencies to an adviser, the adviser can build and run.

## The object-model discipline this carries over

Inherited from `spec/CONSTITUTION.md` and lived through in the library:

- **Protein / Conformation split.** `QtProtein` holds identity, topology,
  bonds, ring-class hierarchy. `QtConformation` is the trajectory (one H5
  file = one conformation). `QtFrame` is one sampled XTC frame. Atoms
  query identity through the protein back-pointer; per-frame data through
  the frame slab.
- **Typed enums, never strings for identity.** Element, AtomRole,
  Hybridisation, BondCategory, BondOrder, RingTypeIndex, DSSP 8-class,
  ProtonationVariant — all enums decoded once at H5 load.
- **Ring class hierarchy reconstructed from `topology/ring_type`.**
  `QtPheBenzeneRing`, `QtHisImidazoleRing`, `QtHidImidazoleRing`,
  `QtHieImidazoleRing`, `QtIndolePerimeterRing` etc. — virtual
  `Intensity()`, `NitrogenCount()`, `JBLobeOffset()`. Same physics the
  library bakes into its ring types.
- **NamingRegistry at the H5 boundary.** CHARMM-flavoured `atom_name`
  entries decode to canonical for display; unknown names rejected, not
  silently passed through.
- **SphericalTensor preserved, never collapsed to T0.** Every per-atom
  shielding dataset is (T, N, 9) decoded to `QtSphericalTensor` with
  T0 (scalar), T1 (3-array), T2 (5-array). Glyphs render T2 as ellipsoids;
  scalars expose T0 where useful.

## Reference implementation bar

- Every diagnostic string states values at the point of failure.
- Every QObject constructor registers with the census (`CENSUS_REGISTER`).
- Every signal/slot connection goes through the auditor (`ACONNECT`).
- Every thread-sensitive method asserts affinity (`ASSERT_THREAD`).
- Every external-library boundary catches and logs, never swallows.
- Every file I/O path goes through `QFile`/`QDir`/`QStandardPaths` —
  no POSIX-only assumptions.
- UDP log on port 9997 is the primary debug channel. `udp_listen.py`
  tails it. When misbehaviour appears, the log is consulted BEFORE
  code changes.
- No `QTimer` on interactive controls. Playback uses a dedicated
  `QtPlaybackController` whose timer is load-bearing; that is the one
  legitimate exception.
