# Rediscover — current state (2026-05-31, EOD)

Freshest current state, for the next session/agent. Read `GUIDANCE.md`
(orientation) + `DESIGN.md` (class model) first; this supersedes the
"Status" section in GUIDANCE.md, which a build agent left mid-build/stale.

## UPDATE 2026-05-31 (late) — multi-scenario surface BUILT, faithful-rebuild gate PASSED

The general surface is implemented and validated (built by **codex** — codex has
full unsandboxed agency via `codex exec --dangerously-bypass-approvals-and-sandbox`
from the lead; Claude Agent-tool subagents are sandbox-blocked from compile, a
wrong-knob misconfig, not a hard limit — see `reference_subagent_build_agency`).
Committed `1bd61aa` on `h5-reader-pysr-spike` (NEVER MERGED).

- **Spine built** (`src/rediscover/`): `Catalog`, `TemporalIndex`, `TypedAtomIndex`
  (scoped `select`/`selectUnique`, IUPAC-trap-safe — the positional `front()` anchor
  fallback is REMOVED), `SpatialIndexSet` (per-cloud KD trees, near + range/annulus),
  `RingGeometryCache`, `ChargeStore` + FF14SB read from `topol.top` inline `[atoms]`
  (typed resnr/order cross-check, no glob), `ResidentIndexes`, `OutputManifest`,
  `AnalysisBody`. Per-relationship schema + `relationship_kind` + T2 sidecar NPYs
  documented in `python/nmr_extract/_catalog.py`. No `PbcCellSeries` (PBC=None).
- **ring_current + mcconnell ported** to the Body/catalog/index surface and
  **reproduce the ORACLE from a fresh rebuild** (not the one-off output):
  ring k=21.1, coupled within-atom R²=0.616, equivariant T2 R²=0.468, |T2| r=0.757,
  basis 4.88e-8, frame rot ~1e-4°; mcconnell scalar R²=0.55, kernel readout
  r=0.918/R²=0.843. `ctest h5reader_rediscover_tests` passes; GUI untouched.
- **The 7 others fail loud** (exit 2, ValidateScenario): buckingham_efield, efg,
  charge_dipole, charge_quadrupole, larsen_hbond, charge_response_gradient,
  aimnet2_embedding. MOPAC charge source = absent→loud (per-frame data lands AM);
  AIMNet2 recognized but multipole reducers intentionally unrunnable.

**Next (build, when data/decisions land):** wire the 7 fail-loud stubs — charge
multipoles once charges flow (FF14SB done, MOPAC AM, AIMNet2 charge), the
per-atom-feature items (efg/efield/CRG/embedding) against the carrier, larsen
once its detection/classifier decisions are made. Equivariant-T2 path is proven on
ring; extend to the new items per the frame resolution. The one emergent issue
codex hit + handled: topol↔model atom-name aliasing during charge load (typed
residue/order cross-checks, not positional).

## UPDATE 2026-05-31 (evening) — canonical re-run done + scalar fit landed

The value-affecting Codex fixes are IN, built clean (lead session), re-run on
1P9J, and verified. Canonical output: `/tmp/rediscover-out-v2/` (same shape:
141000/20500 ring, 812205/26000 mc). Fixes that landed:
- **Ring identity + self/bonded exclusion** — `ring_index` + `is_self_or_bonded`
  per source; aggregated split `sum_dipolar_all` vs `sum_dipolar_producer_valid`
  + `n_sources_valid`. Self/fused detected by identity (own-ring set + shared-atom
  overlap), NOT a distance proxy; self-rings verified at r≈2.49 Å, 21.6% of rows.
- **Aromatic-H frame anchor → typed CG/CD2** (the topology-convention time bomb):
  `typedRingAnchor` picks the unique γ-carbon (or δ-carbon CD2 for TRP-benzene) by
  typed `Locant`, no name strings, no positional first-atom. `frame_anchor_atom_index`
  emitted, 100% populated (13 distinct anchors).
- **Ring normal flipped to canonical traversal** `(v1−v0)×(v2−v0)` (local to the
  rediscover code; `FitRingGeometry`/GUI untouched).
- **McConnell**: `bond_axis_local` (unit axis in local frame, verified |·|=1) +
  `bond_atom_a/b` endpoints; `cutoff_A` recorded (CLI `--mc-cutoff`, default 8.0,
  the conventions' aromatic value — producer's exact MC cutoff still unverified).
- nanoflann `[[nodiscard]]` captured.

**Scalar fit ran and is credible (see `analysis/FINDINGS.md`).** Pipeline:
sum-pooling NN learns the per-source function, PySR distils it. Recovering the
Pople form from the producer kernel `bare_T0` is partly CIRCULAR (reverse-engineers
the producer's own Giessner-Prettre formula) — demoted from the headline. The
NON-circular, instantaneous result (ring current is a static map; the trajectory
is a geometry sampler, not a process): a universal coefficient k≈21 ppm·Å³ predicts
the within-atom shielding of HELD-OUT ATOMS (leave-atoms-out, autocorrelation-free)
at R²=0.62 on the ~7 coupled atoms, against independent DFT, on the identity-clean
`sum_dipolar_producer_valid`. Thin (one protein, ~7 coupled aromatic H) but real.

**T2 Cartesian-frame caveat — RESOLVED (2026-05-31 eve), the disciplined way.**
The check lives in the READER, not a Python h5py hack (the reader owns H5): the
ORCA parser now additively retains the `CARTESIAN COORDINATES (ANGSTROEM)` block
(`DftAtomShielding::orca_coord`), and `ExtractionSupport::CheckDftFrameAlignment`
Kabsch-compares the ORCA-tensor frame to the H5 positions the extractor already
holds. Verdict on 1P9J (500 frames × 846 atoms): rotation mean 8.9e-5°, max
2.4e-4°; RMSD 0.0005 Å. The ORCA tensors ARE in the H5 frame — **T2 components
are valid as emitted, no rotation needed** (the sub-millidegree residual is
float-print precision). The equivariant T2 model (task #23) is UNBLOCKED.

**Still open:** (b) DFT
validation (element-equality, raw total≈dia+para) — DEFERRED (defensive, not
fit-affecting). (c) the literature-coefficient-FIXED test (un-fitted constant →
DFT) — the final de-circularising check, not yet run. (d) McConnell producer-kernel
reconstruction gap (R²≈0.55 — likely fuller anisotropy than one bond-axis angle).

## The extractor works end-to-end (verified in the lead session)

`h5reader-rediscover` (branch `h5-reader-pysr-spike`) **compiles, links,
and runs on 1P9J**, and the basis fixture passes (`√6`).

- Build: `cmake --build build/linux-gcc --target h5reader_extract h5reader_rediscover_tests`
- Test: `ctest --test-dir build/linux-gcc -R h5reader_rediscover_tests`
- Run: `build/linux-gcc/h5reader_extract --run /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft --out <outdir> --case all`
- Verified output: 4 CSVs — `ring_current_{sources,aggregated}.csv`
  (141000 / 20500 rows; 41 aromatic ring-facing H × 500 DFT frames),
  `mcconnell_{sources,aggregated}.csv` (812205 / 26000; 52 backbone HN ×
  500). 846 atoms, 751 frames, 500 DFT rows, 0 gaps. Two row kinds per
  case (un-summed per-source + aggregated summed-feature) sharing the
  identity + DFT-target columns.

**Verified:** builds, runs, well-formed tables, basis `√6`.
**NOT yet verified:** the physics correctness of the *values* (Codex is
checking — see below).

## Fix already applied

`RingCurrentNeighborhood.cpp` + `McConnellNeighborhood.cpp` had the wrong
include `"../model/QtTrajectoryH5.h"` → corrected to `"../io/QtTrajectoryH5.h"`
(that header lives in `src/io/`). This was the only compile error.

## Two confirmed must-fix caveats (Jessica confirmed both real)

1. **Aromatic-H local-frame anchor — topology-convention violation, the
   project's classic time bomb.** The frame x-axis is anchored on the ring's
   FIRST canonical-walk atom (positional) instead of the chemistry-typed
   anchor the conventions require: CG for PHE/TYR/TRP-pyrrole/HID/HIE/HIP,
   CD2 for TRP-benzene (per `spec/substrate_conventions_2026-05-30.md`). This
   is identity-from-position — the banned pattern (`feedback_identity_from_chemistry_not_position`,
   the IUPAC-revert episode). FIX: anchor x on centroid→(typed atom);
   **re-read the exact per-residue anchors, do not guess.** Invariants
   (z = ring normal, r, cosθ = z/r) are unaffected; only the azimuthal x/y
   frame is wrong. File: `src/rediscover/LocalFrameBasis.{h,cpp}`.
2. **T2 Cartesian frame.** ORCA DFT tensors are in the ORCA-input frame; H5
   kernels in the MD frame. If those orientations differ, T2 *components*
   don't correspond — comparison is meaningless. T0 (iso) and |T2|
   (magnitude) are rotation-invariant and safe. RESOLVE empirically:
   compare ORCA-input Cartesian coords (from a `.out`) to H5 positions at
   that frame. Same → T2 stands; rotated → rotate the DFT tensor into the
   H5 frame via the shared atoms, or restrict T2 to invariants + flag
   components unverified. **Do not assume.**

## In flight: Codex critique

A Codex critique of the WORKING code + output is running (background,
started ~16:30). Output log: `/tmp/claude-1000/-shared-2026Thesis-nmr-shielding/37ab6452-233d-4c70-9077-8b3440b92984/tasks/blu78kz6j.output`
(huge log — grep near the end for the verdict). It assesses the physics of
the values, both caveats above (incl. running the frame-match check),
match-to-design, the additive edits, and the nanoflann `radiusSearch`
`[[nodiscard]]` warning. Read it, fold it in, then apply the two fixes +
whatever else it surfaces.

## Agent build-agency (flagged twice)

Subagents' Bash is sandboxed in this environment and denied the compiler/
CMake even with their override; only the lead session's
`dangerouslyDisableSandbox` works. So agents could WRITE but not BUILD —
which is why the include bug shipped blind. Agents have had build agency
before; this is a sandbox/settings restriction to RESTORE so the dozen
incoming agents can compile what they write. NOT yet fixed. Interim:
lead/Codex builds; agents write + review.

## Next actions (rough order)

1. Read Codex's report (`blu78kz6j`); apply its fixes.
2. Fix caveat 1 (anchor → typed CG/CD2, re-reading the convention) — committed.
3. Resolve caveat 2 (frame check → rotate or flag).
4. Restore agent build-agency (settings) for the incoming agents.
5. Run a fitter on the substrate (task #22) ONLY after correctness is
   settled. The fitter is OPEN — ridge / scalar SR / equivariant SR /
   equivariant sum-pooling; NOT decided, NOT PySR-only, GNN/equivariant
   not foreclosed.

## Durable discipline (unchanged)

Substrate is method-agnostic; fitter open. Truthful docs, no hyperbole.
Reuse the reader; additive edits only; GUI untouched. Experimental
one-shot branch (no integration). Guardrail memory:
`feedback_build_inmemory_export_dont_relitigate`.

## Codex critique findings (landed, blu78kz6j)

**Confirmed correct:** builds + test pass; scalar physics right — ring
`cosθ=z/r` + `(3cos²θ−1)/r³`; McConnell bond-axis `(3cos²θ−1)/r³` (not /3r³);
`σ_iso=trace/3`; library basis ordering + `√6` fixture; intensity from
`LiteratureIntensity`; per-source rows sum exactly to the aggregate. No
`inf`; `nan` only in 30 ring rows (`ring_in_plane_angle`, near-axis, expected).

**Real problems to fix (source semantics — the scalar math itself is sound):**
1. **Self/bonded ring included** (biggest new physics finding): ring-current
   uses every H5 ring slot incl. the H's OWN parent ring; the producer
   excludes self/bonded. Fix: emit `ring_index` + an `is_self_or_bonded`
   flag; split `sum_dipolar_all` vs `sum_dipolar_producer_valid`.
2. **Ring source identity missing** — add `ring_index` to `SourceSlot` + CSV
   (McConnell emits `bond_index`; ring doesn't).
3. **Aromatic-frame normal sign not canonicalized** — SVD normal not flipped
   to the ring-traversal convention (library `Ring.cpp:32` does). Scalar OK;
   `disp_local` / local-tensor / azimuthal components can randomly flip.
   Fix: flip normal per `(v1−v0)×(v2−v0)`.
4. **Aromatic-H anchor (caveat 1, confirmed)** — `ring.atomIndices.front()`
   not the typed CG/CD2; works today only by coincidence. Fix: typed anchor
   + emit `frame_anchor_atom_index`.
5. **McConnell source rows lack the bond-axis vector in the local frame** —
   an equivariant/tensor fit can't reconstruct the McConnell tensor from
   scalar+midpoint. Fix: add `bond_axis_local_{x,y,z}` + endpoint indices.
6. **McConnell cutoff hard-coded 8.0 Å** — producer uses a configurable
   cutoff (commonly 10 Å); design says don't hide cutoffs. Fix: record/require
   it; use 10 Å when comparing to the producer.
7. **T2 frame unverified (caveat 2, confirmed)** — fix: mark T2 columns
   frame-unverified, OR Kabsch-verify ORCA-input vs H5 frame and rotate the
   tensors back if ORCA reoriented.
8. **DFT validation weaker than spec** — `DftShieldingLoader` doesn't check
   parsed-element == protein-atom element, nor raw 3×3 `total≈dia+para` (only
   iso). Fix: add both.
9. **nanoflann nodiscard** (`FrameSpatialIndex.cpp:59`) — cosmetic (`hits`
   filled by reference); tidy anyway.

**DECIDED — C–H only for now** (Jessica, 2026-05-31). Honest caveat she
raised: the accumulated narrowing choices — single protein, DFT subset,
C–H only, cutoffs — may starve the fit of statistical depth; flagged, and
we follow it to the end regardless (report effective N alongside any fit).
`IsAromaticRingHydrogen` (HA/H4/H5) = aromatic
**C–H only**; it excludes N-bound aromatic/exchangeable H (TRP `HE1`,
protonated HIS N–H). Right if the stratum is aromatic C–H; if it should be
"all aromatic-ring-attached H," change to parent-heavy-atom-in-aromatic-ring.

Codex bottom line: core formulas, basis, tensor retention, row shape, and
sampled magnitudes are sound; the real risks are source semantics
(self/bonded rings, unverified tensor frame, missing orientation/provenance
columns, normal/anchor conventions).
