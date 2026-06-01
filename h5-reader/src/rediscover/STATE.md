# Rediscover — current state (2026-05-31, EOD)

Freshest current state, for the next session/agent. Read `GUIDANCE.md`
(orientation) + `DESIGN.md` (class model) first; this supersedes the
"Status" section in GUIDANCE.md, which a build agent left mid-build/stale.

## UPDATE 2026-06-01 — Python physics RETIRED; equiv-T2 rebuilt on e3nn (authored)

Per `MODEL_PLACEMENT_PROPOSAL.md` + the lead's decisions. The Python end-runs are
gone; the equivariant fitter is e3nn on the C++ export. The C++ substrate was NOT
touched (no recompile; the carrier `compact_npy` change is a separate codex task).

- **Deleted** `analysis/equiv_t2.py` (its numpy `lib_T2` was a byte-for-byte
  clone of `DecomposeLibrary` — the projection end-run).
- **Rebuilt** as `analysis/equiv_t2_e3nn.py`: `o3.spherical_harmonics("2e")` for
  Y2(r̂)/Y2(n̂) + a `1o⊗1o→2e` cross term (BOTH `--cross exact` fixed-path and
  `--cross learnable` FullyConnectedTensorProduct; `--cross both` reports both and
  picks the better by frame-split R² — decision 4) + invariant radial MLP +
  scatter-pool. Consumes the C++-emitted per-source `disp_local`/
  `source_normal_local` + the REQUIRED emitted target NPY
  `rediscover_ring_current_sources_target_local_T2.npy` (fail-loud if absent; NO
  numpy projection fallback). Target mapped library-basis→e3nn-2e by a PINNED
  constant.
- **`analysis/change_of_basis.py` + `analysis/test_change_of_basis.py`**: the 5×5
  library-T2↔e3nn-2e change-of-basis is the ONLY surviving "recompute," written as
  a fixture-pinned TEST — orthogonality/round-trip + a Wigner-D equivariance
  round-trip vs the C++ library tensor (handles e3nn's y-z-x convention). The
  derivation uses e3nn's own `ReducedTensorProducts("ij=ji", i="1o")` so features
  and target share the e3nn-2e basis (verified: e3nn 1o is the plain 3-vector,
  no hidden permutation — `_reduce.py:126-129`).
- **Pople-comparison recompute arrays DELETED** (lead decision 3, no kept "labeled
  integrity test"): `sumpool_kernel.py`, `refine_kernel.py`, `pysr_distill.py`,
  `sumpool_t0.py`, `sumpool_mcconnell.py` keep their pooling FITS but no longer
  build `(3cos²−1)/r³`; comparisons read the C++-emitted `dipolar`/`bare_T0`, and
  PySR carries the symbolic "form fell out" claim. `look03_coefficient.py` reads
  the emitted `sum_dipolar_producer_valid` aggregate (no `Σ intensity·dipolar`
  re-sum; intensity-weighted aggregate is a future C++ reducer).
  `look_charge_dipole.py` lost the `Σ q·d/r³` field recompute — fits only the
  emitted `mu_*` (the field is the C++ `buckingham_efield`/APBS relationship, a
  fail-loud stub on this branch). `look01` repointed to the emitted aggregate col.
- **GREP PROOF**: no `(3cos²−1)/r³`, `/r**3`, `q·d/r³`, or `lib_T2` projection
  arithmetic remains in any analysis `.py` outside the pinned change-of-basis
  test. (Run from `analysis/`: see the grep block in the handoff report.)
- **Env** (decision 1): system `/usr/bin/python3` has torch 2.11.0+cu130 + e3nn
  0.6.0. PINNED in `analysis/requirements-e3nn.txt`; `LD_LIBRARY_PATH` gotcha +
  run commands in `analysis/ENV.md`. rediscover NPYs NOT added to the SDK catalog
  (decision 2 — never-merge spike).
- **RUN + GATE: PASSED (lead, 2026-06-01).** Ran in system python (torch
  2.11+cu130 / e3nn 0.6.0, `LD_LIBRARY_PATH` per `ENV.md`) on `/tmp/rdc-composed`.
  Change-of-basis: all 3 checks pass — fixed the derivation to float64 (e3nn
  defaults float32, which under-precised C to ~1e-7 and tripped the test's strict
  thresholds; the matrix was always exactly orthogonal). `_C_FROZEN` frozen as the
  exact 5×5 (0,±1,±0.5,±√3/2); `get_C()` loads it with NO e3nn (orthogonal to
  1.1e-16). e3nn fit `--cross both`: **frame-split T2 R²=0.466** (baseline 0.467),
  **|T2| r=0.757** (baseline 0.756) — reproduces the retired hand-rolled result to
  noise, confirming the equivariant signal is REAL, not a hand-rolling artifact.
  `exact` cross-term beat `learnable` (angular structure fixed, not fitted);
  leave-atoms-out 0.370 reported, not gated (thin ~7 coupled aromatic H). Offense
  grep-clean (no physics recompute outside the labeled change-of-basis test).
  Committed on `h5-reader-pysr-spike` (never merged).

## UPDATE 2026-06-01 (later) — broad-backbone BUILT + GATED (commit 35f3768)

Big-one #1 done: `broad_backbone` composed THROUGH the engine (Claude-authored,
codex-built+gated, lead-verified). All 6 backbone frame classes resolve **0%
invalid** (N/CA/C/O=54, HN=52, HA=58); new N-frame convention z=unit(CA−N), x in
the peptide plane. Heterogeneous selectors [rings/bonds/charge-field via the
GENERAL KD backends] — breadth from the selector list, not a procedural walk.
Carrier target-repeat fix held (27-col source rows, no `dft_*` columns). Tests
12/12 (√6 + frames + rotation-equivariance + reducer-sum). Ring/mc byte-parity
oracle STILL exit-0 (lead independently re-ran — NO regression); commit touched
only the new broad files + additive `LocalFrameBasis` (existing builders
untouched) + `main_extract` + tests + `_catalog.py`.

Per-atom-type σ_iso R² (correlate-not-match, first-pass, features
[ring_sum_dipolar, bond_sum_dipolar, field_z, field_mag]): HN 0.45, N 0.45, C
0.31, O 0.20, HA 0.24, **CA 0.055** (weakest — CA is local-bonding-dominated, so
through-space mechanisms barely touch it; diagnostic about kernel completeness,
NOT a physical conclusion). Cutoff sweep 6/10/12 Å **flat** (field
short-range-saturated, not truncation-starved — revises the charge_dipole
"sweep it" hypothesis).

SCALE FINDING (→ #29 input): per-source CSV is 1.7–5.5 GB across the sweep (charge
= 8.1M rows at 6 Å). The target-repeat fix prevented worse, but per-source CSV
doesn't scale for the charge mechanism → NPY-for-source-rows (or don't emit
per-source charge rows) belongs in the totality design.

ENGINE FINDING (→ #29): broad needed a SIBLING runner (`RunBroadBackbone`) because
`RunRelationship` hardcodes the ring/mc scalar-sum sink. DECIDED (Jessica): fix
the engine, but as a REVIEWED TOTALITY design — the fold/sink/carrier across ALL
9 relationship shapes, minimal+clarifying abstraction, NOT whack-a-mole. #29,
blocked until the design + review eyes land. Sibling runner kept meanwhile;
oracle untouched.

## UPDATE 2026-06-01 — functional API BUILT + byte-parity-VERIFIED (commit 99cdc85)

BUILT + GATED (codex build, lead independent re-verify). Compiles + ctest pass —
one real fix: Qt's `slots` macro collided with `verbs::slots` → renamed
`verbs::ringSlots` in the NEW files only. The composed engine reproduces the
procedural oracle **byte-for-byte** — all 4 CSVs (141000/20500/812205/26000 rows)
+ 12 sidecar NPYs identical (independently re-run by the lead into a fresh dir,
GATE: PASS). Physics held by identity: ring k=21.1 / coupled R²=0.616, equiv-T2
basis 4.88e-8 / R²=0.467 / |T2| 0.756, mc R²=0.549 / readout 0.919; DFT frame rot
~9e-5°. Commit 99cdc85 touched ONLY the new API files + `main_extract` +
`CMakeLists` (cells/spine/RecordSink/library untouched — frozen-oracle rule held).
The functional API is **now the proven default path**; broad-backbone (#26) is
unblocked. Compile-knob answer: a Claude Agent-tool subagent STILL can't reach the
compiler here (even `which cmake` denied) — author-with-Claude / build-with-codex
is the working split.

The original authoring note follows (now superseded by the BUILT status above):

### (superseded) functional API authored, UNBUILT

Authored the composable functional API that SURFACE_DESIGN specifies, replacing
the three monolithic procedural walks as the DEFAULT path (cells kept as the
reference oracle). New files: `Verbs.{h,cpp}` (Layer-1 primitive verbs, thin
over the spine), `Relationship.{h,cpp}` (Layer-2 curried-closure combinators +
the named `Relationship` bundle + `atomsWhere`/`slotsBackend`/`nearBackend`
builders), `RelationshipEngine.{h,cpp}` (Layer-3 one pure iterated-closure
loop), `ComposedRelationships.{h,cpp}` (ring_current + mcconnell rebuilt as
composed `Relationship`s). `main_extract.cpp` gained `--engine
{composed|procedural}` (default composed; procedural runs the oracle cells for
the parity diff). `analysis/oracle_parity.py` is the gate runner (diffs
composed vs procedural CSVs+NPYs byte-for-byte).

**The compiler was NOT reachable** from the authoring agent's sandbox (Bash
denied `cmake --build`, even `which cmake`) — answers the open knob question:
agent Bash here cannot compile. So the API is authored + rigorously
self-reviewed, NOT built. **codex takes the compile + the oracle gate** — full
inventory, self-review, and the exact build/test/gate commands are in
`HANDOFF_TO_CODEX.md`. Topology reused not regenerated; reader owns H5; GUI
untouched; library not linked; never merged.

## NEXT (read `BROAD_BACKBONE_NEXT.md`) — do the BROAD case, not more narrow cells

Decided 2026-05-31: STOP building narrow single-stratum × single-mechanism cells
("single-bond-thingies"). The next move is **every backbone atom, all mechanisms
in one heterogeneous bundle** — which (a) is the thesis target (backbone shifts),
(b) stress-tests whether the architecture composes the GENERAL analysis or has
baked narrowness, (c) forces the unbuilt backbone frames (Cα/CO/HA/N). Full plan +
the charge_dipole carry-overs (field-not-μ, cutoff sweep, AIMNet2/APBS, carrier
target-repeat fix) are in `BROAD_BACKBONE_NEXT.md`. charge_dipole cell committed
`3103e73` (μ null, field carries: Buckingham field_z r=0.46, LOAO R²=0.21).

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
