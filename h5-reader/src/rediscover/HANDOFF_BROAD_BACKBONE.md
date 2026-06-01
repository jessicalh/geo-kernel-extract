# HANDOFF — broad_backbone (composed) + backbone frames + frames test

Branch `h5-reader-pysr-spike`. NEVER MERGE. Authored 2026-06-01 by an Agent-tool
session whose sandbox **denied the compiler** (cmake/gcc present but `cmake
--build` blocked — the known Claude-subagent compile block, per
`reference_subagent_build_agency`). So this is written, self-reviewed, and
UNTESTED-AT-COMPILE. Codex takes the build-fix-run loop.

## Build command (the gate's entry point)

```
cd /shared/2026Thesis/nmr-shielding/h5-reader && cmake --build build/linux-gcc \
  --target h5reader_extract h5reader_rediscover_tests \
           h5reader_rediscover_backbone_frame_tests
```

Then ctest the frames fixture:
```
ctest --test-dir build/linux-gcc -R h5reader_rediscover_backbone_frame_tests -V
```

Run (cutoff sweep — record each):
```
build/linux-gcc/h5reader_extract --run <1p9j-calcset> --out <dir> \
  --case broad_backbone --charge-source ff14sb \
  --ring-cutoff 8 --bond-cutoff 8 --charge-cutoff 6     # then 10, then 12
```

## What was written (all additive)

1. **Backbone local frames** — `LocalFrameBasis.{h,cpp}` (NEW functions; the
   existing `BuildHNFrame` / `BuildAromaticHFrame` are UNTOUCHED so the ring/mc
   oracle byte-parity holds). Added `FrameVariant` enumerators BackboneN /
   BackboneN_NTerminus / BackboneCA / BackboneCarbonylC / BackboneCarbonylO /
   BackboneHA, and builders:
   - `BuildBackboneNFrame(n, ca, cRef, c_prev_valid)` — **the DEFINED N frame**:
     z = unit(CA−N) (the N→Cα bond, the dominant covalent axis at N — N has no H
     in PRO so the HN z=N→H choice can't generalise); x = in-plane component of
     the peptide reference ⟂ z, where ref = (C_prev − N) interior, falling back
     to (C_own − N) at the N-terminus. {N, CA, ref} spans the peptide plane N
     sits in. Directly analogous to the HN frame but with N→Cα as z.
   - `BuildBackboneCaFrame(ca, n, c)` — z = unit bisector(Cα→N, Cα→C); x = Cα→N.
     (conventions doc "Cα frame".)
   - `BuildBackboneCarbonylCFrame(c, o, ca)` — z = unit(O−C) (the carbonyl bond,
     == the McConnell reference); x in-plane ref C→CA. (conventions "C=O frame".)
   - `BuildBackboneCarbonylOFrame(o, c, ca)` — z = unit(C−O) (O→C); x in-plane
     ref C→CA. Same peptide plane, referenced from O.
   - `BuildBackboneHaFrame(ha, ca, n)` — z = unit(HA−Cα); x = Cα→N. (conventions
     "HA / Cα chirality frame".)
   All right-handed (y = z×x), all return `is_valid=false` (no NaN) on
   coincident/collinear anchors via the shared `finishFrame` Gram-Schmidt.

2. **Frames unit test** — `tests/rediscover_backbone_frame_tests.cpp`, NEW
   executable `h5reader_rediscover_backbone_frame_tests` (CMakeLists, links only
   LocalFrameBasis + Eigen + Qt::Test — pure synthetic Vec3, no calcset/H5).
   Asserts, per class, on KNOWN coordinates: axes unit-length + mutually
   orthogonal; right-handed (x×y≈z); z along the defined convention; x-azimuth
   toward the typed anchor (in-plane projection has +x, ~0 y); valid==true; plus
   degenerate (coincident & collinear-bisector & in-plane-ref∥z) → valid==false,
   no NaN. (I hand-verified the Cα/CO numeric expectations; codex confirms under
   ctest.) Sibling target, NOT folded into `h5reader_rediscover_tests`, because
   that exe owns a `QTEST_GUILESS_MAIN` already — one main per Qt::Test exe.

3. **broad_backbone composed relationship** — `BroadBackbone.{h,cpp}`:
   - `MakeBroadBackboneRelationship(ring_cut, bond_cut, charge_cut, charge_src,
     exclude_residue)` returns a `BroadRelationship` = a shared `Relationship`
     bundle (the SAME curried closures the engine runs) + a `BroadReducer` +
     recorded cutoffs.
   - **Stratum** = `atomsWhere([](a){ return a.IsBackbone() || a.IsAnyAlphaHydrogen(); })`.
     `IsBackbone()` (typed BackboneRole != None) covers N/CA/C/O/HN + non-GLY HA;
     GLY HA2/HA3 carry Locant::Alpha + BackboneRole::None, so `IsAnyAlphaHydrogen()`
     is OR'd in. Chemistry, not a name scan. (No `IsBackbone` enum was added —
     `QtAtom::IsBackbone()` already exists, QtAtom.h:79.)
   - **frame_fn** = `backboneFrameFn` — dispatches on the typed BackboneRole (+
     the GLY-HA case), looks up anchors via the **resident** `QtResidue`
     N/CA/C/O/HA cache (built at load from BackboneRole+Locant — collision-safe,
     never positional) and the typed `prevResidueIndex` link. No re-parse, no
     name scan, no positional anchor.
   - **selectors** (a LIST) = `nearBackend(RingCenters, ring_cut)`,
     `nearBackend(BondMidpoints, bond_cut)`, `nearBackend(ChargeSites, charge_cut)`
     — the GENERAL KD backends, NOT the aromatic-H `slotsBackend()`.
   - **attachers** (a LIST, each branches on the typed source kind):
     `ringAttacher` (ring via RingGeometryCache + QtRing virtuals — the general
     KD ring backend reading `raw.ref.entity_index`, NOT the slots `raw.ring`),
     `bondAttacher` (mirrors mcAttacher), and `makeChargeAttacher(src, excl)`
     (FF14SB charge via Catalog `value(Ff14sbCharge,…)`; fills q + disp_local;
     rejects with `source_atom_index = -1` sentinel).
   - **source_filter** drops NaN-dipolar ring/bond rejects and charge rejects
     (`source_atom_index < 0`).
   - **broad_reducer** sums per-mechanism dipolar (ring/bond) and builds the
     local Coulomb **FIELD** (not μ): `E = Σ qᵢ (r_atom − rᵢ)/r³` in the local
     frame (= −Σ qᵢ·disp_local/r³ since disp_local = rᵢ − r_atom), with
     `field_z` = E·ẑ (the Buckingham r=+0.46 axis-projected signal) and `|E|`.
     **μ is kept** (`charge_mu_local = Σ qᵢ·disp_local`) for the tensor story.
   - **target_fn** = `BuildTarget` (raw 3×3 + library T0/T1/T2 + total_local;
     T0 = σ_iso). No bare-kernel cross-check (no single producer kernel spans
     all backbone atoms × all mechanisms).
   - `RunBroadBackbone(brel, body, sink)` drives the SAME closure protocol as
     `RelationshipEngine::RunRelationship` (stratum → frame_fn → flatten
     selectors → attachers → source_filter), diverging ONLY in the reducer
     output shape (BroadAggregate) and the sink (below).

4. **Carrier with the target-repeat FIX** — `BroadBackboneSink.{h,cpp}` (NEW;
   does NOT touch `RecordSink`, so ring/mc + the procedural cells keep exact
   byte-parity). Two-kind output keyed by a `row_id`:
   - `broad_backbone_aggregated.csv` — ONE row per (atom,frame): identity +
     frame + per-mechanism summed scalars + the field (E, field_z, |E|) + μ + the
     full DFT target (raw 3×3, σ_iso, T1, T2, local_frame_valid). The target
     lives HERE, once.
   - `broad_backbone_sources.csv` — one row per (atom,frame,source): `row_id` +
     mechanism + source geometry + per-mechanism identity. **NO target columns**
     (the 828 MB charge_dipole bloat fix — the ~50-col target no longer repeats
     per source row; join on `row_id`).
   - NPYs (one ArraySpec each in `python/nmr_extract/_catalog.py`, additive):
     `broad_backbone_aggregated_target_T2.npy` (n,5),
     `broad_backbone_aggregated_target_local_T2.npy` (n,5),
     `broad_backbone_aggregated_field_local.npy` (n,3). Keyed by aggregated-row
     order.

5. **CLI** — `main_extract.cpp`: `--case broad_backbone`, new `--ring-cutoff`
   (8.0) / `--bond-cutoff` (8.0); reuses `--charge-cutoff` (6.0, swept) +
   `--charge-source` (ff14sb only — ValidateScenario fails loud otherwise and
   when FF14SB charges are absent). The broad case runs on its own
   BroadBackboneSink + manifest, then returns. `--case all` does NOT include
   broad (deliberate opt-in).

## Self-review (compiler-as-adversary, since I couldn't run it)

- `QtAtom::IsBackbone()` exists (QtAtom.h:79) — no enum addition needed.
- `QtResidue` exposes N/CA/C/O/H/HA (NONE = -1 sentinel) + prevResidueIndex —
  used for every anchor; all checked against NONE before use.
- `QtBond` fields (bondIndex/atomIndexA/B/order/category) confirmed (QtBond.h).
- `QtRing::parentResidueIndex` (QtRing.h:72) used for the self-ring flag.
- `RingGeometry{center,normal,radius}` from `verbs::ringGeom` (QtRing.h:35).
- `Catalog::value/present/has(ArrayId::Ff14sbCharge,…)` is the charge read path
  (same as ChargeDipoleNeighborhood) — no topol.top re-read.
- `body.run.conformation->atomPosition(frame, atom)` for the source charge lab
  pos (Conformation.h, same as ChargeDipoleNeighborhood).
- NPY writer duplicated locally in BroadBackboneSink.cpp (same bytes as
  RecordSink's) — additive, no coupling.
- `_catalog.py` parses + the 3 entries instantiate (verified — Python WAS
  runnable; C++ was not).

### Things for codex to watch at compile / first run
- `makeChargeAttacher` returns an `Attacher` (std::function) captured by value in
  the `rel.attachers` vector — fine, but confirm no dangling capture (it captures
  `charge_source` QString + bool by value).
- `BroadBackboneSink` source-row writer prints `s.cos_theta`/`s.dipolar` for
  charge rows (both 0 — charge attacher doesn't set them). Harmless; documented.
- The `field_local` NPY units are labelled `e/Angstrom^2` in the catalog (units
  cancel in a correlation fit; the manifest/_catalog label is the recorded
  convention — verify it reads sensibly, not load-bearing for the fit).
- I could not confirm the 1P9J stratum actually resolves a frame for EVERY
  backbone class at runtime (PRO has no HN → backboneFrameFn returns invalid for
  a PRO N's amide-H, which simply doesn't exist; PRO N still gets the N frame).
  This is exactly what the gate must check: every frame class resolves, no
  narrowness blocks the run.

## The gate (broad — no single oracle; run after a GREEN build)

1. Runs across ALL backbone atom types — every frame class resolves (N/CA/C/O/
   HN/HA), no narrowness blocks the run; well-formed two-kind output + manifest;
   **confirm the carrier fix** — `broad_backbone_sources.csv` has NO target
   columns and does NOT bloat (the 828 MB failure mode); the target appears once
   per (atom,frame) on the aggregated row + the target NPYs.
2. The combined mechanism features (ring_sum_dipolar, bond_sum_dipolar,
   field_z/|E|) explain backbone σ_iso **per-element / per-atom-type** in the
   Stage-1 ballpark — correlate-not-match; report per-atom-type R² + effective N
   (analysis venv `analysis/venv`). Sweep charge-cutoff 6/10/12 Å; note whether
   field_z / the fit climb (6 Å truncates the long-range 1/r² field).
   If a frame class breaks or the carrier bloats, that is the latent narrowness
   found — report it precisely.

## Discipline confirmed
Composed-not-procedural (no `BroadBackboneNeighborhood::extract` walk — the same
`Relationship` closures + verbs, run through a sibling runner that shares the
engine's loop shape); topology reused (frames/anchors/selectors/charges off the
resident QtResidue cache + KD clouds + ChargeStore — no re-parse / re-derived
connectivity / name scan / positional anchor); ring/mc oracle parity intact
(RecordSink + RelationshipEngine + ComposedRelationships + the procedural cells
UNTOUCHED); reader owns H5 (charges via Catalog, no Python H5 read); GUI
untouched; library not linked; NEVER MERGE.
