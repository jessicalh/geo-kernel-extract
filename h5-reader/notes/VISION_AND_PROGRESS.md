# h5-reader — vision, goals, and progress

**Living tracker.** A status snapshot plus the working vision for the
current build. It *points to* the living sources (design, contract, code)
rather than restating them — when this doc and a living source disagree,
the living source wins. Last updated **2026-05-26**.

## Vision

The reader is **the player, the interface, and the exposed formal model**
for the rich numeric trajectory the `nmr_shielding` library produces — the
read-side *programmable mirror* of the library's typed object model. Not a
viewer bolted onto an H5: the canonical typed embodiment of a trajectory
that humans inspect and other systems build on.

Six capabilities, all **projections of one typed model** (build the model
first; layer these onto its seams):

1. **Player** — interactive playback + movie export + relationship animations.
2. **Multi-scale navigation** — µs run → 1–10 ns window → sampled-frame detail.
3. **Observability** — Python scripts referenced/triggered by npy→npy or Welford changes.
4. **GNN export** — topology = edges, per-atom groups = node features (irreps preserved).
5. **Network streaming** — post per-frame state / graph / events while it plays.
6. **The formal model itself** — typed, no-strings; the join key back to the library.

**Topology is the spine**: it relates the reader's model to the library's
and makes every numeric buffer interpretable in topology terms. **Labels
are projections, never identity — the reader is not label-driven.**

## Locked design decisions

- **Mirror the library's ConformationResult groups** as typed read-side
  groups (the SDK `Protein` groups) on a `QtConformationSnapshot`, off the
  `QtProtein` topology spine. T2 preserved end-to-end.
- **Two linked layers:** dense H5 animates; sparse per-frame NPY snapshots
  give full-fidelity detail when parked, joined by frame index.
- **Buffering:** lazy load on park + bounded LRU; each snapshot loads
  complete (incl. the AIMNet2 256-d embedding); the LRU bound keeps it safe
  on adviser hardware.
- **No-strings:** a generated C++ field catalog (`QtFieldCatalog.gen.h`,
  from `_catalog.py`) resolves each NPY stem → typed `FieldKind` at the load
  boundary; runtime is enums + typed objects only.
- **Names are selectable projections** (AMBER/IUPAC/BMRB × Verbatim/Derived),
  never identity.
- **Scale:** monolithic H5, eager, first pass (adviser H5s = 15 ns; µs
  "monsters" run on 128 GB batcave/spark/Strix-Halo). µs-tiering deferred.
- **Deployment split:** portable adviser build (display + movies, no embedded
  Python) vs capable build (128 GB; embedded Python + GNN + streaming,
  build-optioned). Python mechanism (embedded CPython vs out-of-process
  streaming reusing the SDK) **deferred** — the observability + serialization
  seam decouples it from the model.
- **Run shapes — trajectory vs single-conformation (decided 2026-05-26):** the
  reader stays a trajectory reader, but a single pose (`--orca`/`--mutant`/`--pdb`
  — no H5, just sidecar+manifest+flat NPYs) is a single-conformation **sibling**
  of the H5-trajectory under a shared Qt base class ("subclass, don't hack" — NOT
  a faked dummy frame). Orca is a sparse per-frame overlay: pre-extracted
  `orca_*.npy` via `QtOrcaGroup` (mutants/poses), plus a **3rd linked layer**
  parsing raw `_nmr.out` for DFT-on-trajectory (positional matchup, provably
  stricter than the library's element-only check). Full spec →
  `notes/SINGLE_POSE_AND_ORCA_DESIGN_2026-05-26.md`.

## Progress

- ✅ **Topology-spine fidelity audit** (read-only agent + codex brief) —
  faithful; 4 low-severity fixes applied (stale comments + one UNK field).
- ✅ **Field catalog** — `scripts/gen_field_catalog.py` (AST-based, faithful,
  stdlib-only) → `QtFieldCatalog.gen.h` (109 fields / 28 groups / 10 axes);
  compiles, spot-checked against `_catalog.py`.
- ✅ **Snapshot + storage** — `QtConformationSnapshot` (per-`FieldKind`
  `NpyColumn` store) + `QtResultBlocks` (PerRingTypeT0/T2, RingCounts).
- ✅ **Selectable name projection** (#7) — verbatim residue struct + loader
  population + `atomLabel`/`residueLabel` selectors; full build clean.
- 🔄 **Result-group mirror** (#3): eight families built + compile-clean +
  reviewed-clear — **ring-kernel** (BS/HM/PQ/Disp/RingChi), **bond/electrostatic**
  (McConnell/Coulomb/HBond), **scalar/vector** (Dssp/Sasa/WaterField/Hydration/
  WaterPolarization/Eeq/Bonded/Gromacs/Apbs), **AIMNet2** (charges/embedding/
  EFG×3/charge-response-gradient + the `AIMNet2Embedding` 256-d view),
  **PlanarGeometry** (multi-axis: per-atom pyramidalization, per-residue ω×3,
  per-aromatic-ring χ₂, per-saturated-ring pucker Q/θ), **Tripeptide** (ProCS15
  backbone lookup + Larsen-Eq-3 neighbour: shielding + residual-vec ML features
  + method-tag), **LarsenHBond** (Eq-5 four-term + water + Cβ diagnostic +
  count), **MOPAC** (3 sub-groups: Core with the Bond-axis `BondOrders` +
  Protein-axis `MopacGlobal`, + Coulomb/McConnell mirrors of their ff14SB
  siblings) and **Orca** (DFT quantum-reference total/diamagnetic/paramagnetic
  shielding; file-loaded, single-pose only). **Remaining:** the identity-group
  `RingContributions` / `RingGeometry` blocks, then the single-pose / dual-purpose
  loader / raw-`_nmr.out` architectural work (→ `SINGLE_POSE_AND_ORCA_DESIGN`).
- ⏳ **Lazy snapshot loader + bounded LRU** (#4) — load a per-atom NPY dir on
  park, **directory-agnostic** (a trajectory `frame_NNNNNN/` OR a single-pose run
  root — see `SINGLE_POSE_AND_ORCA_DESIGN_2026-05-26.md`).
  **Loader contract (accreted across the per-family reviews — the NPY header is
  truth, the catalog is a cross-check):** read BOTH shape AND dtype (`descr`)
  from each NPY header; widen float32→float64 element-wise for `aimnet2_aim` (the
  ONLY `<f4` array — a blind memcpy of f32 into the double buffer corrupts all
  256 cols + halves the stride); set `NpyColumn.cols` from the actual shape with
  the catalog `-1`→1 fixup (charges/scalars rely on `cols==1`, never 0); when the
  actual shape disagrees with the catalog `cols`, trust the NPY and log loud
  (e.g. `gromacs_energy` 43 vs catalog 42); fail loud on truncation. The existing
  `QtNpyReader` is sidecar-only (1-D structured) — the per-frame loader is new
  and is the first place a non-`<f8` array becomes load-bearing.
  **Verified non-`<f8` inventory (full fixture scan, 7 of 91 per-frame arrays):**
  `<f4` → `aimnet2_aim`; `<i4` → `element`, `residue_index`, `residue_type`,
  `larsen_hbond_count`; `|i1` → `omega_is_xpro`, `tripeptide_bb_method_tag`. All
  widen to double in `NpyColumn`; the other 84 are `<f8`. (`element` /
  `residue_index` / `residue_type` are per-frame identity dupes of the QtProtein
  spine — the loader may skip them rather than store redundant columns.)
- ⏳ **Wire snapshots as the second layer** (#5) — align `frame_NNNNNN` ↔ H5
  `/trajectory/frames/original_index`; surface detail when parked.
- ⏳ **Build + verify against the real fixture** (#6).

## Living sources (authoritative; this doc only points)

- `notes/H5_READER_REWRITE_DESIGN_2026-05-23.md` — the typed-mirror design.
- `notes/SINGLE_POSE_AND_ORCA_DESIGN_2026-05-26.md` — run shapes, the
  conformation base class, the dual-purpose loader, the Orca raw-`_nmr.out`
  positional-matchup spec + regression test.
- `../python/nmr_extract/_catalog.py` — the NPY format contract (single source of truth).
- `../GEOMETRIC_KERNEL_CATALOGUE.md` — the physics of each kernel.
- `../OBJECT_MODEL.md` / `../PATTERNS.md` — the library object model + patterns.
- `notes/REVIEW_BRIEF_spine_and_catalog_2026-05-26.md` — the codex review brief.
- `../src/FrameNpyEmitter.{h,cpp}` — what a `--trajectory` run emits per frame.
- Fixture with per-frame NPYs:
  `/shared/2026Thesis/shielding-calcsets/baselines/old-system-baseline/trajectory_1p9j_20260524T172121/`.

## Stale-doc warnings

- `notes/SCOPE.md` and `notes/H5_FIELD_GLOSSARY.md` predate the per-TR-group
  rewrite — they describe the old monolithic `analysis_file.h` schema and the
  `fileformat/` dependency the reader no longer uses. Trust the rewrite
  design + `_catalog.py`, not these, for the current format.

## Review cadence (per step)

After each build step: produce a lead-maximising codex brief
(`notes/REVIEW_BRIEF_*`) + run a read-only audit agent; **verify every
finding deterministically** against the source before acting (reviews
hallucinate ~1 high-confidence false positive per pass — flag generously,
trust nothing unverified). Only verified findings get acted on.

**Review log**

- **Ring-kernel family (2026-05-26): cleared.** Row-decode, FieldKind
  mappings, 8×5 per-type ordering, shell order, units, absent-not-faked all
  verified correct. The one high-confidence agent finding ("T1 should be
  `1o` not `1e`") was a **verified hallucination** — T1 (antisymmetric
  pseudovector) is axial → parity-even → e3nn `1e`, and a definite-parity
  rank-2 tensor is all-even (`0e+1e+2e`); `gen.h` (mirroring `_catalog.py`)
  is correct, not changed. It *did* surface a real **library-side**
  issue: SDK `_tensors.py:39` declares `1x0e+1x1o+1x2e` (parity-inconsistent
  `1o`) — flagged for the team; the reader does not inherit it. Loader (#4)
  carries the actionable follow-ups (column-shape validation, `cols -1→1`,
  bounds), since the views are thin and bounds-check-free by design.
- **Bond/electrostatic family (2026-05-26): cleared.** Zero decode/column-order
  defects — McConnell-category 5-enum, 5×5 `category_T2`, all scalar-block
  field orders, `QtEfg` 5-comp T2, and every FieldKind mapping verified
  three-way (reader ↔ `_tensors.py`/`_types.py` ↔ the `*Result.cpp` writers).
  Two verified fixes: `CoulombScalars::E_backbone_frac` comment corrected (it's
  a SIGNED projection V/Å, not a fraction — `CoulombResult.cpp:287`); field
  `inv_r3` → `inv_d3` to match the SDK/writer name. **(Retro-note 2026-05-26:**
  the MOPAC review later found this family's shared `McConnellScalars` left the
  `nearest_CO/CN_dist` `99.0` no-data sentinel undocumented — a real gap this
  pass missed; fixed then with `hasNearestCO/CN()` guards. See the MOPAC entry.)
- **Scalar/vector family (2026-05-26): cleared.** All 9 array-families (DSSP,
  SASA, water-field, hydration, water-polarization, EEQ, bonded, GROMACS, APBS)
  decode correctly — verified against the WRITERS (writer-is-definitive rule,
  `feedback_writer_is_definitive` memory); a read-only audit agent independently
  re-derived the headline findings. Two writer-vs-catalog drifts caught by
  applying the rule: **(1)** `gromacs_energy` emits **43** columns (fixture NPY
  is `(1,43)`) but `_catalog.py`/`gen.h` say 42 — off-by-one. The GromacsEnergy
  block decodes 43; the loader (#4) MUST size `NpyColumn.cols` from the actual
  NPY header (not the catalog) and log the mismatch. **LIBRARY FLAG:** fix
  `_catalog.py` gromacs_energy 42→43 and regenerate `gen.h`. **(2)** `dssp_chi`
  writes 0.0 (not NaN) for a missing angle though the writer's own comment says
  NaN — code wins. One reader doc fix: `DsspHBondEnergy`'s acceptor/donor comment
  made a confident directional claim that contradicted the library and isn't
  documented in `dssp.hpp` → rewritten to honest provenance + footgun note.
  Cosmetic: dropped QtDsspGroup's unused `#include "Types.h"`.
- **AIMNet2 family (2026-05-26): cleared.** All 7 arrays (charges, aim,
  efg/_aromatic/_backbone, charge_response_gradient + _scalar) decode correctly —
  verified against the writers; the audit agent independently re-derived the
  float32 embedding via byte-size arithmetic (866432 = 128 + 846·256·**4**) and
  the EFG T0=T1=0 structural zeros. EFG total/aromatic/backbone map with no swap;
  the `AIMNet2Embedding` 256-d view is non-owning with a stated lifetime; the
  charge-response gradient is correctly NOT labelled a polarisability. No reader
  change needed. The one HIGH item is downstream + captured in the #4 loader
  contract above: `aimnet2_aim` is the first load-bearing float32 array, so the
  loader MUST read `descr` and widen — never memcpy f32 into the double buffer.
- **PlanarGeometry family (2026-05-26): cleared.** All 7 arrays correct — the
  first MULTI-AXIS group (per-atom pyramidalization, per-residue ω×3, per-
  aromatic-ring χ₂, per-saturated-ring pucker Q/θ); the four distinct fixture
  lengths (846/54/15/1) are a natural axis test. Audit independently confirmed
  the aromatic-first ring-ordinal↔global mapping (vs `TopologySidecar.cpp`) and
  PROVED the writer-over-catalog ω rule: fixture residue 7 is X→Pro with
  `omega_actual=3.027` (real value, NOT NaN) — the catalog's "NaN at X-Pro"
  (`_catalog.py:397`) is stale; the reader follows the writer. `omega_is_xpro` is
  int8. No reader change. Consolidated the complete verified non-`<f8` dtype
  inventory into the #4 loader contract above (7 of 91: 1 `<f4`, 4 `<i4`, 2 `|i1`).
- **Tripeptide + LarsenHBond family (2026-05-26): cleared.** The reference-
  shielding pair — Tripeptide (ProCS15 backbone DB lookup + Larsen 2015 Eq 3
  neighbour) and LarsenHBond (Eq 5 four-term + water + Cβ diagnostic). 15 arrays,
  all decode correct: shieldings 9-col SphericalTensor, residuals 3-col Vec3,
  scalars/tags 1-col; no 1p/2p or HB/HαB transposition — the audit independently
  re-derived the sum identity `total == 1pHB+2pHB+1pHaB+2pHaB` to 5e-15 with the
  Cβ diagnostic correctly EXCLUDED. **Headline verified-fact: the two families
  use OPPOSITE per-atom sentinels** — Tripeptide writes **NaN** where unmatched
  (38 bb atoms; 474/487/488 neighbour sum/prev/next), LarsenHBond writes **0.0**
  (PackFull9 unconditionally; 657 all-zero == 657 `count==0`, zero NaN anywhere).
  `method_tag` {0 NoMatch, 1 OPBE, 2 ORCA-PBE} — the tag==0 set IS the 38 NaN
  rows. Three audit findings, all verified REAL against fixture+writer and fixed
  (comment-only): (1) the ML note now includes the water term (catalog
  `is_feature=true`, not just the 4 decomposition terms); (2) `count()` reworded —
  it's the Table-2 contribution count, NOT all geometric candidates
  (`LarsenHBondShieldingResult.cpp:652`); (3) `hasContribution()` now warns it
  gates the `shielding` TENSOR only — 3 fixture atoms have `count==0` yet
  `waterTerm==2.07`. No decode defect. One new typed block (`TripeptideMethodTag`
  enum); everything else reuses `UnpackSphericalTensor` + `Vec3` + scalars.
- **MOPAC family (2026-05-26): cleared.** Three sub-groups (`QtMopacCoreGroup`,
  `QtMopacCoulombGroup`, `QtMopacMcConnellGroup`) + three new blocks
  (`MopacScalars`, `MopacBondOrder`, `MopacGlobal`); 12 arrays, all decode
  correct. The audit independently confirmed the Coulomb + McConnell groups are
  EXACT structural mirrors of their ff14SB siblings (no field/category swap;
  category order backbone/sidechain/aromatic/CO/CN == `McConnellCategory`) and
  re-derived the fixture invariants: `mopac_scalars[:,0] == mopac_charges`
  (duplicate-charge), `mopac_bond_orders` atomA<atomB for all 896 rows.
  **This family introduces the Bond and Protein axes** to the group views:
  `mopac_bond_orders` is a (B,3) Bond-axis array whose B (896) is MOPAC's OWN
  unique-pair count — NOT the topology bond count, ordinal arbitrary (hash
  order), join via atomA/atomB; `mopac_global` is a (4,) 1-D Protein-axis row
  (loader must map → 1×4, like gromacs_energy). Units pinned to the writer:
  charge e, valency Wiberg, HoF kcal/mol (`MopacResult.cpp:435`), dipole Debye
  (`:212`). Four findings, all verified + fixed: **(F1, the keeper)**
  `McConnellScalars.nearest_CO/CN_dist` carry `NO_DATA_SENTINEL = 99.0`
  (`PhysicalConstants.h:88`) when no partner is within cutoff — a real gap
  MISSED in the Step-2 clearance (shared block), surfaced only by MOPAC's reuse;
  fixed on the block with a sentinel doc + `hasNearestCO()`/`hasNearestCN()`
  guards (mirrors `HydrationShell::hasNearestIon()`), so BOTH McConnell groups
  benefit. (2) `bondOrderCount()` returns 0 for absent AND ran-but-empty (writer
  emits explicit (0,3)) — comment now points to `snapshot.has()`. (3)
  `MopacGlobal.dipole` is origin-dependent for the net-−1 fixture — caveat added.
  (4) dropped QtMopacCoreGroup's unused `Types.h`. No decode defect; F2/F3
  (loader-contingent reshape / cols=-1→1) are the known #4 loader plan.
- **Orca family (2026-05-26): cleared.** `QtOrcaGroup` — 3 accessors
  (total/diamagnetic/paramagnetic), the DFT quantum-reference shielding (ORCA
  r²SCAN/def2-SVP, ppm; the calibration target, NOT a kernel), each a 9-col
  SphericalTensor reusing `UnpackSphericalTensor`; no new block. The audit
  independently verified the field mapping (no total↔dia↔para swap), the
  `UnpackSphericalTensor`↔`PackFull9` byte-inverse, and — against the fixture
  bytes — the Ramsey identity `total == dia + para` (max residual 1.63e-3 ppm,
  print precision); and confirmed by grep that no `Trajectory*`/`Run*`/`*Frame*`
  source computes Orca, so the "never in `--trajectory`" framing holds. One
  low-sev finding, verified + fixed: the `_nmr.out` is an OPTIONAL mode-3/4 input
  (CLAUDE.md 5-mode spec), so a mode-3 run without one has no orca_* — comment
  tightened (the nullopt path was already correct). File-loaded against the
  single-pose baseline (`orca-single_*`, (732,9) `<f8`), absent from the
  trajectory fixture by design.

## Deferred / open

- Python interface mechanism (embedded vs out-of-process) — decide when the
  observability seam lands.
- µs-scale overview tier — not needed for the first pass (monolithic H5).
- Consumer layers (movies + relationship animation, event/scripting, GNN
  export, network streaming) — layer onto the formal model once complete.
