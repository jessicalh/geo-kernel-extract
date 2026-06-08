# EFG / charge — spec-vs-code reality-check (2026-06-07)

**Provenance:** codex gpt-5.5 xhigh, read-only, full source walk (~251k tokens), file:line
throughout. Brief: `/tmp/efg_code_check_brief.txt`. **Run FIRST, before the structured grounding** —
a deliberate reorder of the McConnell cadence: EFG is *confirm-and-refine* (the code already computes
the physics), so the build substrate is written against verified code, not the first-stage grounding's
secondhand read. The grounding (`efg_grounding_agent1.md`) was treated as claims to check, not citations
to repeat. Verdict: mostly confirmed; one stale spot (AIMNet2).

## Verified facts (file:line)

**1. All-pairs vs cutoff — CONFIRMED, no documented reason.**
- FF14SB `CoulombResult`: cutoff + spatial index, default **20 Å** (`CoulombResult.cpp:113,137,144`;
  `CalculatorConfig.cpp:53`; `data/calculator_params.toml:34`).
- MOPAC `MopacCoulombResult`: **full all-pairs N²** (`MopacCoulombResult.cpp:128`); the file states the
  difference out loud ("source set: all pairs (full N²) vs AtomsWithinRadius cutoff", line 41).
- **Intent: the code is *aware* of the difference; there is NO physics rationale, no config switch, no
  comment justifying** why MOPAC is all-pairs. The "reason" is folklore, not documented.
- Range is applied **uniformly to both E and EFG** within each path (`CoulombResult.cpp:162`,
  `MopacCoulombResult.cpp:145`); the code does NOT encode an E-needs-all-pairs / EFG-fine-with-cutoff
  distinction.
- Harmonization seam exists: MOPAC already depends on `SpatialIndexResult` (`MopacCoulombResult.cpp:18`)
  so it could take a cutoff; FF14SB could be opened up (no "infinite" sentinel today).

**2. Source inventory — five EFG sources; the AIMNet2-primary line is stale.**
- FF14SB (`CoulombResult`): E + EFG from `partial_charge`; emits
  `coulomb_{shielding,E,efg_backbone,efg_aromatic,scalars}` (`:153,:395`).
- MOPAC (`MopacCoulombResult`): E + EFG from `mopac_charge`; emits `mopac_coulomb_*` (`:135,:329`).
- AIMNet2 (`AIMNet2Result`): **EFG only — NO E-field path** (`:271,:326,:385`); emits
  `aimnet2_{efg,efg_aromatic,efg_backbone}` (+ `charges`, `aim`).
- APBS (`ApbsFieldResult`): E + EFG from PB grid central differences (`:72,:86`); emits `apbs_E`, `apbs_efg`.
- Water (`WaterFieldResult`): E + EFG from explicit water; `water_efield`, `water_efg` + first-shell.
- So `charge_efg.md`'s "AIMNet2-primary" is partly reflected (EFG) but stale: the `1o` field is
  FF14SB/MOPAC/APBS/water, **not** AIMNet2.

**3. Naming + metadata thin spots.**
- `coulomb_shielding.npy` = `Decompose(EFG_total)` — bare EFG, not calibrated shielding
  (`CoulombResult.cpp:310,378`). Same for `mopac_coulomb_shielding` (`:256,:312`; trajectory header
  says "bare EFG, no gamma"). `pq_shielding.npy` = the geometric π-quad kernel.
- NPYs carry shape/dtype only; the topology manifest attaches only McConnell feature metadata
  (`TopologySidecar.cpp:620`). SDK catalog has units/parity/irreps but SDK-side.
- Catalog marks E-fields `parity="odd"` yet `irreps="1e"` not `1o` (`_catalog.py:369,402,494`) —
  metadata inconsistency (the parked **parity pass**), not a C++ physics bug.

**4. Computed-but-not-exported — true and understated.**
- Sidechain EFG computed + projected, explicitly **"not stored or written"** (FF14SB
  `CoulombResult.cpp:181,206`; MOPAC `MopacCoulombResult.cpp:161,184`). `ConformationAtom` has
  total/backbone/aromatic EFG fields, no sidechain field.
- Partitioned **E-fields** (sidechain/aromatic/backbone) computed internally but **only total E
  emitted**. APBS-delta fields stored-not-emitted (`CoulombResult.cpp:302`).
- AIMNet2 computes total/backbone/aromatic only (sidechain = total minus partitions).
- Re-extraction is cheap; a generous-separable emit needs these channels turned on.

**5. π-quad vs EFG — NOT a simple double-count.**
- `PiQuadrupoleResult` does **not read partial charges**; it builds a point axial quadrupole from
  ring-center / ring-normal / position (Stone's T-tensor, pinned per-ring-type scale; the physical
  `−Θ/2` prefactor left to `Q_type`) (`PiQuadrupoleResult.cpp:51,25`; `.h:3`).
- It overlaps aromatic EFG in electrostatic content but is a **distinct geometric basis** that may carry
  π-cloud anisotropy point charges miss.
- **Verdict: decide by ablation** (aromatic-EFG present/absent × π-quad present/absent), not by dropping
  it as a duplicate on code grounds. The "leave it out if simple" lever does **not** pull cleanly.

**6. FF14SB role — present everywhere via FullFat (NOT retired in practice).**
- The code comments say vacuum Coulomb is "retired from production" and skipped when APBS is active
  (`RunConfiguration.cpp:126`; `OperationRunner.h:53`), re-enabled only in FullFat
  (`RunConfiguration.cpp:294`). **But FullFat is now the universal run shape** (MOPAC runs per-frame on
  everything), so FF14SB-Coulomb is **de-facto always-on** — the "retired" comments are stale relative to
  the current cadence. **FF14SB stays; it is not marginal, and it is not EFG's job to remove or change it.**
- It keeps its **cutoff (20 Å)** — do NOT harmonize it up to MOPAC's all-pairs. The two ranging choices
  coexist; that difference is a caveat of the MOPAC-vs-FF14SB diagnostic, **not an EFG problem to solve**.
- It drives the MOPAC-vs-FF14SB reconciliation diagnostic (signed T2 cosine between the two bare EFG
  kernels, `MopacVsFf14SbReconciliationTrajectoryResult.cpp:56,66,139`), and APBS consumes the assigned
  charges + PB radii as its input (`ApbsFieldResult.cpp:145`). So FF14SB is load-bearing both as a field
  source and as APBS's charge assignment — keep it as-is.

## Still-open design questions (empirical, not code-resolvable)
- The winning/default source (the emergent fork — **map, don't decide**).
- The all-pairs/cutoff range difference is a caveat of the MOPAC-vs-FF14SB diagnostic, **left as-is** —
  FF14SB keeps its cutoff, MOPAC its all-pairs; neither is made to match the other, and harmonizing is
  **not an EFG deliverable**.
- Add an AIMNet2 E-field path or not.
- Partition aromatic EFG vs π-quad — by **ablation**, not a code-grounds drop.

## What this feeds
The EFG structured grounding writes against these facts: sources kept separate (fork stays emergent);
the range difference left as-is (FF14SB keeps its cutoff — not an EFG deliverable); the sidechain + partitioned-E dropped channels turned into separable emits
(cheap re-extract); the π-quad partition flagged as an ablation; the `coulomb_shielding` name fixed; the
catalog `1e`/`1o` parity inconsistency folded into the parity pass.
