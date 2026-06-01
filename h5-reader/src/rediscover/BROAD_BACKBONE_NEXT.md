> **SUPERSEDED (2026-06-01) — historical.** The broad-backbone case is BUILT + GATED
> (commit `35f3768`) and the equivariant-T2 exemplar landed (`901d1df`). Current
> handoff: `NEXT_SESSION_PROMPT.md` + `REDISCOVERY_MAP.md` + `STATE.md`. Ignore the
> below (historical).

# Next step — the BROAD case (every backbone atom, all mechanisms)

Status: the decided next move (2026-05-31, post charge_dipole critique). Read with
`SURFACE_DESIGN.md`, `STATE.md`. Branch `h5-reader-pysr-spike` — NEVER MERGE.

**The two big ones to stand up the functional API (chosen 2026-05-31):** this
broad backbone case (**A**, source_sum, forces the frames) PLUS **efg #4** (APBS
EFG → DFT T2, all atoms) as the orthogonal companion — efg #4 exercises the
`per_atom_feature` carrier branch A never touches and the equivariant-T2 fit in
the lab frame (tumbling-safe, no new frame; reuses the proven ring-T2 machinery,
APBS data present now). Together they cover BOTH relationship_kinds and BOTH
fit-shapes — that is the API genuinely in place, not half-exercised.

## The pivot (the lesson — do NOT slide back into narrow cells)

Every cell so far is **one stratum × one source-type** (aromatic-H←rings,
HN←bonds, HN←charges). That is the "single-bond-thingy" trap. Jessica's concern:
we may have filled the code with narrow concepts and never addressed the GENERAL
target. Honest breakdown:

- **Genuinely narrow / baked:** only TWO local frames coded (HN + aromatic-H) — a
  backbone N/CA/C/O has NO frame, so "every backbone atom" can't run today; and
  the ring selector's `slots()` backend is aromatic-H-specific.
- **Actually general (API designed right, just unexercised):** `Stratum` is a
  predicate (`atomsWhere(IsBackbone)` is a one-liner); selectors are a LIST
  (heterogeneous by design); the KD backends (ring-centres / bond-midpoints /
  charge-sites — all built) serve ANY target atom; DFT covers all 846 atoms; the
  carrier (`relationship_kind`, per-relationship schema) handles a heterogeneous
  bundle.

So the narrowness is mostly in what we've EXERCISED + the unbuilt frames — NOT the
architecture. But the broad case is UNPROVEN, so we don't know it composes.

## The task — prove it by running the broad case

A SINGLE relationship, not another narrow cell:
- **Stratum** = all backbone atoms (N / CA / C / O / HN / HA, every residue) via
  `atomsWhere(IsBackbone)`.
- **Heterogeneous selectors** on each backbone atom: `[rings (KD ring-centres —
  the GENERAL backend, NOT the aromatic slots), anisotropic bonds (KD),
  charge-sites]`; emit each mechanism's summed feature + the un-summed per-source
  rows.
- **Frames for every backbone atom type** — build Cα / CO / HA / N frames (the
  conventions doc `substrate_conventions_2026-05-30.md` "Local frames per atom
  class" specifies HA / Cα / CO; define N). THIS is the forced fix for the baked
  narrowness — the broad run can't start without it.
- **Target** = DFT σ_iso (+ tensor, T1 flagged) per backbone atom.
- **Carrier**: the FIRST heterogeneous bundle — per-relationship schema; array
  payloads → NPYs (`_catalog.py`); apply the carrier fix below.

## Gate (no single oracle — broad)
Runs across ALL backbone atom types (frames resolve, no narrowness blocks);
well-formed two-kind output + manifest; the combined mechanism features explain
backbone σ_iso **per-element / per-atom-type** in the Stage-1 ballpark
(correlate-not-match; report per-atom-type R² + effective N). If it composes, the
surface delivers the GENERAL analysis (the thesis), not a pile of explorers. If it
breaks, the baked narrowness (frames / slots / carrier-under-heterogeneous-load)
is found NOW.

## Carry from the charge_dipole critique (commit 3103e73)
- **Field, not μ.** The dipole moment μ is null vs σ_iso; the derived Coulomb
  FIELD carries it (Buckingham `field_z` r=+0.46 along N–H, leave-atoms-out
  R²=+0.21). Emit the FIELD (E, bond-projected E_z) as the scalar feature; keep μ
  + per-source rows for the tensor story.
- **Cutoff sweep.** 6 Å truncates a long-range (1/r²) field — the signal is
  under-measured. Sweep 6/10/12 Å (and/or all-protein-charges); see if field_z /
  LOAO R² climb.
- **Charge-source axis:** FF14SB (static, geometry-only) now; **AIMNet2** (H5,
  per-frame → captures polarization FF14SB can't) and the **APBS field** (#3,
  PB-solved → includes implicit-solvent screening point-charges miss) are the
  natural stronger-field comparisons.
- **Carrier fix (applies to ALL relationships):** the per-source CSV bloats
  (828 MB) because the ~50-column DFT target repeats on every source row. Key the
  target ONCE per (atom,frame) (aggregated row / a target NPY); source rows carry
  only source fields + an (atom,frame) index.

## USE the resident topology — do NOT regenerate it

The body already HOLDS the topology, resident from load (Layer 0: atoms / bonds /
rings / residues + reverse indices; `TypedAtomIndex`, `RingGeometryCache`,
`ChargeStore` with FF14SB charges, the per-cloud KD trees — ALL built day 1).
Build every frame, anchor, selector, charge, and target off THOSE indexes and the
typed model (`QtProtein`/`QtResidue`/`QtRing`/`QtAtom`). Do **NOT** re-parse the
PDB / `topol.top`, re-derive connectivity/bonds/rings, name-string-scan, or
re-detect aromaticity — that work is done and lives in the resident body. The
forced new work (backbone frames) is the ONLY thing being added; it consumes the
existing topology, it does not regenerate it.

Concretely, for A + efg #4:
- **Backbone frames (N/CA/C/O/HA):** get the anchor atoms via
  `selectUnique(residue, {typed locant})` on `TypedAtomIndex` — the same
  collision-safe typed lookup the aromatic-H CG/CD2 anchor now uses — NOT a name
  string, NOT a positional index, NOT a re-parse. (identity-from-chemistry; the
  IUPAC-revert trap — `feedback_identity_from_chemistry_not_position`.)
- **Ring sources:** the `ring-centres` KD cloud + `RingGeometryCache`; do not
  re-detect rings.
- **Charge field:** read `ChargeStore` (FF14SB already loaded, typed
  resnr/order-checked); do not re-read `topol.top`.
- **efg #4 target/source:** read the APBS EFG time-series already in the H5 via
  the catalog `value(arrayId, atom, frame, …)`; it is a per-atom datum to READ,
  not to recompute.

If a needed datum is not yet in the resident body, the fix is to ADD a Layer-0
catalog entry / day-1 index for it — never an ad-hoc regeneration inside the hot
loop. "No file discovery, no regeneration": documented conventions + resident
typed indexes only.

## How to run it
Plan first (lead, cheap). Build via **codex** (`codex exec
--dangerously-bypass-approvals-and-sandbox` from the lead — full agency; Claude
Agent-tool subagents are sandbox-blocked from compile, see
`reference_subagent_build_agency`) or the lead. Settle the cell, then codex-do,
then gate + independent verify (the proven pattern). NEVER MERGE, full stop.
