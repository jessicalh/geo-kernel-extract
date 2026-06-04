# Codex brief — Build 4: dominance-completion C++ emit (the dominance handle for the law-bearing five)

> **Historical run brief — not current truth (trued 2026-06-04).** Session
> provenance only; current rediscover truth is `NOW.md` and corrected `STATE.md`.

Status: **DRAFT — pending lead plan-vet, not yet fired.**

You own the grind; the lead vets + judges + owns ALL git. This is the C++ half that
COMPLETES the dominance arc. Today the dominance/isolation handle exists for only THREE
mechanisms (ring, charge, mc); the lead's law set is FIVE. The two missing — **field
(Buckingham / FF14SB E)** and **H-bond (Larsen)** — are the two laws we most want to judge
now, and they have NO dominance handle at all. The deliverable: give every law-bearing
mechanism the same dominance/isolation/bin machinery + a CaseHunter habitat, so the
dominance-clean law fit and the navigable examples cover what is on the tin — not a wee bit.

Read `PerAtomSubstrate.cpp` (the existing `dominance` lambda ~`:2755`, `rowPairContributions`,
the Build-1 isolation block ~`:1703`, and the existing **field** + **H-bond** emit blocks) +
`CaseHunter.{h,cpp}` (the `habitat()` switch ~`:186` + the 3-mechanism loop ~`:354`) FIRST.

## SCOPE (exactly this — no scope creep)

IN: dominance/isolation/bins + CaseHunter habitat for **field** and **hbond** (ring/charge/mc
already have them). Move the dom quantile **bin-ids into C++** for all five (retire Build 3's
Python-quantile crutch). Widen the **charge** habitat to include carbons.
OUT (explicitly, this loop): dispersion, π-quadrupole, EFG, water — emitted as features, NOT
law-bearing now, NO dominance handle this loop. ringχ sign-fix is ANALYSIS-side (next stage),
NOT here. The fitter decomposition is the NEXT stage, NOT here.

## HARD RULES

- **C++ spine ONLY. Zero Python computes/derives any emitted value.** Mirror the existing
  dominance/isolation blocks (typed, no string dispatch, ArraySpec per output, present-flags,
  units/mechanism metadata). Producer/extractor/CMake/ctest untouched beyond this rediscover
  code. Read existing H5/NPY + the live model/indexes only; NO re-extraction.
- **Write NO new Python in this loop — don't grow the Python pile.** No analysis scripts, no
  `.py` helpers, no edits to the fitter. BUILD 4 is pure C++; ALL fitting/analysis is stage 2.
- **NO GIT. You do not run any git command (no commit/add/reset/rebase/checkout).** Leave the
  gated-green changes uncommitted; report the gate results + a tight diff-stat. The LEAD
  reviews the gates and commits. Branch `h5-reader-pysr-spike` — NEVER merge/switch/rebase/PR.
- **DISK GUARD:** `df` before any emit; if free < 20 G, ABORT + report. Drop-old by EXPLICIT
  FULL PATH only (NEVER a regex/glob — name each dir you remove, print it). Total rediscover
  output < 15 G. NEVER write/delete under `/shared`.
- **Lean / keyed (atom,frame) only.** No per-source/pairwise dump. `query_frame_slots` stays 1
  (no boa-constrictor). Each dominance/gap value is a lean per-(atom,frame) scalar reduction.
- **FAIL-LOUD-AND-LOCATE (load-bearing — this is the design risk).** For field and hbond, FIRST
  LOCATE the per-source contribution structure in the resident model (how the existing field /
  hbond emit + `rowPairContributions` are built). Define `dominant_fraction` /
  `gap_to_2nd_*_r` CONSISTENTLY with ring/charge/mc:
  `dominant_fraction = max|source contribution| / Σ|source contribution|` over the mechanism's
  source set; `gap_to_2nd = r(2nd-nearest source) − r(nearest source)`. If a mechanism's
  per-source contributions are genuinely NOT available to dominate over (e.g. field is an
  APBS-grid pre-sum with no decomposable per-source term), DO NOT fabricate — STOP that
  mechanism, emit present=0/NaN, and report the obstacle + where you looked in the postmortem.
  "Absent under the obvious name" ≠ absent — locate first (read the existing emit, not glob).
- **Anti-circular:** the new dominance/gap/bins are INPUT-side only; never read the DFT
  target/residual/coef to compute them. The Build-1Fix real selection-guard must still fire.
- **Self/same-residue exclusion** applied consistently (the Build-1Fix charge lesson): field and
  hbond source sets exclude self + same-residue where ring/charge/mc do.

## Piece A — dominance/isolation for field + hbond (mirror ring/charge/mc)

- **field — the MOPAC-Coulomb leg, NOT Amber/FF14SB.** The field signal that tracks DFT is the
  **MOPAC-Coulomb** field, far more than the FF14SB/Amber field (Stage-1 mining: APBS-EFG tiny,
  Coulomb/MOPAC EFG moderate → the field leg is MOPAC; `STATE.md` 2026-06-02 EOD). Define
  `dominant_fraction_field` / `gap_to_2nd_field_r` over the **MOPAC-Coulomb** field's per-source
  contributions — the decomposable `Σ q_mopac·r̂/r²` source set over the resident MOPAC charges +
  geometry. Make the source explicit (column name or manifest `source=mopac_coulomb`). Do **NOT**
  build an FF14SB/Amber-field dominance — that is the weak leg, out of scope this loop. APBS field
  is a grid solve, never per-source-decomposable, never the field-dominance source.
  - **LOCATE then build:** if the MOPAC-Coulomb per-source contributions already live in
    `rowPairContributions` / a charge-source cloud, dominate over them directly. If only the
    pre-summed `mopac_coulomb_E` is emitted, COMPUTE the per-source field contributions C++-side
    from the resident MOPAC charges + positions as a lean in-memory reduction (mirroring the
    charge q/r³ source build) — this is a rediscover-internal spine reduction, **NOT** a
    re-enable of the retired production `coulomb_*` NPY path (CoulombResult stays FullFat-only).
    Only if the MOPAC charges/geometry are genuinely absent for these atoms do you fail-loud:
    emit present=0/NaN + report where you looked. Do not fabricate.
- **hbond:** `dominant_fraction_hbond`, `gap_to_2nd_hbond_r` over the H-bond partner set
  (nearest partner's contribution / Σ over partners; gap = 2nd-nearest − nearest partner
  distance), using the Larsen per-class / `hbond_scalars` partner geometry already emitted.

## Piece B — C++ dominance bin-ids for all FIVE (the deferred-from-Build3 piece)

Emit per-(atom,frame) dom quantile bin-id columns for
`dominant_fraction_{ring,charge,mc,field,hbond}` — compute the quintile edges C++-side in one
pass, record edges in the manifest. Python bins by LOOKUP, never derives. Bin-ids must be
non-degenerate (multi-bin), like the charge H1 fix.

## Piece C — CaseHunter habitats for field + hbond + widen charge

- `HunterMechanism::Field` — habitat = the atoms where the field law lives (polar/charged
  heavy atoms + backbone N/HN), typed via `TypedAtomIndex`.
- `HunterMechanism::Hbond` — habitat = H-bond donors/acceptors, role-resolved (donor→HN/N,
  acceptor→O/C′), preferring the Larsen per-class assignment.
- **Widen `HunterMechanism::Charge`** habitat to include carbons (C/CA) — currently O/N/S/H;
  carbons feel charge fields too.
- Same anti-circular guard (DFT read during selection THROWS); deterministic (same thresholds →
  identical output).

## Piece D — cases manifest

Extend `equations/<mechanism>/cases_manifest.csv` to `field` and `hbond` (same navigable
addressing as ring/charge/mc).

## Gates (ALL pass before you hand back; lead commits on green)

- **Disk guard** passed (free ≥ 20 G; drop-old; total < 15 G).
- **Oracle parity:** every pre-existing NPY/CSV byte-identical (append-only); `backbone_audit`
  byte-identical; DFT target parity exact.
- **Dominance two-path cross-check for the NEW mechanisms** (field/hbond): the new
  per-(atom,frame) `dominant_fraction` equals a 1-frame named-query value (as Build 1 proved for
  ring/charge/mc, to ~5e-13).
- **New-channel coverage:** present flags 0/1; present ⇒ finite; bin-ids non-degenerate; report
  each new column's present-count + range + (if a mechanism was stopped) the located-but-absent
  reason.
- **CaseHunter:** > 0 candidates for the field AND hbond habitats; deterministic re-run
  identical; anti-circular assertion fires on a deliberate DFT leak (regression test).
- **R & uniqueness:** 846 / 660 / 558,360; dense `row_id`; unique `(atom_index, original_frame_index)`.

## Build / run / output

- **Build (unsandboxed):** `cmake --build build/linux-gcc --target h5reader_extract h5reader_rediscover_tests`.
- **Input — locate, do NOT guess/glob:** re-run on the SAME producer calcset that made the
  `…-build1` substrate; read `…-per-atom-substrate-build1/run_audit.json` for its `--run` input
  path. `h5reader_extract` reads the EXISTING producer H5/NPY (where `mopac_coulomb_*` lives) —
  it does NOT re-run `nmr_extract`, which is SACRED.
- **Output:** a fresh dir `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`. **KEEP the
  `…-build1` substrate** (the current good one; ~130 G free, plenty of room) until the LEAD
  verifies build4 — do NOT drop build1 this loop. You MAY drop-old only your OWN superseded
  intermediate dirs, by EXPLICIT FULL PATH (never a glob), printing each.
- **Run command:** `h5reader_extract --run <input from build1 audit> --out
  /tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4 --case all`.

## Resilience (bank progress against the codex context-image bug)

Implement + gate in TWO atomic chunks so a crash never loses banked work: (1) field+hbond
dominance/isolation/bins (Pieces A+B); (2) CaseHunter habitats + charge widen + manifest
(Pieces C+D). Gate each chunk before moving on. Short runs + retry if the context-image bug bites.

## OUTPUT — TINY VERDICT

Write `src/rediscover/POSTMORTEM_BUILD4.md` (≤55 lines): per-chunk gate pass/fail; the field +
hbond dominance DEFINITIONS you settled on (and any mechanism you STOPPED on fail-loud, with
where you looked); the dominance two-path agreement for field/hbond; new-column present-counts +
ranges; CaseHunter candidate counts for field/hbond; disk-guard free before/after; the run dir;
and a one-line "ready for lead commit" (NO git run). Print ONLY a ≤12-line summary + that path;
DO NOT echo diffs/tables to stdout.
