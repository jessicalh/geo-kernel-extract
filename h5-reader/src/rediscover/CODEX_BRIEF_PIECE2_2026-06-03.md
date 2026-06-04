# Codex brief — Piece 2: all-atoms joint fit + partition-by-condition

> **Historical run brief — not current truth (trued 2026-06-04).** Session
> provenance only; current rediscover truth is `NOW.md` and corrected `STATE.md`.

Status: **landed** at `5b5525b45f363e7c4a05cfff6a02a19d2d8c15ab`; retained as the execution brief.

You own the grind; the lead and the project owner vet and judge. This is **layer 2
of the thesis reporting arc** (equation calibrations) and the **engine of the
reframed law-example hunter**: fit across ALL atoms so the between/static axis
becomes determinable (846 atoms ≫ features — the per-54-atom-stratum fits couldn't
resolve it), then partition by input-side conditions so the favourable cases
**emerge** rather than being hunted. Keep the run focused (one fit+partition pass);
there is a codex context-image platform bug on long runs — if you crash, that is
the platform bug, the lead re-fires.

## Authoritative spec

Implement, EXACTLY, three sections of `src/rediscover/ALLATOM_FIT_SPEC_2026-06-03.md`:
"**Piece 2: All-Atoms Joint Fit**", "**Partition By Condition**", and
"**Fit-Stage Checks**". Read them in full. The spec is authoritative for inputs,
join contract, target, tiers, ridge/preprocessing, held-out axes, score columns,
condition families, bins, response-curve columns, favourable-case eligibility, and
artifacts. This brief adds only the owner's overrides, the disciplines, and the
report/commit contract — it does not restate the spec.

## Input (the charge-complete substrate from Loop 1)

```
/tmp/rediscover-runs/2026-06-03-per-atom-substrate-charge-scalars-piece1-final
```

Committed at `4bb9a01`; gates green; `charge_complete_rows = 558360` (100%).
Confirm via its manifest: `relationship = per_atom_substrate`,
`relationship_kind = per_atom_aggregate`, shape 846 × 660. Read THIS directory
only.

## Owner's overrides (these WIN over the spec defaults)

1. **Embedding → training-only PCA, NOT the spec's full-256-d default.** Fit PCA on
   TRAINING rows only (no leakage), apply the fixed transform to held-out rows, and
   record n_components + explained-variance-ratio + the transform provenance in
   `run_audit.json`. Rationale: 256 raw dims on 846 between-atoms would re-break the
   very static axis this fit exists to make determinable. Full-256-d is allowed ONLY
   as a clearly-labelled secondary sensitivity, never the primary.
2. **Same charge-complete row set across ALL three tiers.** All of
   `classical_mechanisms_combined` / `plus_AIMNet2` / `all` fit on the identical
   charge-complete row set, so `delta_R2_vs_previous_tier` is honest apples-to-apples.
   (Here that is all 558,360 rows since coverage is 100% — still enforce the filter
   and report `charge_complete_rows` used.)
3. This IS Loop 2. Do NOT touch the C++ emit / `PerAtomSubstrate.*` / the H5 read
   path. Python edge only.

## Disciplines (hard — these are the project's law)

- **Read the emitted substrate directory ONLY.** NEVER open `trajectory.h5`, ORCA
  outputs, `/tmp/combined-mopac-layer3`, broad-backbone dirs, MOPAC dirs,
  all-atom-equivariant dirs, or any pair/source dump. The audit must state this
  explicitly. Sidecar joins by row id; `driver_modulation_by_atom` joins by
  `atom_index`; feature names/blocks come from `per_atom_substrate_column_specs.json`.
- **No Python physics recompute.** The fit reads emitted feature columns and derives
  no physics. Algebraic conditioners for PARTITIONING are allowed, but every derived
  formula must be listed in `run_audit.json`.
- **Frozen basis.** Use `change_of_basis.get_C()` as-is; do not reimplement it.
  Assert `abs(C.T @ C - I).max() < 1e-12`. Target is DFT total T2:
  `y_2e = target_T2[:, 0:5] @ get_C().T`; score the 5-component 2e vector, NOT only
  `|T2|`. Every T2 feature block in the library basis gets the SAME `get_C()`
  transform; scalars/vectors are not projected.
- **Anti-circularity (load-bearing).** Partition on INPUT-side conditions ONLY. NEVER
  define a bin using the DFT target, residuals, or fitted coefficients. Select the
  habitat by input geometry/identity/modulation; TEST recovery against DFT.
- **Held-out everywhere; report held-out R² only.** Between = leave-atoms-out over the
  846 atom-mean points (this is cheap — do NOT naively refit on 558K rows 846 times;
  use atom means per the spec). Within = held-out contiguous frame block, purge ≥1
  adjacent frame each side, report cross-split lag-1 (target 0). Report AR(1)-aware
  N_eff per T2 component + median lag-1 ρ. Jackknife-over-atoms CIs.
- **SETI — do NOT editorialize.** Produce the numbers, the response-curve SHAPES
  (monotone rise / fall / U / threshold / flat / unstable-thin), and the
  favourable-case shortlist as FACTS. Do NOT grade "law"/"not law", do NOT assert
  physical conclusions, do NOT rank by chemistry intuition. The statistical-position
  grading and the verdict are later steps and the owner's call. Use "expected
  relationship + probability + fit" language only where the spec's own table does.

## Environment

Use the existing analysis venv (`src/rediscover/analysis/venv`); ridge via sklearn
or a numpy closed form. If a required package is genuinely missing, report it rather
than improvising. Write the analysis as a committed script under
`src/rediscover/analysis/`. Write all result artifacts to a fresh drop-old results
dir, e.g. `/tmp/rediscover-runs/2026-06-03-allatom-fit-piece2/` (ephemeral; do not
commit heavy result data).

## Artifacts (per spec)

Score outputs: `allatom_fit_score_table.csv`, `allatom_fit_score_long.csv`,
`allatom_fit_report.md`, `join_coverage.csv`, `feature_blocks_used.json`,
`run_audit.json` — with the exact `allatom_fit_score_table.csv` columns the spec
lists. Partition outputs: `partition_response_curves.csv`,
`partition_response_curves_long.csv`, `partition_favorable_cases.csv`,
`partition_report.md`. Primary atom-type axis = emitted `stratum`; also report
`role` and `ff_atom_type_ord` coverage.

## Fit-stage checks — run ALL, report each

(From the spec's "Fit-Stage Checks".) Input acceptance; **no external merge** (audit
proves substrate-dir-only); basis + target (`C.TC≈I`, target shape (R,5), consistent
T2 transforms); **CV integrity** (between held-out atoms never in train; within
held-out block excluded from all preprocessing/impute/standardize/alpha/centering;
purged frames absent from both; scores are held-out; N_eff + lag-1 ρ per row);
**partition integrity** (every condition emitted or derived-from-emitted input-side;
bin defs written before ranking; no condition uses DFT target/residual/coef; thin
bins flagged + excluded from shortlist). If any check fails: STOP, report it, do not
present scores as valid.

## Commit (only if the fit-stage checks pass)

Branch `h5-reader-pysr-spike`. **NEVER merge / switch / rebase / PR / checkout
another branch.** Commit the analysis SCRIPT(s) under `analysis/` (reproducible code
+ provenance), atomic. Do NOT commit the heavy result data (it lives in the drop-old
results dir; the owner decides what becomes durable). Message: describe the
all-atoms charge-complete joint fit + partition.

## Report back — drift-assessment + postmortem (plain, no editorializing)

- What landed vs what Piece 2 asked. Any drift?
- Each fit-stage check: pass/fail.
- HEADLINE numbers: per-tier × per-stratum (N/CA/C/O/HN/HA …) **between-LOAO** and
  **within-frameblock** held-out R² with CIs; the tier deltas
  (`delta_R2_vs_previous_tier`); the AR(1) N_eff (min/median/max + median ρ).
- The **emergent favourable-case shortlist** (the partitions where recovery is
  strongest) and, factually, the **response-curve shapes per mechanism/conditioner**
  (which rise, which fall, which are flat).
- Commit hash; the results-dir path.
- Anything the spec left to your discretion that you decided.

The owner will read the numbers and steer the next loop. Keep it tight.
