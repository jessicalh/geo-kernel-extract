# Opus adversarial review — SPEC_INTERP_1P9J_SCOPING_2026-06-04

> **Historical review record — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Reviewer: opus, adversarial-but-fair. Scope: docs-only. I did not run, build, fit,
or touch git. I read the spec, its brief, the e3nn exemplars + protocol + fitter,
STATE/NOW, the two postmortems, and the substrate on disk (read-only).

**Bottom line: the GO is honest and the NO-GO is correctly scoped.** The spec is
unusually disciplined about its own feasibility. Three sharpenings below close a
soft spot in the headline number, a label on the re-emission cost, and one
flatter-risk panel. None is a blocker.

---

## 1. Feasibility honesty — is the GO real? YES, with one cost mis-labelled

What I verified actually exists (so the "very close" claim is real, not glossed):

- **Model spine** — `analysis/equiv_t2_e3nn.py` (311 lines) is a genuine e3nn
  source-sum 2e model: `o3.spherical_harmonics`, `o3.TensorProduct` /
  `FullyConnectedTensorProduct` for the `1o x 1o -> 2e` cross, invariant radial MLP
  gate, `index_add_` scatter-pool, `change_of_basis` target map, no numpy SH
  fallback. `equiv_t2_backbone_e3nn.py` (548 lines) is a real heterogeneous
  three-source-kind extension. The spec's "extend, don't author" is feasible: these
  are lean and already the right class.
- **Clean protocol** — `analysis/e3nn_protocol.py` is real: blocked contiguous
  held-out block, `PURGE_NEIGHBOUR_FRAMES=1`, train-only `centred_by_train_atom`,
  train-only `normalised_by_train_rows`, explicit `g_tr`/`g_te`. Both exemplars
  were wired to it (`POSTMORTEM_E3NN_PROTOCOL_FIX`).
- **Basis map** — `change_of_basis.get_C()` is frozen, `|C^T C - I| ~ 1e-16`. The
  spec's preflight `< 1e-12` gate (line 351) is satisfied with margin.
- **Fitter not regrown** — `allatom_fit_common.py` (2431 lines) is the ridge edge
  for the *aggregate* substrate; it already carries the within-frame split, purge,
  alpha-grid, and disk gates (`MIN_DISK_FREE_BYTES`, `MAX_REDISCOVER_BYTES`). The
  spec correctly steers Track A to the e3nn exemplar, not this pile (Step 1, line
  378: "Do not add a second scalar model. Do not add a Python source aggregator.").
- **Disk** — `/tmp/rediscover-runs` is 6.7 GiB of 15; `/tmp` has 160 GiB free. The
  spec's gates (>= 20 GiB free, < 15 GiB in rediscover-runs) hold.

**The one real hidden cost the spec under-labels: Track A's source input does not
exist on disk and must be re-emitted by the C++ producer.** I searched: there is no
`ring_current_sources.csv`, no `*_target_local_T2.npy`, no `efg_aggregated.csv`/
`efg_*T2.npy`, and the historical `/tmp/rdc-composed`, `/tmp/rdc-efg`,
`/tmp/rediscover-out-v2` directories are all gone. The spec DOES flag this — line
36 ("did not find a current `ring_current_sources.csv`... The build must preflight
this") and the Step-0 gate at lines 352-364 ("source-mode input files exist, or
stop for lead approval to re-emit"). Good. But the prose framing repeatedly calls
the delta "packaging" (lines 153-160, 178-184: "graph/rendering: small
analysis-output delta"). For the **primary** Track A that is misleading: the gating
path is a **producer re-emission run**, not packaging. Re-emission is explicitly
not Python physics — it is the sanctioned C++-emits-the-source-vectors step
(`analysis/PATTERNS.md` rule 6) — but it is a build action with its own gate
(byte-parity of the existing oracle, no library change in service of UI, extractor
SACRED). **Sharpening 1** makes this a named milestone, not a footnote, so the lead
green-lights the re-emission with eyes open rather than discovering it mid-loop.

Net: GO is honest. "Very close" is true for the *model and protocol*; it is "one
gated producer re-emission + a T0 head + a figure" for the *runnable Track-A graph*.
That is still close — but say it in those words.

## 2. NO-GO scope — correctly identified, correctly bounded

The NO-GO is exactly right and I confirmed its premise on disk. The live Build1
substrate is **aggregate-only**: every feature NPY is `(558360, N)` 2-D, and
558360 = 846 atoms x 660 frames, one row per (atom, frame). The columns are
summed/counted per atom (`ring_n`, `charge_n`, `bond_n`, `mopac_welford_mean_charge`,
scalar `*_shielding`); `features_ring_paths` is a fixed-width `(558360, 226)`
per-atom encoding, NOT a variable-length per-source row set. There is no
`disp_local`, no `source_normal_local` anywhere in it. So the spec's core claim —
that the source-sum e3nn cannot be fed from this substrate because the per-source
geometry was summed away (lines 116, 188-190) — is literally true, not cautious
hand-waving.

Crucially, the NO-GO is bounded the right way: it does **not** demand the forbidden
terabyte all-atom per-source materialization (558K atoms x all sources). It scopes
the source path to the **dominance-clean exemplar subsets** the exemplars already
target — through-space ring sources (`is_self_or_bonded == 0`, ~5-7 coupled
aromatic H) and the broad-backbone heterogeneous neighbourhood — which are small,
bounded, and the same rows the existing `equiv_t2_*` scripts consume. That is the
right line between "v0 advisor graph" and "the full Stage-3 chewer." The NO-GO
correctly says the latter needs the per-source carrier / chewer engineering and is
Stage-3 build work (lines 11, 188-190, 519).

One precision note for the build: the spec says the source path is "bounded" but
never states the expected row magnitude. Add the bound explicitly in the
re-emission gate (ring through-space ~10^4-10^5 source rows over ~5-7 atoms; broad
backbone ~6 strata) so a reviewer can see at a glance it is the exemplar subset and
not a creeping all-atom dump. Covered in **Sharpening 1**.

## 3. Overclaim guard — well held, with ONE soft spot in the headline

The within-axis / direction-not-destination discipline is held almost everywhere
and is genuinely good: the honest-reading block (lines 51-58), the forbidden list
(lines 224-231: no `trajectory.h5`, no recompute, no Python spatial graph, no
row-selection-by-residual, **no old LOAO numbers as between evidence**), the SETI
avoid-list (lines 472-479: no "validated/transferable/between recovery/BMRB-ready/
process model/solves Stage 3"), and the mandatory caption (lines 298-300). The
between-axis retraction is correctly load-bearing (lines 131-139) and matches
`POSTMORTEM_TRUE_LOAO` (true-LOAO charge 0.036, ring -1.04, unified -105) and NOW.md
("Do NOT quote any between/LOAO number"). Held-atoms-out is correctly demoted to a
within-modulation audit, never the gate (lines 243, 391).

**The soft spot: the headline "0.466 / 0.757" is doing argumentative work it has not
earned in its clean form.** The spec leans on "frame-split T2 R2 about 0.466 and |T2|
r about 0.757" as the reason the GO is strong (lines 9, 72, 125, 515) — but those are
the **pre-protocol-fix** numbers. `POSTMORTEM_E3NN_PROTOCOL_FIX` line 13 and NOW.md
lines 65-68 are explicit: the e3nn re-run under the clean blocked/purged protocol is
**HELD for the lead**, and "the clean-vs-leaky verdict... pends the HELD e3nn re-run."
The old 0.466 came from the all-group-demean + unpurged-random-split path that audit
issue #3 flagged as possibly leaky. The spec is half-aware of this — it does say "a
clean-protocol re-run may move the number; a large collapse must be investigated"
(line 396) and "expect a small metric move" (line 84) — but it still uses 0.466 as a
GO justification and as the sanity anchor (line 395: "T2 R2 should be near the
recorded 0.466"). That is a quiet circularity: justifying the GO with a number the
same spec says is not yet clean.

This does **not** sink the GO — the model/protocol/substrate-path existence is the
real basis for GO, and that stands independently. But two tightenings are needed
(**Sharpening 2**): (a) label 0.466/0.757 in-line every time as "pre-protocol-fix,
leaky-suspect; clean re-run held," and rest the GO on existence-of-machinery, not on
that figure; (b) anchor the interpretation gate (T2 R2 >= 0.25 = "direction") to the
**clean-protocol** number the build itself produces, not to 0.466 — otherwise a clean
re-run that lands at, say, 0.2 could be silently graded against a leaky 0.466
expectation. The "crude graph is still a win if held-out T2 recovery is positive and
the caption is honest" (line 397) is the correct fallback and should be the load-bearing
sentence, not the 0.466 anchor.

Minor: the spec is internally honest that Track B (aggregate-shadow e3nn) is weaker,
but "aggregate-shadow e3nn v0; emitted source sums, not full source chewer" (line 170)
should be the **on-figure** label if Track B is ever graphed, not just a metrics-JSON
field — an advisor reading the PNG must see it is the fallback. The spec implies this
but should state it as a hard caption rule for the Track-B branch.

## 4. The graph — right artifact, honestly captioned, one flatter-risk

The 2x2 is the right advisor-facing object and the panel choices are sound:

- Panel A (5-component held-out T2 scatter, train-atom de-meaned, y=x, R2 annotated)
  — the correct primary; T2 not scalar-collapsed (`feedback_t2_sacred` held).
- Panel B (|T2| modulation recovery, rotation-invariant r) — the right
  advisor-readable invariant; honest that it is magnitude.
- Panel D (caption + flow sketch) — the mandatory within-axis/geometry-sampler/
  correlate-not-match/direction-not-destination text is present and correct.

**Flatter-risk on Panel C (sigma_iso/T0 companion).** The spec already knows the
trap — "absolute sigma_iso can be dominated by per-atom chemical baseline" (line 212)
and "caption says the T0 head is companion, not the T2 thesis target" (line 286). But
it leaves the baseline-restored absolute-sigma_iso inset *optional* and does not make
the modulation-R2 mandatory **on that panel**. An advisor's eye will land on the
tightest scatter in the figure; a baseline-restored absolute-sigma_iso panel will be
near-diagonal almost regardless of model quality because the per-atom chemical
baseline (tens of ppm of spread across atom types) swamps the few-ppm modulation. If
that panel is shown at all, the modulation R2 must be printed **on it, larger than the
restored-R2**, and the restored panel labelled "baseline-dominated; not a model-skill
metric." **Sharpening 3.** (Line 245 already requires this in general — "if a
baseline-restored panel is shown, it must also show the modulation R2" — promote it
from general rule to a hard Panel-C caption requirement, and make the restored panel
the inset, never the main C axes.)

The interpretation gate buckets (lines 410-412) are honest and well-graded against
the project's statistical-position standard. The hard graph gate (lines 400-406:
non-empty held-out rows, no leakage in split audit, finite R2, row-count match) is the
right "what proves it works." Good.

## 5. Build sketch — runnable + gated, ready for the loop

The build sketch is genuinely loop-ready: fresh explicit out dir (no glob), Step-0
preflight gates (disk, basis orthogonality, input-file existence with fail-loud,
no-H5 / no-recompute boundary audits), the reuse-the-exemplar Step 1, the
train/score Step 2 with the clean protocol, the artifact set (predictions CSV +
metrics JSON + graph + run audit), and a SETI report that is numbers-first with an
explicit caveat label. It reuses the exemplar and does not hand-roll a second model.
The boundary list (lines 224-231) and the proof-it-worked list (lines 482-492) are
the right gates.

Two build-readiness gaps to close before firing (both fold into the sharpenings):

- The re-emission, if triggered, needs its own explicit sub-gate: byte-parity of the
  existing emitted oracle for unchanged columns, additive-only, broad-specific so the
  ring/mc oracle stays intact, and the extractor-SACRED rule honoured (the producer
  is not modified in service of a UI/analysis feature — re-emitting an existing source
  sink is the allowed path, a schema change is not). The spec gestures at this (line
  364) but should make it a named gate with the byte-parity check, mirroring how
  `allatom_fit_common`'s Loop-1 emit was gated.
- The clean-protocol expectation anchor (Sharpening 2) must be in the metrics JSON:
  record `protocol = clean_blocked_purged`, `cross_split_lag1_pairs == 0`, and stamp
  the result as the *new clean baseline*, explicitly NOT compared to the leaky 0.466.

---

## Three sharpenings for the build (priority order)

1. **Re-label the Track-A delta as a gated producer re-emission, not "packaging,"
   and state its row bound.** Make "re-emit the bounded ring/backbone source substrate
   (~10^4-10^5 ring through-space rows over ~5-7 atoms; ~6 backbone strata)" a named
   milestone with a byte-parity + additive-only + extractor-SACRED sub-gate. This is
   the one real hidden cost; surface it so the lead approves it deliberately.

2. **Stop resting the GO on 0.466/0.757; anchor the gate to the clean re-run.** Tag
   those numbers "pre-protocol-fix, leaky-suspect, clean re-run held" at every
   mention; rest the GO on machinery-exists; grade the produced clean number on its
   own (>=0.25 = direction), never against the leaky 0.466; make "positive held-out
   T2 recovery + honest caption = a win" the load-bearing fallback sentence.

3. **Make Panel C anti-flatter mandatory.** Modulation R2 printed on Panel C, larger
   than any restored-R2; baseline-restored absolute sigma_iso demoted to an inset and
   labelled "baseline-dominated, not a skill metric"; never the main C axes. Same hard
   on-figure rule for the Track-B "aggregate-shadow, not full chewer" label if Track B
   is ever the graphed path.

---

## Verdict

- **GO honest?** Yes. Model spine, clean protocol, frozen basis map, and the
  exemplar all exist and are lean; the delta is small in code. The one real hidden
  cost (Track-A source CSVs are gone -> a gated producer re-emission) is flagged and
  gated but mis-labelled as "packaging" — fix the label, not the verdict.
- **NO-GO correctly scoped?** Yes, and confirmed on disk: Build1 is strictly
  `(558360, N)` per-(atom,frame) summed features, no per-source `disp_local`/
  `source_normal_local`; the full all-atom chewer genuinely needs the per-source
  carrier. The NO-GO bounds the source path to the dominance-clean exemplar subsets,
  not the forbidden all-atom terabyte. Correct.
- **Overclaim held?** Mostly yes (within-axis, correlate-not-match,
  direction-not-destination, between-axis retraction all held). One soft spot: the
  0.466/0.757 headline is a pre-protocol-fix (leaky-suspect) number used as a GO
  justification while the clean re-run is HELD. Tighten per Sharpening 2.
- **Graph right?** Yes — 2x2 (held-out T2 scatter / |T2| recovery / sigma_iso-T0
  companion / direction caption) is the right honestly-captioned advisor artifact.
  Sole risk: Panel C baseline-restored sigma_iso can flatter; make modulation-R2
  mandatory and the restored panel an inset (Sharpening 3).
