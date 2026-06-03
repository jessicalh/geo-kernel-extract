# Fresh-Look Review — rediscover docs vs the units-audit corrected state (2026-06-02)

Read-only fresh-eyes audit before handoff. Branch `h5-reader-pysr-spike` (NEVER MERGE).
Checked the rediscover docs against the corrected state in `UNITS_AND_ISSUES_AUDIT.md` +
`STATE.md` (top). Banners added to genuinely-stale/misleading docs; no doc body rewritten.
Scientific verdicts reserved for the lead. This doc + the 5 banners are the only changes.

## The one corrected fact the docs must reflect

The literature-fixed de-circularisation "failure" (ring T0 γ=−11.3, McConnell T0 γ=−4.75)
is a **units/labeling bug, not a physics failure**. The test read producer **bare** kernels
(ring = unit-current Biot-Savart, `ppm_T_per_nA`; McConnell = full `M/r^3`, `Angstrom^-3`),
NOT literature-ppm-scaled predictions. γ≠1 because γ silently absorbs the omitted
`LiteratureIntensity` (ring) / Δχ×unit-prefactor (McConnell). The negative ring sign is the
**correct** diamagnetic ring-current convention (`G=−n·B·PPM_FACTOR`); the DFT target is ORCA
σ (`dft_sigma_iso`), so it is NOT a σ-vs-δ flip. Ring is likely a clean recovered law once
re-run in literature-scaled units. Full map + 5-item fix-set: `UNITS_AND_ISSUES_AUDIT.md`.

## Banners added (5 docs)

| doc | why stale/misleading | banner type |
|---|---|---|
| `analysis/LITERATURE_FIXED_DECIRC.md` | line 7 says emitted kernels + DFT "compared in matched ppm units"; the columns read are `bare_T0`/bare-T2 = producer bare kernels (`literature_fixed_decirc.py:32-46,381-403`). This is the central offender the audit flags (`UNITS_AND_ISSUES_AUDIT.md:70-99`). | SUPERSEDED |
| `analysis/BACKBONE_LAW_EVIDENCE.md` | "De-Circularising Correlations / fixed-kernel correlations against DFT T2" (lines 42-53) + "Literature-coefficient-fixed prediction remains a substrate gap" (line 116) read `*_literature_kernel_T2` sidecars that the reducer fills from H5 bare/geometric-fallback kernels, not literature ppm (`BroadBackbone.cpp:65-127,468-477`; `UNITS_AND_ISSUES_AUDIT.md:101-109`). | SUPERSEDED |
| `analysis/STATIC_ENVIRONMENT_CALIBRATION.md` | `mcconnell_T2_fixed` / `ring_current_T2_fixed` rows tagged `calibrated-to-physics`, `literature_value=1.0000`, "fixed emitted literature-kernel multiplier" — script labels bare kernels as fixed-literature coef-1 (`static_environment_calibration.py:493-523`), which the emit does not support (`UNITS_AND_ISSUES_AUDIT.md:106-109,147`). γ values (e.g. McConnell-N 25.96, ring-N −391.7) are unitful, not dimensionless-≈1. | SUPERSEDED |
| `analysis/VARIANCE_DECOMPOSITION_METHOD.md` | §3 (lines 254-306) makes the CPCM "solvation treatment mismatch" a blanket caveat on ALL electrostatic results incl. the within axis; line 305 mandates "labelled as mismatched-treatment correlations." STATE's OPEN guard (lines 7-22) says this is the wrong, too-global excuse for the weak WITHIN field σ-response. Method itself stands; the blanket framing is partially superseded. | PARTIALLY SUPERSEDED |
| `analysis/FINDINGS.md` | (a) "in-plane radial R²=0.67 / point-dipole breaks down in the ring plane" (lines 36-37, 92-93) — J-B recovery (`165ed08`, `JOHNSON_BOVEY_REGION_RECOVERY.md`) shows that was producer-kernel form-reconstruction, NOT DFT signal (in-plane DFT ring-modulation ~0). (b) "literature-coefficient-FIXED test … TODO" (lines 152-155, 164) is the bare-kernel γ=−11.3 units bug. Living doc, already self-dated 05-31. | PARTIALLY SUPERSEDED |

## NOT bannered (checked, found current / correctly-handled)

- `analysis/EFG_ARC_EVIDENCE.md` — the EFG "O 0.34" is correctly presented as the **lab-frame
  rotation confound** (lab γ R²=0.315/0.342 → local ≈0; lines 37-44, 84-91), exactly the
  corrected reading. Local-frame de-circularisation reports honest null. Solvation caveat at
  112-113 is a between/absolute footnote (defensible), not a within-axis excuse. Current.
- `analysis/JOHNSON_BOVEY_REGION_RECOVERY.md` — this IS the corrective doc for the in-plane
  claim. Uses real literature-scaled `jb_T0`/`jb_T2_local_*` cols (audit confirms `jb_T*` are
  literature-ppm, `UNITS_AND_ISSUES_AUDIT.md:57`), oracle-PASS vs H5 `bare_T0`. Current.
- `analysis/PATTERNS.md` — guidance, not results. Lines 51-54 still call the per-cell emit a
  "fixed-literature kernel T2 for the de-circularising test" / "literature-kernel T2 (#34)";
  the audit shows that emit is actually bare-kernel-scale. Not a false RESULT, but the
  next-session brief should note the #34 emit must be literature-scaled (or relabeled) before
  any cell's step-5 de-circularisation is read as un-fitted. Flag only; no banner.
- `STATE.md`, `NEXT_SESSION_PROMPT.md` — canonical corrected docs. STATE top (OPEN guards +
  units RESOLVED, lines 7-54) supersedes its own older 184a1ee section (lines 124-125, which
  still states "UN-FITTED literature-coefficient kernel T2 predicts DFT for N 0.69, O 0.53,
  C 0.51"). NEXT_SESSION line 28-29 + 64 already carry "form-recovered, scale-fitted" + the
  open γ-vs-units thread. Left as-is — these are the docs that DO the correcting.

## Stale roadmap-doc rows the next session should not quote at face value

These are roadmap/status docs, not results docs; their stale rows are corrected by STATE top.
Not bannered (would clutter the roadmaps), but listed so the next session does not re-quote them:

- `REDISCOVERY_MAP.md:31` — broad-backbone "De-circularised (un-fitted literature kernel T2 vs
  DFT, component r): N 0.69, O 0.53, C 0.51 — textbook physics predicts … un-fitted." The
  kernel was bare, not literature; this is bare-kernel-scale, not an un-fitted literature
  prediction. The component-r magnitudes themselves are scale-invariant correlations and survive;
  the "textbook physics predicts un-fitted" gloss does not.
- `WORK_CATALOG.md:18` — same "un-fitted literature kernel predicts N 0.69 / O 0.53 / C 0.51"
  in the Done-arc summary. Same correction.
- `STATE.md:124-125` (older 184a1ee section) — same claim; superseded by STATE top item 1.

## Discipline spot-checks (all clean)

- Per-stratum: every results table breaks out N/CA/C/O/HN/HA(/HA2/HA3); no pooled summary
  stands in for a stratum. Good.
- Within-protein confidence only: docs consistently scope to "this one protein, ~500 frames,
  ~50 atoms/stratum," jackknife/block-bootstrap, AR(1) N_eff; no across-proteins population
  claim. Thin strata (HA2/HA3, 4 atoms) flagged. Good.
- No "treatment mismatch" hand-wave for the within axis after the VARIANCE banner — STATE's
  three falsifiable candidates (units/prefactor, projection axis, genuinely-weak-driver,
  decided by raw ΔField-vs-Δσ) are the concrete framing. Good.
- model-is-spine / correlate-not-match: the `|T2| r` + component-r headline (scale-invariant)
  is the right metric and is unaffected by the units bug; only the γ/"literature" framing is.

## Verdict for the next session — is the record honest + consistent?

Honest where it counts, with one labeling fault now bannered. The corrected canon
(`UNITS_AND_ISSUES_AUDIT.md`, `STATE.md` top, `JOHNSON_BOVEY_REGION_RECOVERY.md`,
`EFG_ARC_EVIDENCE.md`) is accurate: the de-circularisation "failure" is a units/labeling bug
not a physics failure; the EFG O signal is a tumbling confound; the J-B in-plane completion is
a true negative. The fault was that four results docs + three roadmap rows propagated
"fixed-literature ppm" / "calibrated-to-physics" language for what were bare-kernel-scale
correlations — now flagged. After the banners, a fresh session reading top-down (STATE →
UNITS_AUDIT → the bannered docs) will not be misled. The two load-bearing scale-invariant
results — ring k≈21 LOAO R²=0.62 (clean recovered law) and the per-stratum `|T2| r` /
component-r correlations — are untouched by the units bug and stand. The single most important
next-session action implied here: re-run the de-circularisation in literature-scaled units
(`jb_T*` / BS×LiteratureIntensity for ring; named Δχ+prefactor for McConnell) before the γ≈1
de-circularised verdict can be read for any mechanism. Scientific verdict reserved for the lead.

No merge. No commit beyond this doc + the 5 banners.
