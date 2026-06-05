# Data Dictionary v2 — Language Vet (2026-06-05)

Read-only language vet of `notes/DATA_DICTIONARY_v2_2026-06-05.md` (204 metrics) against the two
naming failure modes (unintelligible / too-subjective) and the description writing standard
(plain, grounded, not condescending, not a jargon poem). Spec: `notes/DICTIONARY_NAMING_PIPELINE_2026-06-05.md`.

This vet touches LANGUAGE only. Where a name is the producer's own field key and a chemist
cannot parse it (the Larsen 1p/2p terms, KTB), the fix is to spell the chemistry out in the
NAME so the reader is not sent to the code — exactly the "go further afield when the catalog
name is junk" instruction.

---

## VET SUMMARY

Overall: **strong pass.** The author's register is consistent — plain, declarative, "what it is /
what it is not," with units and tensor-rank discipline carried correctly throughout. The
description standard is met on the large majority of rows; the jargon-poem and condescension
failure modes are essentially absent (the author erred toward dry-but-correct, never toward
flowery). The flags are concentrated in two pockets: (1) a small family of producer-key NAMES
that no chemist can parse without the code (Larsen `1pHB`/`2pHaB`, `KTB`), and (2) a handful of
descriptions that restate the name or lean on an undefined in-house term.

Counts:

- **Names flagged — unintelligible: 8.** (6 Larsen pair-term names × the "1pHB/2pHaB"-style code
  key; 1 "KTB" in the spectral-density name; 1 "iRED" left unexpanded in the name.) Note: many of
  these recur as both a trajectory and a snapshot row, so the *distinct naming patterns* are
  fewer (~3) but each surfaced row is counted.
- **Names flagged — too-subjective: 0.** No name asserts an interpretation over identifying the
  quantity. The author was disciplined here; the closest calls ("contact tensor," "latent
  descriptors") are defensible identifiers, not editorialising.
- **Descriptions flagged — condescending: 0.**
- **Descriptions flagged — jargon-poem: 1.** (One "loses memory over lag" anthropomorphism; mild.)
- **Descriptions flagged — vague / restates-name / undefined-term: 6.**

Total flagged rows: **15 of 204 → ~93% clean pass-rate.**

Honest assessment: this is a publishable-quality first vet target. The author did the hard part
— meanings are grounded, units corrected against the generated catalog, and the "kernel not
shielding" distinction (A^-n vs ppm) is held everywhere it matters. The residue is naming
hygiene on a few inherited code-keys and six prose tightenings, none of which require re-deriving
physics. The 10 still-uncertain meanings (carried in the author's own "Still Uncertain" table)
are correctly fenced and not papered over — that honesty is itself a pass.

---

## CORRECTIONS TABLE

Only FLAGGED rows appear. Corrections stay in the author's grounded register; `id` is the row's
catalog id. For the NAME fixes, the corrected NAME is given; for DESCRIPTION fixes, the corrected
sentence(s).

### Names

| id | Failure mode | Corrected wording |
|---|---|---|
| `h5:larsen_hbond_1pHB_shielding_time_series` | Unintelligible — `1pHB` is a code key; a chemist cannot parse "1p"/"HB" without the source. | NAME → **Larsen primary amide-H backbone H-bond shielding trajectory** (keep `1pHB` only in the Origin column). |
| `h5:larsen_hbond_1pHaB_shielding_time_series` | Unintelligible — `1pHaB`. | NAME → **Larsen primary alpha-H backbone H-bond shielding trajectory**. |
| `h5:larsen_hbond_2pHB_shielding_time_series` | Unintelligible — `2pHB`. | NAME → **Larsen secondary amide-H backbone H-bond shielding trajectory**. |
| `h5:larsen_hbond_2pHaB_shielding_time_series` | Unintelligible — `2pHaB`. | NAME → **Larsen secondary alpha-H backbone H-bond shielding trajectory**. |
| `npy:larsen_hbond_1pHB_shielding` (and `2pHB`, `1pHaB`, `2pHaB` snapshots) | Unintelligible — same `Np`/`HaB` code keys in the snapshot rows. | NAME → mirror the trajectory fix: **Larsen {primary,secondary} {amide-H,alpha-H} backbone H-bond shielding snapshot**. |
| `h5:reorient_spectral_density` | Unintelligible — "KTB" is an in-field acronym (the 5-frequency relaxation sampling 0, ωN, ωH−ωN, ωH, ωH+ωN); a chemist will not know "KTB." | NAME → **Spectral density J(ω) at the five relaxation frequencies** (keep "KTB" out of the name; it can live in the description as the convention label if wanted). |
| `h5:ired_s2` | Unintelligible (borderline) — "iRED" unexpanded in the name. | NAME → **iRED (isotropic reorientational eigenmode) N-H order parameter** — expand on first/only use so the acronym is anchored. (Low severity; iRED is standard in relaxation work — acceptable to leave if the description expands it, which it does not currently.) |

Notes on the Larsen block: "primary/secondary" = the 1-bond / 2-bond pair terms, and "HB" vs
"HaB" = the amide-H (backbone N-H) vs alpha-H (CαH) acceptor families, per
`doc/calculators/larsenhbond-numerics.tex:357-365`. The producer key stays authoritative in
Origin; only the human-facing NAME changes. This is the single highest-value naming fix in the
document — six rows that are currently code-only.

### Descriptions

| id | Failure mode | Corrected wording |
|---|---|---|
| `h5:kernel_dynamics_acf` | Jargon-poem (mild anthropomorphism) — "how quickly each channel loses memory over lag" is cute but imprecise. | "The curves measure how fast each channel decorrelates with increasing time lag." |
| `npy:residue_index` | Vague / restates name — "It maps atom-level descriptors back to residue-level context" is near-circular. | "Integer index of the residue each atom belongs to; used to group atom-level descriptors by residue." |
| `npy:residue_type` | Vague / restates name — "supplies the residue chemistry context used by selectors and descriptors" says little. | "Residue (or residue-variant) name for each atom — e.g. ALA, HID — giving the residue chemistry used by selectors and per-residue descriptors." |
| `topology:residues` | Restates name — "the residue-level identity map used to group atoms" echoes the name without adding content. | "Per-residue identity table from the topology sidecar: residue name, number, and chain, keyed to the atoms it contains." |
| `geometry:atom_displacement` | Vague (acknowledged) — "from the active reference used by the dataset or view" is genuinely unclear, and the prose admits it. | Keep, but tighten the hedge: "Per-atom displacement vector relative to a reference configuration. The reference-frame/alignment convention is not pinned down by the catalog — treat the direction as dataset-defined until confirmed." (See still-needs-a-human.) |
| `npy:aimnet2_aim` | Vague — "atom-wise learned descriptors or AIM-related quantities" hedges two different meanings without landing either. | "AIMNet2 per-atom outputs beyond charge — learned atomic descriptors and/or atoms-in-molecules-style quantities. The exact column packing of this block is not grounded in reader-side docs." (Honest, and flags the gap inline.) |

---

## STILL-NEEDS-A-HUMAN

Items I cannot confidently fix on language grounds alone, because the residual is a *meaning* gap,
not a wording gap. The lead (or a code-grounded pass) should resolve these; a vetter inventing
prose here would be guessing.

1. **The author's own 10 "Still Uncertain" meanings** (their §"Still Uncertain"). I endorse all 10
   as genuinely unresolved from the allowed sources and have NOT papered over them:
   - `geometry:atom_displacement` — reference-frame/alignment convention unknown.
   - `npy:ring_contributions` — 58 columns confirmed, per-column packing not.
   - `npy:ring_geometry` — 10 columns, exact packing not spelled out.
   - `npy:larsen_hbond_diagnostic_CB_shielding` — "diagnostic CB" operational meaning not found.
   - `npy:delta_scalars` — only total-delta-T0 column grounded.
   - `npy:delta_apbs` — 12-column packing + sign convention not grounded.
   - `npy:delta_ring_proximity` — variable-column packing + sign/proximity convention not grounded.
   - `h5:kernel_dynamics_psd` — PSD normalization/window convention unspecified.
   - `npy:aimnet2_aim` — exact column packing of the non-charge block not grounded.
   For each of these the NAME is fine and the DESCRIPTION's hedge is honest; they need a fact, not
   a rewrite.

2. **"KTB" provenance for the spectral-density description.** I moved the acronym out of the NAME
   (the fix above), but if the team wants to *name* the convention in the description, someone
   should confirm what KTB stands for in this codebase before it is written into a human-facing
   doc. (The reader source uses it as a bare label; I did not chase it to a citation.)

3. **`larsen_hbond_diagnostic_CB_shielding` NAME.** I left the name as-authored ("Larsen diagnostic
   CB shielding") because "CB" = carbon-beta is parseable, but whether "diagnostic" is the right
   human label depends on the unresolved meaning in item 1 — flagging so it is not lost.

4. **GROMACS energy snapshot column count (42 vs 43).** The description honestly reports the
   catalog/reader disagreement; that is a producer-side fact to reconcile, not a language fix.
   Noted here only so the contradiction is not mistaken for a prose defect.

---

## Method note

Vet was language-only and read-only. To ground the cryptic-name fixes I consulted
`doc/calculators/larsenhbond-numerics.tex` (primary/secondary + amide-H/alpha-H term naming) and
the h5-reader source comments for the KTB frequency set (`FixedFreqPanel.cpp:165`,
`QtBondVectorBuffers.h:108`). No `nmr_shielding` library files were read; nothing outside
`h5-reader/` was touched; no code or git was modified.
