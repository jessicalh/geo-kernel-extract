# Codex — data dictionary: AUTHOR pass (names + descriptions, sourced beyond the SDK)

## (Preamble prepended. READ-ONLY authoring — NO code, NO git. Cite sources. You MAY read `../doc/calculators/` READ-ONLY for physics meaning; do NOT read the `nmr_shielding` library source (repo-root `src/`) and modify NOTHING outside `h5-reader/`.)

## Why / what
This is the AUTHOR pass of the naming pipeline (spec: `notes/DICTIONARY_NAMING_PIPELINE_2026-06-05.md`).
A first draft exists: `notes/DATA_DICTIONARY_DRAFT_2026-06-05.md` (204 metrics, 84 `[VERIFY]`). Refine it
into a production names+descriptions set. A SEPARATE opus agent will VET your output afterward (for
unintelligible / too-subjective names and for plain-grounded-not-a-jargon-poem descriptions) — so your
job is the best-sourced, honest draft, not the final word.

## Do — for every metric in the draft
1. **NAME** — a meaningful human name. The catalog/SDK label is the FIRST source, but where it is
   junk/cryptic/abbreviation-soup, go FURTHER AFIELD: the `../doc/calculators/` exposition docs, standard
   NMR/MD terminology. LOCATE-BEFORE-ABSENT — exhaust the sources before settling for a poor name.
   Resolve the 84 `[VERIFY]` tags wherever the sources let you.
2. **DESCRIPTION** — 1–3 sentences on what it ACTUALLY means, plain and grounded: not condescending,
   not a jargon poem. What a chemist can act on. (Draft to your best; the vetter polishes the language.)
3. **UNITS / rank-shape** — confirm or fix. T2 stays T2 — never collapse a tensor to a scalar.
4. Keep the **origin/infra** column (descriptor id, concept key, storage path, source, catalog `file:line`).
5. **EXPECTED-BUT-EMPTY** — if a metric is catalogued but you have grounds to think its data is typically
   absent/empty (vs genuinely populated), FLAG it for the reality-check list (this feeds the registry's
   empty-check). Do NOT silently drop it.

## Output
Write `notes/DATA_DICTIONARY_v2_2026-06-05.md`: the refined table (Name | What it is | Units | Rank/shape |
Origin), grouped by physics category; plus a short "still-uncertain" list (what you could not ground) and an
"expected-but-empty candidates" list. State the metric count and how many `[VERIFY]` you resolved vs remain.
No code, no git. Note anything in the catalog that contradicts the draft.
