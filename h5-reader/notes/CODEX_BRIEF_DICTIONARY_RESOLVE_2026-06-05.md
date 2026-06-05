# Codex — resolve the 18 unresolved dictionary items from the PRODUCER source + web

## Scope / license — READ CAREFULLY (this is NOT a standard h5-reader brief)
The reader is downstream and cannot ground these meanings. The PRODUCER — the `nmr_shielding`
library calculators that COMPUTE and EMIT these NPY/H5 fields — is the authority. For THIS task you
have EXPLICIT license to read, READ-ONLY:
- The producer source: `/shared/2026Thesis/nmr-shielding/src/` (the calculators that compute + emit
  each field), `/shared/2026Thesis/nmr-shielding/spec/`, and the SDK catalog
  `/shared/2026Thesis/nmr-shielding/python/nmr_extract/_catalog.py`.
- `/shared/2026Thesis/nmr-shielding/doc/calculators/` (the per-calculator exposition docs).
- The open WEB — cite URLs — for standard NMR/MD terms (KTB spectral density, iRED / Lipari–Szabo,
  Larsen H-bond conventions, etc.).
MODIFY NOTHING anywhere (read-only); no git; no build. **Verify from the producer source before
asserting; cite `file:line` for every resolution.** "I couldn't ground it" beats a guess.

## The 18 to resolve (= 8 cryptic names + 10 fact-gap meanings)
Read `notes/DATA_DICTIONARY_VET_2026-06-05.md` (the vet's flagged + still-needs-a-human lists) and
`notes/DATA_DICTIONARY_v2_2026-06-05.md` (the dictionary). The 7 wording-only flags (1 jargon-poem +
6 vague) are ALREADY fixed by the vet — IGNORE those. Resolve the 18 that need real grounding:
- **8 cryptic NAMES** — the six Larsen pair-terms (`1pHB` / `2pHaB` / `1pHaB` / `2pHB`), `KTB` in the
  spectral-density name, and `iRED`. Find what each ACTUALLY is from the producer calculator source
  (the Larsen H-bond calculator, the kernel-dynamics / spectral-density code, the iRED code) + web,
  and give the chemistry-spelled-out human name.
- **10 still-uncertain MEANINGS** — the fact-gaps the author/vet flagged (incl. KTB provenance and
  the diagnostic-CB term if not already covered). Find the true meaning + units + shape from the
  producer source that emits each.
List EXACTLY which items you took (the vet's two lists define them; the lead's count is ~18).

## For each item produce
field id | RESOLVED human name (chemistry spelled out) | grounded 1–3 sentence meaning (plain, not a
jargon poem) | units | shape — each citing the PRODUCER source `file:line` (+ web URL for a standard
term). If a field STILL can't be resolved even from the producer, say so honestly and state what's
missing — do not guess.

## Output
Write `notes/DATA_DICTIONARY_RESOLVED_2026-06-05.md`: lead with a count (resolved vs still-unresolved),
then the table of the ~18. This feeds the final dictionary merge. No code, no git.
