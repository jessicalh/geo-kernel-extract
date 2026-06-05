# Dictionary / metric-naming production pipeline — spec (2026-06-05)

The first draft (`notes/DATA_DICTIONARY_DRAFT_2026-06-05.md`, 204 metrics, 84 `[VERIFY]`) was a
single SDK-bound codex pass. The lead's three requirements turn it into a rigorous, multi-agent
production pipeline. This doc is the spec; agents run it on her go.

## Why (lead, 2026-06-05)
Metric selection must lead with a **meaningful name + plain description**, infra/origin tucked
out of sight. The draft was too catalog-bound, didn't vet its own names, and treated every
empty thing as cruft. Three fixes:

## 1. NAME — author, then vet
- **Author.** A meaningful, human name for each metric. Source from the catalog/SDK FIRST, but
  **go further afield when the catalog name is junk** — the per-calculator exposition docs
  (`../doc/calculators/`, read-only), standard NMR/MD literature, domain knowledge. Honour
  LOCATE-BEFORE-ABSENT: exhaust the sources before settling for a poor name.
- **Vet.** A separate agent checks every name against two failure modes:
  - **Unintelligible** — cryptic, abbreviation-soup, no one outside the code would parse it.
  - **Too subjective** — arbitrary/opinionated/editorialising; a name that asserts an
    interpretation rather than identifies the quantity.
  Flag + propose a fix for each; the name must survive both.

## 2. EMPTY = REALITY CHECK, not a silent cut
The hollow/inactive will NOT all just "fall out." Distinguish:
- **No renderer / no data ever** → structural cut (the registry handles this; silent is fine).
- **We expected data here and it came back EMPTY** → a genuine reality check: producer gap, or
  are we reading it wrong? This list does NOT vanish silently — it goes to the **log / a
  report** ("expected something for X, found empty"). LOCATE-BEFORE-ABSENT first: check the
  canonical name AND plausible alternates before declaring empty. This is a fine reality-check
  moment — a data-integrity signal, not UI cruft. (Implementation home: the Stage-1 registry
  offerability gate — see `notes/STAGE1_REGISTRY_PLAN_2026-06-05.md`. The dictionary pipeline
  contributes the "what we EXPECTED to be populated" list that the gate checks against.)

## 3. DESCRIPTION — author, then language-vet
- **Author.** 1–3 sentences on what the metric ACTUALLY means, in **plain, grounded language**:
  not condescending, not a jargon poem. Say the literal truth a chemist can act on.
- **Language vet.** A separate agent holds the prose to exactly that bar (the project technical-
  writing standard: exquisite knowledge + taste + intellectual empathy; plain ≠ dry). Rewrites
  the condescending and the jargon-poem alike.

## Shape
Relay/pipeline per metric: name(author → vet) ∥ description(author → vet); the expected-but-empty
set is emitted to the reality-check log. Output: the production dictionary + the reality-check
report. Grounded in: technical-writing-standard, plain-spoken-but-charming,
locate-before-absent, darlings-vs-cruft. Run as its own pass; independent of the build stages,
can go in parallel once the lead green-lights the agent spend.
