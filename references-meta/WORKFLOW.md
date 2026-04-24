# references-meta — Companion files for the PDF corpus

Every PDF in `../references/` gets two companion files here:

- `<basename>-summary.txt` — plain-English summary, 100–500 words body.
- `<basename>-keywords.txt` — one keyword per line, content terms only.

`<basename>` matches the PDF filename so the correspondence is obvious at the shell. When the corpus reaches Zotero, these files become the `Abstract` (summary text) and `Tags` (keywords, comma-joined at upload) fields.

## When to produce

At PDF acquisition time, before the paper joins the citable corpus. Reading without summarising lets comprehension fall behind acquisition and rots the corpus by the time it matters.

## Ingestion pipeline

When a PDF lands in `references/incoming/`, run the ingest tool to produce the text and image layers that the summary workflow depends on:

```bash
scripts/references/ingest_pdf.sh incoming/<messy-pdf-name>.pdf <basename>
```

Produces:
- `references/<basename>.pdf` — PDF moved to final location.
- `references-text/<basename>-text-N.txt` — 3-page text chunks (default).
- `references-images/<basename>-page-N.png` — per-page renders at 150 DPI.

Then produce the summary + keywords files per the discipline below. For dense papers, the chunked text enables the scratch-then-compose pattern (see "Chunk-and-distill for dense papers" below) which avoids hitting Claude's working-memory ceiling on 15+ page reviews.

Details: `scripts/references/README.md`.

## Summary format

Line 1 (header metadata): full citation. Include DOI, open-access status (PMCID if present), page count, reference count. If the paper declares competing interests or funding that shapes its argument, add one short sentence noting this in the metadata block.

Body (100–500 words, no section headers, paragraphs flow):

1. **What the paper is** — scope and stance (research, Perspective, review, textbook chapter); what it argues or measures.
2. **Methods / central results** — for a research paper: methods, dataset, scope. For a review or Perspective: the central theorems, arguments, or organising principles.
3. **Findings / intellectual narrative** — for a research paper: what was observed, what was not observed, magnitudes. For a review: who did what, how the field arrived here.
4. **Relations to our corpus** — which papers referenced here we already have, which the paper points at that we might want. Secondary-cite flagging is useful for reviews.
5. **Implications for our work, widely considered** — across thesis stages, methods chapter, committee questions. This is the most project-specific paragraph; it should not make factual claims beyond what the paper supports.

## Voice

- Plain English first; terms of art defined briefly inline the first time they appear. No reaching for jargon to sound technical.
- Dry. No "pioneering", "crucial", "groundbreaking", "elegant", "profound", "remarkable", "novel", "fundamental", "comprehensive" unless quoting the authors. If the paper uses them, direct quotes are allowed.
- Do not editorialise in the summary voice. "Our kit differs from theirs" framing belongs in the lit-review outline or in eventual evidence cards, not here.
- Cite page numbers when quoting; quotes are verbatim.

## Two-draft discipline

1. **Draft 1 — associative pass.** Fast, vocabulary-rich, short. Makes sure nothing is missed. Internal only.
2. **Draft 2 — ultrathink pass.** Structured, written for the reader, paragraphs in plain English with terms of art defined. Only Draft 2 is saved.

Show both drafts in chat only if the user asks, or on the first summary of a new paper type (research vs. review vs. textbook) as a format check. Otherwise Draft 1 stays internal.

## Reading the chunks

The text chunks in `references-text/<basename>-text-N.txt` are the primary read surface. The PDF itself is for Jessica's reading; the per-page images are on-demand when text references a figure that matters. Every paper gets the same uniform chunk layout regardless of length — the ingest pipeline does not distinguish "dense" from "not dense." That call belongs to the reader.

After reading chunk 1, pick one of two composition modes:

**One-shot mode** (default for short or loosely-structured papers): read the remaining chunks in sequence, hold the content in working context, go straight to the two-draft summary flow. Skip scratch files.

**Scratch-and-compose mode** (default for long / visually-dense / figure-heavy papers, and whenever one-shot would risk hitting the working-memory ceiling):

1. Read chunk 1 → write `references-meta/<basename>-scratch-1.txt` (200–500 words of structured notes: what this chunk covers, key numbers, cited references, relevance hooks). Stop.
2. Repeat for chunks 2 through N. One scratch file per chunk.
3. Compose Draft 2 by reading **only the scratch files + the paper's abstract/first page** — not the full text chunks. Working memory during the compose step is ~500 words per chunk × N chunks instead of full-paper size.
4. Save Draft 2 as `<basename>-summary.txt`. Delete the scratch files once the summary is verified (they are intentionally ephemeral).

The decision is yours after chunk 1 — ingest does not make it. If you start in one-shot mode and find chunk 3 or 4 is pushing working memory hard, you can switch to scratch mode retroactively by writing scratches for chunks already read.

## Snarf at writing time

When you need to pull everything about a topic from the corpus at writing time, spawn a subagent rather than trying to fit the corpus into main context:

1. Main thread spawns an Explore or general-purpose subagent with the topic query.
2. Subagent greps `references-meta/*-keywords.txt` first (high-signal, small payload), then `references-meta/*-summary.txt` (lower-signal, larger payload) for the topic term and close synonyms.
3. Subagent returns a ranked list of papers with one-line why-relevant hooks.
4. Main thread picks 2–5 papers; if deeper extraction is needed, spawn a second subagent to read specific chunks from `references-text/`.
5. Main thread synthesises and writes. Full-paper payloads never enter main context.

For quick scans without a subagent, `scripts/references/regen_index.sh` regenerates `references-meta/INDEX.md`, a keyword → paper cross-reference you can grep directly.

## Citation discipline

Nothing in a summary is a citation to a paper we do not own the PDF for. If a reviewed paper mentions a third paper we do not have, flag it by author-year as a **pointer for acquisition**, not a citation. Downstream argument-writing cannot use pointers as citations — only cards derived from PDFs.

## Keywords format

One keyword per line. Line 1: `# <cite-key>` comment so the file can find its paper if pulled out of context. Keep keywords to paper-content only (methods, models, concepts). Thesis-structural tags (`stage-1`, `trajectory-averaging`, `geometric-kernel`) are a separate axis, not yet settled on format; do not mix into the paper keyword file.

## Supporting Information (SI) files

When a paper's supporting information is available as a separate PDF (conventionally `<basename>-SI.pdf`), ingest it like any other paper but do **not** produce an independent summary or keywords file. The SI content should be covered by the parent paper's summary, typically with a one-sentence note ("Supplementary companion `<basename>-SI.pdf` contains additional proofs / data; covered by this summary per the SI convention"). Rationale: SIs exist to support a paper already in the corpus; splitting them into independent summaries fragments the evidence trail and inflates the INDEX noise-floor. If the SI material matters enough to cite on its own, promote it to its own full summary; that's a conscious choice.

## Length calibration

- Short research paper, 4–6 pages: ~200–300 words.
- Standard research paper, 6–12 pages: ~300–450 words.
- Dense review or Perspective, 10+ pages, many refs: ~450–500 words.
- Textbook chapters and books: may need their own structure; not covered here.

## Paper-type templates

Research paper — Paragraph 2 is methods-heavy; Paragraph 3 is findings/magnitudes; Paragraph 4 shorter (research papers usually reference less widely than reviews).

Review / Perspective — Paragraph 2 is central arguments or organising framework; Paragraph 3 is the intellectual narrative and landscape; Paragraph 4 is the most useful paragraph because reviews exist to map relationships.

Textbook excerpt — rarely worth a summary file; excerpt and store the relevant chapters as PDFs with short abstracts.

## Tentative observations → TICKLERS.md

When a summary surfaces an observation that *might* affect our code, statistics, or argument but is not ready to act on (foreign framing, narrow calibration context, unreplicated claim, etc.), do not bury it in the summary. Add an entry to `TICKLERS.md` in the format documented there. Ticklers are explicitly non-load-bearing — they exist so future coding and stats sessions can review them with full context and decide whether to promote (to memory or evidence card) or reject (with a journal note).

## Companion documents

- `WORKFLOW.md` — this file, rules for producing summary + keywords files.
- `TICKLERS.md` — tentative observations from summaries that may bear on code or argument.
- `INDEX.md` — auto-generated keyword → paper cross-reference, regenerated via `scripts/references/regen_index.sh`.
- `../references/PENDING_ACQUISITIONS.md` — papers we want but do not yet have on disk.
- `../scripts/references/README.md` — tool documentation for the ingest and index scripts.
