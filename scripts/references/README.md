# scripts/references — reference corpus tools

Vetted tools for the reference-corpus pipeline. The pipeline is described in full in `references-meta/WORKFLOW.md`; this README is the tool-side reference.

## Pipeline overview

Four layers, all sharing a common `<basename>` per paper (convention: `firstauthor(-coauthor)-year-short-title`):

```
references/<basename>.pdf                  — PDF for Jessica (human reading)
references-text/<basename>-text-N.txt      — 3-page text chunks (AI reading with bounded working memory)
references-images/<basename>-page-N.png    — per-page renders at 150 DPI (figure-on-demand)
references-meta/<basename>-summary.txt     — 200–500 word summary (two-draft workflow)
references-meta/<basename>-keywords.txt    — one content term per line (Zotero tags + snarf index)
```

Optional: `references-meta/<basename>-scratch-N.txt` — per-chunk scratch notes for dense papers, deletable once the summary lands.

## Tools

### `ingest_pdf.sh`

Run at PDF acquisition time. Produces text chunks, per-page images, and moves the PDF to its final location.

```bash
scripts/references/ingest_pdf.sh <pdf_path> <basename>
```

Example:
```bash
scripts/references/ingest_pdf.sh incoming/1-s2.0-S0079656525000196-main.pdf \
    klukowski-2025-machine-learning-nmr-spectroscopy
```

**Environment overrides:**
- `CHUNK_PAGES` (default `3`) — pages per text chunk. Denser papers benefit from smaller chunks to keep AI working-memory light.
- `IMAGE_DPI` (default `150`) — image resolution. `100` is smaller but still readable; `300` for figure-heavy papers where detail matters.
- `SKIP_IMAGES` (set to any value) — skip image extraction entirely. Useful when storage is a concern and figures aren't load-bearing.

**Idempotency:** re-running is safe. Outputs that already exist are skipped; missing outputs are filled in. Useful after a partial failure or when upgrading the chunk size.

**Dependencies:** `pdftotext`, `pdftoppm`, `pdfinfo` (all from `poppler-utils`; `apt install poppler-utils` on Ubuntu).

### `regen_index.sh`

Regenerate `references-meta/INDEX.md` — the keyword → paper-basename cross-reference.

```bash
scripts/references/regen_index.sh
```

Scans every `*-keywords.txt`, groups papers by keyword (case-folded), and writes a sorted Markdown index. Overwrites the existing INDEX.md. Safe to run any time.

**Caveat:** keywords are free-form, not a controlled vocabulary. Near-synonyms appear as separate entries (`ring current` vs `ring-current`). When snarfing, grep multiple phrasings.

## The snarf pattern

For topic queries at writing time — "let's get everything about __" — the recommended flow is:

1. Main thread spawns a subagent (Explore or general-purpose).
2. Subagent greps `references-meta/*-keywords.txt` and `references-meta/*-summary.txt` for the topic term(s).
3. Subagent returns a ranked list: papers where the topic is a keyword (high relevance) and papers where the topic is only implied in the summary (lower relevance).
4. Main thread picks 2–5 papers; spawn a second subagent to read the specific text chunks from `references-text/` if deeper extraction is needed.
5. Main thread writes prose using what comes back. Full-paper payloads never enter main context.

This keeps writing sessions efficient and preserves the evidence-discipline convention (nothing cited without the PDF on disk).

For quick manual snarfs without a subagent:
```bash
# all papers whose keywords mention "ring current"
grep -li "ring current" references-meta/*-keywords.txt

# summaries mentioning "order parameter"
grep -l "order parameter" references-meta/*-summary.txt
```

## Scale notes

Designed to hold up to 300–500 papers. At that scale:
- Text chunks: ~10–30 files per paper × 500 = ~15,000 files, ~500 MB.
- Per-page images at 150 DPI: ~10 pages × 300 KB × 500 = ~1.5 GB. Worth keeping, can be skipped via `SKIP_IMAGES` for batch re-ingests.
- Keyword files: ~500 × 40 terms ≈ 20,000 term rows, ~1 MB. INDEX.md will become ~5,000 unique keywords; browser-scrollable but probably wants chunking at that scale.

Gitignore policy (adjust as needed):
- `references/` and `references-meta/` are committed (the corpus and its distilled outputs are project-durable).
- `references-text/` and `references-images/` are gitignored candidates — trivially regenerable from the PDFs via `ingest_pdf.sh`.

## Future additions (not yet implemented)

- `snarf_topic.sh` — grep-based wrapper around the snarf pattern, for when spawning a subagent is overkill.
- Figure-selective extraction — pdftoppm renders entire pages; tools like `pdffigures2` could isolate just figures + captions, saving space.
- Bib-entry generator — produce a BibTeX line from a summary file's header, for pasting into LaTeX.
