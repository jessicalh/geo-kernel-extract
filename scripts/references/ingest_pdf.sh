#!/usr/bin/env bash
#
# ingest_pdf.sh — ingest a PDF into the reference pipeline.
#
# Usage:
#   ingest_pdf.sh <pdf_path> <basename>
#
# Produces (all using <basename> as stem):
#   references/<basename>.pdf                       — the PDF (moved from input)
#   references-text/<basename>-text-N.txt           — 3-page text chunks
#   references-images/<basename>-page-N.png         — per-page images at 150 DPI
#
# Conventions:
#   - <basename> follows our naming scheme: firstauthor(-coauthor)-year-short-title
#   - The input PDF is MOVED (not copied) to references/. Save your original elsewhere
#     if you want to preserve it.
#   - Idempotent: output files that already exist are not overwritten. Re-running is
#     safe and will only fill in missing outputs.
#
# Environment overrides:
#   CHUNK_PAGES   (default 3)     — pages per text chunk
#   IMAGE_DPI     (default 150)   — resolution for per-page renders
#   SKIP_IMAGES   (default unset) — set to any value to skip image extraction entirely
#
# After running, produce companion files in references-meta/ per the discipline in
# references-meta/WORKFLOW.md:
#   - <basename>-summary.txt  (200–500 words, two-draft process)
#   - <basename>-keywords.txt (one content term per line)
#
# Requires: pdftotext, pdftoppm, pdfinfo (all from poppler-utils).

set -euo pipefail

# ─── args ──────────────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <pdf_path> <basename>" >&2
    echo "" >&2
    echo "Example: $0 incoming/messy-elsevier-name.pdf klukowski-2025-ml-nmr" >&2
    exit 2
fi

PDF_IN="$1"
BASENAME="$2"
CHUNK_PAGES="${CHUNK_PAGES:-3}"
IMAGE_DPI="${IMAGE_DPI:-150}"

# ─── resolve paths relative to project root ────────────────────────────────
# Script lives at scripts/references/ingest_pdf.sh; project root is two levels up.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

REF_DIR="$PROJECT_ROOT/references"
TEXT_DIR="$PROJECT_ROOT/references-text"
IMG_DIR="$PROJECT_ROOT/references-images"
META_DIR="$PROJECT_ROOT/references-meta"

mkdir -p "$TEXT_DIR" "$IMG_DIR" "$META_DIR"

# ─── sanity checks ─────────────────────────────────────────────────────────
if [[ ! -f "$PDF_IN" ]]; then
    echo "error: input PDF not found: $PDF_IN" >&2
    exit 1
fi

for cmd in pdftotext pdftoppm pdfinfo; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "error: required command not found: $cmd (install poppler-utils)" >&2
        exit 1
    fi
done

# ─── page count ────────────────────────────────────────────────────────────
NPAGES=$(pdfinfo "$PDF_IN" | awk '/^Pages:/ {print $2}')
if [[ -z "$NPAGES" || "$NPAGES" -lt 1 ]]; then
    echo "error: could not determine page count for $PDF_IN" >&2
    exit 1
fi

echo "ingesting: $PDF_IN"
echo "  basename:    $BASENAME"
echo "  pages:       $NPAGES"
echo "  chunk size:  $CHUNK_PAGES pages"
echo "  image dpi:   $IMAGE_DPI"
echo ""

# ─── text chunks ───────────────────────────────────────────────────────────
NCHUNKS=$(( (NPAGES + CHUNK_PAGES - 1) / CHUNK_PAGES ))  # ceil
echo "extracting $NCHUNKS text chunk(s)..."
for (( i=1; i<=NCHUNKS; i++ )); do
    START=$(( (i - 1) * CHUNK_PAGES + 1 ))
    END=$(( i * CHUNK_PAGES ))
    if [[ $END -gt $NPAGES ]]; then END=$NPAGES; fi

    OUT_TXT="$TEXT_DIR/${BASENAME}-text-${i}.txt"
    if [[ -f "$OUT_TXT" ]]; then
        echo "  [skip] $OUT_TXT (exists)"
    else
        pdftotext -f "$START" -l "$END" -layout "$PDF_IN" "$OUT_TXT"
        echo "  [ok]   $OUT_TXT (pages $START-$END)"
    fi
done

# ─── per-page images ───────────────────────────────────────────────────────
if [[ -z "${SKIP_IMAGES:-}" ]]; then
    echo ""
    echo "rendering $NPAGES page image(s) at ${IMAGE_DPI} DPI..."
    # pdftoppm writes files named <prefix>-<pagenum>.png with zero-padded numbers
    # depending on page count. We normalise by rendering one page at a time.
    for (( p=1; p<=NPAGES; p++ )); do
        OUT_PNG="$IMG_DIR/${BASENAME}-page-${p}.png"
        if [[ -f "$OUT_PNG" ]]; then
            echo "  [skip] $OUT_PNG (exists)"
        else
            # pdftoppm appends -N to the prefix; write to a temp prefix then rename
            TMP_PREFIX="$IMG_DIR/${BASENAME}-tmp"
            pdftoppm -f "$p" -l "$p" -r "$IMAGE_DPI" -png "$PDF_IN" "$TMP_PREFIX"
            # find the produced file (pdftoppm pads with zeros if npages requires)
            PRODUCED=$(ls "${TMP_PREFIX}"-*.png 2>/dev/null | head -1)
            if [[ -n "$PRODUCED" && -f "$PRODUCED" ]]; then
                mv "$PRODUCED" "$OUT_PNG"
                echo "  [ok]   $OUT_PNG"
            else
                echo "  [warn] page $p did not produce an image"
            fi
        fi
    done
else
    echo ""
    echo "skipping image extraction (SKIP_IMAGES set)"
fi

# ─── move PDF to references/ ───────────────────────────────────────────────
OUT_PDF="$REF_DIR/${BASENAME}.pdf"
if [[ "$(readlink -f "$PDF_IN")" == "$(readlink -f "$OUT_PDF")" ]]; then
    echo ""
    echo "PDF is already at final location: $OUT_PDF"
elif [[ -f "$OUT_PDF" ]]; then
    echo ""
    echo "warn: $OUT_PDF already exists; input $PDF_IN was NOT moved."
    echo "      review manually to avoid losing the newer copy."
else
    mv "$PDF_IN" "$OUT_PDF"
    echo ""
    echo "moved PDF: $OUT_PDF"
fi

# ─── next steps ────────────────────────────────────────────────────────────
echo ""
echo "────────────────────────────────────────────────────────────────"
echo "ingest complete for $BASENAME."
echo ""
echo "next: produce companion files per references-meta/WORKFLOW.md"
echo "  - references-meta/${BASENAME}-summary.txt"
echo "  - references-meta/${BASENAME}-keywords.txt"
echo ""
echo "scratch notes (one per chunk) if the paper is dense:"
echo "  - references-meta/${BASENAME}-scratch-N.txt  (delete after summary lands)"
echo "────────────────────────────────────────────────────────────────"
