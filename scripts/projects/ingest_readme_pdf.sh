#!/usr/bin/env bash
# ingest_readme_pdf.sh — capture a GitHub README PDF as a project entry.
#
# For when Jessica has print-to-PDF'd a GitHub README and dropped it into
# references/incoming/. This is the project-side counterpart to ingest_pdf.sh:
# no chunking, no image rendering, no Qwen pass — the README text IS the
# durable record. PDFs that are actual papers go through ingest_pdf.sh as
# usual; this is for README PDFs only.
#
# Usage:
#     scripts/projects/ingest_readme_pdf.sh <pdf_path> [<slug>]
#
# Produces:
#     projects/<slug>.md      — extracted README text + YAML frontmatter
#     projects/<slug>.pdf     — original PDF moved here (source preserved)
#     projects/INDEX.md       — appended index entry (file created on first run)
#
# If <slug> is omitted, it's derived from any github.com/<org>/<repo> URL found
# in the PDF text (preferred), or from the PDF filename (fallback). Final
# slug is always lowercased and hyphen-separated.
#
# No GitHub API calls. URL detection is best-effort against the PDF text only.
# Hand-edit projects/<slug>.md frontmatter or projects/INDEX.md to add or
# correct metadata.
#
# Idempotent: skips if projects/<slug>.md already exists. Use --force to overwrite.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <pdf_path> [<slug>] [--force]" >&2
    exit 2
fi

PDF_IN="$1"
SLUG_IN=""
FORCE=0
shift
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        *) [[ -z "$SLUG_IN" ]] && SLUG_IN="$arg" ;;
    esac
done

if [[ ! -f "$PDF_IN" ]]; then
    echo "error: input PDF not found: $PDF_IN" >&2
    exit 1
fi

for cmd in pdftotext jq; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "error: required command not found: $cmd" >&2
        exit 1
    fi
done

# ─── paths ───────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_DIR="$PROJECT_ROOT/projects"
INDEX="$OUT_DIR/INDEX.md"

mkdir -p "$OUT_DIR"

# ─── extract text ────────────────────────────────────────────────────────
TXT_TMP="$(mktemp)"
trap 'rm -f "$TXT_TMP"' EXIT
pdftotext -layout "$PDF_IN" "$TXT_TMP"

# ─── detect github URL ───────────────────────────────────────────────────
# Match github.com/<org>/<repo> with permissive chars; first hit wins.
GITHUB_URL=""
DETECTED_ORG=""
DETECTED_REPO=""
URL_LINE=$(grep -oE 'https?://github\.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+' "$TXT_TMP" | head -1 || true)
if [[ -z "$URL_LINE" ]]; then
    # Fall back to bare github.com/org/repo
    URL_LINE=$(grep -oE 'github\.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+' "$TXT_TMP" | head -1 || true)
    if [[ -n "$URL_LINE" ]]; then
        URL_LINE="https://$URL_LINE"
    fi
fi

if [[ -n "$URL_LINE" ]]; then
    # Strip any trailing punctuation, .git
    GITHUB_URL=$(echo "$URL_LINE" | sed -E 's#\.git$##; s#[.,;:]+$##')
    REPO_PATH=$(echo "$GITHUB_URL" | sed -E 's#^https?://github\.com/##')
    DETECTED_ORG=$(echo "$REPO_PATH" | cut -d/ -f1)
    DETECTED_REPO=$(echo "$REPO_PATH" | cut -d/ -f2)
fi

# ─── determine slug ──────────────────────────────────────────────────────
if [[ -n "$SLUG_IN" ]]; then
    SLUG=$(echo "$SLUG_IN" | tr '[:upper:]' '[:lower:]')
elif [[ -n "$DETECTED_ORG" && -n "$DETECTED_REPO" ]]; then
    SLUG=$(echo "${DETECTED_ORG}-${DETECTED_REPO}" | tr '[:upper:]' '[:lower:]')
else
    # Fallback: PDF filename
    BASE=$(basename "$PDF_IN" .pdf)
    SLUG=$(echo "$BASE" | tr '[:upper:]' '[:lower:]' | tr ' _' '--' | tr -cd 'a-z0-9-')
    echo "warn: no github URL detected in PDF; slug derived from filename: $SLUG" >&2
    echo "      hand-edit projects/$SLUG.md frontmatter to add the URL." >&2
fi

OUT_MD="$OUT_DIR/$SLUG.md"
OUT_PDF="$OUT_DIR/$SLUG.pdf"

# ─── idempotency ─────────────────────────────────────────────────────────
if [[ -f "$OUT_MD" && $FORCE -eq 0 ]]; then
    echo "skip: $OUT_MD exists (use --force to overwrite)" >&2
    exit 0
fi

# ─── write project file ──────────────────────────────────────────────────
{
    cat <<EOF
---
slug: $SLUG
url: ${GITHUB_URL:-(unknown — fill in)}
captured: $(date +%Y-%m-%d)
source_pdf: projects/$SLUG.pdf
---

EOF
    cat "$TXT_TMP"
} > "$OUT_MD"

echo "wrote: $OUT_MD" >&2

# ─── move PDF ────────────────────────────────────────────────────────────
if [[ "$(readlink -f "$PDF_IN")" != "$(readlink -f "$OUT_PDF")" ]]; then
    if [[ -f "$OUT_PDF" && $FORCE -eq 0 ]]; then
        echo "warn: $OUT_PDF already exists; input PDF NOT moved" >&2
    else
        mv "$PDF_IN" "$OUT_PDF"
        echo "moved PDF: $OUT_PDF" >&2
    fi
fi

# ─── append to INDEX.md ──────────────────────────────────────────────────
if [[ ! -f "$INDEX" ]]; then
    cat > "$INDEX" <<'EOF'
# Projects index

GitHub repositories tracked alongside the paper corpus. Each project's README
is captured verbatim (extracted from a print-to-PDF) in `projects/<slug>.md`,
with the source PDF preserved at `projects/<slug>.pdf`.

Projects are NOT run through Qwen summarisation, paper chunking, or image
rendering — the README IS the summary; the project author's framing is the
data.

To grab the source when needed, the URL in each entry is the canonical handle.

## Projects

EOF
fi

cat >> "$INDEX" <<EOF

### $SLUG
- URL: ${GITHUB_URL:-(unknown — fill in)}
- Captured: $(date +%Y-%m-%d)
- README: projects/$SLUG.md
- Source PDF: projects/$SLUG.pdf
EOF

echo "appended index entry to $INDEX" >&2
