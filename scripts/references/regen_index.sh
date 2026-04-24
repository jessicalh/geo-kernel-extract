#!/usr/bin/env bash
#
# regen_index.sh — regenerate references-meta/INDEX.md from keyword files.
#
# Produces a topic → paper-basename cross-reference by scanning every
# <basename>-keywords.txt under references-meta/ and grouping papers by keyword.
# Keywords are case-folded for grouping but displayed in their first-seen form.
#
# Usage:
#   regen_index.sh
#
# Caveats:
#   - Keywords are free-form, not a controlled vocabulary. Near-synonyms
#     ("ring current" vs "ring-current") appear as separate topic entries.
#     Scan the output and consolidate manually if it matters for a given
#     writing task.
#   - This script is idempotent and safe to run any time. Overwrites
#     references-meta/INDEX.md.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
META_DIR="$PROJECT_ROOT/references-meta"
OUT_FILE="$META_DIR/INDEX.md"

if [[ ! -d "$META_DIR" ]]; then
    echo "error: $META_DIR does not exist" >&2
    exit 1
fi

# ─── collect (keyword, basename) pairs ────────────────────────────────────
# Each keyword file: lines starting with '#' are comments (convention is '# <cite-key>').
# Blank lines ignored. Every other line is one keyword.

TMP_PAIRS="$(mktemp)"
trap 'rm -f "$TMP_PAIRS"' EXIT

NPAPERS=0
for kw_file in "$META_DIR"/*-keywords.txt; do
    [[ -e "$kw_file" ]] || continue
    # derive basename from filename
    fname="$(basename "$kw_file")"
    basename="${fname%-keywords.txt}"
    NPAPERS=$((NPAPERS + 1))

    # collect keywords
    grep -v '^#' "$kw_file" | grep -v '^$' | while IFS= read -r kw; do
        # trim whitespace
        kw="$(echo "$kw" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
        [[ -z "$kw" ]] && continue
        # emit lowercased key + original + basename, tab-separated
        printf '%s\t%s\t%s\n' "$(echo "$kw" | tr '[:upper:]' '[:lower:]')" "$kw" "$basename"
    done >> "$TMP_PAIRS"
done

if [[ "$NPAPERS" -eq 0 ]]; then
    echo "error: no keyword files found in $META_DIR" >&2
    exit 1
fi

# ─── emit INDEX.md ────────────────────────────────────────────────────────
{
    echo "# references-meta/INDEX.md"
    echo ""
    echo "Auto-generated keyword → paper index. Regenerate with:"
    echo ""
    echo '```'
    echo "scripts/references/regen_index.sh"
    echo '```'
    echo ""
    echo "Source: every \`*-keywords.txt\` file under \`references-meta/\`."
    echo "Papers: $NPAPERS. Generated: $(date -u +'%Y-%m-%d %H:%M UTC')."
    echo ""
    echo "**Caveat:** keywords are free-form, not a controlled vocabulary. Near-synonyms appear as separate entries. When snarfing a topic, grep more than one phrasing."
    echo ""
    echo "---"
    echo ""

    # Sort by lowercased key, then unique (key, basename) per key.
    # For each unique lowercased keyword, print the original form (first-seen) and the list of papers.
    sort "$TMP_PAIRS" | awk -F'\t' '
    {
        key = $1; orig = $2; paper = $3
        if (key != last_key) {
            if (NR > 1) print ""
            printf("**%s**\n", orig)
            last_key = key
            delete seen
        }
        if (!(paper in seen)) {
            printf("- %s\n", paper)
            seen[paper] = 1
        }
    }
    '

    echo ""
    echo ""
    echo "---"
    echo ""
    echo "## Papers indexed"
    echo ""
    ls "$META_DIR"/*-keywords.txt 2>/dev/null | while read -r f; do
        fname="$(basename "$f")"
        basename="${fname%-keywords.txt}"
        echo "- \`$basename\`"
    done | sort
} > "$OUT_FILE"

echo "wrote $OUT_FILE"
echo "  papers indexed: $NPAPERS"
echo "  unique keywords: $(awk -F'\t' '{print $1}' "$TMP_PAIRS" | sort -u | wc -l | tr -d ' ')"
