#!/usr/bin/env bash
# lit-map crawler — history-first broad hit. Ask every paper, with respect to ONE area:
#   does it DEAL WITH the area, or does its intro/discussion TOUCH ON it? fail early on NO.
# Capture HOW (foregrounding roots/history) + foundational refs it CITES (dead-zone leads).
# Engaged papers -> living .md ; every verdict (incl NO) -> .jsonl (resume ledger + dead-zone math).
#
# Usage:   lit-map/crawl.sh <slug>
#   reads  lit-map/questions/<slug>/area.txt     the broad area (one phrase)
#          lit-map/questions/<slug>/source.txt    summary|fulltext   (optional, default summary)
#   writes lit-map/results/<slug>.{md,jsonl}
#
# Env: LIMIT (pilot cap, 0=all) · PAPERS (space-sep basenames to restrict to) ·
#      SOURCE · MAXTOK (250000) · MAXTIME (900) · ENDPOINT · MODEL
# Resumable: a paper already in the jsonl is skipped. Ctrl-C safe.
set -uo pipefail
ROOT=/shared/2026Thesis/nmr-shielding
META="$ROOT/references-meta"; TEXT="$ROOT/references-text"; LM="$ROOT/lit-map"
ENDPOINT="${ENDPOINT:-http://localhost:8199/v1/chat/completions}"
MODEL="${MODEL:-Qwen3.6-35B-A3B-UD-Q4_K_XL.gguf}"
MAXTOK="${MAXTOK:-250000}"; MAXTIME="${MAXTIME:-900}"; LIMIT="${LIMIT:-0}"

slug="${1:?usage: crawl.sh <slug>}"
QDIR="$LM/questions/$slug"
[ -f "$QDIR/area.txt" ] || { echo "ERROR: no $QDIR/area.txt" >&2; exit 1; }
AREA="$(< "$QDIR/area.txt")"
SOURCE="${SOURCE:-$(cat "$QDIR/source.txt" 2>/dev/null || echo summary)}"

mkdir -p "$LM/results"
OUT_MD="$LM/results/$slug.md"; OUT_JSONL="$LM/results/$slug.jsonl"

SYS="You are mapping a thesis introduction that traces the ROOTS AND HISTORY of: ${AREA}.
Judge the paper below ONLY with respect to that area:
(1) DEALS-WITH - is the area a subject of this paper's own work, theory, or results?
(2) TOUCHES-ON - if not, does the paper's INTRODUCTION or DISCUSSION engage the area (its history, origin, physical mechanism, theory, or application) as background, motivation, comparison, or implication?
If NEITHER, output exactly: NO
Otherwise output exactly these three lines and nothing else:
ENGAGES: one of {deals-with | touches-intro | touches-discussion}
HOW: 1 to 3 sentences, grounded in the paper, on what it says about the area; foreground anything about its ROOTS or HISTORY (who, when, which model or idea)
CITES: foundational references it points to for the area, as author-year, semicolon-separated; or none
Be terse and factual; do not invent. If unsure whether it engages the area, lean NO."

if [ ! -f "$OUT_MD" ]; then
  { printf '# Lit-map (history-first): %s\n\n' "$slug"
    printf '**Area:** %s\n\n' "$AREA"
    printf '**Source:** %s  ·  **Model:** %s  ·  thinking on, max_tokens=%s\n\n_Engaged papers only; full verdicts incl. NO in the .jsonl._\n\n---\n' "$SOURCE" "$MODEL" "$MAXTOK"; } > "$OUT_MD"
fi

declare -A done=()
if [ -f "$OUT_JSONL" ]; then while IFS= read -r bn; do [ -n "$bn" ] && done["$bn"]=1; done < <(jq -r '.basename' "$OUT_JSONL" 2>/dev/null); fi

if [ -n "${PAPERS:-}" ]; then read -r -a papers <<< "$PAPERS"
else mapfile -t papers < <(ls "$META"/*-summary.txt 2>/dev/null | sed 's#.*/##; s#-summary\.txt$##' | sort); fi
total=${#papers[@]}
echo "[$slug] source=$SOURCE papers=$total done=${#done[@]} limit=$LIMIT"

n=0; processed=0
for b in "${papers[@]}"; do
  n=$((n+1)); [ -n "${done[$b]:-}" ] && continue
  if [ "$SOURCE" = fulltext ]; then paperfile="$(mktemp)"; cat $(ls -v "$TEXT/$b"-text-*.txt 2>/dev/null) > "$paperfile" 2>/dev/null
  else paperfile="$META/$b-summary.txt"; fi
  if [ ! -s "$paperfile" ]; then echo "[$n/$total] $b: NO SOURCE, skip"; [ "$SOURCE" = fulltext ] && rm -f "$paperfile"; continue; fi
  payload="$(jq -n --arg m "$MODEL" --arg s "$SYS" --argjson mt "$MAXTOK" --rawfile p "$paperfile" \
    '{model:$m,max_tokens:$mt,stream:false,messages:[{role:"system",content:$s},{role:"user",content:("PAPER:\n"+$p)}]}')"
  [ "$SOURCE" = fulltext ] && rm -f "$paperfile"
  resp="$(curl -sS --max-time "$MAXTIME" "$ENDPOINT" -H 'Content-Type: application/json' -d "$payload" 2>/dev/null)"
  content="$(printf '%s' "$resp" | jq -r '.choices[0].message.content // empty' 2>/dev/null)"
  if [ -z "$content" ]; then resp="$(curl -sS --max-time "$MAXTIME" "$ENDPOINT" -H 'Content-Type: application/json' -d "$payload" 2>/dev/null)"; content="$(printf '%s' "$resp" | jq -r '.choices[0].message.content // empty' 2>/dev/null)"; fi
  finish="$(printf '%s' "$resp" | jq -r '.choices[0].finish_reason // "error"' 2>/dev/null)"
  gentok="$(printf '%s' "$resp" | jq -r '.usage.completion_tokens // 0' 2>/dev/null)"
  if [ -z "$content" ]; then echo "[$n/$total] $b: EMPTY after retry — leaving unrecorded for re-run"; continue; fi
  if printf '%s' "$content" | grep -qiE '^[[:space:]]*NO[[:space:]]*$'; then engages="no"
  else engages="$(printf '%s\n' "$content" | sed -n 's/^ENGAGES:[[:space:]]*//p' | head -1)"; [ -z "$engages" ] && engages="yes"; fi
  jq -nc --arg b "$b" --arg e "$engages" --arg f "$finish" --argjson g "${gentok:-0}" --arg c "$content" \
    '{basename:$b,engages:$e,finish:$f,gen_tokens:$g,content:$c}' >> "$OUT_JSONL"
  [ "$engages" != no ] && { printf '\n## %s\n\n' "$b"; printf '%s\n' "$content"; } >> "$OUT_MD"
  done["$b"]=1; processed=$((processed+1))
  echo "[$n/$total] $b -> $engages (${gentok}tok)"
  if [ "$LIMIT" -gt 0 ] && [ "$processed" -ge "$LIMIT" ]; then echo "LIMIT reached"; break; fi
done
echo "[$slug] processed=$processed done=${#done[@]} engaged=$(grep -c '^## ' "$OUT_MD" 2>/dev/null || echo 0)"
