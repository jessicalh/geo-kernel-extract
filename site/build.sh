#!/usr/bin/env bash
# Build the static documentation site.
# Copies all generated docs into site/ so the result is fully
# self-contained — no symlinks, servable from anywhere.
#
# Usage:  bash site/build.sh          (from repo root)
#         bash build.sh               (from site/)

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SITE="$REPO/site"

echo "Building docs site in $SITE"

# ── Diagrams (PNG thumbnails + SVG full-size) ────────
mkdir -p "$SITE/img/diagrams"
cp "$REPO/doc/diagrams/"*.png "$SITE/img/diagrams/"
cp "$REPO/doc/diagrams/"*.svg "$SITE/img/diagrams/"

# ── C++ library Doxygen ──────────────────────────────
echo "  Copying C++ API docs (doxygen)..."
rm -rf "$SITE/api"
cp -r "$REPO/doc/generated/doxygen" "$SITE/api"

# ── Viewer Doxygen ───────────────────────────────────
echo "  Copying viewer API docs (doxygen)..."
rm -rf "$SITE/viewer"
cp -r "$REPO/ui/doc/generated/doxygen" "$SITE/viewer"

# ── Python SDK Sphinx ────────────────────────────────
echo "  Copying Python SDK docs (sphinx)..."
rm -rf "$SITE/python"
cp -r "$REPO/python/doc/_build" "$SITE/python"

# ── Summary ──────────────────────────────────────────
echo ""
echo "Done. Contents:"
du -sh "$SITE/api" "$SITE/viewer" "$SITE/python" "$SITE/img" 2>/dev/null
echo ""
echo "Serve with:  python3 -m http.server -d $SITE 8000"
