#!/usr/bin/env python3
"""Quote relocated paths in the official Apptainer v1.5.3 installer wrappers.

The engine itself still needs a cache path without whitespace; these fixes do
not remove its separate FUSE session-path limitation. Run only on a fresh tree
created by the pinned installer documented in APPTAINER_ENGINE.txt.
"""
from pathlib import Path
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine", type=Path)
    args = parser.parse_args()
    changes = {
        "bin/apptainer": [
            ('$(/usr/bin/realpath $0)', '$(/usr/bin/realpath "$0")'),
            ('exec $APPTDIR/bin/apptainer "$@"', 'exec "$APPTDIR/bin/apptainer" "$@"'),
        ],
        "x86_64/utils/bin/.wrapper": [('} $REALME "$@"', '} "$REALME" "$@"')],
        "x86_64/libexec/apptainer/bin/.wrapper": [('} $REALME "$@"', '} "$REALME" "$@"')],
    }
    prepared = []
    for relative, replacements in changes.items():
        path = args.engine / relative
        text = path.read_text()
        for old, new in replacements:
            if text.count(old) != 1:
                parser.error(f"Unexpected or already patched wrapper: {relative}: {old}")
            text = text.replace(old, new)
        prepared.append((path, text))
    for path, text in prepared:
        path.write_text(text)


if __name__ == "__main__":
    main()
