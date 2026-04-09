#!/usr/bin/env python3
"""CLI for secondary analysis tools.

Usage:
    cd learn/src
    python -m secondary divergence --config calibration.toml
    python -m secondary ring_strata --config calibration.toml
    python -m secondary kernel_structure --config calibration.toml
    python -m secondary ablation --config calibration.toml
    python -m secondary --all --config calibration.toml
    python -m secondary loader --config calibration.toml  # verify loader
"""

from __future__ import annotations

import argparse
import sys

from mutation_set.config import load_config
from .loader import setup_sdk


TOOLS = {
    "divergence": "Classical vs MOPAC T2 divergence",
    "ring_strata": "Ring-type conditional analysis",
    "kernel_structure": "Inter-kernel T2 structure (DFT-free)",
    "ablation": "Per-ring-type R² ablation",
    "scalar_ablation": "Scalar group interaction analysis (linear→equivariant gap)",
    "viewer_export": "Export per-atom annotations as JSON for Qt/VTK viewer",
    "export_for_r": "Export kernel matrix + scalars + metadata for R analysis",
    "element_physics": "Per-element kernel physics decomposition (Boyd/Sahakyan test)",
    "loader": "Verify loader: print ring-type distribution",
}


def _run_loader(cfg, max_proteins):
    """Diagnostic: iterate all proteins, print ring-type distribution."""
    from .loader import iter_proteins, RING_TYPE_NAMES, N_RING_TYPES

    rt_counts = [0] * N_RING_TYPES
    n_proteins = 0
    n_atoms = 0
    n_with_mopac = 0
    strata_counts = {}

    for rec in iter_proteins(cfg, max_proteins=max_proteins):
        n_proteins += 1
        n_atoms += len(rec.matched_idx)
        if rec.has_mopac:
            n_with_mopac += 1
        for rt in rec.protein_ring_types:
            rt_counts[rt] += 1

        for mask in rec.atom_ring_masks:
            from .loader import assign_stratum
            s = assign_stratum(int(mask), cfg.secondary.strata)
            if s:
                strata_counts[s] = strata_counts.get(s, 0) + 1

    print(f"\nLoader summary:")
    print(f"  Proteins: {n_proteins}")
    print(f"  Matched atoms: {n_atoms}")
    print(f"  With MOPAC: {n_with_mopac}")
    print(f"\n  Ring type distribution (proteins containing type):")
    for i in range(N_RING_TYPES):
        print(f"    {RING_TYPE_NAMES[i]:15s}  {rt_counts[i]:4d} proteins")
    print(f"\n  Atom strata:")
    for s in cfg.secondary.strata:
        print(f"    {s:15s}  {strata_counts.get(s, 0):6d} atoms")


def main():
    parser = argparse.ArgumentParser(
        description="Secondary analysis toolkit for NMR shielding calibration")
    parser.add_argument("tool", nargs="?", choices=list(TOOLS.keys()),
                        help="Which tool to run")
    parser.add_argument("--all", action="store_true",
                        help="Run all 4 analysis tools")
    parser.add_argument("--config", type=str, default="calibration.toml")
    parser.add_argument("--max-proteins", type=int, default=0,
                        help="Limit protein count (0 = all)")
    parser.add_argument("--protein", type=str, default="",
                        help="Single protein ID (for viewer_export)")
    args = parser.parse_args()

    if not args.tool and not args.all:
        parser.print_help()
        print(f"\nAvailable tools:")
        for name, desc in TOOLS.items():
            print(f"  {name:20s}  {desc}")
        sys.exit(1)

    cfg = load_config(args.config)
    if args.max_proteins > 0:
        from dataclasses import replace
        cfg = replace(cfg, data=replace(cfg.data,
                                        max_proteins=args.max_proteins))
    setup_sdk(cfg)

    if args.tool == "loader":
        _run_loader(cfg, args.max_proteins)
        return

    tools_to_run = []
    if args.all:
        tools_to_run = ["divergence", "ring_strata", "kernel_structure",
                        "ablation", "scalar_ablation"]
    elif args.tool:
        tools_to_run = [args.tool]

    for tool_name in tools_to_run:
        print(f"\n{'=' * 60}")
        print(f"  {tool_name}: {TOOLS[tool_name]}")
        print(f"{'=' * 60}\n")

        if tool_name == "divergence":
            from .divergence import run
            run(cfg, args.max_proteins)
        elif tool_name == "ring_strata":
            from .ring_strata import run
            run(cfg, args.max_proteins)
        elif tool_name == "kernel_structure":
            from .kernel_structure import run
            run(cfg, args.max_proteins)
        elif tool_name == "ablation":
            from .ablation import run
            run(cfg, args.max_proteins)
        elif tool_name == "scalar_ablation":
            from .scalar_ablation import run
            run(cfg, args.max_proteins)
        elif tool_name == "viewer_export":
            from .viewer_export import run
            run(cfg, args.max_proteins, args.protein)
        elif tool_name == "export_for_r":
            from .export_for_r import run
            run(cfg, args.max_proteins)
        elif tool_name == "element_physics":
            from .element_physics import run
            run(cfg, args.max_proteins)


if __name__ == "__main__":
    main()
