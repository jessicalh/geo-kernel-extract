#!/usr/bin/env python3
"""
Example: load a trajectory H5, find all carbon atoms near ring 0,
plot their Biot-Savart T0 time series.

Usage:
    python example_ring_proximity.py /path/to/md_analysis.h5
"""

import sys
from pathlib import Path

import numpy as np

# Add parent to path for direct execution
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from trajectory_reader import load_analysis, Element


def main():
    if len(sys.argv) < 2:
        print("Usage: python example_ring_proximity.py <path_to_md_analysis.h5>")
        sys.exit(1)

    h5_path = sys.argv[1]
    protein = load_analysis(h5_path)
    print(protein.summary())
    print()

    if not protein.rings:
        print("No rings in this protein.")
        sys.exit(0)

    ring = protein.rings[0]
    print(f"Target ring: {ring} in residue {ring.residue}")
    print(f"  Ring atom indices: {ring.atom_indices}")
    print()

    # Find carbon atoms within 8 Angstroms of ring center at frame 0
    ring_center_0 = ring.center(0)
    cutoff = 8.0
    carbons_near_ring = []

    for atom in protein.atoms:
        if atom.element != Element.C:
            continue
        # Skip ring atoms themselves
        if atom.index in ring.atom_indices:
            continue
        dist = np.linalg.norm(atom.positions[0] - ring_center_0)
        if dist < cutoff:
            carbons_near_ring.append((atom, dist))

    carbons_near_ring.sort(key=lambda x: x[1])
    print(f"Carbon atoms within {cutoff} A of ring {ring.index} center (frame 0):")
    for atom, dist in carbons_near_ring:
        print(f"  {atom}  dist={dist:.2f} A")
    print()

    if not carbons_near_ring:
        print("No carbon atoms found near ring. Try a larger cutoff.")
        sys.exit(0)

    # Plot BS T0 time series for the PheBenzene ring type (index 0)
    # for the nearest carbon atoms
    ring_type_idx = ring.type_index.value  # index into the 8-type axis

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(12, 6))

        n_to_plot = min(6, len(carbons_near_ring))
        times_ns = protein.frame_times / 1000.0  # ps -> ns

        for atom, dist in carbons_near_ring[:n_to_plot]:
            bs_t0 = atom.ring_current.bs_T0_per_type[:, ring_type_idx]
            label = f"{atom.pdb_name} ({atom.residue.name}{atom.residue.sequence_number}) d={dist:.1f}A"
            ax.plot(times_ns, bs_t0, label=label, alpha=0.8, linewidth=0.5)

        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(f"BS T0 ({ring.type_name} component, ppm)")
        ax.set_title(
            f"{protein.protein_id}: Biot-Savart T0 for carbons near "
            f"{ring.type_name} ring (residue {ring.residue.name}{ring.residue.sequence_number})"
        )
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        out_path = Path(h5_path).parent / "example_bs_t0_near_ring.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Plot saved to: {out_path}")
        plt.close(fig)

    except ImportError:
        print("matplotlib not available. Printing numerical summary instead.")
        print()
        for atom, dist in carbons_near_ring[:6]:
            bs_t0 = atom.ring_current.bs_T0_per_type[:, ring_type_idx]
            print(
                f"  {atom.pdb_name} ({atom.residue.name}{atom.residue.sequence_number}): "
                f"mean={np.mean(bs_t0):.6f}, std={np.std(bs_t0):.6f}, "
                f"range=[{np.min(bs_t0):.6f}, {np.max(bs_t0):.6f}]"
            )

    # Also demonstrate frame view access
    print()
    print("Frame view demonstration:")
    fv = protein.frame(0)
    print(f"  {fv}")
    for atom, dist in carbons_near_ring[:3]:
        fav = fv[atom.index]
        print(f"  Atom {atom.index} ({atom.pdb_name}) position at t=0: {fav.position}")


if __name__ == "__main__":
    main()
