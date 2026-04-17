#!/usr/bin/env python3
"""
Example: load an analysis H5, find carbon atoms near ring 0,
plot their Biot-Savart T0 time series. Also demonstrate per-ring data access.

Usage:
    python -m trajectory_reader_v2.example /path/to/analysis.h5

    Or from the parent directory:
    python -m trajectory_reader_v2.example \
        /shared/2026Thesis/fleet_calibration-working/1B1V_4292/analysis_output/1B1V_4292_analysis.h5
"""

import sys
import numpy as np

# Add parent directory to path so we can import the package
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from trajectory_reader_v2 import load_analysis, Element, AtomRole, RingTypeIndex
from trajectory_reader_v2.model import SphericalTensorView


def main():
    if len(sys.argv) < 2:
        print("Usage: python example.py <analysis.h5>")
        sys.exit(1)

    h5_path = sys.argv[1]

    with load_analysis(h5_path) as protein:
        print(protein)
        print()

        # ----------------------------------------------------------------
        # Part 1: Protein identity -- constant topology
        # ----------------------------------------------------------------

        print("=== Topology ===")
        print(f"  Atoms:    {protein.n_atoms}")
        print(f"  Residues: {protein.n_residues}")
        print(f"  Rings:    {protein.n_rings}")
        print(f"  Bonds:    {protein.n_bonds}")
        print(f"  Frames:   {protein.n_frames}")
        print(f"  Time:     {protein.frame_times[0]:.1f} - {protein.frame_times[-1]:.1f} ps")
        print()

        # Show residue sequence
        print("=== Sequence ===")
        for r in protein.residues:
            print(f"  {r}")
        print()

        # Show rings
        print("=== Aromatic rings ===")
        for ring in protein.rings:
            print(f"  {ring}")
        print()

        # Element breakdown
        for elem in Element:
            count = len(protein.atoms_by_element(elem))
            if count > 0:
                print(f"  {elem.symbol}: {count} atoms")
        print()

        # ----------------------------------------------------------------
        # Part 2: Find carbon atoms near ring 0
        # ----------------------------------------------------------------

        if protein.n_rings == 0:
            print("No aromatic rings in this protein.")
            return

        ring0 = protein.ring(0)
        print(f"=== Ring 0: {ring0} ===")
        print(f"  Ring atoms: {[protein.atom(i).pdb_name for i in ring0.atom_indices]}")
        print(f"  Residue: {ring0.residue}")
        print()

        # Find carbons: compute mean distance to ring center over all frames
        # using the ring's atom positions to determine center.
        # h5py fancy indexing requires sorted indices, so sort them.
        sorted_ring_atoms = sorted(ring0.atom_indices)
        ring_atom_positions = protein.positions.xyz[:, sorted_ring_atoms, :]  # (T, n_ring_atoms, 3)
        ring_centers = ring_atom_positions.mean(axis=1)  # (T, 3)

        # All carbon atoms
        carbons = protein.atoms_by_element(Element.C)
        print(f"=== Carbon atoms near ring 0 (within 8A mean distance) ===")

        near_carbons = []
        for c in carbons:
            # Skip ring's own atoms
            if c.index in ring0.atom_indices:
                continue
            c_pos = protein.positions.xyz[:, c.index, :]  # (T, 3)
            dists = np.linalg.norm(c_pos - ring_centers, axis=1)  # (T,)
            mean_dist = dists.mean()
            if mean_dist < 8.0:
                near_carbons.append((c, mean_dist))

        near_carbons.sort(key=lambda x: x[1])
        for c, d in near_carbons:
            res = c.residue
            print(f"  atom {c.index:3d} {c.pdb_name:4s} ({res.amino_acid.name} {res.sequence_number}) "
                  f"mean_dist={d:.2f} A")
        print()

        # ----------------------------------------------------------------
        # Part 3: Biot-Savart T0 time series for the nearest carbons
        # ----------------------------------------------------------------

        if near_carbons:
            print("=== BS T0 time series for nearest 5 carbons ===")
            bs_T0 = protein.ring_current.bs_T0_per_type  # (T, N, 8) dataset

            for c, d in near_carbons[:5]:
                # BS T0 for ring type 0 (PheBenzene)
                ring_type_idx = ring0.type_index.value
                ts = bs_T0[:, c.index, ring_type_idx]  # (T,)
                print(f"  atom {c.index} {c.pdb_name}: "
                      f"mean={ts.mean():.6f}, std={ts.std():.6f}, "
                      f"min={ts.min():.6f}, max={ts.max():.6f}")
            print()

            # ----------------------------------------------------------------
            # Part 4: BS shielding as SphericalTensor
            # ----------------------------------------------------------------

            c0 = near_carbons[0][0]
            stv = protein.spherical_tensor("ring_current/bs_shielding", atom_index=c0.index)
            print(f"=== BS shielding SphericalTensor for atom {c0.index} {c0.pdb_name} ===")
            print(f"  T0 time series: mean={stv.T0.mean():.6f}, std={stv.T0.std():.6f}")
            print(f"  T2 magnitude time series: mean={stv.T2_magnitude.mean():.6f}")
            print()

        # ----------------------------------------------------------------
        # Part 5: Per-ring K=6 data access
        # ----------------------------------------------------------------

        if protein.per_ring is not None:
            print(f"=== Per-ring data (K={protein.per_ring.K}) ===")
            print(f"  Geometry fields: {protein.per_ring.geometry_fields}")
            print()

            # Pick a backbone amide H atom to show per-ring data
            amide_atoms = protein.find_atoms(role=AtomRole.AmideH)
            if not amide_atoms:
                amide_atoms = protein.find_atoms(element=Element.H)

            if amide_atoms:
                test_atom = amide_atoms[0]
                print(f"  Per-ring data for atom {test_atom.index} ({test_atom.pdb_name}, "
                      f"{test_atom.residue.amino_acid.name} {test_atom.residue.sequence_number}):")
                print()

                # Show frame 0 per-ring breakdown
                t = 0
                ring_types_at_t = protein.per_ring.ring_type[t, test_atom.index, :]  # (K,)
                geometry_at_t = protein.per_ring.geometry[t, test_atom.index, :, :]   # (K, F)
                bs_T2_at_t = protein.per_ring.bs_T2[t, test_atom.index, :, :]         # (K, 5)

                for k in range(protein.per_ring.K):
                    rt = int(ring_types_at_t[k])
                    if rt < 0:
                        print(f"    ring slot {k}: empty (no ring)")
                        continue
                    try:
                        ring_type = RingTypeIndex(rt)
                        rt_name = ring_type.short_name
                    except ValueError:
                        rt_name = f"type={rt}"

                    geom_vals = geometry_at_t[k, :]
                    bs_t2 = bs_T2_at_t[k, :]
                    t2_mag = np.sqrt(np.sum(bs_t2 ** 2))

                    # Show geometry field values
                    geom_str = ", ".join(
                        f"{name}={geom_vals[j]:.3f}"
                        for j, name in enumerate(protein.per_ring.geometry_fields)
                    )
                    print(f"    ring slot {k}: {rt_name} -- {geom_str}")
                    print(f"      BS T2 = [{', '.join(f'{v:.6f}' for v in bs_t2)}], |T2|={t2_mag:.6f}")
                print()

                # Show T2 magnitude time series for nearest ring
                bs_T2_k0 = protein.per_ring.bs_T2[:, test_atom.index, 0, :]  # (T, 5)
                t2_mag_ts = np.sqrt(np.sum(bs_T2_k0 ** 2, axis=1))  # (T,)
                print(f"    Nearest ring BS |T2| time series: "
                      f"mean={t2_mag_ts.mean():.6f}, std={t2_mag_ts.std():.6f}")
                print()

        # ----------------------------------------------------------------
        # Part 6: Ring geometry time series
        # ----------------------------------------------------------------

        if protein.ring_geometry is not None and protein.n_rings > 0:
            print(f"=== Ring geometry ===")
            print(f"  Fields: {protein.ring_geometry.fields}")
            print(f"  Data shape: {protein.ring_geometry.data.shape}")

            # Show ring 0 geometry at frame 0
            r0_geom = protein.ring_geometry.data[0, 0, :]
            for j, fname in enumerate(protein.ring_geometry.fields):
                print(f"    ring 0, frame 0: {fname} = {r0_geom[j]:.4f}")
            print()

        # ----------------------------------------------------------------
        # Part 7: Charge time series
        # ----------------------------------------------------------------

        if protein.charges is not None:
            c_atoms = protein.atoms_by_element(Element.C)[:3]
            print("=== AIMNet2 charge time series (first 3 carbons) ===")
            for c in c_atoms:
                q = protein.charges.aimnet2_charge[:, c.index]  # (T,)
                print(f"  atom {c.index} {c.pdb_name}: "
                      f"mean={q.mean():.4f}, std={q.std():.4f}")
            print()

        print("Done.")


if __name__ == "__main__":
    main()
