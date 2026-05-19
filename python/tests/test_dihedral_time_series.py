"""SDK round-trip tests for /trajectory/dihedral_time_series/.

Verifies that:
  - load_trajectory() exposes the dihedral group as
    TrajectoryData.dihedrals (DihedralTimeSeriesGroup);
  - per-frame channels (phi, psi, omega, omega_deviation, chi,
    rama_region) round-trip with correct shape and dtype;
  - static per-residue masks (chi_exists, omega_is_xpro, is_glycine,
    is_proline, is_pre_proline, residue_terminal_state) round-trip;
  - chain_id_per_residue variable-length strings decode correctly;
  - residue_index_per_atom (N,) atom→residue lookup is the
    SDK-side broadcast bridge for the movie viewer;
  - convention-pin attributes (angle_units, periodicity, value_range,
    angle_convention, chain_break_policy, rama_region_legend,
    rama_region_boundaries, chi_symmetry_caveats, source_attached_policy,
    chunking_policy) propagate verbatim;
  - missing group (no dihedral TR ran) returns None.

Synthesises trajectory.h5 from scratch with h5py — mirrors
DihedralTimeSeriesTrajectoryResult::WriteH5Group exactly.
"""

import h5py
import numpy as np
import pytest

from nmr_extract import (
    load_trajectory,
    DihedralTimeSeriesGroup,
)


N_ATOMS    = 12
N_RESIDUES = 4
N_FRAMES   = 5


def _write_required_traj_root(f: h5py.File) -> None:
    f.attrs["protein_id"] = "TEST_PROTEIN"
    f.attrs["n_atoms"]    = N_ATOMS
    f.attrs["finalized"]  = True

    traj = f.create_group("trajectory")
    frames = traj.create_group("frames")
    frames.attrs["n_frames"] = N_FRAMES
    times = np.arange(N_FRAMES, dtype=np.float64) * 0.1
    frames.create_dataset("time_ps", data=times)
    frames.create_dataset("original_index",
                          data=np.arange(N_FRAMES, dtype=np.uint64))

    pos = traj.create_group("positions")
    pos.attrs["n_atoms"]     = N_ATOMS
    pos.attrs["n_frames"]    = N_FRAMES
    pos.attrs["result_name"] = "PositionsTimeSeriesTrajectoryResult"
    pos.attrs["finalized"]   = True
    pos.create_dataset("xyz",
                       data=np.zeros((N_ATOMS * N_FRAMES, 3), dtype=np.float64))


def _write_dihedral_time_series(f: h5py.File) -> dict:
    """Emit /trajectory/dihedral_time_series with a recognisable pattern.

    Convention is residue-major flat: residue i, frames 0..T-1 contiguous.
    Returns a dict of expected payloads for round-trip assertion.
    """
    grp = f.create_group("/trajectory/dihedral_time_series")

    # Provenance + convention attrs (must match the C++ emit verbatim)
    grp.attrs["result_name"] = "DihedralTimeSeriesTrajectoryResult"
    grp.attrs["n_residues"]  = N_RESIDUES
    grp.attrs["n_atoms"]     = N_ATOMS
    grp.attrs["n_frames"]    = N_FRAMES
    grp.attrs["finalized"]   = True
    grp.attrs["angle_units"]   = "radians"
    grp.attrs["periodicity"]   = "2pi"
    grp.attrs["value_range"]   = "[-pi, pi] (atan2 closed range)"
    grp.attrs["angle_convention"] = (
        "IUPAC signed dihedral atan2(y,x). Note: DSSP uses NEGATED "
        "convention (phi_DSSP = -phi_IUPAC).")
    grp.attrs["chain_break_policy"] = (
        "BackboneConnected queries the LegacyAmber bond graph directly")
    grp.attrs["omega_deviation_policy"] = (
        "WrapPi(omega - pi); emitted for X->Pro too")
    grp.attrs["rama_region_legend"] = (
        "0=unassigned, 1=alphaR, 2=beta, 3=alphaL, 4=PPII, 5=other")
    grp.attrs["rama_region_boundaries"] = (
        "Resolution order: alphaR -> alphaL -> PPII -> beta -> other")
    grp.attrs["chi_symmetry_caveats"] = (
        "PHE/TYR chi2 ring-flip; ASP/GLU carboxylate; "
        "TRP/HIS chi2 NOT symmetric")
    grp.attrs["residue_terminal_state_legend"] = (
        "0=internal, 1=n_terminus, 2=c_terminus, 3=n_and_c_terminus, 4=unknown")
    grp.attrs["source_attached_policy"] = "always_attached"
    grp.attrs["chunking_policy"] = (
        "Per-residue datasets chunked as {R, min(T, 64)}")
    grp.attrs["residue_axis"] = "protein_residue_index"
    grp.attrs["atom_axis"]    = "protein_atom_index"
    grp.attrs["source"]       = "positions + AminoAcidType.chi_angles"

    # Per-frame (R, T) datasets — fill with recognisable patterns.
    # phi[ri, t] = ri/10 + t/100  (small distinct values)
    phi = np.fromfunction(lambda ri, t: ri / 10.0 + t / 100.0,
                          (N_RESIDUES, N_FRAMES), dtype=np.float64)
    psi = -phi
    omega = np.full((N_RESIDUES, N_FRAMES), np.pi, dtype=np.float64)
    omega_deviation = np.zeros((N_RESIDUES, N_FRAMES), dtype=np.float64)
    # N-terminus phi NaN; C-terminus psi+omega NaN
    phi[0, :]                = np.nan
    psi[N_RESIDUES - 1, :]   = np.nan
    omega[N_RESIDUES - 1, :] = np.nan
    omega_deviation[N_RESIDUES - 1, :] = np.nan

    grp.create_dataset("phi", data=phi)
    grp.create_dataset("psi", data=psi)
    grp.create_dataset("omega", data=omega)
    grp.create_dataset("omega_deviation", data=omega_deviation)

    rama_region = np.full((N_RESIDUES, N_FRAMES), 5, dtype=np.uint8)
    rama_region[0, :] = 0  # N-terminus unassigned (NaN phi)
    grp.create_dataset("rama_region", data=rama_region)

    chi = np.full((N_RESIDUES, N_FRAMES, 4), np.nan, dtype=np.float64)
    chi[1, :, 0] = 1.5  # residue 1 has chi1 only
    chi[2, :, :2] = -0.5  # residue 2 has chi1+chi2
    grp.create_dataset("chi", data=chi)

    # Static masks
    chi_exists = np.zeros((N_RESIDUES, 4), dtype=np.uint8)
    chi_exists[1, 0] = 1
    chi_exists[2, :2] = 1
    grp.create_dataset("chi_exists", data=chi_exists)

    omega_is_xpro = np.array([0, 1, 0, 0], dtype=np.uint8)
    is_glycine    = np.array([1, 0, 0, 0], dtype=np.uint8)
    is_proline    = np.array([0, 0, 1, 0], dtype=np.uint8)
    is_pre_proline= np.array([0, 1, 0, 0], dtype=np.uint8)
    residue_terminal_state = np.array([1, 0, 0, 2], dtype=np.uint8)
    grp.create_dataset("omega_is_xpro", data=omega_is_xpro)
    grp.create_dataset("is_glycine",   data=is_glycine)
    grp.create_dataset("is_proline",   data=is_proline)
    grp.create_dataset("is_pre_proline", data=is_pre_proline)
    grp.create_dataset("residue_terminal_state",
                       data=residue_terminal_state)

    # Variable-length strings for chain_id_per_residue
    chain_ids = np.array(["A", "A", "A", "A"], dtype=h5py.string_dtype())
    grp.create_dataset("chain_id_per_residue", data=chain_ids)

    # Per-atom lookup: 3 atoms per residue
    ria = np.repeat(np.arange(N_RESIDUES, dtype=np.int32), 3)
    grp.create_dataset("residue_index_per_atom", data=ria)

    # Per-frame metadata
    frame_indices = np.arange(N_FRAMES, dtype=np.uint64)
    frame_times = np.arange(N_FRAMES, dtype=np.float64) * 0.1
    source_attached = np.ones((N_FRAMES,), dtype=np.uint8)
    grp.create_dataset("frame_indices", data=frame_indices)
    grp.create_dataset("frame_times", data=frame_times)
    grp.create_dataset("source_attached_per_frame", data=source_attached)

    return {
        "phi": phi, "psi": psi, "omega": omega,
        "omega_deviation": omega_deviation, "chi": chi,
        "rama_region": rama_region,
        "chi_exists": chi_exists,
        "omega_is_xpro": omega_is_xpro,
        "is_glycine": is_glycine, "is_proline": is_proline,
        "is_pre_proline": is_pre_proline,
        "residue_terminal_state": residue_terminal_state,
        "chain_ids": ["A", "A", "A", "A"],
        "ria": ria,
        "frame_indices": frame_indices,
        "frame_times": frame_times,
        "source_attached": source_attached,
    }


def test_dihedral_round_trip(tmp_path):
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        expected = _write_dihedral_time_series(f)

    traj = load_trajectory(h5_path)
    assert traj.dihedrals is not None
    g = traj.dihedrals
    assert isinstance(g, DihedralTimeSeriesGroup)

    # Per-frame channel shape + dtype + payload
    assert g.phi.shape == (N_RESIDUES, N_FRAMES)
    assert g.phi.dtype == np.float64
    np.testing.assert_array_equal(g.phi, expected["phi"], strict=False)
    np.testing.assert_array_equal(g.psi, expected["psi"], strict=False)
    np.testing.assert_array_equal(g.omega, expected["omega"], strict=False)
    np.testing.assert_array_equal(g.omega_deviation,
                                  expected["omega_deviation"], strict=False)

    assert g.chi.shape == (N_RESIDUES, N_FRAMES, 4)
    np.testing.assert_array_equal(g.chi, expected["chi"], strict=False)

    assert g.rama_region.shape == (N_RESIDUES, N_FRAMES)
    assert g.rama_region.dtype == np.uint8
    np.testing.assert_array_equal(g.rama_region, expected["rama_region"])


def test_dihedral_static_masks(tmp_path):
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        expected = _write_dihedral_time_series(f)
    traj = load_trajectory(h5_path)
    g = traj.dihedrals

    assert g.chi_exists.shape == (N_RESIDUES, 4)
    np.testing.assert_array_equal(g.chi_exists, expected["chi_exists"])
    np.testing.assert_array_equal(g.omega_is_xpro, expected["omega_is_xpro"])
    np.testing.assert_array_equal(g.is_glycine, expected["is_glycine"])
    np.testing.assert_array_equal(g.is_proline, expected["is_proline"])
    np.testing.assert_array_equal(g.is_pre_proline, expected["is_pre_proline"])
    np.testing.assert_array_equal(g.residue_terminal_state,
                                  expected["residue_terminal_state"])


def test_dihedral_chain_id_decodes(tmp_path):
    """chain_id_per_residue is variable-length strings — must decode to
    Python str for downstream readers."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        _write_dihedral_time_series(f)
    traj = load_trajectory(h5_path)
    g = traj.dihedrals
    assert g.chain_id_per_residue.shape == (N_RESIDUES,)
    for ri in range(N_RESIDUES):
        assert isinstance(g.chain_id_per_residue[ri], str)
        assert g.chain_id_per_residue[ri] == "A"


def test_dihedral_atom_axis_broadcast(tmp_path):
    """residue_index_per_atom (N,) is the SDK-side bridge that lets the
    movie viewer broadcast residue-axis data to atom axis at render time.
    Verify the lookup is correct and the broadcast pattern works."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        expected = _write_dihedral_time_series(f)
    traj = load_trajectory(h5_path)
    g = traj.dihedrals
    assert g.residue_index_per_atom.shape == (N_ATOMS,)
    np.testing.assert_array_equal(g.residue_index_per_atom, expected["ria"])

    # The movie viewer pattern: phi_per_atom[a, t] = phi[ria[a], t]
    phi_per_atom = g.phi[g.residue_index_per_atom, :]
    assert phi_per_atom.shape == (N_ATOMS, N_FRAMES)
    # 3 atoms per residue → all 3 atoms of residue 0 share residue-0 phi
    np.testing.assert_array_equal(phi_per_atom[0:3], g.phi[0:1].repeat(3, axis=0))


def test_dihedral_convention_attrs(tmp_path):
    """Convention-pin attrs must propagate verbatim — downstream
    consumers / advisor read these to understand sign conventions
    and DSSP divergence."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        _write_dihedral_time_series(f)
    traj = load_trajectory(h5_path)
    g = traj.dihedrals
    assert g.angle_units == "radians"
    assert g.periodicity == "2pi"
    assert "[-pi, pi]" in g.value_range
    assert "IUPAC" in g.angle_convention
    assert "DSSP" in g.angle_convention
    assert "bond graph" in g.chain_break_policy
    assert "X->Pro" in g.omega_deviation_policy
    assert "alphaR" in g.rama_region_legend
    assert "Resolution order" in g.rama_region_boundaries
    assert "ring-flip" in g.chi_symmetry_caveats
    assert "NOT symmetric" in g.chi_symmetry_caveats
    assert g.source_attached_policy == "always_attached"
    assert "min(T, 64)" in g.chunking_policy


def test_dihedral_missing_group(tmp_path):
    """A trajectory.h5 without the dihedral TR group → dihedrals is None."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        # NO dihedral_time_series group written.
    traj = load_trajectory(h5_path)
    assert traj.dihedrals is None


def test_dihedral_frame_metadata(tmp_path):
    """frame_indices, frame_times, source_attached_per_frame all (T,)."""
    h5_path = tmp_path / "trajectory.h5"
    with h5py.File(h5_path, "w") as f:
        _write_required_traj_root(f)
        expected = _write_dihedral_time_series(f)
    traj = load_trajectory(h5_path)
    g = traj.dihedrals
    np.testing.assert_array_equal(g.frame_indices, expected["frame_indices"])
    np.testing.assert_array_equal(g.frame_times, expected["frame_times"])
    np.testing.assert_array_equal(g.source_attached_per_frame,
                                  expected["source_attached"])
    # source_attached_per_frame should be trivially all-1
    assert (g.source_attached_per_frame == 1).all()
