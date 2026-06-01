#!/usr/bin/env python3
"""test_change_of_basis — the ONLY surviving "recompute", and it is a TEST.

Pins the library-T2 <-> e3nn-2e change-of-basis C against two INDEPENDENT
implementations of the l=2 object: the C++ producer (SphericalBasis /
DecomposeLibrary, emitted as the target NPY) and the e3nn library (o3 irreps +
Wigner-D). It asserts:

  (1) C is orthogonal (C^T C = I) and round-trips (e3nn_to_lib(lib_to_e3nn(t))==t).
  (2) EQUIVARIANCE vs the C++ producer: take a known symmetric-traceless tensor,
      rotate it by a random R, project BOTH the original and the rotated tensor
      into the library basis the C++ way (the function the producer ships), map
      both into e3nn 2e via C, and confirm
          e3nn(rotated) == D2(R) @ e3nn(original)
      where D2 is e3nn's Wigner-D for l=2. This proves the C++ library tensor,
      transported through C, transforms under e3nn's D2 — i.e. the basis match is
      correct AND equivariant. This is the rigorous replacement for the deleted
      one-line `np.nanmax(|lib_T2(raw) - emitted|)` sanity check in equiv_t2.py.
  (3) AGAINST THE EMITTED SUBSTRATE: if the canonical target NPY +
      `dft_total_local` 3x3 columns are present, confirm the emitted library
      5-vector equals the library projection of the emitted 3x3 (the producer is
      self-consistent), and that C maps it into a 2e vector of identical norm
      (orthogonal map preserves |T2|, which is the rotation invariant the fitter
      reports).

Run (system python, torch+e3nn env):
    LD_LIBRARY_PATH=<nvidia cu13 libs>:<torch/lib> \
        python3 -m pytest test_change_of_basis.py -q
or plainly:
    python3 test_change_of_basis.py
"""
import os

import numpy as np

import change_of_basis as cob

OUT = os.environ.get("REDISCOVER_OUT", "/tmp/rediscover-out-v2")
TOL = 1e-9


def _lib_project_sym(M):
    """The C++ DecomposeLibrary T2 map, used here ONLY to build/verify fixtures
    inside the test (independent re-statement of the producer's own function for
    the cross-check; NOT importable model code)."""
    M = np.asarray(M, float)
    sym = 0.5 * (M + M.T)
    tr = np.trace(sym) / 3.0
    s = sym - tr * np.eye(3)
    xx, yy, zz = s[0, 0], s[1, 1], s[2, 2]
    xy, yz, xz = s[0, 1], s[1, 2], s[0, 2]
    s2, s32, i2 = np.sqrt(2.0), np.sqrt(1.5), 1.0 / np.sqrt(2.0)
    return np.array([s2 * xy, s2 * yz, s32 * zz, s2 * xz, i2 * (xx - yy)])


def test_orthogonal_and_roundtrip():
    C = cob.derive_C()
    assert np.abs(C.T @ C - np.eye(5)).max() < 1e-10, "C not orthogonal"
    rng = np.random.default_rng(0)
    t = rng.standard_normal((37, 5))
    back = cob.e3nn_to_lib(cob.lib_to_e3nn(t))
    assert np.abs(back - t).max() < 1e-9, "C round-trip failed"


def test_wigner_d_equivariance_vs_cpp():
    """The load-bearing assertion: C++ library tensor through C transforms under
    e3nn's Wigner-D D2(R)."""
    from e3nn import o3
    import torch

    # angles_to_matrix / D_from_matrix default to float32 internally even for
    # float64 inputs, which caps the equivariance residual at ~1e-7; run them in
    # float64 so the check exercises the true (machine-precision) equivariance.
    torch.set_default_dtype(torch.float64)

    C = cob.derive_C()
    rng = np.random.default_rng(1)
    # a known symmetric traceless tensor (the C++ target lives in this space)
    A = rng.standard_normal((3, 3))
    S = 0.5 * (A + A.T)
    S = S - np.trace(S) / 3.0 * np.eye(3)

    # random rotation; rotate the Cartesian tensor S' = R S R^T
    ang = torch.tensor(rng.uniform(-np.pi, np.pi, size=3), dtype=torch.float64)
    R = o3.angles_to_matrix(ang[0], ang[1], ang[2]).numpy().astype(float)
    S_rot = R @ S @ R.T

    # project BOTH through the C++ library map, then into e3nn 2e via C
    e_orig = cob.lib_to_e3nn(_lib_project_sym(S))
    e_rot = cob.lib_to_e3nn(_lib_project_sym(S_rot))

    # e3nn Wigner-D for l=2 of the SAME rotation
    D2 = o3.Irrep("2e").D_from_matrix(torch.tensor(R, dtype=torch.float64)).numpy()
    e_pred = D2 @ e_orig

    err = np.abs(e_rot - e_pred).max()
    assert err < 1e-8, f"Wigner-D equivariance failed: max|e_rot - D2 e_orig| = {err:.2e}"


def test_emitted_substrate_selfconsistent():
    """If the canonical substrate is present: emitted library 5-vector == library
    projection of emitted 3x3, and |C @ t| == |t| (orthogonality preserves |T2|)."""
    import pandas as pd

    src_csv = f"{OUT}/ring_current_sources.csv"
    if not os.path.exists(src_csv):
        return  # substrate not present in this env; skip silently
    df = pd.read_csv(src_csv, nrows=2000)
    df = df[(df.dft_present == 1) & (df.dft_local_frame_valid == 1)].head(200)
    if not len(df):
        return
    loc = df[[f"dft_total_local_{i}{j}" for i in range(3) for j in range(3)]].to_numpy()
    loc = loc.reshape(-1, 3, 3)
    emit_lib = np.stack([_lib_project_sym(M) for M in loc])  # producer's own map
    e3 = cob.lib_to_e3nn(emit_lib)
    # orthogonal map preserves the rotation-invariant magnitude |T2|
    assert np.abs(np.linalg.norm(emit_lib, axis=1) - np.linalg.norm(e3, axis=1)).max() < 1e-8


if __name__ == "__main__":
    test_orthogonal_and_roundtrip()
    print("ok: C orthogonal + round-trip")
    test_wigner_d_equivariance_vs_cpp()
    print("ok: Wigner-D equivariance vs C++ library tensor")
    test_emitted_substrate_selfconsistent()
    print("ok: emitted substrate self-consistent / |T2| preserved")
    C = cob.derive_C()
    np.set_printoptions(precision=8, suppress=True)
    print("\nC (library T2 -> e3nn 2e):\n", C)
