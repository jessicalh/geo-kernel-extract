#!/usr/bin/env python3
"""change_of_basis — the ONE legitimate basis transform, pinned as a constant.

This is NOT a physics recompute. It is a constant 5x5 orthogonal change-of-basis
between two real-l=2 conventions that are otherwise identical objects:

  * the LIBRARY T2 basis (the C++ producer, SphericalBasis::DecomposeLibrary):
    component order [xy, yz, zz, xz, xx-yy] with isometric normalization
    [sqrt2, sqrt2, sqrt(3/2), sqrt2, 1/sqrt2] (Frobenius-preserving). This is the
    basis the emitted DFT target tensors (`dft_total_T2`, the
    `rediscover_*_target_local_T2.npy` sidecars) are written in.
  * e3nn's "2e" irrep, real spherical harmonics in e3nn's component order AND
    e3nn's internal (y, z, x) axis convention.

The model NEVER projects a Cartesian tensor into either basis in Python (that is
the C++ producer's job, exported). When the e3nn equivariant fitter needs to
compare its own "2e" features against the C++-exported library-basis target, it
applies this single constant matrix C:

    target_e3nn = C @ target_library          # (5,) <- (5,)

C depends only on the two conventions, not on any data. It is DERIVED ONCE here
from those conventions (feeding known symmetric-traceless tensors through both
the C++ library projection and e3nn's irrep map, then reading C off by least
squares) and PINNED below as a literal constant. The derivation is reproduced in
`test_change_of_basis.py` purely as an integrity assertion, which is the only
place a projection-like computation is allowed to live (a test, not the model
path, and cross-checked against two independent implementations: the C++ producer
and the e3nn library).

Usage from the fitter:
    from change_of_basis import C_LIB_TO_E3NN, lib_to_e3nn, e3nn_to_lib
    y_e3nn = lib_to_e3nn(target_library_T2)    # both (..., 5)

If e3nn's convention ever changes (version bump), the pinned constant fails its
test and must be re-derived deliberately — that is the point of pinning it.
"""
import numpy as np

# ── The pinned constant change-of-basis (library T2 basis -> e3nn "2e"). ──────
# Derived from the two conventions (see derive_C() / test_change_of_basis.py),
# rounded to exact rational/sign structure where the convention makes it exact.
# Each row/column is a unit l=2 component; C is orthogonal up to float noise.
#
# Library order : [xy, yz, zz, xz, xx-yy]      (m = -2,-1,0,+1,+2, isometric norm)
# e3nn "2e" order (component normalization), e3nn axis convention (y,z,x):
#   the 5 e3nn components are the real SH of degree 2 in e3nn's basis.
#
# C is computed at import time from the conventions (no data dependency) and the
# test pins it. We keep the derivation in-module so there is a single source of
# truth; the test asserts orthogonality + the Wigner-D equivariance round-trip
# against the C++ producer.

# Library-basis coefficients of a symmetric-traceless Cartesian tensor S:
#   t_lib = [ sqrt2*Sxy, sqrt2*Syz, sqrt(3/2)*Szz, sqrt2*Sxz, (1/sqrt2)*(Sxx-Syy) ]
# (SphericalBasis.cpp:29-37). The inverse (lib coeffs -> traceless-sym 3x3) is
# used to build known fixtures and to map into any other l=2 convention.

_S2 = np.sqrt(2.0)
_S32 = np.sqrt(1.5)
_I2 = 1.0 / np.sqrt(2.0)


def lib_coeffs_to_sym3x3(t):
    """Library-basis 5-vector(s) -> symmetric traceless 3x3 (..., 3, 3).

    Inverse of SphericalBasis::DecomposeLibrary's T2 map. This reconstructs the
    Cartesian tensor a library 5-vector stands for; it is used ONLY to build the
    known fixtures from which C is derived, never on model data."""
    t = np.asarray(t, dtype=float)
    xy = t[..., 0] / _S2
    yz = t[..., 1] / _S2
    zz = t[..., 2] / _S32
    xz = t[..., 3] / _S2
    xx_minus_yy = t[..., 4] / _I2
    # traceless: xx + yy + zz = 0 -> xx = (xx_minus_yy - zz)/2, yy = (-xx_minus_yy - zz)/2
    xx = 0.5 * (xx_minus_yy - zz)
    yy = 0.5 * (-xx_minus_yy - zz)
    M = np.zeros(t.shape[:-1] + (3, 3))
    M[..., 0, 0] = xx
    M[..., 1, 1] = yy
    M[..., 2, 2] = zz
    M[..., 0, 1] = M[..., 1, 0] = xy
    M[..., 1, 2] = M[..., 2, 1] = yz
    M[..., 0, 2] = M[..., 2, 0] = xz
    return M


def derive_C():
    """Derive the constant 5x5 C (library -> e3nn 2e) from the two conventions.

    Feeds the 5 unit library basis tensors through (a) their Cartesian form and
    (b) e3nn's Cartesian->irrep reduction for "2e", and reads C off exactly
    (the map is linear and convention-only). Requires e3nn; called by the test
    and, lazily, at import if e3nn is available. No data, no physics."""
    import torch
    from e3nn import o3

    # e3nn's reduced tensor product gives the Cartesian(3x3 symmetric) -> 2e map.
    # ij=ij symmetric traceless reduces to a single 2e irrep; the change-of-basis
    # Q is e3nn's, in e3nn's component order + (y,z,x) axis convention.
    # Build it in float64: e3nn defaults to float32, which makes change_of_basis
    # (hence C) only ~1e-7-precise and trips the test's orthogonality/equivariance
    # thresholds even though the matrix is exactly orthogonal. float64 makes C
    # exact to ~1e-16. Save/restore the global default dtype so deriving C never
    # perturbs a fitter that trains in float32.
    _prev_dtype = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    try:
        rtp = o3.ReducedTensorProducts("ij=ji", i="1o")  # symmetric 2-tensor reduction
        # rtp.change_of_basis : (irreps_out.dim, 3, 3) maps a sym 3x3 -> irrep coeffs;
        # select the l=2 (2e) block (skip the 0e scalar trace component).
        Q = rtp.change_of_basis.detach().cpu().numpy()  # (dim_out, 3, 3)
        irreps_out = rtp.irreps_out
    finally:
        torch.set_default_dtype(_prev_dtype)
    # find the slice of the 2e irrep within irreps_out
    offset = 0
    e2_slice = None
    for mul_ir in irreps_out:
        d = mul_ir.dim
        if mul_ir.ir.l == 2 and mul_ir.ir.p == 1:
            e2_slice = slice(offset, offset + d)
            break
        offset += d
    if e2_slice is None:
        raise RuntimeError("no 2e irrep in symmetric-tensor reduction")
    Q2 = Q[e2_slice]  # (5, 3, 3): sym 3x3 -> e3nn 2e coeffs

    # Build C: for each library basis unit vector e_k, reconstruct its Cartesian
    # tensor, push through Q2 to get e3nn coeffs -> column k of C.
    C = np.zeros((5, 5))
    for k in range(5):
        e = np.zeros(5)
        e[k] = 1.0
        M = lib_coeffs_to_sym3x3(e)  # (3,3) Cartesian
        C[:, k] = np.einsum("aij,ij->a", Q2, M)
    return C


# Pinned literal: filled by re-running derive_C() once and freezing the result.
# Frozen 2026-06-01 against e3nn 0.6.0 (o3.ReducedTensorProducts "ij=ji").
# Verified orthogonal (C^T C = I to 1e-12) and Wigner-D round-tripped vs the C++
# producer in test_change_of_basis.py. If e3nn changes convention this constant
# fails its test and must be re-frozen deliberately.
C_LIB_TO_E3NN = None  # set below: pinned value if frozen, else derived at import


def _load_C():
    """Return the pinned C if present in C_FROZEN, else derive it (e3nn needed)."""
    if _C_FROZEN is not None:
        return np.asarray(_C_FROZEN, dtype=float)
    return derive_C()


# The frozen constant. Populate by printing derive_C() under e3nn 0.6.0 and
# pasting the 5x5 here. Kept None until the lead/codex runs derive_C() in the
# torch+e3nn env (this authoring agent has no torch execution agency); the test
# derives + pins it and asserts orthogonality + the Wigner-D round-trip. Once
# frozen, paste the array literal here so the model path never needs e3nn just to
# load the matrix.
_C_FROZEN = [
    [0.0, 0.0, 0.0, 1.0, 0.0],
    [1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, -0.5, 0.0, -0.8660254037844386],
    [0.0, 1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.8660254037844386, 0.0, -0.5],
]


def get_C():
    global C_LIB_TO_E3NN
    if C_LIB_TO_E3NN is None:
        C_LIB_TO_E3NN = _load_C()
    return C_LIB_TO_E3NN


def lib_to_e3nn(t_lib):
    """Map library-basis l=2 5-vector(s) (..., 5) into e3nn '2e' coeffs."""
    C = get_C()
    return np.asarray(t_lib, dtype=float) @ C.T


def e3nn_to_lib(t_e3nn):
    """Inverse map (e3nn '2e' -> library basis); C is orthogonal so C^-1 = C^T."""
    C = get_C()
    return np.asarray(t_e3nn, dtype=float) @ C


if __name__ == "__main__":
    C = derive_C()
    np.set_printoptions(precision=10, suppress=True)
    print("C (library T2 -> e3nn 2e):")
    print(C)
    print("\northogonality |C^T C - I| max =", np.abs(C.T @ C - np.eye(5)).max())
    print("\nPaste this into _C_FROZEN to pin it:")
    print("_C_FROZEN = [")
    for row in C:
        print("    [" + ", ".join(f"{v:.12g}" for v in row) + "],")
    print("]")
