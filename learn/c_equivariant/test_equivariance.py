#!/usr/bin/env python3
"""
Equivariance rotation test for ShieldingT2Model.

Verifies that rotating L=2 inputs produces the same rotation in outputs.
A wrong e3nn irreps string would be a silent numerical bug — this catches it.

Usage:
    python learn/c_equivariant/test_equivariance.py
"""

import sys
from pathlib import Path

import torch
import numpy as np
from e3nn import o3

sys.path.insert(0, str(Path(__file__).parent.parent))
from c_equivariant.model import make_model, KernelMixingHead, EquivariantCorrectionHead
from c_equivariant.dataset import N_KERNELS, N_SCALAR_FEATURES


def rotate_t2(t2: torch.Tensor, D2: torch.Tensor) -> torch.Tensor:
    """Apply Wigner-D rotation to L=2 vectors.

    t2: (..., 5) — L=2 components in m=-2..+2 order
    D2: (5, 5)  — Wigner D-matrix for L=2
    returns: (..., 5)
    """
    return t2 @ D2.T


def test_mixing_head():
    """KernelMixingHead: scalar(MLP) × L=2 = L=2.  Trivially equivariant."""
    torch.manual_seed(42)
    batch = 16
    head = KernelMixingHead(N_SCALAR_FEATURES, N_KERNELS)
    head.eval()

    scalars = torch.randn(batch, N_SCALAR_FEATURES)
    kernels = torch.randn(batch, N_KERNELS, 5)

    # Random rotation
    angles = torch.randn(3)
    D2 = o3.wigner_D(2, *angles)  # (5, 5)

    # Forward, then rotate output
    out_orig = head(scalars, kernels)
    out_rotated = rotate_t2(out_orig, D2)

    # Rotate inputs, then forward
    kernels_rot = torch.stack(
        [rotate_t2(kernels[:, k, :], D2) for k in range(N_KERNELS)], dim=1)
    out_from_rot = head(scalars, kernels_rot)

    err = (out_rotated - out_from_rot).abs().max().item()
    print(f"KernelMixingHead equivariance error: {err:.2e}")
    assert err < 1e-5, f"FAILED: error {err} too large"
    print("  PASSED")


def test_correction_head():
    """EquivariantCorrectionHead: e3nn tensor products with mixed irreps."""
    torch.manual_seed(42)
    batch = 16
    head = EquivariantCorrectionHead(N_SCALAR_FEATURES, N_KERNELS)
    head.eval()

    scalars = torch.randn(batch, N_SCALAR_FEATURES)
    kernels = torch.randn(batch, N_KERNELS, 5)

    # Random rotation
    angles = torch.randn(3)
    D2 = o3.wigner_D(2, *angles)

    # Forward, then rotate output
    out_orig = head(scalars, kernels)
    out_rotated = rotate_t2(out_orig, D2)

    # Rotate ALL L=2 kernel inputs, then forward
    kernels_rot = torch.stack(
        [rotate_t2(kernels[:, k, :], D2) for k in range(N_KERNELS)], dim=1)
    out_from_rot = head(scalars, kernels_rot)

    err = (out_rotated - out_from_rot).abs().max().item()
    print(f"EquivariantCorrectionHead equivariance error: {err:.2e}")
    assert err < 1e-4, f"FAILED: error {err} too large"
    print("  PASSED")


def test_full_model():
    """ShieldingT2Model: mixing + correction combined."""
    torch.manual_seed(42)
    batch = 16
    model = make_model(use_correction=True)
    model.eval()

    scalars = torch.randn(batch, N_SCALAR_FEATURES)
    kernels = torch.randn(batch, N_KERNELS, 5)

    angles = torch.randn(3)
    D2 = o3.wigner_D(2, *angles)

    out_orig = model(scalars, kernels)
    out_rotated = rotate_t2(out_orig, D2)

    kernels_rot = torch.stack(
        [rotate_t2(kernels[:, k, :], D2) for k in range(N_KERNELS)], dim=1)
    out_from_rot = model(scalars, kernels_rot)

    err = (out_rotated - out_from_rot).abs().max().item()
    print(f"Full model equivariance error: {err:.2e}")
    assert err < 1e-4, f"FAILED: error {err} too large"
    print("  PASSED")


def test_multiple_rotations():
    """Test equivariance across 10 random rotations for robustness."""
    torch.manual_seed(123)
    batch = 8
    model = make_model(use_correction=True)
    model.eval()

    scalars = torch.randn(batch, N_SCALAR_FEATURES)
    kernels = torch.randn(batch, N_KERNELS, 5)

    max_err = 0.0
    for _ in range(10):
        angles = torch.randn(3)
        D2 = o3.wigner_D(2, *angles)

        out_orig = model(scalars, kernels)
        out_rotated = rotate_t2(out_orig, D2)

        kernels_rot = torch.stack(
            [rotate_t2(kernels[:, k, :], D2) for k in range(N_KERNELS)], dim=1)
        out_from_rot = model(scalars, kernels_rot)

        err = (out_rotated - out_from_rot).abs().max().item()
        max_err = max(max_err, err)

    print(f"10-rotation max equivariance error: {max_err:.2e}")
    assert max_err < 1e-4, f"FAILED: max error {max_err} too large"
    print("  PASSED")


if __name__ == "__main__":
    print("=" * 60)
    print("Equivariance rotation tests")
    print("=" * 60)
    test_mixing_head()
    test_correction_head()
    test_full_model()
    test_multiple_rotations()
    print("\nAll equivariance tests passed.")
