"""Pre-flight check: does data/models/aimnet2_wb97m_0.jpt support
requires_grad on the coordinate input tensor?

Gates whether the planned AIMNet2PolarisabilityResult slice (per
spec/PLANNED_CALCULATORS_2026-04-22.md Amendment 2026-05-08(b)) can
be implemented as a single autograd backward pass through the
TorchScript model, or whether the model needs re-exporting from the
original PyTorch with grad tracking enabled.

Inputs to the AIMNet2 model (per src/AIMNet2Result.cpp lines 200-275):

  coord     (N+1, 3) float32  — atom positions + sentinel padding row
  numbers   (N+1,)   int64    — atomic numbers + sentinel
  charge    (1,)     float32  — net molecular charge
  mol_idx   (N+1,)   int64    — molecule index per atom (all zeros for
                                 single molecule + sentinel zero)
  nbmat     (N+1, max_nb)   int32 — short-range half-neighbour list
  nbmat_lr  (N+1, max_nb_lr) int32 — long-range half-neighbour list
  cutoff_lr (1,)     float32  — long-range cutoff

This script builds a minimal 3-atom system (water — H, O, H), runs
forward with `coord.requires_grad_()`, calls backward on the sum of
predicted charges, and asserts coord.grad propagates.

Outcome decides the calculator slice plan:

  - PASSES → the AIMNet2PolarisabilityResult calculator is ~50 lines
    of C++: forward without NoGradGuard, single backward, copy
    gradients off CPU. Tractable next session.

  - FAILS → the .jpt was exported with grad-tracking disabled. Need
    to regenerate from the original PyTorch with explicit gradient
    flow enabled. Upstream model-asset task before any C++.
"""

from pathlib import Path
import sys

import torch


JPT_PATH = Path("/shared/2026Thesis/nmr-shielding/data/models/aimnet2_wb97m_0.jpt")


def build_minimal_input(device: str = "cpu"):
    """Build a 3-atom water input: H, O, H. Sentinel padding row at index 3."""
    N = 3
    N1 = N + 1

    # Water-like geometry, Angstroms.
    coord_data = torch.tensor(
        [
            [0.000,  0.760, 0.000],   # H1
            [0.000,  0.000, 0.000],   # O
            [0.000, -0.760, 0.000],   # H2
            [0.000,  0.000, 0.000],   # sentinel
        ],
        dtype=torch.float32,
        device=device,
    )

    numbers = torch.tensor([1, 8, 1, 0], dtype=torch.int64, device=device)
    charge = torch.zeros(1, dtype=torch.float32, device=device)
    mol_idx = torch.zeros(N1, dtype=torch.int64, device=device)

    # Neighbour matrices: each atom's neighbours within cutoff. For 3
    # atoms we hand-build pairs. AIMNet2 uses half-neighbour lists
    # (each pair listed once). max_nb = 4 is plenty.
    max_nb = 4
    # Initialise with sentinel index (N) — "no neighbour."
    nbmat = torch.full((N1, max_nb), N, dtype=torch.int32, device=device)
    nbmat_lr = torch.full((N1, max_nb), N, dtype=torch.int32, device=device)

    # Half-neighbour pairs for water (within ~1.0 A): O-H1, O-H2.
    # Atom 0 (H1) → atom 1 (O)
    nbmat[0, 0] = 1
    # Atom 1 (O) → atom 2 (H2)  (and atom 0 by symmetry of the half-list,
    # but AIMNet2's half-list convention only stores one direction)
    nbmat[1, 0] = 2

    # Long-range list — same content for this small test.
    nbmat_lr[0, 0] = 1
    nbmat_lr[1, 0] = 2

    cutoff_lr = torch.tensor([10.0], dtype=torch.float32, device=device)

    return coord_data, numbers, charge, mol_idx, nbmat, nbmat_lr, cutoff_lr


def main() -> int:
    if not JPT_PATH.exists():
        print(f"FAIL: model not found at {JPT_PATH}")
        return 2

    # Stick to CPU so we don't fight codex for the GPU.
    device = "cpu"
    print(f"Loading {JPT_PATH} on {device} ...")
    try:
        model = torch.jit.load(str(JPT_PATH), map_location=device)
    except Exception as exc:
        print(f"FAIL: torch.jit.load raised: {exc}")
        return 2
    model.eval()

    coord, numbers, charge, mol_idx, nbmat, nbmat_lr, cutoff_lr = (
        build_minimal_input(device=device)
    )

    # Enable grad on coords. This is the question the test asks: does
    # the TorchScript model propagate gradient through the coord input?
    coord = coord.detach().clone().requires_grad_(True)

    inputs = {
        "coord":     coord,
        "numbers":   numbers,
        "charge":    charge,
        "mol_idx":   mol_idx,
        "nbmat":     nbmat,
        "nbmat_lr":  nbmat_lr,
        "cutoff_lr": cutoff_lr,
    }

    print("Running forward pass ...")
    try:
        output = model(inputs)
    except Exception as exc:
        print(f"FAIL: forward pass raised: {exc}")
        return 3

    if "charges" not in output:
        print(f"FAIL: output dict missing 'charges' key. Keys: {list(output.keys())}")
        return 3

    charges = output["charges"]
    print(f"Charges (first 3): {charges[:3].detach().tolist()}")
    print(f"Charges grad_fn:   {charges.grad_fn}")

    # Test 1: charges should have a grad_fn (i.e., the model didn't
    # detach them internally).
    if charges.grad_fn is None:
        print("FAIL: charges tensor has no grad_fn — the .jpt model "
              "appears to have been exported with grad tracking disabled "
              "(NoGradGuard, .detach(), or torch.no_grad context inside "
              "the scripted forward). AIMNet2PolarisabilityResult cannot "
              "be implemented from this .jpt; need to re-export the "
              "model from PyTorch with grad tracking enabled.")
        return 1

    # Test 2: backward pass should populate coord.grad.
    print("Running backward pass on sum(charges) ...")
    loss = charges[:3].sum()  # sum non-sentinel atoms only
    loss.backward()

    if coord.grad is None:
        print("FAIL: coord.grad is None after backward. Gradient flow "
              "broke between coord input and charges output.")
        return 1

    grad_norm = coord.grad.norm().item()
    print(f"coord.grad.norm() = {grad_norm:.6e}")

    if grad_norm < 1e-12:
        print("FAIL: coord.grad has zero norm. Gradient flowed but "
              "yielded zero everywhere, which is wrong for a model "
              "where charges depend on coordinates.")
        return 1

    print("")
    print("PASS: requires_grad propagates through the .jpt model. The "
          "AIMNet2PolarisabilityResult calculator slice is unblocked.")
    print("Per-atom polarisability is the gradient of charges with "
          "respect to coordinates, which this run confirmed is "
          "computable in a single backward pass.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
