"""Sanity check: AIMNet2 architecturally enforces charge conservation
(Σ q_j = q_input).

Motivation (science review S7, 2026-05-20): the AIMNet2ChargeResponseGradientResult
calculator differentiates L = Σ q_j² with respect to atomic coordinates,
producing the charge-response gradient channel. Note d(Σq_j²)/dr_i =
2 Σ_j q_j · J_ji where J = ∂q/∂r. If AIMNet2 enforces charge
conservation EXACTLY (architectural projection on the charge head),
then `∂(Σ q_j)/∂r_i ≡ 0` for all i, and the per-atom gradient of
Σ q² reduces to the centred charge-weighted sum
    d(Σq_j²)/dr_i = 2 Σ_j (q_j − q̄)·J_ji  (since Σ_j J_ji = 0).
If AIMNet2 enforces conservation only APPROXIMATELY (soft loss during
training), there is a residual all-ones bias term in the gradient
that should either be projected out at backward time or flagged in
the docstring.

This script:
  1. Loads `data/models/aimnet2_wb97m_0.jpt`.
  2. Builds a minimal water input (same as
     aimnet2_requires_grad_check.py for direct comparison).
  3. Runs forward with coord.requires_grad_().
  4. Computes L_sum = Σ_j q_j (NOT L = Σ q_j²) and calls backward.
  5. Asserts the resulting coord.grad has |grad|_max < 1e-5.

Outcome:
  - PASSES → AIMNet2 architecturally enforces charge conservation;
    the charge-response gradient channel has no spurious all-ones bias.
  - FAILS → conservation is approximate. Either (a) the charge-response
    gradient needs explicit projection at backward time, or (b) the
    AIMNet2ChargeResponseGradientResult docstring needs to flag the residual
    bias. Calibration-ridge would absorb the bias, but flagging it
    keeps the audit trail honest.
"""

from pathlib import Path
import os
import sys

import torch


REPO_ROOT = Path(__file__).resolve().parents[3]
JPT_PATH = Path(os.environ.get(
    "NMR_AIMNET2_MODEL",
    REPO_ROOT / "data/models/aimnet2_wb97m_0.jpt",
))
TOL_MAX_ABS_GRAD = 1e-5


def build_minimal_input(device: str = "cpu"):
    """Build a 3-atom water input: H, O, H. Sentinel padding row at index 3.
    Identical to aimnet2_requires_grad_check.py build."""
    N = 3
    N1 = N + 1

    coord_data = torch.tensor(
        [
            [0.000,  0.760, 0.000],
            [0.000,  0.000, 0.000],
            [0.000, -0.760, 0.000],
            [0.000,  0.000, 0.000],
        ],
        dtype=torch.float32,
        device=device,
    )

    numbers = torch.tensor([1, 8, 1, 0], dtype=torch.int64, device=device)
    charge = torch.zeros(1, dtype=torch.float32, device=device)
    mol_idx = torch.zeros(N1, dtype=torch.int64, device=device)

    max_nb = 4
    nbmat = torch.full((N1, max_nb), N, dtype=torch.int32, device=device)
    nbmat_lr = torch.full((N1, max_nb), N, dtype=torch.int32, device=device)
    nbmat[0, 0] = 1
    nbmat[1, 0] = 2

    cutoff_lr = torch.tensor([15.0], dtype=torch.float32, device=device)

    return coord_data, numbers, charge, mol_idx, nbmat, nbmat_lr, cutoff_lr, N


def main():
    if not JPT_PATH.exists():
        print(f"SKIP: AIMNet2 .jpt not found at {JPT_PATH}")
        sys.exit(0)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = torch.jit.load(str(JPT_PATH), map_location=device).eval()

    coord, numbers, charge, mol_idx, nbmat, nbmat_lr, cutoff_lr, N = (
        build_minimal_input(device=device)
    )
    coord = coord.detach().clone().requires_grad_(True)

    output_dict = model(
        coord=coord, numbers=numbers, charge=charge,
        mol_idx=mol_idx, nbmat=nbmat, nbmat_lr=nbmat_lr,
        cutoff_lr=cutoff_lr,
    )
    if "charges" not in output_dict:
        print("FAIL: model output missing 'charges' key")
        sys.exit(1)

    q = output_dict["charges"].to(torch.float64).squeeze()
    # Drop sentinel row N (per AIMNet2 padding convention).
    q_real = q[:N]
    L_sum = q_real.sum()

    # Backward Σ q → ∂(Σ q)/∂r.
    coord.grad = None
    L_sum.backward()
    grad_max = coord.grad.abs().max().item()

    print(f"  AIMNet2 charge-conservation check ({device}):")
    print(f"    Σ q (3-atom water) = {L_sum.item():.6e}")
    print(f"    |∂(Σq)/∂r|_max     = {grad_max:.6e}")
    print(f"    tolerance          = {TOL_MAX_ABS_GRAD:.1e}")

    if grad_max < TOL_MAX_ABS_GRAD:
        print("PASS — AIMNet2 architecturally enforces charge conservation.")
        print("       ∂(Σq²)/∂r in AIMNet2ChargeResponseGradientResult is the centred")
        print("       charge-weighted Jacobian sum; no spurious all-ones bias.")
        sys.exit(0)
    else:
        print("FAIL — ∂(Σq)/∂r exceeds tolerance; conservation is approximate.")
        print("       AIMNet2ChargeResponseGradientResult's gradient has a residual")
        print("       all-ones bias of magnitude ~{:.2e}. Calibration-ridge".format(grad_max))
        print("       will absorb it but the charge-response gradient should be")
        print("       projected to remove the bias, OR the docstring should flag")
        print("       this so downstream readers don't treat the gradient as the")
        print("       physically-pure centred-charge sum.")
        sys.exit(1)


if __name__ == "__main__":
    main()
