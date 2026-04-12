#!/usr/bin/env python3
"""
Test whether AIMNet2 TorchScript model supports autograd for
d(charges)/d(positions) via backpropagation.

Goal: Replace 10 finite-difference forward passes with 1 backward pass
for charge sensitivity (the Jacobian of charges w.r.t. coordinates).

Key insight from reading calculator.py:
  - AIMNet2Calculator already uses torch.autograd.grad for forces (d(energy)/d(coord))
  - And even computes the full Hessian via autograd
  - So the model graph IS differentiable w.r.t. coord
  - But forces differentiate ENERGY, not CHARGES. We need to verify charges
    are also differentiable w.r.t. coord.
"""

import torch
import time
import sys
import os

sys.path.insert(0, '/tmp/aimnet2_repo')

MODEL_PATH = '/tmp/aimnet2_repo/aimnet2calc/assets/aimnet2/aimnet2_wb97m_0.jpt'
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"Device: {DEVICE}")
print(f"PyTorch: {torch.__version__}")

# ============================================================
# Part 1: Load model and build a small test molecule
# ============================================================
print("\n=== Loading model ===")
model = torch.jit.load(MODEL_PATH, map_location=DEVICE)
print(f"Model cutoff: {model.cutoff}")
print(f"Model cutoff_lr: {getattr(model, 'cutoff_lr', 'N/A')}")

N = 10  # real atoms
numbers_list = [6, 6, 6, 7, 8, 8, 1, 1, 1, 1]  # C3 N O2 H4

torch.manual_seed(42)
coord_np = torch.randn(N, 3, device=DEVICE) * 2.0

coord = torch.zeros(N + 1, 3, device=DEVICE, dtype=torch.float32)
coord[:N] = coord_np.float()

numbers = torch.zeros(N + 1, device=DEVICE, dtype=torch.long)
numbers[:N] = torch.tensor(numbers_list, device=DEVICE, dtype=torch.long)

charge = torch.tensor([0.0], device=DEVICE, dtype=torch.float32)
mol_idx = torch.zeros(N + 1, device=DEVICE, dtype=torch.long)

max_nb = N
nbmat = torch.full((N + 1, max_nb), N, device=DEVICE, dtype=torch.int32)
for i in range(N):
    nb_idx = 0
    for j in range(N):
        if i != j:
            nbmat[i, nb_idx] = j
            nb_idx += 1

max_nb_lr = N
nbmat_lr = nbmat.clone()
cutoff_lr = torch.tensor([15.0], device=DEVICE, dtype=torch.float32)

print(f"Test molecule: {N} atoms, types = {numbers_list}")
print(f"Coord shape: {coord.shape}, nbmat shape: {nbmat.shape}")

# ============================================================
# Part 2: Test autograd through charges
# ============================================================
print("\n=== Test 1: backward() on sum of charges ===")

coord_grad = coord.clone().detach().requires_grad_(True)
input_dict = {
    'coord': coord_grad,
    'numbers': numbers,
    'charge': charge,
    'mol_idx': mol_idx,
    'nbmat': nbmat,
    'nbmat_lr': nbmat_lr,
    'cutoff_lr': cutoff_lr,
}

try:
    with torch.jit.optimized_execution(False):
        output = model(input_dict)

    charges = output['charges']
    print(f"Output charges shape: {charges.shape}")
    print(f"Output charges: {charges}")
    print(f"charges.requires_grad: {charges.requires_grad}")

    if not charges.requires_grad:
        print("\nWARNING: charges tensor does not require grad!")
        print("This means the charge computation path may not be connected to coord")
        print("in the autograd graph.")

        if 'energy' in output:
            energy = output['energy']
            print(f"\nenergy.requires_grad: {energy.requires_grad}")
            if energy.requires_grad:
                grad_e = torch.autograd.grad(energy.sum(), coord_grad, retain_graph=True)
                print(f"energy grad shape: {grad_e[0].shape}")
                print("Energy IS differentiable (expected from calculator.py)")
                print("\nThis means charges are computed in a non-differentiable branch")
                print("(e.g., detach() or no_grad() inside the model for charge prediction)")
    else:
        charge_sum = charges[:N].sum()
        charge_sum.backward(retain_graph=True)

        if coord_grad.grad is not None:
            print(f"\nSUCCESS: coord.grad populated!")
            print(f"grad shape: {coord_grad.grad.shape}")
            print(f"grad norm: {coord_grad.grad.norm().item():.6e}")
            print(f"grad (first 3 atoms):\n{coord_grad.grad[:3]}")
            print(f"grad (padding atom): {coord_grad.grad[N]}")
        else:
            print("\nFAILED: coord.grad is None after backward()")

except Exception as e:
    import traceback
    print(f"\nERROR in backward test: {type(e).__name__}: {e}")
    traceback.print_exc()

# ============================================================
# Part 3: Test full Jacobian d(charges)/d(coord) via autograd
# ============================================================
print("\n=== Test 2: Full Jacobian d(charges)/d(coord) ===")

coord_grad2 = coord.clone().detach().requires_grad_(True)
input_dict2 = {
    'coord': coord_grad2,
    'numbers': numbers.clone(),
    'charge': charge.clone(),
    'mol_idx': mol_idx.clone(),
    'nbmat': nbmat.clone(),
    'nbmat_lr': nbmat_lr.clone(),
    'cutoff_lr': cutoff_lr.clone(),
}

try:
    with torch.jit.optimized_execution(False):
        output2 = model(input_dict2)

    charges2 = output2['charges'][:N]

    if charges2.requires_grad:
        jacobian_rows = []
        for i in range(N):
            grad_i = torch.autograd.grad(charges2[i], coord_grad2, retain_graph=True)[0]
            jacobian_rows.append(grad_i[:N])

        jacobian = torch.stack(jacobian_rows)  # (N, N, 3)
        print(f"Jacobian shape: {jacobian.shape}")
        print(f"Jacobian norm: {jacobian.norm().item():.6e}")
        print(f"Jacobian[0,0,:] (self-sensitivity of atom 0): {jacobian[0,0,:]}")
        print(f"Jacobian[0,1,:] (cross-sensitivity atom 0 from atom 1): {jacobian[0,1,:]}")

        # Check charge conservation: sum_i d(q_i)/d(r_j) should be ~0
        charge_conservation = jacobian.sum(dim=0)
        print(f"\nCharge conservation check (should be ~0):")
        print(f"  sum_i d(q_i)/d(r_j) norm: {charge_conservation.norm().item():.6e}")
        print(f"  max abs: {charge_conservation.abs().max().item():.6e}")

        # Diagonal vs off-diagonal magnitude
        diag_norms = torch.stack([jacobian[i, i, :] for i in range(N)]).norm()
        off_diag_norms = jacobian.norm() - diag_norms
        print(f"\nDiag (self) norm: {diag_norms.item():.6e}")
        print(f"Total norm: {jacobian.norm().item():.6e}")

        print("\nFull Jacobian computed successfully via autograd!")
    else:
        print("charges do not require grad -- trying torch.autograd.functional.jacobian")

        def charge_fn(c):
            d = {
                'coord': c,
                'numbers': numbers.clone(),
                'charge': charge.clone(),
                'mol_idx': mol_idx.clone(),
                'nbmat': nbmat.clone(),
                'nbmat_lr': nbmat_lr.clone(),
                'cutoff_lr': cutoff_lr.clone(),
            }
            with torch.jit.optimized_execution(False):
                out = model(d)
            return out['charges'][:N]

        try:
            J = torch.autograd.functional.jacobian(charge_fn, coord.clone().requires_grad_(True))
            print(f"Jacobian via functional API shape: {J.shape}")
            print(f"Jacobian norm: {J.norm().item():.6e}")
        except Exception as e2:
            print(f"functional.jacobian also failed: {type(e2).__name__}: {e2}")

except Exception as e:
    import traceback
    print(f"\nERROR in Jacobian test: {type(e).__name__}: {e}")
    traceback.print_exc()

# ============================================================
# Part 4: Timing comparison
# ============================================================
print("\n=== Test 3: Timing comparison ===")

def run_forward_only():
    c = coord.clone().detach()
    d = {
        'coord': c,
        'numbers': numbers.clone(),
        'charge': charge.clone(),
        'mol_idx': mol_idx.clone(),
        'nbmat': nbmat.clone(),
        'nbmat_lr': nbmat_lr.clone(),
        'cutoff_lr': cutoff_lr.clone(),
    }
    with torch.jit.optimized_execution(False), torch.no_grad():
        out = model(d)
    return out['charges'][:N]

def run_finite_diff(n_perturbations=10):
    results = []
    for _ in range(n_perturbations):
        c = coord.clone().detach()
        c[:N] += torch.randn_like(c[:N]) * 0.001
        d = {
            'coord': c,
            'numbers': numbers.clone(),
            'charge': charge.clone(),
            'mol_idx': mol_idx.clone(),
            'nbmat': nbmat.clone(),
            'nbmat_lr': nbmat_lr.clone(),
            'cutoff_lr': cutoff_lr.clone(),
        }
        with torch.jit.optimized_execution(False), torch.no_grad():
            out = model(d)
        results.append(out['charges'][:N])
    return results

def run_autograd_jacobian():
    c = coord.clone().detach().requires_grad_(True)
    d = {
        'coord': c,
        'numbers': numbers.clone(),
        'charge': charge.clone(),
        'mol_idx': mol_idx.clone(),
        'nbmat': nbmat.clone(),
        'nbmat_lr': nbmat_lr.clone(),
        'cutoff_lr': cutoff_lr.clone(),
    }
    with torch.jit.optimized_execution(False):
        out = model(d)
    charges_out = out['charges'][:N]

    if not charges_out.requires_grad:
        return None

    rows = []
    for i in range(N):
        g = torch.autograd.grad(charges_out[i], c, retain_graph=True)[0]
        rows.append(g[:N])
    return torch.stack(rows)

# Warmup
print("Warming up...")
for _ in range(3):
    run_forward_only()
if DEVICE == 'cuda':
    torch.cuda.synchronize()

N_REPEATS = 20

# Time forward only
if DEVICE == 'cuda':
    torch.cuda.synchronize()
t0 = time.perf_counter()
for _ in range(N_REPEATS):
    run_forward_only()
if DEVICE == 'cuda':
    torch.cuda.synchronize()
t_forward = (time.perf_counter() - t0) / N_REPEATS
print(f"Single forward pass: {t_forward*1000:.2f} ms")

# Time 10 finite-diff forward passes
if DEVICE == 'cuda':
    torch.cuda.synchronize()
t0 = time.perf_counter()
for _ in range(N_REPEATS):
    run_finite_diff(10)
if DEVICE == 'cuda':
    torch.cuda.synchronize()
t_fd = (time.perf_counter() - t0) / N_REPEATS
print(f"10 finite-diff forward passes: {t_fd*1000:.2f} ms")

# Time autograd Jacobian
try:
    # Warmup autograd
    for _ in range(3):
        run_autograd_jacobian()
    if DEVICE == 'cuda':
        torch.cuda.synchronize()

    t0 = time.perf_counter()
    for _ in range(N_REPEATS):
        J = run_autograd_jacobian()
    if DEVICE == 'cuda':
        torch.cuda.synchronize()
    t_autograd = (time.perf_counter() - t0) / N_REPEATS

    if J is not None:
        print(f"Autograd Jacobian ({N} backward passes): {t_autograd*1000:.2f} ms")
        print(f"\nSpeedup vs finite-diff: {t_fd/t_autograd:.2f}x")
        print(f"Autograd gives EXACT Jacobian; finite-diff is approximate.")
    else:
        print("Autograd path returned None (charges not differentiable)")
except Exception as e:
    print(f"Autograd timing failed: {type(e).__name__}: {e}")

# ============================================================
# Part 5: Output inspection
# ============================================================
print("\n=== Output inspection ===")
c = coord.clone().detach().requires_grad_(True)
d = {
    'coord': c,
    'numbers': numbers.clone(),
    'charge': charge.clone(),
    'mol_idx': mol_idx.clone(),
    'nbmat': nbmat.clone(),
    'nbmat_lr': nbmat_lr.clone(),
    'cutoff_lr': cutoff_lr.clone(),
}
with torch.jit.optimized_execution(False):
    out = model(d)

for k, v in out.items():
    if isinstance(v, torch.Tensor):
        print(f"  {k}: shape={v.shape}, dtype={v.dtype}, requires_grad={v.requires_grad}")

print("\n=== DONE ===")
