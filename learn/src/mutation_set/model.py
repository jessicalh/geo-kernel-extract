"""Equivariant T2 model: gated kernel mixing + optional correction head.

The MLP learns per-atom relative importance of each kernel.  Each kernel
gates its own weight by its T2 magnitude — a kernel that reports zero
at an atom contributes nothing.  Distance dependence emerges from the
kernels' own physics (1/r^3 falloff, cutoffs) rather than from a
training hyperparameter.

Equivariance: scalar_function(scalars) x L=2_tensor = L=2_tensor.

Distance dependence emerges from the kernels' own physics: 1/r^3 dipolar
(BS, HM, McConnell, ring susceptibility, H-bond), 1/r^5 quadrupolar
(pi-quadrupole), and short-range switching (dispersion).
"""

from __future__ import annotations

import torch
import torch.nn as nn

try:
    import e3nn
    from e3nn import o3
    from e3nn.nn import FullyConnectedNet
    e3nn.set_optimization_defaults(jit_script_fx=False)
except AttributeError:
    pass

from .config import Config


class KernelMixingHead(nn.Module):
    """Scalar MLP → per-kernel weights, gated by kernel magnitude."""

    def __init__(self, n_scalars: int, n_kernels: int, hidden: int):
        super().__init__()
        self.n_kernels = n_kernels
        self.mlp = FullyConnectedNet(
            [n_scalars, hidden, hidden, n_kernels],
            act=torch.nn.functional.silu,
        )
        self.register_buffer("gate_threshold", torch.ones(n_kernels))

    def set_gate_thresholds(self, thresholds: torch.Tensor):
        self.gate_threshold.copy_(thresholds)

    def forward(self, scalars, kernel_t2s):
        import math
        relative = self.mlp(scalars)
        magnitude = kernel_t2s.norm(dim=-1)
        gate = magnitude / (magnitude + self.gate_threshold)
        # Gate sparsifies inactive kernels; sqrt(n_kernels) normalizes output
        # scale so initial predictions match target std. Both are needed.
        weights = relative * gate / math.sqrt(self.n_kernels)
        return torch.einsum("bk,bkm->bm", weights, kernel_t2s)


class EquivariantCorrectionHead(nn.Module):
    """L=2 correction via e3nn tensor products of selected kernels."""

    def __init__(self, n_scalars: int, l2_indices: list[int], hidden: int):
        super().__init__()
        self.l2_indices = l2_indices
        n_l2 = len(l2_indices) + 1  # +1 for kernel_sum

        self.tp1 = o3.FullyConnectedTensorProduct(
            o3.Irreps(f"{n_scalars}x0e + {n_l2}x2e"),
            o3.Irreps(f"{n_scalars}x0e + {n_l2}x2e"),
            o3.Irreps(f"{hidden}x0e + {hidden}x2e"),
        )
        self.tp2 = o3.FullyConnectedTensorProduct(
            o3.Irreps(f"{hidden}x0e + {hidden}x2e"),
            o3.Irreps(f"{hidden}x0e + {hidden}x2e"),
            o3.Irreps("1x2e"),
        )

    def forward(self, scalars, kernel_t2s):
        parts = [kernel_t2s[:, i, :] for i in self.l2_indices]
        parts.append(kernel_t2s.sum(dim=1))
        x = torch.cat([scalars] + parts, dim=-1)
        return self.tp2(self.tp1(x, x), self.tp1(x, x))


class ShieldingT2Model(nn.Module):
    """Gated kernel mixing + optional equivariant correction."""

    def __init__(self, cfg: Config, n_scalars: int, n_kernels: int):
        super().__init__()
        self.mixing = KernelMixingHead(
            n_scalars, n_kernels, cfg.model.hidden_mix)
        self.use_correction = cfg.model.use_correction
        if self.use_correction:
            self.correction = EquivariantCorrectionHead(
                n_scalars, cfg.model.correction_l2_kernels,
                cfg.model.hidden_corr)
            self.correction_scale = nn.Parameter(
                torch.tensor(cfg.model.correction_scale_init))

    def forward(self, scalars, kernel_t2s):
        mixed = self.mixing(scalars, kernel_t2s)
        if self.use_correction:
            corr = self.correction(scalars, kernel_t2s)
            return mixed + self.correction_scale * corr
        return mixed

    def parameter_count(self):
        return sum(p.numel() for p in self.parameters())
