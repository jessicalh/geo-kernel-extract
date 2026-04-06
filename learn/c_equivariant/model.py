"""
Level C: Equivariant T2 model.

Architecture:
    Input per atom:
        - Scalar features (L=0): ring counts, element, residue type,
          distances to nearest ring/bond/hbond
        - Calculator T2 kernels (L=2): one 5-vector per kernel type

    Output per atom:
        - Predicted DFT T2 (L=2, 5 components)

    The model learns position-dependent mixing weights for the calculator
    kernels, plus a correction term from the equivariant layers.

    Equivariance: e3nn tensor products guarantee that rotating the
    input rotates the output correctly. The L=2 kernels from the
    calculators already encode the angular physics.
"""

import torch
import torch.nn as nn
import e3nn
from e3nn import o3
from e3nn.nn import FullyConnectedNet

# Blackwell (sm_121) workaround: NVRTC JIT doesn't support sm_121.
# Eager mode is slightly slower but works.
try:
    e3nn.set_optimization_defaults(jit_script_fx=False)
except AttributeError:
    pass


class KernelMixingHead(nn.Module):
    """Learn scalar mixing weights from scalar features, apply to L=2 kernels.

    This is the equivariant generalisation of Level A: instead of global
    scalar weights, the weights depend on local environment (ring proximity,
    element type, etc.) through a small MLP on scalar features.

    Equivariant because: scalar_function(scalars) × L=2_tensor = L=2_tensor.
    """

    def __init__(self, n_scalar_features: int, n_kernels: int, hidden: int = 64):
        super().__init__()
        self.n_kernels = n_kernels
        # MLP: scalar features → one weight per kernel
        self.mlp = FullyConnectedNet(
            [n_scalar_features, hidden, hidden, n_kernels],
            act=torch.nn.functional.silu,
        )

    def forward(self, scalars, kernel_t2s):
        """
        scalars: (batch, n_scalar_features) — L=0 input
        kernel_t2s: (batch, n_kernels, 5) — L=2 input from each calculator
        returns: (batch, 5) — weighted sum of kernel T2s
        """
        weights = self.mlp(scalars)  # (batch, n_kernels)
        # Scale weights so initial output std ≈ 1 (matches target std after normalization).
        # Without this, summing n_kernels terms with std~1 gives output std~sqrt(n_kernels).
        import math
        weights = weights / math.sqrt(self.n_kernels)
        # Weighted sum: each weight is a scalar multiplying an L=2 vector
        out = torch.einsum("bk,bkm->bm", weights, kernel_t2s)
        return out


class EquivariantCorrectionHead(nn.Module):
    """Learn an L=2 correction from scalar and L=2 inputs via tensor products.

    This captures physics that NO linear combination of existing kernels
    can represent — genuinely new angular structure from the local environment.

    Architecture:
        1. Tensor product of L=2 kernels with themselves → new L=2 features
        2. Scalar gating from environment features
        3. Output: L=2 correction tensor
    """

    def __init__(self, n_scalar_features: int, n_kernels: int, hidden: int = 32):
        super().__init__()
        # Input irreps: n_kernels × L=2 tensors + scalars
        # Sum of all kernel T2s as one L=2 input, plus top individual kernels:
        # MC, Coulomb, BS, MopacCoulomb (QM charge EFG — independent angular info)
        self.input_irreps = o3.Irreps(f"{n_scalar_features}x0e + 5x2e")
        self.hidden_irreps = o3.Irreps(f"{hidden}x0e + {hidden}x2e")
        self.output_irreps = o3.Irreps("1x2e")

        self.tp1 = o3.FullyConnectedTensorProduct(
            self.input_irreps, self.input_irreps, self.hidden_irreps
        )
        self.tp2 = o3.FullyConnectedTensorProduct(
            self.hidden_irreps, self.hidden_irreps, self.output_irreps
        )

    def forward(self, scalars, kernel_t2_sum, mc_t2, coulomb_t2, bs_t2, mopac_coulomb_t2):
        """
        scalars: (batch, n_scalar_features)
        kernel_t2_sum: (batch, 5) — sum of all calculator T2s
        mc_t2, coulomb_t2, bs_t2: (batch, 5) — top 3 classical T2s
        mopac_coulomb_t2: (batch, 5) — QM charge EFG T2
        returns: (batch, 5) — L=2 correction
        """
        x = torch.cat([scalars, kernel_t2_sum, mc_t2, coulomb_t2, bs_t2,
                        mopac_coulomb_t2], dim=-1)
        h = self.tp1(x, x)
        out = self.tp2(h, h)
        return out


class ShieldingT2Model(nn.Module):
    """Combined model: kernel mixing + equivariant correction.

    prediction = mixing_head(scalars, kernels) + correction_head(scalars, kernels)

    The mixing head recovers Level A (linear) as a special case when
    the MLP learns constant weights. The correction head adds what
    linear mixing cannot capture.
    """

    def __init__(self, n_scalar_features: int, n_kernels: int,
                 hidden_mix: int = 64, hidden_corr: int = 32,
                 use_correction: bool = True):
        super().__init__()
        self.mixing = KernelMixingHead(n_scalar_features, n_kernels, hidden_mix)
        self.use_correction = use_correction
        if use_correction:
            self.correction = EquivariantCorrectionHead(
                n_scalar_features, n_kernels, hidden_corr
            )
            # Scaling gate: starts small so correction doesn't destabilize
            # the working mixing head. Its magnitude at convergence measures
            # how much L=2 structure is beyond linear kernel mixing.
            self.correction_scale = nn.Parameter(torch.tensor(0.1))

    def forward(self, scalars, kernel_t2s, mc_t2, coulomb_t2, bs_t2, mopac_coulomb_t2):
        """
        scalars: (batch, n_scalar_features)
        kernel_t2s: (batch, n_kernels, 5) — all kernel T2 vectors
        mc_t2, coulomb_t2, bs_t2: (batch, 5) — individual top classical kernels
        mopac_coulomb_t2: (batch, 5) — QM charge EFG kernel
        returns: (batch, 5) — predicted DFT T2
        """
        mixed = self.mixing(scalars, kernel_t2s)

        if self.use_correction:
            kernel_sum = kernel_t2s.sum(dim=1)  # (batch, 5)
            corr = self.correction(scalars, kernel_sum, mc_t2, coulomb_t2, bs_t2,
                                   mopac_coulomb_t2)
            return mixed + self.correction_scale * corr
        return mixed

    def parameter_count(self):
        return sum(p.numel() for p in self.parameters())


def make_model(n_scalar_features=13, n_kernels=46,
               use_correction=True) -> ShieldingT2Model:
    """Create model with sensible defaults."""
    model = ShieldingT2Model(
        n_scalar_features=n_scalar_features,
        n_kernels=n_kernels,
        use_correction=use_correction,
    )
    return model


if __name__ == "__main__":
    # Smoke test
    batch = 32
    scalars = torch.randn(batch, 14)
    kernels = torch.randn(batch, 46, 5)
    mc = torch.randn(batch, 5)
    co = torch.randn(batch, 5)
    bs = torch.randn(batch, 5)
    mopac_co = torch.randn(batch, 5)

    model = make_model(use_correction=False)
    print(f"Mixing-only model: {model.parameter_count()} parameters")
    out = model(scalars, kernels, mc, co, bs, mopac_co)
    print(f"Output shape: {out.shape}  (should be [{batch}, 5])")

    # With correction
    model_full = make_model(use_correction=True)
    print(f"Full model: {model_full.parameter_count()} parameters")
    out_full = model_full(scalars, kernels, mc, co, bs, mopac_co)
    print(f"Output shape: {out_full.shape}")

    # GPU test
    if torch.cuda.is_available():
        model_full = model_full.cuda()
        out_gpu = model_full(
            scalars.cuda(), kernels.cuda(),
            mc.cuda(), co.cuda(), bs.cuda(), mopac_co.cuda()
        )
        print(f"GPU output: {out_gpu.shape}, device={out_gpu.device}")
