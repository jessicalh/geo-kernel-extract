"""nmr_extract — typed reader for NMR shielding extraction data.

    from nmr_extract import load
    p = load("path/to/extraction")

    p.biot_savart.shielding.T2          # (N, 5)
    p.biot_savart.shielding.irreps      # Irreps("1x0e + 1x1o + 1x2e")
    p.biot_savart.shielding.torch()     # torch.Tensor (N, 9)
    p.ring_contributions.bs.T2          # (P, 5), per-ring
"""

from e3nn.o3 import Irreps

from ._protein import load, Protein
from ._trajectory import load_trajectory, TrajectoryData, TrajectoryRollup, BondRollup
from ._types import RingType, BondCategory, N_RING_TYPES, N_BOND_CATEGORIES
from ._tensors import (
    SphericalTensor,
    ShieldingTensor,
    EFGTensor,
    VectorField,
    PerRingTypeT0,
    PerRingTypeT2,
    PerBondCategoryT2,
    MopacScalars,
    MopacGlobal,
    AIMNet2Charges,
    AIMNet2AimEmbedding,
    AIMNet2ChargeSensitivity,
)
from ._ring import RingContributions, RingGeometry
from ._catalog import CATALOG, ArraySpec
from ._protein import AIMNet2Group, WaterPolarizationGroup, EeqGroup

__all__ = [
    "load",
    "Protein",
    "Irreps",
    "RingType",
    "BondCategory",
    "SphericalTensor",
    "ShieldingTensor",
    "EFGTensor",
    "VectorField",
    "PerRingTypeT0",
    "PerRingTypeT2",
    "PerBondCategoryT2",
    "RingContributions",
    "RingGeometry",
    "CATALOG",
    "ArraySpec",
    "AIMNet2Charges",
    "AIMNet2AimEmbedding",
    "AIMNet2ChargeSensitivity",
    "MopacScalars",
    "MopacGlobal",
    "AIMNet2Group",
    "WaterPolarizationGroup",
    "EeqGroup",
    "load_trajectory",
    "TrajectoryData",
    "TrajectoryRollup",
    "BondRollup",
]
