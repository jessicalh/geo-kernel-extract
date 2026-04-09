"""Load calibration.toml into a frozen Config dataclass.

Every number in the pipeline traces here.  CLI args override TOML values.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomllib
    except ModuleNotFoundError:
        import tomli as tomllib  # type: ignore[no-redef]


@dataclass(frozen=True)
class PathsConfig:
    features: Path
    runs: Path
    sdk: Path


@dataclass(frozen=True)
class DataConfig:
    min_matched_atoms: int
    val_fraction: float
    seed: int
    max_proteins: int


@dataclass(frozen=True)
class KernelsConfig:
    ring_calculators: list[str]
    bond_calculators: list[str]
    bond_calculators_mopac: list[str]
    total_calculators: list[str]
    efg_calculators: list[str]
    per_ring_calculators: list[str]
    per_ring_k: int
    rbf_centers: list[float] = None
    rbf_sigma: float = 2.0
    ang_magic_theta: float = 0.9553


@dataclass(frozen=True)
class ScalarsConfig:
    top_k_rings: int
    sentinel_distance: float
    inverse_epsilon: float
    bond_order_floor: float
    bond_order_pi: float


@dataclass(frozen=True)
class NormConfig:
    kernel_std_floor: float
    scalar_std_floor: float
    gate_threshold_min_atoms: int


@dataclass(frozen=True)
class TrainingConfig:
    epochs: int
    batch_size: int
    lr: float
    weight_decay: float
    scheduler: str


@dataclass(frozen=True)
class ModelConfig:
    hidden_mix: int
    hidden_corr: int
    use_correction: bool
    correction_scale_init: float
    correction_l2_kernels: list[int]


@dataclass(frozen=True)
class AnalysisConfig:
    ridge_lambda: float
    distance_bands: list[tuple[float, float]]
    max_forward_steps: int
    min_atoms_per_band: int


@dataclass(frozen=True)
class SecondaryConfig:
    output_dir: Path
    strata: list[str]
    redundancy_threshold: float
    ridge_lambda: float


@dataclass(frozen=True)
class Config:
    paths: PathsConfig
    data: DataConfig
    kernels: KernelsConfig
    scalars: ScalarsConfig
    normalization: NormConfig
    training: TrainingConfig
    model: ModelConfig
    analysis: AnalysisConfig
    secondary: SecondaryConfig


_SECONDARY_DEFAULTS = {
    "output_dir": "output/secondary",
    "strata": ["hie_only", "phe_only", "tyr_only", "trp_only", "no_hie", "all"],
    "redundancy_threshold": 0.95,
}


def _load_secondary(raw_sec: dict, raw_all: dict) -> SecondaryConfig:
    """Build SecondaryConfig with defaults for missing keys."""
    merged = {**_SECONDARY_DEFAULTS, **raw_sec}
    return SecondaryConfig(
        output_dir=Path(merged["output_dir"]),
        strata=merged["strata"],
        redundancy_threshold=merged["redundancy_threshold"],
        ridge_lambda=raw_all.get("analysis", {}).get("ridge_lambda", 0.01),
    )


def load_config(path: str | Path) -> Config:
    """Load calibration.toml and return a frozen Config."""
    path = Path(path)
    with open(path, "rb") as f:
        raw = tomllib.load(f)

    return Config(
        paths=PathsConfig(
            features=Path(raw["paths"]["features"]),
            runs=Path(raw["paths"]["runs"]),
            sdk=Path(raw["paths"]["sdk"]),
        ),
        data=DataConfig(**raw["data"]),
        kernels=KernelsConfig(**raw["kernels"]),
        scalars=ScalarsConfig(**raw["scalars"]),
        normalization=NormConfig(**raw["normalization"]),
        training=TrainingConfig(**raw["training"]),
        model=ModelConfig(
            **{k: v for k, v in raw["model"].items()
               if k != "correction_l2_kernels"},
            correction_l2_kernels=raw["model"]["correction_l2_kernels"],
        ),
        analysis=AnalysisConfig(
            ridge_lambda=raw["analysis"]["ridge_lambda"],
            distance_bands=[tuple(b) for b in raw["analysis"]["distance_bands"]],
            max_forward_steps=raw["analysis"]["max_forward_steps"],
            min_atoms_per_band=raw["analysis"]["min_atoms_per_band"],
        ),
        secondary=_load_secondary(raw.get("secondary", {}), raw),
    )
