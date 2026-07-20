#!/usr/bin/env python3
"""Export the sealed F003/R004 model for the Reader's native LibTorch helper."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import sys
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch
from e3nn import o3
from e3nn.math import soft_one_hot_linspace


MODEL_FILE = "model.ts"
MANIFEST_FILE = "manifest.json"
LABEL_KEYS = (
    "element",
    "iupac_atom",
    "iupac_residue",
    "amber_residue",
    "terminal_state",
    "residue_variant",
    "formal_charge",
    "polar_h_kind",
    "backbone_role",
    "planar_group",
    "locant",
    "ring_position",
    "aromatic",
    "exchangeable",
    "iupac_pair",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def import_package(package: Path):
    sys.path.insert(0, str(package))
    try:
        infer = importlib.import_module("infer_f003")
        loader = importlib.import_module("f003_loader")
        model_module = importlib.import_module("f003_model")
    finally:
        sys.path.pop(0)
    return infer, loader, model_module


class ReaderModelWrapper(torch.nn.Module):
    """F003 with frozen normalization and a precomputed graph boundary."""

    def __init__(
        self,
        core: torch.nn.Module,
        normalizer: Mapping[str, np.ndarray],
        project_to_e3nn: np.ndarray,
    ) -> None:
        super().__init__()
        self.core = core
        self.register_buffer(
            "scalar_mean",
            torch.as_tensor(normalizer["scalar_mean"], dtype=torch.float32),
        )
        self.register_buffer(
            "scalar_scale",
            torch.as_tensor(normalizer["scalar_scale"], dtype=torch.float32),
        )
        self.register_buffer(
            "l1_scale",
            torch.as_tensor(normalizer["l1_scale"], dtype=torch.float32),
        )
        self.register_buffer(
            "l2_scale",
            torch.as_tensor(normalizer["l2_scale"], dtype=torch.float32),
        )
        self.register_buffer(
            "pair_iso_mean",
            torch.as_tensor(normalizer["pair_iso_mean"], dtype=torch.float32),
        )
        self.register_buffer(
            "pair_iso_scale",
            torch.as_tensor(normalizer["pair_iso_scale"], dtype=torch.float32),
        )
        self.register_buffer(
            "pair_l2_scale",
            torch.as_tensor(normalizer["pair_l2_scale"], dtype=torch.float32),
        )
        self.register_buffer(
            "project_to_e3nn",
            torch.as_tensor(project_to_e3nn, dtype=torch.float32),
        )
        self.iupac_pair_column = LABEL_KEYS.index("iupac_pair")

    def forward(
        self,
        pos: torch.Tensor,
        l1: torch.Tensor,
        l1_valid: torch.Tensor,
        l2: torch.Tensor,
        l2_valid: torch.Tensor,
        scalars: torch.Tensor,
        scalar_valid: torch.Tensor,
        applicability: torch.Tensor,
        label_ids: torch.Tensor,
        edge_src: torch.Tensor,
        edge_dst: torch.Tensor,
        radial: torch.Tensor,
    ) -> torch.Tensor:
        scalars_normalized = (scalars - self.scalar_mean[None, :]) / self.scalar_scale[None, :]
        scalars_normalized = torch.where(
            scalar_valid > 0.5,
            scalars_normalized,
            torch.zeros_like(scalars_normalized),
        )
        l1_normalized = l1 / self.l1_scale[None, :, None]
        l1_normalized = torch.where(
            l1_valid[:, :, None] > 0.5,
            l1_normalized,
            torch.zeros_like(l1_normalized),
        )
        l2_normalized = l2 / self.l2_scale[None, :, None]
        l2_normalized = torch.where(
            l2_valid[:, :, None] > 0.5,
            l2_normalized,
            torch.zeros_like(l2_normalized),
        )

        batch: Dict[str, torch.Tensor] = {
            "pos": pos,
            "scalars": scalars_normalized,
            "applicability": applicability,
            "l1": l1_normalized,
            "l2": l2_normalized,
        }
        for column, key in enumerate(LABEL_KEYS):
            batch[f"label:{key}"] = label_ids[:, column]

        node = self.core._pack_nodes(batch)
        relative = pos[edge_dst] - pos[edge_src]
        spherical_harmonics = o3.spherical_harmonics(
            self.core.irreps_sh,
            relative,
            normalize=True,
            normalization="component",
        )

        for convolution, gate in zip(self.core.convolutions, self.core.gates):
            weights = convolution.radial(radial)
            messages = convolution.tp(
                node[edge_src],
                spherical_harmonics,
                weights,
            )
            aggregate = node.new_zeros((node.shape[0], convolution.irreps_out.dim))
            aggregate.index_add_(0, edge_dst, messages)
            degree = node.new_zeros((node.shape[0], 1))
            degree.index_add_(
                0,
                edge_dst,
                torch.ones_like(edge_dst, dtype=node.dtype).unsqueeze(1),
            )
            node = gate(
                convolution.self_interaction(node)
                + aggregate / degree.clamp_min(1.0)
            )

        pair_ids = label_ids[:, self.iupac_pair_column]
        prediction = self.core.shared_readout(node)
        prediction = prediction + self.core.pair_delta_readout(node, pair_ids)
        iso = (
            prediction[:, 0:1] * self.pair_iso_scale[pair_ids, None]
            + self.pair_iso_mean[pair_ids, None]
        )
        l2_project = (
            prediction[:, 1:6]
            * self.pair_l2_scale[pair_ids, None]
        ) @ self.project_to_e3nn
        return torch.cat((iso, l2_project), dim=1)


def build_edges(
    model: torch.nn.Module, pos: torch.Tensor
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    src, dst, _relative, length = model._edges(
        pos,
        [(0, int(pos.shape[0]))],
    )
    radial = soft_one_hot_linspace(
        length,
        0.0,
        model.radius,
        model.radial_dim,
        basis="smooth_finite",
        cutoff=True,
    )
    cutoff = 0.5 * (
        torch.cos(
            math.pi * torch.clamp(length / model.radius, max=1.0)
        )
        + 1.0
    )
    return src, dst, radial * cutoff[:, None]


def graph_tensors(
    graph,
    model: torch.nn.Module,
    indices: Optional[np.ndarray] = None,
) -> Tuple[torch.Tensor, ...]:
    def selected(values: np.ndarray) -> np.ndarray:
        return values if indices is None else values[indices]

    pos = torch.as_tensor(selected(graph.pos), dtype=torch.float32)
    src, dst, radial = build_edges(model, pos)
    labels = np.column_stack(
        [selected(graph.label_ids[key]) for key in LABEL_KEYS]
    ).astype(np.int64)
    return (
        pos,
        torch.as_tensor(selected(graph.l1), dtype=torch.float32),
        torch.as_tensor(selected(graph.l1_valid), dtype=torch.float32),
        torch.as_tensor(selected(graph.l2), dtype=torch.float32),
        torch.as_tensor(selected(graph.l2_valid), dtype=torch.float32),
        torch.as_tensor(selected(graph.scalars), dtype=torch.float32),
        torch.as_tensor(selected(graph.scalar_valid), dtype=torch.float32),
        torch.as_tensor(selected(graph.applicability), dtype=torch.float32),
        torch.as_tensor(labels, dtype=torch.long),
        src,
        dst,
        radial,
    )


def trace_indices(graph, count: int = 24) -> np.ndarray:
    positions = np.asarray(graph.pos, dtype=np.float32)
    center = positions.mean(axis=0, keepdims=True)
    distance = np.linalg.norm(positions - center, axis=1)
    return np.argsort(distance)[: min(count, len(positions))]


def write_raw(path: Path, tensor: torch.Tensor) -> None:
    array = tensor.detach().cpu().numpy()
    path.parent.mkdir(parents=True, exist_ok=True)
    array.tofile(path)


def write_smoke(
    root: Path,
    name: str,
    tensors: Tuple[torch.Tensor, ...],
    prediction: torch.Tensor,
) -> Dict[str, int]:
    names = (
        "pos",
        "l1",
        "l1_valid",
        "l2",
        "l2_valid",
        "scalars",
        "scalar_valid",
        "applicability",
        "label_ids",
        "edge_src",
        "edge_dst",
        "radial",
    )
    input_dir = root / "smoke" / name / "input"
    for field, tensor in zip(names, tensors):
        write_raw(input_dir / f"{field}.bin", tensor)

    n = int(tensors[0].shape[0])
    e = int(tensors[9].numel())
    shapes = {
        "N": n,
        "E": e,
        "C1": int(tensors[1].shape[1]),
        "C2": int(tensors[3].shape[1]),
        "C0": int(tensors[5].shape[1]),
        "A": int(tensors[7].shape[1]),
        "label_count": int(tensors[8].shape[1]),
        "radial_dim": int(tensors[11].shape[1]),
    }
    with (input_dir / "shapes.txt").open("w", encoding="ascii") as handle:
        for key, value in shapes.items():
            handle.write(f"{key} {value}\n")
    write_raw(root / "smoke" / name / "reference_output.bin", prediction)
    return {
        "atoms": n,
        "edges": e,
        "scalar_channels": shapes["C0"],
        "applicability_channels": shapes["A"],
        "l1_channels": shapes["C1"],
        "l2_channels": shapes["C2"],
        "label_columns": shapes["label_count"],
    }


def max_abs(left: torch.Tensor, right: torch.Tensor) -> float:
    return float(torch.max(torch.abs(left - right)).item())


def export(package: Path, output: Path) -> Dict[str, object]:
    infer, loader, model_module = import_package(package)
    bundle = infer.Bundle(package)
    device = torch.device("cpu")
    model = bundle.build_model(
        device,
        edge_chunk_atoms=-1,
        message_chunk_edges=-1,
    )
    normalizer = infer._normalization_arrays(bundle.normalizer)
    project_to_e3nn = np.asarray(
        bundle.metadata["project_to_e3nn"],
        dtype=np.float32,
    )
    wrapper = ReaderModelWrapper(
        model,
        normalizer,
        project_to_e3nn,
    ).eval()

    fixture_member = loader.discover_members(
        package / "fixture" / "input"
    )[0]
    graph = loader.load_graph(
        fixture_member,
        fixture_member.frame_dirs[0],
        bundle.contract,
        bundle.metadata,
        bundle.vocabs,
    )
    full_tensors = graph_tensors(graph, model)
    small_tensors = graph_tensors(graph, model, trace_indices(graph))

    with torch.inference_mode():
        original_prediction, _ = infer.infer_graph(
            model,
            graph,
            normalizer,
            project_to_e3nn,
            device,
        )
        original = torch.as_tensor(original_prediction, dtype=torch.float32)
        eager_full = wrapper(*full_tensors)
        traced = torch.jit.trace(
            wrapper,
            small_tensors,
            strict=False,
            check_trace=False,
        )
        traced = torch.jit.freeze(traced.eval())
        traced_small = traced(*small_tensors)
        eager_small = wrapper(*small_tensors)
        traced_full = traced(*full_tensors)

    validation = {
        "original_vs_wrapper_max_abs": max_abs(original, eager_full),
        "small_eager_vs_trace_max_abs": max_abs(eager_small, traced_small),
        "full_eager_vs_trace_max_abs": max_abs(eager_full, traced_full),
    }
    if max(validation.values()) > 5.0e-4:
        raise RuntimeError(f"TorchScript export parity failed: {validation}")

    output.mkdir(parents=True, exist_ok=True)
    model_path = output / MODEL_FILE
    traced.save(str(model_path))
    smoke = {
        "small": write_smoke(
            output,
            "small",
            small_tensors,
            traced_small,
        ),
        "full": write_smoke(
            output,
            "full",
            full_tensors,
            traced_full,
        ),
    }

    metadata = bundle.metadata
    contract = bundle.contract
    manifest = {
        "name": "Experimental Shielding ML",
        "bundle_version": "F003-R004-v1-reader-runtime",
        "bundle_date": "2026-07-19",
        "source_bundle_id": bundle.bundle_id,
        "source_model_sha256": metadata["model_state"]["sha256"],
        "model_file": MODEL_FILE,
        "model_sha256": sha256_file(model_path),
        "target": metadata["target"],
        "output_contract": {
            "dtype": "float32",
            "shape": "N x 6",
            "columns": metadata["output_columns"],
            "tensor_basis": "project 0e plus traceless 2e",
            "t1": "not predicted and not synthesized",
        },
        "models": [
            {
                "id": "f003_r004",
                "label": "Experimental Shielding ML",
                "role": "July 2026 no-MOPAC trajectory-compatible model",
                "model_file": MODEL_FILE,
                "input_preset": "f003_no_mopac_common_sense",
                "selection": "all loaded July-contract runs",
                "training": {
                    "run": metadata["training_run"],
                    "best_epoch": metadata["best_epoch"],
                    "producer_commit": metadata["training_producer_commit"],
                },
            }
        ],
        "inference_schema": {
            "version": 2,
            "model_id": "f003_r004",
            "graph": {
                "radius_angstrom": metadata["architecture"]["radius_A"],
                "max_neighbors": metadata["architecture"]["max_neighbors"],
                "radial_dim": metadata["architecture"]["radial_dim"],
                "radial_basis": "e3nn smooth_finite with cosine cutoff",
            },
            "label_keys": list(LABEL_KEYS),
            "label_vocabs": bundle.vocabs,
            "numeric_features": contract["numeric_features"],
            "expected_channels": {
                "scalars": len(metadata["channels"]["scalar"]),
                "applicability": len(metadata["channels"]["applicability"]),
                "l1": len(metadata["channels"]["l1"]),
                "l2": len(metadata["channels"]["l2"]),
            },
            "expected_channel_names": {
                key: metadata["channels"][key]
                for key in ("scalar", "applicability", "l1", "l2")
            },
            "ring_type_order": metadata["ring_type_order"],
            "ring_active": metadata["ring_active"],
            "ring_intensity_nA_per_T": metadata["ring_intensity_nA_per_T"],
            "source_contract_schema": contract["schema_version"],
        },
        "runtime_contract": {
            "python_required": False,
            "helper_protocol": 2,
            "input_files": [
                "pos.bin",
                "l1.bin",
                "l1_valid.bin",
                "l2.bin",
                "l2_valid.bin",
                "scalars.bin",
                "scalar_valid.bin",
                "applicability.bin",
                "label_ids.bin",
                "edge_src.bin",
                "edge_dst.bin",
                "radial.bin",
                "shapes.txt",
            ],
        },
        "export_validation": validation,
        "smoke_fixture": smoke,
    }
    manifest_path = output / MANIFEST_FILE
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "model": str(model_path),
        "manifest": str(manifest_path),
        "model_bytes": model_path.stat().st_size,
        "validation": validation,
        "smoke": smoke,
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    package = Path(args.package).resolve()
    output = Path(args.output).resolve()
    if package == output or package in output.parents or output in package.parents:
        raise RuntimeError("source package and Reader runtime may not overlap")
    print(json.dumps(export(package, output), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
