"""Capture hidden e3nn Gate outputs for Reader, without changing inference.

This is an optional model-side diagnostic, not a Reader runtime dependency.
Call capture_frame around the model's ordinary evaluation-mode forward call.
The supplied rotation maps extraction-frame column vectors into model space.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import torch
from e3nn.io import CartesianTensor
from e3nn.nn import Gate


def capture_frame(
    model: torch.nn.Module,
    forward: Callable[[], Any],
    *,
    frame_index: int,
    raw_positions: np.ndarray,
    atom_indices: Sequence[int],
    extraction_to_model_rotation: np.ndarray,
) -> tuple[Any, dict]:
    """Return the unchanged forward result and one frame of hidden 2e channels."""
    if any(module.training for module in model.modules()):
        raise ValueError("capture requires evaluation mode")
    positions = np.asarray(raw_positions, dtype=np.float64)
    atoms = np.asarray(atom_indices, dtype=np.int64)
    rotation = np.asarray(extraction_to_model_rotation, dtype=np.float64)
    if (
        positions.ndim != 2
        or positions.shape[1] != 3
        or not np.isfinite(positions).all()
    ):
        raise ValueError("raw_positions must be finite N x 3 extraction coordinates")
    if atoms.ndim != 1 or not len(atoms) or len(set(atoms)) != len(atoms):
        raise ValueError("atom_indices must be distinct and nonempty")
    if atoms.min() < 0 or atoms.max() >= len(positions):
        raise ValueError("atom index outside the loaded molecule")
    if rotation.shape != (3, 3) or not np.isfinite(rotation).all():
        raise ValueError("expected a finite 3 x 3 rotation")
    np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-8)
    if np.linalg.det(rotation) < 0:
        raise ValueError("expected a proper rotation")
    if not isinstance(frame_index, int) or frame_index < 0:
        raise ValueError("frame_index must be a nonnegative extraction frame index")

    cartesian = CartesianTensor("ij=ji")
    if str(cartesian) != "1x0e+1x2e":
        raise RuntimeError(f"unexpected symmetric Cartesian basis: {cartesian}")
    channels: list[dict] = []
    captured_names: set[str] = set()

    def hook(name: str, gate: Gate):
        def record(_module, _inputs, output):
            if output.ndim != 2 or output.shape[0] != len(positions):
                raise ValueError(f"{name}: hidden rows do not match the full molecule")
            channel_number = 0
            for (multiplicity, irrep), block in zip(
                gate.irreps_out, gate.irreps_out.slices(), strict=True
            ):
                if irrep.l != 2 or irrep.p != 1:
                    continue
                coefficients = output.detach()[atoms.tolist(), block].cpu().double()
                coefficients = coefficients.reshape(len(atoms), multiplicity, 5)
                for index in range(multiplicity):
                    label = f"{name} / 2e {channel_number}"
                    if label in captured_names:
                        raise ValueError(
                            f"{name} ran twice: capture one model frame at a time"
                        )
                    captured_names.add(label)
                    spherical = torch.cat(
                        (
                            torch.zeros((len(atoms), 1), dtype=torch.float64),
                            coefficients[:, index],
                        ),
                        dim=1,
                    )
                    in_model = cartesian.to_cartesian(spherical).numpy()
                    # T_model = R T_extraction R^T, so invert only this known
                    # transform. Reader will apply its own current display R.
                    tensors = rotation.T @ in_model @ rotation
                    if not np.isfinite(tensors).all():
                        raise ValueError(f"nonfinite activity in {label}")
                    channels.append(
                        {
                            "name": label,
                            "tensors": tensors[
                                :, [0, 0, 0, 1, 1, 2], [0, 1, 2, 1, 2, 2]
                            ].tolist(),
                        }
                    )
                    channel_number += 1

        return record

    handles = []
    try:
        for name, module in model.named_modules():
            if isinstance(module, Gate):
                handles.append(module.register_forward_hook(hook(name, module)))
        if not handles:
            raise ValueError("model has no e3nn Gate modules to capture")
        with torch.inference_mode():
            prediction = forward()
        if not channels:
            raise ValueError("forward call produced no hidden 2e channels")
    finally:
        for handle in handles:
            handle.remove()
    return prediction, {
        "frame_index": frame_index,
        "positions": positions[atoms].tolist(),
        "channels": channels,
    }


def write_capture(
    path: Path,
    *,
    model_name: str,
    atom_count: int,
    atom_indices: Sequence[int],
    frames: Sequence[dict],
) -> None:
    """Write the small display capture; no weights or input data are changed."""
    payload = {
        "schema_version": 1,
        "model": model_name,
        "coordinate_frame": "extraction",
        "atom_count": atom_count,
        "atom_indices": list(atom_indices),
        "frames": list(frames),
    }
    Path(path).write_text(
        json.dumps(payload, allow_nan=False, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
