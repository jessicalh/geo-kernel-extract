"""Compare original and portable trajectory loads through a running Reader."""

import argparse
import json
import logging
import math
import time
from pathlib import Path

import httpx

LOG = logging.getLogger("reader-publication-check")


def prediction(client):
    client.post("/frame/set", json={"frame": 0}).raise_for_status()
    response = client.post(
        "/dashboard/metric",
        json={
            "descriptor_id": "ml:experimental_shielding_iso",
            "anchor": {"atom": 16},
            "modes": ["strip.scalar"],
        },
    )
    response.raise_for_status()
    signal = response.json()["id"]
    try:
        deadline = time.monotonic() + 180
        while time.monotonic() < deadline:
            response = client.get("/dashboard/display")
            response.raise_for_status()
            for track in response.json()["strip_tracks"]:
                if (
                    track["signal_id"] == signal
                    and track["valid"]
                    and track["valid"][0]
                ):
                    value = track["values"][0]
                    if not math.isfinite(value):
                        raise AssertionError("Non-finite shielding prediction")
                    return value
            time.sleep(0.2)
        raise TimeoutError("Shielding prediction did not finish within 180 seconds")
    finally:
        client.post("/dashboard/metric/remove", json={"id": signal}).raise_for_status()


def snapshot(client, path):
    response = client.post("/api/run/load", json={"path": str(path.resolve())})
    response.raise_for_status()
    response = client.post("/transform", json={"kind": "all_atom_fit"})
    response.raise_for_status()
    response = client.get("/protein/atoms")
    response.raise_for_status()
    atom_count = response.json()["count"]
    response = client.get("/frame/current")
    response.raise_for_status()
    frame_count = response.json()["count"]
    frames = []
    for frame in sorted({0, frame_count // 2, frame_count - 1}):
        client.post("/frame/set", json={"frame": frame}).raise_for_status()
        positions = client.post(
            "/positions", json={"frame": frame, "atoms": list(range(atom_count))}
        )
        positions.raise_for_status()
        catalog = client.get("/catalog")
        catalog.raise_for_status()
        tensors = []
        for atom in sorted(
            {0, min(16, atom_count - 1), atom_count // 2, atom_count - 1}
        ):
            response = client.get("/csa", params={"atom": atom})
            response.raise_for_status()
            tensors.append(response.json())
        timing = client.get("/frame/current")
        timing.raise_for_status()
        frames.append(
            {
                "frame": timing.json(),
                "positions": positions.json(),
                "tensors": tensors,
                "fields": {
                    item["id"]: item["availability"]
                    for item in catalog.json()["descriptors"]
                },
            }
        )
        LOG.info(
            "%s: frame %d, %d atoms, %d catalog fields",
            path.name,
            frame,
            atom_count,
            len(frames[-1]["fields"]),
        )
    state = client.get("/ui/state")
    state.raise_for_status()
    model = state.json()["experimentalShieldingMl"]
    return {
        "atom_count": atom_count,
        "frame_count": frame_count,
        "frames": frames,
        "model_input_ready": model["inferenceReady"],
        "model_input_error": model.get("inferenceError"),
        "prediction": prediction(client) if model["inferenceReady"] else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--portable", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    with httpx.Client(base_url=args.url, timeout=600) as client:
        original = snapshot(client, args.original)
        portable = snapshot(client, args.portable)
    original_prediction = original["prediction"]
    portable_prediction = portable["prediction"]
    prediction_matches = (
        original_prediction is None and portable_prediction is None
    ) or (
        original_prediction is not None
        and portable_prediction is not None
        and math.isclose(
            original_prediction, portable_prediction, rel_tol=1e-7, abs_tol=1e-6
        )
    )
    exact_keys = original.keys() - {"prediction"}
    equal = prediction_matches and all(
        original[key] == portable[key] for key in exact_keys
    )
    args.report.write_text(
        json.dumps(
            {"original": original, "portable": portable, "equal": equal},
            indent=2,
        )
        + "\n"
    )
    if not equal:
        raise AssertionError(f"Reader results differ; inspect {args.report}")
    LOG.info(
        "PASS: coordinates, times, field availability, sampled ORCA tensors and model results agree"
    )
    if original_prediction is None:
        LOG.warning(
            "Neither load has complete model inputs; no prediction was compared"
        )


if __name__ == "__main__":
    main()
