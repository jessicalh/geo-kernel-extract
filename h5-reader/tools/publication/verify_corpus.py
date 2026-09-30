"""Load a list of trajectories and check their displayed data through Reader REST."""

import argparse
import json
import logging
import math
import time
from pathlib import Path

import httpx

from verify_reader import snapshot

LOG = logging.getLogger("reader-corpus-check")


def check_trajectory(client, path, require_prediction):
    started = time.monotonic()
    result = snapshot(client, path)
    if result["atom_count"] < 1 or result["frame_count"] < 1:
        raise AssertionError("Trajectory has no atoms or frames")
    for frame in result["frames"]:
        positions = frame["positions"]["positions"]
        if len(positions) != result["atom_count"]:
            raise AssertionError("Displayed atom count differs from loaded topology")
        if not all(
            math.isfinite(value) for atom in positions for value in atom["position"]
        ):
            raise AssertionError("Displayed coordinates contain non-finite values")
    if require_prediction and result["prediction"] is None:
        raise AssertionError(
            f"Shielding prediction unavailable: {result['model_input_error']}"
        )
    response = client.get("/health")
    response.raise_for_status()
    if not response.json()["ok"]:
        raise AssertionError("Reader health check failed after trajectory inspection")
    return {
        "path": str(path),
        "atoms": result["atom_count"],
        "frames": result["frame_count"],
        "sampled_frames": [frame["frame"]["frame"] for frame in result["frames"]],
        "prediction": result["prediction"],
        "seconds": round(time.monotonic() - started, 3),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument(
        "--paths", type=Path, required=True, help="One absolute LGS path per line"
    )
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--expected-count", type=int)
    parser.add_argument("--require-prediction", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    paths = [
        Path(line)
        for line in args.paths.read_text(encoding="utf-8-sig").splitlines()
        if line
    ]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("Trajectory list must be nonempty and contain no duplicates")
    if args.expected_count is not None and len(paths) != args.expected_count:
        raise ValueError(
            f"Expected {args.expected_count} trajectories, received {len(paths)}"
        )
    for path in paths:
        if not path.is_absolute() or not path.is_file():
            raise FileNotFoundError(f"Not an accessible absolute LGS path: {path}")

    report = {"requested": len(paths), "passed": [], "failure": None}
    try:
        with httpx.Client(base_url=args.url, timeout=600) as client:
            for index, path in enumerate(paths, 1):
                LOG.info("Checking %d/%d: %s", index, len(paths), path)
                try:
                    checked = check_trajectory(client, path, args.require_prediction)
                except Exception as error:
                    report["failure"] = {"path": str(path), "error": str(error)}
                    LOG.exception("Trajectory check failed: %s", path)
                    raise
                report["passed"].append(checked)
                LOG.info(
                    "PASS %d/%d: %d atoms, %d frames, prediction=%s, %.1fs",
                    index,
                    len(paths),
                    checked["atoms"],
                    checked["frames"],
                    checked["prediction"],
                    checked["seconds"],
                )
    finally:
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    LOG.info("All %d trajectories passed", len(paths))


if __name__ == "__main__":
    main()
