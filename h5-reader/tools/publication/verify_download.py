"""Exercise a published bundle through Reader's catalog and real HTTPS endpoint."""

import argparse
import json
import logging
from pathlib import Path
import time

import httpx

LOG = logging.getLogger("reader-download-check")


def wait_for(description, check, timeout=1800):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        result = check()
        if result:
            LOG.info("PASS: %s", description)
            return result
        time.sleep(0.25)
    raise TimeoutError(description)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    logging.getLogger("httpx").setLevel(logging.WARNING)
    with httpx.Client(base_url=args.url, timeout=600) as client:

        def state():
            response = client.get("/api/trajectories")
            response.raise_for_status()
            value = response.json()
            if value["error"]:
                raise RuntimeError(value["error"])
            return value

        def post(action):
            response = client.post(
                f"/api/trajectories/{action}", json={"key": args.key}
            )
            response.raise_for_status()
            return response

        initial = state()
        entry = next(row for row in initial["entries"] if row["key"] == args.key)
        if entry["included"] or entry["downloaded"] or initial["busy"]:
            raise ValueError(
                "Use an idle Reader with this bundle neither installed nor cached"
            )
        cache_root = Path(initial["cache_root"])
        post("show")
        post("open")

        def transfer_started():
            current = state()
            if not current["busy"]:
                raise AssertionError("Transfer ended before the cancellation check")
            return any(
                path.stat().st_size >= 1024 * 1024
                for path in cache_root.glob(".partial-*/download.archive*")
            )

        wait_for("at least 1 MiB transferred and UI still responds", transfer_started)
        post("cancel")
        wait_for("download cancelled", lambda: not state()["busy"])
        if list(cache_root.glob(".partial-*")):
            raise AssertionError("Temporary files remain after cancellation")
        if (cache_root / args.key).exists():
            raise AssertionError("Cancelled download published an incomplete bundle")
        post("open")
        phases = set()

        def finished():
            current = state()
            phases.add(current["phase"])
            return current if not current["busy"] else None

        final = wait_for("retry downloaded, expanded and opened the bundle", finished)
        if 2 not in phases:
            raise AssertionError("Did not observe the extraction phase")
        row = next(row for row in final["entries"] if row["key"] == args.key)
        if not row["downloaded"]:
            raise AssertionError("Completed bundle is not in the catalog cache")
        response = client.get("/frame/current")
        response.raise_for_status()
        if response.json()["count"] != entry["frames"]:
            raise AssertionError("Downloaded trajectory has the wrong frame count")
        response = client.post("/api/trajectories/clear", json={"key": args.key})
        if response.status_code != 409:
            raise AssertionError(
                "Currently loaded data was not protected from cache clearing"
            )
        args.report.write_text(
            json.dumps(
                {"initial": initial, "final": final, "phases": sorted(phases)}, indent=2
            )
            + "\n",
            encoding="utf-8",
        )
        LOG.info(
            "PASS: complete real bundle retained; currently loaded data cannot be cleared"
        )


if __name__ == "__main__":
    main()
