"""Build the Reader-ready packages beside a provenance drive's source datasets."""

import argparse
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter

from make_reader_collection import one_lgs, source_runs


PUBLISHER = Path(__file__).with_name("package_trajectory.py")
LOG = logging.getLogger("reader-drive")
GIB = 1024**3


def frame_count(run):
    with os.scandir(run / "npys") as items:
        count = sum(item.name.startswith("frame_") and item.is_dir() for item in items)
    if count < 1:
        raise ValueError(f"No frame directories in {run / 'npys'}")
    return count


def complete_package(packages, key, frames):
    archive = packages / f"{key}.tar.xz"
    metadata = packages / f"{key}.json"
    if not archive.exists() and not metadata.exists():
        return False
    if not archive.is_file() or not metadata.is_file():
        raise ValueError(f"Incomplete package for {key}: {archive} / {metadata}")
    record = json.loads(metadata.read_text(encoding="utf-8"))
    if (record.get("key") != key or record.get("archive") != archive.name
            or record.get("entry_point") != "run.LGS"
            or record.get("frames") != frames
            or record.get("archive_bytes") != archive.stat().st_size
            or record.get("expanded_bytes", 0) < 1):
        raise ValueError(f"Package metadata does not match {archive}")
    return True


def require_space(directory, needed_bytes, reserve_gib):
    free_bytes = shutil.disk_usage(directory).free
    if free_bytes < needed_bytes + reserve_gib * GIB:
        raise OSError(
            f"Insufficient space at {directory}: {free_bytes / GIB:.1f} GiB free; "
            f"need {needed_bytes / GIB:.1f} GiB and retain {reserve_gib:g} GiB free"
        )
    LOG.info("Disk space at %s: %.1f GiB free", directory, free_bytes / GIB)


def copy_verified(source, destination):
    """Create a new output, then compare its bytes before publishing metadata."""
    with source.open("rb") as incoming, destination.open("xb") as outgoing:
        shutil.copyfileobj(incoming, outgoing, length=1024**2)
        outgoing.flush()
        os.fsync(outgoing.fileno())
    verify_copy(source, destination)
    LOG.info("Copied and verified %s: %d bytes", destination, source.stat().st_size)


def verify_copy(source, destination):
    with source.open("rb") as original, destination.open("rb") as copied:
        while True:
            block = original.read(1024**2)
            if copied.read(1024**2) != block:
                raise OSError(f"Copy verification failed: {destination}")
            if not block:
                break


def remove_local_package(working, key):
    # Never recurse: only these two outputs in the private local work directory.
    (working / f"{key}.tar.xz").unlink()
    (working / f"{key}.json").unlink()
    working.rmdir()


def publish_one(run, title, group, packages, work_root, fields, minimum_free_gib):
    key = f"{run.name.split('_', 1)[0].lower()}-reader-v1"
    frames = frame_count(run)
    working = work_root / key
    if complete_package(packages, key, frames):
        LOG.info("%s: already complete", key)
        if working.exists():
            if not complete_package(working, key, frames):
                raise ValueError(f"Incomplete retained local work at {working}")
            for filename in (f"{key}.tar.xz", f"{key}.json"):
                verify_copy(working / filename, packages / filename)
            remove_local_package(working, key)
            LOG.info("%s: matching retained local outputs removed", key)
        return

    if not working.exists():
        require_space(work_root, 0, minimum_free_gib)
        working.mkdir()
        subprocess.run([
            sys.executable, "-B", "-u", str(PUBLISHER), str(one_lgs(run)),
            "--fields", str(fields),
            "--output", str(working),
            "--key", key,
            "--title", title,
            "--description", f"{frames}-frame {group.lower()} trajectory",
            "--frames", str(frames),
            "--local",
            "--minimum-free-gib", str(minimum_free_gib),
        ], check=True, stderr=subprocess.STDOUT)
    if not complete_package(working, key, frames):
        raise RuntimeError(f"Incomplete local work at {working}; see the batch log")

    archive = working / f"{key}.tar.xz"
    metadata = working / f"{key}.json"
    require_space(packages, archive.stat().st_size + metadata.stat().st_size, 10)
    copy_verified(archive, packages / archive.name)
    copy_verified(metadata, packages / metadata.name)
    if not complete_package(packages, key, frames):
        raise RuntimeError(f"Publication did not create {key}")

    remove_local_package(working, key)
    LOG.info("%s: published; local package files removed", key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reader_root", type=Path)
    parser.add_argument("--fields", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True,
                        help="Local directory for one package at a time")
    parser.add_argument("--minimum-free-gib", type=float, default=100)
    parser.add_argument("--only", help="Publish one run name, such as bmr68 or trpcage")
    parser.add_argument("--group", choices=("Small proteins", "MD trajectories"),
                        help="Publish only this part of the source inventory")
    args = parser.parse_args()
    root = args.reader_root.resolve()
    packages = root / "packages"
    work_root = args.work_dir.resolve()
    if work_root.is_relative_to(root) or root.is_relative_to(work_root):
        parser.error("the work directory must be outside the provenance tree")
    if root.drive and root.drive.casefold() == work_root.drive.casefold():
        parser.error("local work must be on a different drive from provenance")
    if args.minimum_free_gib < 0:
        parser.error("minimum free disk space must be non-negative")
    fields = args.fields.resolve()
    if not fields.is_file():
        parser.error(f"field list not found: {fields}")
    runs = list(source_runs(root))
    if (sum(group == "Small proteins" for _, _, group in runs) != 3
            or sum(group == "MD trajectories" for _, _, group in runs) != 173):
        raise ValueError("Expected all 3 small proteins and 173 MD trajectories")
    if args.group:
        runs = [item for item in runs if item[2] == args.group]
    if args.only:
        runs = [item for item in runs if item[0].name.split("_", 1)[0].lower() == args.only.lower()]
        if len(runs) != 1:
            parser.error(f"Expected one run named {args.only}, found {len(runs)}")
    packages.mkdir(exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)
    LOG.info("Publishing %d Reader packages to %s; local work at %s", len(runs), packages, work_root)
    for number, (run, title, group) in enumerate(runs, 1):
        started = perf_counter()
        LOG.info("[%d/%d] %s", number, len(runs), run.name)
        publish_one(run, title, group, packages, work_root, fields, args.minimum_free_gib)
        LOG.info("[%d/%d] finished in %.1f seconds", number, len(runs), perf_counter() - started)
    LOG.info("All %d requested packages are complete", len(runs))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                        format="%(asctime)s %(levelname)s %(message)s")
    try:
        main()
    except Exception:
        LOG.exception("Batch stopped; published packages and source data are left in place")
        sys.exit(1)
