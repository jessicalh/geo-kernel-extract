"""Publish a solid tar.xz containing the Reader's ordinary files and relative LGS."""

import argparse
import io
import json
import logging
import lzma
import os
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
from time import perf_counter

LOG = logging.getLogger("reader-package")
XZ_DICTIONARY_BYTES = 1024 * 1024 * 1024


def resolve(root, value):
    path = Path(value)
    return path if path.is_absolute() else root / path


def read_arrays(extraction, fields, expected_frames):
    with os.scandir(extraction / "npys") as entries:
        frames = sorted(
            Path(entry.path)
            for entry in entries
            if entry.name.startswith("frame_") and entry.is_dir()
        )
    if len(frames) != expected_frames:
        raise ValueError(f"Expected {expected_frames} NPY frames, found {len(frames)}")
    arrays = []
    byte_count = 0
    for index, frame in enumerate(frames, 1):
        number = frame.name.removeprefix("frame_")
        if not number.isascii() or not number.isdecimal():
            raise ValueError(f"Invalid NPY frame directory: {frame}")
        has_positions = False
        with os.scandir(frame) as entries:
            for entry in entries:
                if entry.name not in fields or not entry.is_file():
                    continue
                data = Path(entry.path).read_bytes()
                arrays.append((f"extraction/npys/{frame.name}/{entry.name}", data))
                byte_count += len(data)
                has_positions |= entry.name == "pos.npy"
        if not has_positions:
            raise FileNotFoundError(frame / "pos.npy")
        if index == 1 or index % 100 == 0 or index == len(frames):
            LOG.info(
                "Read %d / %d frames into memory: %d arrays, %d bytes",
                index,
                len(frames),
                len(arrays),
                byte_count,
            )
    # Reorder bytes already in memory, not source-file visits, for solid compression.
    arrays.sort(key=lambda item: (item[0].rsplit("/", 1)[1], item[0]))
    return arrays


def write_archive(path, members):
    """Compress the collected bytes without reopening or querying source files."""
    filters = [
        {
            "id": lzma.FILTER_LZMA2,
            "preset": 9 | lzma.PRESET_EXTREME,
            "dict_size": XZ_DICTIONARY_BYTES,
            "lc": 1,
            "lp": 3,
            "pb": 3,
        }
    ]
    LOG.info(
        "Solid XZ: preset 9 extreme, dictionary %d bytes, lc=1 lp=3 pb=3",
        XZ_DICTIONARY_BYTES,
    )
    with lzma.open(path, "wb", filters=filters) as compressed:
        with tarfile.open(fileobj=compressed, mode="w|") as bundled:
            for index, (name, data) in enumerate(members, 1):
                if index == 1 or index % 1000 == 0 or index == len(members):
                    LOG.info("tar.xz %d / %d files: %s", index, len(members), name)
                info = tarfile.TarInfo(name)
                info.size = len(data)
                bundled.addfile(info, io.BytesIO(data))


def package(args):
    started = perf_counter()
    if (
        not args.key
        or args.key in {".", ".."}
        or any(character in args.key for character in "/\\:")
    ):
        raise ValueError("Publication key must be a single directory name")
    if args.frames <= 0:
        raise ValueError("Expected frame count must be positive")
    source_lgs = args.lgs.resolve()
    original_bytes = source_lgs.read_bytes()
    manifest = json.loads(original_bytes.decode("utf-8-sig"))
    if manifest["kind"] != "trajectory":
        raise ValueError("This publisher accepts trajectory LGS files only")
    source = manifest["trajectory"]
    root = source_lgs.parent
    extraction = resolve(root, source["extraction_dir"])
    fields = {
        line.strip()
        for line in args.fields.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    }
    if "pos.npy" not in fields:
        raise ValueError("Reader field list is missing pos.npy")
    if any(
        not name.endswith(".npy") or any(character in name for character in "/\\:")
        for name in fields
    ):
        raise ValueError("Reader fields must be job-local .npy filenames")
    archive = args.output / f"{args.key}.tar.xz"
    metadata = args.output / f"{args.key}.json"
    if archive.exists() or metadata.exists():
        raise FileExistsError(f"Publication already exists: {archive}")
    if args.install_root and (args.install_root / args.key).exists():
        raise FileExistsError(
            f"Installed example already exists: {args.install_root / args.key}"
        )

    protected = {root, extraction.resolve()}
    members = []
    LOG.info("Reading publication contents into memory: %s", args.key)
    with os.scandir(extraction) as entries:
        for entry in entries:
            path = Path(entry.path)
            if path.suffix.lower() in {".npy", ".json", ".h5"} and entry.is_file():
                members.append((f"extraction/{entry.name}", path.read_bytes()))
    members.sort(key=lambda item: item[0])
    arrays = read_arrays(extraction, fields, args.frames)

    # Only publication metadata changes; arrays and ORCA outputs are copied as bytes.
    for name in ("md_dir", "topology_top"):
        source.pop(name, None)
    for name in ("trajectory_h5", "extraction_manifest"):
        original = resolve(root, source[name])
        expected = extraction / original.name
        if original.resolve() != expected.resolve() or not expected.is_file():
            raise ValueError(
                f"{name} is not a file in the declared extraction: {original}"
            )
        source[name] = f"extraction/{original.name}"
    source["extraction_dir"] = "extraction"
    if source.get("reference_pdb"):
        reference = resolve(root, source["reference_pdb"])
        members.append(("reference.pdb", reference.read_bytes()))
        protected.add(reference.resolve().parent)
        source["reference_pdb"] = "reference.pdb"
    for frame in manifest.get("dft", {}).get("frames", []):
        meta = resolve(root, frame["meta_json"])
        protected.add(meta.resolve().parent)
        meta_bytes = meta.read_bytes()
        record = json.loads(meta_bytes.decode("utf-8-sig"))
        output = record["files"]["out_primary"]
        if (
            not output
            or output in {".", ".."}
            or any(character in output for character in "/\\:")
        ):
            raise ValueError(
                f"DFT output is not a job-local filename: {meta}: {output}"
            )
        destination = f"dft/frame_{frame['frame_index']:06d}"
        members.extend(
            [
                (f"{destination}/{meta.name}", meta_bytes),
                (f"{destination}/{output}", (meta.parent / output).read_bytes()),
            ]
        )
        frame["meta_json"] = f"{destination}/{meta.name}"

    destinations = {"run.LGS", "provenance/source.LGS"}
    for destination, _ in members:
        if destination in destinations:
            raise ValueError(f"Duplicate cache path: {destination}")
        destinations.add(destination)
    for target in (args.output, args.install_root):
        if target and any(target.resolve().is_relative_to(path) for path in protected):
            raise ValueError(
                f"Publication output/install root must be outside source directories: {target}"
            )

    manifest_bytes = (json.dumps(manifest, indent=2) + "\n").encode()
    members = (
        [("run.LGS", manifest_bytes), ("provenance/source.LGS", original_bytes)]
        + members
        + arrays
    )
    expanded_bytes = sum(len(data) for _, data in members)
    collected = perf_counter()
    LOG.info(
        "Collected %d files, %d bytes in %.3f seconds; compression reads only memory",
        len(members),
        expanded_bytes,
        collected - started,
    )
    LOG.info(
        "Packing %s: %d frames, %d selected fields, %d NPY arrays, %d source files",
        args.key,
        args.frames,
        len(fields),
        len(arrays),
        len(members) - 1,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(
        prefix=f".{args.key}-", suffix=".partial", dir=args.output
    ) as temporary:
        staging = Path(temporary)
        pending_archive = staging / archive.name
        write_archive(pending_archive, members)
        LOG.info("Compression finished in %.3f seconds", perf_counter() - collected)
        archive_bytes = pending_archive.stat().st_size
        entry = {
            "key": args.key,
            "title": args.title,
            "description": args.description,
            "frames": args.frames,
            "entry_point": "run.LGS",
            "url": f"https://semantic.construction/files/{archive.name}",
            "expanded_bytes": expanded_bytes,
            "archive_bytes": archive_bytes,
        }
        pending_metadata = staging / metadata.name
        pending_metadata.write_text(
            json.dumps(entry, indent=2) + "\n", encoding="utf-8"
        )
        pending_archive.replace(archive)
        try:
            pending_metadata.replace(metadata)
        except BaseException:
            archive.unlink()
            raise
    LOG.info(
        "Archive ready: %s; tar.xz %d bytes; expanded %d bytes; %d cache files; %d NPY arrays",
        archive,
        archive_bytes,
        expanded_bytes,
        len(members),
        len(arrays),
    )
    if args.install_root:
        destination = args.install_root / args.key
        args.install_root.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(
            prefix=f".{args.key}-", suffix=".partial", dir=args.install_root
        ) as temporary:
            pending = Path(temporary) / args.key
            pending.mkdir()
            LOG.info(
                "Installing %s: decompressing tar.xz once to %s", args.key, destination
            )
            with tarfile.open(archive, "r:xz") as bundled:
                bundled.extractall(pending, filter="data")
            pending.rename(destination)
        LOG.info("Installed example prepared: %s", destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lgs", type=Path)
    parser.add_argument("--fields", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--description", required=True)
    parser.add_argument("--frames", type=int, required=True)
    parser.add_argument("--install-root", type=Path)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    package(args)


if __name__ == "__main__":
    main()
