"""Write a Reader collection LGS for packages on a provenance drive."""

import argparse
import json
from pathlib import Path


SMALL_TITLES = {
    "1p9j": "1P9J",
    "chignolin": "Chignolin",
    "trpcage": "Trp-cage",
}


def one_lgs(run):
    files = [path for path in run.iterdir() if path.is_file() and path.suffix.lower() == ".lgs"]
    if len(files) != 1:
        raise ValueError(f"Expected one .lgs in {run}, found {len(files)}")
    return files[0]


def source_runs(root):
    datasets = root / "datasets"
    for run in sorted((datasets / "Small_proteins").iterdir()):
        if run.is_dir():
            lgs = one_lgs(run)
            source = json.loads(lgs.read_text(encoding="utf-8-sig"))
            short_name = run.name.split("_", 1)[0].lower()
            yield run, SMALL_TITLES.get(short_name, source["human_name"]), "Small proteins"
    md_runs = (run for run in (datasets / "MD_173").iterdir()
               if run.is_dir() and run.name.startswith("bmr"))
    for run in sorted(md_runs, key=lambda path: int(path.name[3:])):
        one_lgs(run)
        yield run, f"BMRB {run.name[3:]}", "MD trajectories"


def make_entries(reader_root):
    entries = []
    for run, title, group in source_runs(reader_root):
        key = f"{run.name.split('_', 1)[0].lower()}-reader-v1"
        metadata_path = reader_root / "packages" / f"{key}.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        archive_name = f"{key}.tar.xz"
        archive = metadata_path.parent / archive_name
        if (metadata["key"] != key or metadata.get("archive") != archive_name
                or metadata["entry_point"] != "run.LGS"
                or metadata["archive_bytes"] != archive.stat().st_size
                or metadata["expanded_bytes"] < 1 or metadata["frames"] < 1):
            raise ValueError(f"Package metadata does not match {archive}")
        entries.append({
            "key": key,
            "title": title,
            "group": group,
            "description": metadata["description"],
            "frames": metadata["frames"],
            "archive": archive.relative_to(reader_root).as_posix(),
            "entry_point": "run.LGS",
            "archive_bytes": metadata["archive_bytes"],
            "expanded_bytes": metadata["expanded_bytes"],
        })
    if not entries:
        raise ValueError(f"No Reader runs found under {reader_root / 'datasets'}")
    return entries


def packaged_entries(reader_root, descriptions):
    packages = {}
    for path in sorted((reader_root / "packages").iterdir()):
        if not path.is_file() or path.suffix != ".json":
            continue
        metadata = json.loads(path.read_text(encoding="utf-8-sig"))
        key = metadata["key"]
        archive_name = f"{key}.tar.xz"
        archive = path.parent / archive_name
        if (path.stem != key or metadata["archive"] != archive_name
                or metadata["entry_point"] != "run.LGS"
                or metadata["archive_bytes"] != archive.stat().st_size
                or metadata["expanded_bytes"] < 1 or metadata["frames"] < 1):
            raise ValueError(f"Package metadata does not match {archive}")
        member = key.split("-", 1)[0]
        if member in packages:
            raise ValueError(f"More than one package for {member}")
        packages[member] = metadata

    entries = []
    for description in descriptions["entries"]:
        member = description["key"].split("-", 1)[0]
        if member not in packages:
            raise ValueError(f"No unambiguous package for catalog entry {member}")
        metadata = packages.pop(member)
        entry = dict(description)
        for field in ("key", "entry_point", "frames", "archive_bytes", "expanded_bytes"):
            entry[field] = metadata[field]
        entry["archive"] = "packages/" + metadata["archive"]
        entries.append(entry)
    if packages or not entries:
        raise ValueError(f"Catalog does not cover the packages: {', '.join(sorted(packages))}")
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reader_root", type=Path, help="Directory containing packages")
    parser.add_argument("output", type=Path, help="Collection .lgs to write")
    parser.add_argument("--descriptions", type=Path,
                        help="Existing catalog with deposited names; no datasets directory is needed")
    args = parser.parse_args()
    root = args.reader_root.resolve()
    output = args.output.resolve()
    if output.parent != root:
        parser.error("the collection .lgs must sit beside the packages directory")
    if args.descriptions:
        descriptions = json.loads(args.descriptions.read_text(encoding="utf-8-sig"))
        entries = packaged_entries(root, descriptions)
    else:
        entries = make_entries(root)
        if (sum(entry["group"] == "Small proteins" for entry in entries) != 3
                or sum(entry["group"] == "MD trajectories" for entry in entries) != 173):
            raise ValueError("Expected all 3 small proteins and 173 MD trajectories")
    document = {
        "schema_version": 1,
        "kind": "collection",
        "title": "Reader datasets",
        "entries": entries,
    }
    with output.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(document, indent=2, ensure_ascii=True) + "\n")
    print(f"Wrote {len(entries)} packaged runs to {output}")


if __name__ == "__main__":
    main()
