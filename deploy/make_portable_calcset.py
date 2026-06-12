#!/usr/bin/env python3
"""make_portable_calcset.py — assemble a self-contained, relocatable calcset.

A trajectory calcset's .LGS references external paths the viewer does not bundle:
the raw MD dir + topology (provenance only) and the DFT job tree of ORCA `.out`
files (the shielding overlay). On another machine those paths do not exist, so
the calcset will not open as-is. This tool produces ONE directory that opens
anywhere:

  * all viewing data (trajectory.h5, sidecars, per-frame npys, pdbs) copied in;
  * each DFT frame's meta.json AND its `out_primary` ORCA .out copied into
    dft/<job>/ so the shielding overlay works offline;
  * topology_top bundled if it resolves; md_dir left as provenance (the reader's
    loader treats md_dir/topology_top as optional — see CalcsetManifest.cpp);
  * every path in the .LGS rewritten relative to the bundle.

Usage:
    make_portable_calcset.py --src <calcset_dir> --out <bundle_dir> [--no-dft]

The result is portable: copy it to the target box (or bind it into the
h5reader Singularity image) and open the .LGS inside it.
"""
import argparse
import json
import os
import shutil
import sys
from pathlib import Path


def find_lgs(src: Path) -> Path:
    lgs = sorted(src.glob("*.lgs")) + sorted(src.glob("*.LGS"))
    if len(lgs) != 1:
        sys.exit(f"expected exactly one .LGS in {src}, found {len(lgs)}: {lgs}")
    return lgs[0]


def resolve(base: Path, p: str) -> Path:
    return Path(p) if os.path.isabs(p) else (base / p).resolve()


def main() -> None:
    ap = argparse.ArgumentParser(description="Assemble a self-contained portable calcset.")
    ap.add_argument("--src", required=True, type=Path, help="source calcset directory")
    ap.add_argument("--out", required=True, type=Path, help="output bundle directory (must be empty/new)")
    ap.add_argument("--no-dft", action="store_true", help="omit the DFT .out overlay (smaller bundle)")
    args = ap.parse_args()

    src = args.src.resolve()
    out = args.out.resolve()
    if not src.is_dir():
        sys.exit(f"--src is not a directory: {src}")
    if out.exists() and any(out.iterdir()):
        sys.exit(f"--out exists and is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)

    lgs_path = find_lgs(src)
    manifest = json.loads(lgs_path.read_text())

    # 1. Copy all self-contained viewing data (everything except the .LGS, which
    #    we rewrite). extraction_dir="." / trajectory_h5 / extraction_manifest are
    #    already relative, so they keep resolving inside the bundle unchanged.
    print(f"[1/4] copying viewing data {src} -> {out}  (the bulk; minutes for the npys)")
    for entry in sorted(src.iterdir()):
        if entry.suffix.lower() == ".lgs":
            continue
        dst = out / entry.name
        if entry.is_dir():
            shutil.copytree(entry, dst, symlinks=False)
        else:
            shutil.copy2(entry, dst)
        print(f"      + {entry.name}")

    traj = manifest.get("trajectory", {})

    # 2. topology_top: bundle if it resolves (small); md_dir stays provenance.
    top = traj.get("topology_top")
    if top:
        top_abs = resolve(src, top)
        if top_abs.is_file():
            shutil.copy2(top_abs, out / top_abs.name)
            traj["topology_top"] = top_abs.name
            print(f"[2/4] bundled topology_top -> {top_abs.name}")
        else:
            print(f"[2/4] topology_top not found ({top_abs}); left as-is (non-fatal at load)")
    if traj.get("md_dir"):
        print(f"      md_dir left as provenance (raw MD not bundled; non-fatal at load)")

    # 3. DFT overlay: each frame's meta.json + its out_primary .out -> dft/<job>/.
    dft = manifest.get("dft")
    if dft and not args.no_dft:
        frames = dft.get("frames", [])
        kept = skipped = 0
        print(f"[3/4] bundling {len(frames)} DFT frames (meta.json + ORCA .out)")
        for fr in frames:
            mj = fr.get("meta_json")
            if not mj:
                skipped += 1
                continue
            mj_abs = resolve(src, mj)
            if not mj_abs.is_file():
                print(f"      ! meta.json missing, skipping frame {fr.get('frame_index')}: {mj_abs}")
                skipped += 1
                continue
            job = mj_abs.parent.name
            jobdst = out / "dft" / job
            jobdst.mkdir(parents=True, exist_ok=True)
            shutil.copy2(mj_abs, jobdst / mj_abs.name)
            meta = json.loads(mj_abs.read_text())
            out_primary = meta.get("files", {}).get("out_primary")
            if out_primary:
                out_abs = mj_abs.parent / out_primary
                if out_abs.is_file():
                    shutil.copy2(out_abs, jobdst / out_primary)
                else:
                    print(f"      ! out_primary missing for {job}: {out_abs}")
            fr["meta_json"] = os.path.join("dft", job, mj_abs.name)
            kept += 1
        print(f"      DFT: {kept} bundled, {skipped} skipped")
    elif args.no_dft and "dft" in manifest:
        del manifest["dft"]
        print("[3/4] DFT overlay omitted (--no-dft)")

    # 4. Write the rewritten, fully-relative .LGS into the bundle.
    (out / lgs_path.name).write_text(json.dumps(manifest, indent=2))
    total_gb = sum(f.stat().st_size for f in out.rglob("*") if f.is_file()) / 1e9
    print(f"[4/4] wrote portable .LGS -> {out / lgs_path.name}")
    print(f"\nBundle ready: {out}  ({total_gb:.1f} GB)")
    print(f"Open with:    apptainer run --bind {out.parent} h5reader.sif {out / lgs_path.name}")


if __name__ == "__main__":
    main()
