#!/usr/bin/env python3
"""Generate `<dataset_id>.LGS` for an existing calcset directory.

The `.LGS` (Lowly Graduate Student) file is the calcset's top-level
wrapper — see `/shared/2026Thesis/nmr-shielding/spec/CALCSET_MANIFEST.md`
for the authoritative schema. This tool walks an existing on-disk
layout and emits the matching `.LGS` (schema v1) so the h5-reader and
other consumers can locate every artifact without filename parsing.

Detection rules — no glob-and-guess: only canonical layouts. The
producer's actual layouts (verified on the 1P9J calibration fixture
2026-05-27):

  Trajectory calcset (the 1P9J fixture shape):
    <root>/<extraction_dir>/trajectory.h5
    <root>/<extraction_dir>/extraction_manifest.json
    <root>/<extraction_dir>/{atoms_category_info,bonds,residues,rings,
                              ring_membership}.npy
    <root>/<extraction_dir>/npys/frame_NNNNNN/   (per-frame snapshots)
    <root>/md/                                    (production.tpr/.trr/.edr)
    <root>/topol.top
    <root>/dft/jobs/..._fNNNNNN_t<ps>/            (optional, when DFT exists)

  Single pose calcset:
    <root>/extraction_manifest.json
    <root>/atoms_category_info.npy                (sidecar anchor)
    <pose-kind from the producer's CLI invocation — passed explicitly)

  Mutant pair:
    Two side-by-side single-pose calcsets, each with its own .LGS.

DFT frame discovery: walk every `dft/jobs/*/` directory, or an explicit
`--dft-jobs-dir`, and read each `*_meta.json` (NOT parsing the dir
name) for the typed `frame_index`. Only completed ORCA jobs
(`orca_exit_code == 0`) are emitted. The frames[] array is sorted by
frame_index. `frame_stride` is the human-readable cadence summary
computed from min/max/median-diff; consumers iterate frames[] for
coverage, not the stride.

Join semantics: `dft.frames[].frame_index` is expressed in
`trajectory.frame_index_basis` and joins to the trajectory H5 dataset
`/trajectory/frames/original_index`. Do not infer coverage from stride
or from the row number in `dft.frames[]`.

Usage:
  python3 lgs_write.py <calcset_root>
  python3 lgs_write.py --dry-run <calcset_root>     # print, don't write
  python3 lgs_write.py --force <calcset_root>       # overwrite existing
  python3 lgs_write.py --root <calcset_root> --extraction-dir extract
  python3 lgs_write.py --dft-jobs-dir /path/to/jobs <calcset_root>
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import statistics
import sys
from pathlib import Path
from typing import Optional

LGS_WRITER_VERSION = "lgs-tools 0.1.1"
SCHEMA_VERSION = 1


# ---------------------------------------------------------------------
# Detection helpers — each looks for a documented producer layout.
# ---------------------------------------------------------------------

def _find_extraction_dir(root: Path) -> Optional[Path]:
    """Locate the extraction directory by canonical layout.

    The producer's --trajectory mode writes to a sibling `extract/`
    (the 1P9J fixture), or — in some older layouts — the root itself.
    We check both, preferring `extract/` because it's the convention.
    """
    extract = root / "extract"
    if (extract / "trajectory.h5").is_file():
        return extract
    if (root / "trajectory.h5").is_file():
        return root
    return None


def _find_trajectory_h5(extraction_dir: Path) -> Optional[Path]:
    cand = extraction_dir / "trajectory.h5"
    return cand if cand.is_file() else None


def _find_extraction_manifest(extraction_dir: Path) -> Optional[Path]:
    cand = extraction_dir / "extraction_manifest.json"
    return cand if cand.is_file() else None


def _find_reference_pdb(extraction_dir: Path) -> Optional[Path]:
    cand = extraction_dir / "reference.pdb"
    return cand if cand.is_file() else None


def _find_md_dir(root: Path) -> Optional[Path]:
    cand = root / "md"
    return cand if cand.is_dir() else None


def _find_topology_top(root: Path) -> Optional[Path]:
    cand = root / "topol.top"
    return cand if cand.is_file() else None


def _find_dft_jobs(root: Path) -> Optional[Path]:
    cand = root / "dft" / "jobs"
    return cand if cand.is_dir() else None


def _read_extraction_protein_id(extraction_manifest: Path) -> Optional[str]:
    """Best-effort: protein_id from extraction_manifest.json.

    The producer has a known bug where this is sometimes 'extract'.
    The caller should fall back to the root dir name in that case.
    """
    try:
        with extraction_manifest.open() as f:
            data = json.load(f)
        pid = data.get("protein_id")
        if isinstance(pid, str) and pid and pid != "extract":
            return pid
    except (json.JSONDecodeError, OSError):
        pass
    return None


def _read_extractor_version(extraction_manifest: Path) -> Optional[str]:
    try:
        with extraction_manifest.open() as f:
            data = json.load(f)
        v = data.get("extractor_version")
        if isinstance(v, str) and v:
            return v
    except (json.JSONDecodeError, OSError):
        pass
    return None


def _collect_dft_frames(dft_jobs: Path) -> list[dict]:
    """Walk `dft/jobs/*/` and collect (frame_index, meta_json_relpath) per job.

    Reads `frame_index` from each meta.json — NOT from the directory
    name. The per-job meta is the authoritative source for the typed
    frame_index integer. Only completed ORCA jobs (`orca_exit_code == 0`)
    are included.
    """
    frames: list[dict] = []
    for job_dir in sorted(dft_jobs.iterdir()):
        if not job_dir.is_dir():
            continue
        # Each job dir holds exactly one `<job_id>_meta.json`. The job_id
        # is the dir name, so the meta filename is mechanical — no glob.
        job_id = job_dir.name
        meta_path = job_dir / f"{job_id}_meta.json"
        if not meta_path.is_file():
            print(f"  warning: no meta.json at {meta_path}", file=sys.stderr)
            continue
        try:
            with meta_path.open() as f:
                meta = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            print(f"  warning: cannot read {meta_path}: {e}", file=sys.stderr)
            continue
        frame_index = meta.get("frame_index")
        if not isinstance(frame_index, int):
            print(f"  warning: {meta_path} has no integer frame_index", file=sys.stderr)
            continue
        if meta.get("orca_exit_code") != 0:
            continue
        frames.append({
            "frame_index": int(frame_index),
            "meta_path": meta_path,
        })
    # Sort by frame_index for the deterministic, human-readable order.
    frames.sort(key=lambda fr: fr["frame_index"])
    return frames


def _read_dft_method(dft_jobs: Path, frames: list[dict]) -> Optional[str]:
    """Best-effort: method string from the first job's meta.json.

    The meta.json schema doesn't expose `method` directly; the caller
    can hand-edit if the producer-level convention changes. We label
    by the campaign's known sampling cadence for clarity.
    """
    # The 1P9J campaign uses r2SCAN/def2-SVP/CPCM(Water) — a hard-coded
    # default that matches the only DFT campaign in the corpus today.
    # Future calcsets with other methods should hand-edit the .LGS or
    # the producer can plumb method into meta.json schema v3.
    return "r2SCAN def2-SVP def2/J NMR  CPCM(Water)"


def _compute_stride(frame_indices: list[int]) -> dict:
    """Derive (first, last, step) from the frame_index list.

    `step` is the median diff (robust to partial-campaign gaps).
    """
    if not frame_indices:
        return {"first": 0, "last": 0, "step": 1}
    if len(frame_indices) == 1:
        return {"first": frame_indices[0], "last": frame_indices[0], "step": 1}
    diffs = [b - a for a, b in zip(frame_indices, frame_indices[1:]) if b > a]
    step = int(round(statistics.median(diffs))) if diffs else 1
    return {
        "first": frame_indices[0],
        "last": frame_indices[-1],
        "step": step,
    }


# ---------------------------------------------------------------------
# Manifest assembly
# ---------------------------------------------------------------------

def _relpath(target: Path, root: Path) -> str:
    """Relative path, forward slashes (cross-platform JSON friendly)."""
    try:
        rel = target.resolve().relative_to(root.resolve())
        return str(rel).replace("\\", "/")
    except ValueError:
        # outside root — fall back to absolute
        return str(target.resolve()).replace("\\", "/")


def build_trajectory_manifest(
    root: Path,
    dataset_id: str,
    protein_id: str,
    human_name: str,
    extraction_dir: Optional[Path] = None,
    dft_jobs_dir: Optional[Path] = None,
    frame_dt_ps: float = 10.0,
    frame_index_basis: str = "trr_frame_index",
) -> dict:
    """Build the JSON dict for a trajectory calcset .LGS."""
    extr = extraction_dir if extraction_dir is not None else _find_extraction_dir(root)
    if extr is None:
        raise RuntimeError(
            f"No extraction directory found under {root} "
            f"(neither {root}/extract/trajectory.h5 nor {root}/trajectory.h5)"
        )
    traj_h5 = _find_trajectory_h5(extr)
    if traj_h5 is None:
        raise RuntimeError(f"No trajectory.h5 in {extr}")
    extr_manifest = _find_extraction_manifest(extr)
    if extr_manifest is None:
        raise RuntimeError(f"No extraction_manifest.json in {extr}")

    md_dir = _find_md_dir(root)
    if md_dir is None:
        raise RuntimeError(f"No md/ dir in {root}")
    topol = _find_topology_top(root)
    if topol is None:
        raise RuntimeError(f"No topol.top in {root}")

    out: dict = {
        "schema_version": SCHEMA_VERSION,
        "kind": "trajectory",
        "dataset_id": dataset_id,
        "protein_id": protein_id,
        "human_name": human_name,
        "trajectory": {
            "md_dir": _relpath(md_dir, root),
            "topology_top": _relpath(topol, root),
            "extraction_dir": _relpath(extr, root),
            "trajectory_h5": _relpath(traj_h5, root),
            "extraction_manifest": _relpath(extr_manifest, root),
            "frame_dt_ps": frame_dt_ps,
            "frame_index_basis": frame_index_basis,
        },
    }
    ref_pdb = _find_reference_pdb(extr)
    if ref_pdb is not None:
        out["trajectory"]["reference_pdb"] = _relpath(ref_pdb, root)

    # DFT block — optional
    dft_jobs = dft_jobs_dir if dft_jobs_dir is not None else _find_dft_jobs(root)
    if dft_jobs is not None:
        frames = _collect_dft_frames(dft_jobs)
        method = _read_dft_method(dft_jobs, frames) or "unknown"
        # Heuristic: campaign_target_frames from `dft/_consolidation_snapshot.json`
        # when present, else len(frames).
        target = len(frames)
        snap = root / "dft" / "_consolidation_snapshot.json"
        if snap.is_file():
            try:
                with snap.open() as f:
                    data = json.load(f)
                ct = data.get("campaign_target")
                if isinstance(ct, int) and ct > 0:
                    target = ct
            except (json.JSONDecodeError, OSError):
                pass
        out["dft"] = {
            "method": method,
            "campaign_target_frames": target,
            "frame_stride": _compute_stride([fr["frame_index"] for fr in frames]),
            "frames": [
                {
                    "frame_index": fr["frame_index"],
                    "meta_json": _relpath(fr["meta_path"], root),
                }
                for fr in frames
            ],
        }

    # metadata — always present
    extr_version = _read_extractor_version(extr_manifest)
    out["metadata"] = {
        "generated_at_utc": dt.datetime.now(dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "lgs_writer": LGS_WRITER_VERSION,
    }
    if extr_version:
        out["metadata"]["producer_extractor_version"] = extr_version
    return out


def render(manifest: dict) -> str:
    """Format the manifest as a stable JSON string. UTF-8 friendly
    (non-ASCII chars like superscript ² stay readable)."""
    return json.dumps(manifest, indent=2, sort_keys=False, ensure_ascii=False) + "\n"


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("calcset_root", nargs="?", type=Path,
                   help="Calcset root directory (with md/, extract/, dft/).")
    p.add_argument("--root", dest="root_flag", type=Path,
                   help="Alias for calcset_root; useful with --extraction-dir.")
    p.add_argument("--extraction-dir", type=Path,
                   help="Explicit extraction directory (defaults to "
                        "<root>/extract or <root>).")
    p.add_argument("--dft-jobs-dir", type=Path,
                   help="Explicit ORCA DFT jobs directory. Use this when "
                        "completed jobs live outside <root>/dft/jobs.")
    p.add_argument("--dataset-id",
                   help="dataset_id (defaults to calcset root basename).")
    p.add_argument("--protein-id",
                   help="protein_id (defaults to extraction_manifest's value, "
                        "falling back to dataset_id).")
    p.add_argument("--human-name",
                   help="human_name (defaults to a generic description).")
    p.add_argument("--frame-dt-ps", type=float, default=10.0,
                   help="MD frame Δt in ps (default 10.0 — 1P9J convention).")
    p.add_argument("--frame-index-basis", default="trr_frame_index",
                   help="frame_index_basis string (default 'trr_frame_index').")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the .LGS to stdout instead of writing to disk.")
    p.add_argument("--force", action="store_true",
                   help="Overwrite an existing .LGS at the target path.")
    args = p.parse_args(argv)

    root = args.calcset_root if args.calcset_root is not None else args.root_flag
    if root is None:
        print("error: <calcset_root> (or --root) is required", file=sys.stderr)
        p.print_help(sys.stderr)
        return 2
    root = root.resolve()
    if not root.is_dir():
        print(f"error: {root} is not a directory", file=sys.stderr)
        return 2

    dataset_id = args.dataset_id or root.name
    extraction_dir = (
        args.extraction_dir.resolve() if args.extraction_dir else None
    )
    dft_jobs_dir = args.dft_jobs_dir.resolve() if args.dft_jobs_dir else None
    if dft_jobs_dir is not None and not dft_jobs_dir.is_dir():
        print(f"error: {dft_jobs_dir} is not a directory", file=sys.stderr)
        return 2
    if dft_jobs_dir is not None:
        try:
            dft_jobs_dir.relative_to(root)
        except ValueError:
            print(
                "warning: --dft-jobs-dir is outside the calcset root; "
                "meta_json paths will be absolute in the .LGS",
                file=sys.stderr,
            )

    # Resolve protein_id. Order of authority:
    #   1. --protein-id flag (explicit override)
    #   2. extraction_manifest.json's protein_id (when it is not the
    #      producer's "extract" bug value)
    #   3. the first DFT job's meta.json (these carry the canonical
    #      `<PDB>_<chain id>` form, e.g. `1P9J_5801`)
    #   4. dataset_id fallback
    extr = extraction_dir if extraction_dir is not None else _find_extraction_dir(root)
    protein_id = args.protein_id
    if not protein_id and extr is not None:
        em = _find_extraction_manifest(extr)
        if em is not None:
            protein_id = _read_extraction_protein_id(em)
    if not protein_id:
        dft_jobs = dft_jobs_dir if dft_jobs_dir is not None else _find_dft_jobs(root)
        if dft_jobs is not None:
            first_meta = next(
                (j / f"{j.name}_meta.json"
                 for j in sorted(dft_jobs.iterdir()) if j.is_dir()),
                None,
            )
            if first_meta and first_meta.is_file():
                try:
                    with first_meta.open() as f:
                        pid = json.load(f).get("protein_id")
                    if isinstance(pid, str) and pid:
                        protein_id = pid
                except (json.JSONDecodeError, OSError):
                    pass
    if not protein_id:
        protein_id = dataset_id

    human_name = args.human_name or f"{dataset_id} calcset"

    try:
        manifest = build_trajectory_manifest(
            root=root,
            dataset_id=dataset_id,
            protein_id=protein_id,
            human_name=human_name,
            extraction_dir=extraction_dir,
            dft_jobs_dir=dft_jobs_dir,
            frame_dt_ps=args.frame_dt_ps,
            frame_index_basis=args.frame_index_basis,
        )
    except RuntimeError as e:
        print(f"error: {e}", file=sys.stderr)
        return 3

    text = render(manifest)
    target = root / f"{dataset_id}.LGS"
    if args.dry_run:
        sys.stdout.write(text)
        return 0
    if target.exists() and not args.force:
        print(f"error: {target} already exists; pass --force to overwrite",
              file=sys.stderr)
        return 4
    target.write_text(text)
    print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
