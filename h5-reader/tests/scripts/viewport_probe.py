#!/usr/bin/env python3
"""Viewport probe — establish the marker-drift floor under plane lock
and the bug magnitude without it.

Drives experiments against a running reader started with
``--rest <port>``. The reader is launched separately (e.g. by
``HARNESS_BASELINE_2026-05-30.md``'s recipe); this script only talks
to the REST surface.

Strategy (per notes/VIEWPORT_OBSERVATIONS_2026-05-30.md §5b):

* Pixel-perfect snapshot comparison is doomed for our renderer (FXAA,
  imposter shaders, GPU driver edge noise). We use centroid-of-blob
  analysis instead: ``POST /selection/instrument`` switches the
  ``MeasurementOverlay`` spheres to a CPK-distinct fixed palette
  (magenta / spring green / deep pink / vivid violet) at opacity 1.0
  and radius 1.5 Å. The marker hue falls outside every CPK element
  colour, so a hue threshold isolates the marker against the rendered
  molecule. ``scipy.ndimage.label`` extracts the largest connected
  component and ``scipy.ndimage.center_of_mass`` returns the sub-pixel
  centroid.

* Drift metric: per-frame centroid distance from the frame-0 centroid.
  Mean, std, max, fraction-within-5px, fraction-within-10px.

* Two experiments per run:

    1. **with_lock** — plane lock enabled on the 3 atoms. The marker
       (slot 0, magenta) tracks atom 0 of the selection. Drift here is
       the atom-vibration floor: how much the locked-camera-stable
       view of a thermal-vibrating atom moves on screen.

    2. **no_lock** — plane lock disabled, same 3 atoms still in the
       selection. Drift here is floor + the centroid-delta camera bug.
       If (mean-no-lock) >> (mean-with-lock), the bug is measurable
       and the next session's refactor has a baseline to beat.

Outputs:

* Per-frame PNG dumps to ``<out-dir>/with_lock/`` and
  ``<out-dir>/no_lock/`` for visual review.
* Per-frame CSV (frame, centroid_x, centroid_y, area, hash) to
  ``<out-dir>/with_lock.csv`` and ``<out-dir>/no_lock.csv``.
* Summary JSON to stdout.

Usage:

    python tests/scripts/viewport_probe.py \\
        --base http://127.0.0.1:9988 \\
        --atoms 12 13 14 \\
        --frames 200 \\
        --out-dir /tmp/viewport_probe_run1

The reader binary is NOT launched by this script; start it separately
with ``./h5reader --rest <port> <fixture-dir>``. The script verifies
``/health`` before each experiment.

Marker hue thresholds (slot 0, magenta) — derived empirically from
a baseline snapshot on the 1P9J fixture with FXAA on:

    H in [285°, 315°], S >= 0.6, V >= 0.6

Adjust if the GPU/driver shifts the rendered hue (e.g. ANGLE on
Windows desaturates differently from Mesa llvmpipe).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import requests
from PIL import Image
from scipy import ndimage


# ----------------------------------------------------------------------------
# Marker hue presets — one entry per AtomSelection slot. The instrument-mode
# palette is hard-coded in MeasurementOverlay.cpp (kInstrumentRgb); keep these
# in sync with that table or the harness will not find the blob.
# ----------------------------------------------------------------------------

@dataclass(frozen=True)
class MarkerHue:
    name: str
    hue_lo_deg: float
    hue_hi_deg: float
    sat_min: float = 0.6
    val_min: float = 0.6


MARKER_HUES: dict[int, MarkerHue] = {
    0: MarkerHue("magenta",      hue_lo_deg=285.0, hue_hi_deg=315.0),
    1: MarkerHue("spring green", hue_lo_deg=130.0, hue_hi_deg=170.0),
    2: MarkerHue("deep pink",    hue_lo_deg=320.0, hue_hi_deg=350.0),
    3: MarkerHue("vivid violet", hue_lo_deg=270.0, hue_hi_deg=285.0),
}


# ----------------------------------------------------------------------------
# REST helpers — thin wrappers, no retry. The reader is a long-running test
# fixture; if a call times out something is wrong and we want a hard failure.
# ----------------------------------------------------------------------------

def _post(base: str, path: str, body: dict | None = None, timeout: float = 10.0) -> requests.Response:
    r = requests.post(f"{base}{path}", json=body or {}, timeout=timeout)
    r.raise_for_status()
    return r


def _get(base: str, path: str, timeout: float = 5.0) -> requests.Response:
    r = requests.get(f"{base}{path}", timeout=timeout)
    r.raise_for_status()
    return r


def health_check(base: str) -> None:
    r = _get(base, "/health")
    j = r.json()
    if not j.get("ok"):
        raise RuntimeError(f"reader health check failed: {j}")


def atom_count(base: str) -> int:
    return int(_get(base, "/protein/atoms").json()["count"])


def set_atoms(base: str, atoms: Iterable[int]) -> None:
    _post(base, "/selection/atoms", {"atoms": list(atoms)})


def clear_selection(base: str) -> None:
    _post(base, "/selection/clear")


def set_instrument(base: str, enabled: bool, focus_only: bool = False) -> None:
    _post(base, "/selection/instrument",
          {"enabled": enabled, "focus_only": focus_only})


def enable_plane_lock(base: str, atoms: Iterable[int]) -> None:
    _post(base, "/plane-lock/enable", {"atoms": list(atoms)})


def disable_plane_lock(base: str) -> None:
    _post(base, "/plane-lock/disable")


def set_docks_visible(base: str, visible: bool) -> None:
    """POST /docks/visible — hides or restores the dock widgets to give the
    central VTK viewport more pixels. The reader stashes pre-hide
    visibility so restore is safe even if a dock was already hidden."""
    _post(base, "/docks/visible", {"visible": bool(visible)})


def set_transform(base: str, kind: str, **kwargs) -> None:
    """POST /transform — sets the wrapped Conformation's transform mode.

    kind: "identity" | "center_com" | "fit_reference" | "fit_subset"

    kwargs are forwarded as JSON fields:
      reference_frame: int (default 0; used by fit_reference / fit_subset)
      subset_atoms:    [int, ...] (fit_subset only)
      backbone_only:   bool (fit_subset shorthand — selects backbone via
                       QtAtom::IsBackbone, no atom-name string parsing)
    """
    body: dict = {"kind": kind}
    body.update(kwargs)
    _post(base, "/transform", body)


def get_transform(base: str) -> dict:
    return _get(base, "/transform").json()


def set_camera_mode(base: str, mode: str, **kwargs) -> None:
    """POST /camera/mode — sets the CameraComposer mode.

    mode: "free" | "atom" | "bond" | "dihedral" | "plane" | "subset"

    kwargs are forwarded as JSON fields. Examples:
      set_camera_mode(base, "atom", atom=14)
      set_camera_mode(base, "bond", a=12, b=13)
      set_camera_mode(base, "dihedral", a=12, b=13, c=14, d=15)
      set_camera_mode(base, "plane", a=12, b=13, c=14)
      set_camera_mode(base, "subset", backbone_only=True)
      set_camera_mode(base, "subset", atoms=[12, 13, 14, 15, 16])
    """
    body: dict = {"mode": mode}
    body.update(kwargs)
    _post(base, "/camera/mode", body)


def clear_camera_mode(base: str) -> None:
    _post(base, "/camera/clear")


def get_camera_mode(base: str) -> dict:
    return _get(base, "/camera/mode").json()


def set_camera_focus_atom(base: str, atom: int, kind: str = "plane") -> None:
    """POST /camera/focus_atom — derive a typed CameraMode + policy from a
    focus atom by reaching into its residue's typed backbone-atom-index
    cache.

    kind: "plane" | "dihedral" | "dihedral_phi" | "dihedral_psi"

    Plane is the canonical "focus atom + local neighborhood coherent"
    recipe — a 3-atom plane lock on the residue's N/CA/C backbone.
    Median drift should sit at the plane-lock floor (≤5 px).
    """
    _post(base, "/camera/focus_atom", {"atom": int(atom), "kind": kind})


def set_log_mask(base: str, mask: int) -> dict:
    """POST /log/mask — bitmask gate for the structured logger."""
    return _post(base, "/log/mask", {"mask": int(mask)}).json()


def set_frame(base: str, frame: int) -> None:
    _post(base, "/frame/set", {"frame": int(frame)})


def frame_count(base: str) -> int:
    return int(_get(base, "/frame/current").json()["count"])


def screenshot(base: str, force_render: bool = True) -> bytes:
    r = requests.post(
        f"{base}/screenshot",
        json={"target": "scene", "force_render": bool(force_render)},
        timeout=15.0,
    )
    r.raise_for_status()
    return r.content


# ----------------------------------------------------------------------------
# Blob detection — RGB → HSV → hue/sat/val threshold → largest connected
# component → sub-pixel centroid via scipy.ndimage.center_of_mass.
# ----------------------------------------------------------------------------

def _rgb_to_hsv(rgb: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised RGB→HSV, RGB in [0,1]. Returns (H deg, S, V) each in
    the same (H, W) shape as the input's first two axes."""
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    maxc = np.maximum(np.maximum(r, g), b)
    minc = np.minimum(np.minimum(r, g), b)
    v = maxc
    delta = maxc - minc
    s = np.where(maxc > 0, delta / np.where(maxc == 0, 1.0, maxc), 0.0)
    # avoid div-by-zero on perfectly grey pixels by guarding delta in the
    # denominator; the np.select branch already masks those out.
    safe_delta = np.where(delta == 0, 1.0, delta)
    rc = (maxc - r) / safe_delta
    gc = (maxc - g) / safe_delta
    bc = (maxc - b) / safe_delta
    h = np.select(
        [maxc == r, maxc == g, maxc == b],
        [bc - gc, 2.0 + rc - bc, 4.0 + gc - rc],
        default=0.0,
    )
    h = (h / 6.0) % 1.0
    return h * 360.0, s, v


def find_marker_centroid(
    png_bytes: bytes,
    hue: MarkerHue,
    min_area_px: int = 5,
) -> tuple[float, float, int] | None:
    """Returns (x, y, area_pixels) of the largest matching blob, or None.

    Hue mask handles wrap-around: ``hue_hi_deg < hue_lo_deg`` is treated
    as ``[lo, 360) ∪ [0, hi]``.

    Sub-pixel accurate (``center_of_mass`` returns floats). Y-axis is
    image-space (top-down); convert to GL-style if needed downstream.
    """
    img = np.asarray(Image.open(io.BytesIO(png_bytes)).convert("RGB"),
                     dtype=np.float32) / 255.0
    h_deg, s, v = _rgb_to_hsv(img)
    if hue.hue_hi_deg < hue.hue_lo_deg:
        hue_mask = (h_deg >= hue.hue_lo_deg) | (h_deg <= hue.hue_hi_deg)
    else:
        hue_mask = (h_deg >= hue.hue_lo_deg) & (h_deg <= hue.hue_hi_deg)
    mask = hue_mask & (s >= hue.sat_min) & (v >= hue.val_min)
    if not mask.any():
        return None
    labels, n = ndimage.label(mask)
    if n == 0:
        return None
    # Largest component
    sizes = ndimage.sum(mask, labels, index=range(1, n + 1)).astype(np.int64)
    largest = int(np.argmax(sizes)) + 1
    area = int(sizes[largest - 1])
    if area < min_area_px:
        return None
    cy, cx = ndimage.center_of_mass(mask, labels, largest)
    return float(cx), float(cy), area


# ----------------------------------------------------------------------------
# Experiment driver — sweeps frames, captures, finds marker, accumulates
# per-frame records. PNG dumps are kept for visual review.
# ----------------------------------------------------------------------------

@dataclass
class FrameRecord:
    frame: int
    centroid_x: float | None
    centroid_y: float | None
    area_px: int | None
    png_sha1: str


def _png_sha1(b: bytes) -> str:
    return hashlib.sha1(b).hexdigest()[:12]


def run_drift_experiment(
    base: str,
    atoms: list[int],
    *,
    plane_lock: bool,
    transform_kind: str = "identity",
    transform_kwargs: dict | None = None,
    frames: list[int],
    hue: MarkerHue,
    out_dir: Path,
    keep_pngs: bool = True,
    focus_only: bool = True,
    camera_mode: str | None = None,
    camera_mode_kwargs: dict | None = None,
    focus_atom: int | None = None,
    focus_atom_kind: str = "plane",
) -> list[FrameRecord]:
    """Configure the reader for the experiment, sweep frames, return
    per-frame records.

    Setup order (matters):
      1. clear selection
      2. bulk-set the atoms
      3. set transform mode (identity / center_com / fit_reference / fit_subset)
      4. instrument mode ON with focus_only=True by default — only the
         focus-slot magenta sphere renders, eliminating the slot-1
         eclipses-slot-0 problem the prior no-lock baseline run hit
      5. one of (precedence top-down):
           a. focus_atom is set → POST /camera/focus_atom (derives a
              typed CameraMode from the focus atom's residue backbone)
           b. camera_mode is set → POST /camera/mode (manual mode)
           c. plane_lock True → POST /plane-lock/enable (legacy shim)
           d. otherwise → clear all camera state
      6. for each frame: set_frame → screenshot → blob → record

    The instrument-mode toggle is applied after bulkSet because the
    overlay applies CPK-distinct colours to whichever spheres are visible
    at that moment; doing it in this order matches the manual workflow a
    user would follow (pick atoms, then enable instrument mode). The
    transform is applied before instrument mode so the first marker render
    is in the transformed coordinate frame.

    focus_atom=N + focus_atom_kind="plane" exercises the
    /camera/focus_atom endpoint: the helper looks up atom N's residue's
    N/CA/C backbone atoms and applies them as a Plane camera mode.
    Median drift should match the plane-lock floor.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Setup
    clear_selection(base)
    disable_plane_lock(base)  # known-clean start
    clear_camera_mode(base)
    set_atoms(base, atoms)
    # Transform applies upstream of instrument-mode marker placement; the
    # wrapped Conformation re-emits position queries through the rigid-body
    # transform so MeasurementOverlay's per-frame SetCenter call lands at
    # the stabilised atom position. Default kind="identity" → no-op.
    set_transform(base, transform_kind, **(transform_kwargs or {}))
    set_instrument(base, True, focus_only=focus_only)
    if focus_atom is not None:
        set_camera_focus_atom(base, focus_atom, focus_atom_kind)
    elif camera_mode is not None:
        set_camera_mode(base, camera_mode, **(camera_mode_kwargs or {}))
    elif plane_lock:
        enable_plane_lock(base, atoms)
    else:
        disable_plane_lock(base)

    records: list[FrameRecord] = []
    for f in frames:
        set_frame(base, f)
        png = screenshot(base, force_render=True)
        result = find_marker_centroid(png, hue)
        sha = _png_sha1(png)
        if result is None:
            records.append(FrameRecord(
                frame=f, centroid_x=None, centroid_y=None,
                area_px=None, png_sha1=sha,
            ))
        else:
            cx, cy, area = result
            records.append(FrameRecord(
                frame=f, centroid_x=cx, centroid_y=cy,
                area_px=area, png_sha1=sha,
            ))
        if keep_pngs:
            (out_dir / f"frame_{f:05d}.png").write_bytes(png)

    return records


# ----------------------------------------------------------------------------
# Summary statistics — drift relative to frame-0 centroid.
# ----------------------------------------------------------------------------

def summarize(records: list[FrameRecord]) -> dict:
    n_total = len(records)
    found = [r for r in records if r.centroid_x is not None]
    n_found = len(found)
    out: dict = {
        "n_total": n_total,
        "n_marker_found": n_found,
        "n_marker_missing": n_total - n_found,
    }
    if n_found < 2:
        # Need at least two valid points to compute a drift.
        out["marker_drift"] = None
        out["marker_area"] = None
        return out

    # Reference centroid = the first valid record's centroid (the "zero"
    # against which we measure drift). All other valid records get a
    # pixel-distance reading relative to that.
    cx0 = found[0].centroid_x
    cy0 = found[0].centroid_y
    drifts = np.array([
        float(np.hypot(r.centroid_x - cx0, r.centroid_y - cy0))
        for r in found
    ])
    areas = np.array([r.area_px for r in found], dtype=np.float64)

    out["reference_frame"] = found[0].frame
    out["reference_centroid_px"] = {"x": cx0, "y": cy0}
    out["marker_drift"] = {
        "mean_px":     float(drifts.mean()),
        "std_px":      float(drifts.std(ddof=0)),
        "max_px":      float(drifts.max()),
        "median_px":   float(np.median(drifts)),
        "p95_px":      float(np.percentile(drifts, 95)),
        "fraction_within_2px":  float((drifts <= 2.0).mean()),
        "fraction_within_5px":  float((drifts <= 5.0).mean()),
        "fraction_within_10px": float((drifts <= 10.0).mean()),
    }
    out["marker_area"] = {
        "mean_px": float(areas.mean()),
        "std_px":  float(areas.std(ddof=0)),
        "min_px":  int(areas.min()),
        "max_px":  int(areas.max()),
    }
    return out


def write_csv(records: list[FrameRecord], path: Path) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "centroid_x", "centroid_y", "area_px", "png_sha1"])
        for r in records:
            w.writerow([
                r.frame,
                "" if r.centroid_x is None else f"{r.centroid_x:.4f}",
                "" if r.centroid_y is None else f"{r.centroid_y:.4f}",
                "" if r.area_px    is None else r.area_px,
                r.png_sha1,
            ])


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", default="http://127.0.0.1:9988",
                        help="REST base URL (default: %(default)s)")
    parser.add_argument("--atoms", type=int, nargs=3, required=True,
                        help="three atom indices defining the plane")
    parser.add_argument("--frames", type=int, default=200,
                        help="number of frames to sweep (default: %(default)s)")
    parser.add_argument("--frame-stride", type=int, default=1,
                        help="frame stride; the script samples 0, S, 2S, ... up to "
                             "--frames frames (default: %(default)s)")
    parser.add_argument("--slot", type=int, default=0, choices=sorted(MARKER_HUES),
                        help="AtomSelection slot whose marker to track (default: 0)")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="directory for PNG dumps and CSVs")
    parser.add_argument("--no-pngs", action="store_true",
                        help="don't save per-frame PNGs (CSVs + summary only)")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    hue = MARKER_HUES[args.slot]

    health_check(args.base)

    n_atoms = atom_count(args.base)
    for a in args.atoms:
        if a < 0 or a >= n_atoms:
            print(f"atom {a} out of range; protein has {n_atoms} atoms",
                  file=sys.stderr)
            return 2

    n_frames_avail = frame_count(args.base)
    sampled = list(range(0, n_frames_avail, args.frame_stride))[:args.frames]
    if len(sampled) < 2:
        print(f"not enough frames to compute drift (sampled {len(sampled)})",
              file=sys.stderr)
        return 2

    # Hide the docks to expand the central viewport: the marker is more
    # findable in a larger pixel area, AND the drift floor measurement is
    # less affected by dock-induced viewport compression that varies
    # slightly between reader sessions. Restored in a try/finally below.
    set_docks_visible(args.base, False)

    try:
        t0 = time.monotonic()
        floor = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=True, transform_kind="identity",
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "with_lock",
            keep_pngs=not args.no_pngs,
        )
        write_csv(floor, args.out_dir / "with_lock.csv")
        t1 = time.monotonic()

        no_lock = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=False, transform_kind="identity",
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "no_lock",
            keep_pngs=not args.no_pngs,
        )
        write_csv(no_lock, args.out_dir / "no_lock.csv")
        t2 = time.monotonic()

        # Experiment 3: transform-only (fit_subset on backbone) — no
        # plane lock. Tests the upstream layer in isolation. Hypothesis:
        # if rigid-body motion was driving the no-lock drift, this should
        # bring the marker well within the with-lock floor.
        transform_only = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=False,
            transform_kind="fit_subset",
            transform_kwargs={"backbone_only": True, "reference_frame": 0},
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "transform_only",
            keep_pngs=not args.no_pngs,
        )
        write_csv(transform_only, args.out_dir / "transform_only.csv")
        t3 = time.monotonic()

        # Experiment 4: transform + plane lock. Belt+suspenders — if the
        # transform fully stabilises rigid-body motion, this should
        # roughly match transform_only (the lock then has nothing extra
        # to do; in fact the lock could add a small bias if the lock
        # atoms vibrate inside the now-stabilised molecule).
        transform_plus_lock = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=True,
            transform_kind="fit_subset",
            transform_kwargs={"backbone_only": True, "reference_frame": 0},
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "transform_plus_lock",
            keep_pngs=not args.no_pngs,
        )
        write_csv(transform_plus_lock, args.out_dir / "transform_plus_lock.csv")
        t4 = time.monotonic()

        # Experiment 5: CameraMode::Atom on the focus atom (atoms[-1]
        # which the harness has been treating as the focus). Identity
        # transform; the camera absolute-writes the focal to atoms[-1]
        # every frame, so drift should be at the atom-vibration floor.
        focus_atom = args.atoms[-1]
        mode_atom = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=False, transform_kind="identity",
            camera_mode="atom", camera_mode_kwargs={"atom": focus_atom},
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "mode_atom",
            keep_pngs=not args.no_pngs,
        )
        write_csv(mode_atom, args.out_dir / "mode_atom.csv")
        t5 = time.monotonic()

        # Experiment 6: CameraMode::Dihedral on a hard-coded mid-protein
        # backbone psi dihedral. For the 1P9J fixture, residue 27 (LEU,
        # index 26 zero-based) has backbone atom indices N=395, CA=397,
        # C=412, and the next residue's N=414 (residue 28, ASP).
        # Psi(27) = N(27)-CA(27)-C(27)-N(28) = (395, 397, 412, 414) —
        # a well-conditioned mid-protein dihedral that should sit
        # sub-pixel-stable under the new explicit sign-continuity guard
        # (Codex finding #1).
        #
        # The harness selection still tracks the user-supplied
        # args.atoms for marker colouring; only the camera mode's
        # dihedral atoms come from this hard-coded tuple. The selection
        # focus atom remains atoms[-1] from the CLI.
        #
        # If a different fixture is used, update DIHEDRAL_TUPLE to a
        # mid-protein psi or phi dihedral on that fixture (verify via
        # atoms_category_info.npy: pick a non-terminal residue's
        # backbone N/CA/C and the adjacent residue's N).
        DIHEDRAL_TUPLE = (395, 397, 412, 414)  # 1P9J residue 27 psi
        mode_dihedral = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=False, transform_kind="identity",
            camera_mode="dihedral",
            camera_mode_kwargs={
                "a": DIHEDRAL_TUPLE[0], "b": DIHEDRAL_TUPLE[1],
                "c": DIHEDRAL_TUPLE[2], "d": DIHEDRAL_TUPLE[3],
            },
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "mode_dihedral",
            keep_pngs=not args.no_pngs,
        )
        write_csv(mode_dihedral, args.out_dir / "mode_dihedral.csv")
        t6 = time.monotonic()

        # Experiment 7: CameraMode::Subset on backbone — equivalent to
        # today's transform_only's Kabsch but applied to the camera
        # instead of positions. One Kabsch, written to camera; no
        # transform-vs-lock interference. With the writeSubset rotation
        # half implemented, this should match transform_only (~67 px) —
        # before, centroid-only-follow produced ~307 px.
        mode_subset = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=False, transform_kind="identity",
            camera_mode="subset", camera_mode_kwargs={"backbone_only": True},
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "mode_subset_backbone",
            keep_pngs=not args.no_pngs,
        )
        write_csv(mode_subset, args.out_dir / "mode_subset_backbone.csv")
        t7 = time.monotonic()

        # Experiment 8: /camera/focus_atom on the focus atom with kind=plane.
        # The helper derives a 3-atom plane lock from the focus atom's
        # residue's N/CA/C backbone. Should produce median ≤5 px (matches
        # with_plane_lock floor) — it IS a plane lock, with the atoms
        # picked from chemistry-typed residue topology instead of typed
        # by the harness.
        focus_atom = args.atoms[-1]
        focus_atom_local = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=False, transform_kind="identity",
            focus_atom=focus_atom, focus_atom_kind="plane",
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "focus_atom_local",
            keep_pngs=not args.no_pngs,
        )
        write_csv(focus_atom_local, args.out_dir / "focus_atom_local.csv")
        t8 = time.monotonic()

        # Experiment 9: the canonical two-step "focus atom + whole protein
        # steady" recipe documented in HARNESS_BASELINE_PIPELINE_2026-05-30.
        # Stage 1 stabilises positions via TransformedConformation::
        # FitSubset(backbone); Stage 2 keeps the focal on the focus atom
        # via CameraMode::Atom. The two layers compose additively (Atom
        # is translation-only, no rotation interference). Should match
        # transform_plus_lock (~30 px).
        focus_atom_global = run_drift_experiment(
            args.base, args.atoms,
            plane_lock=False,
            transform_kind="fit_subset",
            transform_kwargs={"backbone_only": True, "reference_frame": 0},
            camera_mode="atom", camera_mode_kwargs={"atom": focus_atom},
            frames=sampled, hue=hue,
            out_dir=args.out_dir / "focus_atom_global_compose",
            keep_pngs=not args.no_pngs,
        )
        write_csv(focus_atom_global, args.out_dir / "focus_atom_global_compose.csv")
        t9 = time.monotonic()
    finally:
        # Restore the reader to a clean state for any follow-up runs.
        set_instrument(args.base, False)
        disable_plane_lock(args.base)
        clear_camera_mode(args.base)
        set_transform(args.base, "identity")
        clear_selection(args.base)
        set_docks_visible(args.base, True)

    summary = {
        "atoms": list(args.atoms),
        "slot_tracked": args.slot,
        "marker": asdict(hue),
        "frames_sampled": len(sampled),
        "frame_stride": args.frame_stride,
        "frames_available": n_frames_avail,
        "atoms_total": n_atoms,
        "with_plane_lock": summarize(floor),
        "without_plane_lock": summarize(no_lock),
        "transform_only": summarize(transform_only),
        "transform_plus_lock": summarize(transform_plus_lock),
        "mode_atom": summarize(mode_atom),
        "mode_subset_backbone": summarize(mode_subset),
        "focus_atom_local": summarize(focus_atom_local),
        "focus_atom_global_compose": summarize(focus_atom_global),
        "timing_s": {
            "with_lock":                  round(t1 - t0, 3),
            "no_lock":                    round(t2 - t1, 3),
            "transform_only":             round(t3 - t2, 3),
            "transform_plus_lock":        round(t4 - t3, 3),
            "mode_atom":                  round(t5 - t4, 3),
            "mode_dihedral":              round(t6 - t5, 3),
            "mode_subset_backbone":       round(t7 - t6, 3),
            "focus_atom_local":           round(t8 - t7, 3),
            "focus_atom_global_compose":  round(t9 - t8, 3),
        },
        "mode_dihedral": summarize(mode_dihedral),
    }
    print(json.dumps(summary, indent=2))

    summary_path = args.out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
