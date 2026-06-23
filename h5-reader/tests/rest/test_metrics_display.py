"""Read-to-display coverage for EVERY metric descriptor.

The dashboard catalog carries ~188 signal descriptors. By hand, nobody can
verify that each one actually reaches the UI in a sane form -- there are too
many, and "it's in the catalog" says nothing about whether data survives the
descriptor -> sampler -> SignalSample -> panel -> pixels pipeline. This test
drives that pipeline over REST for every descriptor and classifies the result,
so a regression (a metric that used to display stops displaying) fails loudly
and the known gaps stay documented.

What it does, per descriptor (from GET /catalog):
  1. Choose the display path from the descriptor's `axis` (anchor kind) and
     `modes` (a strip.* mode -> temporal track; a static.* mode -> panel).
  2. Add it via POST /dashboard/metric with a valid anchor (iterating a few
     candidate indices, since data presence varies per atom/ring/bond).
  3. Jump to a bounded frame WINDOW so the freshly-added strip backfills that
     many frames (the controller samples frames 0..current on add), then read
     GET /dashboard/display -- the unified manifest of strip tracks + static
     panels (the latter closed the static-panel REST blind spot).
  4. Classify the verdict AND measure variation: a track that never changes
     across the window is "flat"; one with no data is "empty".

The window is bounded on purpose: jumping to the LAST frame backfills the whole
~1500-frame trajectory on EVERY add, which (times the candidate-anchor retries)
runs for many minutes. A few dozen consecutive frames is enough to tell a
genuinely constant track from one that moves frame-to-frame -- "flat" here means
"did not change at all across the window," a strong (window-honest) bupkis
signal, not a whole-trajectory variance claim.

Behaviours this test deliberately treats as EXPECTED, not failures (each was
confirmed against the running app, not assumed):
  - Multi-value-as-one-track: vector3/EFG/spherical/tensor-component metrics
    collapse to ONE scalar per track (a magnitude/invariant trend); the full
    tensor lives only in the focus-driven CSA glyph. Flagged, not failed.
  - static.tensor is a no-op in the dashboard for ALL descriptors that offer it
    (the scene-glyph trigger is deferred); only the SOLE-mode case
    (reorient_orientation_tensor) yields no display at all.
  - Absent-in-this-dataset descriptors are refused with 409 (correct).
  - embedding + the 5 topology tables are non-displayable by policy.

KNOWN_DEBT (calibrated to the default 1P9J fixture) is the small set of genuine
gaps the sweep surfaces; the inventory test asserts it EXACTLY. Flatness is
REPORTED (not asserted): a flat track is still displayed, and many metrics are
legitimately constant -- the report lets a human spot the suspicious ones.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest


# --- verdict vocabulary -------------------------------------------------

CLEAN = {"PASS", "PASS_FLAG_MULTIVALUE", "ABSENT_REFUSED", "NON_DISPLAYABLE"}

# Genuine display debt on the default 1P9J fixture. id -> expected verdict.
#   DEFERRED_GLYPH: the only metric with no dashboard path at all -- its sole
#     mode is static.tensor, whose scene-glyph trigger is deferred
#     (DashboardDisplayController ~1769-1787). Tensors are shown via the
#     focus-driven CSA glyph; wiring a dashboard-driven scene glyph is the
#     deferred A1/A3 renderer work.
#
# (The two empty welfords -- water_field / aimnet2_charge_response_gradient --
# that used to be EMPTY_AVAILABLE are gone from here: the whole rollup-moment
# family was de-stripped (DisplayPolicy) because those summaries only ever drew a
# flat line; they now classify NON_DISPLAYABLE. When a static mean+/-std readout
# ships for rollups, expect them to return as PASS and update accordingly.)
KNOWN_DEBT = {
    "h5:reorient_orientation_tensor": "DEFERRED_GLYPH",
}

SOFT_VERDICTS = {"EMPTY_AVAILABLE", "DEFERRED_GLYPH"}


# --- anchor + mode helpers (ground truth from parseDashboardAnchor) -----

def build_anchor(axis: str, idx: int) -> dict:
    return {
        "atom": {"atom": idx},
        "residue": {"residue": idx},
        "bond": {"bond": idx},
        "bond vector": {"kind": "bond_vector", "residue": idx, "kind_id": 0},
        "ring": {"kind": "ring", "ring": idx},
        "aromatic ring": {"kind": "aromatic_ring", "ring": idx},
        "saturated ring": {"kind": "saturated_ring", "ring": idx},
        "ring contribution pair": {"kind": "ring_contribution_pair", "pair": idx},
        "ring membership": {"kind": "ring_membership", "membership": idx},
        "mutation match pair": {"kind": "mutation_match_pair", "pair": idx},
        "system": {"kind": "system"},
        "event": {"kind": "event"},
        "atom tuple": {"follows_focus": True},
    }.get(axis, {"kind": "none"})


# Candidate anchors to try (data presence varies per anchor). Kept short: each
# attempt backfills a frame window, so retries multiply cost.
def candidate_indices(axis: str) -> list[int]:
    return {
        "atom": [16, 0, 100],
        "residue": [0, 5, 10],
        "bond": [0, 16, 50],
        "bond vector": [0, 5, 10],
        "ring": [0, 1, 2],
        "aromatic ring": [0, 1, 2],
        "saturated ring": [0, 1, 2],
        "ring contribution pair": [0, 1],
        "ring membership": [0, 1],
        "mutation match pair": [0],
    }.get(axis, [0])


# Frames to backfill per add when measuring variation. Bounded for speed.
SAMPLE_WINDOW = 24


STRIP_PREF = [
    "strip.scalar", "strip.count", "strip.category", "strip.per-class",
    "strip.rollup", "strip.event", "strip.system", "strip.vector.magnitude",
    "strip.tensor.T0", "strip.tensor.T2", "strip.tensor.component",
    "strip.vector.component", "strip.tensor.T1",
]
STATIC_VISIBLE = [
    "static.chord.coupling", "static.curve.lag.animated",
    "static.spectrum.power", "static.fixed_freq", "static.bar.sequence",
]
MULTI_SHAPES = {"vector3", "EFG/T2 tensor", "spherical tensor", "tensor components"}


def _pick_mode(modes: list[str], pref: list[str]) -> str | None:
    for m in pref:
        if m in modes:
            return m
    return None


def _err(resp) -> str:
    try:
        return str(resp.json().get("error", resp.text))
    except Exception:
        return resp.text


def _track_variation(values: list) -> tuple[int, float, bool]:
    """(non-null count, span = max-min over non-null, varies) for one track."""
    nums = [v for v in values if v is not None]
    if not nums:
        return 0, 0.0, False
    lo, hi = min(nums), max(nums)
    span = hi - lo
    # "varies" = a visible change over time; tolerate fp noise relative to scale.
    varies = span > 1e-9 * max(1.0, abs(hi), abs(lo))
    return len(nums), span, varies


def _probe_add(client, d: dict, mode: str, idx: int, frame: int) -> dict:
    """Add one metric at (axis, idx), jump to `frame` so the strip backfills that
    many frames, read the manifest, remove it. Never raises on a 4xx."""
    body = {"descriptor_id": d["id"], "anchor": build_anchor(d["axis"], idx), "modes": [mode]}
    client.post("/frame/set", json={"frame": frame})
    r = client.post("/dashboard/metric", json=body)
    if r.status_code != 200:
        return {"ok": False, "code": r.status_code, "msg": _err(r)}
    payload = r.json()
    rid, refs = payload.get("id"), payload.get("added_refs", -1)
    disp = client.get("/dashboard/display").json()
    tracks = disp.get("strip_tracks", []) or []
    panels = disp.get("panels", []) or []
    max_non, vlen, metric_span, any_varies = 0, 0, 0.0, False
    for t in tracks:
        vals = t.get("values") or []
        vlen = max(vlen, len(vals))
        nn, span, varies = _track_variation(vals)
        max_non = max(max_non, nn)
        metric_span = max(metric_span, span)
        any_varies = any_varies or varies
    res = {"ok": True, "code": 200, "refs": refs, "tracks": len(tracks),
           "vlen": vlen, "nonNull": max_non, "span": metric_span, "varies": any_varies,
           "panels": len(panels), "panelKind": "", "panelPts": 0,
           "panelFin": 0, "panelNan": 0, "msg": ""}
    if panels:
        p = panels[0]
        res.update(panelKind=p.get("kind", ""), panelPts=p.get("point_count", 0),
                   panelFin=p.get("finite_count", 0), panelNan=p.get("nan_count", 0))
    if rid:
        client.post("/dashboard/metric/remove", json={"id": rid})
    return res


def _probe_until(client, d: dict, mode: str, want: str, frame: int) -> dict:
    """Try candidate anchors until one yields data (want='track'|'panel')."""
    best = None
    for idx in candidate_indices(d["axis"]):
        r = _probe_add(client, d, mode, idx, frame)
        if r["ok"] and want == "track" and r["tracks"] > 0 and r["nonNull"] > 0:
            return r
        # A panel with points but zero FINITE values renders blank -- keep
        # trying anchors for one with real data before settling.
        if r["ok"] and want == "panel" and r["panelPts"] > 0 and r["panelFin"] > 0:
            return r
        if best is None or (r["ok"] and not best["ok"]):
            best = r
    return best


def _classify(client, d: dict, displayable: bool, frame: int) -> dict:
    """Drive the pipeline for one descriptor and return a verdict record."""
    modes = d.get("modes") or []
    strip_mode = _pick_mode(modes, STRIP_PREF)
    static_mode = _pick_mode(modes, STATIC_VISIBLE)
    only_tensor = bool(modes) and not strip_mode and not static_mode and "static.tensor" in modes
    rec = {"id": d["id"], "value_shape": d.get("value_shape", ""), "axis": d.get("axis", ""),
           "avail": d.get("availability", ""), "path": "", "mode": "", "verdict": "",
           "shape": "", "detail": ""}

    if not modes:
        rec.update(path="none", verdict="NON_DISPLAYABLE", detail=f"displayable={displayable}")
    elif d.get("availability") == "Absent":
        r = _probe_add(client, d, modes[0], 0, frame)
        rec["path"] = "absent"
        if r.get("code") == 409:
            rec.update(verdict="ABSENT_REFUSED", detail="409 not available (expected for this dataset)")
        else:
            rec.update(verdict="ABSENT_UNEXPECTED", detail=f"add code={r.get('code')} {r.get('msg','')}")
    elif only_tensor:
        r = _probe_add(client, d, "static.tensor", 0, frame)
        rec["path"] = "deferred_glyph"
        if r["ok"] and r["panels"] == 0 and r["tracks"] == 0:
            rec.update(verdict="DEFERRED_GLYPH",
                       detail=f"static.tensor refs={r['refs']} -> no dashboard element (scene-glyph trigger deferred)")
        else:
            rec.update(verdict="UNEXPECTED_TENSOR", detail=f"tracks={r['tracks']} panels={r['panels']}")
    elif d.get("axis") == "atom tuple":
        r = _probe_add(client, d, strip_mode, 0, frame)
        rec["path"] = "tuple"
        if r["ok"] and (r["tracks"] > 0 or r["panels"] > 0):
            rec.update(verdict="PASS", shape=("varying" if r["varies"] else "flat"),
                       detail=f"via follows_focus tracks={r['tracks']} span={r['span']:.3g}")
        else:
            rec.update(verdict="NEEDS_SELECTION", detail=f"needs multi-atom pick; code={r.get('code')} {r.get('msg','')}")
    elif strip_mode:
        rec.update(path="temporal", mode=strip_mode)
        r = _probe_until(client, d, strip_mode, "track", frame)
        if not r["ok"]:
            rec.update(verdict="BIND_FAIL", detail=f"code={r.get('code')} {r.get('msg','')}")
        elif r["tracks"] == 0:
            rec.update(verdict="NO_DATA", detail=f"bound but 0 tracks (refs={r['refs']})")
        elif r["nonNull"] == 0:
            rec.update(verdict="EMPTY_AVAILABLE", shape="empty",
                       detail=f"Available + binds + {r['tracks']} track(s) but ALL-NULL across every candidate anchor")
        else:
            rec.update(verdict="PASS", shape=("varying" if r["varies"] else "flat"),
                       detail=f"tracks={r['tracks']} vlen={r['vlen']} nonNull={r['nonNull']} span={r['span']:.3g}")
            if d.get("value_shape") in MULTI_SHAPES:
                rec["verdict"] = "PASS_FLAG_MULTIVALUE"
                rec["detail"] += " (multi-value -> per-track scalar; full tensor in CSA glyph)"
    elif static_mode:
        rec.update(path="static", mode=static_mode)
        r = _probe_until(client, d, static_mode, "panel", frame)
        if not r["ok"]:
            rec.update(verdict="BIND_FAIL", detail=f"code={r.get('code')} {r.get('msg','')}")
        elif r["panels"] == 0:
            rec.update(verdict="NO_DATA", detail=f"bound but 0 panels (refs={r['refs']})")
        elif r["panelPts"] == 0:
            rec.update(verdict="EMPTY_AVAILABLE", shape="empty", detail=f"panel {r['panelKind']} empty")
        elif r["panelFin"] == 0:
            rec.update(verdict="EMPTY_AVAILABLE", shape="blank",
                       detail=f"panel {r['panelKind']} pts={r['panelPts']} but ALL non-finite -> renders BLANK")
        else:
            rec.update(verdict="PASS", shape="static",
                       detail=f"{r['panelKind']} pts={r['panelPts']} fin={r['panelFin']} nan={r['panelNan']}")
    else:
        rec.update(path="unknown", verdict="NO_PATH", detail="modes=" + ",".join(modes))
    return rec


# --- the sweep (run once per session) -----------------------------------

@pytest.fixture(scope="session")
def metrics_sweep(h5reader_session):
    """Classify every catalog descriptor once over the WHOLE trajectory; write a
    TSV report; return the records. Session-scoped to amortize the REST calls."""
    client = h5reader_session.client
    # Jump past the end so playback clamps to the true last frame; read it back.
    client.post("/frame/set", json={"frame": 1_000_000_000})
    last_frame = int(client.get("/frame/current").json().get("frame", 0))
    window = min(SAMPLE_WINDOW, last_frame) if last_frame else SAMPLE_WINDOW

    catalog = client.get("/catalog").json()
    descriptors = catalog.get("descriptors", [])
    assert descriptors, "GET /catalog returned no descriptors -- nothing loaded?"
    audit = client.get("/catalog/display-audit").json()
    displayable = {f["id"]: f.get("displayable", False) for f in audit.get("fields", [])}

    records = [_classify(client, d, displayable.get(d["id"], False), window) for d in descriptors]

    # Report artifact: TSV + a printed summary (pytest -s surfaces it).
    out = Path(os.environ.get("H5READER_METRIC_REPORT")
               or (Path(tempfile.gettempdir()) / "h5reader_metric_report.tsv"))
    header = ["id", "value_shape", "axis", "avail", "path", "mode", "verdict", "shape", "detail"]
    lines = ["\t".join(header)]
    lines += ["\t".join([r["id"], r["value_shape"], r["axis"], r["avail"], r["path"],
                         r["mode"], r["verdict"], r["shape"], r["detail"]]) for r in records]
    out.write_text("\n".join(lines), encoding="ascii", errors="replace")

    summary: dict[str, int] = {}
    for r in records:
        summary[r["verdict"]] = summary.get(r["verdict"], 0) + 1
    print(f"\n=== metrics read-to-display sweep: {len(records)} descriptors "
          f"({window}-frame window, last_frame={last_frame}) ===")
    for v in sorted(summary, key=lambda k: -summary[k]):
        print(f"  {v:<22} {summary[v]}")

    # Flat/empty quantification (the "flat and empty as it flashes by" view):
    temporal = [r for r in records if r["path"] in ("temporal", "tuple")]
    flat = [r for r in temporal if r["shape"] == "flat"]
    varying = [r for r in temporal if r["shape"] == "varying"]
    empty = [r for r in records if r["shape"] == "empty"]
    print(f"  --- temporal display quality: {len(varying)} varying, {len(flat)} FLAT, "
          f"{len(empty)} empty (of {len(temporal)} temporal) ---")
    if flat:
        print(f"  FLAT (constant across the {window}-frame window):")
        for r in flat:
            print(f"      {r['id']:<46} {r['value_shape']}")
    print(f"  report -> {out}")
    return records


# --- assertions ---------------------------------------------------------

def test_coverage_is_total(metrics_sweep, h5reader_session):
    """The sweep classified every catalog descriptor (no silent drop)."""
    total = h5reader_session.client.get("/catalog").json().get("count")
    assert total is not None
    assert len(metrics_sweep) == total, f"swept {len(metrics_sweep)} but catalog has {total}"
    assert all(r["verdict"] for r in metrics_sweep), "some descriptor went unclassified"


def test_no_pipeline_breaks(metrics_sweep):
    """No descriptor lands in a hard-fail bucket. These verdicts mean the
    read-to-display pipeline (or the anchor map) is genuinely broken, as
    opposed to the documented soft gaps in KNOWN_DEBT."""
    hard = {"BIND_FAIL", "NO_DATA", "NO_PATH", "ABSENT_UNEXPECTED",
            "UNEXPECTED_TENSOR", "NEEDS_SELECTION"}
    broken = [(r["id"], r["verdict"], r["detail"]) for r in metrics_sweep if r["verdict"] in hard]
    assert not broken, "pipeline breaks:\n" + "\n".join(f"  {i}  [{v}]  {d}" for i, v, d in broken)


def test_clean_or_known_debt(metrics_sweep):
    """Every descriptor is either cleanly displayed/excluded, or a documented
    debt item with the expected verdict."""
    unexpected = []
    for r in metrics_sweep:
        if r["verdict"] in CLEAN:
            continue
        if KNOWN_DEBT.get(r["id"]) == r["verdict"]:
            continue
        unexpected.append((r["id"], r["verdict"], r["detail"]))
    assert not unexpected, ("descriptors neither clean nor known-debt:\n"
                            + "\n".join(f"  {i}  [{v}]  {d}" for i, v, d in unexpected))


def test_known_debt_inventory(metrics_sweep):
    """The set of soft-gap descriptors EXACTLY matches KNOWN_DEBT. Forces the
    allowlist to track reality: fixing a gap (or a new one appearing) fails
    until KNOWN_DEBT is updated. Calibrated to the default 1P9J fixture."""
    observed = {r["id"]: r["verdict"] for r in metrics_sweep if r["verdict"] in SOFT_VERDICTS}
    assert observed == KNOWN_DEBT, (
        f"known-debt drift:\n  observed={observed}\n  expected={KNOWN_DEBT}\n"
        "  (a gap was fixed, or a new one appeared -- update KNOWN_DEBT)")
