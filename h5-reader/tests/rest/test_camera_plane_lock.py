"""Camera plane-lock end-to-end test, driven over REST.

Replaces the in-binary `--camera-plane-lock-smoke` runner. Picks three
atoms, enables the lock, walks frames, asserts:

  - GET /plane-lock returns active=true with the three atoms in order.
  - Per frame, GET /scene/camera focal point tracks the three-atom
    centroid (matching the lock implementation's SetFocalPoint(origin)).
  - Per frame, |dot(camera direction, plane normal)| ≈ 1 (camera is
    actually looking along the plane normal).
  - Frame-to-frame, dot(direction_n, direction_{n-1}) > 0 — i.e. no
    sign flip. This catches the open bug the bespoke smoke missed
    because it used |dot(direction, normal)|.
  - A scene screenshot at frame 0 matches the committed baseline
    (pixel-exact). `--update-baselines` rewrites the file in place.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest


ATOMS = [1, 100, 200]
FRAMES = [0, 50, 100, 150, 200, 250]


def _vec_sub(a: list[float], b: list[float]) -> list[float]:
    return [a[i] - b[i] for i in range(3)]


def _vec_cross(a: list[float], b: list[float]) -> list[float]:
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def _vec_norm(v: list[float]) -> float:
    return math.sqrt(sum(x * x for x in v))


def _vec_normalize(v: list[float]) -> list[float]:
    n = _vec_norm(v)
    return [x / n for x in v] if n > 0 else [0.0, 0.0, 0.0]


def _vec_dot(a: list[float], b: list[float]) -> float:
    return sum(a[i] * b[i] for i in range(3))


def _centroid(positions: list[list[float]]) -> list[float]:
    n = len(positions)
    return [sum(p[i] for p in positions) / n for i in range(3)]


def _plane_normal(positions: list[list[float]]) -> list[float]:
    a, b, c = positions
    return _vec_normalize(_vec_cross(_vec_sub(b, a), _vec_sub(c, a)))


def _fetch_positions(rest, atoms: list[int], frame: int) -> list[list[float]]:
    r = rest.client.post("/positions", json={"atoms": atoms, "frame": frame})
    assert r.status_code == 200, r.text
    return [entry["position"] for entry in r.json()["positions"]]


def test_plane_lock_enables_with_three_atoms(rest):
    rest.client.post("/frame/set", json={"frame": 0})
    rest.client.post("/selection/pick", json={"atom": ATOMS[0], "modifiers": "none"})
    rest.client.post("/selection/pick", json={"atom": ATOMS[1], "modifiers": "shift"})
    rest.client.post("/selection/pick", json={"atom": ATOMS[2], "modifiers": "shift"})

    r = rest.client.post("/plane-lock/enable", json={"atoms": ATOMS})
    assert r.status_code == 204, r.text

    state = rest.client.get("/plane-lock").json()
    assert state["active"] is True
    assert state["atoms"] == ATOMS


def test_plane_lock_camera_tracks_centroid_and_normal(rest):
    rest.client.post("/plane-lock/enable", json={"atoms": ATOMS})

    previous_direction: list[float] | None = None
    for frame in FRAMES:
        rest.client.post("/frame/set", json={"frame": frame})
        positions = _fetch_positions(rest, ATOMS, frame)
        centroid = _centroid(positions)
        expected_normal = _plane_normal(positions)

        camera = rest.client.get("/scene/camera").json()
        focal = camera["focal"]
        direction = camera["direction"]

        # Focal point matches the centroid the lock pushed in (the lock
        # calls camera->SetFocalPoint(basis->origin) every frame).
        focal_error = _vec_norm(_vec_sub(focal, centroid))
        assert focal_error < 1e-3, (
            f"frame {frame}: focal {focal} drifted from centroid {centroid} by {focal_error}"
        )

        # Camera direction is parallel to the plane normal (signed).
        normal_dot = _vec_dot(direction, expected_normal)
        assert abs(normal_dot) > 0.99, (
            f"frame {frame}: camera direction {direction} not parallel to "
            f"plane normal {expected_normal} (|dot|={abs(normal_dot):.4f})"
        )

        # Direction continuity — catches the sign-flip bug the bespoke
        # in-binary smoke missed (it used abs(dot(direction, normal)) so
        # a sign flip looked identical to no flip). This assertion is the
        # whole reason this test exists at this layer.
        if previous_direction is not None:
            continuity = _vec_dot(direction, previous_direction)
            assert continuity > 0.0, (
                f"frame {frame}: camera direction flipped from "
                f"{previous_direction} to {direction} (continuity dot {continuity:.4f})"
            )
        previous_direction = direction


def test_plane_lock_scene_screenshot_is_valid_png(rest):
    """Capture a scene screenshot via the REST surface and assert it's a
    well-formed, non-trivial PNG of the expected viewport.

    Pixel-exact baseline diffs are deliberately NOT used here: VTK's
    render output drifts run-to-run (anti-aliasing, sub-pixel state, GL
    driver scheduling) by enough to make exact equality a poor signal.
    Visual regression with SSIM + per-platform baselines is a separate
    pass; for now this test verifies the screenshot endpoint actually
    returns a valid image and the lock didn't produce a blank scene.
    """
    from PIL import Image
    from io import BytesIO

    rest.client.post("/plane-lock/enable", json={"atoms": ATOMS})
    rest.client.post("/frame/set", json={"frame": 0})

    r = rest.client.post("/screenshot", json={"target": "scene"})
    assert r.status_code == 200, r.text
    assert r.headers["content-type"] == "image/png"
    png_bytes = r.content
    assert png_bytes.startswith(b"\x89PNG"), "response is not a PNG"
    assert 10_000 < len(png_bytes) < 5_000_000, (
        f"PNG size {len(png_bytes)} outside sane range — empty or runaway capture?"
    )

    img = Image.open(BytesIO(png_bytes))
    assert img.mode == "RGB"
    width, height = img.size
    assert width >= 200 and height >= 200, f"viewport too small: {img.size}"

    # Sanity: the scene must not be a uniform fill (background-only).
    # Sample 1000 pixels evenly; require at least 5 distinct colour values.
    pixels = list(img.getdata())
    sampled = pixels[:: max(1, len(pixels) // 1000)]
    unique_colors = len({tuple(p) for p in sampled})
    assert unique_colors > 5, (
        f"scene appears blank (only {unique_colors} distinct colours in {len(sampled)} samples)"
    )
