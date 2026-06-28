"""REST coverage for high-resolution (poster / print) scene capture.

`POST /api/screenshot {target:"scene", scale:N}` renders the scene at N x N tiles
and stitches them (vtkWindowToImageFilter SetScale + FixBoundary) -> a high-DPI
image whose extra pixels are effective supersampled antialiasing. Used to export
figure-quality stills for posters/thesis. These tests pin the contract:

  - scale=2 yields EXACTLY twice the pixel dimensions of scale=1 (the tiled
    render stitched without dropping/duplicating rows),
  - the result is still a valid, non-blank PNG (depth-peeled translucency +
    overlay layer survive the tiled capture),
  - omitting scale is identical to scale=1.
"""

from __future__ import annotations

from io import BytesIO


def _scene_png(rest, scale=None):
    body = {"target": "scene"}
    if scale is not None:
        body["scale"] = scale
    r = rest.client.post("/api/screenshot", json=body)
    assert r.status_code == 200, r.text
    assert r.headers["content-type"] == "image/png"
    assert r.content.startswith(b"\x89PNG"), "response is not a PNG"
    return r.content


def test_scale_two_doubles_dimensions(rest):
    """scale=2 is exactly 2x scale=1 in each axis, and still a real scene."""
    from PIL import Image

    rest.client.post("/frame/set", json={"frame": 0})

    base = Image.open(BytesIO(_scene_png(rest, scale=1)))
    assert base.mode == "RGB"
    w1, h1 = base.size
    assert w1 >= 200 and h1 >= 200, f"viewport too small: {base.size}"

    big = Image.open(BytesIO(_scene_png(rest, scale=2)))
    assert big.size == (w1 * 2, h1 * 2), (
        f"scale=2 should be 2x {(w1, h1)} -> {(w1 * 2, h1 * 2)}, got {big.size}"
    )
    # The supersampled capture must not be a blank tile-stitch.
    pixels = list(big.getdata())
    sampled = pixels[:: max(1, len(pixels) // 1000)]
    assert len({tuple(p) for p in sampled}) > 5, "supersampled scene appears blank"


def test_scale_omitted_matches_scale_one(rest):
    """No scale field == scale 1 (same dimensions)."""
    from PIL import Image

    rest.client.post("/frame/set", json={"frame": 0})
    default = Image.open(BytesIO(_scene_png(rest)))
    one = Image.open(BytesIO(_scene_png(rest, scale=1)))
    assert default.size == one.size, (
        f"omitting scale ({default.size}) must equal scale=1 ({one.size})"
    )
