"""pytest fixtures for the h5reader REST end-to-end suite.

The session-scoped `h5reader_session` fixture launches `h5reader --rest 0
<fixture-dataset>` as a subprocess, scrapes the chosen port from the
`H5READER_REST_PORT=NNNNN` stderr handshake, and yields an httpx.Client
bound to that port. Tests share the one running binary across the session
to amortize the ~3-5s startup cost (Qt + VTK + topology load).

Env contract:
- `H5READER_BINARY`   — absolute path to the built h5reader executable.
                        CTest sets this via set_tests_properties.
- `H5READER_REST_FIXTURE` — absolute .LGS or single-.LGS directory loadable by
                            h5reader. CTest defaults it to the 1P9J dataset.
- `H5READER_REST_RELOAD_FIXTURE` — optional second load path used to exercise a
                                   real run change through `POST /api/run/load`.
- `H5READER_REST_ADDRESS` — optional literal bind address passed through to
                            `--rest-address`; omitted means Reader's loopback default.

Headless: VTK needs a real GL FBO, which `QT_QPA_PLATFORM=offscreen` does
not provide. On non-Windows hosts with no DISPLAY, the fixture wraps the
binary in `xvfb-run -a` (which provisions a Xvfb-backed display on demand).
On Windows and on dev machines with an X display already set, the binary
runs directly.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Generator

import httpx
import pytest


PORT_RE = re.compile(rb"H5READER_REST_PORT=(\d+)")


def _env_path(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        pytest.skip(f"env {name} not set — REST suite needs an h5reader binary and a fixture trajectory")
    path = Path(value)
    if not path.exists():
        pytest.fail(f"env {name}={value} does not exist on disk")
    return path


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--update-baselines",
        action="store_true",
        default=False,
        help="Rewrite PNG / JSON baseline files under tests/rest/baselines/ instead of asserting.",
    )


@dataclass
class RestSession:
    process: subprocess.Popen
    port: int
    base_url: str
    client: httpx.Client

    def sampled_frames(self, intervals: int) -> list[int]:
        """Return evenly spaced valid frames, including both trajectory ends."""
        if intervals < 1:
            raise ValueError("intervals must be positive")
        response = self.client.post(
            "/frame/set", json={"frame": 1_000_000_000}
        )
        if response.status_code != 204:
            raise RuntimeError(response.text)
        last_frame = self.client.get("/frame/current").json()["frame"]
        self.client.post("/frame/set", json={"frame": 0})
        return sorted(
            {round(last_frame * step / intervals) for step in range(intervals + 1)}
        )

    def reset(self) -> None:
        """Per-test reset: clear plane lock, clear selection, return to frame 0."""
        try:
            self.client.post("/plane-lock/disable")
        except httpx.HTTPError:
            pass
        try:
            self.client.post("/selection/clear")
        except httpx.HTTPError:
            pass
        try:
            self.client.post("/frame/set", json={"frame": 0})
        except httpx.HTTPError:
            pass


@pytest.fixture(scope="session")
def h5reader_session() -> Generator[RestSession, None, None]:
    binary = _env_path("H5READER_BINARY")
    fixture_dir = _env_path("H5READER_REST_FIXTURE")

    env = {**os.environ}
    cmd: list[str] = []
    if os.name != "nt" and not env.get("DISPLAY"):
        xvfb = shutil.which("xvfb-run")
        if not xvfb:
            pytest.skip("no DISPLAY and xvfb-run not on PATH; install xvfb or run under an X session")
        # Xvfb's default depth lacks GLX; VTK's QOpenGLWidget needs both
        # the GLX extension and a 24-bit visual to obtain an FBO.
        cmd.extend([xvfb, "-a", "-s", "-screen 0 1280x720x24 +extension GLX +extension RANDR +render -noreset"])
        # Force Mesa software OpenGL so CI hosts without a GPU still render.
        env.setdefault("LIBGL_ALWAYS_SOFTWARE", "1")
        env.setdefault("GALLIUM_DRIVER", "llvmpipe")
    cmd.extend([str(binary), "--rest", "0"])
    rest_address = env.get("H5READER_REST_ADDRESS")
    if rest_address:
        cmd.extend(["--rest-address", rest_address])
    cmd.append(str(fixture_dir))

    # Redirect both streams to an on-disk log so the kernel pipe buffer
    # never fills (h5reader logs heavily via the Qt structured logger; a
    # PIPE backed by an unread Python end deadlocks the binary within a
    # few hundred ms of dashboard activity).
    log_file = tempfile.NamedTemporaryFile(
        prefix="h5reader_rest_", suffix=".log", delete=False
    )
    log_path = Path(log_file.name)
    log_file.close()
    log_handle = open(log_path, "wb")
    proc = subprocess.Popen(
        cmd,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        env=env,
    )

    # Scrape H5READER_REST_PORT from the log file. Large authoritative H5
    # fixtures loaded over SMB can opt into a longer safety deadline; readiness
    # still comes only from the Reader's real port-handshake event.
    port: int | None = None
    startup_timeout = float(os.environ.get("H5READER_REST_STARTUP_TIMEOUT", "30"))
    if not 1.0 <= startup_timeout <= 1800.0:
        raise RuntimeError("H5READER_REST_STARTUP_TIMEOUT must be between 1 and 1800 seconds")
    deadline = time.monotonic() + startup_timeout
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            log_handle.close()
            tail = log_path.read_bytes()[-2000:].decode(errors="replace")
            raise RuntimeError(
                f"h5reader exited before REST handshake (rc={proc.returncode}); "
                f"log tail (also at {log_path}):\n{tail}"
            )
        try:
            buf = log_path.read_bytes()
        except OSError:
            buf = b""
        m = PORT_RE.search(buf)
        if m:
            port = int(m.group(1))
            break
        time.sleep(0.05)
    if port is None:
        proc.terminate()
        proc.wait(timeout=5)
        log_handle.close()
        tail = log_path.read_bytes()[-2000:].decode(errors="replace")
        raise RuntimeError(
            f"timed out after {startup_timeout:g}s waiting for "
            f"H5READER_REST_PORT handshake; "
            f"log tail (also at {log_path}):\n{tail}"
        )

    base_url = f"http://127.0.0.1:{port}"
    # The all-metric display sweep intentionally backfills a 24-frame window.
    # On Windows that may cold-load thousands of NPY arrays for a single request.
    client = httpx.Client(base_url=base_url, timeout=60.0)

    # One-shot health check — fail fast if the server isn't actually serving.
    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline:
        try:
            r = client.get("/health")
            if r.status_code == 200 and r.json().get("ok"):
                break
        except httpx.HTTPError:
            pass
        time.sleep(0.1)
    else:
        proc.terminate()
        proc.wait(timeout=5)
        raise RuntimeError(f"REST server failed health check on {base_url}")

    session = RestSession(process=proc, port=port, base_url=base_url, client=client)
    yield session

    # Teardown: exercise /shutdown first; terminate/kill are fallbacks.
    if proc.poll() is None:
        try:
            client.post("/shutdown", timeout=5.0)
        except httpx.HTTPError:
            pass
    client.close()
    if proc.poll() is None:
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=5)
    try:
        log_handle.close()
    except Exception:
        pass
    # Leave log_path on disk on failure so the user can inspect it;
    # remove on clean exit to avoid /tmp clutter.
    if proc.returncode == 0:
        try:
            log_path.unlink()
        except OSError:
            pass


@pytest.fixture
def rest(h5reader_session: RestSession) -> Generator[RestSession, None, None]:
    """Per-test fixture that resets state before yielding."""
    h5reader_session.reset()
    yield h5reader_session


@pytest.fixture
def baselines_dir() -> Path:
    """Per-platform baseline directory. Linux x86_64 only for now."""
    base = Path(__file__).parent / "baselines" / "linux-x64"
    base.mkdir(parents=True, exist_ok=True)
    return base


@pytest.fixture
def update_baselines(request: pytest.FixtureRequest) -> bool:
    return bool(request.config.getoption("--update-baselines"))
