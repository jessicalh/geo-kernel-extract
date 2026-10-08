"""Exercise real frame reads and replacement while probing the GUI event loop."""

from concurrent.futures import ThreadPoolExecutor
import os
import time

import httpx
import pytest


pytestmark = pytest.mark.skipif(
    os.environ.get("H5READER_TEST_RESPONSIVENESS") != "1",
    reason="opt-in responsiveness and shutdown checks on local trajectory fixtures",
)


def probe_during(rest, operation):
    delays = []
    with ThreadPoolExecutor(max_workers=1) as executor:
        pending = executor.submit(operation)
        # Use a separate connection: HTTP/1 serializes requests on one socket.
        with httpx.Client(base_url=rest.base_url, timeout=30) as probe:
            until = time.monotonic() + 2
            while not pending.done() or time.monotonic() < until:
                start = time.monotonic()
                probe.get("/health").raise_for_status()
                delays.append(time.monotonic() - start)
                time.sleep(0.03)
        result = pending.result()
    print(f"GUI response: worst {max(delays):.3f}s over {len(delays)} probes")
    assert max(delays) < 1.5, delays
    return result


def test_frame_reads_and_run_replacement_stay_responsive(rest):
    rest.client.post("/selection/pick", json={"atom": 16}).raise_for_status()
    probe_during(rest, lambda: rest.client.post("/frame/set", json={"frame": 12})).raise_for_status()
    # A failed replacement must leave the current trajectory usable.
    failed = probe_during(rest, lambda: rest.client.post(
        "/api/run/load", json={"path": "missing-reader-responsiveness-fixture.LGS"}))
    assert failed.status_code == 409
    assert rest.client.get("/ui/state").json()["loaded"]

    replacement = os.environ["H5READER_REST_RELOAD_FIXTURE"]
    rest.client.post("/frame/set", json={"frame": 26}).raise_for_status()
    result = probe_during(rest, lambda: rest.client.post(
        "/api/run/load", json={"path": replacement}))
    result.raise_for_status()
    assert result.json()["ok"]
    assert rest.client.get("/frame/current").json()["frame"] == 0
    rest.client.post("/selection/pick", json={"atom": 16}).raise_for_status()
    probe_during(rest, lambda: rest.client.post("/frame/set", json={"frame": 30})).raise_for_status()


def test_z_close_during_load_finishes_cleanly(rest):
    with ThreadPoolExecutor(max_workers=1) as executor:
        loading = executor.submit(rest.client.post, "/api/run/load",
                                  json={"path": os.environ["H5READER_REST_FIXTURE"]})
        with httpx.Client(base_url=rest.base_url, timeout=30) as control:
            # Observe admission without sending a second competing request.
            for _ in range(100):
                if loading.done():
                    pytest.fail("fixture loaded before the shutdown overlap could be tested")
                if control.get("/ui/state").json()["runLoading"]:
                    break
                time.sleep(0.01)
            else:
                pytest.fail("did not observe an active load")
            for route, body in (
                ("/frame/set", {"frame": 10}),
                ("/selection/pick", {"atom": 16}),
                ("/api/run/load", {"path": os.environ["H5READER_REST_FIXTURE"]}),
            ):
                rejected = control.post(route, json=body)
                assert rejected.status_code == 409, (route, rejected.text)
            control.post("/shutdown").raise_for_status()
        response = loading.result(timeout=60)
        assert response.status_code == 409, response.text
    assert rest.process.wait(timeout=60) == 0
