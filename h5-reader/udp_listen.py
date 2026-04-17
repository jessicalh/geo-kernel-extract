#!/usr/bin/env python3
"""Tail h5reader's UDP log stream on port 9997.

The reader emits structured JSON datagrams via StructuredLogger (see
src/diagnostics/StructuredLogger.cpp). This script binds the same port
and prints each datagram one per line.

WARNING: Linux unicast UDP delivers each datagram to exactly one socket.
If the reader is open, it likely binds 9997 for its Operations Log dock,
and this listener will receive nothing. Run this script ONLY when the
reader is not running — it is for batch/CLI debugging, not live viewer
sessions.

SO_REUSEADDR + SO_REUSEPORT are set so neither side errors on "address
already in use", but reception is still single-socket.

Usage:
    python3 udp_listen.py                     # foreground, ^C to stop
    python3 udp_listen.py > log.jsonl         # capture to file
    python3 udp_listen.py --pretty            # format JSON human-readably
"""

from __future__ import annotations

import argparse
import json
import socket
import sys
from datetime import datetime


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=9997)
    parser.add_argument(
        "--pretty", action="store_true",
        help="Decode JSON and render one field per line, human-readable.",
    )
    args = parser.parse_args()

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEPORT, 1)
    sock.bind((args.host, args.port))
    print(f"Listening on UDP {args.host}:{args.port} ...", file=sys.stderr)

    try:
        while True:
            data, _src = sock.recvfrom(4096)
            raw = data.decode("utf-8", errors="replace")
            if args.pretty:
                try:
                    obj = json.loads(raw)
                    ts = obj.get("ts", "?")
                    sev = obj.get("severity", "?")
                    cat = obj.get("category", "?")
                    thr = obj.get("thread", "?")
                    msg = obj.get("message", "")
                    print(f"{ts}  {sev:<7}  {thr:<10}  {cat:<28}  {msg}",
                          flush=True)
                except json.JSONDecodeError:
                    print(raw, flush=True)
            else:
                print(raw, flush=True)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
