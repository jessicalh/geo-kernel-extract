#!/usr/bin/env python3
"""Tail h5reader's UDP log stream on port 9997.

The reader emits structured JSON datagrams via StructuredLogger (see
src/diagnostics/StructuredLogger.cpp). This script binds the same port
and prints each datagram one per line.

When the configured host is an IPv4 multicast address (224.0.0.0/4,
typically 239.x.y.z), multiple listeners can join the group at the
same port — the reader's Operations Log dock and this script can
both receive simultaneously, which Linux unicast UDP forbids. For
unicast destinations the old "one socket per port" rule still
applies.

SO_REUSEADDR + SO_REUSEPORT are set so neither side errors on
"address already in use".

Usage:
    python3 udp_listen.py                     # foreground, ^C to stop
    python3 udp_listen.py > log.jsonl         # capture to file
    python3 udp_listen.py --pretty            # format JSON human-readably
    python3 udp_listen.py --host 239.255.0.1  # multicast group
"""

from __future__ import annotations

import argparse
import json
import socket
import struct
import sys


def is_ipv4_multicast(host: str) -> bool:
    """True iff host is an IPv4 multicast address (224.0.0.0/4)."""
    try:
        first = int(host.split(".")[0])
    except (ValueError, IndexError):
        return False
    return 224 <= first <= 239


def open_socket(host: str, port: int) -> socket.socket:
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEPORT, 1)
    if is_ipv4_multicast(host):
        # Multicast receive: bind ANY (so the OS hands us datagrams
        # destined for the group), then join the group on all
        # interfaces (INADDR_ANY).
        sock.bind(("", port))
        mreq = struct.pack(
            "4sl", socket.inet_aton(host), socket.INADDR_ANY
        )
        sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)
        print(
            f"Listening on UDP multicast {host}:{port} (joined group on INADDR_ANY) ...",
            file=sys.stderr,
        )
    else:
        sock.bind((host, port))
        print(
            f"Listening on UDP unicast {host}:{port} ...", file=sys.stderr
        )
    return sock


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=9997)
    parser.add_argument(
        "--pretty", action="store_true",
        help="Decode JSON and render one field per line, human-readable.",
    )
    args = parser.parse_args()

    sock = open_socket(args.host, args.port)

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
