#!/usr/bin/env python3
"""Listen for OperationLog UDP messages on port 9998.

When the configured host is an IPv4 multicast address (224.0.0.0/4,
typically 239.x.y.z), multiple listeners can join the group at the
same port — the viewer's Log tab and this script can both receive
simultaneously, which Linux unicast UDP forbids. For unicast
destinations the old "one socket per port" rule still applies.

For BATCH / CLI sessions, or alongside the viewer when the channel is
configured as multicast in ~/.nmr_tools.toml.

SO_REUSEADDR + SO_REUSEPORT are set so neither side errors on
"address already in use".

Usage:
    python3 ui/udp_listen.py                     # foreground, ^C to stop
    python3 ui/udp_listen.py > log.txt           # capture to file
    python3 ui/udp_listen.py --host 239.255.0.1  # multicast group
"""
from __future__ import annotations

import argparse
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
    parser.add_argument("--port", type=int, default=9998)
    args = parser.parse_args()

    sock = open_socket(args.host, args.port)
    try:
        while True:
            data, _ = sock.recvfrom(4096)
            print(data.decode("utf-8", errors="replace"), end="", flush=True)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
