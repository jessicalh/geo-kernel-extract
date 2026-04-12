#!/usr/bin/env python3
"""Listen for OperationLog UDP messages on port 9998.

For BATCH / CLI sessions only (nmr_extract runs, no viewer open).
Do NOT run alongside the viewer — both bind port 9998, and on Linux
unicast UDP delivers each datagram to exactly one socket (typically
the last to bind). If this script runs first, the viewer's log tab
will be starved. If the viewer runs first, this script will be starved.

SO_REUSEADDR + SO_REUSEPORT are set so this listener is always a
"safe" co-binder: it never gives "address already in use" errors and
never blocks other processes from binding 9998. But it cannot receive
datagrams that go to another socket on the same port.

Usage:
    python3 ui/udp_listen.py            # foreground, ^C to stop
    python3 ui/udp_listen.py > log.txt  # capture to file
"""
import socket, sys

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEPORT, 1)
sock.bind(("127.0.0.1", 9998))
print("Listening on UDP 127.0.0.1:9998 ...", file=sys.stderr)

try:
    while True:
        data, _ = sock.recvfrom(4096)
        print(data.decode("utf-8", errors="replace"), end="", flush=True)
except KeyboardInterrupt:
    pass
