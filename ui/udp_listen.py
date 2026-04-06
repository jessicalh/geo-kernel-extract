#!/usr/bin/env python3
"""Listen for OperationLog UDP messages on port 9998.

Usage:
    python3 ui/udp_listen.py &          # background, prints to stdout
    python3 ui/udp_listen.py > log.txt & # capture to file
"""
import socket, sys

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.bind(("127.0.0.1", 9998))
print("Listening on UDP 127.0.0.1:9998 ...", file=sys.stderr)

try:
    while True:
        data, _ = sock.recvfrom(4096)
        print(data.decode("utf-8", errors="replace"), end="", flush=True)
except KeyboardInterrupt:
    pass
