import socket, sys
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
s.bind(("127.0.0.1", 9998))
s.settimeout(50)
n = 0
try:
    while True:
        d, _ = s.recvfrom(8192)
        print(d.decode("utf-8", "replace").strip(), flush=True)
        n += 1
except socket.timeout:
    pass
print("=== %d msgs ===" % n, flush=True)
