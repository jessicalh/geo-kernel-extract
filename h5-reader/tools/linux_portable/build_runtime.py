#!/usr/bin/env python3
"""Build an offline CPU userspace payload from an installed Linux Reader bundle.

Run on the Ubuntu 24.04 packaging host. No mount, root, Docker or GPU is used.
The output is staging material; this tool never installs onto a source drive.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

BASE_SHA256 = "c1e67ef7b17a6300e136118bd1dc04725009cb376c1aad10abcf8cd453628d58"


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def stage_bundled_apptainer(source: Path, output: Path) -> dict:
    """Stage a previously accepted relocatable engine without modifying the SIF."""
    if str(source).startswith("/mnt/provenance") or str(output).startswith("/mnt/provenance"):
        raise ValueError("Use writable staging; the Provenance mount is excluded.")
    if not (source / "bin/apptainer").is_file():
        raise ValueError("Bundled engine is missing bin/apptainer")
    entries = {}
    for path in sorted(source.rglob("*")):
        if path.is_symlink():
            entries[str(path.relative_to(source))] = {"symlink": os.readlink(path)}
        elif path.is_file():
            entries[str(path.relative_to(source))] = {"sha256": digest(path), "mode": path.stat().st_mode & 0o7777}
    identity = hashlib.sha256(json.dumps(entries, sort_keys=True).encode()).hexdigest()
    target = output / "runtime/apptainer"
    shutil.copytree(source, target, symlinks=True)
    (output / "runtime/apptainer.identity").write_text(identity + "\n")
    record = {"identity": identity, "files": entries,
              "requirements": "glibc 2.28+, permitted user namespaces, local executable workspace without spaces"}
    (output / "runtime/apptainer-manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    return {"identity": identity, "manifest": "runtime/apptainer-manifest.json"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True,
                        help="Verified ubuntu-base-24.04.4-base-amd64.tar.gz")
    parser.add_argument("--proot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bundled-apptainer", type=Path,
                        help="Optional tested relocatable Apptainer tree, copied separately from the SIF")
    parser.add_argument("--build-root", type=Path,
                        help="Separate temporary userspace directory (outside the deliverable)")
    parser.add_argument("--apptainer", default=shutil.which("apptainer") or shutil.which("singularity"),
                        help="Apptainer/Singularity executable for the primary SIF image")
    args = parser.parse_args()
    if not args.apptainer:
        parser.error("Apptainer or Singularity is required to build the primary image")
    # Until installation is explicitly coordinated, this builder must not touch
    # the actual archive mount, including to resolve a path or inspect metadata.
    for path in (args.bundle, args.base, args.proot, args.output, args.build_root, args.bundled_apptainer):
        if str(path).startswith("/mnt/provenance"):
            parser.error("Use writable staging; the Provenance mount is excluded.")
    if digest(args.base) != BASE_SHA256:
        parser.error("Ubuntu base checksum does not match the pinned release")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    root = args.build_root.resolve() if args.build_root else output.parent / (output.name + "-rootfs-build")
    if root.exists():
        parser.error(f"Build directory already exists: {root}; use a fresh output")
    root.mkdir(parents=True)
    subprocess.run(["tar", "--extract", "--gzip", "--file", str(args.base),
                    "--directory", str(root), "--no-same-owner"], check=True)
    shutil.copytree(args.bundle, root / "opt/h5reader", symlinks=True)
    clean_env = dict(os.environ)
    clean_env.pop("LD_LIBRARY_PATH", None)
    clean_env.pop("LD_PRELOAD", None)
    copied: dict[str, str] = {}

    def copy_library(path: Path) -> None:
        if str(path) in copied:
            return
        destination = root / str(path).lstrip("/")
        destination.parent.mkdir(parents=True, exist_ok=True)
        # Dereference at packaging time: the payload cannot depend on host links.
        shutil.copy2(path, destination, follow_symlinks=True)
        copied[str(path)] = digest(destination)

    # Mesa loads the DRI/Gallium driver dynamically, so it is not discovered by
    # the application's ELF dependency walk. Include the software renderer and
    # its complete dependency set; no NVIDIA/CUDA device or driver is bundled.
    seeds = [Path("/usr/lib/x86_64-linux-gnu/libGLX_mesa.so.0"),
             Path("/usr/lib/x86_64-linux-gnu/libEGL_mesa.so.0"),
             Path("/usr/lib/x86_64-linux-gnu/dri/swrast_dri.so")]
    seeds += list(Path("/usr/lib/x86_64-linux-gnu").glob("libgallium-*.so"))
    for seed in seeds:
        if not seed.is_file():
            raise RuntimeError(f"Software renderer dependency is absent: {seed}")
        copy_library(seed)
        report = subprocess.check_output(["ldd", str(seed)], env=clean_env, text=True)
        if "not found" in report:
            raise RuntimeError(report)
        for dependency in re.findall(r"(?:=>\s+)?(/[^\s]+)\s+\(", report):
            copy_library(Path(dependency))
    for directory in ("usr/share/fonts/truetype/dejavu", "usr/share/fonts/truetype/liberation",
                      "etc/fonts", "usr/share/fontconfig", "usr/share/X11/xkb",
                      "etc/ssl/certs"):
        source = Path("/") / directory
        if source.is_dir():
            shutil.copytree(source, root / directory, dirs_exist_ok=True, symlinks=False)
    for directory in ("provenance", "workspace", "dev/shm", "proc", "tmp/.X11-unix"):
        (root / directory).mkdir(parents=True, exist_ok=True)
    templates = Path(__file__).resolve().parent
    shutil.copy2(templates / "guest-start.sh", root / "opt/h5reader/guest-start.sh")
    (root / "opt/h5reader/guest-start.sh").chmod(0o755)
    payload = output / "payload"
    payload.mkdir()
    shutil.copy2(args.proot, payload / "proot-x86_64")
    (payload / "proot-x86_64").chmod(0o755)
    archive = payload / "reader-userspace.tar.gz"
    with archive.open("wb") as destination:
        tar = subprocess.Popen(["tar", "--create", "--directory", str(root), "."], stdout=subprocess.PIPE)
        gzip = subprocess.run(["gzip", "-1", "-n"], stdin=tar.stdout, stdout=destination, check=True)
        tar.stdout.close()
        if tar.wait():
            raise RuntimeError("Could not create userspace archive")
    sif = payload / "reader.sif"
    build_tmp = root.parent / (root.name + "-apptainer-tmp")
    build_tmp.mkdir(exist_ok=True)
    build_env = dict(clean_env, APPTAINER_TMPDIR=str(build_tmp),
                     APPTAINER_CACHEDIR=str(build_tmp / "cache"),
                     CUDA_VISIBLE_DEVICES="-1", HIP_VISIBLE_DEVICES="-1", ROCR_VISIBLE_DEVICES="-1")
    with (root.parent / (root.name + "-sif-build.log")).open("wb") as build_log:
        subprocess.run([args.apptainer, "build", "--mksquashfs-args",
                        "-processors 4 -mem 512M -comp gzip", str(sif), str(root)],
                       env=build_env, stdout=build_log, stderr=subprocess.STDOUT, check=True)
    manifest = {"schema_version": 1, "architecture": "x86_64", "cpu_only": True,
                "reader_sha256": digest(root / "opt/h5reader/h5reader"),
                "local_catalog_sha256": digest(root / "opt/h5reader/linux-local-library.json"),
                "local_catalog_entries": len(json.loads((root / "opt/h5reader/linux-local-library.json").read_text())["datasets"]),
                "base": args.base.name, "base_sha256": BASE_SHA256,
                "runtime_sha256": digest(archive), "proot_sha256": digest(args.proot),
                "sif_sha256": digest(sif), "preferred_backend": "Apptainer/Singularity",
                "proot_version": subprocess.check_output([str(args.proot), "--version"], text=True).strip(),
                "mesa_and_dependencies": copied,
                "requirements": ["x86_64 Linux", "X11 or XWayland", "Apptainer/Singularity",
                                 "PRoot fallback additionally needs ptrace and executable local storage"],
                "security": "Apptainer mounts source data read-only. PRoot fallback is not a security sandbox."}
    if args.bundled_apptainer:
        manifest["bundled_apptainer"] = stage_bundled_apptainer(args.bundled_apptainer, output)
    (output / "runtime-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (payload / "SHA256SUMS").write_text(
        f"{manifest['runtime_sha256']}  reader-userspace.tar.gz\n"
        f"{manifest['proot_sha256']}  proot-x86_64\n")
    (payload / "reader.sif.sha256").write_text(f"{manifest['sif_sha256']}  reader.sif\n")
    for filename in ("start-reader.sh", "Start Reader.desktop", "check-compatibility.sh",
                     "Check compatibility.desktop", "Start Reader - compatibility mode.desktop",
                     "Start Reader - Singularity.desktop", "Start Reader - Apptainer.desktop",
                     "Start Reader - bundled Apptainer.desktop",
                     "prepare-bundled-apptainer.sh",
                     "README.txt", "TRY_NEXT.txt", "ON_THE_DESK.txt", "APPTAINER_ENGINE.txt"):
        shutil.copy2(templates / filename, output / filename)
        (output / filename).chmod(0o644 if filename.endswith(".txt") else 0o755)
    print(json.dumps({"output": str(output), "archive_bytes": archive.stat().st_size,
                      "runtime_sha256": manifest["runtime_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
