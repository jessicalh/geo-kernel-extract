import copy
from collections import Counter
import io
import json
import lzma
import os
from pathlib import Path
import struct
import subprocess
import sys
import tarfile
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from package_trajectory import package, write_archive


def npy_bytes(value):
    header = b"{'descr': '<f8', 'fortran_order': False, 'shape': (1,), }"
    header += b" " * (-(10 + len(header) + 1) % 64) + b"\n"
    return (
        b"\x93NUMPY\x01\x00"
        + struct.pack("<H", len(header))
        + header
        + struct.pack("<d", value)
    )


class PublicationTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "source"
        self.extraction = self.source / "extraction"
        self.frames = [
            self.extraction / "npys" / f"frame_{index:06d}" for index in (15, 17)
        ]
        self.fields = {
            "pos.npy",
            "aimnet2_charges.npy",
            "mopac_atom_s_population.npy",
            "bond_length.npy",
        }
        for index, frame in enumerate(self.frames):
            frame.mkdir(parents=True)
            for field in sorted(self.fields | {"unused.npy", "mopac_restart.npy"}):
                if field != "bond_length.npy" or index == 0:
                    (frame / field).write_bytes(npy_bytes(index + len(field)))
        (self.extraction / "npys" / "frame_not_a_directory").write_bytes(b"ignored")
        root_files = {
            "trajectory.h5": b"hdf5-data",
            "extraction_manifest.json": b'{"original": true}\r\n',
            "atom_info.npy": npy_bytes(1),
            "topology.npy": npy_bytes(2),
            "topology.json": b'{"bonds": []}\n',
        }
        self.expected_copies = {}
        for name, data in root_files.items():
            (self.extraction / name).write_bytes(data)
            self.expected_copies[f"extraction/{name}"] = data
        (self.extraction / "unused.txt").write_bytes(b"not published")
        reference = self.source / "reference-input.pdb"
        reference.write_bytes(b"HEADER    reference PDB\r\nEND\r\n")
        self.expected_copies["reference.pdb"] = reference.read_bytes()
        dft_frames = []
        for index in (15, 17):
            job = self.source / f"job-{index}"
            job.mkdir()
            meta = job / "meta.json"
            meta.write_bytes(
                b"\xef\xbb\xbf"
                + json.dumps(
                    {"files": {"out_primary": "orca.out"}, "original": index}
                ).encode()
            )
            output = job / "orca.out"
            output.write_bytes(f"Original ORCA output, frame {index}\r\n".encode())
            (job / "restart.gbw").write_bytes(b"not published")
            for path in (meta, output):
                self.expected_copies[f"dft/frame_{index:06d}/{path.name}"] = (
                    path.read_bytes()
                )
            dft_frames.append(
                {
                    "frame_index": index,
                    "meta_json": f"job-{index}/meta.json",
                    "energy": -index,
                }
            )
        self.lgs = self.source / "source.LGS"
        self.manifest = {
            "schema_version": 1,
            "kind": "trajectory",
            "dataset_id": "test-trajectory",
            "protein_id": "example",
            "human_name": "Unchanged identity",
            "annotations": {"keep": [1, 2]},
            "trajectory": {
                "extraction_dir": "extraction",
                "trajectory_h5": str(self.extraction / "trajectory.h5"),
                "extraction_manifest": "extraction/extraction_manifest.json",
                "md_dir": "unneeded-md",
                "topology_top": "unneeded-topology",
                "reference_pdb": "reference-input.pdb",
                "frame_dt_ps": 2.5,
                "frame_indices": [15, 17],
            },
            "dft": {"frames": dft_frames, "method": "keep-original"},
        }
        self.write_manifest()
        fields_file = self.root / "fields.txt"
        fields_file.write_text(
            "\n".join(sorted(self.fields)) + "\n\n pos.npy \n", encoding="utf-8-sig"
        )
        self.args = SimpleNamespace(
            lgs=self.lgs,
            fields=fields_file,
            output=self.root / "published",
            key="example-v2",
            title="Example",
            description="Two sparse frames",
            frames=2,
            install_root=self.root / "installed",
        )
        self.archive = self.args.output / "example-v2.tar.xz"
        self.metadata = self.args.output / "example-v2.json"
        self.original = self.snapshot()

    def write_manifest(self):
        self.lgs.write_text(
            json.dumps(self.manifest, indent=2) + "\n", encoding="utf-8-sig"
        )

    def snapshot(self):
        return {
            path.relative_to(self.source): (path.read_bytes(), path.stat().st_mtime_ns)
            for path in self.source.rglob("*")
            if path.is_file()
        }

    def read_bundle(self):
        compressed = self.archive.read_bytes()
        self.assertTrue(compressed.startswith(b"\xfd7zXZ\x00"))
        decoder = lzma.LZMADecompressor()
        uncompressed = decoder.decompress(compressed)
        self.assertTrue(decoder.eof)
        self.assertEqual(decoder.unused_data, b"")
        with tarfile.open(fileobj=io.BytesIO(uncompressed), mode="r:") as bundled:
            members = bundled.getmembers()
            self.assertTrue(all(item.isfile() for item in members))
            self.assertEqual(len(members), len({item.name for item in members}))
            files = {item.name: bundled.extractfile(item).read() for item in members}
            entry = json.loads(self.metadata.read_text(encoding="utf-8"))
            self.assertEqual(
                entry["expanded_bytes"], sum(item.size for item in members)
            )
            self.assertEqual(entry["archive_bytes"], len(compressed))
            self.assertEqual(entry["frames"], 2)
            self.assertEqual(entry["entry_point"], "run.LGS")
            self.assertEqual(
                entry["url"], "https://semantic.construction/files/example-v2.tar.xz"
            )
        return files

    def expected_arrays(self):
        return {
            f"extraction/npys/{frame.name}/{path.name}": path.read_bytes()
            for frame in self.frames
            for path in frame.glob("*.npy")
            if path.name in self.fields
        }

    def assert_arrays(self, files):
        arrays = {
            name: data
            for name, data in files.items()
            if name.startswith("extraction/npys/")
        }
        self.assertEqual(arrays, self.expected_arrays())
        self.assertEqual(
            list(arrays), sorted(arrays, key=lambda name: (name.split("/")[-1], name))
        )

    def assert_no_publication(self):
        self.assertFalse(self.archive.exists())
        self.assertFalse(self.metadata.exists())
        if self.args.output.exists():
            self.assertEqual(list(self.args.output.iterdir()), [])

    def test_portable_copy_keeps_source_bytes_paths_dfts_and_counts(self):
        with self.assertLogs("reader-package", level="INFO") as logs:
            package(self.args)
        self.assertEqual(self.snapshot(), self.original)
        files = self.read_bundle()
        self.assertEqual(
            set(files),
            set(self.expected_copies)
            | {"run.LGS", "provenance/source.LGS"}
            | set(self.expected_arrays()),
        )
        for name, data in self.expected_copies.items():
            self.assertEqual(files[name], data)
        self.assertEqual(files["provenance/source.LGS"], self.lgs.read_bytes())
        self.assert_arrays(files)
        expected_manifest = copy.deepcopy(self.manifest)
        trajectory = expected_manifest["trajectory"]
        trajectory.pop("md_dir")
        trajectory.pop("topology_top")
        trajectory.update(
            trajectory_h5="extraction/trajectory.h5",
            reference_pdb="reference.pdb",
        )
        for frame in expected_manifest["dft"]["frames"]:
            frame["meta_json"] = f"dft/frame_{frame['frame_index']:06d}/meta.json"
        self.assertEqual(json.loads(files["run.LGS"]), expected_manifest)
        installed = self.args.install_root / self.args.key
        self.assertEqual(
            {
                path.relative_to(installed).as_posix(): path.read_bytes()
                for path in installed.rglob("*")
                if path.is_file()
            },
            files,
        )
        self.assertTrue((installed / "extraction" / "npys").is_dir())
        entry = json.loads(self.metadata.read_text())
        self.assertEqual(
            entry["expanded_bytes"],
            sum(path.stat().st_size for path in installed.rglob("*") if path.is_file()),
        )
        self.assertIn(
            "2 frames, 4 selected fields, 7 NPY arrays, 18 source files",
            "\n".join(logs.output),
        )
        self.assertIn("19 cache files; 7 NPY arrays", "\n".join(logs.output))
        self.assertEqual(set(self.args.output.iterdir()), {self.archive, self.metadata})
        self.assertEqual(list(self.args.install_root.iterdir()), [installed])

    def test_npys_are_readable_with_numpy(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("NumPy is optional for this format interoperability check")
        self.args.install_root = None
        package(self.args)
        files = self.read_bundle()
        for index, frame in enumerate(self.frames):
            array = np.load(
                io.BytesIO(files[f"extraction/npys/{frame.name}/pos.npy"]),
                allow_pickle=False,
            )
            self.assertEqual(array.tolist(), [index + len("pos.npy")])

    def test_source_bytes_are_read_once_and_compression_never_reopens_them(self):
        self.args.install_root = None
        reads = Counter()
        read_bytes = Path.read_bytes
        open_path = Path.open

        def record_read(path):
            reads[path] += 1
            return read_bytes(path)

        def compress_from_memory(path, members):
            def no_source_open(candidate, *args, **kwargs):
                if candidate.is_relative_to(self.source):
                    self.fail(f"Compression reopened source: {candidate}")
                return open_path(candidate, *args, **kwargs)

            with (
                patch.object(Path, "open", no_source_open),
                patch.object(
                    tarfile.TarFile,
                    "gettarinfo",
                    side_effect=AssertionError("Source stat during compression"),
                ),
            ):
                write_archive(path, members)

        with (
            patch.object(Path, "read_bytes", record_read),
            patch("package_trajectory.write_archive", compress_from_memory),
        ):
            package(self.args)
        self.assertTrue(reads)
        self.assertEqual(set(reads.values()), {1})
        self.assert_arrays(self.read_bundle())

    def test_source_read_failure_is_not_silently_omitted(self):
        missing = self.frames[1] / "pos.npy"
        read_bytes = Path.read_bytes

        def fail_read(path):
            if path == missing:
                raise OSError(f"Source read failed: {path}")
            return read_bytes(path)

        with patch.object(Path, "read_bytes", fail_read):
            with self.assertRaisesRegex(OSError, "Source read failed"):
                package(self.args)
        self.assert_no_publication()

    def test_actual_reader_policy_fields(self):
        if "H5READER_PUBLICATION_FIELDS" not in os.environ:
            self.skipTest("Set H5READER_PUBLICATION_FIELDS to test the Reader field policy")
        executable = Path(os.environ["H5READER_PUBLICATION_FIELDS"]).resolve()
        self.assertTrue(
            executable.is_file(),
            f"Reader field-export executable does not exist: {executable}",
        )
        result = subprocess.run(
            [str(executable)], capture_output=True, text=True, check=True
        )
        self.fields = set(result.stdout.splitlines())
        self.assertEqual(len(self.fields), 281)
        self.args.fields.write_text(result.stdout, encoding="utf-8")
        for index, frame in enumerate(self.frames):
            for field in self.fields:
                (frame / field).write_bytes(npy_bytes(index + len(field)))
        original = self.snapshot()
        self.args.install_root = None
        with self.assertLogs("reader-package", level="INFO") as logs:
            package(self.args)
        self.assert_arrays(self.read_bundle())
        self.assertIn(
            "2 frames, 281 selected fields, 562 NPY arrays", "\n".join(logs.output)
        )
        self.assertEqual(self.snapshot(), original)

    def test_no_dft_or_reference_and_no_install(self):
        self.manifest.pop("dft")
        self.manifest["trajectory"].pop("reference_pdb")
        self.write_manifest()
        self.args.install_root = None
        package(self.args)
        files = self.read_bundle()
        self.assertFalse(
            any(name.startswith("dft/") or name == "reference.pdb" for name in files)
        )
        self.assertFalse((self.root / "installed").exists())

    def test_existing_publication_is_not_overwritten(self):
        package(self.args)
        before = {path: path.read_bytes() for path in self.args.output.iterdir()}
        with self.assertRaises(FileExistsError):
            package(self.args)
        self.assertEqual(before, {path: path.read_bytes() for path in before})

    def test_existing_metadata_or_install_is_not_overwritten(self):
        for destination in (self.metadata, self.args.install_root / self.args.key):
            with self.subTest(destination=destination):
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(b"keep")
                with self.assertRaises(FileExistsError):
                    package(self.args)
                self.assertEqual(destination.read_bytes(), b"keep")
                self.assertFalse(self.archive.exists())
                destination.unlink()

    def test_invalid_frame_count_fields_or_key(self):
        for name, value, message in (
            ("frames", 3, "Expected 3"),
            ("frames", 0, "positive"),
            ("key", "../escape", "single directory"),
        ):
            with self.subTest(name=name, value=value):
                original = getattr(self.args, name)
                setattr(self.args, name, value)
                with self.assertRaisesRegex(ValueError, message):
                    package(self.args)
                setattr(self.args, name, original)
                self.assert_no_publication()
        for fields, message in (
            ("unused.npy\n", "missing pos.npy"),
            ("pos.npy\n../escape.npy\n", "filenames"),
        ):
            self.args.fields.write_text(fields)
            with self.assertRaisesRegex(ValueError, message):
                package(self.args)
            self.assert_no_publication()

    def test_missing_position_or_dft_output(self):
        for path in (self.frames[0] / "pos.npy", self.source / "job-15" / "orca.out"):
            with self.subTest(path=path):
                data = path.read_bytes()
                path.unlink()
                with self.assertRaises(FileNotFoundError):
                    package(self.args)
                self.assert_no_publication()
                path.write_bytes(data)

    def test_rejects_nonlocal_dft_output_and_duplicate_destinations(self):
        meta = self.source / "job-15" / "meta.json"
        for output in ("../orca.out", "..\\orca.out", "C:orca.out", ""):
            meta.write_text(json.dumps({"files": {"out_primary": output}}))
            with self.assertRaisesRegex(ValueError, "job-local filename"):
                package(self.args)
            self.assert_no_publication()
        meta.write_text(json.dumps({"files": {"out_primary": "orca.out"}}))
        self.manifest["dft"]["frames"].append(self.manifest["dft"]["frames"][0])
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "Duplicate cache path"):
            package(self.args)
        self.assert_no_publication()

    def test_rejects_wrong_kind_and_h5_outside_extraction(self):
        self.manifest["kind"] = "single_pose"
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "trajectory LGS"):
            package(self.args)
        self.manifest["kind"] = "trajectory"
        self.manifest["trajectory"]["trajectory_h5"] = "job-15/orca.out"
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "declared extraction"):
            package(self.args)
        self.assert_no_publication()

    def test_output_and_install_never_write_inside_source(self):
        for option in ("output", "install_root"):
            for parent in (self.source, self.extraction, self.source / "job-15"):
                with self.subTest(option=option, parent=parent):
                    original = getattr(self.args, option)
                    setattr(self.args, option, parent / "new-publication")
                    with self.assertRaisesRegex(ValueError, "outside source"):
                        package(self.args)
                    self.assertFalse((parent / "new-publication").exists())
                    setattr(self.args, option, original)
        self.assertEqual(self.snapshot(), self.original)

    def test_failed_tar_or_metadata_write_leaves_no_final_or_partial_output(self):
        for target in (
            "tarfile.TarFile.addfile",
            "pathlib.Path.write_text",
        ):
            with self.subTest(target=target):
                with patch(target, side_effect=OSError("simulated write failure")):
                    with self.assertRaisesRegex(OSError, "simulated write failure"):
                        package(self.args)
                self.assert_no_publication()
                self.assertEqual(self.snapshot(), self.original)

    def test_failed_metadata_publish_removes_archive(self):
        replace = Path.replace

        def fail_metadata(path, target):
            if target == self.metadata:
                raise OSError("simulated metadata rename failure")
            return replace(path, target)

        with patch.object(Path, "replace", fail_metadata):
            with self.assertRaisesRegex(OSError, "metadata rename"):
                package(self.args)
        self.assert_no_publication()

    def test_failed_install_keeps_complete_bundle_but_no_partial_cache(self):
        with patch.object(
            tarfile.TarFile,
            "extractall",
            side_effect=OSError("simulated install failure"),
        ):
            with self.assertRaisesRegex(OSError, "install failure"):
                package(self.args)
        self.assert_arrays(self.read_bundle())
        self.assertEqual(list(self.args.install_root.iterdir()), [])
        self.assertEqual(self.snapshot(), self.original)

    def test_existing_cli(self):
        script = Path(__file__).with_name("package_trajectory.py")
        command = [sys.executable, "-B", str(script), str(self.args.lgs)]
        for name in (
            "fields",
            "output",
            "key",
            "title",
            "description",
            "frames",
            "install_root",
        ):
            command.extend(
                ["--" + name.replace("_", "-"), str(getattr(self.args, name))]
            )
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        self.assertIn("Archive ready:", result.stderr)
        self.assert_arrays(self.read_bundle())
        self.assertTrue((self.args.install_root / self.args.key / "run.LGS").is_file())


if __name__ == "__main__":
    unittest.main()
