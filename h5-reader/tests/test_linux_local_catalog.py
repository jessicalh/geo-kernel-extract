"""Offline inventory checks; no scientific payload or mounted drive is needed."""

import copy
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "linux"))
import build_local_catalog as catalog


def tsv(rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


class LocalCatalogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.datasets = []
        cls.pairs = []
        for index in range(176):
            member = f"protein_{index}"
            frames, step = (100, 1) if index < 173 else [(751, 2), (1501, 1), (5001, 1)][index - 173]
            group = "MD_173" if index < 173 else "Small_proteins"
            directory = f"01_Reader/datasets/{group}/{member}"
            source = f"{directory}/source.lgs"
            cls.datasets.append({
                "member": member, "dataset_id": member, "title": f"Protein {index}",
                "frames": str(frames), "dft_frames": "0" if index < 173 else str(frames),
                "disk_member_directory": f"/mnt/provenance/{directory}",
                "disk_lgs": f"/mnt/provenance/{source}", "disk_lgs_relative_to_drive": source,
                "staged_lgs": f"/unused/ssd/{source}",
            })
            if index >= 173:
                for frame in range(0, frames * step, step):
                    base = f"frame_{frame:06d}"
                    cls.pairs.append({
                        "member": member, "frame_index": str(frame),
                        "disk_metadata": f"/mnt/provenance/04_MD_records/Small_proteins/{member}/{base}.json",
                        "disk_primary_output": f"/mnt/provenance/04_MD_records/Small_proteins/{member}/{base}.out",
                        "portable_metadata": f"dft/{base}/{base}.json",
                        "portable_primary_output": f"dft/{base}/{base}.out",
                    })
        cls.inventory_bytes = tsv(cls.datasets)
        cls.pair_bytes = tsv(cls.pairs)
        cls.contract_bytes = json.dumps({
            "records": catalog.EXPECTED,
            "data_state": {"verified_counts": {"frame_directories_checked": 24553}},
        }).encode()

    def build(self, datasets=None, pairs=None, **kwargs):
        return catalog.build_catalog(
            self.inventory_bytes if datasets is None else tsv(datasets),
            self.pair_bytes if pairs is None else tsv(pairs), self.contract_bytes, **kwargs,
        )

    def test_complete_catalog_without_any_filesystem_access(self):
        with patch.object(Path, "open", side_effect=AssertionError("No source files may be opened")):
            document, roster = self.build()
        self.assertEqual(document["records"], {**catalog.EXPECTED, "total_frames": 24553})
        self.assertEqual(len(roster["pairs"]), 7253)
        self.assertNotIn("/mnt/provenance", json.dumps((document, roster)))
        self.assertNotIn("url", document["datasets"][0])
        self.assertEqual(document, self.build()[0])

    def test_legitimate_reference_to_raw_records_stays_within_drive(self):
        self.assertEqual(
            catalog.relative_path("../../../../04_MD_records/Small_proteins/x/meta.json",
                                  "01_Reader/datasets/Small_proteins/x"),
            "04_MD_records/Small_proteins/x/meta.json",
        )

    def test_bmrb_title_is_readable_and_source_identity_is_retained(self):
        datasets = copy.deepcopy(self.datasets)
        row = datasets[0]
        original_id = row["dataset_id"]
        original_title = "run_20260611T134525Z_991811 calcset"
        row["title"] = original_title
        row["member"] = "bmr10005"
        for field in ("disk_member_directory", "disk_lgs", "disk_lgs_relative_to_drive"):
            row[field] = row[field].replace("protein_0", "bmr10005")
        document, _ = self.build(datasets=datasets)
        generated = next(entry for entry in document["datasets"] if entry["member"] == "bmr10005")
        self.assertEqual(generated["title"], "BMRB 10005")
        self.assertEqual(generated["source_title"], original_title)
        self.assertEqual(generated["dataset_id"], original_id)
        self.assertEqual(generated["key"], "local-bmr10005")
        self.assertEqual(generated["source_lgs"], row["disk_lgs_relative_to_drive"])

    def test_named_small_proteins_keep_source_titles_and_dataset_identity(self):
        datasets = copy.deepcopy(self.datasets)
        pairs = copy.deepcopy(self.pairs)
        for index, (member, title) in enumerate(catalog.SMALL_PROTEIN_TITLES.items(), start=173):
            row = datasets[index]
            old_member = row["member"]
            row["member"] = member
            for field in ("disk_member_directory", "disk_lgs", "disk_lgs_relative_to_drive"):
                row[field] = row[field].replace(old_member, member)
            for pair in pairs:
                if pair["member"] == old_member:
                    pair["member"] = member
                    for field in ("disk_metadata", "disk_primary_output"):
                        pair[field] = pair[field].replace(old_member, member)
        document, _ = self.build(datasets=datasets, pairs=pairs)
        for index, (member, title) in enumerate(catalog.SMALL_PROTEIN_TITLES.items(), start=173):
            generated = next(entry for entry in document["datasets"] if entry["member"] == member)
            self.assertEqual(generated["title"], title)
            self.assertEqual(generated["source_title"], self.datasets[index]["title"])
            self.assertEqual(generated["dataset_id"], self.datasets[index]["dataset_id"])

    def test_paths_cannot_escape_drive_or_be_urls(self):
        for path in ("../outside.lgs", "/absolute.lgs", "file://bad", "C:\\bad.lgs", "x/../../bad", "bad\nname"):
            with self.subTest(path=path), self.assertRaises(catalog.ContractError):
                catalog.relative_path(path)

    def test_missing_dataset_rejected(self):
        with self.assertRaisesRegex(catalog.ContractError, "Incomplete corpus"):
            self.build(datasets=self.datasets[1:])

    def test_duplicate_members_rejected(self):
        with self.assertRaisesRegex(catalog.ContractError, "Duplicate member"):
            self.build(datasets=self.datasets + [self.datasets[0]])

    def test_inconsistent_source_paths_rejected(self):
        datasets = copy.deepcopy(self.datasets)
        datasets[0]["disk_lgs"] = "/mnt/provenance/outside.lgs"
        with self.assertRaisesRegex(catalog.ContractError, "Inconsistent LGS"):
            self.build(datasets=datasets)

    def test_missing_dft_frame_rejected(self):
        with self.assertRaisesRegex(catalog.ContractError, "Incomplete DFT roster"):
            self.build(pairs=self.pairs[:-1])

    def test_duplicate_dft_frame_rejected(self):
        with self.assertRaisesRegex(catalog.ContractError, "Duplicate DFT frame"):
            self.build(pairs=self.pairs + [self.pairs[0]])

    def test_dft_reference_outside_drive_rejected(self):
        pairs = copy.deepcopy(self.pairs)
        pairs[0]["disk_metadata"] = "/mnt/provenance/../unrelated/meta.json"
        with self.assertRaisesRegex(catalog.ContractError, "escapes the drive"):
            self.build(pairs=pairs)

    def test_dft_destination_cannot_move_frames(self):
        pairs = copy.deepcopy(self.pairs)
        pairs[0]["portable_metadata"] = "dft/frame_000001/frame_000000.json"
        with self.assertRaisesRegex(catalog.ContractError, "Unexpected portable DFT"):
            self.build(pairs=pairs)

    def test_staged_lgs_exact_roster_including_noncontiguous_frame_indices(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rosters = {}
            for row in self.pairs:
                rosters.setdefault(row["member"], []).append({
                    "frame_index": int(row["frame_index"]),
                    "meta_json": "../../../../" + row["disk_metadata"].removeprefix("/mnt/provenance/"),
                })
            for row in self.datasets:
                path = root / row["disk_lgs_relative_to_drive"]
                path.parent.mkdir(parents=True)
                path.write_text(json.dumps({
                    "kind": "trajectory", "dataset_id": row["dataset_id"], "protein_id": "identity-" + row["member"],
                    "dft": {"frames": rosters.get(row["member"], [])},
                }))
            document, roster = self.build(staged_lgs_root=root)
            entry = next(entry for entry in document["datasets"] if entry["member"] == "protein_173")
            self.assertEqual(entry["frames"], 751)
            self.assertEqual(entry["protein_id"], "identity-protein_173")
            self.assertEqual(max(pair["frame_index"] for pair in roster["pairs"] if pair["member"] == "protein_173"), 1500)
            path = root / self.datasets[173]["disk_lgs_relative_to_drive"]
            document = json.loads(path.read_text())
            document["dft"]["frames"][1]["frame_index"] = 1
            path.write_text(json.dumps(document))
            with self.assertRaisesRegex(catalog.ContractError, "Staged LGS DFT roster differs"):
                self.build(staged_lgs_root=root)

    def test_staged_symlink_rejected_without_following_target(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "01_Reader").symlink_to("/does/not/exist")
            with self.assertRaisesRegex(catalog.ContractError, "contains a symlink"):
                self.build(staged_lgs_root=root)


if __name__ == "__main__":
    unittest.main()
