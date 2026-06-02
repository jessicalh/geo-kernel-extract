from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_force_clean_refuses_nmr_extract_signature_under_run_root(tmp_path: Path) -> None:
    root = tmp_path / "runs"
    run = root / "old"
    extract = run / "extract"
    (extract / "npys" / "frame_000000").mkdir(parents=True)
    (extract / "pdbs").mkdir()
    (extract / "trajectory.h5").write_text("", encoding="utf-8")
    (extract / "extraction_manifest.json").write_text("{}", encoding="utf-8")
    (run / ".rediscover-run.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "created_at": "2020-01-01T00:00:00+00:00",
                "completed_at": "2020-01-01T00:00:00+00:00",
                "substrate_files": ["extract/trajectory.h5"],
            }
        ),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "tools" / "rediscover_run.py"),
            "--root",
            str(root),
            "clean",
            "--keep-substrate",
            "0",
            "--active-minutes",
            "0",
            "--force",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert (extract / "trajectory.h5").exists()
    assert "REFUSE delete" in result.stderr
    assert "nmr_extract extraction signature" in result.stderr
