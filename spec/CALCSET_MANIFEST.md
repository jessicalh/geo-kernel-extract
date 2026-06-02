# Calcset manifest — `.LGS` (Lowly Graduate Student)

**Status:** settled 2026-05-31. Schema version 1.

The `.LGS` file is the project-wide top-level wrapper for a *calcset*:
the directory tree containing one trajectory (or one pose, or one
mutant pair) plus its extraction output and optional DFT ground truth.
It is the consumer's entry point — open the directory, read the
`.LGS`, follow its pointers to find every artifact.

## Motivation

Every calcset today carries naming conventions consumers must reverse-
engineer: DFT job dirs named
`batcave_local_..._f{NNNNNN}_t{ps}.0/` with the frame index parsed out
of the substring; extracted frame dirs named
`extract/npys/frame_NNNNNN/` with the same trick; primary ORCA outputs
matched by `*_nmr.out` glob. Each new consumer re-implements the
parsing. The conventions drift; consumers break silently.

The `.LGS` retires the naming tricks by carrying the typed
`frame_index → artifact_path` map explicitly. Consumers stop globbing
and parsing; they read the manifest and chase its pointers.

`.LGS` is *the wrapper*, not the definitions. It points at the
producer-written sidecars (`extraction_manifest.json`,
`{job_id}_meta.json`) for axis, enum, charge, MD5, and per-frame
provenance — it does not duplicate them. The `.LGS` owns: top-level
identity (dataset_id, protein_id, kind), the artifact-pointer table,
and the frame-index → DFT-output map.

## File location and naming

`<calcset_root>/<dataset_id>.LGS`

`<calcset_root>` is the directory a consumer opens (the directory that
already holds `md/`, `extract/`, `dft/`, `topol.top`, etc. for a
trajectory calcset). The `.LGS` sits at the root, named after the
dataset.

`<dataset_id>` is the calcset's canonical name (typically the
calcset_root's basename), e.g. `1p9j-calibration-with-dft`.
The same string appears in the file body under `dataset_id`.

Convention: lowercase, hyphen-separated, no spaces. The wrapper
classes accept the directory path and locate the unique `.LGS` inside.
If a directory contains zero or more than one `.LGS` it is a hard
error — no glob-and-pick.

## Schema (v1)

```jsonc
{
  "schema_version": 1,
  "kind": "trajectory",                  // or "single_pose" or "mutant_pair"
  "dataset_id": "1p9j-calibration-with-dft",
  "protein_id": "1P9J_5801",
  "human_name": "1P9J calibration trajectory (15 ns r²SCAN/def2-SVP DFTs)",

  // present iff kind == "trajectory"
  "trajectory": {
    "md_dir":                "md",
    "topology_top":          "topol.top",
    "extraction_dir":        "extraction",
    "trajectory_h5":         "extraction/trajectory.h5",
    "extraction_manifest":   "extraction/extraction_manifest.json",
    "frame_dt_ps":           10.0,
    "frame_index_basis":     "trr_frame_index",
    "reference_pdb":         "extraction/reference.pdb"   // optional
  },

  // present iff kind == "single_pose"
  "single_pose": {
    "pose_kind":            "orca",      // or "pdb" or "protonated_pdb"
    "pose_dir":             ".",
    "extraction_manifest":  "extraction_manifest.json"
  },

  // present iff kind == "mutant_pair"
  "mutant_pair": {
    "wt_lgs":  "wt/wt.LGS",
    "ala_lgs": "ala/ala.LGS"
  },

  // present iff DFT data exists for this calcset
  "dft": {
    "method": "r2SCAN def2-SVP def2/J NMR  CPCM(Water)",
    "campaign_target_frames": 751,
    "frame_stride": { "first": 0, "last": 1500, "step": 2 },
    "frames": [
      { "frame_index": 0,    "meta_json": "dft/jobs/.../foo_meta.json" },
      { "frame_index": 2,    "meta_json": "dft/jobs/.../bar_meta.json" },
      { "frame_index": 1002, "meta_json": "dft/jobs/.../baz_meta.json" }
    ]
  },

  "metadata": {
    "generated_at_utc":           "2026-05-31T22:00:00Z",
    "lgs_writer":                 "lgs-tools 0.1.1",
    "producer_extractor_version": "0.2.0"
  }
}
```

### Required vs optional by `kind`

| Key                                   | Required             | Type                |
|---------------------------------------|----------------------|---------------------|
| `schema_version`                      | always               | integer (== 1)      |
| `kind`                                | always               | enum string         |
| `dataset_id`                          | always               | string              |
| `protein_id`                          | always               | string              |
| `human_name`                          | always               | string              |
| `trajectory.md_dir`                   | when trajectory      | dir path            |
| `trajectory.topology_top`             | when trajectory      | file path           |
| `trajectory.extraction_dir`           | when trajectory      | dir path            |
| `trajectory.trajectory_h5`            | when trajectory      | file path           |
| `trajectory.extraction_manifest`      | when trajectory      | file path           |
| `trajectory.frame_dt_ps`              | when trajectory      | number (positive)   |
| `trajectory.frame_index_basis`        | when trajectory      | string              |
| `trajectory.reference_pdb`            | optional             | file path           |
| `single_pose.pose_kind`               | when single_pose     | enum string         |
| `single_pose.pose_dir`                | when single_pose     | dir path            |
| `single_pose.extraction_manifest`     | when single_pose     | file path           |
| `mutant_pair.wt_lgs`                  | when mutant_pair     | file path (`.LGS`)  |
| `mutant_pair.ala_lgs`                 | when mutant_pair     | file path (`.LGS`)  |
| `dft`                                 | optional (any kind)  | object              |
| `dft.method`                          | when dft present     | string              |
| `dft.frames`                          | when dft present     | array (authoritative)|
| `dft.campaign_target_frames`          | when dft present     | integer             |
| `dft.frame_stride`                    | when dft present     | object              |
| `dft.frames[].frame_index`            | per entry            | integer (≥ 0)       |
| `dft.frames[].meta_json`              | per entry            | file path           |
| `metadata.generated_at_utc`           | always               | ISO 8601 timestamp  |
| `metadata.lgs_writer`                 | always               | string              |
| `metadata.producer_extractor_version` | when known           | string              |

`kind`-specific sub-objects are mutually exclusive: a `kind:
"trajectory"` manifest has a `trajectory` block and no `single_pose`
or `mutant_pair` block; ditto for the other two.

### `dft.frames` is authoritative; `frame_stride` is decoration

`frame_stride` is the human-readable cadence summary. It tells a
reader at a glance "this campaign sampled every 2nd TRR frame from 0
to 1500". Consumers never trust it for coverage — they iterate
`frames[]`. The two can disagree (partial campaign: `frames.length <
((last - first) / step) + 1`); that's fine and visible at glance.

### `dft.frames[].frame_index` joins through `trajectory.frame_index_basis`

For trajectory calcsets, `dft.frames[].frame_index` is expressed in
the same index basis declared by `trajectory.frame_index_basis`. With
the current 1P9J trajectory convention, that basis is
`trr_frame_index`.

Consumers join DFT rows to a trajectory extraction by matching each
`dft.frames[].frame_index` value to the H5
`/trajectory/frames/original_index` value. This is the declared join
for sparse or irregular DFT coverage over a denser extraction. Do not
derive coverage from `dft.frame_stride`, from the row number in
`dft.frames[]`, or from a stride factor between DFT and extraction
rows.

### `dft.frames[].meta_json` is the per-frame pointer

One field per entry: the relative path to the DFT job's `meta_json`.
The meta.json carries `frame_index` (cross-checkable against the
`.LGS` row), `frame_ps`, `frame_ns`, `files.out_primary` (the ORCA
output filename), `files.out_files[]` (with MD5 + role), `charge`,
`atom_count`, `elapsed_secs`, `orca_exit_code`, etc. Consumers read
the meta.json to get the actual ORCA output path; the wrapper exposes
a `frame.orca_out()` accessor that lazy-loads it.

Rationale: the meta.json IS the per-frame truth. Duplicating its
fields in the `.LGS` would force two-place updates. The `.LGS`'s job
is the *index* (which frames are present, where to find each), not
the per-frame *content*.

## Loader behaviour

`CalcsetManifest::Load(calcset_root_or_lgs_path)`:

1. Resolve `lgs_path`:
   * if path is a directory, look for a single `*.LGS` file inside.
     Zero or > 1 matches → hard error.
   * if path is a file ending in `.LGS`, use directly.
2. Parse JSON. Schema version mismatch → hard error (forward-
   compatible: v1 loaders refuse v2 manifests rather than guessing).
3. Validate `kind`-specific blocks are present; required keys exist;
   declared paths exist on disk; file/dir distinction holds.
4. Optional paths declared but missing on disk → hard error (manifest
   lies). Optional sub-blocks absent → silently absent.

`.LGS` does **not** validate the contents of the artifacts it points
at — that's the consumer's job (e.g. the reader validates trajectory.h5
shape on its own). `.LGS` is the index, not a parser.

## Relationship to producer-written sidecars

Pointed at, never duplicated:

| Producer artifact                  | `.LGS` carries                       |
|------------------------------------|--------------------------------------|
| `extraction_manifest.json`         | a path to it                         |
| `{job_id}_meta.json` (per DFT)     | the typed `frame_index → meta_json` map |
| `extraction/atoms_category_info.npy`, `bonds.npy`, etc. (the 5-NPY topology sidecar) | implicit inside `extraction_dir`     |
| `extraction/npys/frame_NNNNNN/`    | implicit inside `extraction_dir`     |
| `extraction/pdbs/`                 | implicit inside `extraction_dir`     |
| `extraction/trajectory.h5`         | explicit path (may live outside extraction_dir in future fixtures) |
| `extraction/run.log`               | not referenced (producer-internal)   |
| `extraction.pid`, `extraction.log` | not referenced (producer-internal)   |

Where the producer fixes its `protein_id` bug (currently emits "extract"
in `extraction_manifest.json`), `.LGS`'s `protein_id` is the
authoritative one for downstream consumers.

## Scope

`.LGS` is consumer-facing: it points at artifacts produced by an
extraction run (trajectory.h5, NPYs, PDBs, DFT job outputs) so the
reader and downstream tools can locate them without filename parsing.

Out of scope for v1:

* **Raw producer inputs.** `.xyz`, `.prmtop`, `_nmr.out` triples
  consumed by `--orca`/`--mutant` modes are not in `.LGS`. Producer
  workflows continue to derive these from `--root NAME` conventions.
* **Per-MD-file pointers.** `production.tpr` / `.trr` / `.edr` are
  implicit inside `trajectory.md_dir` per producer convention.
* **Per-frame NPY/PDB pointers.** `frame_NNNNNN/` and
  `md_f{NNNNNN}_t{ps}.pdb` are implicit inside
  `trajectory.extraction_dir` per producer naming. Consumers honour
  the convention; they do not enumerate via globbing.

These will gain explicit fields if and when a producer-side `.LGS`
writer (eventual home: the extractor library itself) needs them.
That is a separate, later piece of work — not part of v1.

## Relationship to `h5reader_manifest.toml`

`.LGS` **retires** `h5reader_manifest.toml`. The TOML loader was a
single-week scaffold; the calcset-manifest design supersedes it
before it accreted real fixtures. No deprecation window:

* `h5-reader/src/io/ReaderInputManifest.{h,cpp}` is deleted.
* `h5-reader/tools/generate_manifest.py` is replaced by
  `h5-reader/tools/lgs_write.py` in the same change.
* `h5-reader/extern/tomlplusplus/` is removed (vendored only for the
  TOML loader; JSON parsing uses the C++ standard library / Qt's
  `QJsonDocument`).
* `h5-reader/spec/INPUT_DIRECTORY.md` is marked superseded (header
  note pointing at this spec) but kept as historical record.
* Existing fixtures get `.LGS` generated by `lgs_write.py` in the
  same change so the reader continues to open them.

## C++ wrapper

`h5-reader/src/io/CalcsetManifest.{h,cpp}`. Typed surface:

```cpp
namespace h5reader::io {

struct DftFrame {
    int32_t frame_index = 0;
    QString meta_json_relpath;       // as written in .LGS
    QString meta_json_abspath;       // resolved against calcset_root
    // Accessors that lazy-load the meta.json:
    int32_t framePs() const;
    QString orcaOutAbspath() const;
    int     orcaExitCode() const;
};

struct CalcsetManifest {
    enum class Kind { Trajectory, SinglePose, MutantPair };
    enum class PoseKind { Pdb, ProtonatedPdb, Orca };

    int     schema_version = 0;
    Kind    kind = Kind::Trajectory;
    QString dataset_id;
    QString protein_id;
    QString human_name;

    struct Trajectory {
        QString md_dir_abspath;
        QString topology_top_abspath;
        QString extraction_dir_abspath;
        QString trajectory_h5_abspath;
        QString extraction_manifest_abspath;
        double  frame_dt_ps = 0.0;
        QString frame_index_basis;
    };
    struct SinglePose {
        PoseKind pose_kind = PoseKind::Pdb;
        QString  pose_dir_abspath;
        QString  extraction_manifest_abspath;
    };
    struct MutantPair {
        QString wt_lgs_abspath;
        QString ala_lgs_abspath;
    };
    struct Dft {
        QString method;
        int     campaign_target_frames = 0;
        struct Stride { int first = 0; int last = 0; int step = 1; };
        Stride  frame_stride;
        std::vector<DftFrame> frames;
    };

    std::optional<Trajectory> trajectory;
    std::optional<SinglePose> single_pose;
    std::optional<MutantPair> mutant_pair;
    std::optional<Dft>        dft;

    QString calcset_root_abspath;  // resolved at load time

    // Constants
    static constexpr int kSupportedSchemaVersion = 1;
    static constexpr const char* kExtension = ".LGS";

    // Loader entry point.
    static std::optional<CalcsetManifest> Load(const QString& root_or_lgs_path,
                                                QString* err_out = nullptr);
};

}  // namespace h5reader::io
```

Constructor-style validation; no exceptions across boundary (loader
returns `nullopt` + writes a message to `err_out`).

## Python wrapper

`python/nmr_extract/_calcset.py`. Typed surface (Python dataclasses):

```python
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, List

class Kind(Enum):
    TRAJECTORY = "trajectory"
    SINGLE_POSE = "single_pose"
    MUTANT_PAIR = "mutant_pair"

class PoseKind(Enum):
    PDB = "pdb"
    PROTONATED_PDB = "protonated_pdb"
    ORCA = "orca"

@dataclass(frozen=True)
class DftFrame:
    frame_index: int
    meta_json: Path                 # absolute, resolved
    # accessors that lazy-load the meta.json:
    def frame_ps(self) -> float: ...
    def orca_out(self) -> Path: ...
    def orca_exit_code(self) -> int: ...

@dataclass(frozen=True)
class Trajectory:
    md_dir: Path
    topology_top: Path
    extraction_dir: Path
    trajectory_h5: Path
    extraction_manifest: Path
    frame_dt_ps: float
    frame_index_basis: str

@dataclass(frozen=True)
class SinglePose:
    pose_kind: PoseKind
    pose_dir: Path
    extraction_manifest: Path

@dataclass(frozen=True)
class MutantPair:
    wt_lgs: Path
    ala_lgs: Path

@dataclass(frozen=True)
class Dft:
    method: str
    campaign_target_frames: int
    frame_stride: tuple                  # (first, last, step)
    frames: List[DftFrame]

@dataclass(frozen=True)
class CalcsetManifest:
    schema_version: int
    kind: Kind
    dataset_id: str
    protein_id: str
    human_name: str
    trajectory: Optional[Trajectory]
    single_pose: Optional[SinglePose]
    mutant_pair: Optional[MutantPair]
    dft: Optional[Dft]
    calcset_root: Path

    SCHEMA_VERSION: int = 1
    EXTENSION: str = ".LGS"

    @classmethod
    def load(cls, root_or_lgs_path: Path) -> "CalcsetManifest": ...
```

## Schema evolution

`schema_version` is an integer. Bumps:

* `1` (current) — initial schema, three `kind`s, the `dft` block.
* `≥ 2` (future) — add fields only; old fields keep their meaning.
  Loaders strictly reject unknown versions rather than guess.

If a field's semantics need to change, bump the version and update
both the spec, the parsers, and a migration tool in the same commit.

## Generator tool

`scripts/lgs_write.py` (or `tools/lgs_write.py` — TBD when written):

* Detects existing layouts (mirrors the trajectory / single-pose /
  mutant-pair detection of the loader).
* Writes the `.LGS` with all relative paths resolved and the DFT
  frame list populated from `dft/jobs/*/meta.json` (reading typed
  `frame_index` from each, NOT parsing the dir name).
* `--dry-run` / `--force` / `--root <path>`.

## See also

* `h5-reader/spec/INPUT_DIRECTORY.md` — superseded predecessor (TOML,
  reader-only)
* `python/nmr_extract/_catalog.py` — the per-NPY contract (one
  `ArraySpec` per artifact); `.LGS` is the per-*calcset* analogue
* `memory: project_proper_manifest_design` — the directive that
  motivated this design
* `memory: feedback_data_format_discuss_first` — why JSON was chosen
  through discussion rather than guessed
* `memory: feedback_no_file_discovery` — the discipline this implements
