# Input directory contract — `h5reader_manifest.toml`

**Status**: settled 2026-05-30. Schema version 1. Implementation:
`src/io/ReaderInputManifest.{h,cpp}`, wired into
`src/io/QtProteinLoader.cpp::LoadRunPath`.

## Motivation

h5-reader used to discover the layout of an input directory by trying a
hard-coded list of relative paths (the loader's "did `trajectory.h5`
exist there? at `extract/trajectory.h5`? was a `npys/` dir alongside?
was there a `dft/jobs` sibling?"). That worked for the canonical 1P9J
fixture and broke silently for variants. Per the project's
`feedback_no_file_discovery` discipline — "no try-and-fail, no glob, no
regex on paths; documented conventions only; if an expected file is
missing, log it and stop" — we replace the heuristics with a
declarative manifest at the directory root.

The manifest IS the contract. If it says a file is there, the file is
there. If it doesn't mention a thing, the thing isn't part of this run.

## File location

`<run_root>/h5reader_manifest.toml`

`<run_root>` is the directory the user passes to `h5reader` on the CLI
(or chooses via File ▸ Open Directory). All paths declared in the
manifest are resolved relative to `<run_root>`; absolute paths pass
through unchanged.

## Schema (v1)

```toml
schema_version = 1
run_kind       = "trajectory"      # or "single_pose" or "mutant_pair"
protein_id     = "1P9J"
human_name     = "1P9J calibration trajectory"   # optional; informational
reference_pdb  = "extract/reference.pdb"          # optional, top-level; file

# populated only when run_kind = "trajectory"
[trajectory]
trajectory_h5        = "extract/trajectory.h5"   # required, file
topology_sidecar_dir = "extract"                 # required, dir
per_frame_npys_dir   = "extract/npys"            # optional, dir

# populated only when run_kind = "single_pose"
[single_pose]
pose_dir = "./"                                   # required, dir
# pose_dir holds the 5-NPY topology sidecar AND the flat per-atom
# calculator NPYs (per the producer's --pdb / --orca single-pose layout).

# populated only when run_kind = "mutant_pair"
[mutant_pair]
wt_pose_dir  = "wt/"                              # required, dir
ala_pose_dir = "ala/"                             # required, dir
# h5-reader auto-opens WT in the active window and exposes "Open mutant
# alternate (ALA)…" as a File menu action that launches a fresh process
# on the ALA dir. Single-protein-per-window is preserved (see CLAUDE.md).

# optional, any run_kind
[dft]
jobs_dir = "dft/jobs"                             # required-within-section, dir

# optional informational pass-through; not validated
[metadata]
extracted_at    = 2026-05-28T14:33:21Z            # TOML datetime, or ISO 8601 string
extract_version = "0.5.0"
fixture_kind    = "calibration"                   # calibration | production | test
```

### Required vs optional

| Key                                   | Required           | Type      |
|---------------------------------------|--------------------|-----------|
| `schema_version`                      | yes                | integer   |
| `run_kind`                            | yes                | string    |
| `protein_id`                          | yes                | string    |
| `human_name`                          | no                 | string    |
| `reference_pdb`                       | no                 | file path |
| `[trajectory].trajectory_h5`          | when trajectory    | file path |
| `[trajectory].topology_sidecar_dir`   | when trajectory    | dir path  |
| `[trajectory].per_frame_npys_dir`     | no                 | dir path  |
| `[single_pose].pose_dir`              | when single_pose   | dir path  |
| `[mutant_pair].wt_pose_dir`           | when mutant_pair   | dir path  |
| `[mutant_pair].ala_pose_dir`          | when mutant_pair   | dir path  |
| `[dft].jobs_dir`                      | if `[dft]` present | dir path  |
| `[metadata].*`                        | no                 | any       |

### Validation

The loader (`ReaderInputManifest::Load`) enforces:

1. `schema_version` must equal `kSupportedSchemaVersion` (currently 1).
   Reject with a clear error otherwise — supports forward-compatibility
   when a future schema bump lands.
2. `run_kind` must parse to one of `Trajectory` / `SinglePose` /
   `MutantPair`. Unknown values are rejected.
3. The per-`run_kind` table must be present with its required keys.
   A `run_kind = "trajectory"` manifest with no `[trajectory]` block
   is an error.
4. Every declared path is resolved to an absolute file-system path and
   existence-checked. Required paths missing on disk → hard error.
   Optional paths declared but missing → hard error (manifest lies).
   Optional paths absent → silently empty.
5. The directory/file distinction is enforced: `trajectory_h5` and
   `reference_pdb` must be files; sidecar dirs and pose dirs must be
   directories.

`metadata` is informational and never fails validation.

## Loader behaviour

`QtProteinLoader::LoadRunPath(path)`:

1. Resolve `manifestDir` (= `path` if dir; else `path`'s parent).
2. If `<manifestDir>/h5reader_manifest.toml` exists:
   * Parse + validate. If invalid → fail the load with the parse error.
   * If valid → dispatch by `runKind`:
     - `Trajectory` → `LoadTrajectory(manifest.trajectoryH5, manifest.topologySidecarDir, manifest.perFrameNpysDir)`
     - `SinglePose` → `LoadPose(manifest.poseDir)`
     - `MutantPair` → `LoadPose(manifest.wtPoseDir)` (auto-open WT)
3. If the manifest is absent:
   * Emit a CRITICAL log line naming the missing manifest and the
     `tools/generate_manifest.py` bridge tool.
   * Fall back to the legacy bounded-convention check (search for
     `trajectory.h5` and `extract/trajectory.h5` in `path` if it's a
     directory; else treat the path as the H5 itself; else single pose).
   * `[dft].jobs_dir` is replaced by the legacy `locateDftJobsDir`
     heuristic.

The fallback is **transitional**. Once existing fixtures have manifests
(via `tools/generate_manifest.py` or hand-editing), the fallback becomes
a hard error. See the deprecation timeline below.

## Cross-checks against extractor data

When the load completes, two protein-identifier sources are compared:

1. The input manifest's `protein_id` (this file).
2. The extractor's sidecar `extraction_manifest.json` `protein_id`
   (already cross-checked against the H5 root attribute in
   `QtProteinLoader::LoadTrajectory`).

Mismatches log a warning but do not fail the load. The motivating
scenario: someone copies a fixture and updates the input manifest but
forgets to update the sidecar. The warning surfaces the inconsistency
without blocking the session.

## Mutant pair semantics

The producer's `--mutant` mode emits two parallel single-pose runs (WT
and ALA), each with its own sidecar + per-atom NPYs. h5-reader is
single-protein-per-window per `CLAUDE.md`, so:

* The mutant manifest declares both pose dirs.
* `LoadFromManifest` for `MutantPair` opens the WT pose
  (`manifest.wtPoseDir`).
* `ReaderMainWindow` adds a File menu action "Open mutant alternate
  (ALA)…" that launches a fresh `h5reader` process on the ALA dir,
  using the same `QProcess::startDetached` path as recent files.

The alternate dir typically itself contains an `h5reader_manifest.toml`
(a single-pose manifest pointing at the ALA-only data). When opened
that way, the alternate window doesn't know it's "the other half" of a
mutant pair — it's just a single pose in its own right. The pairing is
a top-level fact at the parent directory; each pose is independently
addressable.

## Generating manifests for existing fixtures

`tools/generate_manifest.py` walks an existing layout and writes the
TOML. Detection rules match the loader's legacy fallback exactly:

```
python3 tools/generate_manifest.py <run_root>
python3 tools/generate_manifest.py --mutant-pair WT_DIR ALA_DIR --root <run_root>
python3 tools/generate_manifest.py --dry-run <run_root>     # print only
python3 tools/generate_manifest.py --force <run_root>       # overwrite existing
```

The script reads `protein_id` from `extraction_manifest.json` when
present; falls back to the root directory name. The user should
hand-edit the generated file if the extractor's `protein_id` is wrong
(some older fixtures have `"extract"` as the value, a producer bug).

## Deprecation timeline

| Phase | Reader behaviour |
|-------|------------------|
| **0.5.x (current)** | Manifest is preferred; legacy fallback emits CRITICAL log but still loads. |
| **0.6.0** | Manifest is required for new fixtures; existing pre-0.6 fixtures get a one-release grace period (loader prints the migration command and continues). |
| **0.7.0** | Legacy fallback removed. `LoadRunPath` requires a valid manifest at the run root. |

Bumps the loader version constant in `main_reader.cpp` (currently
`H5READER_VERSION="0.5.0-rewrite"`).

## Schema evolution

`schema_version` is an integer. Bumps:

* `1` (current) — initial schema, three run kinds.
* `2` (proposed) — would add e.g. `[replicas]` for multi-trajectory
  comparison, or a `[trajectory.frames]` block for stride/start/end
  filtering. The parser must accept v1 manifests forever; new fields
  are added under a v2 schema and the loader picks the right code path.

The parser explicitly rejects unknown `schema_version` values today
rather than silently ignoring unknown keys. If you find yourself
wanting to add fields, bump the version and update both the parser
and this doc in the same commit.

## See also

* `src/io/ReaderInputManifest.h` — typed surface
* `src/io/ReaderInputManifest.cpp` — toml++ parser
* `src/io/QtProteinLoader.cpp::LoadFromManifest` — dispatch
* `tools/generate_manifest.py` — bridge for existing fixtures
* `extern/tomlplusplus/README.md` — vendored single-header parser
* memory: `feedback_no_file_discovery` — the principle this implements
