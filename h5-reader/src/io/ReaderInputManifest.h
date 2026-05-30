// ReaderInputManifest — declarative descriptor for the layout of an
// input directory the user opens.
//
// Stored at `<run_root>/h5reader_manifest.toml` and parsed via the
// vendored toml++ single header in `extern/tomlplusplus/include/`.
// Replaces the try-and-fail / conventional-path discovery the loader
// used to do, per feedback_no_file_discovery: "no glob, no regex on
// paths. Documented conventions only. If an expected file is missing,
// log it and stop."
//
// Schema v1 covers the three input shapes h5-reader actually loads:
//
//   schema_version = 1
//   run_kind       = "trajectory" | "single_pose" | "mutant_pair"
//   protein_id     = "1P9J"
//   human_name     = "1P9J calibration"             # optional
//   reference_pdb  = "extract/reference.pdb"        # optional, top-level
//
//   # populated only when run_kind = "trajectory"
//   [trajectory]
//   trajectory_h5        = "extract/trajectory.h5"  # required, file
//   topology_sidecar_dir = "extract"                # required, dir
//   per_frame_npys_dir   = "extract/npys"           # optional, dir
//
//   # populated only when run_kind = "single_pose"
//   [single_pose]
//   pose_dir = "./"                                  # required, dir; holds sidecar + per-atom NPYs together
//
//   # populated only when run_kind = "mutant_pair"  (h5-reader auto-opens WT)
//   [mutant_pair]
//   wt_pose_dir  = "wt/"                             # required, dir
//   ala_pose_dir = "ala/"                            # required, dir
//
//   # optional, any run_kind
//   [dft]
//   jobs_dir = "dft/jobs"                            # required-within-section, dir
//
//   # optional informational pass-through; not validated
//   [metadata]
//   extracted_at    = 2026-05-28T14:33:21Z           # TOML datetime or ISO string
//   extract_version = "0.5.0"
//   fixture_kind    = "calibration"                  # calibration | production | test
//
// Validation:
//   * Required key missing                -> hard error (ok=false).
//   * Required path missing on disk       -> hard error.
//   * Optional path declared but missing  -> hard error.
//   * Optional path absent                -> silently empty.
//   * schema_version != kSupportedSchemaVersion -> hard error.
//   * run_kind not in {trajectory, single_pose, mutant_pair} -> hard error.
//
// For mutant_pair: the loader uses `primaryPoseDir()` (WT) for the initial
// load and `alternatePoseDir()` (ALA) for the toolbar/menu switch action.

#pragma once

#include <QString>

namespace h5reader::io {

struct ReaderInputManifest {
    enum class RunKind { Trajectory, SinglePose, MutantPair };

    bool ok = false;
    QString error;          // empty on success

    QString manifestPath;   // absolute path to h5reader_manifest.toml
    QString rootDir;        // absolute path to the directory containing the manifest

    int schemaVersion = 0;
    RunKind runKind = RunKind::Trajectory;  // value only valid when ok=true
    QString proteinId;
    QString humanName;

    // Per-run-kind paths — only the table for the active runKind is
    // populated/validated. The others remain empty.
    QString trajectoryH5;         // trajectory: required, file
    QString topologySidecarDir;   // trajectory: required, dir
    QString perFrameNpysDir;      // trajectory: optional, dir (empty = none)
    QString poseDir;              // single_pose: required, dir
    QString wtPoseDir;            // mutant_pair: required, dir
    QString alaPoseDir;           // mutant_pair: required, dir

    // Optional, any run_kind:
    QString dftJobsDir;           // dir; empty if [dft] omitted
    QString referencePdb;         // file; empty if reference_pdb omitted

    QString extractedAt;
    QString extractVersion;
    QString fixtureKind;

    static constexpr const char* kFileName = "h5reader_manifest.toml";
    static constexpr int kSupportedSchemaVersion = 1;

    // Parse + validate the manifest at `<dir>/h5reader_manifest.toml`.
    // Returns ok=false with `error` set on any failure; the caller
    // decides whether to fall back to convention-based loading or
    // refuse to open the directory.
    static ReaderInputManifest Load(const QString& dir);

    // Cheap existence check (no parse, no validation).
    static bool ExistsAt(const QString& dir);

    // For mutant_pair runs: primaryPoseDir is the WT pose (auto-opened
    // by the loader), alternatePoseDir is the ALA pose (offered as a
    // "Switch to ALA" menu action that spawns a fresh process on the
    // alternate dir). Both empty for non-mutant kinds.
    QString primaryPoseDir() const;
    QString alternatePoseDir() const;

    static const char* NameForRunKind(RunKind k);
};

}  // namespace h5reader::io
