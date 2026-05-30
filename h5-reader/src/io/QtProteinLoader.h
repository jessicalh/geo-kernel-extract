// QtProteinLoader — top-level orchestrator: a run directory into a typed
// QtProtein + a Conformation (the base; concrete subclass depends on run shape).
//
// Two entry points, one per run shape (decided 2026-05-26):
//   Load(h5_path)     — a --trajectory run: sidecar + trajectory.h5 ->
//                       QtProtein + TrajectoryConformation (N frames, dense H5,
//                       lazy per-frame snapshots from per_frame_npys/ or npys/).
//   LoadPose(run_dir) — a single-pose run (--orca / --mutant / --pdb): sidecar +
//                       the flat run-root per-atom NPYs -> QtProtein +
//                       SingleConformation (one frame, no H5).
//
// Failure modes log via ErrorBus and return ok=false with error populated.
// Per-TR-group reads inside QtTrajectoryH5 log Warn and skip; they don't fail
// the whole load.

#pragma once

#include "ReaderInputManifest.h"
#include "../model/Conformation.h"
#include "../model/QtProtein.h"

#include <QString>
#include <memory>

namespace h5reader::io {

struct QtLoadResult {
    std::unique_ptr<h5reader::model::QtProtein> protein;
    std::unique_ptr<h5reader::model::Conformation> conformation;  // TrajectoryConformation or SingleConformation

    QString proteinId;
    bool ok = false;
    QString error;
    int decodeWarnings = 0;

    // The path that was opened (the trajectory.h5 file or the run dir). Used by
    // the window to locate sibling artifacts when no manifest is present (the
    // legacy bounded convention path).
    QString runPath;

    // The input-directory manifest that drove this load. manifest.ok is true
    // when the run was opened via a h5reader_manifest.toml; consumers should
    // prefer manifest paths over runPath-derived heuristics when ok.
    ReaderInputManifest manifest;
};

class QtProteinLoader {
public:
    // h5_path is the absolute path to trajectory.h5. The sidecar (5 NPYs +
    // extraction_manifest.json) and frame snapshot directory are found in the
    // same directory by bounded convention: per_frame_npys/ first, then npys/.
    static QtLoadResult Load(const QString& h5_path);

    // Trajectory load with explicit sidecar + per-frame snapshot dirs — the
    // path the manifest-driven flow uses (the manifest declares each
    // explicitly rather than letting the loader guess). `perFrameNpysDir`
    // may be empty (no per-frame snapshots; dense H5 still drives playback).
    static QtLoadResult LoadTrajectory(const QString& h5_path,
                                       const QString& topologySidecarDir,
                                       const QString& perFrameNpysDir);

    // run_dir is a single-pose run root holding the 5 sidecar NPYs +
    // extraction_manifest.json + the flat per-atom calculator NPYs (no
    // trajectory.h5). Builds the QtProtein spine and one snapshot.
    static QtLoadResult LoadPose(const QString& run_dir);

    // Manifest-driven entry. Dispatches by manifest.runKind:
    //   * Trajectory  -> LoadTrajectory(manifest.trajectoryH5, ...)
    //   * SinglePose  -> LoadPose(manifest.poseDir)
    //   * MutantPair  -> LoadPose(manifest.wtPoseDir)  [auto-open WT]
    // Returns ok=false with error set on any failure. The manifest is
    // copied into the result for downstream consumers (DFT path lookup,
    // ALA-switch action for mutant pairs).
    static QtLoadResult LoadFromManifest(const ReaderInputManifest& manifest);

    // Sniff a run path and dispatch:
    //   1. If a h5reader_manifest.toml is at the path's directory, parse it
    //      and call LoadFromManifest. Manifest is the documented contract.
    //   2. Otherwise fall back to the legacy bounded-convention check:
    //         - a dir containing trajectory.h5             -> Load(that h5)
    //         - a dataset root containing extract/trajectory.h5 -> Load(that h5)
    //         - the trajectory.h5 file itself              -> Load(path)
    //         - a single-pose run-root directory           -> LoadPose(path)
    //      The fallback emits a CRITICAL log naming the missing manifest so
    //      the user knows to generate one (see tools/generate_manifest.py).
    static QtLoadResult LoadRunPath(const QString& path);
};

}  // namespace h5reader::io
