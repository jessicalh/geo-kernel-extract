// QtProteinLoader — top-level orchestrator: a run directory into a typed
// QtProtein + a Conformation (the base; concrete subclass depends on run shape).
//
// Two entry points, one per run shape (decided 2026-05-26):
//   Load(h5_path)     — a --trajectory run: sidecar + trajectory.h5 ->
//                       QtProtein + TrajectoryConformation (N frames, dense H5,
//                       lazy per-frame snapshots from per_frame_npys/).
//   LoadPose(run_dir) — a single-pose run (--orca / --mutant / --pdb): sidecar +
//                       the flat run-root per-atom NPYs -> QtProtein +
//                       SingleConformation (one frame, no H5).
//
// Failure modes log via ErrorBus and return ok=false with error populated.
// Per-TR-group reads inside QtTrajectoryH5 log Warn and skip; they don't fail
// the whole load.

#pragma once

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

    // The path that was opened (the trajectory.h5 file or the run dir). Lets the
    // window locate sibling artifacts by documented convention — e.g. the DFT
    // campaign at <dataset-root>/dft/jobs. A head-of-directory TOML descriptor
    // will formalise this later (#25); for now it is a bounded convention check.
    QString runPath;
};

class QtProteinLoader {
public:
    // h5_path is the absolute path to trajectory.h5. The sidecar (5 NPYs +
    // extraction_manifest.json) and per_frame_npys/ are found in the same
    // directory by canonical name — per feedback_no_file_discovery, no globbing.
    static QtLoadResult Load(const QString& h5_path);

    // run_dir is a single-pose run root holding the 5 sidecar NPYs +
    // extraction_manifest.json + the flat per-atom calculator NPYs (no
    // trajectory.h5). Builds the QtProtein spine and one snapshot.
    static QtLoadResult LoadPose(const QString& run_dir);

    // Sniff a run path and dispatch by documented convention (NOT globbing):
    //   - a directory containing trajectory.h5  -> Load(that h5)
    //   - the trajectory.h5 file itself          -> Load(path)
    //   - a single-pose run-root directory       -> LoadPose(path)
    // This is the "open a directory" entry point; main() and the window's
    // Open-Directory action both route through it.
    static QtLoadResult LoadRunPath(const QString& path);
};

}  // namespace h5reader::io
