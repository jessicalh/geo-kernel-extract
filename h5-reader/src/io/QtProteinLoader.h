// QtProteinLoader — top-level orchestrator: a calcset directory into a typed
// QtProtein + a Conformation (the base; concrete subclass depends on run shape).
//
// Reads the calcset's `.LGS` (CalcsetManifest, spec/CALCSET_MANIFEST.md),
// then dispatches by `kind`:
//   * Trajectory  -> sidecar + trajectory.h5 -> QtProtein + TrajectoryConformation
//   * SinglePose  -> sidecar -> QtProtein + SingleConformation
//   * MutantPair  -> opens the WT side (the ALA `.LGS` is exposed for a
//                    separate-process "Open mutant alternate" action).
//
// There is exactly ONE entry point — LoadRunPath. The pre-`.LGS` legacy
// convention-based fallback (Load(h5_path), bounded-convention sniffing)
// was deleted along with ReaderInputManifest in the 2026-05-31 SIMPLIFY
// pass. Missing/malformed `.LGS` is a hard error.
//
// Failure modes log via ErrorBus and return ok=false with error populated.

#pragma once

#include "CalcsetManifest.h"
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

    // The path that was opened (the calcset directory). Used by the window
    // for status bar / recent-files / window title display.
    QString runPath;

    // The .LGS that drove this load. ok-set sub-blocks expose the typed
    // pointers downstream consumers need (dft frames for the strip chart,
    // alternate-LGS for mutant pair switch). When ok=false on the result,
    // manifest will be a default-constructed `nullopt`-everywhere instance.
    CalcsetManifest manifest;
};

class QtProteinLoader {
public:
    // Trajectory load — used by LoadFromManifest for the Trajectory kind.
    //   h5_path: the trajectory.h5 file
    //   extraction_dir: the directory holding the 5-NPY topology sidecar
    //                   AND the per-frame `npys/frame_NNNNNN/` snapshots
    //                   (producer convention; not enumerated by globbing
    //                   in this function, only at scan time)
    //   extraction_manifest_path: explicit path to extraction_manifest.json
    //                   (the `.LGS` carries it as
    //                   trajectory.extraction_manifest)
    static QtLoadResult LoadTrajectory(const QString& h5_path,
                                       const QString& extraction_dir,
                                       const QString& extraction_manifest_path);

    // Single-pose load — used by LoadFromManifest for SinglePose / the WT
    // half of MutantPair.
    //   pose_dir: directory holding the 5-NPY topology sidecar AND the
    //             flat per-atom calculator NPYs (no trajectory.h5)
    //   extraction_manifest_path: explicit path to extraction_manifest.json
    static QtLoadResult LoadPose(const QString& pose_dir,
                                 const QString& extraction_manifest_path);

    // Dispatch by `manifest.kind`:
    //   * Trajectory  -> LoadTrajectory(manifest.trajectory)
    //   * SinglePose  -> LoadPose(manifest.single_pose)
    //   * MutantPair  -> recursively Load() the WT `.LGS`; the parent
    //                    manifest is preserved on the result so the
    //                    window can offer "Open mutant alternate (ALA)…".
    // Returns ok=false with error set on any failure. The manifest is
    // copied into the result for downstream consumers.
    static QtLoadResult LoadFromManifest(const CalcsetManifest& manifest);

    // The single user-facing entry point.
    // Argument is either a calcset directory holding the single `.LGS`,
    // or a `.LGS` file directly. CalcsetManifest::Load resolves and
    // dispatches; any failure stops the load with a clear error.
    static QtLoadResult LoadRunPath(const QString& path);
};

}  // namespace h5reader::io
