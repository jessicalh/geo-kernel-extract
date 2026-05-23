// QtProteinLoader — top-level orchestrator: sidecar + trajectory.h5
// into typed QtProtein + QtConformation.
//
// Sequence:
//   1. Resolve sidecar dir = dirname(h5_path)
//   2. QtTopologySidecar::Load(sidecar_dir) → 5 NPYs + manifest
//      decoded into typed Qt model objects
//   3. Construct QtTrajectoryH5(h5_path) — eager-load all 56 TR groups
//   4. Cross-check atomCount() / proteinId between sidecar and H5
//   5. Build QtProtein with atoms/names/residues/topology
//   6. Build QtConformation with h5 + protein back-ptr
//   7. Return QtLoadResult
//
// Failure modes log via ErrorBus and return ok=false with error
// populated. Per-TR-group reads inside QtTrajectoryH5 log Warn and
// skip; they don't fail the whole load.

#pragma once

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QString>
#include <memory>

namespace h5reader::io {

struct QtLoadResult {
    std::unique_ptr<h5reader::model::QtProtein> protein;
    std::unique_ptr<h5reader::model::QtConformation> conformation;

    QString proteinId;
    bool ok = false;
    QString error;
    int decodeWarnings = 0;
};

class QtProteinLoader {
public:
    // h5_path is the absolute path to trajectory.h5. Sidecar (the 5
    // NPYs + extraction_manifest.json) is found in the same directory
    // by canonical filename — per `feedback_no_file_discovery`, no
    // globbing.
    static QtLoadResult Load(const QString& h5_path);
};

}  // namespace h5reader::io
