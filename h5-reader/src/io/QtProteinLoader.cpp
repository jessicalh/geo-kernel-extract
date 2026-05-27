// QtProteinLoader implementation.
//
// String-discipline note: this file is the ONE place where the
// new format's string surfaces (atom_nom names, chain_id, mangled
// selection type_index names) cross the boundary from H5/NPY into
// the typed Qt model. The QtAtom typed substrate is populated from
// atoms_category_info.npy int8 columns at QtTopologySidecar::Load();
// the parallel QtAtomNames is the projection layer for display.
// Past this file no code in the reader compares strings for chemistry.

#include "QtProteinLoader.h"

#include "FrameNpyLoader.h"
#include "QtTopologySidecar.h"
#include "QtTrajectoryH5.h"

#include "../diagnostics/ErrorBus.h"
#include "../model/QtTopology.h"
#include "../model/SingleConformation.h"
#include "../model/TrajectoryConformation.h"

#include <QDir>
#include <QFileInfo>
#include <QLoggingCategory>
#include <stdexcept>

Q_LOGGING_CATEGORY(cLoader, "h5reader.loader")

namespace h5reader::io {

QtLoadResult QtProteinLoader::Load(const QString& h5_path) {
    QtLoadResult result;
    result.runPath = h5_path;
    QFileInfo fi(h5_path);
    if (!fi.exists()) {
        result.error = QStringLiteral("QtProteinLoader: H5 path does not exist: %1").arg(h5_path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                h5_path);
        return result;
    }

    // Derive a display ID from the H5 file's parent dir basename — the
    // sidecar manifest also carries protein_id, but the dir name is
    // typically the canonical run identifier (e.g.
    // "run_01_baseline_2026-05-21_sha7d8dbe9").
    const QString sidecar_dir = fi.absolutePath();
    result.proteinId = QFileInfo(sidecar_dir).fileName();

    qCInfo(cLoader).noquote() << "Loading" << result.proteinId << "from" << h5_path;

    // Step 1+2: sidecar
    auto sidecar = QtTopologySidecar::Load(sidecar_dir);
    if (!sidecar.ok) {
        result.error = sidecar.error;
        return result;
    }
    result.decodeWarnings += sidecar.warningCount;

    // Step 3: trajectory.h5 (eager-load all TR groups; throws on
    // structural failure).
    std::unique_ptr<QtTrajectoryH5> traj;
    try {
        traj = std::make_unique<QtTrajectoryH5>(h5_path);
    } catch (const std::exception& e) {
        result.error = QStringLiteral("QtTrajectoryH5: %1").arg(QString::fromUtf8(e.what()));
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                h5_path);
        return result;
    }

    // Step 4: cross-check atom count.
    if (sidecar.atoms.size() != traj->atomCount()) {
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                QStringLiteral("QtProteinLoader"),
                                                QStringLiteral("sidecar atom count %1 != H5 /atoms count %2; "
                                                               "typed substrate may not align with per-atom TRs")
                                                    .arg(sidecar.atoms.size())
                                                    .arg(traj->atomCount()),
                                                h5_path);
        ++result.decodeWarnings;
    }

    // Optional cross-check: manifest.protein_id vs H5 root attr
    // protein_id. They typically agree; log Warn on mismatch but do
    // not fail.
    if (!sidecar.manifest.proteinId.isEmpty() && !traj->proteinId().isEmpty()
        && sidecar.manifest.proteinId != traj->proteinId()) {
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                QStringLiteral("QtProteinLoader"),
                                                QStringLiteral("manifest.protein_id (%1) disagrees with "
                                                               "H5 root protein_id (%2)")
                                                    .arg(sidecar.manifest.proteinId, traj->proteinId()),
                                                h5_path);
        ++result.decodeWarnings;
    }

    // Step 5: build QtProtein.
    auto protein = std::make_unique<h5reader::model::QtProtein>();
    protein->proteinId_ = result.proteinId;
    protein->atoms_ = std::move(sidecar.atoms);
    protein->atomNames_ = std::move(sidecar.atomNames);
    protein->residueNames_ = std::move(sidecar.residueNames);
    protein->residues_ = std::move(sidecar.residues);
    protein->topology_ = std::make_unique<h5reader::model::QtTopology>(protein->atoms_.size(),
                                                                       std::move(sidecar.bonds),
                                                                       std::move(sidecar.rings),
                                                                       std::move(sidecar.ringMemberships),
                                                                       sidecar.aromaticRingCount,
                                                                       sidecar.saturatedRingCount);

    // Step 6: build the trajectory conformation. The per-frame snapshot detail
    // is lazily loaded from per_frame_npys/ (a sibling of trajectory.h5, present
    // only when the run used --emit-frame-npys); absent → no detail, playback
    // still runs off the dense H5.
    QString perFrameNpysDir;
    {
        const QString candidate = QStringLiteral("%1/per_frame_npys").arg(sidecar_dir);
        if (QFileInfo::exists(candidate))
            perFrameNpysDir = candidate;
    }
    auto conformation = std::make_unique<h5reader::model::TrajectoryConformation>(protein.get(), std::move(traj),
                                                                                  perFrameNpysDir);

    qCInfo(cLoader).noquote() << "Loaded" << result.proteinId << "| atoms=" << protein->atomCount()
                              << "| residues=" << protein->residueCount() << "| bonds=" << protein->bondCount()
                              << "| rings=" << protein->ringCount() << "| frames=" << conformation->frameCount()
                              << "| per_frame_npys=" << (perFrameNpysDir.isEmpty() ? "absent" : "present")
                              << "| warnings=" << result.decodeWarnings;

    result.protein = std::move(protein);
    result.conformation = std::move(conformation);
    result.ok = true;
    return result;
}

QtLoadResult QtProteinLoader::LoadPose(const QString& run_dir) {
    QtLoadResult result;
    result.runPath = run_dir;
    QFileInfo fi(run_dir);
    if (!fi.exists() || !fi.isDir()) {
        result.error = QStringLiteral("QtProteinLoader::LoadPose: run dir does not exist: %1").arg(run_dir);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                run_dir);
        return result;
    }
    result.proteinId = QFileInfo(run_dir).fileName();
    qCInfo(cLoader).noquote() << "Loading pose" << result.proteinId << "from" << run_dir;

    // Steps 1-2 + 5: sidecar -> QtProtein (no H5). Duplicated from Load() per
    // the project's duplication-over-chaining preference; both are friends of
    // QtProtein.
    auto sidecar = QtTopologySidecar::Load(run_dir);
    if (!sidecar.ok) {
        result.error = sidecar.error;
        return result;
    }
    result.decodeWarnings += sidecar.warningCount;

    auto protein = std::make_unique<h5reader::model::QtProtein>();
    protein->proteinId_ = result.proteinId;
    protein->atoms_ = std::move(sidecar.atoms);
    protein->atomNames_ = std::move(sidecar.atomNames);
    protein->residueNames_ = std::move(sidecar.residueNames);
    protein->residues_ = std::move(sidecar.residues);
    protein->topology_ = std::make_unique<h5reader::model::QtTopology>(protein->atoms_.size(),
                                                                       std::move(sidecar.bonds),
                                                                       std::move(sidecar.rings),
                                                                       std::move(sidecar.ringMemberships),
                                                                       sidecar.aromaticRingCount,
                                                                       sidecar.saturatedRingCount);

    // The single pose's per-atom calculator NPYs sit flat in the run root —
    // the same FrameNpyLoader that reads a trajectory frame dir reads it.
    auto pose = FrameNpyLoader::LoadSnapshotDir(run_dir, protein.get(), /*frameIndex=*/0, /*timePs=*/0.0);
    if (!pose) {
        result.error =
            QStringLiteral("QtProteinLoader::LoadPose: no per-atom NPYs loaded from %1").arg(run_dir);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                run_dir);
        return result;
    }

    auto conformation = std::make_unique<h5reader::model::SingleConformation>(protein.get(), std::move(pose));

    qCInfo(cLoader).noquote() << "Loaded pose" << result.proteinId << "| atoms=" << protein->atomCount()
                              << "| residues=" << protein->residueCount() << "| bonds=" << protein->bondCount()
                              << "| rings=" << protein->ringCount() << "| warnings=" << result.decodeWarnings;

    result.protein = std::move(protein);
    result.conformation = std::move(conformation);
    result.ok = true;
    return result;
}

QtLoadResult QtProteinLoader::LoadRunPath(const QString& path) {
    QFileInfo fi(path);
    if (fi.isDir()) {
        // A directory: a trajectory run has trajectory.h5 beside its sidecar;
        // a single-pose run root does not. This is a documented-convention
        // presence check (the same QFileInfo::exists Load() uses for the
        // sidecar), not file discovery / globbing.
        const QString h5 = QStringLiteral("%1/trajectory.h5").arg(QDir(path).absolutePath());
        if (QFileInfo::exists(h5))
            return Load(h5);
        return LoadPose(path);
    }
    // A file path: treat it as the trajectory.h5 to load.
    return Load(path);
}

}  // namespace h5reader::io
