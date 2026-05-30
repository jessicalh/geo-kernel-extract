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
    // Legacy convention-based entry — derive sidecar and per_frame_npys
    // dirs by name from the H5's parent, then call the explicit-path
    // form. Manifest-driven loads bypass this entirely.
    QFileInfo fi(h5_path);
    if (!fi.exists()) {
        QtLoadResult result;
        result.runPath = h5_path;
        result.error = QStringLiteral("QtProteinLoader: H5 path does not exist: %1").arg(h5_path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                h5_path);
        return result;
    }
    const QString sidecar_dir = fi.absolutePath();
    QString perFrameNpysDir;
    for (const QString& name : {QStringLiteral("per_frame_npys"), QStringLiteral("npys")}) {
        const QString candidate = QStringLiteral("%1/%2").arg(sidecar_dir, name);
        if (QFileInfo::exists(candidate)) {
            perFrameNpysDir = candidate;
            break;
        }
    }
    return LoadTrajectory(h5_path, sidecar_dir, perFrameNpysDir);
}

QtLoadResult QtProteinLoader::LoadTrajectory(const QString& h5_path,
                                              const QString& topologySidecarDir,
                                              const QString& perFrameNpysDir) {
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
    if (!QFileInfo(topologySidecarDir).isDir()) {
        result.error = QStringLiteral("QtProteinLoader: sidecar dir does not exist: %1").arg(topologySidecarDir);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                topologySidecarDir);
        return result;
    }

    // Display ID — the sidecar dir's basename (typically a canonical run id).
    result.proteinId = QFileInfo(topologySidecarDir).fileName();

    qCInfo(cLoader).noquote() << "Loading" << result.proteinId << "| h5=" << h5_path
                              << "| sidecar=" << topologySidecarDir
                              << "| frame_npys=" << (perFrameNpysDir.isEmpty() ? "absent" : perFrameNpysDir);

    auto sidecar = QtTopologySidecar::Load(topologySidecarDir);
    if (!sidecar.ok) {
        result.error = sidecar.error;
        return result;
    }
    result.decodeWarnings += sidecar.warningCount;

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

    // Per-frame snapshot dir is whatever the caller passed — empty means
    // no per-frame detail (dense H5 still drives playback).
    QString resolvedPerFrame = perFrameNpysDir;
    if (!resolvedPerFrame.isEmpty() && !QFileInfo(resolvedPerFrame).isDir()) {
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                QStringLiteral("QtProteinLoader"),
                                                QStringLiteral("per_frame_npys_dir does not exist: %1")
                                                    .arg(resolvedPerFrame),
                                                resolvedPerFrame);
        resolvedPerFrame.clear();
        ++result.decodeWarnings;
    }
    auto conformation = std::make_unique<h5reader::model::TrajectoryConformation>(protein.get(), std::move(traj),
                                                                                  resolvedPerFrame);

    qCInfo(cLoader).noquote() << "Loaded" << result.proteinId << "| atoms=" << protein->atomCount()
                              << "| residues=" << protein->residueCount() << "| bonds=" << protein->bondCount()
                              << "| rings=" << protein->ringCount() << "| frames=" << conformation->frameCount()
                              << "| warnings=" << result.decodeWarnings;

    result.protein = std::move(protein);
    result.conformation = std::move(conformation);
    result.ok = true;
    return result;
}

QtLoadResult QtProteinLoader::LoadFromManifest(const ReaderInputManifest& manifest) {
    QtLoadResult result;
    result.runPath = manifest.rootDir;
    result.manifest = manifest;

    if (!manifest.ok) {
        result.error = QStringLiteral("QtProteinLoader::LoadFromManifest: manifest not valid: %1")
                          .arg(manifest.error);
        return result;
    }

    qCInfo(cLoader).noquote() << "Loading via manifest |" << manifest.manifestPath
                              << "| run_kind=" << ReaderInputManifest::NameForRunKind(manifest.runKind)
                              << "| protein_id=" << manifest.proteinId;

    switch (manifest.runKind) {
        case ReaderInputManifest::RunKind::Trajectory: {
            QtLoadResult inner = LoadTrajectory(manifest.trajectoryH5,
                                                 manifest.topologySidecarDir,
                                                 manifest.perFrameNpysDir);
            inner.manifest = manifest;
            inner.runPath = manifest.rootDir;
            return inner;
        }
        case ReaderInputManifest::RunKind::SinglePose: {
            QtLoadResult inner = LoadPose(manifest.poseDir);
            inner.manifest = manifest;
            inner.runPath = manifest.rootDir;
            return inner;
        }
        case ReaderInputManifest::RunKind::MutantPair: {
            // Auto-open WT — per project decision (the ALA pose is exposed
            // through ReaderInputManifest::alternatePoseDir() for the
            // toolbar/menu switch action).
            QtLoadResult inner = LoadPose(manifest.wtPoseDir);
            inner.manifest = manifest;
            inner.runPath = manifest.rootDir;
            qCInfo(cLoader).noquote() << "Mutant pair: opened WT |" << manifest.wtPoseDir
                                      << "| alternate ALA available |" << manifest.alaPoseDir;
            return inner;
        }
    }
    result.error = QStringLiteral("QtProteinLoader::LoadFromManifest: unhandled run_kind");
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

    // Manifest-driven entry — preferred. The directory the user opens is
    // the "run root"; the manifest lives at the conventional file name
    // at that root.
    const QString manifestDir = fi.isDir() ? path : fi.absolutePath();
    if (ReaderInputManifest::ExistsAt(manifestDir)) {
        const auto manifest = ReaderInputManifest::Load(manifestDir);
        if (!manifest.ok) {
            QtLoadResult result;
            result.runPath = path;
            result.error = manifest.error;
            return result;
        }
        return LoadFromManifest(manifest);
    }

    // Legacy bounded-convention fallback. Log CRITICAL so the user knows
    // the manifest is missing and that this path will become a hard
    // error in a future release. See spec/INPUT_DIRECTORY.md.
    qCCritical(cLoader).noquote()
        << "no h5reader_manifest.toml at" << manifestDir
        << "— falling back to legacy convention; generate one via tools/generate_manifest.py";

    if (fi.isDir()) {
        const QDir dir(path);
        const QString directH5 = QStringLiteral("%1/trajectory.h5").arg(dir.absolutePath());
        if (QFileInfo::exists(directH5))
            return Load(directH5);
        const QString extractH5 = QStringLiteral("%1/extract/trajectory.h5").arg(dir.absolutePath());
        if (QFileInfo::exists(extractH5))
            return Load(extractH5);
        return LoadPose(path);
    }
    return Load(path);
}

}  // namespace h5reader::io
