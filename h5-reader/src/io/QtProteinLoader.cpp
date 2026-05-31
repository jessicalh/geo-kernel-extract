// QtProteinLoader implementation.
//
// String-discipline note: this file is the ONE place where the new
// format's string surfaces (atom_nom names, chain_id, mangled selection
// type_index names) cross the boundary from H5/NPY into the typed Qt
// model. The QtAtom typed substrate is populated from
// atoms_category_info.npy int8 columns at QtTopologySidecar::Load(); the
// parallel QtAtomNames is the projection layer for display. Past this
// file no code in the reader compares strings for chemistry.
//
// `.LGS` discipline (2026-05-31 SIMPLIFY): the loader has ONE entry
// path. The argument lands in CalcsetManifest::Load which either yields
// a typed manifest or a hard error. There is no fallback chain, no
// convention sniffing, no "if missing try these other names". Per
// feedback_no_file_discovery.

#include "QtProteinLoader.h"

#include "CalcsetManifest.h"
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

namespace {

// `<extraction_dir>/npys/frame_NNNNNN/` is the producer-documented
// per-frame snapshot layout. The reader honours this convention without
// enumerating it via glob; FrameNpyLoader receives one resolved frame
// dir per request. Returns empty if the convention dir does not exist.
QString DerivePerFrameNpysDir(const QString& extraction_dir) {
    if (extraction_dir.isEmpty()) return {};
    const QString cand = QStringLiteral("%1/npys").arg(extraction_dir);
    return QFileInfo(cand).isDir() ? cand : QString();
}

}  // namespace

QtLoadResult QtProteinLoader::LoadTrajectory(const QString& h5_path,
                                              const QString& extraction_dir,
                                              const QString& extraction_manifest_path) {
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
    if (!QFileInfo(extraction_dir).isDir()) {
        result.error = QStringLiteral("QtProteinLoader: extraction dir does not exist: %1").arg(extraction_dir);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                extraction_dir);
        return result;
    }

    // Display ID — the extraction dir's basename (typically a canonical run id).
    result.proteinId = QFileInfo(extraction_dir).fileName();

    const QString perFrameNpysDir = DerivePerFrameNpysDir(extraction_dir);
    qCInfo(cLoader).noquote() << "Loading" << result.proteinId << "| h5=" << h5_path
                              << "| extraction=" << extraction_dir
                              << "| extraction_manifest=" << extraction_manifest_path
                              << "| frame_npys=" << (perFrameNpysDir.isEmpty() ? "absent" : perFrameNpysDir);

    auto sidecar = QtTopologySidecar::Load(extraction_dir, extraction_manifest_path);
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

    // No producer-side protein_id consistency check here: the .LGS
    // protein_id is authoritative (spec/CALCSET_MANIFEST.md) and is
    // stamped onto the result in LoadFromManifest. Comparing the two
    // producer-written values (extraction_manifest.json vs H5 root)
    // against each other says nothing about which (if either) is right.

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

    auto conformation = std::make_unique<h5reader::model::TrajectoryConformation>(protein.get(), std::move(traj),
                                                                                  perFrameNpysDir);

    qCInfo(cLoader).noquote() << "Loaded" << result.proteinId << "| atoms=" << protein->atomCount()
                              << "| residues=" << protein->residueCount() << "| bonds=" << protein->bondCount()
                              << "| rings=" << protein->ringCount() << "| frames=" << conformation->frameCount()
                              << "| warnings=" << result.decodeWarnings;

    result.protein = std::move(protein);
    result.conformation = std::move(conformation);
    result.ok = true;
    return result;
}

QtLoadResult QtProteinLoader::LoadPose(const QString& pose_dir,
                                        const QString& extraction_manifest_path) {
    QtLoadResult result;
    result.runPath = pose_dir;
    QFileInfo fi(pose_dir);
    if (!fi.exists() || !fi.isDir()) {
        result.error = QStringLiteral("QtProteinLoader::LoadPose: pose dir does not exist: %1").arg(pose_dir);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                pose_dir);
        return result;
    }
    result.proteinId = QFileInfo(pose_dir).fileName();
    qCInfo(cLoader).noquote() << "Loading pose" << result.proteinId << "from" << pose_dir
                              << "| extraction_manifest=" << extraction_manifest_path;

    auto sidecar = QtTopologySidecar::Load(pose_dir, extraction_manifest_path);
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
    auto pose = FrameNpyLoader::LoadSnapshotDir(pose_dir, protein.get(), /*frameIndex=*/0, /*timePs=*/0.0);
    if (!pose) {
        result.error =
            QStringLiteral("QtProteinLoader::LoadPose: no per-atom NPYs loaded from %1").arg(pose_dir);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtProteinLoader"),
                                                result.error,
                                                pose_dir);
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

QtLoadResult QtProteinLoader::LoadFromManifest(const CalcsetManifest& manifest) {
    QtLoadResult result;
    result.runPath = manifest.calcset_root_abspath;
    result.manifest = manifest;

    qCInfo(cLoader).noquote() << "Loading via .LGS |" << manifest.lgs_path_abspath
                              << "| kind=" << CalcsetManifest::NameForKind(manifest.kind)
                              << "| protein_id=" << manifest.protein_id;

    switch (manifest.kind) {
        case CalcsetManifest::Kind::Trajectory: {
            const auto& t = *manifest.trajectory;
            QtLoadResult inner = LoadTrajectory(t.trajectory_h5_abspath,
                                                 t.extraction_dir_abspath,
                                                 t.extraction_manifest_abspath);
            // .LGS protein_id is authoritative; overrides producer-side
            // values stamped in LoadTrajectory from the extraction dir
            // basename and from the H5 root attribute.
            if (!manifest.protein_id.isEmpty() && inner.protein) {
                inner.proteinId = manifest.protein_id;
                inner.protein->proteinId_ = manifest.protein_id;
            }
            inner.manifest = manifest;
            inner.runPath = manifest.calcset_root_abspath;
            return inner;
        }
        case CalcsetManifest::Kind::SinglePose: {
            const auto& s = *manifest.single_pose;
            QtLoadResult inner = LoadPose(s.pose_dir_abspath, s.extraction_manifest_abspath);
            if (!manifest.protein_id.isEmpty() && inner.protein) {
                inner.proteinId = manifest.protein_id;
                inner.protein->proteinId_ = manifest.protein_id;
            }
            inner.manifest = manifest;
            inner.runPath = manifest.calcset_root_abspath;
            return inner;
        }
        case CalcsetManifest::Kind::MutantPair: {
            // Auto-open WT by recursively loading its `.LGS`. The parent
            // manifest (with both WT/ALA paths) is preserved on the
            // result so the window can offer the ALA-switch action.
            const auto& mp = *manifest.mutant_pair;
            qCInfo(cLoader).noquote() << "Mutant pair: opening WT |" << mp.wt_lgs_abspath
                                      << "| alternate ALA available |" << mp.ala_lgs_abspath;
            QString err;
            const auto wtManifest = CalcsetManifest::Load(mp.wt_lgs_abspath, &err);
            if (!wtManifest) {
                result.error = QStringLiteral(
                    "QtProteinLoader::LoadFromManifest: failed to load WT .LGS %1: %2")
                    .arg(mp.wt_lgs_abspath, err);
                return result;
            }
            QtLoadResult inner = LoadFromManifest(*wtManifest);
            // Preserve the PARENT mutant_pair manifest (not the nested
            // WT one) so the ALA-alternate action can find the ALA path.
            inner.manifest = manifest;
            inner.runPath = manifest.calcset_root_abspath;
            return inner;
        }
    }
    result.error = QStringLiteral("QtProteinLoader::LoadFromManifest: unhandled kind");
    return result;
}

QtLoadResult QtProteinLoader::LoadRunPath(const QString& path) {
    QtLoadResult result;
    result.runPath = path;
    QString err;
    auto manifest = CalcsetManifest::Load(path, &err);
    if (!manifest) {
        result.error = err.isEmpty()
            ? QStringLiteral("QtProteinLoader::LoadRunPath: failed to load .LGS at %1").arg(path)
            : err;
        // ErrorBus already reported in CalcsetManifest::Load.
        return result;
    }
    return LoadFromManifest(*manifest);
}

}  // namespace h5reader::io
