#include "RunData.h"

#include "ChargeStore.h"

#include "../io/DftShieldingLoader.h"
#include "../io/QtProteinLoader.h"
#include "../model/QtBond.h"

#include <QLoggingCategory>

#include <algorithm>
#include <cmath>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cRun, "h5reader.rediscover.run")

bool proteinLooksWhole(const RunData& run, QString* err_out) {
    if (!run.protein || !run.conformation) return true;
    const model::QtTopology& topo = run.protein->topology();
    constexpr double kWrappedBondThresholdA = 5.0;
    for (std::size_t frame = 0; frame < run.conformation->frameCount(); ++frame) {
        for (std::size_t bi = 0; bi < topo.bondCount(); ++bi) {
            const model::QtBond& b = topo.bondAt(bi);
            if (b.atomIndexA < 0 || b.atomIndexB < 0) continue;
            const model::Vec3 a = run.conformation->atomPosition(frame, static_cast<std::size_t>(b.atomIndexA));
            const model::Vec3 c = run.conformation->atomPosition(frame, static_cast<std::size_t>(b.atomIndexB));
            const double d = (c - a).norm();
            if (d > kWrappedBondThresholdA) {
                if (err_out) {
                    *err_out = QStringLiteral("protein appears wrapped: frame %1 bond %2 length %3 A "
                                              "(PBC mode is None; pbc_whole must run upstream)")
                                   .arg(frame)
                                   .arg(b.bondIndex)
                                   .arg(d);
                }
                return false;
            }
        }
    }
    return true;
}
}

std::optional<FrameMap> FrameMap::Build(const model::TrajectoryConformation& traj,
                                        const DftFrameSet& dft,
                                        const QString& frame_index_basis,
                                        QString* err_out) {
    // Fail loud on basis mismatch — a different basis would silently
    // mis-key the DFT lookup.
    if (frame_index_basis != QStringLiteral("trr_frame_index")) {
        if (err_out)
            *err_out = QStringLiteral("unexpected frame_index_basis '%1' (expected 'trr_frame_index')")
                           .arg(frame_index_basis);
        return std::nullopt;
    }

    FrameMap fm;
    const std::size_t n = traj.frameCount();
    fm.originalByRow_.resize(n);
    for (std::size_t row = 0; row < n; ++row) {
        const std::size_t orig = traj.originalFrameIndex(row);
        fm.originalByRow_[row] = orig;
        if (dft.Has(orig)) fm.dftRows_.push_back(row);
    }
    std::sort(fm.dftRows_.begin(), fm.dftRows_.end());

    if (fm.originalByRow_.size() != n) {
        if (err_out) *err_out = QStringLiteral("frame map row count disagreement");
        return std::nullopt;
    }
    return fm;
}

std::optional<RunData> RunLoader::Load(const QString& calcset_path, QString* err_out) {
    // 1. Protein + conformation via the existing loader (resolves the `.LGS`,
    //    sidecar, trajectory.h5). No file discovery here.
    io::QtLoadResult loaded = io::QtProteinLoader::LoadRunPath(calcset_path);
    if (!loaded.ok) {
        if (err_out) *err_out = QStringLiteral("protein load failed: %1").arg(loaded.error);
        return std::nullopt;
    }
    if (loaded.manifest.kind != io::CalcsetManifest::Kind::Trajectory || !loaded.manifest.trajectory) {
        if (err_out) *err_out = QStringLiteral("rediscover requires a trajectory calcset");
        return std::nullopt;
    }
    const model::TrajectoryConformation* traj =
        loaded.conformation ? loaded.conformation->asTrajectory() : nullptr;
    if (!traj) {
        if (err_out) *err_out = QStringLiteral("loaded conformation is not H5-backed trajectory");
        return std::nullopt;
    }

    RunData run;
    run.protein = std::move(loaded.protein);
    run.conformation = std::move(loaded.conformation);
    run.manifest = loaded.manifest;

    if (run.manifest.trajectory) {
        QString chargeErr;
        if (!LoadFf14sbChargesFromTopol(run.manifest.trajectory->topology_top_abspath,
                                        *run.protein, &chargeErr)) {
            if (err_out) *err_out = chargeErr;
            return std::nullopt;
        }
        qCInfo(cRun).noquote() << "FF14SB charges loaded from"
                               << run.manifest.trajectory->topology_top_abspath;
    }

    QString pbcErr;
    if (!proteinLooksWhole(run, &pbcErr)) {
        if (err_out) *err_out = pbcErr;
        return std::nullopt;
    }

    // 2. DFT frames via DftShieldingLoader, keyed by manifest frame_index.
    //    A partial campaign is fine — gaps are honest, not failures.
    std::size_t dftLoaded = 0, dftGap = 0;
    if (run.manifest.dft) {
        for (const io::DftFrame& f : run.manifest.dft->frames) {
            auto frame = io::DftShieldingLoader::LoadAndValidate(f.meta_json_abspath, run.protein.get());
            if (frame) {
                run.dft.Insert(static_cast<std::size_t>(f.frame_index), std::move(frame));
                ++dftLoaded;
            } else {
                ++dftGap;
            }
        }
    }
    qCInfo(cRun).noquote() << "DFT frames loaded=" << dftLoaded << "gaps=" << dftGap;

    // 3. Frame map (basis check + row→original + dft rows).
    QString fmErr;
    auto fm = FrameMap::Build(*run.trajectory(), run.dft,
                              run.manifest.trajectory->frame_index_basis, &fmErr);
    if (!fm) {
        if (err_out) *err_out = fmErr;
        return std::nullopt;
    }
    run.frameMap = std::move(*fm);

    qCInfo(cRun).noquote() << "RunData ready | atoms=" << run.protein->atomCount()
                           << "| frames=" << run.trajectory()->frameCount()
                           << "| dft_rows=" << run.frameMap.dftRows().size();
    return run;
}

}  // namespace h5reader::rediscover
