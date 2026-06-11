#include "RunData.h"

#include "CanonicalSpineGuard.h"
#include "ChargeStore.h"
#include "RowDesignCatalogCoverage.h"

#include "../io/DftShieldingLoader.h"
#include "../io/QtNpyReader.h"
#include "../io/QtProteinLoader.h"
#include "../model/QtBond.h"

#include <QDir>
#include <QFileInfo>
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

QString qsv(std::string_view s) {
    return QString::fromUtf8(s.data(), static_cast<qsizetype>(s.size()));
}

QString fieldStem(const io::FieldSpec& spec) { return qsv(spec.stem); }

bool structuredField(const io::FieldSpec& spec) {
    return spec.kind == io::FieldKind::AtomsCategoryInfo
           || spec.group == io::FieldGroup::Topology;
}

bool dftField(io::FieldKind kind) {
    return kind == io::FieldKind::OrcaTotal
           || kind == io::FieldKind::OrcaDiamagnetic
           || kind == io::FieldKind::OrcaParamagnetic;
}

void normalizeProducerArray(const io::FieldSpec& spec, io::QtNpyReader::WidenedArray* a) {
    if (!a || !a->ok) return;
    if (spec.axis == io::NativeAxis::Protein
        && spec.cols > 1
        && a->cols == 1
        && a->rows == static_cast<std::size_t>(spec.cols)) {
        a->rows = 1;
        a->cols = static_cast<std::size_t>(spec.cols);
    }
}

QString trajectoryFrameNpyPath(const RunData& run, std::size_t row, const io::FieldSpec& spec) {
    if (!run.manifest.trajectory) return {};
    const QString npysDir = QDir(run.manifest.trajectory->extraction_dir_abspath)
                                .filePath(QStringLiteral("npys"));
    const std::size_t original = run.frameMap.originalIndex(row);
    const QString frameDir = QDir(npysDir).filePath(
        QStringLiteral("frame_%1").arg(original, 6, 10, QLatin1Char('0')));
    return QDir(frameDir).filePath(fieldStem(spec) + QStringLiteral(".npy"));
}

bool validateFixedNativeRows(const io::FieldSpec& spec,
                             const model::QtProtein& protein,
                             std::size_t rows,
                             const QString& path,
                             QString* err_out) {
    std::optional<std::size_t> expected;
    switch (spec.axis) {
    case io::NativeAxis::Atom:
        expected = protein.atomCount();
        break;
    case io::NativeAxis::Residue:
        expected = protein.residueCount();
        break;
    case io::NativeAxis::Protein:
        expected = 1;
        break;
    case io::NativeAxis::AromaticRing:
        expected = protein.topology().aromaticRingCount();
        break;
    case io::NativeAxis::SaturatedRing:
        expected = protein.topology().saturatedRingCount();
        break;
    case io::NativeAxis::Ring:
        expected = protein.topology().ringCount();
        break;
    case io::NativeAxis::RingMembership:
        expected = protein.topology().ringMembershipCount();
        break;
    default:
        break;
    }
    if (expected && rows != *expected) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 native rows; %3 expects %4")
                           .arg(path)
                           .arg(rows)
                           .arg(fieldStem(spec))
                           .arg(*expected);
        }
        return false;
    }
    return true;
}

bool loadTrajectoryProducerArrays(RunData& run, QString* err_out) {
    if (!run.manifest.trajectory || !run.protein) return true;
    const std::size_t frames = run.frameMap.frameCount();
    for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
        if (!spec || structuredField(*spec) || dftField(spec->kind)) continue;

        StaticNpyArray resident;
        resident.stem = fieldStem(*spec);
        resident.path = QDir(run.manifest.trajectory->extraction_dir_abspath)
                            .filePath(QStringLiteral("npys"));
        resident.frameVarying = true;
        resident.frameCount = frames;
        resident.frameRowOffsets.reserve(frames);
        resident.frameRowCounts.reserve(frames);
        const bool keepFloat = spec->kind == io::FieldKind::AIMNet2Aim;
        std::size_t cols = 0;
        std::size_t totalRows = 0;

        for (std::size_t row = 0; row < frames; ++row) {
            const QString path = trajectoryFrameNpyPath(run, row, *spec);
            io::QtNpyReader::WidenedArray arr = io::QtNpyReader::ReadArrayWidened(path);
            if (!arr.ok) {
                if (err_out) *err_out = arr.error;
                return false;
            }
            normalizeProducerArray(*spec, &arr);
            if (spec->cols != -1 && arr.cols != static_cast<std::size_t>(spec->cols)) {
                if (err_out) {
                    *err_out = QStringLiteral("%1 has %2 columns; catalog %3 expects %4")
                                   .arg(path)
                                   .arg(arr.cols)
                                   .arg(fieldStem(*spec))
                                   .arg(spec->cols);
                }
                return false;
            }
            if (!validateFixedNativeRows(*spec, *run.protein, arr.rows, path, err_out))
                return false;
            if (row == 0) {
                cols = arr.cols;
                resident.dtype_descr = QString::fromStdString(arr.descr);
                resident.atomsPerFrame = arr.rows;
            } else if (arr.cols != cols) {
                if (err_out) {
                    *err_out = QStringLiteral("%1 column drift at frame row %2: expected %3, saw %4")
                                   .arg(fieldStem(*spec))
                                   .arg(row)
                                   .arg(cols)
                                   .arg(arr.cols);
                }
                return false;
            }

            resident.frameRowOffsets.push_back(totalRows);
            resident.frameRowCounts.push_back(arr.rows);
            totalRows += arr.rows;
            if (keepFloat) {
                resident.floatValues.reserve(resident.floatValues.size() + arr.data.size());
                for (double v : arr.data)
                    resident.floatValues.push_back(static_cast<float>(v));
            } else {
                resident.values.insert(resident.values.end(), arr.data.begin(), arr.data.end());
            }
        }
        resident.rows = totalRows;
        resident.cols = cols;
        run.producerArrays[static_cast<int>(spec->kind)] = std::move(resident);
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

FrameMap FrameMap::Static(std::size_t originalIndex, bool hasDft) {
    FrameMap fm;
    fm.originalByRow_.push_back(originalIndex);
    if (hasDft) fm.dftRows_.push_back(0);
    return fm;
}

std::optional<RunData> RunLoader::Load(const QString& calcset_path, QString* err_out) {
    if (!ValidateCanonical1p9jRoot(calcset_path, err_out)) return std::nullopt;

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

    QString residentErr;
    if (!loadTrajectoryProducerArrays(run, &residentErr)) {
        if (err_out) *err_out = residentErr;
        return std::nullopt;
    }

    qCInfo(cRun).noquote() << "RunData ready | atoms=" << run.protein->atomCount()
                           << "| frames=" << run.trajectory()->frameCount()
                           << "| dft_rows=" << run.frameMap.dftRows().size()
                           << "| producer_arrays=" << run.producerArrays.size();
    return run;
}

}  // namespace h5reader::rediscover
