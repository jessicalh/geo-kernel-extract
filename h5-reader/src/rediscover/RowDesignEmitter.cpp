#include "RowDesignEmitter.h"

#include "Catalog.h"
#include "ExtractionSupport.h"
#include "RamaRegion.h"
#include "RelationshipEngine.h"
#include "SphericalBasis.h"
#include "Verbs.h"

#include "../io/QtFieldCatalog.gen.h"
#include "../io/QtNpyReader.h"
#include "../io/QtTrajectoryH5.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtPerResidueBuffers.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtResultBlocks.h"
#include "../model/QtRing.h"
#include "../model/QtTopology.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <QDir>
#include <QFileInfo>

namespace h5reader::rediscover {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kRingAzimuthNormFloor = 1e-10;
constexpr std::size_t kRingTensorSidecarDims = 162;

static_assert(static_cast<int>(model::RingTypeIndex::PheBenzene) == 0, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::TyrPhenol) == 1, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::TrpBenzene) == 2, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::TrpPyrrole) == 3, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::TrpPerimeter) == 4, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::HisImidazole) == 5, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::HidImidazole) == 6, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::HieImidazole) == 7, "ring intensity slot mismatch");
static_assert(static_cast<int>(model::RingTypeIndex::ProPyrrolidine)
              == model::kAromaticRingTypeCount,
              "Pro must remain outside the per-aromatic ring arrays");

constexpr std::array<double, model::kAromaticRingTypeCount> kRingLiteratureIntensity = {
    -12.0,   // PheBenzene
    -11.28,  // TyrPhenol
    -12.48,  // TrpBenzene
    -6.72,   // TrpPyrrole
    -19.2,   // Trp indole perimeter
    -5.16,   // HisImidazole
    -5.16,   // HidImidazole
    -5.16,   // HieImidazole
};

struct RingArrayRequirement {
    ArrayId id;
    bool required = true;
};

struct TrajectoryNpyRequirement {
    ArrayId id;
    bool required = true;
};

constexpr RingArrayRequirement kRequiredRingArrays[] = {
    {ArrayId::BSPerTypeT0, true},
    {ArrayId::BSPerTypeT1, true},
    {ArrayId::BSPerTypeT2, true},
    {ArrayId::HMPerTypeT0, true},
    {ArrayId::HMPerTypeT1, true},
    {ArrayId::HMPerTypeT2, true},
};

constexpr TrajectoryNpyRequirement kTrajectorySourceArrays[] = {
    {ArrayId::MopacCharge, true},
    {ArrayId::MopacCoulombEfield, true},
    {ArrayId::Aimnet2Efg, true},
    {ArrayId::EeqCoordinationNumber, true},
    {ArrayId::LarsenHBondShielding, true},
    {ArrayId::Pyramidalization, true},
    {ArrayId::AromaticChi2, true},
    {ArrayId::PuckerQ, true},
    {ArrayId::PuckerTheta, true},
};

QString fieldStem(const io::FieldSpec& spec) {
    return QString::fromUtf8(spec.stem.data(), static_cast<qsizetype>(spec.stem.size()));
}

std::optional<io::FieldKind> producerKind(ArrayId id, QString* err_out) {
    const std::optional<io::FieldKind> kind = ProducerFieldFor(id);
    if (!kind && err_out) {
        *err_out = QStringLiteral("row_design ArrayId has no producer FieldKind");
    }
    return kind;
}

bool validateResidentArray(const RunData& run,
                           const RingArrayRequirement& req,
                           std::size_t atomCount,
                           QString* err_out) {
    const std::optional<io::FieldKind> kind = producerKind(req.id, err_out);
    if (!kind) return false;
    const io::FieldSpec& spec = io::FieldSpecFor(*kind);
    if (spec.cols < 0) {
        if (err_out) *err_out = QStringLiteral("required ring array %1 has variable catalog width")
                                    .arg(fieldStem(spec));
        return false;
    }
    const StaticNpyArray* a = run.producerArray(*kind);
    if (!a) {
        if (err_out) *err_out = QStringLiteral("row_design ring restore requires %1.npy")
                                    .arg(fieldStem(spec));
        return false;
    }
    const std::size_t expectedRows =
        a->frameVarying ? atomCount * run.frameMap.frameCount() : atomCount;
    const std::size_t expectedCols = static_cast<std::size_t>(spec.cols);
    if (a->rows != expectedRows || a->cols != expectedCols) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has shape (%2,%3); expected (%4,%5)%6")
                           .arg(a->path)
                           .arg(a->rows)
                           .arg(a->cols)
                           .arg(expectedRows)
                           .arg(expectedCols)
                           .arg(a->frameVarying ? QStringLiteral(" flattened frame-major")
                                                : QString());
        }
        return false;
    }
    return true;
}

std::size_t rowsPerFrameForAxis(const RunData& run,
                                io::NativeAxis axis,
                                std::size_t atomCount) {
    if (!run.protein) return 0;
    switch (axis) {
    case io::NativeAxis::Atom:
        return atomCount;
    case io::NativeAxis::AromaticRing:
        return run.protein->topology().aromaticRingCount();
    case io::NativeAxis::SaturatedRing:
        return run.protein->topology().saturatedRingCount();
    default:
        return 0;
    }
}

QString trajectoryFrameNpyPath(const RunData& run,
                               const QString& npysDir,
                               std::size_t row,
                               io::FieldKind kind) {
    const std::size_t original = run.frameMap.originalIndex(row);
    const QString frameDir = QDir(npysDir).filePath(
        QStringLiteral("frame_%1").arg(original, 6, 10, QLatin1Char('0')));
    return QDir(frameDir).filePath(fieldStem(io::FieldSpecFor(kind)) + QStringLiteral(".npy"));
}

bool loadTrajectorySourceArray(RunData& run,
                               const TrajectoryNpyRequirement& req,
                               std::size_t atomCount,
                               QString* err_out) {
    const std::optional<io::FieldKind> kind = producerKind(req.id, err_out);
    if (!kind) return false;
    if (run.producerArray(*kind)) return true;
    const io::FieldSpec& spec = io::FieldSpecFor(*kind);
    if (!run.manifest.trajectory) {
        if (err_out) *err_out = QStringLiteral("trajectory source arrays require trajectory manifest");
        return false;
    }
    const QString npysDir = QDir(run.manifest.trajectory->extraction_dir_abspath)
                                .filePath(QStringLiteral("npys"));
    if (!QFileInfo(npysDir).isDir()) {
        if (err_out) *err_out = QStringLiteral("trajectory per-frame NPY dir is missing: %1").arg(npysDir);
        return false;
    }

    const std::size_t rowsPerFrame = rowsPerFrameForAxis(run, spec.axis, atomCount);
    if (rowsPerFrame == 0) return true;
    const QString firstPath = trajectoryFrameNpyPath(run, npysDir, 0, *kind);
    if (!QFileInfo::exists(firstPath)) {
        if (req.required && err_out) {
            *err_out = QStringLiteral("required trajectory producer NPY is missing: %1")
                           .arg(firstPath);
        }
        return !req.required;
    }

    auto first = io::QtNpyReader::ReadArrayWidened(firstPath);
    if (!first.ok) {
        if (err_out) *err_out = first.error;
        return false;
    }
    if (first.rows != rowsPerFrame) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; catalog axis %3 expects %4")
                           .arg(firstPath)
                           .arg(first.rows)
                           .arg(static_cast<int>(spec.axis))
                           .arg(rowsPerFrame);
        }
        return false;
    }
    if (spec.cols != -1 && first.cols != static_cast<std::size_t>(spec.cols)) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 columns; catalog %3 expects %4")
                           .arg(firstPath)
                           .arg(first.cols)
                           .arg(fieldStem(spec))
                           .arg(spec.cols);
        }
        return false;
    }
    const std::size_t cols = first.cols;

    StaticNpyArray resident;
    resident.stem = fieldStem(spec);
    resident.path = npysDir;
    resident.rows = rowsPerFrame * run.frameMap.frameCount();
    resident.cols = cols;
    resident.frameVarying = true;
    resident.atomsPerFrame = rowsPerFrame;
    resident.frameCount = run.frameMap.frameCount();
    resident.dtype_descr = QStringLiteral("flattened_frame_series");
    resident.values.resize(resident.rows * resident.cols);

    const std::size_t valuesPerFrame = rowsPerFrame * cols;
    std::copy(first.data.begin(), first.data.end(), resident.values.begin());
    for (std::size_t frame = 0; frame < run.frameMap.frameCount(); ++frame) {
        if (frame == 0) continue;
        const QString path = trajectoryFrameNpyPath(run, npysDir, frame, *kind);
        if (!QFileInfo::exists(path)) {
            if (err_out) *err_out = QStringLiteral("trajectory source array is missing: %1").arg(path);
            return false;
        }
        auto arr = io::QtNpyReader::ReadArrayWidened(path);
        if (!arr.ok) {
            if (err_out) *err_out = arr.error;
            return false;
        }
        if (arr.rows != rowsPerFrame || arr.cols != cols) {
            if (err_out) {
                *err_out = QStringLiteral("%1 has shape (%2,%3); expected (%4,%5)")
                               .arg(path)
                               .arg(arr.rows)
                               .arg(arr.cols)
                               .arg(rowsPerFrame)
                               .arg(cols);
            }
            return false;
        }
        const std::size_t offset = frame * valuesPerFrame;
        std::copy(arr.data.begin(), arr.data.end(), resident.values.begin() + offset);
    }

    run.producerArrays[static_cast<int>(*kind)] = std::move(resident);
    return true;
}

bool reconstructTrajectoryPerTypeT1FromRingContributions(const QString& frameDir,
                                                         io::FieldKind kind,
                                                         std::size_t atomCount,
                                                         double* frameOut,
                                                         QString* err_out);

bool loadTrajectoryRingArray(RunData& run,
                             const RingArrayRequirement& req,
                             std::size_t atomCount,
                             QString* err_out) {
    const std::optional<io::FieldKind> kind = producerKind(req.id, err_out);
    if (!kind) return false;
    const io::FieldSpec& spec = io::FieldSpecFor(*kind);
    if (spec.cols < 0) {
        if (err_out) *err_out = QStringLiteral("required ring array %1 has variable catalog width")
                                    .arg(fieldStem(spec));
        return false;
    }
    const std::size_t cols = static_cast<std::size_t>(spec.cols);
    if (!run.manifest.trajectory) {
        if (err_out) *err_out = QStringLiteral("trajectory ring arrays require trajectory manifest");
        return false;
    }
    const QString npysDir = QDir(run.manifest.trajectory->extraction_dir_abspath)
                                .filePath(QStringLiteral("npys"));
    if (!QFileInfo(npysDir).isDir()) {
        if (err_out) *err_out = QStringLiteral("trajectory per-frame NPY dir is missing: %1").arg(npysDir);
        return false;
    }

    StaticNpyArray resident;
    resident.stem = fieldStem(spec);
    resident.path = npysDir;
    resident.rows = atomCount * run.frameMap.frameCount();
    resident.cols = cols;
    resident.frameVarying = true;
    resident.atomsPerFrame = atomCount;
    resident.frameCount = run.frameMap.frameCount();
    resident.dtype_descr = QStringLiteral("flattened_frame_series");
    resident.values.resize(resident.rows * resident.cols);

    const std::size_t valuesPerFrame = atomCount * cols;
    for (std::size_t frame = 0; frame < run.frameMap.frameCount(); ++frame) {
        const std::size_t original = run.frameMap.originalIndex(frame);
        const QString frameDir = QDir(npysDir).filePath(
            QStringLiteral("frame_%1").arg(original, 6, 10, QLatin1Char('0')));
        const QString path = QDir(frameDir).filePath(fieldStem(spec) + QStringLiteral(".npy"));
        if (!QFileInfo::exists(path)) {
            if ((*kind == io::FieldKind::BSPerTypeT1 || *kind == io::FieldKind::HMPerTypeT1)
                && reconstructTrajectoryPerTypeT1FromRingContributions(
                    frameDir, *kind, atomCount, resident.values.data() + frame * valuesPerFrame, err_out)) {
                continue;
            }
            if (err_out) *err_out = QStringLiteral("required ring array is missing: %1").arg(path);
            return false;
        }
        auto arr = io::QtNpyReader::ReadArrayWidened(path);
        if (!arr.ok) {
            if (err_out) *err_out = arr.error;
            return false;
        }
        if (arr.rows != atomCount || arr.cols != cols) {
            if (err_out) {
                *err_out = QStringLiteral("%1 has shape (%2,%3); expected (%4,%5)")
                               .arg(path)
                               .arg(arr.rows)
                               .arg(arr.cols)
                               .arg(atomCount)
                               .arg(cols);
            }
            return false;
        }
        const std::size_t offset = frame * valuesPerFrame;
        std::copy(arr.data.begin(), arr.data.end(), resident.values.begin() + offset);
    }

    run.producerArrays[static_cast<int>(*kind)] = std::move(resident);
    return true;
}

bool reconstructTrajectoryPerTypeT1FromRingContributions(const QString& frameDir,
                                                         io::FieldKind kind,
                                                         std::size_t atomCount,
                                                         double* frameOut,
                                                         QString* err_out) {
    const bool bs = kind == io::FieldKind::BSPerTypeT1;
    const bool hm = kind == io::FieldKind::HMPerTypeT1;
    if (!bs && !hm) return false;

    const io::FieldSpec& outSpec = io::FieldSpecFor(kind);
    const io::FieldSpec& sourceSpec = io::FieldSpecFor(io::FieldKind::RingContributions);
    const std::size_t outCols = static_cast<std::size_t>(outSpec.cols);
    const QString path = QDir(frameDir).filePath(fieldStem(sourceSpec) + QStringLiteral(".npy"));
    if (!QFileInfo::exists(path)) {
        if (err_out) *err_out = QStringLiteral("cannot reconstruct %1: missing %2")
                                    .arg(fieldStem(outSpec), path);
        return false;
    }
    auto arr = io::QtNpyReader::ReadArrayWidened(path);
    if (!arr.ok) {
        if (err_out) *err_out = arr.error;
        return false;
    }
    if (arr.cols < 36) {
        if (err_out) *err_out = QStringLiteral("%1 has %2 columns; expected at least 36")
                                    .arg(path)
                                    .arg(arr.cols);
        return false;
    }

    std::fill(frameOut, frameOut + atomCount * outCols, 0.0);
    const std::size_t sourceBase = bs ? 9 : 27;
    for (std::size_t r = 0; r < arr.rows; ++r) {
        const std::size_t base = r * arr.cols;
        const int atom = static_cast<int>(std::llround(arr.data[base + 0]));
        const int ringType = static_cast<int>(std::llround(arr.data[base + 2]));
        if (atom < 0 || static_cast<std::size_t>(atom) >= atomCount
            || ringType < 0 || ringType >= model::kAromaticRingTypeCount) {
            continue;
        }
        for (int c = 0; c < 3; ++c) {
            frameOut[static_cast<std::size_t>(atom) * outCols
                     + static_cast<std::size_t>(ringType * 3 + c)] +=
                arr.data[base + sourceBase + 1 + static_cast<std::size_t>(c)];
        }
    }
    return true;
}

QString num(double v) {
    if (std::isnan(v)) return QStringLiteral("NaN");
    if (!std::isfinite(v)) return v > 0 ? QStringLiteral("Inf") : QStringLiteral("-Inf");
    return QString::number(v, 'g', 12);
}
QString intv(qint64 v) { return QString::number(v); }
QString boolv(bool v) { return v ? QStringLiteral("1") : QStringLiteral("0"); }

std::vector<QString> blankValues() {
    std::vector<QString> values;
    values.reserve(RowDesignSchema().size());
    for (const RowColumnSpec& c : RowDesignSchema()) {
        switch (c.type) {
        case RowColType::String: values.push_back(QString()); break;
        case RowColType::Bool: values.push_back(QStringLiteral("0")); break;
        case RowColType::Int: values.push_back(QStringLiteral("-1")); break;
        case RowColType::Double: values.push_back(QStringLiteral("NaN")); break;
        }
    }
    return values;
}

struct RowWriteContext {
    std::vector<QString>& values;
    RowDesignStats* stats = nullptr;
    std::vector<unsigned char> counted;

    RowWriteContext(std::vector<QString>& rowValues, RowDesignStats* rowStats)
        : values(rowValues), stats(rowStats), counted(RowDesignSchema().size(), 0) {
        if (stats && stats->populatedCounts.size() != RowDesignSchema().size())
            stats->populatedCounts.assign(RowDesignSchema().size(), 0);
    }

    void count(std::size_t idx, bool populated) {
        if (!stats || !populated || counted[idx]) return;
        ++stats->populatedCounts[idx];
        counted[idx] = 1;
    }

    void set(const QString& name, const QString& v) {
        const int idx = RowDesignColumnIndex(name);
        if (idx < 0) return;
        const std::size_t pos = static_cast<std::size_t>(idx);
        values[pos] = v;
        const RowColType type = RowDesignSchema()[pos].type;
        count(pos, type == RowColType::String ? !v.isEmpty() : !v.isEmpty());
    }
    void set(const char* name, const QString& v) { set(QString::fromLatin1(name), v); }

    void set(const QString& name, double v) {
        const int idx = RowDesignColumnIndex(name);
        if (idx < 0) return;
        const std::size_t pos = static_cast<std::size_t>(idx);
        values[pos] = num(v);
        count(pos, std::isfinite(v));
    }
    void set(const char* name, double v) { set(QString::fromLatin1(name), v); }

    void setInt(const char* name, qint64 v) {
        const int idx = RowDesignColumnIndex(QString::fromLatin1(name));
        if (idx < 0) return;
        const std::size_t pos = static_cast<std::size_t>(idx);
        values[pos] = intv(v);
        count(pos, true);
    }

    void setBool(const char* name, bool v) {
        const int idx = RowDesignColumnIndex(QString::fromLatin1(name));
        if (idx < 0) return;
        const std::size_t pos = static_cast<std::size_t>(idx);
        values[pos] = boolv(v);
        count(pos, true);
    }
};

struct PredecessorChi1Region {
    int region = -1;
    bool present = false;
    QString missingReason;
};

PredecessorChi1Region predecessorChi1Region(const Body& body,
                                            const model::QtProtein& protein,
                                            const model::QtResidue* residue,
                                            std::size_t frame) {
    PredecessorChi1Region out;
    if (!residue) {
        out.missingReason = QStringLiteral("no_current_residue");
        return out;
    }
    if (residue->prevResidueIndex < 0
        || static_cast<std::size_t>(residue->prevResidueIndex) >= protein.residueCount()) {
        out.missingReason = QStringLiteral("no_predecessor");
        return out;
    }
    const model::QtResidue& prev =
        protein.residue(static_cast<std::size_t>(residue->prevResidueIndex));
    if (prev.atomIndices.empty()) {
        out.missingReason = QStringLiteral("predecessor_empty");
        return out;
    }
    for (int32_t atomIndex : prev.atomIndices) {
        if (atomIndex < 0 || static_cast<std::size_t>(atomIndex) >= protein.atomCount()) continue;
        const DihedralState chi =
            body.idx.dihedrals.state(DihedralKind::Chi1,
                                     static_cast<std::size_t>(atomIndex),
                                     frame);
        if (!chi.present) continue;
        out.region = chi.fixed_bin;
        out.present = true;
        return out;
    }
    out.missingReason = QStringLiteral("predecessor_chi1_absent");
    return out;
}

double t2Magnitude(const std::array<double, 5>& t2) {
    double s = 0.0;
    for (double v : t2) {
        if (!std::isfinite(v)) return kNaN;
        s += v * v;
    }
    return std::sqrt(s);
}

bool finiteTensor(const model::SphericalTensor& t) {
    if (!std::isfinite(t.T0)) return false;
    for (double v : t.T1)
        if (!std::isfinite(v)) return false;
    for (double v : t.T2)
        if (!std::isfinite(v)) return false;
    return true;
}

void addTensor(model::SphericalTensor& dst, const model::SphericalTensor& src) {
    dst.T0 += src.T0;
    for (std::size_t i = 0; i < dst.T1.size(); ++i) dst.T1[i] += src.T1[i];
    for (std::size_t i = 0; i < dst.T2.size(); ++i) dst.T2[i] += src.T2[i];
}

bool addPresentTensor(const Body& body, ArrayId id, std::size_t atom, std::size_t row,
                      model::SphericalTensor& dst) {
    if (!body.catalog.present(body, id, atom, row)) return false;
    const model::SphericalTensor t = body.catalog.valueTensor(body, id, atom, row);
    if (!finiteTensor(t)) return false;
    addTensor(dst, t);
    return true;
}

double requiredComponent(const Body& body, ArrayId id, std::size_t atom, std::size_t row, int comp) {
    if (!body.catalog.present(body, id, atom, row)) {
        throw std::runtime_error(QStringLiteral("row_design required ring array is absent: %1")
                                     .arg(body.catalog.spec(id).name)
                                     .toStdString());
    }
    const double v = body.catalog.value(body, id, atom, row, -1, comp);
    if (!std::isfinite(v)) {
        throw std::runtime_error(QStringLiteral("row_design required ring array has non-finite value: %1")
                                     .arg(body.catalog.spec(id).name)
                                     .toStdString());
    }
    return v;
}

model::SphericalTensor scaledPerTypeTensor(const Body& body,
                                           ArrayId t0Id,
                                           ArrayId t1Id,
                                           ArrayId t2Id,
                                           std::size_t atom,
                                           std::size_t row,
                                           std::array<double, model::kAromaticRingTypeCount>* scaledT0ByType = nullptr) {
    model::SphericalTensor out;
    if (scaledT0ByType) scaledT0ByType->fill(0.0);
    for (int t = 0; t < model::kAromaticRingTypeCount; ++t) {
        const double intensity = kRingLiteratureIntensity[static_cast<std::size_t>(t)];
        const double scaledT0 = intensity * requiredComponent(body, t0Id, atom, row, t);
        out.T0 += scaledT0;
        if (scaledT0ByType) (*scaledT0ByType)[static_cast<std::size_t>(t)] = scaledT0;
        for (int c = 0; c < 3; ++c)
            out.T1[static_cast<std::size_t>(c)] +=
                intensity * requiredComponent(body, t1Id, atom, row, t * 3 + c);
        for (int c = 0; c < 5; ++c)
            out.T2[static_cast<std::size_t>(c)] +=
                intensity * requiredComponent(body, t2Id, atom, row, t * 5 + c);
    }
    return out;
}

void setTensorColumns(RowWriteContext& row,
                      const char* prefix,
                      const model::SphericalTensor& t) {
    row.set(QStringLiteral("%1_T0").arg(QString::fromLatin1(prefix)), t.T0);
    for (int c = 0; c < 3; ++c) {
        row.set(QStringLiteral("%1_T1_%2").arg(QString::fromLatin1(prefix)).arg(c),
                t.T1[static_cast<std::size_t>(c)]);
    }
    for (int c = 0; c < 5; ++c) {
        row.set(QStringLiteral("%1_T2_%2").arg(QString::fromLatin1(prefix)).arg(c),
                t.T2[static_cast<std::size_t>(c)]);
    }
    row.set(QStringLiteral("%1_absT2").arg(QString::fromLatin1(prefix)), t2Magnitude(t.T2));
}

void appendTensor9(std::array<double, kRingTensorSidecarDims>& out,
                   std::size_t* idx,
                   const model::SphericalTensor& t) {
    out[(*idx)++] = t.T0;
    for (double v : t.T1) out[(*idx)++] = v;
    for (double v : t.T2) out[(*idx)++] = v;
}

void fillRingTensorSidecar(const Body& body,
                           std::size_t atom,
                           std::size_t row,
                           std::array<double, kRingTensorSidecarDims>* out) {
    std::size_t idx = 0;
    for (int c = 0; c < 8; ++c)
        (*out)[idx++] = requiredComponent(body, ArrayId::BSPerTypeT0, atom, row, c);
    for (int c = 0; c < 24; ++c)
        (*out)[idx++] = requiredComponent(body, ArrayId::BSPerTypeT1, atom, row, c);
    for (int c = 0; c < 40; ++c)
        (*out)[idx++] = requiredComponent(body, ArrayId::BSPerTypeT2, atom, row, c);
    for (int c = 0; c < 8; ++c)
        (*out)[idx++] = requiredComponent(body, ArrayId::HMPerTypeT0, atom, row, c);
    for (int c = 0; c < 24; ++c)
        (*out)[idx++] = requiredComponent(body, ArrayId::HMPerTypeT1, atom, row, c);
    for (int c = 0; c < 40; ++c)
        (*out)[idx++] = requiredComponent(body, ArrayId::HMPerTypeT2, atom, row, c);

    if (!body.catalog.present(body, ArrayId::KernelBs, atom, row)
        || !body.catalog.present(body, ArrayId::HmShielding, atom, row)) {
        throw std::runtime_error("row_design required summed ring tensors are absent");
    }
    appendTensor9(*out, &idx, body.catalog.valueTensor(body, ArrayId::KernelBs, atom, row));
    appendTensor9(*out, &idx, body.catalog.valueTensor(body, ArrayId::HmShielding, atom, row));
    if (idx != kRingTensorSidecarDims)
        throw std::runtime_error("row_design internal ring sidecar width mismatch");
}

bool mcProxyTensor(const Body& body, std::size_t atom, std::size_t row,
                   model::SphericalTensor* out) {
    model::SphericalTensor total;

    if (addPresentTensor(body, ArrayId::MopacMcShielding, atom, row, total)) {
        *out = total;
        return true;
    }

    const ArrayId bondOrderMc[] = {
        ArrayId::McPeptideCoBo,
        ArrayId::McPeptideCnBo,
        ArrayId::McBackboneOtherBo,
        ArrayId::McSidechainCoBo,
        ArrayId::McSidechainOtherBo,
        ArrayId::McDisulfideBo,
        ArrayId::McAromaticZeroedBo,
    };
    bool any = false;
    total = {};
    for (ArrayId id : bondOrderMc)
        any = addPresentTensor(body, id, atom, row, total) || any;
    if (any) {
        *out = total;
        return true;
    }

    const ArrayId nearestMc[] = {ArrayId::McNearestCoT2, ArrayId::McNearestCnT2};
    total = {};
    any = false;
    for (ArrayId id : nearestMc)
        any = addPresentTensor(body, id, atom, row, total) || any;
    if (any) {
        *out = total;
        return true;
    }

    if (addPresentTensor(body, ArrayId::KernelMc, atom, row, total)) {
        *out = total;
        return true;
    }

    const ArrayId fixedMc[] = {
        ArrayId::McPeptideCoFixed,
        ArrayId::McPeptideCnFixed,
        ArrayId::McBackboneOtherFixed,
        ArrayId::McSidechainCoFixed,
        ArrayId::McSidechainOtherFixed,
        ArrayId::McDisulfideFixed,
        ArrayId::McAromaticZeroedFixed,
    };
    total = {};
    any = false;
    for (ArrayId id : fixedMc)
        any = addPresentTensor(body, id, atom, row, total) || any;
    if (any) {
        *out = total;
        return true;
    }

    return false;
}

bool validResidue(const model::QtProtein& p, int32_t res) {
    return res >= 0 && static_cast<std::size_t>(res) < p.residueCount();
}

const StaticNpyArray* staticArray(const RunData& run, ArrayId id) {
    const std::optional<io::FieldKind> kind = ProducerFieldFor(id);
    return kind ? run.producerArray(*kind) : nullptr;
}

double staticByResidue(const RunData& run, ArrayId id, int32_t residueIndex) {
    const StaticNpyArray* a = staticArray(run, id);
    if (!a || residueIndex < 0 || static_cast<std::size_t>(residueIndex) >= a->rows) return kNaN;
    return a->value(static_cast<std::size_t>(residueIndex));
}

double staticAxisValue(const RunData& run,
                       ArrayId id,
                       std::size_t nativeIndex,
                       std::size_t row) {
    const StaticNpyArray* a = staticArray(run, id);
    if (!a || a->cols == 0) return kNaN;
    if (a->frameVarying) {
        if (nativeIndex >= a->atomsPerFrame || row >= a->frameCount) return kNaN;
    }
    const std::size_t sourceRow = a->rowFor(nativeIndex, row);
    if (sourceRow >= a->rows) return kNaN;
    const double v = a->value(sourceRow);
    return std::isfinite(v) ? v : kNaN;
}

int ringNativeAxisForAtom(const model::QtProtein& p,
                          std::size_t atom,
                          model::RingKind kind) {
    const model::QtTopology& topo = p.topology();
    for (int32_t membershipIndex : topo.ringMembershipsForAtom(atom)) {
        if (membershipIndex < 0
            || static_cast<std::size_t>(membershipIndex) >= topo.ringMembershipCount()) {
            continue;
        }
        const model::QtRingMembership& membership =
            topo.ringMembershipAt(static_cast<std::size_t>(membershipIndex));
        if (membership.ringId < 0
            || static_cast<std::size_t>(membership.ringId) >= topo.ringCount()) {
            continue;
        }
        const model::QtRing& ring = topo.ringAt(static_cast<std::size_t>(membership.ringId));
        if (ring.ringKind == kind && ring.nativeAxisIndex >= 0) return ring.nativeAxisIndex;
    }
    return -1;
}

int ss3FromDssp(model::DsspCode c) {
    switch (c) {
    case model::DsspCode::AlphaHelix:
    case model::DsspCode::Helix310:
    case model::DsspCode::PiHelix:
        return 0;
    case model::DsspCode::ExtendedStrand:
    case model::DsspCode::BetaBridge:
        return 1;
    default:
        return 2;
    }
}

struct AngleBlock {
    double phi = kNaN;
    double psi = kNaN;
    double omega = kNaN;
    std::array<double, 4> chi = {kNaN, kNaN, kNaN, kNaN};
    std::array<bool, 4> chiExists = {};
    bool phiPsiPresent = false;
};

AngleBlock readAngles(const Body& body, const model::QtAtom& atom, std::size_t row) {
    AngleBlock out;
    const int32_t resIdx = atom.residueIndex;
    if (resIdx < 0) return out;
    if (const io::QtTrajectoryH5* h5 = body.run.h5()) {
        const model::QtDihedralTimeSeries* dih = h5->dihedrals();
        if (!dih) throw std::runtime_error("row_design requires trajectory dihedrals");
        if (body.run.conformation && dih->n_frames != body.run.conformation->frameCount())
            throw std::runtime_error("row_design dihedral frame count does not match trajectory");
        const std::size_t r = static_cast<std::size_t>(resIdx);
        // The current 1P9J dihedral_time_series H5 channel is stored with
        // phi/psi sign opposite to standard IUPAC despite the producer attr.
        // Normalize here so the row-design contract is one convention.
        out.phi = -dih->phiAt(r, row);
        out.psi = -dih->psiAt(r, row);
        out.omega = dih->omegaAt(r, row);
        out.phiPsiPresent = std::isfinite(out.phi) && std::isfinite(out.psi);
        for (int k = 0; k < 4; ++k) {
            const std::size_t mask = r * 4 + static_cast<std::size_t>(k);
            out.chiExists[k] = mask < dih->chi_exists.size() && dih->chi_exists[mask] != 0;
            out.chi[k] = out.chiExists[k] ? dih->chiAt(r, row, k) : kNaN;
        }
        return out;
    }

    const StaticNpyArray* bb = staticArray(body.run, ArrayId::DsspBackbone);
    if (!bb) throw std::runtime_error("row_design static run requires dssp_backbone.npy");
    const std::size_t atomIndex = static_cast<std::size_t>(atom.atomIndex);
    if (atomIndex < bb->rows && bb->cols >= 2) {
        // The full720 static dssp_backbone arrays are already in standard
        // IUPAC backbone convention; do not apply the trajectory H5 flip.
        out.phi = bb->value(atomIndex, 0);
        out.psi = bb->value(atomIndex, 1);
        out.phiPsiPresent = std::isfinite(out.phi) && std::isfinite(out.psi);
    }
    out.omega = staticByResidue(body.run, ArrayId::OmegaActual, resIdx);
    const StaticNpyArray* chi = staticArray(body.run, ArrayId::DsspChi);
    if (chi && atomIndex < chi->rows && chi->cols >= 12) {
        for (int k = 0; k < 4; ++k) {
            const std::size_t base = static_cast<std::size_t>(k) * 3;
            const double exists = chi->value(atomIndex, base + 2);
            out.chiExists[k] = exists != 0.0;
            if (out.chiExists[k]) {
                const double c = chi->value(atomIndex, base);
                const double s = chi->value(atomIndex, base + 1);
                out.chi[k] = std::atan2(s, c);
            }
        }
    }
    return out;
}

struct DsspBlock {
    int ss8 = static_cast<int>(model::DsspCode::Unknown);
    int ss3 = 2;
    double energy = kNaN;
    bool donor = false;
    bool acceptor = false;
    bool present = false;
};

DsspBlock readDssp(const Body& body, const model::QtAtom& atom, std::size_t row) {
    DsspBlock out;
    if (atom.residueIndex < 0) return out;
    const std::size_t res = static_cast<std::size_t>(atom.residueIndex);
    if (const io::QtTrajectoryH5* h5 = body.run.h5()) {
        const model::QtDssp8TimeSeries* dssp = h5->dssp8();
        if (!dssp || res >= dssp->n_residues || row >= dssp->n_frames || !dssp->sourceAttachedAt(row))
            return out;
        const model::DsspCode code = dssp->codeAt(res, row);
        out.ss8 = static_cast<int>(code);
        out.ss3 = ss3FromDssp(code);
        const std::size_t base = (res * dssp->n_frames + row) * 2;
        std::vector<double> energies;
        if (base + 1 < dssp->hbond_acceptor_energy.size()) {
            out.acceptor = dssp->hbond_acceptor_partner[base] >= 0
                           || dssp->hbond_acceptor_partner[base + 1] >= 0;
            energies.push_back(dssp->hbond_acceptor_energy[base]);
            energies.push_back(dssp->hbond_acceptor_energy[base + 1]);
        }
        if (base + 1 < dssp->hbond_donor_energy.size()) {
            out.donor = dssp->hbond_donor_partner[base] >= 0
                        || dssp->hbond_donor_partner[base + 1] >= 0;
            energies.push_back(dssp->hbond_donor_energy[base]);
            energies.push_back(dssp->hbond_donor_energy[base + 1]);
        }
        for (double e : energies)
            if (std::isfinite(e) && e != 0.0)
                out.energy = std::isnan(out.energy) ? e : std::min(out.energy, e);
        out.present = true;
        return out;
    }

    const std::size_t atomIndex = static_cast<std::size_t>(atom.atomIndex);
    if (const StaticNpyArray* ss = staticArray(body.run, ArrayId::DsspSs8)) {
        if (atomIndex < ss->rows && ss->cols >= 8) {
            std::size_t best = 0;
            for (std::size_t i = 1; i < 8; ++i)
                if (ss->value(atomIndex, i) > ss->value(atomIndex, best)) best = i;
            out.ss8 = static_cast<int>(best);
            out.ss3 = ss3FromDssp(static_cast<model::DsspCode>(best));
            out.present = true;
        }
    }
    if (const StaticNpyArray* e = staticArray(body.run, ArrayId::DsspHBondEnergy)) {
        if (atomIndex < e->rows && e->cols >= 4) {
            for (std::size_t i = 0; i < 4; ++i) {
                const double v = e->value(atomIndex, i);
                if (std::isfinite(v) && v != 0.0)
                    out.energy = std::isnan(out.energy) ? v : std::min(out.energy, v);
            }
            out.acceptor = (e->value(atomIndex, 0) != 0.0 || e->value(atomIndex, 1) != 0.0);
            out.donor = (e->value(atomIndex, 2) != 0.0 || e->value(atomIndex, 3) != 0.0);
        }
    }
    return out;
}

std::vector<std::size_t> allAtoms(const Body& body) {
    std::vector<std::size_t> atoms;
    if (!body.run.protein) return atoms;
    atoms.reserve(body.run.protein->atomCount());
    for (std::size_t i = 0; i < body.run.protein->atomCount(); ++i) atoms.push_back(i);
    return atoms;
}

FrameResult labFrame(const Body&, std::size_t, std::size_t) { return {}; }

std::size_t countNear(const Body& body, CloudKind kind, std::size_t atom, std::size_t row, double cutoff) {
    std::size_t n = 0;
    for (const SourceRef& ref : verbs::near(body, kind, atom, row, cutoff)) {
        if (kind == CloudKind::Atoms && ref.entity_index == static_cast<int32_t>(atom)) continue;
        ++n;
    }
    return n;
}

double nearestDist(const Body& body, CloudKind kind, std::size_t atom, std::size_t row,
                   int32_t* entityOut = nullptr) {
    const Vec3 p = verbs::pos(body, atom, row);
    double best = kNaN;
    int32_t bestEntity = -1;
    for (const SourceRef& ref : verbs::near(body, kind, atom, row, 99.0)) {
        if (kind == CloudKind::Atoms && ref.entity_index == static_cast<int32_t>(atom)) continue;
        const Vec3 q = body.idx.spatial.tree(kind, row).pointAt(static_cast<std::size_t>(ref.cloud_index));
        const double d = (q - p).norm();
        if (std::isnan(best) || d < best) {
            best = d;
            bestEntity = ref.entity_index;
        }
    }
    if (entityOut) *entityOut = bestEntity;
    return best;
}

RowDesignRow buildRow(const Body& body, std::size_t atom, std::size_t row, std::size_t original,
                      RowDesignStats* stats) {
    RowDesignRow out;
    out.values = blankValues();
    RowWriteContext rowCtx(out.values, stats);

    const model::QtProtein& p = *body.run.protein;
    const model::QtAtom& a = p.atom(atom);
    const model::QtResidue* r = validResidue(p, a.residueIndex)
                                    ? &p.residue(static_cast<std::size_t>(a.residueIndex))
                                    : nullptr;

    const QString dataset = body.run.manifest.dataset_id.isEmpty()
                                ? (body.run.poseKind() == PoseKind::Static ? QStringLiteral("720_static")
                                                                           : QStringLiteral("1p9j"))
                                : body.run.manifest.dataset_id;
    rowCtx.set("dataset_id", dataset);
    rowCtx.set("protein_id", body.run.protein->proteinId());
    rowCtx.set("case_id", QStringLiteral("row_design"));
    rowCtx.set("pose_kind", body.run.poseKind() == PoseKind::Static ? QStringLiteral("static")
                                                                         : QStringLiteral("trajectory"));
    rowCtx.setInt("frame_slot", static_cast<qint64>(row));
    rowCtx.set("split_group_id", body.run.protein->proteinId());
    rowCtx.setInt("atom_index", a.atomIndex);
    rowCtx.setInt("h5_row", static_cast<qint64>(row));
    rowCtx.setInt("original_index", static_cast<qint64>(original));
    rowCtx.set("time_ps", body.run.timePs(row));
    rowCtx.setInt("element", a.AtomicNumber());
    rowCtx.setInt("residue_index", a.residueIndex);
    rowCtx.setInt("residue_type", r ? static_cast<int>(r->aminoAcid) : -1);
    rowCtx.setInt("residue_number", r ? r->address.residueNumber : 0);
    rowCtx.set("atom_name", p.atomLabel(atom, model::NamingConvention::Bmrb));
    rowCtx.setInt("iupac_role", static_cast<int>(a.pseudoatomKind));
    rowCtx.setInt("backbone_role", static_cast<int>(a.backboneRole));
    rowCtx.setInt("locant", static_cast<int>(a.locant));
    rowCtx.setInt("branch_outer", a.branch.outer);
    rowCtx.setInt("branch_inner", a.branch.inner);
    rowCtx.setInt("di_index", static_cast<int>(a.diIndex));
    rowCtx.setInt("prochiral", static_cast<int>(a.prochiral));
    rowCtx.setInt("ring_position_primary", static_cast<int>(a.ringPositionPrimary));
    rowCtx.setInt("ring_position_secondary", static_cast<int>(a.ringPositionSecondary));
    rowCtx.setInt("equivalence_class", a.equivalenceClass);
    rowCtx.setInt("formal_charge", a.formalCharge);
    rowCtx.setBool("is_exchangeable", a.isExchangeable);
    rowCtx.setInt("polar_h_kind", static_cast<int>(a.polarH));
    rowCtx.setBool("aromatic", a.aromatic);
    const ResidueClass rc = r ? ClassifyResidue(r->aminoAcid) : ResidueClass::Other;
    rowCtx.setInt("residue_class", static_cast<int>(rc));
    rowCtx.setInt("prev_class", r ? static_cast<int>(ClassifyResidue(r->prevResidueType)) : -1);
    rowCtx.setInt("next_class", r ? static_cast<int>(ClassifyResidue(r->nextResidueType)) : -1);
    rowCtx.setInt("prev_restype", r ? static_cast<int>(r->prevResidueType) : -1);
    rowCtx.setInt("next_restype", r ? static_cast<int>(r->nextResidueType) : -1);
    const PredecessorChi1Region prevChi1 = predecessorChi1Region(body, p, r, row);
    rowCtx.setInt("prev_chi1_region", prevChi1.region);
    rowCtx.setBool("prev_chi1_present", prevChi1.present);
    rowCtx.set("prev_chi1_missing_reason", prevChi1.missingReason);
    rowCtx.setBool("pre_proline", r && r->nextResidueType == model::AminoAcid::PRO);
    rowCtx.setBool("is_pro", r && r->isProline);
    rowCtx.setBool("is_gly", r && r->aminoAcid == model::AminoAcid::GLY);
    rowCtx.setInt("terminal_state", r ? static_cast<int>(r->terminalState) : -1);

    const AngleBlock angles = readAngles(body, a, row);
    rowCtx.set("phi", angles.phi);
    rowCtx.set("psi", angles.psi);
    rowCtx.set("sin_phi", std::sin(angles.phi));
    rowCtx.set("cos_phi", std::cos(angles.phi));
    rowCtx.set("sin_psi", std::sin(angles.psi));
    rowCtx.set("cos_psi", std::cos(angles.psi));
    const RowRamaRegion rama = ClassifyRowRama(angles.phi, angles.psi);
    rowCtx.setInt("rama_region", static_cast<int>(rama));
    for (int k = 0; k < 4; ++k) {
        const QString prefix = QStringLiteral("chi%1").arg(k + 1);
        rowCtx.set(QStringLiteral("sin_%1").arg(prefix).toLatin1().constData(),
            angles.chiExists[k] ? std::sin(angles.chi[k]) : kNaN);
        rowCtx.set(QStringLiteral("cos_%1").arg(prefix).toLatin1().constData(),
            angles.chiExists[k] ? std::cos(angles.chi[k]) : kNaN);
        rowCtx.setBool(QStringLiteral("%1_exists").arg(prefix).toLatin1().constData(), angles.chiExists[k]);
    }
    rowCtx.set("omega", angles.omega);
    rowCtx.set("sin_omega", std::sin(angles.omega));
    rowCtx.set("cos_omega", std::cos(angles.omega));
    if (angles.phiPsiPresent) ++stats->phiPsiPresent;
    const bool phiPsiGateEligible =
        a.IsBackbone() && r && !r->IsTerminalN() && !r->IsTerminalC();
    if (phiPsiGateEligible) {
        ++stats->phiPsiEligible;
        if (angles.phiPsiPresent) ++stats->phiPsiFiniteEligible;
    }

    const int saturatedNative =
        ringNativeAxisForAtom(p, atom, model::RingKind::Saturated);
    if (saturatedNative >= 0) {
        const double q = staticAxisValue(body.run, ArrayId::PuckerQ,
                                         static_cast<std::size_t>(saturatedNative), row);
        const double theta = staticAxisValue(body.run, ArrayId::PuckerTheta,
                                             static_cast<std::size_t>(saturatedNative), row);
        if (std::isfinite(q)) rowCtx.set("pucker_Q", q);
        if (std::isfinite(theta)) rowCtx.set("pucker_theta", theta);
    }
    if (body.catalog.present(body, ArrayId::Pyramidalization, atom, row))
        rowCtx.set("pyramidalization",
            body.catalog.value(body, ArrayId::Pyramidalization, atom, row));
    const int aromaticNative =
        ringNativeAxisForAtom(p, atom, model::RingKind::Aromatic);
    if (aromaticNative >= 0) {
        const double chi2 = staticAxisValue(body.run, ArrayId::AromaticChi2,
                                            static_cast<std::size_t>(aromaticNative), row);
        if (std::isfinite(chi2)) rowCtx.set("aromatic_ringflip_state", chi2);
    }

    const DsspBlock dssp = readDssp(body, a, row);
    rowCtx.setInt("dssp_ss8", dssp.ss8);
    rowCtx.setInt("dssp_ss3", dssp.ss3);
    rowCtx.set("dssp_hbond_energy", dssp.energy);
    rowCtx.setBool("hbond_donor", dssp.donor);
    rowCtx.setBool("hbond_acceptor", dssp.acceptor);
    if (dssp.present) ++stats->dsspPresent;

    for (int cutoff : {4, 6, 8, 10}) {
        rowCtx.setInt(QStringLiteral("n_atoms_%1A").arg(cutoff).toLatin1().constData(),
               static_cast<qint64>(countNear(body, CloudKind::Atoms, atom, row, cutoff)));
        rowCtx.setInt(QStringLiteral("n_rings_%1A").arg(cutoff).toLatin1().constData(),
               static_cast<qint64>(countNear(body, CloudKind::RingCenters, atom, row, cutoff)));
        rowCtx.setInt(QStringLiteral("n_charges_%1A").arg(cutoff).toLatin1().constData(),
               static_cast<qint64>(countNear(body, CloudKind::ChargeSites, atom, row, cutoff)));
        rowCtx.setInt(QStringLiteral("n_bonds_%1A").arg(cutoff).toLatin1().constData(),
               static_cast<qint64>(countNear(body, CloudKind::AllBondMidpoints, atom, row, cutoff)));
    }

    int32_t nearestRing = -1, nearestCharge = -1, nearestBond = -1, nearestAtom = -1;
    rowCtx.set("nearest_ring_dist", nearestDist(body, CloudKind::RingCenters, atom, row, &nearestRing));
    rowCtx.setInt("nearest_ring_identity_ord", nearestRing);
    rowCtx.set("nearest_charge_dist", nearestDist(body, CloudKind::ChargeSites, atom, row, &nearestCharge));
    rowCtx.setInt("nearest_charge_identity_ord", nearestCharge);
    rowCtx.set("nearest_charge_sign",
        nearestCharge >= 0 && body.catalog.present(body, ArrayId::Ff14sbCharge, nearestCharge, row)
            ? body.catalog.value(body, ArrayId::Ff14sbCharge, nearestCharge, row)
            : kNaN);
    rowCtx.set("nearest_bond_dist", nearestDist(body, CloudKind::AllBondMidpoints, atom, row, &nearestBond));
    rowCtx.setInt("nearest_bond_identity_ord", nearestBond);
    rowCtx.set("nearest_atom_dist", nearestDist(body, CloudKind::Atoms, atom, row, &nearestAtom));
    rowCtx.setInt("nearest_atom_identity_ord", nearestAtom);
    if (nearestRing >= 0) {
        const model::RingGeometry& g = body.idx.ringGeometry.at(static_cast<std::size_t>(nearestRing), row);
        const Vec3 disp = verbs::pos(body, atom, row) - g.center;
        const double z = disp.dot(g.normal);
        const Vec3 dPlane = disp - z * g.normal;
        const double rho = dPlane.norm();
        rowCtx.set("ring_cyl_z", z);
        rowCtx.set("ring_cyl_rho", rho);
        rowCtx.set("ring_angle_to_normal", disp.norm() > 1e-12 ? std::acos(std::clamp(std::abs(z) / disp.norm(), 0.0, 1.0)) : kNaN);
        const model::QtRing& ring = p.ring(static_cast<std::size_t>(nearestRing));
        if (!ring.atomIndices.empty() && body.run.conformation) {
            const Vec3 vertex0 =
                body.run.conformation->atomPosition(row, static_cast<std::size_t>(ring.atomIndices.front()));
            const Vec3 ref = vertex0 - g.center;
            const Vec3 refPlane = ref - ref.dot(g.normal) * g.normal;
            const double refNorm = refPlane.norm();
            if (rho > kRingAzimuthNormFloor && refNorm > kRingAzimuthNormFloor) {
                const Vec3 dHat = dPlane / rho;
                const Vec3 refHat = refPlane / refNorm;
                rowCtx.set("ring_cos_phi", dHat.dot(refHat));
                rowCtx.set("ring_sin_phi", dHat.cross(refHat).dot(g.normal));
            }
        }
    }

    if (body.catalog.present(body, ArrayId::KernelBs, atom, row)) {
        const model::SphericalTensor bs = body.catalog.valueTensor(body, ArrayId::KernelBs, atom, row);
        rowCtx.set("ring_bs_T0", bs.T0);
        rowCtx.set("ring_bs_absT2", t2Magnitude(bs.T2));
    }
    if (body.catalog.present(body, ArrayId::HmShielding, atom, row)) {
        const model::SphericalTensor hm = body.catalog.valueTensor(body, ArrayId::HmShielding, atom, row);
        rowCtx.set("ring_hm_T0", hm.T0);
        rowCtx.set("ring_hm_absT2", t2Magnitude(hm.T2));
    }
    std::array<double, model::kAromaticRingTypeCount> bsScaledT0ByType = {};
    const model::SphericalTensor ringJb =
        scaledPerTypeTensor(body, ArrayId::BSPerTypeT0, ArrayId::BSPerTypeT1, ArrayId::BSPerTypeT2,
                            atom, row, &bsScaledT0ByType);
    const model::SphericalTensor ringHmJb =
        scaledPerTypeTensor(body, ArrayId::HMPerTypeT0, ArrayId::HMPerTypeT1, ArrayId::HMPerTypeT2,
                            atom, row);
    setTensorColumns(rowCtx, "ring_jb", ringJb);
    setTensorColumns(rowCtx, "ring_hm_jb", ringHmJb);
    rowCtx.set("ring_jb_T0_phe", bsScaledT0ByType[0]);
    rowCtx.set("ring_jb_T0_tyr", bsScaledT0ByType[1]);
    rowCtx.set("ring_jb_T0_trp6", bsScaledT0ByType[2]);
    rowCtx.set("ring_jb_T0_trp5", bsScaledT0ByType[3]);
    rowCtx.set("ring_jb_T0_trp9", bsScaledT0ByType[4]);
    rowCtx.set("ring_jb_T0_his", bsScaledT0ByType[5]);
    rowCtx.set("ring_jb_T0_hid", bsScaledT0ByType[6]);
    rowCtx.set("ring_jb_T0_hie", bsScaledT0ByType[7]);
    fillRingTensorSidecar(body, atom, row, &out.ringTensors);
    out.ringTensorsPresent = true;
    model::SphericalTensor mc;
    if (mcProxyTensor(body, atom, row, &mc)) {
        rowCtx.set("mc_lit_T0", mc.T0);
        rowCtx.set("mc_lit_absT2", t2Magnitude(mc.T2));
    }
    if (body.catalog.present(body, ArrayId::MopacCoulombShielding, atom, row))
        rowCtx.set("mopac_bare_efg_kernel_absT2",
            t2Magnitude(body.catalog.valueT2(body, ArrayId::MopacCoulombShielding, atom, row)));
    if (body.catalog.present(body, ArrayId::ApbsEfg, atom, row))
        rowCtx.set("apbs_bare_efg_kernel_absT2",
            t2Magnitude(body.catalog.valueT2(body, ArrayId::ApbsEfg, atom, row)));
    if (body.catalog.present(body, ArrayId::Aimnet2Efg, atom, row))
        rowCtx.set("aimnet2_bare_efg_kernel_absT2",
            t2Magnitude(body.catalog.valueT2(body, ArrayId::Aimnet2Efg, atom, row)));
    if (body.catalog.present(body, ArrayId::MopacCoulombEfield, atom, row))
        rowCtx.set("mopac_efield_mag", body.catalog.valueVec3(body, ArrayId::MopacCoulombEfield, atom, row).norm());
    if (body.catalog.present(body, ArrayId::ApbsEfield, atom, row))
        rowCtx.set("apbs_efield_mag", body.catalog.valueVec3(body, ArrayId::ApbsEfield, atom, row).norm());
    if (body.catalog.present(body, ArrayId::Aimnet2Charge, atom, row))
        rowCtx.set("aimnet2_charge", body.catalog.value(body, ArrayId::Aimnet2Charge, atom, row));
    if (body.catalog.present(body, ArrayId::Ff14sbCharge, atom, row))
        rowCtx.set("ff14sb_charge", body.catalog.value(body, ArrayId::Ff14sbCharge, atom, row));
    if (body.catalog.present(body, ArrayId::MopacCharge, atom, row))
        rowCtx.set("mopac_charge", body.catalog.value(body, ArrayId::MopacCharge, atom, row));
    else if (body.catalog.present(body, ArrayId::MopacChargeWelfordMean, atom, row))
        rowCtx.set("mopac_charge", body.catalog.value(body, ArrayId::MopacChargeWelfordMean, atom, row));
    if (body.catalog.present(body, ArrayId::EeqChargeMean, atom, row))
        rowCtx.set("eeq_charge", body.catalog.value(body, ArrayId::EeqChargeMean, atom, row));
    if (body.catalog.present(body, ArrayId::EeqCoordinationNumber, atom, row))
        rowCtx.set("eeq_cn", body.catalog.value(body, ArrayId::EeqCoordinationNumber, atom, row));
    if (body.catalog.present(body, ArrayId::Sasa, atom, row))
        rowCtx.set("sasa", body.catalog.value(body, ArrayId::Sasa, atom, row));
    if (body.catalog.present(body, ArrayId::HbondCount, atom, row))
        rowCtx.set("larsen_hbond_count", body.catalog.value(body, ArrayId::HbondCount, atom, row));
    if (body.catalog.present(body, ArrayId::LarsenHBondShielding, atom, row))
        rowCtx.set("larsen_hbond_absT2",
            t2Magnitude(body.catalog.valueT2(body, ArrayId::LarsenHBondShielding, atom, row)));

    const DftTarget target = BuildTarget(body.run, atom, original, LocalFrame{});
    if (target.present) {
        rowCtx.set("target_T0", target.total_decomp.T0);
        out.targetT2 = target.total_decomp.T2;
        out.targetT2Present = true;
        ++stats->dftPresent;
    }
    rowCtx.set("tensor_frame", QStringLiteral("molecular_lab"));
    rowCtx.setBool("valid_for_T2_model", target.present);
    rowCtx.set("region_def_id", QStringLiteral("row_design_v1"));
    rowCtx.set("rama_region_hdr", QString::fromLatin1(NameForRowRama(rama)));
    rowCtx.set("rotamer_id", QString());
    std::size_t dims = 0;
    out.embedding = body.catalog.valueEmbedding(body, ArrayId::Aimnet2Embedding, atom, row, dims);
    out.embeddingDims = dims;
    if (out.embedding && dims == 256) ++stats->embeddingPresent;
    return out;
}

}  // namespace

bool EnsureRowDesignRingArrays(RunData& run, QString* err_out) {
    const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
    if (atomCount == 0) {
        if (err_out) *err_out = QStringLiteral("row_design ring restore requires a non-empty protein");
        return false;
    }

    if (run.poseKind() == PoseKind::Trajectory) {
        const io::QtTrajectoryH5* h5 = run.h5();
        if (!h5 || !h5->bsShielding() || !h5->hmShielding()) {
            if (err_out) *err_out = QStringLiteral("row_design ring restore requires bs_shielding and hm_shielding in the trajectory H5");
            return false;
        }
        for (const RingArrayRequirement& req : kRequiredRingArrays) {
            if (!staticArray(run, req.id)
                && !loadTrajectoryRingArray(run, req, atomCount, err_out)) {
                return false;
            }
            if (!validateResidentArray(run, req, atomCount, err_out)) return false;
        }
        for (const TrajectoryNpyRequirement& req : kTrajectorySourceArrays) {
            if (!loadTrajectorySourceArray(run, req, atomCount, err_out)) return false;
        }
        return true;
    }

    const StaticNpyArray* bs = staticArray(run, ArrayId::KernelBs);
    const StaticNpyArray* hm = staticArray(run, ArrayId::HmShielding);
    if (!bs || bs->rows != atomCount || bs->cols < 9) {
        if (err_out) {
            *err_out = QStringLiteral("row_design ring restore requires static %1.npy with shape (N,9)")
                           .arg(fieldStem(io::FieldSpecFor(io::FieldKind::BSShielding)));
        }
        return false;
    }
    if (!hm || hm->rows != atomCount || hm->cols < 9) {
        if (err_out) {
            *err_out = QStringLiteral("row_design ring restore requires static %1.npy with shape (N,9)")
                           .arg(fieldStem(io::FieldSpecFor(io::FieldKind::HMShielding)));
        }
        return false;
    }
    for (const RingArrayRequirement& req : kRequiredRingArrays)
        if (!validateResidentArray(run, req, atomCount, err_out)) return false;
    return true;
}

RowDesignStats RunRowDesignEmit(const Body& body,
                                RowDesignSink& sink,
                                const ConditioningSpec&) {
    RowDesignStats stats;
    stats.atoms = body.run.protein ? body.run.protein->atomCount() : 0;
    stats.dftRows = body.run.frameMap.dftRows().size();

    auto recordFn = [&](const Body& b, std::size_t atom, std::size_t row,
                        std::size_t orig, const FrameResult&) {
        return buildRow(b, atom, row, orig, &stats);
    };
    auto sourceFn = [](const Body&, const AtomState&, const FrameResult&, const RowDesignRow&) {
        return std::vector<int>{};
    };
    auto carrier = [&](std::size_t, std::size_t, std::size_t, const FrameResult&,
                       const RowDesignRow& row, const std::vector<int>&) {
        if (!sink.WriteRow(row)) throw std::runtime_error("row_design sink write failed");
        ++stats.rows;
    };

    RunTraversal(body, allAtoms, labFrame, recordFn, sourceFn, carrier);
    return stats;
}

}  // namespace h5reader::rediscover
