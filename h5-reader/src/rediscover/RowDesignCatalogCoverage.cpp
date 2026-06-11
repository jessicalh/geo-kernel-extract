#include "RowDesignCatalogCoverage.h"

#include "RowDesignEmitter.h"

#include "../io/QtNpyReader.h"
#include "../model/DftShielding.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QTextStream>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace h5reader::rediscover {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

QString qsv(std::string_view s) {
    return QString::fromUtf8(s.data(), static_cast<qsizetype>(s.size()));
}

QString fieldStem(const io::FieldSpec& spec) { return qsv(spec.stem); }

QString fieldStem(io::FieldKind kind) { return fieldStem(io::FieldSpecFor(kind)); }

QString groupName(io::FieldGroup g) {
    switch (g) {
    case io::FieldGroup::Identity: return QStringLiteral("identity");
    case io::FieldGroup::Enrichment: return QStringLiteral("enrichment");
    case io::FieldGroup::ChargeAssignment: return QStringLiteral("charge_assignment");
    case io::FieldGroup::BiotSavart: return QStringLiteral("biot_savart");
    case io::FieldGroup::HaighMallion: return QStringLiteral("haigh_mallion");
    case io::FieldGroup::PiQuadrupole: return QStringLiteral("pi_quadrupole");
    case io::FieldGroup::Dispersion: return QStringLiteral("dispersion");
    case io::FieldGroup::McConnell: return QStringLiteral("mcconnell");
    case io::FieldGroup::Coulomb: return QStringLiteral("coulomb");
    case io::FieldGroup::HBond: return QStringLiteral("hbond");
    case io::FieldGroup::DSSP: return QStringLiteral("dssp");
    case io::FieldGroup::SASA: return QStringLiteral("sasa");
    case io::FieldGroup::WaterField: return QStringLiteral("water_field");
    case io::FieldGroup::Hydration: return QStringLiteral("hydration");
    case io::FieldGroup::WaterPolarization: return QStringLiteral("water_polarization");
    case io::FieldGroup::EEQ: return QStringLiteral("eeq");
    case io::FieldGroup::Gromacs: return QStringLiteral("gromacs");
    case io::FieldGroup::Bonded: return QStringLiteral("bonded");
    case io::FieldGroup::MOPACCore: return QStringLiteral("mopac_core");
    case io::FieldGroup::MOPACCoulomb: return QStringLiteral("mopac_coulomb");
    case io::FieldGroup::APBS: return QStringLiteral("apbs");
    case io::FieldGroup::Orca: return QStringLiteral("orca");
    case io::FieldGroup::AIMNet2: return QStringLiteral("aimnet2");
    case io::FieldGroup::PlanarGeometry: return QStringLiteral("planar_geometry");
    case io::FieldGroup::Topology: return QStringLiteral("topology");
    default: return QStringLiteral("excluded");
    }
}

QString axisName(io::NativeAxis a) {
    switch (a) {
    case io::NativeAxis::Atom: return QStringLiteral("atom");
    case io::NativeAxis::RingContributionPair: return QStringLiteral("ring_contribution_pair");
    case io::NativeAxis::AromaticRing: return QStringLiteral("aromatic_ring");
    case io::NativeAxis::Protein: return QStringLiteral("protein");
    case io::NativeAxis::Bond: return QStringLiteral("bond");
    case io::NativeAxis::MOPACBondNeighborPair: return QStringLiteral("mopac_bond_neighbor_pair");
    case io::NativeAxis::MOPACUniquePair: return QStringLiteral("mopac_unique_pair");
    case io::NativeAxis::Residue: return QStringLiteral("residue");
    case io::NativeAxis::SaturatedRing: return QStringLiteral("saturated_ring");
    case io::NativeAxis::Ring: return QStringLiteral("ring");
    case io::NativeAxis::RingMembership: return QStringLiteral("ring_membership");
    default: return QStringLiteral("excluded");
    }
}

bool excludedGroup(io::FieldGroup g) {
    switch (g) {
    case io::FieldGroup::LarsenHBond:
    case io::FieldGroup::Tripeptide:
    case io::FieldGroup::Delta:
    case io::FieldGroup::Rediscover:
    case io::FieldGroup::McConnellLegacy:
    case io::FieldGroup::MOPACMcConnellLegacy:
    case io::FieldGroup::MOPACCoulombLegacy:
    case io::FieldGroup::CoulombLegacy:
        return true;
    default:
        return false;
    }
}

bool scoped(const io::FieldSpec& spec) {
    if (excludedGroup(spec.group)) return false;
    return true;
}

bool structuredField(const io::FieldSpec& spec) {
    return spec.kind == io::FieldKind::AtomsCategoryInfo
           || spec.group == io::FieldGroup::Topology;
}

std::size_t nominalComponentCount(const io::FieldSpec& spec) {
    return spec.cols > 0 ? static_cast<std::size_t>(spec.cols) : 1;
}

bool scalarAtomRowField(const io::FieldSpec& spec) {
    if (!scoped(spec) || structuredField(spec)) return false;
    if (spec.axis != io::NativeAxis::Atom) return false;
    if (spec.tensor_rank != 0) return false;
    const std::size_t cols = nominalComponentCount(spec);
    return cols <= 12;
}

QString rowColumnName(const io::FieldSpec& spec, std::size_t component) {
    const QString stem = fieldStem(spec);
    if (nominalComponentCount(spec) == 1) return QStringLiteral("catalog_%1").arg(stem);
    return QStringLiteral("catalog_%1_%2").arg(stem).arg(component);
}

RowNativeAxis rowAxis(io::NativeAxis a) {
    switch (a) {
    case io::NativeAxis::Atom: return RowNativeAxis::Atom;
    case io::NativeAxis::Residue: return RowNativeAxis::Residue;
    case io::NativeAxis::Protein: return RowNativeAxis::Protein;
    case io::NativeAxis::Bond: return RowNativeAxis::Bond;
    case io::NativeAxis::Ring:
    case io::NativeAxis::AromaticRing:
    case io::NativeAxis::SaturatedRing:
    case io::NativeAxis::RingContributionPair:
    case io::NativeAxis::RingMembership:
        return RowNativeAxis::Ring;
    default:
        return RowNativeAxis::Row;
    }
}

QByteArray npyHeaderLiteral(const QByteArray& descr, const std::vector<std::size_t>& shape) {
    QByteArray header;
    header += "{'descr': ";
    header += descr;
    header += ", 'fortran_order': False, 'shape': (";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (i) header += ", ";
        header += QByteArray::number(static_cast<qulonglong>(shape[i]));
    }
    if (shape.size() == 1) header += ",";
    header += "), }";
    constexpr int kPreambleBytes = 10;
    const int newlineBytes = 1;
    int pad = 16 - ((kPreambleBytes + header.size() + newlineBytes) % 16);
    if (pad == 16) pad = 0;
    header += QByteArray(pad, ' ');
    header += '\n';

    QByteArray prefix;
    prefix.append("\x93NUMPY", 6);
    prefix.append(char(1));
    prefix.append(char(0));
    const quint16 headerLen = static_cast<quint16>(header.size());
    prefix.append(char(headerLen & 0xff));
    prefix.append(char((headerLen >> 8) & 0xff));
    return prefix + header;
}

QByteArray npyHeader(const QByteArray& descr, const std::vector<std::size_t>& shape) {
    QByteArray literal;
    literal += "'";
    literal += descr;
    literal += "'";
    return npyHeaderLiteral(literal, shape);
}

class NpyDoubleWriter {
public:
    NpyDoubleWriter(const QString& path, std::size_t rows, std::size_t cols)
        : file_(path), cols_(cols) {
        ok_ = file_.open(QIODevice::WriteOnly);
        if (!ok_) return;
        const QByteArray header = npyHeader(QByteArrayLiteral("<f8"), {rows, cols});
        ok_ = file_.write(header) == header.size();
    }

    bool writeRow(const std::vector<double>& row) {
        if (!ok_ || row.size() != cols_) return false;
        const qsizetype n = static_cast<qsizetype>(row.size() * sizeof(double));
        return file_.write(reinterpret_cast<const char*>(row.data()), n) == n;
    }

    bool commit() { return ok_ && file_.commit(); }

private:
    QSaveFile file_;
    std::size_t cols_ = 0;
    bool ok_ = false;
};

bool writeText(const QString& path, const QByteArray& bytes, QString* err_out) {
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write %1").arg(path);
        return false;
    }
    if (f.write(bytes) != bytes.size()) {
        if (err_out) *err_out = QStringLiteral("short write to %1").arg(path);
        return false;
    }
    return true;
}

QString trajectoryFrameNpyPath(const RunData& run, std::size_t row, io::FieldKind kind) {
    if (!run.manifest.trajectory) return {};
    const QString npysDir = QDir(run.manifest.trajectory->extraction_dir_abspath)
                                .filePath(QStringLiteral("npys"));
    const std::size_t original = run.frameMap.originalIndex(row);
    const QString frameDir = QDir(npysDir).filePath(
        QStringLiteral("frame_%1").arg(original, 6, 10, QLatin1Char('0')));
    return QDir(frameDir).filePath(fieldStem(kind) + QStringLiteral(".npy"));
}

QString staticFieldPath(const RunData& run, io::FieldKind kind) {
    if (run.manifest.single_pose)
        return QDir(run.manifest.single_pose->pose_dir_abspath)
            .filePath(fieldStem(kind) + QStringLiteral(".npy"));
    if (run.manifest.trajectory)
        return QDir(run.manifest.trajectory->extraction_dir_abspath)
            .filePath(fieldStem(kind) + QStringLiteral(".npy"));
    return {};
}

bool readTrajectoryArray(const RunData& run,
                         std::size_t row,
                         io::FieldKind kind,
                         io::QtNpyReader::WidenedArray* out) {
    const QString path = trajectoryFrameNpyPath(run, row, kind);
    if (path.isEmpty() || !QFileInfo::exists(path)) return false;
    *out = io::QtNpyReader::ReadArrayWidened(path);
    return out->ok;
}

bool readRunArray(const RunData& run,
                  std::size_t row,
                  io::FieldKind kind,
                  io::QtNpyReader::WidenedArray* out) {
    if (run.poseKind() == PoseKind::Trajectory) return readTrajectoryArray(run, row, kind, out);
    const StaticNpyArray* a = run.producerArray(kind);
    if (!a) return false;
    out->ok = true;
    out->rows = a->rows;
    out->cols = a->cols;
    out->descr = a->dtype_descr.toStdString();
    out->data = a->values;
    return true;
}

void normalizeNativeArray(const io::FieldSpec& spec, io::QtNpyReader::WidenedArray* a) {
    if (!a || !a->ok) return;
    if (spec.axis == io::NativeAxis::Protein
        && spec.cols > 1
        && a->cols == 1
        && a->rows == static_cast<std::size_t>(spec.cols)) {
        a->rows = 1;
        a->cols = static_cast<std::size_t>(spec.cols);
    }
}

bool readRunArrayForSpec(const RunData& run,
                         std::size_t row,
                         const io::FieldSpec& spec,
                         io::QtNpyReader::WidenedArray* out) {
    if (!readRunArray(run, row, spec.kind, out)) return false;
    normalizeNativeArray(spec, out);
    return true;
}

bool hasNumericArray(const RunData& run, io::FieldKind kind) {
    if (run.poseKind() == PoseKind::Static) return run.producerArray(kind) != nullptr;
    for (std::size_t row : run.frameMap.dftRows()) {
        const QString path = trajectoryFrameNpyPath(run, row, kind);
        if (QFileInfo::exists(path)) return true;
    }
    return false;
}

bool isDftKind(io::FieldKind kind) {
    return kind == io::FieldKind::OrcaTotal
           || kind == io::FieldKind::OrcaDiamagnetic
           || kind == io::FieldKind::OrcaParamagnetic;
}

const model::Mat3* dftMatrix(const model::DftAtomShielding& atom, io::FieldKind kind) {
    if (kind == io::FieldKind::OrcaTotal) return &atom.total_raw;
    if (kind == io::FieldKind::OrcaDiamagnetic) return &atom.dia_raw;
    if (kind == io::FieldKind::OrcaParamagnetic) return &atom.para_raw;
    return nullptr;
}

void countFinite(const std::vector<double>& row, std::vector<std::size_t>* counts, std::size_t offset = 0) {
    for (std::size_t i = offset; i < row.size(); ++i) {
        if (std::isfinite(row[i])) ++(*counts)[i - offset];
    }
}

std::size_t emittedRowsForRun(const RunData& run) {
    return (run.protein ? run.protein->atomCount() : 0) * run.frameMap.dftRows().size();
}

QJsonObject baseFieldJson(const io::FieldSpec& spec) {
    QJsonObject o;
    o.insert(QStringLiteral("kind"), QString::fromLatin1("FieldKind::") + QString::number(static_cast<int>(spec.kind)));
    o.insert(QStringLiteral("stem"), fieldStem(spec));
    o.insert(QStringLiteral("group"), groupName(spec.group));
    o.insert(QStringLiteral("native_axis"), axisName(spec.axis));
    o.insert(QStringLiteral("catalog_cols"), spec.cols);
    o.insert(QStringLiteral("tensor_rank"), spec.tensor_rank);
    o.insert(QStringLiteral("units"), qsv(spec.units));
    o.insert(QStringLiteral("irreps"), qsv(spec.irreps));
    return o;
}

QJsonArray countsJson(const std::vector<std::size_t>& counts) {
    QJsonArray a;
    for (std::size_t c : counts) a.push_back(static_cast<double>(c));
    return a;
}

QJsonArray stringJson(const std::vector<QString>& values) {
    QJsonArray a;
    for (const QString& v : values) a.push_back(v);
    return a;
}

std::size_t totalEmittedRows(const std::vector<RunData>& runs) {
    std::size_t n = 0;
    for (const RunData& run : runs) n += emittedRowsForRun(run);
    return n;
}

int primaryNativeRingForAtom(const RunData& run, std::size_t atom, io::NativeAxis axis) {
    if (!run.protein) return -1;
    const model::QtTopology& topo = run.protein->topology();
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
        const bool wanted =
            (axis == io::NativeAxis::AromaticRing && ring.ringKind == model::RingKind::Aromatic)
            || (axis == io::NativeAxis::SaturatedRing && ring.ringKind == model::RingKind::Saturated)
            || (axis == io::NativeAxis::Ring);
        if (wanted && ring.nativeAxisIndex >= 0) return ring.nativeAxisIndex;
    }
    return -1;
}

std::vector<QString> directReductionNames(std::size_t cols, const QString& prefix) {
    std::vector<QString> names;
    names.reserve(cols);
    for (std::size_t c = 0; c < cols; ++c)
        names.push_back(QStringLiteral("%1_%2").arg(prefix).arg(c));
    return names;
}

std::vector<QString> ringPairReductionNames(std::size_t cols) {
    std::vector<QString> names;
    names.reserve(1 + cols * 10);
    names.push_back(QStringLiteral("pair_count"));
    for (std::size_t c = 0; c < cols; ++c) {
        names.push_back(QStringLiteral("sum_%1").arg(c));
        names.push_back(QStringLiteral("nearest_%1").arg(c));
        for (int t = 0; t < 8; ++t)
            names.push_back(QStringLiteral("type%1_sum_%2").arg(t).arg(c));
    }
    return names;
}

std::vector<QString> incidentReductionNames(std::size_t cols, bool includeMeanMax) {
    std::vector<QString> names;
    names.reserve(1 + cols * (includeMeanMax ? 3 : 1));
    names.push_back(QStringLiteral("incident_count"));
    for (std::size_t c = 0; c < cols; ++c) {
        names.push_back(QStringLiteral("sum_%1").arg(c));
        if (includeMeanMax) {
            names.push_back(QStringLiteral("mean_%1").arg(c));
            names.push_back(QStringLiteral("max_%1").arg(c));
        }
    }
    return names;
}

void countReductionRow(const std::vector<double>& row, std::vector<std::size_t>* counts) {
    for (std::size_t i = 0; i < row.size(); ++i)
        if (std::isfinite(row[i])) ++(*counts)[i];
}

struct ReductionSpec {
    QString policy;
    std::vector<QString> names;
    bool flaggedForReview = false;
};

ReductionSpec reductionSpecFor(const io::FieldSpec& spec, std::size_t cols) {
    switch (spec.axis) {
    case io::NativeAxis::Protein:
        return {QStringLiteral("broadcast protein-axis values to every emitted atom row"), directReductionNames(cols, QStringLiteral("broadcast")), false};
    case io::NativeAxis::Residue:
        return {QStringLiteral("broadcast residue-axis values by atom.residue_index"), directReductionNames(cols, QStringLiteral("residue_broadcast")), false};
    case io::NativeAxis::AromaticRing:
    case io::NativeAxis::SaturatedRing:
    case io::NativeAxis::Ring:
        return {QStringLiteral("primary ring-membership policy: first topology membership matching the catalog native axis"), directReductionNames(cols, QStringLiteral("primary_membership")), false};
    case io::NativeAxis::RingContributionPair:
        return {QStringLiteral("per-(atom,ring): pair count, sum over rings, nearest-ring value by ring_contributions distance, and ring-type sum bins"), ringPairReductionNames(cols), false};
    case io::NativeAxis::MOPACBondNeighborPair:
        return {QStringLiteral("mopac_bond_neighbors incident aggregate: count, sum, mean, max by atom endpoint"), incidentReductionNames(cols, true), false};
    case io::NativeAxis::Bond:
    case io::NativeAxis::MOPACUniquePair:
        return {QStringLiteral("generic incident endpoint sum for bond-like native axis; review downstream semantics before modelling"), incidentReductionNames(cols, false), true};
    default:
        return {QStringLiteral("generic native-axis sum reduction; review downstream semantics before modelling"), incidentReductionNames(cols, false), true};
    }
}

bool writeStructuredIndex(const QString& outDir,
                          const io::FieldSpec& spec,
                          const std::vector<RunData>& runs,
                          QJsonObject* field,
                          QString* err_out) {
    const QString indexName = QStringLiteral("row_design_field_%1_index.csv").arg(fieldStem(spec));
    const QString sidecarName = QStringLiteral("row_design_field_%1_native.npy").arg(fieldStem(spec));

    struct Chunk {
        QString path;
        QString proteinId;
        std::size_t rowStart = 0;
        std::size_t rowCount = 0;
    };

    std::vector<Chunk> chunks;
    chunks.reserve(runs.size());
    std::size_t total = 0;
    std::size_t recordSize = 0;
    std::string descr;

    for (const RunData& run : runs) {
        const QString path = staticFieldPath(run, spec.kind);
        std::vector<unsigned char> bytes;
        std::size_t rs = 0;
        auto r = io::QtNpyReader::ReadRawBytes(path, bytes, rs);
        if (!r.ok) {
            if (err_out) *err_out = QStringLiteral("required structured field %1 is missing or unreadable: %2")
                                      .arg(fieldStem(spec), path);
            return false;
        }
        if (recordSize == 0) {
            recordSize = rs;
            descr = r.dtype_descr;
        } else if (recordSize != rs || descr != r.dtype_descr) {
            if (err_out) *err_out = QStringLiteral("structured field %1 dtype drift at %2")
                                      .arg(fieldStem(spec), path);
            return false;
        }
        chunks.push_back(Chunk{path, run.protein ? run.protein->proteinId() : QString(), total, r.row_count});
        total += r.row_count;
    }

    QSaveFile rawFile(QDir(outDir).filePath(sidecarName));
    if (!rawFile.open(QIODevice::WriteOnly)) {
        if (err_out) *err_out = QStringLiteral("cannot write %1").arg(sidecarName);
        return false;
    }
    const QByteArray header = npyHeaderLiteral(QByteArray::fromStdString(descr), {total});
    if (rawFile.write(header) != header.size()) {
        if (err_out) *err_out = QStringLiteral("short write to %1").arg(sidecarName);
        return false;
    }
    for (const Chunk& chunk : chunks) {
        std::vector<unsigned char> bytes;
        std::size_t rs = 0;
        auto r = io::QtNpyReader::ReadRawBytes(chunk.path, bytes, rs);
        if (!r.ok || rs != recordSize || r.dtype_descr != descr) {
            if (err_out) *err_out = QStringLiteral("structured field %1 changed while writing %2")
                                      .arg(fieldStem(spec), chunk.path);
            return false;
        }
        const qsizetype n = static_cast<qsizetype>(bytes.size());
        if (rawFile.write(reinterpret_cast<const char*>(bytes.data()), n) != n) {
            if (err_out) *err_out = QStringLiteral("short raw write to %1").arg(sidecarName);
            return false;
        }
    }
    if (!rawFile.commit()) {
        if (err_out) *err_out = QStringLiteral("structured sidecar commit failed for %1").arg(sidecarName);
        return false;
    }

    QFile indexFile(QDir(outDir).filePath(indexName));
    if (!indexFile.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write %1").arg(indexName);
        return false;
    }
    QTextStream ts(&indexFile);
    ts << "protein_ordinal,protein_id,source_path,row_start,row_count\n";

    for (std::size_t i = 0; i < chunks.size(); ++i) {
        const Chunk& chunk = chunks[i];
        ts << i << "," << chunk.proteinId
           << "," << chunk.path << "," << chunk.rowStart << "," << chunk.rowCount << "\n";
    }
    (*field).insert(QStringLiteral("representation"), QStringLiteral("structured_native_sidecar"));
    (*field).insert(QStringLiteral("sidecar"), sidecarName);
    (*field).insert(QStringLiteral("index_sidecar"), indexName);
    (*field).insert(QStringLiteral("shape"), QStringLiteral("(%1,)").arg(total));
    (*field).insert(QStringLiteral("component_names"), QJsonArray{QStringLiteral("record")});
    (*field).insert(QStringLiteral("populated_counts"), QJsonArray{static_cast<double>(total)});
    (*field).insert(QStringLiteral("row_count"), static_cast<double>(total));
    return true;
}

bool writeAtomSidecar(const QString& outDir,
                      const io::FieldSpec& spec,
                      const std::vector<RunData>& runs,
                      QJsonObject* field,
                      QString* err_out) {
    const std::size_t cols = nominalComponentCount(spec);
    bool any = isDftKind(spec.kind);
    for (const RunData& run : runs) any = any || hasNumericArray(run, spec.kind);
    const std::size_t totalRows = [&]() {
        std::size_t n = 0;
        for (const RunData& run : runs) n += emittedRowsForRun(run);
        return n;
    }();
    std::vector<std::size_t> counts(cols, 0);
    if (!any) {
        (*field).insert(QStringLiteral("representation"), QStringLiteral("atom_row_sidecar_absent"));
        (*field).insert(QStringLiteral("sidecar"), QJsonValue());
        (*field).insert(QStringLiteral("shape"), QStringLiteral("(rows,%1)").arg(cols));
        (*field).insert(QStringLiteral("populated_counts"), countsJson(counts));
        (*field).insert(QStringLiteral("row_count"), static_cast<double>(totalRows));
        return true;
    }

    const QString name = QStringLiteral("row_design_field_%1.npy").arg(fieldStem(spec));
    NpyDoubleWriter writer(QDir(outDir).filePath(name), totalRows, cols);
    std::vector<double> values(cols, kNaN);
    for (const RunData& run : runs) {
        const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
        for (std::size_t row : run.frameMap.dftRows()) {
            io::QtNpyReader::WidenedArray arr;
            const bool arrOk = !isDftKind(spec.kind) && readRunArray(run, row, spec.kind, &arr);
            for (std::size_t atom = 0; atom < atomCount; ++atom) {
                std::fill(values.begin(), values.end(), kNaN);
                if (isDftKind(spec.kind)) {
                    const model::DftAtomShielding* dft =
                        run.dft.AtomShielding(atom, run.frameMap.originalIndex(row));
                    if (dft) {
                        if (const model::Mat3* m = dftMatrix(*dft, spec.kind)) {
                            for (int r = 0; r < 3; ++r)
                                for (int c = 0; c < 3; ++c)
                                    values[static_cast<std::size_t>(r * 3 + c)] = (*m)(r, c);
                        }
                    }
                } else if (arrOk && atom < arr.rows) {
                    const std::size_t copyCols = std::min(cols, arr.cols);
                    for (std::size_t c = 0; c < copyCols; ++c)
                        values[c] = arr.data[atom * arr.cols + c];
                }
                countFinite(values, &counts);
                if (!writer.writeRow(values)) {
                    if (err_out) *err_out = QStringLiteral("sidecar write failed for %1").arg(name);
                    return false;
                }
            }
        }
    }
    if (!writer.commit()) {
        if (err_out) *err_out = QStringLiteral("sidecar commit failed for %1").arg(name);
        return false;
    }
    (*field).insert(QStringLiteral("representation"), QStringLiteral("atom_row_sidecar"));
    (*field).insert(QStringLiteral("sidecar"), name);
    (*field).insert(QStringLiteral("shape"), QStringLiteral("(%1,%2)").arg(totalRows).arg(cols));
    (*field).insert(QStringLiteral("populated_counts"), countsJson(counts));
    (*field).insert(QStringLiteral("row_count"), static_cast<double>(totalRows));
    return true;
}

bool writeNativeReductionSidecar(const QString& outDir,
                                  const io::FieldSpec& spec,
                                  std::size_t dataCols,
                                  const std::vector<RunData>& runs,
                                  QJsonObject* field,
                                  QString* err_out) {
    const std::size_t totalRows = totalEmittedRows(runs);
    const ReductionSpec red = reductionSpecFor(spec, dataCols);
    const QString name = QStringLiteral("row_design_field_%1_reduction.npy").arg(fieldStem(spec));
    NpyDoubleWriter writer(QDir(outDir).filePath(name), totalRows, red.names.size());
    std::vector<std::size_t> counts(red.names.size(), 0);
    std::vector<double> out(red.names.size(), kNaN);

    auto writeOut = [&]() -> bool {
        countReductionRow(out, &counts);
        if (!writer.writeRow(out)) {
            if (err_out) *err_out = QStringLiteral("native reduction write failed for %1").arg(name);
            return false;
        }
        return true;
    };

    auto writeDirectRows = [&](const RunData& run,
                               std::size_t frame,
                               const io::QtNpyReader::WidenedArray& arr) -> bool {
        const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
        for (std::size_t atom = 0; atom < atomCount; ++atom) {
            std::fill(out.begin(), out.end(), kNaN);
            std::size_t native = 0;
            bool ok = false;
            if (spec.axis == io::NativeAxis::Protein) {
                native = 0;
                ok = arr.rows > 0;
            } else if (spec.axis == io::NativeAxis::Residue && run.protein) {
                const model::QtAtom& a = run.protein->atom(atom);
                ok = a.residueIndex >= 0;
                native = ok ? static_cast<std::size_t>(a.residueIndex) : 0;
            } else if (spec.axis == io::NativeAxis::AromaticRing
                       || spec.axis == io::NativeAxis::SaturatedRing
                       || spec.axis == io::NativeAxis::Ring) {
                const int idx = primaryNativeRingForAtom(run, atom, spec.axis);
                ok = idx >= 0;
                native = ok ? static_cast<std::size_t>(idx) : 0;
            }
            if (ok && native < arr.rows) {
                for (std::size_t c = 0; c < dataCols && c < arr.cols; ++c)
                    out[c] = arr.data[native * arr.cols + c];
            }
            Q_UNUSED(frame);
            if (!writeOut()) return false;
        }
        return true;
    };

    auto writeRingPairRows = [&](const RunData& run,
                                 const io::QtNpyReader::WidenedArray& arr,
                                 const io::QtNpyReader::WidenedArray& map) -> bool {
        const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
        std::vector<std::vector<std::size_t>> byAtom(atomCount);
        for (std::size_t native = 0; native < map.rows; ++native) {
            const double atomV = map.data[native * map.cols + 0];
            if (!std::isfinite(atomV)) continue;
            const long atomIdx = std::lround(atomV);
            if (atomIdx >= 0 && static_cast<std::size_t>(atomIdx) < atomCount)
                byAtom[static_cast<std::size_t>(atomIdx)].push_back(native);
        }
        for (std::size_t atom = 0; atom < atomCount; ++atom) {
            std::fill(out.begin(), out.end(), kNaN);
            out[0] = static_cast<double>(byAtom[atom].size());
            std::size_t bestNative = std::numeric_limits<std::size_t>::max();
            double bestDist = kNaN;
            for (std::size_t native : byAtom[atom]) {
                const double dist = map.cols > 3 ? map.data[native * map.cols + 3] : kNaN;
                if (std::isfinite(dist) && (std::isnan(bestDist) || dist < bestDist)) {
                    bestDist = dist;
                    bestNative = native;
                }
            }
            for (std::size_t c = 0; c < dataCols; ++c) {
                const std::size_t base = 1 + c * 10;
                double sum = 0.0;
                std::array<double, 8> typeSums = {};
                bool any = false;
                std::array<bool, 8> typeAny = {};
                for (std::size_t native : byAtom[atom]) {
                    if (native >= arr.rows || c >= arr.cols) continue;
                    const double v = arr.data[native * arr.cols + c];
                    if (!std::isfinite(v)) continue;
                    any = true;
                    sum += v;
                    const double typeV = map.cols > 2 ? map.data[native * map.cols + 2] : kNaN;
                    if (std::isfinite(typeV)) {
                        const long t = std::lround(typeV);
                        if (t >= 0 && t < 8) {
                            typeSums[static_cast<std::size_t>(t)] += v;
                            typeAny[static_cast<std::size_t>(t)] = true;
                        }
                    }
                }
                if (any) out[base] = sum;
                if (bestNative != std::numeric_limits<std::size_t>::max()
                    && bestNative < arr.rows
                    && c < arr.cols) {
                    const double v = arr.data[bestNative * arr.cols + c];
                    if (std::isfinite(v)) out[base + 1] = v;
                }
                for (std::size_t t = 0; t < typeSums.size(); ++t)
                    if (typeAny[t]) out[base + 2 + t] = typeSums[t];
            }
            if (!writeOut()) return false;
        }
        return true;
    };

    auto writeIncidentRows = [&](const RunData& run,
                                 const io::QtNpyReader::WidenedArray& arr,
                                 bool includeMeanMax) -> bool {
        const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
        std::vector<std::vector<std::size_t>> byAtom(atomCount);
        for (std::size_t native = 0; native < arr.rows; ++native) {
            auto addEndpoint = [&](std::size_t comp) {
                if (comp >= arr.cols) return;
                const double atomV = arr.data[native * arr.cols + comp];
                if (!std::isfinite(atomV)) return;
                const long atomIdx = std::lround(atomV);
                if (atomIdx >= 0 && static_cast<std::size_t>(atomIdx) < atomCount)
                    byAtom[static_cast<std::size_t>(atomIdx)].push_back(native);
            };
            addEndpoint(0);
            if (spec.axis != io::NativeAxis::MOPACBondNeighborPair) addEndpoint(1);
        }
        for (std::size_t atom = 0; atom < atomCount; ++atom) {
            std::fill(out.begin(), out.end(), kNaN);
            out[0] = static_cast<double>(byAtom[atom].size());
            for (std::size_t c = 0; c < dataCols; ++c) {
                const std::size_t base = 1 + c * (includeMeanMax ? 3 : 1);
                double sum = 0.0;
                double maxv = kNaN;
                std::size_t n = 0;
                for (std::size_t native : byAtom[atom]) {
                    if (native >= arr.rows || c >= arr.cols) continue;
                    const double v = arr.data[native * arr.cols + c];
                    if (!std::isfinite(v)) continue;
                    sum += v;
                    maxv = std::isnan(maxv) ? v : std::max(maxv, v);
                    ++n;
                }
                if (n > 0) {
                    out[base] = sum;
                    if (includeMeanMax) {
                        out[base + 1] = sum / static_cast<double>(n);
                        out[base + 2] = maxv;
                    }
                }
            }
            if (!writeOut()) return false;
        }
        return true;
    };

    for (const RunData& run : runs) {
        for (std::size_t frame : run.frameMap.dftRows()) {
            io::QtNpyReader::WidenedArray arr;
            if (!readRunArrayForSpec(run, frame, spec, &arr)) {
                const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
                std::fill(out.begin(), out.end(), kNaN);
                for (std::size_t atom = 0; atom < atomCount; ++atom)
                    if (!writeOut()) return false;
                continue;
            }
            if (arr.cols != dataCols) {
                if (err_out) {
                    *err_out = QStringLiteral("%1 reduction width drift: expected %2 cols, saw %3")
                                   .arg(fieldStem(spec))
                                   .arg(dataCols)
                                   .arg(arr.cols);
                }
                return false;
            }
            if (spec.axis == io::NativeAxis::Protein
                || spec.axis == io::NativeAxis::Residue
                || spec.axis == io::NativeAxis::AromaticRing
                || spec.axis == io::NativeAxis::SaturatedRing
                || spec.axis == io::NativeAxis::Ring) {
                if (!writeDirectRows(run, frame, arr)) return false;
            } else if (spec.axis == io::NativeAxis::RingContributionPair) {
                io::QtNpyReader::WidenedArray map;
                if (!readRunArrayForSpec(run, frame, io::FieldSpecFor(io::FieldKind::RingContributions), &map)
                    || map.cols < 3
                    || map.rows != arr.rows) {
                    if (err_out) {
                        *err_out = QStringLiteral("%1 reduction requires ring_contributions map with matching rows")
                                       .arg(fieldStem(spec));
                    }
                    return false;
                }
                if (!writeRingPairRows(run, arr, map)) return false;
            } else if (spec.axis == io::NativeAxis::MOPACBondNeighborPair) {
                if (!writeIncidentRows(run, arr, true)) return false;
            } else if (spec.axis == io::NativeAxis::Bond
                       || spec.axis == io::NativeAxis::MOPACUniquePair) {
                if (!writeIncidentRows(run, arr, false)) return false;
            } else {
                if (!writeIncidentRows(run, arr, false)) return false;
            }
        }
    }

    if (!writer.commit()) {
        if (err_out) *err_out = QStringLiteral("native reduction commit failed for %1").arg(name);
        return false;
    }

    QJsonObject r;
    r.insert(QStringLiteral("representation"), QStringLiteral("atom_row_reduction_sidecar"));
    r.insert(QStringLiteral("sidecar"), name);
    r.insert(QStringLiteral("shape"), QStringLiteral("(%1,%2)").arg(totalRows).arg(red.names.size()));
    r.insert(QStringLiteral("policy"), red.policy);
    r.insert(QStringLiteral("component_names"), stringJson(red.names));
    r.insert(QStringLiteral("populated_counts"), countsJson(counts));
    r.insert(QStringLiteral("row_count"), static_cast<double>(totalRows));
    if (red.flaggedForReview) r.insert(QStringLiteral("flagged_for_review"), true);
    (*field).insert(QStringLiteral("row_reduction"), r);
    return true;
}

bool writeNativeNumericSidecar(const QString& outDir,
                               const io::FieldSpec& spec,
                               const std::vector<RunData>& runs,
                               QJsonObject* field,
                               QString* err_out) {
    std::size_t dataCols = nominalComponentCount(spec);
    std::size_t totalRows = 0;
    bool any = false;
    for (const RunData& run : runs) {
        for (std::size_t row : run.frameMap.dftRows()) {
            io::QtNpyReader::WidenedArray arr;
            if (!readRunArrayForSpec(run, row, spec, &arr)) continue;
            any = true;
            dataCols = arr.cols;
            totalRows += arr.rows;
        }
    }
    std::vector<std::size_t> counts(dataCols, 0);
    if (!any) {
        (*field).insert(QStringLiteral("representation"), QStringLiteral("native_axis_sidecar_absent"));
        (*field).insert(QStringLiteral("sidecar"), QJsonValue());
        (*field).insert(QStringLiteral("shape"), QStringLiteral("(native_rows,%1)").arg(4 + dataCols));
        (*field).insert(QStringLiteral("populated_counts"), countsJson(counts));
        (*field).insert(QStringLiteral("row_count"), 0);
        return true;
    }

    const QString name = QStringLiteral("row_design_field_%1_native.npy").arg(fieldStem(spec));
    NpyDoubleWriter writer(QDir(outDir).filePath(name), totalRows, 4 + dataCols);
    std::vector<double> out(4 + dataCols, kNaN);
    int proteinOrdinal = 0;
    for (const RunData& run : runs) {
        for (std::size_t row : run.frameMap.dftRows()) {
            io::QtNpyReader::WidenedArray arr;
            if (!readRunArrayForSpec(run, row, spec, &arr)) continue;
            for (std::size_t native = 0; native < arr.rows; ++native) {
                out[0] = static_cast<double>(proteinOrdinal);
                out[1] = static_cast<double>(row);
                out[2] = static_cast<double>(run.frameMap.originalIndex(row));
                out[3] = static_cast<double>(native);
                for (std::size_t c = 0; c < dataCols; ++c)
                    out[4 + c] = c < arr.cols ? arr.data[native * arr.cols + c] : kNaN;
                countFinite(out, &counts, 4);
                if (!writer.writeRow(out)) {
                    if (err_out) *err_out = QStringLiteral("native sidecar write failed for %1").arg(name);
                    return false;
                }
            }
        }
        ++proteinOrdinal;
    }
    if (!writer.commit()) {
        if (err_out) *err_out = QStringLiteral("native sidecar commit failed for %1").arg(name);
        return false;
    }
    (*field).insert(QStringLiteral("representation"), QStringLiteral("native_axis_sidecar"));
    (*field).insert(QStringLiteral("sidecar"), name);
    (*field).insert(QStringLiteral("shape"), QStringLiteral("(%1,%2)").arg(totalRows).arg(4 + dataCols));
    (*field).insert(QStringLiteral("context_columns"),
                    QJsonArray{QStringLiteral("protein_ordinal"), QStringLiteral("frame_slot"),
                               QStringLiteral("original_index"), QStringLiteral("native_index")});
    (*field).insert(QStringLiteral("populated_counts"), countsJson(counts));
    (*field).insert(QStringLiteral("row_count"), static_cast<double>(totalRows));
    if (!writeNativeReductionSidecar(outDir, spec, dataCols, runs, field, err_out)) return false;
    return true;
}

}  // namespace

const std::vector<const io::FieldSpec*>& ScopedProducerCatalog() {
    static const std::vector<const io::FieldSpec*> fields = [] {
        std::vector<const io::FieldSpec*> out;
        for (const io::FieldSpec& spec : io::kFieldCatalog)
            if (scoped(spec)) out.push_back(&spec);
        return out;
    }();
    return fields;
}

const std::vector<CatalogRowColumn>& CatalogRowColumns() {
    static const std::vector<CatalogRowColumn> cols = [] {
        std::vector<CatalogRowColumn> out;
        for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
            if (!scalarAtomRowField(*spec)) continue;
            const std::size_t n = nominalComponentCount(*spec);
            for (std::size_t c = 0; c < n; ++c) {
                RowColumnSpec rc;
                rc.name = rowColumnName(*spec, c);
                rc.type = RowColType::Double;
                rc.unit = qsv(spec->units);
                rc.irrep = qsv(spec->irreps);
                rc.nativeAxis = rowAxis(spec->axis);
                rc.timeVarying = true;
                rc.isFeature = spec->is_feature;
                out.push_back(CatalogRowColumn{spec->kind, static_cast<int>(c), rc});
            }
        }
        return out;
    }();
    return cols;
}

bool IsCatalogRowColumn(io::FieldKind kind) {
    const io::FieldSpec& spec = io::FieldSpecFor(kind);
    return scalarAtomRowField(spec);
}

bool EnsureCatalogRowColumnArrays(RunData& run, QString* err_out) {
    if (run.poseKind() != PoseKind::Trajectory) return true;
    const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
    if (atomCount == 0) return true;
    for (const CatalogRowColumn& col : CatalogRowColumns()) {
        if (run.producerArray(col.kind)) continue;
        const io::FieldSpec& spec = io::FieldSpecFor(col.kind);
        io::QtNpyReader::WidenedArray first;
        bool found = false;
        for (std::size_t row : run.frameMap.dftRows()) {
            if (readTrajectoryArray(run, row, col.kind, &first)) {
                found = true;
                break;
            }
        }
        if (!found) continue;
        if (first.rows != atomCount) {
            if (err_out) {
                *err_out = QStringLiteral("%1 row-column source has %2 rows; topology has %3 atoms")
                               .arg(fieldStem(spec))
                               .arg(first.rows)
                               .arg(atomCount);
            }
            return false;
        }
        const std::size_t cols = first.cols;
        StaticNpyArray resident;
        resident.stem = fieldStem(spec);
        resident.path = run.manifest.trajectory
                            ? QDir(run.manifest.trajectory->extraction_dir_abspath).filePath(QStringLiteral("npys"))
                            : QString();
        resident.rows = atomCount * run.frameMap.frameCount();
        resident.cols = cols;
        resident.frameVarying = true;
        resident.atomsPerFrame = atomCount;
        resident.frameCount = run.frameMap.frameCount();
        resident.dtype_descr = QStringLiteral("catalog_row_columns_dft_frame_series");
        resident.values.assign(resident.rows * resident.cols, kNaN);
        for (std::size_t row : run.frameMap.dftRows()) {
            io::QtNpyReader::WidenedArray arr;
            if (!readTrajectoryArray(run, row, col.kind, &arr)) continue;
            if (arr.rows != atomCount || arr.cols != cols) {
                if (err_out) {
                    *err_out = QStringLiteral("%1 shape drift in emitted DFT frame").arg(fieldStem(spec));
                }
                return false;
            }
            const std::size_t dst = row * atomCount * cols;
            std::copy(arr.data.begin(), arr.data.end(), resident.values.begin() + dst);
        }
        run.producerArrays[static_cast<int>(col.kind)] = std::move(resident);
    }
    return true;
}

std::optional<double> CatalogRowValue(const RunData& run,
                                      io::FieldKind kind,
                                      int component,
                                      std::size_t atom,
                                      std::size_t frame,
                                      std::size_t) {
    const StaticNpyArray* a = run.producerArray(kind);
    if (!a) return std::nullopt;
    if (component < 0 || static_cast<std::size_t>(component) >= a->cols) return std::nullopt;
    const std::size_t sourceRow = a->frameVarying ? a->rowFor(atom, frame) : atom;
    if (sourceRow >= a->rows) return std::nullopt;
    const double v = a->value(sourceRow, static_cast<std::size_t>(component));
    if (!std::isfinite(v)) return std::nullopt;
    return v;
}

bool WriteCatalogCoverageArtifacts(const QString& outDir,
                                   const std::vector<RunData>& runs,
                                   const std::vector<RowColumnSpec>&,
                                   const RowDesignStats& stats,
                                   CatalogCoverageArtifacts* artifacts,
                                   QString* err_out) {
    QDir().mkpath(outDir);
    QJsonArray fields;
    QStringList sidecars;

    for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
        QJsonObject f = baseFieldJson(*spec);
        if (scalarAtomRowField(*spec)) {
            const std::size_t n = nominalComponentCount(*spec);
            QJsonArray cols;
            QJsonArray counts;
            for (std::size_t c = 0; c < n; ++c) {
                const QString name = rowColumnName(*spec, c);
                cols.push_back(name);
                const int idx = RowDesignColumnIndex(name);
                counts.push_back(idx >= 0 && static_cast<std::size_t>(idx) < stats.populatedCounts.size()
                                     ? static_cast<double>(stats.populatedCounts[static_cast<std::size_t>(idx)])
                                     : 0.0);
            }
            f.insert(QStringLiteral("representation"), QStringLiteral("row_columns"));
            f.insert(QStringLiteral("row_columns"), cols);
            f.insert(QStringLiteral("populated_counts"), counts);
            f.insert(QStringLiteral("row_count"), static_cast<double>(stats.rows));
        } else if (structuredField(*spec)) {
            if (!writeStructuredIndex(outDir, *spec, runs, &f, err_out)) return false;
            sidecars << f.value(QStringLiteral("sidecar")).toString();
        } else if (spec->axis == io::NativeAxis::Atom) {
            if (!writeAtomSidecar(outDir, *spec, runs, &f, err_out)) return false;
            if (f.value(QStringLiteral("sidecar")).isString()) sidecars << f.value(QStringLiteral("sidecar")).toString();
        } else {
            if (!writeNativeNumericSidecar(outDir, *spec, runs, &f, err_out)) return false;
            if (f.value(QStringLiteral("sidecar")).isString()) sidecars << f.value(QStringLiteral("sidecar")).toString();
        }
        fields.push_back(f);
    }

    QJsonObject doc;
    doc.insert(QStringLiteral("schema_version"), 1);
    doc.insert(QStringLiteral("scope"), QStringLiteral("producer catalog minus Larsen, tripeptide, delta, rediscover, legacy groups"));
    doc.insert(QStringLiteral("field_count"), fields.size());
    doc.insert(QStringLiteral("fields"), fields);
    if (!writeText(QDir(outDir).filePath(QStringLiteral("catalog_coverage.json")),
                   QJsonDocument(doc).toJson(QJsonDocument::Indented), err_out)) {
        return false;
    }

    QFile csv(QDir(outDir).filePath(QStringLiteral("catalog_sidecar_support.csv")));
    if (!csv.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write catalog_sidecar_support.csv");
        return false;
    }
    QTextStream ts(&csv);
    ts << "stem,representation,component,populated_count,row_count,sidecar\n";
    for (const QJsonValue& v : fields) {
        const QJsonObject o = v.toObject();
        const QJsonArray counts = o.value(QStringLiteral("populated_counts")).toArray();
        for (int i = 0; i < counts.size(); ++i) {
            ts << o.value(QStringLiteral("stem")).toString() << ","
               << o.value(QStringLiteral("representation")).toString() << ","
               << i << ","
               << QString::number(counts[i].toDouble(), 'f', 0) << ","
               << QString::number(o.value(QStringLiteral("row_count")).toDouble(), 'f', 0) << ","
               << o.value(QStringLiteral("sidecar")).toString() << "\n";
        }
        const QJsonObject reduction = o.value(QStringLiteral("row_reduction")).toObject();
        if (!reduction.isEmpty()) {
            const QJsonArray rCounts = reduction.value(QStringLiteral("populated_counts")).toArray();
            const QJsonArray names = reduction.value(QStringLiteral("component_names")).toArray();
            for (int i = 0; i < rCounts.size(); ++i) {
                const QString component = i < names.size()
                                              ? QStringLiteral("reduction:%1").arg(names[i].toString())
                                              : QStringLiteral("reduction:%1").arg(i);
                ts << o.value(QStringLiteral("stem")).toString() << ","
                   << reduction.value(QStringLiteral("representation")).toString() << ","
                   << component << ","
                   << QString::number(rCounts[i].toDouble(), 'f', 0) << ","
                   << QString::number(reduction.value(QStringLiteral("row_count")).toDouble(), 'f', 0) << ","
                   << reduction.value(QStringLiteral("sidecar")).toString() << "\n";
            }
        }
    }

    if (artifacts) {
        artifacts->fields = fields;
        artifacts->sidecarFiles = sidecars;
    }
    return true;
}

}  // namespace h5reader::rediscover
