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

QByteArray npyHeader(const QByteArray& descr, const std::vector<std::size_t>& shape) {
    QByteArray header;
    header += "{'descr': '";
    header += descr;
    header += "', 'fortran_order': False, 'shape': (";
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

bool writeStructuredIndex(const QString& outDir,
                          const io::FieldSpec& spec,
                          const std::vector<RunData>& runs,
                          QJsonObject* field,
                          QString* err_out) {
    const QString indexName = QStringLiteral("row_design_field_%1_index.csv").arg(fieldStem(spec));
    QFile indexFile(QDir(outDir).filePath(indexName));
    if (!indexFile.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write %1").arg(indexName);
        return false;
    }
    QTextStream ts(&indexFile);
    ts << "protein_ordinal,protein_id,source_path,row_start,row_count\n";

    std::size_t total = 0;
    int proteinOrdinal = 0;
    for (const RunData& run : runs) {
        const QString path = staticFieldPath(run, spec.kind);
        std::vector<unsigned char> bytes;
        std::size_t recordSize = 0;
        std::size_t rows = 0;
        if (QFileInfo::exists(path)) {
            auto r = io::QtNpyReader::ReadRawBytes(path, bytes, recordSize);
            if (r.ok) rows = r.row_count;
        }
        ts << proteinOrdinal << "," << (run.protein ? run.protein->proteinId() : QString())
           << "," << path << "," << total << "," << rows << "\n";
        total += rows;
        ++proteinOrdinal;
    }
    (*field).insert(QStringLiteral("representation"), QStringLiteral("structured_native_sidecar_index"));
    (*field).insert(QStringLiteral("sidecar"), indexName);
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
            if (!readRunArray(run, row, spec.kind, &arr)) continue;
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
            if (!readRunArray(run, row, spec.kind, &arr)) continue;
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
    }

    if (artifacts) {
        artifacts->fields = fields;
        artifacts->sidecarFiles = sidecars;
    }
    return true;
}

}  // namespace h5reader::rediscover
