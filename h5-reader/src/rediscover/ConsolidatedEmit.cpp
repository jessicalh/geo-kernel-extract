#include "ConsolidatedEmit.h"

#include "AnalysisBody.h"
#include "CanonicalSpineGuard.h"
#include "Catalog.h"
#include "ExtractionSupport.h"
#include "ResidentIndexes.h"
#include "RowDesign.h"
#include "RunData.h"
#include "ScopedProducerCatalog.h"
#include "SphericalBasis.h"
#include "StaticRunData.h"

#include "../io/QtNpyReader.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtResidueNames.h"
#include "../model/QtRing.h"
#include "../model/QtRingMembership.h"
#include "../model/SphericalDecomposition.h"

#include <QCryptographicHash>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QTextStream>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace h5reader::rediscover {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

QString qsv(std::string_view s) {
    return QString::fromUtf8(s.data(), static_cast<qsizetype>(s.size()));
}

QString fieldStem(const io::FieldSpec& spec) { return qsv(spec.stem); }

QString num(double v) {
    return std::isfinite(v) ? QString::number(v, 'g', 15) : QString();
}

QString csv(QString s) {
    if (s.contains(QLatin1Char('"'))) s.replace(QStringLiteral("\""), QStringLiteral("\"\""));
    if (s.contains(QLatin1Char(',')) || s.contains(QLatin1Char('\n'))
        || s.contains(QLatin1Char('"'))) {
        return QStringLiteral("\"%1\"").arg(s);
    }
    return s;
}

QString boolText(bool v) { return v ? QStringLiteral("1") : QStringLiteral("0"); }

QString axisName(io::NativeAxis axis) {
    switch (axis) {
    case io::NativeAxis::Atom: return QStringLiteral("Atom");
    case io::NativeAxis::RingContributionPair: return QStringLiteral("RingContributionPair");
    case io::NativeAxis::AromaticRing: return QStringLiteral("AromaticRing");
    case io::NativeAxis::RediscoverSourceRow: return QStringLiteral("RediscoverSourceRow");
    case io::NativeAxis::RediscoverAggregatedRow: return QStringLiteral("RediscoverAggregatedRow");
    case io::NativeAxis::RediscoverTargetRow: return QStringLiteral("RediscoverTargetRow");
    case io::NativeAxis::Protein: return QStringLiteral("Protein");
    case io::NativeAxis::Bond: return QStringLiteral("Bond");
    case io::NativeAxis::MOPACBondNeighborPair: return QStringLiteral("MOPACBondNeighborPair");
    case io::NativeAxis::MOPACUniquePair: return QStringLiteral("MOPACUniquePair");
    case io::NativeAxis::MutationMatchPair: return QStringLiteral("MutationMatchPair");
    case io::NativeAxis::Residue: return QStringLiteral("Residue");
    case io::NativeAxis::SaturatedRing: return QStringLiteral("SaturatedRing");
    case io::NativeAxis::Ring: return QStringLiteral("Ring");
    case io::NativeAxis::RingMembership: return QStringLiteral("RingMembership");
    }
    return QStringLiteral("Unknown");
}

QString ss8Label(SecondaryStructure8 ss8) {
    switch (ss8) {
    case SecondaryStructure8::H: return QStringLiteral("H");
    case SecondaryStructure8::G: return QStringLiteral("G");
    case SecondaryStructure8::I: return QStringLiteral("I");
    case SecondaryStructure8::E: return QStringLiteral("E");
    case SecondaryStructure8::B: return QStringLiteral("B");
    case SecondaryStructure8::T: return QStringLiteral("T");
    case SecondaryStructure8::S: return QStringLiteral("S");
    case SecondaryStructure8::C: return QStringLiteral("C");
    case SecondaryStructure8::Unknown: return QString();
    }
    return QString();
}

QString ss3Label(SecondaryStructure3 ss3) {
    switch (ss3) {
    case SecondaryStructure3::Helix: return QStringLiteral("helix");
    case SecondaryStructure3::Sheet: return QStringLiteral("sheet");
    case SecondaryStructure3::Coil: return QStringLiteral("coil");
    case SecondaryStructure3::Unknown: return QString();
    }
    return QString();
}

QString residueClassName(ResidueClass c) {
    return QString::fromLatin1(NameForResidueClass(c));
}

QJsonArray stringArray(const std::vector<QString>& values) {
    QJsonArray out;
    for (const QString& v : values) out.append(v);
    return out;
}

std::size_t product(const std::vector<std::size_t>& shape) {
    std::size_t out = 1;
    for (std::size_t dim : shape) out *= dim;
    return out;
}

QByteArray npyHeader(const QByteArray& dtype, const std::vector<std::size_t>& shape) {
    QByteArray header;
    header += "{'descr': '";
    header += dtype;
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
    return header;
}

bool writeNpyPrefix(QSaveFile& f, const QByteArray& dtype,
                    const std::vector<std::size_t>& shape) {
    const QByteArray header = npyHeader(dtype, shape);
    if (header.size() > 65535) return false;
    QByteArray prefix;
    prefix.append("\x93NUMPY", 6);
    prefix.append(char(1));
    prefix.append(char(0));
    const quint16 headerLen = static_cast<quint16>(header.size());
    prefix.append(char(headerLen & 0xff));
    prefix.append(char((headerLen >> 8) & 0xff));
    return f.write(prefix) == prefix.size() && f.write(header) == header.size();
}

class RawArrayWriter {
public:
    RawArrayWriter() = default;
    RawArrayWriter(QString finalPath, QByteArray dtype, std::size_t elemSize)
        : finalPath_(std::move(finalPath)),
          tempPath_(finalPath_ + QStringLiteral(".rawtmp")),
          dtype_(std::move(dtype)),
          elemSize_(elemSize),
          temp_(tempPath_) {}

    bool open(QString* err_out) {
        QDir().mkpath(QFileInfo(finalPath_).dir().absolutePath());
        if (!temp_.open(QIODevice::WriteOnly)) {
            if (err_out) *err_out = QStringLiteral("cannot open raw temp: %1").arg(tempPath_);
            return false;
        }
        return true;
    }

    template <typename T>
    bool writeOne(const T& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        if (sizeof(T) != elemSize_) return false;
        if (temp_.write(reinterpret_cast<const char*>(&value), sizeof(T))
            != static_cast<qint64>(sizeof(T))) {
            return false;
        }
        ++count_;
        return true;
    }

    template <typename T>
    bool writeMany(const T* values, std::size_t n) {
        static_assert(std::is_trivially_copyable_v<T>);
        if (sizeof(T) != elemSize_) return false;
        const qint64 bytes = static_cast<qint64>(n * sizeof(T));
        if (bytes > 0
            && temp_.write(reinterpret_cast<const char*>(values), bytes) != bytes) {
            return false;
        }
        count_ += n;
        return true;
    }

    bool commit(const std::vector<std::size_t>& shape, QString* err_out) {
        temp_.flush();
        temp_.close();
        if (product(shape) != count_) {
            if (err_out) {
                *err_out = QStringLiteral("shape/count mismatch for %1: shape product %2, wrote %3")
                               .arg(finalPath_)
                               .arg(product(shape))
                               .arg(count_);
            }
            return false;
        }

        QFile in(tempPath_);
        if (!in.open(QIODevice::ReadOnly)) {
            if (err_out) *err_out = QStringLiteral("cannot reopen raw temp: %1").arg(tempPath_);
            return false;
        }
        QSaveFile out(finalPath_);
        if (!out.open(QIODevice::WriteOnly)) {
            if (err_out) *err_out = QStringLiteral("cannot open npy: %1").arg(finalPath_);
            return false;
        }
        if (!writeNpyPrefix(out, dtype_, shape)) {
            if (err_out) *err_out = QStringLiteral("cannot write npy header: %1").arg(finalPath_);
            return false;
        }
        QByteArray buf(1 << 20, Qt::Uninitialized);
        while (!in.atEnd()) {
            const qint64 n = in.read(buf.data(), buf.size());
            if (n < 0 || out.write(buf.constData(), n) != n) {
                if (err_out) *err_out = QStringLiteral("cannot copy raw payload: %1").arg(finalPath_);
                return false;
            }
        }
        in.close();
        if (!out.commit()) {
            if (err_out) *err_out = QStringLiteral("cannot commit npy: %1").arg(finalPath_);
            return false;
        }
        QFile::remove(tempPath_);
        return true;
    }

    std::size_t count() const { return count_; }
    const QString& path() const { return finalPath_; }

private:
    QString finalPath_;
    QString tempPath_;
    QByteArray dtype_;
    std::size_t elemSize_ = 0;
    QFile temp_;
    std::size_t count_ = 0;
};

struct ColumnSupport {
    QString name;
    QString kind;
    QString source;
    bool required = false;
    std::size_t rowCount = 0;
    std::size_t populated = 0;
    std::map<QString, std::size_t> absence;
};

struct DatasetSupport {
    std::size_t rows = 0;
    std::size_t populated = 0;
    std::map<QString, std::size_t> absence;
};

struct FieldSupport {
    const io::FieldSpec* spec = nullptr;
    QString route;
    QString sidecar;
    QString shape;
    QString dtype;
    QString componentLayout;
    std::size_t count = 0;
    std::size_t populated = 0;
    std::vector<std::size_t> componentPopulated;
    std::map<QString, DatasetSupport> byDataset;
    std::map<QString, std::size_t> absence;
};

void addAbsence(std::map<QString, std::size_t>& m, const QString& reason, std::size_t n = 1) {
    m[reason.isEmpty() ? QStringLiteral("unknown") : reason] += n;
}

QString absenceSummary(const std::map<QString, std::size_t>& m) {
    QStringList parts;
    for (const auto& item : m)
        parts << QStringLiteral("%1:%2").arg(item.first).arg(item.second);
    return parts.join(QLatin1Char(';'));
}

std::size_t componentCountFor(const io::FieldSpec& spec) {
    return spec.cols > 0 ? static_cast<std::size_t>(spec.cols) : 1;
}

bool isStructuredField(const io::FieldSpec& spec) {
    return spec.kind == io::FieldKind::AtomsCategoryInfo
           || spec.group == io::FieldGroup::Topology;
}

bool isTargetField(const io::FieldSpec& spec) {
    return spec.kind == io::FieldKind::OrcaTotal
           || spec.kind == io::FieldKind::OrcaDiamagnetic
           || spec.kind == io::FieldKind::OrcaParamagnetic;
}

bool isRingContributionPair(const io::FieldSpec& spec) {
    return spec.axis == io::NativeAxis::RingContributionPair;
}

bool isIncidenceAxis(const io::FieldSpec& spec) {
    return spec.axis == io::NativeAxis::Bond
           || spec.axis == io::NativeAxis::MOPACBondNeighborPair
           || spec.axis == io::NativeAxis::MOPACUniquePair;
}

QString componentLayoutFor(const io::FieldSpec& spec) {
    if (spec.cols == 1 || spec.cols < 0) return QStringLiteral("scalar");
    if (spec.cols == 3) return QStringLiteral("vec3");
    if (spec.cols == 5) return QStringLiteral("T2_5");
    if (spec.cols == 9) return QStringLiteral("tensor9_raw");
    return QStringLiteral("cols:%1").arg(spec.cols);
}

enum class FieldRoute {
    Measure,
    NativePair,
    Incidence,
    Target,
    Topology,
};

struct RowSchema {
    QStringList columns;
    std::map<QString, std::size_t> index;
};

struct RowBuilder {
    const RowSchema& schema;
    std::vector<QString> values;

    explicit RowBuilder(const RowSchema& s) : schema(s), values(s.columns.size()) {}

    void set(const QString& name, const QString& value) {
        auto it = schema.index.find(name);
        if (it != schema.index.end()) values[it->second] = value;
    }
    void setInt(const QString& name, qint64 value) { set(name, QString::number(value)); }
    void setBool(const QString& name, bool value) { set(name, boolText(value)); }
    void setDouble(const QString& name, double value) { set(name, num(value)); }
};

struct RunContext {
    QString datasetId;
    QString proteinId;
    QString poseId;
    QString conformationId;
    QString poseKind;
    QString splitGroupId;
    std::size_t rowBase = 0;
};

struct RegionBasin {
    QString id;
    QString label;
    double phi = kNaN;
    double psi = kNaN;
};

struct RegionApplySpec {
    bool supplied = false;
    QString regionDefId;
    QString sourcePath;
    QString sha256;
    std::vector<RegionBasin> ramaBasins;
};

struct FieldSink {
    const io::FieldSpec* spec = nullptr;
    FieldRoute route = FieldRoute::Measure;
    std::size_t cols = 1;
    std::unique_ptr<RawArrayWriter> values;
    std::unique_ptr<RawArrayWriter> presence;
    std::unique_ptr<RawArrayWriter> byAtom;
    std::unique_ptr<RawArrayWriter> nativeValues;
    std::unique_ptr<QSaveFile> absenceFile;
    std::unique_ptr<QTextStream> absenceOut;
    std::unique_ptr<QSaveFile> nativeMapFile;
    std::unique_ptr<QTextStream> nativeMapOut;
    FieldSupport support;
    std::size_t rowAlignedRows = 0;
    std::size_t nativeRows = 0;
};

struct TargetSinks {
    RawArrayWriter totalT0;
    RawArrayWriter totalT1;
    RawArrayWriter totalT2;
    RawArrayWriter totalRaw;
    RawArrayWriter diaT0;
    RawArrayWriter diaT1;
    RawArrayWriter diaT2;
    RawArrayWriter diaRaw;
    RawArrayWriter paraT0;
    RawArrayWriter paraT1;
    RawArrayWriter paraT2;
    RawArrayWriter paraRaw;
    RawArrayWriter present;

    explicit TargetSinks(const QString& outDir);
    bool open(QString* err_out);
    bool commit(std::size_t nRows, QString* err_out);
};

struct E3nnSinks {
    QSaveFile conformationsFile;
    std::unique_ptr<QTextStream> conformationsOut;
    RawArrayWriter graphOffsets;
    RawArrayWriter atomIndex;
    RawArrayWriter z;
    RawArrayWriter positions;
    RawArrayWriter aimnet2Embedding;
    RawArrayWriter rowToGraphAtom;
    RawArrayWriter sigmaTotalRaw;
    RawArrayWriter sigmaDiaRaw;
    RawArrayWriter sigmaParaRaw;
    std::size_t graphCount = 0;
    std::size_t graphAtoms = 0;

    explicit E3nnSinks(const QString& outDir);
    bool open(QString* err_out);
    bool commit(std::size_t nRows, QString* err_out);
};

QString routeName(FieldRoute route);
FieldRoute routeFor(const io::FieldSpec& spec);
RegionApplySpec parseRegionSpec(const ConditioningSpec& spec);
std::optional<RegionBasin> nearestRamaBasin(const RegionApplySpec& spec, double phi, double psi);
std::vector<int> configuredCutoffs(const ConditioningSpec& spec);
QString residueTypeLabel(const model::QtProtein& p, const model::QtResidue* r);
QString residueTypeLabelForAa(model::AminoAcid aa);
const model::QtResidue* residueForAtom(const model::QtProtein& p, const model::QtAtom& a);
std::set<int32_t> bondedAtoms(const model::QtProtein& p, std::size_t atom);
std::size_t countNear(const Body& body, CloudKind kind, std::size_t atom, std::size_t frame,
                      double cutoff, const std::set<int32_t>& bonded);
double nearestDist(const Body& body, CloudKind kind, std::size_t atom, std::size_t frame,
                   const std::set<int32_t>& bonded, int32_t* entityOut);
QString ringMembershipId(const model::QtProtein& p, std::size_t atom);
QString ringTypeForAtom(const model::QtProtein& p, std::size_t atom);
std::optional<std::size_t> nativeRowForAtom(const Body& body, const io::FieldSpec& spec,
                                            std::size_t atom, QString* reason_out);
std::array<double, 9> matArray(const model::Mat3& m);
void writeTensorTarget(RawArrayWriter& raw, RawArrayWriter& t0, RawArrayWriter& t1,
                       RawArrayWriter& t2, const model::Mat3& m, bool present);
std::optional<std::pair<int32_t, int32_t>> nativeEndpoints(const Body& body,
                                                           const io::FieldSpec& spec,
                                                           std::size_t nativeRow,
                                                           std::size_t frame);
void emitMeasureField(FieldSink& sink, const Body& body, const RunContext& ctx,
                      std::size_t rowId, std::size_t atom, std::size_t frame);
std::vector<std::vector<std::size_t>> buildIncidentRows(const Body& body,
                                                        const io::FieldSpec& spec,
                                                        std::size_t frame,
                                                        std::size_t nativeCount);
void emitIncidenceNative(FieldSink& sink, const Body& body, const RunContext& ctx,
                         std::size_t frame, std::size_t nativeCount,
                         const std::vector<std::vector<std::size_t>>& incidentRows);
void emitNativePairField(FieldSink& sink, const Body& body, const RunContext& ctx,
                         std::size_t frame, std::size_t nativeCount);
void emitIncidenceRow(FieldSink& sink, const Body& body, const RunContext& ctx,
                      std::size_t rowId, std::size_t atom, std::size_t frame,
                      const std::vector<std::vector<std::size_t>>& incidentRows);
void writeTargets(TargetSinks& targets, E3nnSinks& e3nn, const RunData& run,
                  std::size_t atom, std::size_t frame);
void writeE3nnAtom(E3nnSinks& e3nn, const Body& body, std::size_t graphId,
                   std::size_t graphAtomOffset, std::size_t rowId,
                   std::size_t atom, std::size_t frame);
RowSchema buildRowSchema(const std::vector<int>& cutoffs);
std::vector<QString> enumerate720Poses(const QString& root720);
std::size_t preflightRowCount(const QString& root720, const std::vector<QString>& poses,
                              const RunData& onep9j, QString* err_out);
std::vector<std::unique_ptr<FieldSink>> openFieldSinks(const QString& outDir,
                                                       std::size_t nRows,
                                                       QString* err_out);
void emitRun(RunData& run, const RunContext& ctx, const QString& outDir,
             const RowSchema& rowSchema, const RegionApplySpec& region,
             const std::vector<int>& cutoffs, QTextStream& rowsOut,
             std::map<QString, ColumnSupport>& columnSupport,
             std::vector<std::unique_ptr<FieldSink>>& fieldSinks,
             TargetSinks& targets, E3nnSinks& e3nn,
             std::map<QString, FieldSupport*>& topologySupport,
             std::size_t* nextRowId);
bool targetPresent(const RunData& run, std::size_t atom, std::size_t frame);
void fillRowLabels(RowBuilder& row, const Body& body, const RunContext& ctx,
                   const RegionApplySpec& region, const std::vector<int>& cutoffs,
                   std::size_t rowId, std::size_t atom, std::size_t frame);
void writeRow(QTextStream& out, const RowSchema& schema, const RowBuilder& row,
              std::map<QString, ColumnSupport>& support);
bool commitFieldSinks(std::vector<std::unique_ptr<FieldSink>>& sinks,
                      std::size_t nRows, QString* err_out);
bool writeColumnSupport(const QString& outDir, const RowSchema& schema,
                        const std::map<QString, ColumnSupport>& support,
                        QString* err_out);
bool writeSidecarSupport(const QString& outDir,
                         const std::vector<std::unique_ptr<FieldSink>>& sinks,
                         std::size_t nRows, QString* err_out);
bool writeManifests(const QString& outDir, const ConsolidatedEmitOptions& options,
                    const RowSchema& rowSchema,
                    const std::vector<std::unique_ptr<FieldSink>>& sinks,
                    const RegionApplySpec& region,
                    const std::vector<int>& cutoffs,
                    std::size_t nRows,
                    QString* err_out);
bool survivalGate(const RowSchema& rowSchema,
                  const std::vector<std::unique_ptr<FieldSink>>& sinks,
                  std::size_t nRows,
                  QString* err_out);

}  // namespace
}  // namespace h5reader::rediscover

namespace h5reader::rediscover {

bool RunConsolidatedEmit(const ConsolidatedEmitOptions& options,
                         ConsolidatedEmitStats* stats,
                         QString* err_out) {
    auto fail = [&](const QString& msg) {
        if (err_out) *err_out = msg;
        return false;
    };

    if (options.root720.isEmpty()) return fail(QStringLiteral("consolidated_emit requires --root720"));
    if (options.run1p9j.isEmpty()) return fail(QStringLiteral("consolidated_emit requires --run"));
    if (options.outDir.isEmpty()) return fail(QStringLiteral("consolidated_emit requires --out"));

    QString guardErr;
    if (!ValidateCanonical720Root(options.root720, &guardErr))
        return fail(QStringLiteral("720 canonical guard failed: %1").arg(guardErr));
    if (!ValidateCanonical1p9jRoot(options.run1p9j, &guardErr))
        return fail(QStringLiteral("1P9J canonical guard failed: %1").arg(guardErr));

    QDir outRoot(options.outDir);
    if (outRoot.exists()) {
        const QStringList entries = outRoot.entryList(QDir::NoDotAndDotDot | QDir::AllEntries);
        if (!entries.isEmpty())
            return fail(QStringLiteral("output directory is not empty: %1").arg(options.outDir));
    } else if (!QDir().mkpath(options.outDir)) {
        return fail(QStringLiteral("cannot create output directory: %1").arg(options.outDir));
    }
    for (const QString& sub : {QStringLiteral("measures"),
                               QStringLiteral("native"),
                               QStringLiteral("targets"),
                               QStringLiteral("topology"),
                               QStringLiteral("e3nn")}) {
        if (!outRoot.mkpath(sub)) return fail(QStringLiteral("cannot create output subdir: %1").arg(sub));
    }

    QString loadErr;
    auto onep9j = RunLoader::Load(options.run1p9j, &loadErr);
    if (!onep9j) return fail(QStringLiteral("1P9J load failed: %1").arg(loadErr));
    const std::vector<QString> poses = enumerate720Poses(options.root720);
    if (poses.empty()) return fail(QStringLiteral("720 root has no pose directories"));

    const std::size_t nRows = preflightRowCount(options.root720, poses, *onep9j, &loadErr);
    if (!loadErr.isEmpty()) return fail(QStringLiteral("row preflight failed: %1").arg(loadErr));
    if (nRows == 0) return fail(QStringLiteral("row preflight found zero rows"));

    RegionApplySpec region;
    std::vector<int> cutoffs;
    try {
        region = parseRegionSpec(options.conditioningSpec);
        cutoffs = configuredCutoffs(options.conditioningSpec);
    } catch (const std::exception& e) {
        return fail(QStringLiteral("conditioning spec failed: %1").arg(e.what()));
    }
    const RowSchema rowSchema = buildRowSchema(cutoffs);

    QSaveFile rowsFile(outRoot.filePath(QStringLiteral("rows.csv")));
    if (!rowsFile.open(QIODevice::WriteOnly | QIODevice::Text))
        return fail(QStringLiteral("cannot open rows.csv"));
    QTextStream rowsOut(&rowsFile);
    rowsOut << rowSchema.columns.join(QLatin1Char(',')) << '\n';

    QString sinkErr;
    std::vector<std::unique_ptr<FieldSink>> fieldSinks = openFieldSinks(options.outDir, nRows, &sinkErr);
    if (!sinkErr.isEmpty()) return fail(QStringLiteral("sidecar open failed: %1").arg(sinkErr));
    if (fieldSinks.size() != ScopedProducerCatalog().size())
        return fail(QStringLiteral("sidecar open produced %1 sinks, expected %2")
                        .arg(fieldSinks.size())
                        .arg(ScopedProducerCatalog().size()));

    TargetSinks targets(options.outDir);
    if (!targets.open(&sinkErr)) return fail(QStringLiteral("target sink open failed: %1").arg(sinkErr));
    E3nnSinks e3nn(options.outDir);
    if (!e3nn.open(&sinkErr)) return fail(QStringLiteral("e3nn sink open failed: %1").arg(sinkErr));

    std::map<QString, FieldSupport*> topologySupport;
    for (const auto& sink : fieldSinks) {
        if (sink->route == FieldRoute::Topology)
            topologySupport[fieldStem(*sink->spec)] = &sink->support;
    }
    std::map<QString, ColumnSupport> columnSupport;
    std::size_t nextRowId = 0;

    try {
        QDir root(options.root720);
        for (const QString& pose : poses) {
            const QString poseDir = root.filePath(pose);
            if (!ValidateCanonical720Pose(poseDir, &guardErr)) {
                return fail(QStringLiteral("720 pose canonical guard failed for %1: %2")
                                .arg(pose, guardErr));
            }
            QString staticErr;
            auto staticRun = StaticRunData::Load(poseDir, &staticErr);
            if (!staticRun) {
                return fail(QStringLiteral("720 static load failed for %1: %2")
                                .arg(pose, staticErr));
            }
            RunContext ctx;
            ctx.datasetId = QStringLiteral("720_static");
            ctx.proteinId = staticRun->protein ? staticRun->protein->proteinId()
                                                : staticRun->manifest.protein_id;
            ctx.poseId = pose;
            ctx.conformationId = pose;
            ctx.poseKind = QStringLiteral("static");
            ctx.splitGroupId = ctx.proteinId;
            ctx.rowBase = nextRowId;
            emitRun(*staticRun, ctx, options.outDir, rowSchema, region, cutoffs, rowsOut,
                    columnSupport, fieldSinks, targets, e3nn, topologySupport, &nextRowId);
        }

        RunContext ctx;
        ctx.datasetId = QStringLiteral("1p9j");
        ctx.proteinId = onep9j->protein ? onep9j->protein->proteinId()
                                         : onep9j->manifest.protein_id;
        ctx.poseId = onep9j->manifest.dataset_id.isEmpty()
                         ? QStringLiteral("md_20260610T103128Z")
                         : onep9j->manifest.dataset_id;
        ctx.conformationId = ctx.poseId;
        ctx.poseKind = QStringLiteral("trajectory");
        ctx.splitGroupId = ctx.proteinId;
        ctx.rowBase = nextRowId;
        emitRun(*onep9j, ctx, options.outDir, rowSchema, region, cutoffs, rowsOut,
                columnSupport, fieldSinks, targets, e3nn, topologySupport, &nextRowId);
    } catch (const std::exception& e) {
        return fail(QStringLiteral("emit failed: %1").arg(e.what()));
    }

    if (nextRowId != nRows)
        return fail(QStringLiteral("row-count mismatch: preflight %1, emitted %2").arg(nRows).arg(nextRowId));

    const std::int64_t finalGraphOffset = static_cast<std::int64_t>(e3nn.graphAtoms);
    e3nn.graphOffsets.writeOne(finalGraphOffset);
    rowsOut.flush();
    if (!rowsFile.commit()) return fail(QStringLiteral("cannot commit rows.csv"));
    if (!targets.commit(nRows, &sinkErr)) return fail(QStringLiteral("target sink commit failed: %1").arg(sinkErr));
    if (!e3nn.commit(nRows, &sinkErr)) return fail(QStringLiteral("e3nn sink commit failed: %1").arg(sinkErr));
    if (!commitFieldSinks(fieldSinks, nRows, &sinkErr))
        return fail(QStringLiteral("field sink commit failed: %1").arg(sinkErr));
    if (!writeColumnSupport(options.outDir, rowSchema, columnSupport, &sinkErr))
        return fail(QStringLiteral("column_support.csv failed: %1").arg(sinkErr));
    if (!writeSidecarSupport(options.outDir, fieldSinks, nRows, &sinkErr))
        return fail(QStringLiteral("sidecar_support.csv failed: %1").arg(sinkErr));
    if (!writeManifests(options.outDir, options, rowSchema, fieldSinks, region, cutoffs, nRows, &sinkErr))
        return fail(QStringLiteral("manifest write failed: %1").arg(sinkErr));
    if (!survivalGate(rowSchema, fieldSinks, nRows, &sinkErr))
        return fail(sinkErr);

    if (stats) {
        stats->rows = nRows;
        stats->scopedFieldCount = ScopedProducerCatalog().size();
    }
    return true;
}

}  // namespace h5reader::rediscover

namespace h5reader::rediscover {
namespace {

bool commitFieldSinks(std::vector<std::unique_ptr<FieldSink>>& sinks,
                      std::size_t nRows,
                      QString* err_out) {
    for (auto& sink : sinks) {
        if (sink->absenceOut) sink->absenceOut->flush();
        if (sink->nativeMapOut) sink->nativeMapOut->flush();
        switch (sink->route) {
        case FieldRoute::Measure:
            if (!sink->values->commit(sink->cols == 1 ? std::vector<std::size_t>{nRows}
                                                       : std::vector<std::size_t>{nRows, sink->cols},
                                      err_out)
                || !sink->presence->commit({nRows}, err_out)) {
                return false;
            }
            if (sink->absenceFile && !sink->absenceFile->commit()) {
                if (err_out) *err_out = QStringLiteral("cannot commit absence CSV for %1").arg(fieldStem(*sink->spec));
                return false;
            }
            break;
        case FieldRoute::Incidence:
            if (!sink->byAtom->commit({nRows, 4, sink->cols}, err_out)
                || !sink->presence->commit({nRows}, err_out)
                || !sink->nativeValues->commit({sink->nativeRows, sink->cols}, err_out)) {
                return false;
            }
            if (sink->absenceFile && !sink->absenceFile->commit()) {
                if (err_out) *err_out = QStringLiteral("cannot commit absence CSV for %1").arg(fieldStem(*sink->spec));
                return false;
            }
            if (sink->nativeMapFile && !sink->nativeMapFile->commit()) {
                if (err_out) *err_out = QStringLiteral("cannot commit native map for %1").arg(fieldStem(*sink->spec));
                return false;
            }
            sink->support.shape = QStringLiteral("[%1,4,%2]").arg(nRows).arg(sink->cols);
            break;
        case FieldRoute::NativePair:
            if (!sink->nativeValues->commit({sink->nativeRows, sink->cols}, err_out)) return false;
            if (sink->nativeMapFile && !sink->nativeMapFile->commit()) {
                if (err_out) *err_out = QStringLiteral("cannot commit native map for %1").arg(fieldStem(*sink->spec));
                return false;
            }
            sink->support.shape = QStringLiteral("[%1,%2]").arg(sink->nativeRows).arg(sink->cols);
            break;
        case FieldRoute::Target:
        case FieldRoute::Topology:
            break;
        }
    }
    return true;
}

bool writeColumnSupport(const QString& outDir,
                        const RowSchema& schema,
                        const std::map<QString, ColumnSupport>& support,
                        QString* err_out) {
    QSaveFile f(QDir(outDir).filePath(QStringLiteral("column_support.csv")));
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (err_out) *err_out = QStringLiteral("cannot open column_support.csv");
        return false;
    }
    QTextStream out(&f);
    out << "name,kind,source_file,component,row_count,populated_count,nonnull_fraction,required,absence_reason_summary\n";
    for (const QString& name : schema.columns) {
        auto it = support.find(name);
        const ColumnSupport empty{name, QStringLiteral("selector_label"), QStringLiteral("ConsolidatedEmit")};
        const ColumnSupport& st = it == support.end() ? empty : it->second;
        const double frac = st.rowCount ? static_cast<double>(st.populated) / static_cast<double>(st.rowCount) : 0.0;
        out << csv(name) << ',' << csv(st.kind) << ',' << csv(st.source) << ",,"
            << st.rowCount << ',' << st.populated << ',' << num(frac) << ','
            << boolText(name == QStringLiteral("row_id")
                        || name == QStringLiteral("global_row_key")
                        || name == QStringLiteral("dataset_id")
                        || name == QStringLiteral("protein_id")
                        || name == QStringLiteral("pose_id")
                        || name == QStringLiteral("frame_slot")
                        || name == QStringLiteral("atom_index"))
            << ',' << csv(absenceSummary(st.absence)) << '\n';
    }
    return f.commit();
}

bool writeSidecarSupport(const QString& outDir,
                         const std::vector<std::unique_ptr<FieldSink>>& sinks,
                         std::size_t nRows,
                         QString* err_out) {
    QSaveFile f(QDir(outDir).filePath(QStringLiteral("sidecar_support.csv")));
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (err_out) *err_out = QStringLiteral("cannot open sidecar_support.csv");
        return false;
    }
    QTextStream out(&f);
    out << "sidecar,stem,axis,shape,dtype,component_layout,row_count_or_native_count,populated_count,required,scope\n";
    for (const auto& sink : sinks) {
        const io::FieldSpec& spec = *sink->spec;
        std::size_t rowOrNative = 0;
        if (sink->route == FieldRoute::Measure || sink->route == FieldRoute::Incidence)
            rowOrNative = nRows;
        else if (sink->route == FieldRoute::NativePair)
            rowOrNative = sink->nativeRows;
        else if (sink->route == FieldRoute::Target)
            rowOrNative = nRows;
        else
            rowOrNative = sink->support.count;
        out << csv(sink->support.sidecar) << ','
            << csv(fieldStem(spec)) << ','
            << csv(axisName(spec.axis)) << ','
            << csv(sink->support.shape) << ','
            << csv(sink->support.dtype) << ','
            << csv(sink->support.componentLayout) << ','
            << rowOrNative << ','
            << sink->support.populated << ','
            << boolText(spec.required) << ','
            << csv(sink->support.route) << '\n';

        for (const auto& ds : sink->support.byDataset) {
            out << csv(QStringLiteral("%1#dataset:%2").arg(sink->support.sidecar, ds.first)) << ','
                << csv(fieldStem(spec)) << ','
                << csv(axisName(spec.axis)) << ','
                << csv(sink->support.shape) << ','
                << csv(sink->support.dtype) << ','
                << csv(sink->support.componentLayout) << ','
                << ds.second.rows << ','
                << ds.second.populated << ','
                << boolText(spec.required) << ','
                << csv(absenceSummary(ds.second.absence)) << '\n';
        }
    }
    return f.commit();
}

QJsonObject fieldJson(const io::FieldSpec& spec, const FieldSink& sink) {
    QJsonObject o;
    o.insert(QStringLiteral("stem"), fieldStem(spec));
    o.insert(QStringLiteral("axis"), axisName(spec.axis));
    o.insert(QStringLiteral("route"), routeName(sink.route));
    o.insert(QStringLiteral("sidecar"), sink.support.sidecar);
    o.insert(QStringLiteral("shape"), sink.support.shape);
    o.insert(QStringLiteral("dtype"), sink.support.dtype);
    o.insert(QStringLiteral("component_layout"), sink.support.componentLayout);
    o.insert(QStringLiteral("required"), spec.required);
    o.insert(QStringLiteral("is_feature"), spec.is_feature);
    o.insert(QStringLiteral("units"), qsv(spec.units));
    o.insert(QStringLiteral("mechanism"), qsv(spec.mechanism));
    return o;
}

bool writeJsonFile(const QString& path, const QJsonObject& object, QString* err_out) {
    QSaveFile f(path);
    QDir().mkpath(QFileInfo(path).dir().absolutePath());
    if (!f.open(QIODevice::WriteOnly)) {
        if (err_out) *err_out = QStringLiteral("cannot open %1").arg(path);
        return false;
    }
    f.write(QJsonDocument(object).toJson(QJsonDocument::Indented));
    if (!f.commit()) {
        if (err_out) *err_out = QStringLiteral("cannot commit %1").arg(path);
        return false;
    }
    return true;
}

std::vector<QString> aggregateOfFor(const QString& stem) {
    static const std::map<QString, std::vector<QString>> k = {
        {QStringLiteral("bs_shielding"), {QStringLiteral("bs_per_type_T0"), QStringLiteral("bs_per_type_T1"), QStringLiteral("bs_per_type_T2")}},
        {QStringLiteral("hm_shielding"), {QStringLiteral("hm_per_type_T0"), QStringLiteral("hm_per_type_T1"), QStringLiteral("hm_per_type_T2")}},
        {QStringLiteral("larsen_hbond_shielding"), {QStringLiteral("larsen_hbond_1pHB_shielding"), QStringLiteral("larsen_hbond_2pHB_shielding"), QStringLiteral("larsen_hbond_1pHaB_shielding"), QStringLiteral("larsen_hbond_2pHaB_shielding")}},
        {QStringLiteral("coulomb_efg"), {QStringLiteral("coulomb_efg_backbone"), QStringLiteral("coulomb_efg_sidechain"), QStringLiteral("coulomb_efg_aromatic")}},
        {QStringLiteral("coulomb_E"), {QStringLiteral("coulomb_E_backbone"), QStringLiteral("coulomb_E_sidechain"), QStringLiteral("coulomb_E_aromatic")}},
        {QStringLiteral("mopac_coulomb_efg"), {QStringLiteral("mopac_coulomb_efg_backbone"), QStringLiteral("mopac_coulomb_efg_sidechain"), QStringLiteral("mopac_coulomb_efg_aromatic")}},
        {QStringLiteral("mopac_coulomb_E"), {QStringLiteral("mopac_coulomb_E_backbone"), QStringLiteral("mopac_coulomb_E_sidechain"), QStringLiteral("mopac_coulomb_E_aromatic")}},
        {QStringLiteral("orca_total"), {QStringLiteral("orca_diamagnetic"), QStringLiteral("orca_paramagnetic")}},
    };
    auto it = k.find(stem);
    return it == k.end() ? std::vector<QString>{} : it->second;
}

std::pair<QString, QString> familyAndBlockFor(const QString& stem) {
    static const std::map<QString, std::pair<QString, QString>> k = {
        {QStringLiteral("bs_shielding"), {QStringLiteral("ring_current_kernel_collinear"), QStringLiteral("ring_current_kernel_collinear")}},
        {QStringLiteral("hm_shielding"), {QStringLiteral("ring_current_kernel_collinear"), QStringLiteral("ring_current_kernel_collinear")}},
        {QStringLiteral("bs_per_type_T0"), {QStringLiteral("ring_current_bs_total_parts"), QStringLiteral("ring_current_bs_total_or_parts")}},
        {QStringLiteral("bs_per_type_T1"), {QStringLiteral("ring_current_bs_total_parts"), QStringLiteral("ring_current_bs_total_or_parts")}},
        {QStringLiteral("bs_per_type_T2"), {QStringLiteral("ring_current_bs_total_parts"), QStringLiteral("ring_current_bs_total_or_parts")}},
        {QStringLiteral("hm_per_type_T0"), {QStringLiteral("ring_current_hm_total_parts"), QStringLiteral("ring_current_hm_total_or_parts")}},
        {QStringLiteral("hm_per_type_T1"), {QStringLiteral("ring_current_hm_total_parts"), QStringLiteral("ring_current_hm_total_or_parts")}},
        {QStringLiteral("hm_per_type_T2"), {QStringLiteral("ring_current_hm_total_parts"), QStringLiteral("ring_current_hm_total_or_parts")}},
        {QStringLiteral("coulomb_efg"), {QStringLiteral("ff_coulomb_efg_total_parts"), QStringLiteral("ff_coulomb_efg_total_or_parts")}},
        {QStringLiteral("coulomb_efg_backbone"), {QStringLiteral("ff_coulomb_efg_total_parts"), QStringLiteral("ff_coulomb_efg_total_or_parts")}},
        {QStringLiteral("coulomb_efg_sidechain"), {QStringLiteral("ff_coulomb_efg_total_parts"), QStringLiteral("ff_coulomb_efg_total_or_parts")}},
        {QStringLiteral("coulomb_efg_aromatic"), {QStringLiteral("ff_coulomb_efg_total_parts"), QStringLiteral("ff_coulomb_efg_total_or_parts")}},
        {QStringLiteral("coulomb_E"), {QStringLiteral("ff_coulomb_E_total_parts"), QStringLiteral("ff_coulomb_E_total_or_parts")}},
        {QStringLiteral("coulomb_E_backbone"), {QStringLiteral("ff_coulomb_E_total_parts"), QStringLiteral("ff_coulomb_E_total_or_parts")}},
        {QStringLiteral("coulomb_E_sidechain"), {QStringLiteral("ff_coulomb_E_total_parts"), QStringLiteral("ff_coulomb_E_total_or_parts")}},
        {QStringLiteral("coulomb_E_aromatic"), {QStringLiteral("ff_coulomb_E_total_parts"), QStringLiteral("ff_coulomb_E_total_or_parts")}},
        {QStringLiteral("mopac_coulomb_efg"), {QStringLiteral("mopac_coulomb_efg_total_parts"), QStringLiteral("mopac_coulomb_efg_total_or_parts")}},
        {QStringLiteral("mopac_coulomb_efg_backbone"), {QStringLiteral("mopac_coulomb_efg_total_parts"), QStringLiteral("mopac_coulomb_efg_total_or_parts")}},
        {QStringLiteral("mopac_coulomb_efg_sidechain"), {QStringLiteral("mopac_coulomb_efg_total_parts"), QStringLiteral("mopac_coulomb_efg_total_or_parts")}},
        {QStringLiteral("mopac_coulomb_efg_aromatic"), {QStringLiteral("mopac_coulomb_efg_total_parts"), QStringLiteral("mopac_coulomb_efg_total_or_parts")}},
        {QStringLiteral("mopac_coulomb_E"), {QStringLiteral("mopac_coulomb_E_total_parts"), QStringLiteral("mopac_coulomb_E_total_or_parts")}},
        {QStringLiteral("mopac_coulomb_E_backbone"), {QStringLiteral("mopac_coulomb_E_total_parts"), QStringLiteral("mopac_coulomb_E_total_or_parts")}},
        {QStringLiteral("mopac_coulomb_E_sidechain"), {QStringLiteral("mopac_coulomb_E_total_parts"), QStringLiteral("mopac_coulomb_E_total_or_parts")}},
        {QStringLiteral("mopac_coulomb_E_aromatic"), {QStringLiteral("mopac_coulomb_E_total_parts"), QStringLiteral("mopac_coulomb_E_total_or_parts")}},
        {QStringLiteral("orca_total"), {QStringLiteral("orca_total_parts"), QStringLiteral("orca_total_or_parts")}},
        {QStringLiteral("orca_diamagnetic"), {QStringLiteral("orca_total_parts"), QStringLiteral("orca_total_or_parts")}},
        {QStringLiteral("orca_paramagnetic"), {QStringLiteral("orca_total_parts"), QStringLiteral("orca_total_or_parts")}},
    };
    auto it = k.find(stem);
    return it == k.end() ? std::make_pair(QString(), QString()) : it->second;
}

bool writeContributorFamilyManifest(const QString& outDir, QString* err_out) {
    QSaveFile f(QDir(outDir).filePath(QStringLiteral("contributor_family_manifest.json")));
    if (!f.open(QIODevice::WriteOnly)) {
        if (err_out) *err_out = QStringLiteral("cannot open contributor_family_manifest.json");
        return false;
    }
    QJsonArray rows;
    for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
        const QString stem = fieldStem(*spec);
        const auto fb = familyAndBlockFor(stem);
        QJsonObject r;
        r.insert(QStringLiteral("contributor_stem"), stem);
        r.insert(QStringLiteral("family_id"), fb.first);
        r.insert(QStringLiteral("aggregate_of"), stringArray(aggregateOfFor(stem)));
        r.insert(QStringLiteral("mutually_exclusive_with"), QJsonArray{});
        r.insert(QStringLiteral("legal_permutation_block"), fb.second);
        rows.append(r);
    }
    QJsonObject root;
    root.insert(QStringLiteral("schema"),
                QStringLiteral("contributor_stem,family_id,aggregate_of,mutually_exclusive_with,legal_permutation_block"));
    root.insert(QStringLiteral("source"),
                QStringLiteral("/shared/2026Thesis/nmr-shielding/h5-reader/notes/CONTRIBUTOR_FAMILY_MANIFEST.md"));
    root.insert(QStringLiteral("row_count"), static_cast<int>(rows.size()));
    root.insert(QStringLiteral("rows"), rows);
    f.write(QJsonDocument(root).toJson(QJsonDocument::Indented));
    return f.commit();
}

bool writeManifests(const QString& outDir,
                    const ConsolidatedEmitOptions& options,
                    const RowSchema& rowSchema,
                    const std::vector<std::unique_ptr<FieldSink>>& sinks,
                    const RegionApplySpec& region,
                    const std::vector<int>& cutoffs,
                    std::size_t nRows,
                    QString* err_out) {
    QJsonArray columns;
    for (const QString& name : rowSchema.columns) {
        QJsonObject c;
        c.insert(QStringLiteral("name"), name);
        c.insert(QStringLiteral("source"), QStringLiteral("Selector::label"));
        columns.append(c);
    }
    QJsonObject schema;
    schema.insert(QStringLiteral("row_count"), static_cast<int>(nRows));
    schema.insert(QStringLiteral("columns"), columns);
    if (!writeJsonFile(QDir(outDir).filePath(QStringLiteral("schema_manifest.json")), schema, err_out))
        return false;

    QJsonArray sidecars;
    QJsonArray irreps;
    QJsonArray nativeAxes;
    for (const auto& sink : sinks) {
        sidecars.append(fieldJson(*sink->spec, *sink));
        QJsonObject ir;
        ir.insert(QStringLiteral("stem"), fieldStem(*sink->spec));
        ir.insert(QStringLiteral("irreps"), qsv(sink->spec->irreps));
        ir.insert(QStringLiteral("tensor_rank"), sink->spec->tensor_rank);
        ir.insert(QStringLiteral("parity_odd"), sink->spec->parity_odd);
        ir.insert(QStringLiteral("units"), qsv(sink->spec->units));
        ir.insert(QStringLiteral("mechanism"), qsv(sink->spec->mechanism));
        ir.insert(QStringLiteral("tensor_frame"), QStringLiteral("lab"));
        irreps.append(ir);
        QJsonObject na;
        na.insert(QStringLiteral("stem"), fieldStem(*sink->spec));
        na.insert(QStringLiteral("native_axis"), axisName(sink->spec->axis));
        na.insert(QStringLiteral("route"), routeName(sink->route));
        na.insert(QStringLiteral("row_aligned_sidecar"),
                  sink->route == FieldRoute::Measure || sink->route == FieldRoute::Incidence
                      ? sink->support.sidecar
                      : QString());
        na.insert(QStringLiteral("native_sidecar"),
                  sink->route == FieldRoute::NativePair || sink->route == FieldRoute::Incidence
                      ? QStringLiteral("native/%1.npy").arg(fieldStem(*sink->spec))
                      : QString());
        na.insert(QStringLiteral("map_file"),
                  sink->route == FieldRoute::NativePair || sink->route == FieldRoute::Incidence
                      ? QStringLiteral("native/%1_map.csv").arg(fieldStem(*sink->spec))
                      : QString());
        na.insert(QStringLiteral("reducer_order"),
                  sink->route == FieldRoute::Incidence
                      ? QStringLiteral("sum,mean,max_abs,count_nonnull")
                      : QString());
        nativeAxes.append(na);
    }
    QJsonObject sidecarManifest;
    sidecarManifest.insert(QStringLiteral("row_join"),
                           QStringLiteral("row_id append order; no arithmetic row-id contract"));
    sidecarManifest.insert(QStringLiteral("sidecars"), sidecars);
    if (!writeJsonFile(QDir(outDir).filePath(QStringLiteral("sidecar_manifest.json")),
                       sidecarManifest, err_out)) return false;
    QJsonObject irrep;
    irrep.insert(QStringLiteral("fields"), irreps);
    if (!writeJsonFile(QDir(outDir).filePath(QStringLiteral("irrep_manifest.json")), irrep, err_out))
        return false;
    QJsonObject native;
    native.insert(QStringLiteral("fields"), nativeAxes);
    native.insert(QStringLiteral("ring_contribution_pair_authority"),
                  QStringLiteral("producer ring_contributions columns 0-2 checked against every pair-axis row count and atom/ring bounds"));
    if (!writeJsonFile(QDir(outDir).filePath(QStringLiteral("native_axis_manifest.json")), native, err_out))
        return false;

    QJsonObject provenance;
    provenance.insert(QStringLiteral("root720"), QFileInfo(options.root720).absoluteFilePath());
    provenance.insert(QStringLiteral("run1p9j"), QFileInfo(options.run1p9j).absoluteFilePath());
    provenance.insert(QStringLiteral("scoped_count"), static_cast<int>(ScopedProducerCatalog().size()));
    QJsonArray cutoffJson;
    for (int c : cutoffs) cutoffJson.append(c);
    provenance.insert(QStringLiteral("spatial_cutoffs_A"), cutoffJson);
    provenance.insert(QStringLiteral("exclusion_semantics"),
                      QStringLiteral("exclude target atom, directly bonded atoms, target-endpoint bonds, and self/member rings"));
    if (!writeJsonFile(QDir(outDir).filePath(QStringLiteral("provenance_manifest.json")),
                       provenance, err_out)) return false;

    QJsonObject regionManifest;
    regionManifest.insert(QStringLiteral("supplied"), region.supplied);
    regionManifest.insert(QStringLiteral("source_path"), region.sourcePath);
    regionManifest.insert(QStringLiteral("sha256"), region.sha256);
    regionManifest.insert(QStringLiteral("region_def_id"), region.regionDefId);
    regionManifest.insert(QStringLiteral("absence_reason"),
                          region.supplied ? QString() : QStringLiteral("region-spec-not-supplied"));
    if (!writeJsonFile(QDir(outDir).filePath(QStringLiteral("region_spec_manifest.json")),
                       regionManifest, err_out)) return false;

    QJsonObject nullSpec;
    nullSpec.insert(QStringLiteral("double_null"), QStringLiteral("NaN"));
    nullSpec.insert(QStringLiteral("int_null"), -1);
    nullSpec.insert(QStringLiteral("bool_null"), 0);
    nullSpec.insert(QStringLiteral("string_null"), QString());
    nullSpec.insert(QStringLiteral("sidecar_sanctioned_absence"), QStringLiteral("NaN with presence=0 and absence reason"));
    if (!writeJsonFile(QDir(outDir).filePath(QStringLiteral("null_spec.json")), nullSpec, err_out))
        return false;

    return writeContributorFamilyManifest(outDir, err_out);
}

bool survivalGate(const RowSchema& rowSchema,
                  const std::vector<std::unique_ptr<FieldSink>>& sinks,
                  std::size_t nRows,
                  QString* err_out) {
    static const std::set<QString> oldNames = {
        QStringLiteral("atom_name"),
        QStringLiteral("dssp_ss8"),
        QStringLiteral("dssp_hbond_energy"),
        QStringLiteral("prev_restype"),
        QStringLiteral("n_atoms_4A"),
        QStringLiteral("nearest_ring_dist"),
        QStringLiteral("self_or_bonded_atom_count"),
        QStringLiteral("tensor_frame"),
        QStringLiteral("valid_for_T2_model"),
    };
    for (const QString& c : rowSchema.columns) {
        if (oldNames.count(c)) {
            if (err_out) *err_out = QStringLiteral("survival gate failed: old row-design column reached rows.csv: %1").arg(c);
            return false;
        }
    }
    if (ScopedProducerCatalog().size() != 131) {
        if (err_out) *err_out = QStringLiteral("survival gate failed: scoped catalog count is %1, expected 131")
                                    .arg(ScopedProducerCatalog().size());
        return false;
    }
    std::size_t onep9jOnly = 0;
    for (const auto& sink : sinks) {
        const FieldSupport& st = sink->support;
        if (st.populated == 0 && sink->route != FieldRoute::Topology) {
            bool sanctioned = false;
            for (const auto& ds : st.byDataset) {
                if (ds.second.populated == 0
                    && ds.second.absence.size() == 1
                    && ds.second.absence.count(QStringLiteral("not-produced-in-dataset"))) {
                    sanctioned = true;
                }
            }
            if (!sanctioned) {
                if (err_out) *err_out = QStringLiteral("survival gate failed: scoped field has zero populated rows without sanctioned absence: %1")
                                            .arg(fieldStem(*sink->spec));
                return false;
            }
        }
        auto d720 = st.byDataset.find(QStringLiteral("720_static"));
        auto d1p9j = st.byDataset.find(QStringLiteral("1p9j"));
        if (d720 != st.byDataset.end() && d1p9j != st.byDataset.end()
            && d720->second.populated == 0 && d1p9j->second.populated > 0
            && d720->second.absence.size() == 1
            && d720->second.absence.count(QStringLiteral("not-produced-in-dataset"))) {
            ++onep9jOnly;
        }
    }
    if (onep9jOnly != 20) {
        if (err_out) *err_out = QStringLiteral("survival gate failed: expected 20 1P9J-only fields absent on 720 as not-produced-in-dataset, saw %1")
                                    .arg(onep9jOnly);
        return false;
    }
    for (const auto& sink : sinks) {
        if ((sink->route == FieldRoute::Measure || sink->route == FieldRoute::Incidence)
            && sink->rowAlignedRows != nRows) {
            if (err_out) *err_out = QStringLiteral("survival gate failed: row-aligned sidecar row mismatch for %1")
                                        .arg(fieldStem(*sink->spec));
            return false;
        }
    }
    return true;
}

}  // namespace
}  // namespace h5reader::rediscover

namespace h5reader::rediscover {
namespace {

std::vector<QString> enumerate720Poses(const QString& root720) {
    QDir root(root720);
    const QStringList names = root.entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
    return std::vector<QString>(names.begin(), names.end());
}

std::size_t countStaticRows(const QString& poseDir, QString* err_out) {
    const QString posPath = QDir(poseDir).filePath(QStringLiteral("pos.npy"));
    const io::QtNpyReader::WidenedArray pos = io::QtNpyReader::ReadArrayWidened(posPath);
    if (!pos.ok) {
        if (err_out) *err_out = pos.error;
        return 0;
    }
    return pos.rows;
}

std::size_t preflightRowCount(const QString& root720,
                              const std::vector<QString>& poses,
                              const RunData& onep9j,
                              QString* err_out) {
    std::size_t rows = 0;
    for (const QString& pose : poses) {
        const QString poseDir = QDir(root720).filePath(pose);
        rows += countStaticRows(poseDir, err_out);
        if (err_out && !err_out->isEmpty()) return 0;
    }
    if (onep9j.protein)
        rows += onep9j.protein->atomCount() * onep9j.frameMap.frameCount();
    return rows;
}

void writeTopologyForRun(const Body& body,
                         const RunContext& ctx,
                         const QString& outDir,
                         std::map<QString, FieldSupport*>& topologySupport) {
    const model::QtProtein& p = *body.run.protein;
    const QString prefix = QStringLiteral("%1,%2,%3")
                               .arg(csv(ctx.datasetId), csv(ctx.proteinId), csv(ctx.poseId));
    auto append = [&](const QString& stem, const QString& header, auto writer) {
        const QString path = QDir(outDir).filePath(QStringLiteral("topology/%1.csv").arg(stem));
        const bool exists = QFileInfo::exists(path);
        QFile f(path);
        if (!f.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append))
            throw std::runtime_error(QStringLiteral("cannot open topology CSV: %1").arg(path).toStdString());
        QTextStream out(&f);
        if (!exists) out << header << '\n';
        writer(out);
    };

    append(QStringLiteral("atoms_category_info"),
           QStringLiteral("dataset_id,protein_id,pose_id,atom_index,residue_index,element,iupac_name,bmrb_name,formal_charge,is_exchangeable"),
           [&](QTextStream& out) {
               for (std::size_t atom = 0; atom < p.atomCount(); ++atom) {
                   const model::QtAtom& a = p.atom(atom);
                   out << prefix << ',' << a.atomIndex << ',' << a.residueIndex << ','
                       << a.AtomicNumber() << ',' << csv(p.atomLabel(atom, model::NamingConvention::Iupac))
                       << ',' << csv(p.atomLabel(atom, model::NamingConvention::Bmrb)) << ','
                       << a.formalCharge << ',' << boolText(a.isExchangeable) << '\n';
               }
               if (auto it = topologySupport.find(QStringLiteral("atoms_category_info")); it != topologySupport.end()) {
                   it->second->count += p.atomCount();
                   it->second->populated += p.atomCount();
                   it->second->byDataset[ctx.datasetId].rows += p.atomCount();
                   it->second->byDataset[ctx.datasetId].populated += p.atomCount();
               }
           });
    append(QStringLiteral("residues"),
           QStringLiteral("dataset_id,protein_id,pose_id,residue_index,residue_number,residue_type,prev_residue_index,next_residue_index"),
           [&](QTextStream& out) {
               for (std::size_t ri = 0; ri < p.residueCount(); ++ri) {
                   const model::QtResidue& r = p.residue(ri);
                   out << prefix << ',' << r.residueIndex << ',' << r.address.residueNumber << ','
                       << csv(residueTypeLabel(p, &r)) << ',' << r.prevResidueIndex << ','
                       << r.nextResidueIndex << '\n';
               }
               if (auto it = topologySupport.find(QStringLiteral("residues")); it != topologySupport.end()) {
                   it->second->count += p.residueCount();
                   it->second->populated += p.residueCount();
                   it->second->byDataset[ctx.datasetId].rows += p.residueCount();
                   it->second->byDataset[ctx.datasetId].populated += p.residueCount();
               }
           });
    append(QStringLiteral("bonds"),
           QStringLiteral("dataset_id,protein_id,pose_id,bond_index,atom_a,atom_b,order,category,is_aromatic,is_backbone"),
           [&](QTextStream& out) {
               for (std::size_t bi = 0; bi < p.topology().bondCount(); ++bi) {
                   const model::QtBond& b = p.topology().bondAt(bi);
                   out << prefix << ',' << b.bondIndex << ',' << b.atomIndexA << ','
                       << b.atomIndexB << ',' << static_cast<int>(b.order) << ','
                       << static_cast<int>(b.category) << ',' << boolText(b.isAromatic) << ','
                       << boolText(b.isBackbone) << '\n';
               }
               if (auto it = topologySupport.find(QStringLiteral("bonds")); it != topologySupport.end()) {
                   it->second->count += p.topology().bondCount();
                   it->second->populated += p.topology().bondCount();
                   it->second->byDataset[ctx.datasetId].rows += p.topology().bondCount();
                   it->second->byDataset[ctx.datasetId].populated += p.topology().bondCount();
               }
           });
    append(QStringLiteral("rings"),
           QStringLiteral("dataset_id,protein_id,pose_id,ring_id,native_axis_index,ring_type,parent_residue_index,ring_size,is_aromatic,fused_partner_ring_id"),
           [&](QTextStream& out) {
               for (std::size_t ri = 0; ri < p.ringCount(); ++ri) {
                   const model::QtRing& r = p.ring(ri);
                   out << prefix << ',' << r.ringId << ',' << r.nativeAxisIndex << ','
                       << csv(QString::fromLatin1(r.TypeName())) << ',' << r.parentResidueIndex
                       << ',' << r.RingSizeValue() << ',' << boolText(r.IsAromatic()) << ','
                       << r.fusedPartnerRingId << '\n';
               }
               if (auto it = topologySupport.find(QStringLiteral("rings")); it != topologySupport.end()) {
                   it->second->count += p.ringCount();
                   it->second->populated += p.ringCount();
                   it->second->byDataset[ctx.datasetId].rows += p.ringCount();
                   it->second->byDataset[ctx.datasetId].populated += p.ringCount();
               }
           });
    append(QStringLiteral("ring_membership"),
           QStringLiteral("dataset_id,protein_id,pose_id,membership_row,ring_id,atom_index,ring_atom_order,is_vertex,is_substituent"),
           [&](QTextStream& out) {
               for (std::size_t mi = 0; mi < p.topology().ringMembershipCount(); ++mi) {
                   const model::QtRingMembership& m = p.topology().ringMembershipAt(mi);
                   out << prefix << ',' << mi << ',' << m.ringId << ',' << m.atomIndex << ','
                       << static_cast<int>(m.ringAtomOrder) << ',' << boolText(m.isVertex) << ','
                       << boolText(m.isSubstituent) << '\n';
               }
               if (auto it = topologySupport.find(QStringLiteral("ring_membership")); it != topologySupport.end()) {
                   it->second->count += p.topology().ringMembershipCount();
                   it->second->populated += p.topology().ringMembershipCount();
                   it->second->byDataset[ctx.datasetId].rows += p.topology().ringMembershipCount();
                   it->second->byDataset[ctx.datasetId].populated += p.topology().ringMembershipCount();
               }
           });
}

void emitRun(RunData& run,
             const RunContext& ctx,
             const QString& outDir,
             const RowSchema& rowSchema,
             const RegionApplySpec& region,
             const std::vector<int>& cutoffs,
             QTextStream& rowsOut,
             std::map<QString, ColumnSupport>& columnSupport,
             std::vector<std::unique_ptr<FieldSink>>& fieldSinks,
             TargetSinks& targets,
             E3nnSinks& e3nn,
             std::map<QString, FieldSupport*>& topologySupport,
             std::size_t* nextRowId) {
    Catalog catalog(run);
    ResidentIndexes indexes = BuildResidentIndexes(run);
    Body body{run, indexes, catalog};
    writeTopologyForRun(body, ctx, outDir, topologySupport);

    const std::size_t atomCount = run.protein->atomCount();
    for (std::size_t frame = 0; frame < run.frameMap.frameCount(); ++frame) {
        const std::int64_t graphOffset = static_cast<std::int64_t>(e3nn.graphAtoms);
        e3nn.graphOffsets.writeOne(graphOffset);
        const std::size_t graphId = e3nn.graphCount++;
        *e3nn.conformationsOut << graphId << ','
                                << csv(ctx.datasetId) << ','
                                << csv(ctx.proteinId) << ','
                                << csv(ctx.poseId) << ','
                                << frame << ','
                                << frame << ','
                                << run.frameMap.originalIndex(frame) << ','
                                << num(run.timePs(frame)) << '\n';

        std::unordered_map<const io::FieldSpec*, std::vector<std::vector<std::size_t>>> incidence;
        std::unordered_map<const io::FieldSpec*, std::size_t> nativeCounts;
        for (const auto& sink : fieldSinks) {
            if (sink->route != FieldRoute::Incidence && sink->route != FieldRoute::NativePair) continue;
            const std::size_t nativeCount = catalog.nativeRowCount(body, sink->spec->kind, frame);
            nativeCounts[sink->spec] = nativeCount;
            if (sink->route == FieldRoute::Incidence) {
                auto incident = buildIncidentRows(body, *sink->spec, frame, nativeCount);
                emitIncidenceNative(*sink, body, ctx, frame, nativeCount, incident);
                incidence[sink->spec] = std::move(incident);
            } else {
                emitNativePairField(*sink, body, ctx, frame, nativeCount);
            }
        }

        for (std::size_t atom = 0; atom < atomCount; ++atom) {
            const std::size_t rowId = (*nextRowId)++;
            RowBuilder row(rowSchema);
            fillRowLabels(row, body, ctx, region, cutoffs, rowId, atom, frame);
            writeRow(rowsOut, rowSchema, row, columnSupport);

            const bool dft = targetPresent(run, atom, frame);
            writeTargets(targets, e3nn, run, atom, frame);
            writeE3nnAtom(e3nn, body, graphId, atom, rowId, atom, frame);

            for (const auto& sink : fieldSinks) {
                switch (sink->route) {
                case FieldRoute::Measure:
                    emitMeasureField(*sink, body, ctx, rowId, atom, frame);
                    break;
                case FieldRoute::Incidence:
                    emitIncidenceRow(*sink, body, ctx, rowId, atom, frame, incidence[sink->spec]);
                    break;
                case FieldRoute::Target: {
                    DatasetSupport& ds = sink->support.byDataset[ctx.datasetId];
                    ds.rows++;
                    sink->support.count++;
                    if (dft) {
                        ds.populated++;
                        sink->support.populated++;
                    } else {
                        addAbsence(ds.absence, QStringLiteral("target-not-present"));
                        addAbsence(sink->support.absence, QStringLiteral("target-not-present"));
                    }
                    break;
                }
                case FieldRoute::NativePair:
                case FieldRoute::Topology:
                    break;
                }
            }
        }
    }
}

}  // namespace
}  // namespace h5reader::rediscover

namespace h5reader::rediscover {
namespace {

void recordDatasetAbsence(FieldSink& sink,
                          const QString& datasetId,
                          const QString& reason,
                          std::size_t n = 1) {
    DatasetSupport& ds = sink.support.byDataset[datasetId];
    ds.rows += n;
    addAbsence(ds.absence, reason, n);
    addAbsence(sink.support.absence, reason, n);
}

void recordDatasetPresent(FieldSink& sink,
                          const QString& datasetId,
                          const std::vector<double>& values) {
    DatasetSupport& ds = sink.support.byDataset[datasetId];
    ds.rows++;
    ds.populated++;
    sink.support.count++;
    sink.support.populated++;
    if (sink.support.componentPopulated.size() < values.size())
        sink.support.componentPopulated.resize(values.size(), 0);
    for (std::size_t i = 0; i < values.size(); ++i)
        if (std::isfinite(values[i])) sink.support.componentPopulated[i]++;
}

void writeMeasureAbsent(FieldSink& sink,
                        const QString& datasetId,
                        std::size_t rowId,
                        const QString& reason) {
    const std::vector<double> nan(sink.cols, kNaN);
    sink.values->writeMany(nan.data(), nan.size());
    const std::uint8_t p = 0;
    sink.presence->writeOne(p);
    if (sink.absenceOut) *sink.absenceOut << rowId << ',' << csv(reason) << ",1\n";
    ++sink.rowAlignedRows;
    recordDatasetAbsence(sink, datasetId, reason);
}

void writeMeasurePresent(FieldSink& sink,
                         const QString& datasetId,
                         const std::vector<double>& values) {
    sink.values->writeMany(values.data(), values.size());
    const std::uint8_t p = 1;
    sink.presence->writeOne(p);
    ++sink.rowAlignedRows;
    recordDatasetPresent(sink, datasetId, values);
}

void emitMeasureField(FieldSink& sink,
                      const Body& body,
                      const RunContext& ctx,
                      std::size_t rowId,
                      std::size_t atom,
                      std::size_t frame) {
    QString providerReason;
    if (!body.catalog.has(body.run, sink.spec->kind)) {
        body.catalog.provider(body.run, sink.spec->kind, &providerReason);
        writeMeasureAbsent(sink, ctx.datasetId, rowId,
                           providerReason.isEmpty() ? QStringLiteral("not-produced-in-dataset")
                                                    : providerReason);
        return;
    }
    QString reason;
    const auto native = nativeRowForAtom(body, *sink.spec, atom, &reason);
    if (!native) {
        writeMeasureAbsent(sink, ctx.datasetId, rowId, reason);
        return;
    }
    std::vector<double> values;
    values.reserve(sink.cols);
    for (std::size_t c = 0; c < sink.cols; ++c) {
        QString valueReason;
        const std::optional<double> v =
            body.catalog.value(body, sink.spec->kind, *native, frame, static_cast<int>(c), &valueReason);
        if (!v) {
            writeMeasureAbsent(sink, ctx.datasetId, rowId, valueReason);
            return;
        }
        values.push_back(*v);
    }
    writeMeasurePresent(sink, ctx.datasetId, values);
}

void writeIncidenceAbsent(FieldSink& sink,
                          const QString& datasetId,
                          std::size_t rowId,
                          const QString& reason) {
    std::vector<double> out(4 * sink.cols, kNaN);
    for (std::size_t c = 0; c < sink.cols; ++c) out[3 * sink.cols + c] = 0.0;
    sink.byAtom->writeMany(out.data(), out.size());
    const std::uint8_t p = 0;
    sink.presence->writeOne(p);
    if (sink.absenceOut) *sink.absenceOut << rowId << ',' << csv(reason) << ",1\n";
    ++sink.rowAlignedRows;
    recordDatasetAbsence(sink, datasetId, reason);
}

void writeIncidencePresent(FieldSink& sink,
                           const QString& datasetId,
                           const std::vector<double>& out) {
    sink.byAtom->writeMany(out.data(), out.size());
    const std::uint8_t p = 1;
    sink.presence->writeOne(p);
    ++sink.rowAlignedRows;
    DatasetSupport& ds = sink.support.byDataset[datasetId];
    ds.rows++;
    ds.populated++;
    sink.support.count++;
    sink.support.populated++;
}

void emitIncidenceRow(FieldSink& sink,
                      const Body& body,
                      const RunContext& ctx,
                      std::size_t rowId,
                      std::size_t atom,
                      std::size_t frame,
                      const std::vector<std::vector<std::size_t>>& incidentRows) {
    QString providerReason;
    if (!body.catalog.has(body.run, sink.spec->kind)) {
        body.catalog.provider(body.run, sink.spec->kind, &providerReason);
        writeIncidenceAbsent(sink, ctx.datasetId, rowId,
                             providerReason.isEmpty() ? QStringLiteral("not-produced-in-dataset")
                                                      : providerReason);
        return;
    }
    if (atom >= incidentRows.size() || incidentRows[atom].empty()) {
        writeIncidenceAbsent(sink, ctx.datasetId, rowId, QStringLiteral("no-incident-native-row"));
        return;
    }

    std::vector<double> sum(sink.cols, 0.0);
    std::vector<double> maxAbs(sink.cols, 0.0);
    std::vector<std::size_t> count(sink.cols, 0);
    for (std::size_t native : incidentRows[atom]) {
        for (std::size_t c = 0; c < sink.cols; ++c) {
            QString reason;
            const std::optional<double> v =
                body.catalog.value(body, sink.spec->kind, native, frame, static_cast<int>(c), &reason);
            if (!v || !std::isfinite(*v)) continue;
            sum[c] += *v;
            maxAbs[c] = std::max(maxAbs[c], std::abs(*v));
            count[c]++;
        }
    }
    bool any = false;
    for (std::size_t n : count)
        if (n > 0) any = true;
    if (!any) {
        writeIncidenceAbsent(sink, ctx.datasetId, rowId, QStringLiteral("incident-values-absent"));
        return;
    }
    std::vector<double> out(4 * sink.cols, kNaN);
    for (std::size_t c = 0; c < sink.cols; ++c) {
        if (count[c] > 0) {
            out[c] = sum[c];
            out[sink.cols + c] = sum[c] / static_cast<double>(count[c]);
            out[2 * sink.cols + c] = maxAbs[c];
        }
        out[3 * sink.cols + c] = static_cast<double>(count[c]);
    }
    writeIncidencePresent(sink, ctx.datasetId, out);
}

std::vector<std::vector<std::size_t>>
buildIncidentRows(const Body& body,
                  const io::FieldSpec& spec,
                  std::size_t frame,
                  std::size_t nativeCount) {
    const std::size_t atoms = body.run.protein ? body.run.protein->atomCount() : 0;
    std::vector<std::vector<std::size_t>> out(atoms);
    for (std::size_t native = 0; native < nativeCount; ++native) {
        const auto endpoints = nativeEndpoints(body, spec, native, frame);
        if (!endpoints) continue;
        const int32_t a = endpoints->first;
        const int32_t b = endpoints->second;
        if (a >= 0 && static_cast<std::size_t>(a) < atoms) out[static_cast<std::size_t>(a)].push_back(native);
        if (b >= 0 && static_cast<std::size_t>(b) < atoms && b != a)
            out[static_cast<std::size_t>(b)].push_back(native);
    }
    return out;
}

void emitIncidenceNative(FieldSink& sink,
                         const Body& body,
                         const RunContext& ctx,
                         std::size_t frame,
                         std::size_t nativeCount,
                         const std::vector<std::vector<std::size_t>>& incidentRows) {
    (void)incidentRows;
    for (std::size_t native = 0; native < nativeCount; ++native) {
        std::vector<double> values;
        values.reserve(sink.cols);
        bool present = true;
        QString reason;
        for (std::size_t c = 0; c < sink.cols; ++c) {
            const std::optional<double> v =
                body.catalog.value(body, sink.spec->kind, native, frame, static_cast<int>(c), &reason);
            if (!v) {
                present = false;
                break;
            }
            values.push_back(*v);
        }
        if (!present) continue;
        const std::size_t globalNative = sink.nativeRows++;
        sink.nativeValues->writeMany(values.data(), values.size());
        const auto endpoints = nativeEndpoints(body, *sink.spec, native, frame);
        if (endpoints && sink.nativeMapOut) {
            const std::array<std::pair<int32_t, const char*>, 2> roles = {
                std::make_pair(endpoints->first, "endpoint_a"),
                std::make_pair(endpoints->second, "endpoint_b"),
            };
            for (const auto& role : roles) {
                if (role.first < 0
                    || static_cast<std::size_t>(role.first) >= body.run.protein->atomCount()) {
                    continue;
                }
                const std::size_t rowId =
                    ctx.rowBase + frame * body.run.protein->atomCount()
                    + static_cast<std::size_t>(role.first);
                *sink.nativeMapOut << globalNative << ','
                                   << csv(ctx.datasetId) << ','
                                   << csv(ctx.proteinId) << ','
                                   << csv(ctx.poseId) << ','
                                   << frame << ','
                                   << role.first << ','
                                   << rowId << ','
                                   << csv(axisName(sink.spec->axis)) << ','
                                   << native << ','
                                   << role.second << '\n';
            }
        }
    }
}

void emitNativePairField(FieldSink& sink,
                         const Body& body,
                         const RunContext& ctx,
                         std::size_t frame,
                         std::size_t nativeCount) {
    const StaticNpyArray* pairAnchor = body.run.producerArray(io::FieldKind::RingContributions);
    const StaticNpyArray* fieldArray = body.run.producerArray(sink.spec->kind);
    if (!pairAnchor || !fieldArray) {
        recordDatasetAbsence(sink, ctx.datasetId, QStringLiteral("not-produced-in-dataset"),
                             body.run.protein->atomCount());
        return;
    }
    if (pairAnchor->rowsForFrame(frame) != fieldArray->rowsForFrame(frame)) {
        throw std::runtime_error(QStringLiteral("RingContributionPair native order check failed for %1: row-count mismatch")
                                     .arg(fieldStem(*sink.spec)).toStdString());
    }
    for (std::size_t native = 0; native < nativeCount; ++native) {
        const std::size_t anchorRow = pairAnchor->rowFor(native, frame);
        if (pairAnchor->cols < 3 || anchorRow >= pairAnchor->rows)
            throw std::runtime_error("RingContributionPair native order check failed: malformed ring_contributions");
        const int32_t atomIndex =
            static_cast<int32_t>(std::llround(pairAnchor->value(anchorRow, 0)));
        const int32_t ringId =
            static_cast<int32_t>(std::llround(pairAnchor->value(anchorRow, 1)));
        const int32_t ringType =
            static_cast<int32_t>(std::llround(pairAnchor->value(anchorRow, 2)));
        if (atomIndex < 0 || static_cast<std::size_t>(atomIndex) >= body.run.protein->atomCount()
            || ringId < 0 || static_cast<std::size_t>(ringId) >= body.run.protein->ringCount()) {
            throw std::runtime_error(QStringLiteral("RingContributionPair native order check failed for %1: invalid atom/ring ids")
                                         .arg(fieldStem(*sink.spec)).toStdString());
        }

        std::vector<double> values;
        values.reserve(sink.cols);
        QString reason;
        for (std::size_t c = 0; c < sink.cols; ++c) {
            const std::optional<double> v =
                body.catalog.value(body, sink.spec->kind, native, frame, static_cast<int>(c), &reason);
            if (!v) {
                throw std::runtime_error(QStringLiteral("RingContributionPair value read failed for %1: %2")
                                             .arg(fieldStem(*sink.spec), reason).toStdString());
            }
            values.push_back(*v);
        }
        const std::size_t globalNative = sink.nativeRows++;
        sink.nativeValues->writeMany(values.data(), values.size());
        sink.support.count++;
        sink.support.populated++;
        sink.support.byDataset[ctx.datasetId].rows++;
        sink.support.byDataset[ctx.datasetId].populated++;
        if (sink.nativeMapOut) {
            *sink.nativeMapOut << globalNative << ','
                               << csv(ctx.datasetId) << ','
                               << csv(ctx.proteinId) << ','
                               << csv(ctx.poseId) << ','
                               << frame << ','
                               << atomIndex << ','
                               << ringId << ','
                               << ringType << ','
                               << "atom_ring_pair" << '\n';
        }
    }
}

void writeTargets(TargetSinks& targets,
                  E3nnSinks& e3nn,
                  const RunData& run,
                  std::size_t atom,
                  std::size_t frame) {
    const std::size_t original = run.frameMap.originalIndex(frame);
    const model::DftAtomShielding* dft = run.dft.AtomShielding(atom, original);
    const bool present = dft != nullptr;
    writeTensorTarget(targets.totalRaw, targets.totalT0, targets.totalT1, targets.totalT2,
                      present ? dft->total_raw : model::Mat3::Zero(), present);
    writeTensorTarget(targets.diaRaw, targets.diaT0, targets.diaT1, targets.diaT2,
                      present ? dft->dia_raw : model::Mat3::Zero(), present);
    writeTensorTarget(targets.paraRaw, targets.paraT0, targets.paraT1, targets.paraT2,
                      present ? dft->para_raw : model::Mat3::Zero(), present);
    const std::uint8_t p = present ? 1 : 0;
    targets.present.writeOne(p);

    const std::array<double, 9> nan9{kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, kNaN};
    if (present) {
        const auto total = matArray(dft->total_raw);
        const auto dia = matArray(dft->dia_raw);
        const auto para = matArray(dft->para_raw);
        e3nn.sigmaTotalRaw.writeMany(total.data(), total.size());
        e3nn.sigmaDiaRaw.writeMany(dia.data(), dia.size());
        e3nn.sigmaParaRaw.writeMany(para.data(), para.size());
    } else {
        e3nn.sigmaTotalRaw.writeMany(nan9.data(), nan9.size());
        e3nn.sigmaDiaRaw.writeMany(nan9.data(), nan9.size());
        e3nn.sigmaParaRaw.writeMany(nan9.data(), nan9.size());
    }
}

void writeE3nnAtom(E3nnSinks& e3nn,
                   const Body& body,
                   std::size_t graphId,
                   std::size_t graphAtomOffset,
                   std::size_t rowId,
                   std::size_t atom,
                   std::size_t frame) {
    const model::QtAtom& a = body.run.protein->atom(atom);
    const std::int32_t atomIndex = a.atomIndex;
    const std::int16_t z = static_cast<std::int16_t>(a.AtomicNumber());
    const model::Vec3 pos = body.run.conformation->atomPosition(frame, atom);
    const std::array<double, 3> p{pos.x(), pos.y(), pos.z()};
    e3nn.atomIndex.writeOne(atomIndex);
    e3nn.z.writeOne(z);
    e3nn.positions.writeMany(p.data(), p.size());
    QString reason;
    bool embeddingPresent = true;
    for (std::size_t c = 0; c < 256; ++c) {
        const std::optional<double> v =
            body.catalog.value(body, io::FieldKind::AIMNet2Aim, atom, frame, static_cast<int>(c), &reason);
        const float f = v ? static_cast<float>(*v)
                          : std::numeric_limits<float>::quiet_NaN();
        if (!v) embeddingPresent = false;
        e3nn.aimnet2Embedding.writeOne(f);
    }
    (void)embeddingPresent;
    const std::array<std::int64_t, 2> map{
        static_cast<std::int64_t>(graphId),
        static_cast<std::int64_t>(graphAtomOffset)};
    e3nn.rowToGraphAtom.writeMany(map.data(), map.size());
    e3nn.graphAtoms++;
    (void)rowId;
}

}  // namespace
}  // namespace h5reader::rediscover

namespace h5reader::rediscover {
namespace {

RowSchema buildRowSchema(const std::vector<int>& cutoffs) {
    RowSchema s;
    auto add = [&](const QString& name) {
        s.index[name] = static_cast<std::size_t>(s.columns.size());
        s.columns << name;
    };
    const QStringList fixed = {
        QStringLiteral("row_id"),
        QStringLiteral("global_row_key"),
        QStringLiteral("dataset_id"),
        QStringLiteral("protein_id"),
        QStringLiteral("pose_id"),
        QStringLiteral("conformation_id"),
        QStringLiteral("pose_kind"),
        QStringLiteral("split_group_id"),
        QStringLiteral("frame_slot"),
        QStringLiteral("h5_row"),
        QStringLiteral("original_frame_index"),
        QStringLiteral("time_ps"),
        QStringLiteral("atom_index"),
        QStringLiteral("element"),
        QStringLiteral("residue_index"),
        QStringLiteral("residue_number"),
        QStringLiteral("residue_type"),
        QStringLiteral("iupac_name"),
        QStringLiteral("bmrb_name"),
        QStringLiteral("chemical_category"),
        QStringLiteral("ring_membership_id"),
        QStringLiteral("ring_type"),
        QStringLiteral("in_aromatic_ring"),
        QStringLiteral("formal_charge"),
        QStringLiteral("is_exchangeable"),
        QStringLiteral("polar_h_kind"),
        QStringLiteral("backbone_role"),
        QStringLiteral("locant"),
        QStringLiteral("branch_outer"),
        QStringLiteral("branch_inner"),
        QStringLiteral("di_index"),
        QStringLiteral("prochiral"),
        QStringLiteral("equivalence_class"),
        QStringLiteral("prev_residue_type"),
        QStringLiteral("next_residue_type"),
        QStringLiteral("prev_residue_class"),
        QStringLiteral("next_residue_class"),
        QStringLiteral("pre_proline"),
        QStringLiteral("is_pro"),
        QStringLiteral("is_gly"),
        QStringLiteral("terminal_state"),
        QStringLiteral("ss8"),
        QStringLiteral("ss3"),
        QStringLiteral("dssp_hbond_energy_min"),
        QStringLiteral("phi"),
        QStringLiteral("psi"),
        QStringLiteral("omega"),
        QStringLiteral("phi_bin"),
        QStringLiteral("psi_bin"),
        QStringLiteral("omega_bin"),
        QStringLiteral("chi1_sin"),
        QStringLiteral("chi1_cos"),
        QStringLiteral("chi1_exists"),
        QStringLiteral("chi2_sin"),
        QStringLiteral("chi2_cos"),
        QStringLiteral("chi2_exists"),
        QStringLiteral("chi3_sin"),
        QStringLiteral("chi3_cos"),
        QStringLiteral("chi3_exists"),
        QStringLiteral("chi4_sin"),
        QStringLiteral("chi4_cos"),
        QStringLiteral("chi4_exists"),
        QStringLiteral("region_def_id"),
        QStringLiteral("rama_region_id"),
        QStringLiteral("rama_region_label"),
        QStringLiteral("rama_region_prob"),
        QStringLiteral("rotamer_id"),
        QStringLiteral("rotamer_label"),
        QStringLiteral("rotamer_prob"),
        QStringLiteral("folding_rule_id"),
        QStringLiteral("prop_cutoff_set_id"),
        QStringLiteral("prop_self_or_bonded_atom_count"),
        QStringLiteral("prop_self_or_bonded_bond_count"),
        QStringLiteral("prop_has_self_or_bonded_driver"),
        QStringLiteral("prop_nearest_ring_dist"),
        QStringLiteral("prop_nearest_ring_id"),
        QStringLiteral("prop_nearest_charge_dist"),
        QStringLiteral("prop_nearest_charge_atom"),
        QStringLiteral("prop_nearest_charge_sign"),
        QStringLiteral("prop_nearest_bond_dist"),
        QStringLiteral("prop_nearest_bond_id"),
        QStringLiteral("prop_nearest_atom_dist"),
        QStringLiteral("prop_nearest_atom_index"),
        QStringLiteral("prop_ring_cyl_z"),
        QStringLiteral("prop_ring_cyl_rho"),
        QStringLiteral("prop_ring_angle_to_normal"),
        QStringLiteral("prop_ring_cos_phi"),
        QStringLiteral("prop_ring_sin_phi"),
        QStringLiteral("dft_present"),
        QStringLiteral("target_total_present"),
        QStringLiteral("target_dia_present"),
        QStringLiteral("target_para_present"),
        QStringLiteral("aimnet2_embedding_present"),
    };
    for (const QString& name : fixed) add(name);
    for (int cutoff : cutoffs) {
        add(QStringLiteral("prop_n_atoms_%1A").arg(cutoff));
        add(QStringLiteral("prop_n_rings_%1A").arg(cutoff));
        add(QStringLiteral("prop_n_charges_%1A").arg(cutoff));
        add(QStringLiteral("prop_n_bonds_%1A").arg(cutoff));
    }
    return s;
}

void writeRow(QTextStream& out,
              const RowSchema& schema,
              const RowBuilder& row,
              std::map<QString, ColumnSupport>& support) {
    for (int i = 0; i < schema.columns.size(); ++i) {
        if (i) out << ',';
        out << csv(row.values[static_cast<std::size_t>(i)]);
        ColumnSupport& st = support[schema.columns[i]];
        st.name = schema.columns[i];
        st.kind = QStringLiteral("selector_label");
        st.source = QStringLiteral("ConsolidatedEmit");
        st.rowCount++;
        if (!row.values[static_cast<std::size_t>(i)].isEmpty()) st.populated++;
        else addAbsence(st.absence, QStringLiteral("null-or-sanctioned-empty"));
    }
    out << '\n';
}

std::optional<double> fieldValue(const Body& body,
                                 io::FieldKind kind,
                                 std::size_t nativeRow,
                                 std::size_t frame,
                                 int component) {
    QString reason;
    return body.catalog.value(body, kind, nativeRow, frame, component, &reason);
}

double dsspMinEnergy(const Body& body, std::size_t atom, std::size_t frame) {
    double out = kNaN;
    for (int c = 0; c < 4; ++c) {
        const std::optional<double> v =
            fieldValue(body, io::FieldKind::DSSPHBondEnergy, atom, frame, c);
        if (v && std::isfinite(*v) && *v != 0.0)
            out = std::isfinite(out) ? std::min(out, *v) : *v;
    }
    return out;
}

bool targetPresent(const RunData& run, std::size_t atom, std::size_t frame) {
    return run.dft.AtomShielding(atom, run.frameMap.originalIndex(frame)) != nullptr;
}

bool aimnet2EmbeddingPresent(const Body& body, std::size_t atom, std::size_t frame) {
    QString reason;
    return body.catalog.present(body, io::FieldKind::AIMNet2Aim, atom, frame, 0, &reason);
}

void fillRowLabels(RowBuilder& row,
                   const Body& body,
                   const RunContext& ctx,
                   const RegionApplySpec& region,
                   const std::vector<int>& cutoffs,
                   std::size_t rowId,
                   std::size_t atom,
                   std::size_t frame) {
    const model::QtProtein& p = *body.run.protein;
    const model::QtAtom& a = p.atom(atom);
    const model::QtResidue* r = residueForAtom(p, a);
    const std::size_t original = body.run.frameMap.originalIndex(frame);
    const QString globalKey = QStringLiteral("%1|%2|%3|%4|%5")
                                  .arg(ctx.datasetId, ctx.proteinId, ctx.poseId)
                                  .arg(frame)
                                  .arg(a.atomIndex);
    const QString iupac = p.atomLabel(atom, model::NamingConvention::Iupac);
    const QString bmrb = p.atomLabel(atom, model::NamingConvention::Bmrb);
    const std::set<int32_t> bonded = bondedAtoms(p, atom);

    row.setInt(QStringLiteral("row_id"), static_cast<qint64>(rowId));
    row.set(QStringLiteral("global_row_key"), globalKey);
    row.set(QStringLiteral("dataset_id"), ctx.datasetId);
    row.set(QStringLiteral("protein_id"), ctx.proteinId);
    row.set(QStringLiteral("pose_id"), ctx.poseId);
    row.set(QStringLiteral("conformation_id"), ctx.conformationId);
    row.set(QStringLiteral("pose_kind"), ctx.poseKind);
    row.set(QStringLiteral("split_group_id"), ctx.splitGroupId);
    row.setInt(QStringLiteral("frame_slot"), static_cast<qint64>(frame));
    row.setInt(QStringLiteral("h5_row"), static_cast<qint64>(frame));
    row.setInt(QStringLiteral("original_frame_index"), static_cast<qint64>(original));
    row.setDouble(QStringLiteral("time_ps"), body.run.timePs(frame));
    row.setInt(QStringLiteral("atom_index"), a.atomIndex);
    row.setInt(QStringLiteral("element"), a.AtomicNumber());
    row.setInt(QStringLiteral("residue_index"), a.residueIndex);
    row.setInt(QStringLiteral("residue_number"), r ? r->address.residueNumber : 0);
    row.set(QStringLiteral("residue_type"), residueTypeLabel(p, r));
    row.set(QStringLiteral("iupac_name"), iupac);
    row.set(QStringLiteral("bmrb_name"), bmrb);
    row.set(QStringLiteral("chemical_category"), body.idx.chemicalCategories.categoryForAtom(atom));
    row.set(QStringLiteral("ring_membership_id"), ringMembershipId(p, atom));
    row.set(QStringLiteral("ring_type"), ringTypeForAtom(p, atom));
    row.setBool(QStringLiteral("in_aromatic_ring"), a.aromatic);
    row.setInt(QStringLiteral("formal_charge"), a.formalCharge);
    row.setBool(QStringLiteral("is_exchangeable"), a.isExchangeable);
    row.setInt(QStringLiteral("polar_h_kind"), static_cast<int>(a.polarH));
    row.setInt(QStringLiteral("backbone_role"), static_cast<int>(a.backboneRole));
    row.setInt(QStringLiteral("locant"), static_cast<int>(a.locant));
    row.setInt(QStringLiteral("branch_outer"), a.branch.outer);
    row.setInt(QStringLiteral("branch_inner"), a.branch.inner);
    row.setInt(QStringLiteral("di_index"), static_cast<int>(a.diIndex));
    row.setInt(QStringLiteral("prochiral"), static_cast<int>(a.prochiral));
    row.setInt(QStringLiteral("equivalence_class"), a.equivalenceClass);

    row.set(QStringLiteral("prev_residue_type"), r ? residueTypeLabelForAa(r->prevResidueType) : QString());
    row.set(QStringLiteral("next_residue_type"), r ? residueTypeLabelForAa(r->nextResidueType) : QString());
    row.set(QStringLiteral("prev_residue_class"),
            r ? residueClassName(ClassifyResidue(r->prevResidueType)) : QString());
    row.set(QStringLiteral("next_residue_class"),
            r ? residueClassName(ClassifyResidue(r->nextResidueType)) : QString());
    row.setBool(QStringLiteral("pre_proline"), r && r->nextResidueType == model::AminoAcid::PRO);
    row.setBool(QStringLiteral("is_pro"), r && r->isProline);
    row.setBool(QStringLiteral("is_gly"), r && r->aminoAcid == model::AminoAcid::GLY);
    row.setInt(QStringLiteral("terminal_state"), r ? static_cast<int>(r->terminalState) : -1);

    const SecondaryStructureState ss = body.idx.secondaryStructure.state(atom, frame);
    row.set(QStringLiteral("ss8"), ss8Label(ss.ss8));
    row.set(QStringLiteral("ss3"), ss3Label(ss.ss3));
    row.setDouble(QStringLiteral("dssp_hbond_energy_min"), dsspMinEnergy(body, atom, frame));

    const DihedralState phi = body.idx.dihedrals.state(DihedralKind::Phi, atom, frame);
    const DihedralState psi = body.idx.dihedrals.state(DihedralKind::Psi, atom, frame);
    const DihedralState omega = body.idx.dihedrals.state(DihedralKind::Omega, atom, frame);
    row.setDouble(QStringLiteral("phi"), phi.present ? phi.radians : kNaN);
    row.setDouble(QStringLiteral("psi"), psi.present ? psi.radians : kNaN);
    row.setDouble(QStringLiteral("omega"), omega.present ? omega.radians : kNaN);
    row.setInt(QStringLiteral("phi_bin"), phi.present ? phi.fixed_bin : -1);
    row.setInt(QStringLiteral("psi_bin"), psi.present ? psi.fixed_bin : -1);
    row.setInt(QStringLiteral("omega_bin"), omega.present ? omega.fixed_bin : -1);
    for (int k = 0; k < 4; ++k) {
        const auto kind = static_cast<DihedralKind>(static_cast<int>(DihedralKind::Chi1) + k);
        const DihedralState chi = body.idx.dihedrals.state(kind, atom, frame);
        row.setDouble(QStringLiteral("chi%1_sin").arg(k + 1),
                      chi.present ? std::sin(chi.radians) : kNaN);
        row.setDouble(QStringLiteral("chi%1_cos").arg(k + 1),
                      chi.present ? std::cos(chi.radians) : kNaN);
        row.setBool(QStringLiteral("chi%1_exists").arg(k + 1), chi.present);
    }

    if (region.supplied) {
        row.set(QStringLiteral("region_def_id"), region.regionDefId);
        const auto basin = nearestRamaBasin(region,
                                            phi.present ? phi.radians : kNaN,
                                            psi.present ? psi.radians : kNaN);
        if (basin) {
            row.set(QStringLiteral("rama_region_id"), basin->id);
            row.set(QStringLiteral("rama_region_label"), basin->label);
            row.setDouble(QStringLiteral("rama_region_prob"), 1.0);
        }
    }

    QStringList cutoffLabels;
    for (int cutoff : cutoffs) cutoffLabels << QString::number(cutoff);
    row.set(QStringLiteral("prop_cutoff_set_id"),
            QStringLiteral("spatial_cutoffs_A:%1").arg(cutoffLabels.join(QLatin1Char('|'))));
    row.setInt(QStringLiteral("prop_self_or_bonded_atom_count"),
               static_cast<qint64>(bonded.size() + 1));
    row.setInt(QStringLiteral("prop_self_or_bonded_bond_count"),
               static_cast<qint64>(p.topology().bondIndicesForAtom(atom).size()));
    row.setBool(QStringLiteral("prop_has_self_or_bonded_driver"), !bonded.empty());
    for (int cutoff : cutoffs) {
        row.setInt(QStringLiteral("prop_n_atoms_%1A").arg(cutoff),
                   static_cast<qint64>(countNear(body, CloudKind::Atoms, atom, frame, cutoff, bonded)));
        row.setInt(QStringLiteral("prop_n_rings_%1A").arg(cutoff),
                   static_cast<qint64>(countNear(body, CloudKind::RingCenters, atom, frame, cutoff, bonded)));
        row.setInt(QStringLiteral("prop_n_charges_%1A").arg(cutoff),
                   static_cast<qint64>(countNear(body, CloudKind::ChargeSites, atom, frame, cutoff, bonded)));
        row.setInt(QStringLiteral("prop_n_bonds_%1A").arg(cutoff),
                   static_cast<qint64>(countNear(body, CloudKind::AllBondMidpoints, atom, frame, cutoff, bonded)));
    }

    int32_t nearestRing = -1;
    int32_t nearestCharge = -1;
    int32_t nearestBond = -1;
    int32_t nearestAtom = -1;
    row.setDouble(QStringLiteral("prop_nearest_ring_dist"),
                  nearestDist(body, CloudKind::RingCenters, atom, frame, bonded, &nearestRing));
    row.setInt(QStringLiteral("prop_nearest_ring_id"), nearestRing);
    row.setDouble(QStringLiteral("prop_nearest_charge_dist"),
                  nearestDist(body, CloudKind::ChargeSites, atom, frame, bonded, &nearestCharge));
    row.setInt(QStringLiteral("prop_nearest_charge_atom"), nearestCharge);
    if (nearestCharge >= 0) {
        QString reason;
        const std::optional<double> q =
            body.catalog.value(body, io::FieldKind::FfPartialCharge,
                               static_cast<std::size_t>(nearestCharge), frame, 0, &reason);
        row.setDouble(QStringLiteral("prop_nearest_charge_sign"), q ? *q : kNaN);
    }
    row.setDouble(QStringLiteral("prop_nearest_bond_dist"),
                  nearestDist(body, CloudKind::AllBondMidpoints, atom, frame, bonded, &nearestBond));
    row.setInt(QStringLiteral("prop_nearest_bond_id"), nearestBond);
    row.setDouble(QStringLiteral("prop_nearest_atom_dist"),
                  nearestDist(body, CloudKind::Atoms, atom, frame, bonded, &nearestAtom));
    row.setInt(QStringLiteral("prop_nearest_atom_index"), nearestAtom);
    if (nearestRing >= 0 && static_cast<std::size_t>(nearestRing) < body.idx.ringGeometry.ringCount()) {
        const model::RingGeometry& g =
            body.idx.ringGeometry.at(static_cast<std::size_t>(nearestRing), frame);
        const model::Vec3 pos = body.run.conformation->atomPosition(frame, atom);
        const model::Vec3 disp = pos - g.center;
        const double z = disp.dot(g.normal);
        const model::Vec3 plane = disp - z * g.normal;
        const double rho = plane.norm();
        row.setDouble(QStringLiteral("prop_ring_cyl_z"), z);
        row.setDouble(QStringLiteral("prop_ring_cyl_rho"), rho);
        row.setDouble(QStringLiteral("prop_ring_angle_to_normal"),
                      disp.norm() > 1.0e-12
                          ? std::acos(std::clamp(std::abs(z) / disp.norm(), 0.0, 1.0))
                          : kNaN);
        const model::QtRing& ring = p.ring(static_cast<std::size_t>(nearestRing));
        if (!ring.atomIndices.empty() && rho > 1.0e-12) {
            const model::Vec3 vertex0 =
                body.run.conformation->atomPosition(frame,
                    static_cast<std::size_t>(ring.atomIndices.front()));
            const model::Vec3 ref = vertex0 - g.center;
            const model::Vec3 refPlane = ref - ref.dot(g.normal) * g.normal;
            const double refNorm = refPlane.norm();
            if (refNorm > 1.0e-12) {
                const model::Vec3 dHat = plane / rho;
                const model::Vec3 rHat = refPlane / refNorm;
                row.setDouble(QStringLiteral("prop_ring_cos_phi"), dHat.dot(rHat));
                row.setDouble(QStringLiteral("prop_ring_sin_phi"), dHat.cross(rHat).dot(g.normal));
            }
        }
    }

    const bool dft = targetPresent(body.run, atom, frame);
    row.setBool(QStringLiteral("dft_present"), dft);
    row.setBool(QStringLiteral("target_total_present"), dft);
    row.setBool(QStringLiteral("target_dia_present"), dft);
    row.setBool(QStringLiteral("target_para_present"), dft);
    row.setBool(QStringLiteral("aimnet2_embedding_present"),
                aimnet2EmbeddingPresent(body, atom, frame));
}

}  // namespace
}  // namespace h5reader::rediscover

namespace h5reader::rediscover {
namespace {

QString routeName(FieldRoute route) {
    switch (route) {
    case FieldRoute::Measure: return QStringLiteral("measure");
    case FieldRoute::NativePair: return QStringLiteral("native_pair");
    case FieldRoute::Incidence: return QStringLiteral("incidence_reduction");
    case FieldRoute::Target: return QStringLiteral("target");
    case FieldRoute::Topology: return QStringLiteral("topology");
    }
    return QStringLiteral("unknown");
}

FieldRoute routeFor(const io::FieldSpec& spec) {
    if (isTargetField(spec)) return FieldRoute::Target;
    if (isStructuredField(spec)) return FieldRoute::Topology;
    if (isRingContributionPair(spec)) return FieldRoute::NativePair;
    if (isIncidenceAxis(spec)) return FieldRoute::Incidence;
    return FieldRoute::Measure;
}

bool openCsv(QSaveFile& file, QTextStream*& stream, const QString& header, QString* err_out) {
    QDir().mkpath(QFileInfo(file.fileName()).dir().absolutePath());
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (err_out) *err_out = QStringLiteral("cannot open CSV: %1").arg(file.fileName());
        return false;
    }
    stream = new QTextStream(&file);
    *stream << header << '\n';
    return true;
}

std::vector<std::unique_ptr<FieldSink>>
openFieldSinks(const QString& outDir, std::size_t nRows, QString* err_out) {
    std::vector<std::unique_ptr<FieldSink>> sinks;
    const auto& scoped = ScopedProducerCatalog();
    sinks.reserve(scoped.size());
    for (const io::FieldSpec* spec : scoped) {
        auto sink = std::make_unique<FieldSink>();
        sink->spec = spec;
        sink->route = routeFor(*spec);
        sink->cols = componentCountFor(*spec);
        sink->support.spec = spec;
        sink->support.route = routeName(sink->route);
        sink->support.componentLayout = componentLayoutFor(*spec);
        sink->support.dtype = QStringLiteral("float64");
        sink->support.componentPopulated.assign(sink->cols, 0);

        const QString stem = fieldStem(*spec);
        switch (sink->route) {
        case FieldRoute::Measure: {
            sink->support.sidecar = QStringLiteral("measures/%1.npy").arg(stem);
            sink->support.shape = sink->cols == 1
                                      ? QStringLiteral("[%1]").arg(nRows)
                                      : QStringLiteral("[%1,%2]").arg(nRows).arg(sink->cols);
            sink->values = std::make_unique<RawArrayWriter>(
                QDir(outDir).filePath(QStringLiteral("measures/%1.npy").arg(stem)),
                QByteArray("<f8"), sizeof(double));
            sink->presence = std::make_unique<RawArrayWriter>(
                QDir(outDir).filePath(QStringLiteral("measures/%1_presence.npy").arg(stem)),
                QByteArray("|u1"), sizeof(std::uint8_t));
            sink->absenceFile = std::make_unique<QSaveFile>(
                QDir(outDir).filePath(QStringLiteral("measures/%1_absence_reason.csv").arg(stem)));
            break;
        }
        case FieldRoute::Incidence: {
            sink->support.sidecar = QStringLiteral("measures/%1_by_atom.npy").arg(stem);
            sink->support.shape = QStringLiteral("[%1,4,%2]").arg(nRows).arg(sink->cols);
            sink->byAtom = std::make_unique<RawArrayWriter>(
                QDir(outDir).filePath(QStringLiteral("measures/%1_by_atom.npy").arg(stem)),
                QByteArray("<f8"), sizeof(double));
            sink->presence = std::make_unique<RawArrayWriter>(
                QDir(outDir).filePath(QStringLiteral("measures/%1_presence.npy").arg(stem)),
                QByteArray("|u1"), sizeof(std::uint8_t));
            sink->nativeValues = std::make_unique<RawArrayWriter>(
                QDir(outDir).filePath(QStringLiteral("native/%1.npy").arg(stem)),
                QByteArray("<f8"), sizeof(double));
            sink->absenceFile = std::make_unique<QSaveFile>(
                QDir(outDir).filePath(QStringLiteral("measures/%1_absence_reason.csv").arg(stem)));
            sink->nativeMapFile = std::make_unique<QSaveFile>(
                QDir(outDir).filePath(QStringLiteral("native/%1_map.csv").arg(stem)));
            break;
        }
        case FieldRoute::NativePair: {
            sink->support.sidecar = QStringLiteral("native/%1.npy").arg(stem);
            sink->support.shape = QStringLiteral("[native,%1]").arg(sink->cols);
            sink->nativeValues = std::make_unique<RawArrayWriter>(
                QDir(outDir).filePath(QStringLiteral("native/%1.npy").arg(stem)),
                QByteArray("<f8"), sizeof(double));
            sink->nativeMapFile = std::make_unique<QSaveFile>(
                QDir(outDir).filePath(QStringLiteral("native/%1_map.csv").arg(stem)));
            break;
        }
        case FieldRoute::Target: {
            sink->support.sidecar = QStringLiteral("targets/%1_*.npy").arg(stem);
            sink->support.shape = QStringLiteral("[%1,{T0,T1,T2,raw}]").arg(nRows);
            break;
        }
        case FieldRoute::Topology: {
            sink->support.dtype = QStringLiteral("int32/csv");
            sink->support.sidecar = QStringLiteral("topology/%1.csv").arg(stem);
            sink->support.shape = QStringLiteral("[native,*]");
            break;
        }
        }

        if (sink->values && !sink->values->open(err_out)) return {};
        if (sink->presence && !sink->presence->open(err_out)) return {};
        if (sink->byAtom && !sink->byAtom->open(err_out)) return {};
        if (sink->nativeValues && !sink->nativeValues->open(err_out)) return {};
        if (sink->absenceFile) {
            QTextStream* raw = nullptr;
            if (!openCsv(*sink->absenceFile, raw, QStringLiteral("row_id,reason,count"), err_out))
                return {};
            sink->absenceOut.reset(raw);
        }
        if (sink->nativeMapFile) {
            QTextStream* raw = nullptr;
            const QString header =
                sink->route == FieldRoute::NativePair
                    ? QStringLiteral("native_row,dataset_id,protein_id,pose_id,frame_slot,atom_index,ring_id,ring_type,role")
                    : QStringLiteral("native_row,dataset_id,protein_id,pose_id,frame_slot,atom_index,row_id,source_kind,source_id,source_role");
            if (!openCsv(*sink->nativeMapFile, raw, header, err_out)) return {};
            sink->nativeMapOut.reset(raw);
        }
        sinks.push_back(std::move(sink));
    }
    return sinks;
}

TargetSinks::TargetSinks(const QString& outDir)
    : totalT0(QDir(outDir).filePath(QStringLiteral("targets/orca_total_T0.npy")), QByteArray("<f8"), sizeof(double)),
      totalT1(QDir(outDir).filePath(QStringLiteral("targets/orca_total_T1.npy")), QByteArray("<f8"), sizeof(double)),
      totalT2(QDir(outDir).filePath(QStringLiteral("targets/orca_total_T2.npy")), QByteArray("<f8"), sizeof(double)),
      totalRaw(QDir(outDir).filePath(QStringLiteral("targets/orca_total_raw.npy")), QByteArray("<f8"), sizeof(double)),
      diaT0(QDir(outDir).filePath(QStringLiteral("targets/orca_diamagnetic_T0.npy")), QByteArray("<f8"), sizeof(double)),
      diaT1(QDir(outDir).filePath(QStringLiteral("targets/orca_diamagnetic_T1.npy")), QByteArray("<f8"), sizeof(double)),
      diaT2(QDir(outDir).filePath(QStringLiteral("targets/orca_diamagnetic_T2.npy")), QByteArray("<f8"), sizeof(double)),
      diaRaw(QDir(outDir).filePath(QStringLiteral("targets/orca_diamagnetic_raw.npy")), QByteArray("<f8"), sizeof(double)),
      paraT0(QDir(outDir).filePath(QStringLiteral("targets/orca_paramagnetic_T0.npy")), QByteArray("<f8"), sizeof(double)),
      paraT1(QDir(outDir).filePath(QStringLiteral("targets/orca_paramagnetic_T1.npy")), QByteArray("<f8"), sizeof(double)),
      paraT2(QDir(outDir).filePath(QStringLiteral("targets/orca_paramagnetic_T2.npy")), QByteArray("<f8"), sizeof(double)),
      paraRaw(QDir(outDir).filePath(QStringLiteral("targets/orca_paramagnetic_raw.npy")), QByteArray("<f8"), sizeof(double)),
      present(QDir(outDir).filePath(QStringLiteral("targets/target_present.npy")), QByteArray("|u1"), sizeof(std::uint8_t)) {}

bool TargetSinks::open(QString* err_out) {
    return totalT0.open(err_out) && totalT1.open(err_out) && totalT2.open(err_out)
           && totalRaw.open(err_out) && diaT0.open(err_out) && diaT1.open(err_out)
           && diaT2.open(err_out) && diaRaw.open(err_out) && paraT0.open(err_out)
           && paraT1.open(err_out) && paraT2.open(err_out) && paraRaw.open(err_out)
           && present.open(err_out);
}

bool TargetSinks::commit(std::size_t nRows, QString* err_out) {
    return totalT0.commit({nRows}, err_out) && totalT1.commit({nRows, 3}, err_out)
           && totalT2.commit({nRows, 5}, err_out)
           && totalRaw.commit({nRows, 3, 3}, err_out)
           && diaT0.commit({nRows}, err_out) && diaT1.commit({nRows, 3}, err_out)
           && diaT2.commit({nRows, 5}, err_out)
           && diaRaw.commit({nRows, 3, 3}, err_out)
           && paraT0.commit({nRows}, err_out) && paraT1.commit({nRows, 3}, err_out)
           && paraT2.commit({nRows, 5}, err_out)
           && paraRaw.commit({nRows, 3, 3}, err_out)
           && present.commit({nRows}, err_out);
}

E3nnSinks::E3nnSinks(const QString& outDir)
    : conformationsFile(QDir(outDir).filePath(QStringLiteral("e3nn/conformations.csv"))),
      graphOffsets(QDir(outDir).filePath(QStringLiteral("e3nn/graph_offsets.npy")), QByteArray("<i8"), sizeof(std::int64_t)),
      atomIndex(QDir(outDir).filePath(QStringLiteral("e3nn/atom_index.npy")), QByteArray("<i4"), sizeof(std::int32_t)),
      z(QDir(outDir).filePath(QStringLiteral("e3nn/Z.npy")), QByteArray("<i2"), sizeof(std::int16_t)),
      positions(QDir(outDir).filePath(QStringLiteral("e3nn/positions.npy")), QByteArray("<f8"), sizeof(double)),
      aimnet2Embedding(QDir(outDir).filePath(QStringLiteral("e3nn/aimnet2_embedding.npy")), QByteArray("<f4"), sizeof(float)),
      rowToGraphAtom(QDir(outDir).filePath(QStringLiteral("e3nn/target_row_to_graph_atom.npy")), QByteArray("<i8"), sizeof(std::int64_t)),
      sigmaTotalRaw(QDir(outDir).filePath(QStringLiteral("e3nn/sigma_total_raw.npy")), QByteArray("<f8"), sizeof(double)),
      sigmaDiaRaw(QDir(outDir).filePath(QStringLiteral("e3nn/sigma_diamagnetic_raw.npy")), QByteArray("<f8"), sizeof(double)),
      sigmaParaRaw(QDir(outDir).filePath(QStringLiteral("e3nn/sigma_paramagnetic_raw.npy")), QByteArray("<f8"), sizeof(double)) {}

bool E3nnSinks::open(QString* err_out) {
    QTextStream* raw = nullptr;
    if (!openCsv(conformationsFile, raw,
                 QStringLiteral("graph_id,dataset_id,protein_id,pose_id,frame_slot,h5_row,original_frame_index,time_ps"),
                 err_out)) {
        return false;
    }
    conformationsOut.reset(raw);
    return graphOffsets.open(err_out) && atomIndex.open(err_out) && z.open(err_out)
           && positions.open(err_out) && aimnet2Embedding.open(err_out)
           && rowToGraphAtom.open(err_out) && sigmaTotalRaw.open(err_out)
           && sigmaDiaRaw.open(err_out) && sigmaParaRaw.open(err_out);
}

bool E3nnSinks::commit(std::size_t nRows, QString* err_out) {
    if (conformationsOut) conformationsOut->flush();
    if (!conformationsFile.commit()) {
        if (err_out) *err_out = QStringLiteral("cannot commit e3nn/conformations.csv");
        return false;
    }
    return graphOffsets.commit({graphCount + 1}, err_out)
           && atomIndex.commit({graphAtoms}, err_out)
           && z.commit({graphAtoms}, err_out)
           && positions.commit({graphAtoms, 3}, err_out)
           && aimnet2Embedding.commit({graphAtoms, 256}, err_out)
           && rowToGraphAtom.commit({nRows, 2}, err_out)
           && sigmaTotalRaw.commit({nRows, 3, 3}, err_out)
           && sigmaDiaRaw.commit({nRows, 3, 3}, err_out)
           && sigmaParaRaw.commit({nRows, 3, 3}, err_out);
}

}  // namespace
}  // namespace h5reader::rediscover

namespace h5reader::rediscover {
namespace {

RegionApplySpec parseRegionSpec(const ConditioningSpec& spec) {
    RegionApplySpec out;
    if (spec.sourcePath.isEmpty()) return out;
    out.supplied = true;
    out.sourcePath = spec.sourcePath;
    out.sha256 = spec.configHash();
    out.regionDefId =
        spec.root.value(QStringLiteral("region_def_id")).toString(
            QStringLiteral("fold:unknown:%1").arg(out.sha256.left(12)));
    const QJsonArray basins =
        spec.root.value(QStringLiteral("rama")).toObject()
            .value(QStringLiteral("basins")).toArray();
    for (const QJsonValue& v : basins) {
        const QJsonObject o = v.toObject();
        RegionBasin b;
        b.id = o.value(QStringLiteral("id")).toString();
        b.label = o.value(QStringLiteral("label")).toString(b.id);
        b.phi = o.value(QStringLiteral("center_phi")).toDouble(kNaN);
        b.psi = o.value(QStringLiteral("center_psi")).toDouble(kNaN);
        if (!b.id.isEmpty() && std::isfinite(b.phi) && std::isfinite(b.psi))
            out.ramaBasins.push_back(std::move(b));
    }
    return out;
}

std::optional<RegionBasin> nearestRamaBasin(const RegionApplySpec& spec,
                                            double phi,
                                            double psi) {
    if (!spec.supplied || spec.ramaBasins.empty()
        || !std::isfinite(phi) || !std::isfinite(psi)) {
        return std::nullopt;
    }
    const RegionBasin* best = nullptr;
    double bestD = std::numeric_limits<double>::infinity();
    for (const RegionBasin& b : spec.ramaBasins) {
        const double dPhi = AngularDistance(phi, b.phi);
        const double dPsi = AngularDistance(psi, b.psi);
        const double d = dPhi * dPhi + dPsi * dPsi;
        if (d < bestD) {
            bestD = d;
            best = &b;
        }
    }
    return best ? std::optional<RegionBasin>(*best) : std::nullopt;
}

std::vector<int> configuredCutoffs(const ConditioningSpec& spec) {
    std::vector<int> out;
    const QJsonArray arr = spec.root.value(QStringLiteral("spatial_cutoffs_A")).toArray();
    for (const QJsonValue& v : arr) {
        const int x = v.toInt(-1);
        if (x > 0) out.push_back(x);
    }
    if (out.empty()) throw std::runtime_error("conditioning spec missing spatial_cutoffs_A");
    return out;
}

bool validResidue(const model::QtProtein& p, int32_t residueIndex) {
    return residueIndex >= 0 && static_cast<std::size_t>(residueIndex) < p.residueCount();
}

const model::QtResidue* residueForAtom(const model::QtProtein& p, const model::QtAtom& a) {
    return validResidue(p, a.residueIndex)
               ? &p.residue(static_cast<std::size_t>(a.residueIndex))
               : nullptr;
}

QString residueTypeLabel(const model::QtProtein& p, const model::QtResidue* r) {
    if (!r) return {};
    return p.residueLabel(static_cast<std::size_t>(r->residueIndex),
                          model::NamingConvention::Iupac,
                          model::NamingSource::Derived);
}

QString residueTypeLabelForAa(model::AminoAcid aa) {
    return QString::fromLatin1(model::IupacResidue3LetterFor(aa));
}

std::set<int32_t> bondedAtoms(const model::QtProtein& p, std::size_t atom) {
    std::set<int32_t> out;
    const model::QtTopology& topo = p.topology();
    for (int32_t bi : topo.bondIndicesForAtom(atom)) {
        if (bi < 0) continue;
        const model::QtBond& b = topo.bondAt(static_cast<std::size_t>(bi));
        if (b.atomIndexA == static_cast<int32_t>(atom) && b.atomIndexB >= 0)
            out.insert(b.atomIndexB);
        if (b.atomIndexB == static_cast<int32_t>(atom) && b.atomIndexA >= 0)
            out.insert(b.atomIndexA);
    }
    return out;
}

bool bondTouchesAtom(const model::QtBond& b, std::size_t atom) {
    return b.atomIndexA == static_cast<int32_t>(atom)
           || b.atomIndexB == static_cast<int32_t>(atom);
}

bool atomInRing(const model::QtProtein& p, std::size_t atom, std::size_t ring) {
    for (int32_t mi : p.topology().ringMembershipsForAtom(atom)) {
        if (mi < 0) continue;
        const model::QtRingMembership& m =
            p.topology().ringMembershipAt(static_cast<std::size_t>(mi));
        if (m.ringId == static_cast<int32_t>(ring)) return true;
    }
    return false;
}

bool excludedSource(const Body& body,
                    CloudKind kind,
                    const SourceRef& ref,
                    std::size_t atom,
                    const std::set<int32_t>& bonded) {
    const model::QtProtein& p = *body.run.protein;
    switch (kind) {
    case CloudKind::Atoms:
    case CloudKind::ChargeSites:
        return ref.entity_index == static_cast<int32_t>(atom)
               || bonded.count(ref.entity_index) != 0;
    case CloudKind::AllBondMidpoints:
    case CloudKind::BondMidpoints:
        if (ref.entity_index < 0
            || static_cast<std::size_t>(ref.entity_index) >= p.topology().bondCount()) {
            return true;
        }
        return bondTouchesAtom(p.topology().bondAt(static_cast<std::size_t>(ref.entity_index)), atom);
    case CloudKind::RingCenters:
        if (ref.entity_index < 0
            || static_cast<std::size_t>(ref.entity_index) >= p.ringCount()) {
            return true;
        }
        return atomInRing(p, atom, static_cast<std::size_t>(ref.entity_index));
    }
    return false;
}

std::size_t countNear(const Body& body,
                      CloudKind kind,
                      std::size_t atom,
                      std::size_t frame,
                      double cutoff,
                      const std::set<int32_t>& bonded) {
    const model::Vec3 query = body.run.conformation->atomPosition(frame, atom);
    std::size_t out = 0;
    for (const SourceRef& ref : body.idx.spatial.near(kind, frame, query, cutoff)) {
        if (excludedSource(body, kind, ref, atom, bonded)) continue;
        ++out;
    }
    return out;
}

double nearestDist(const Body& body,
                   CloudKind kind,
                   std::size_t atom,
                   std::size_t frame,
                   const std::set<int32_t>& bonded,
                   int32_t* entityOut) {
    const model::Vec3 query = body.run.conformation->atomPosition(frame, atom);
    double best = kNaN;
    int32_t bestEntity = -1;
    for (const SourceRef& ref : body.idx.spatial.near(kind, frame, query, 99.0)) {
        if (excludedSource(body, kind, ref, atom, bonded)) continue;
        const model::Vec3 q =
            body.idx.spatial.tree(kind, frame).pointAt(static_cast<std::size_t>(ref.cloud_index));
        const double d = (q - query).norm();
        if (!std::isfinite(best) || d < best) {
            best = d;
            bestEntity = ref.entity_index;
        }
    }
    if (entityOut) *entityOut = bestEntity;
    return best;
}

QString ringMembershipId(const model::QtProtein& p, std::size_t atom) {
    QStringList ids;
    for (int32_t mi : p.topology().ringMembershipsForAtom(atom))
        ids << QString::number(mi);
    return ids.join(QLatin1Char('|'));
}

QString ringTypeForAtom(const model::QtProtein& p, std::size_t atom) {
    QStringList types;
    for (int32_t mi : p.topology().ringMembershipsForAtom(atom)) {
        if (mi < 0) continue;
        const model::QtRingMembership& m =
            p.topology().ringMembershipAt(static_cast<std::size_t>(mi));
        if (m.ringId >= 0 && static_cast<std::size_t>(m.ringId) < p.ringCount())
            types << QString::fromLatin1(p.ring(static_cast<std::size_t>(m.ringId)).TypeName());
    }
    types.removeDuplicates();
    return types.join(QLatin1Char('|'));
}

std::optional<std::size_t> uniqueRingNativeRow(const model::QtProtein& p,
                                               std::size_t atom,
                                               io::NativeAxis axis,
                                               QString* reason_out) {
    std::vector<std::size_t> rows;
    for (int32_t mi : p.topology().ringMembershipsForAtom(atom)) {
        if (mi < 0) continue;
        const model::QtRingMembership& m =
            p.topology().ringMembershipAt(static_cast<std::size_t>(mi));
        if (m.ringId < 0 || static_cast<std::size_t>(m.ringId) >= p.ringCount()) continue;
        const model::QtRing& ring = p.ring(static_cast<std::size_t>(m.ringId));
        switch (axis) {
        case io::NativeAxis::AromaticRing:
            if (ring.IsAromatic() && ring.nativeAxisIndex >= 0)
                rows.push_back(static_cast<std::size_t>(ring.nativeAxisIndex));
            break;
        case io::NativeAxis::SaturatedRing:
            if (!ring.IsAromatic() && ring.nativeAxisIndex >= 0)
                rows.push_back(static_cast<std::size_t>(ring.nativeAxisIndex));
            break;
        case io::NativeAxis::Ring:
            rows.push_back(static_cast<std::size_t>(ring.ringId));
            break;
        case io::NativeAxis::RingMembership:
            rows.push_back(static_cast<std::size_t>(mi));
            break;
        default:
            break;
        }
    }
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    if (rows.size() == 1) return rows.front();
    if (reason_out) {
        *reason_out = rows.empty() ? QStringLiteral("no-ring-membership")
                                   : QStringLiteral("non-unique-ring-membership");
    }
    return std::nullopt;
}

std::optional<std::size_t> nativeRowForAtom(const Body& body,
                                            const io::FieldSpec& spec,
                                            std::size_t atom,
                                            QString* reason_out) {
    const model::QtProtein& p = *body.run.protein;
    const model::QtAtom& a = p.atom(atom);
    switch (spec.axis) {
    case io::NativeAxis::Atom:
        return atom;
    case io::NativeAxis::Residue:
        if (a.residueIndex >= 0) return static_cast<std::size_t>(a.residueIndex);
        if (reason_out) *reason_out = QStringLiteral("atom-without-residue");
        return std::nullopt;
    case io::NativeAxis::Protein:
        return 0;
    case io::NativeAxis::AromaticRing:
    case io::NativeAxis::SaturatedRing:
    case io::NativeAxis::Ring:
    case io::NativeAxis::RingMembership:
        return uniqueRingNativeRow(p, atom, spec.axis, reason_out);
    default:
        if (reason_out) *reason_out = QStringLiteral("native-axis-not-row-broadcast");
        return std::nullopt;
    }
}

std::array<double, 9> matArray(const model::Mat3& m) {
    return {m(0, 0), m(0, 1), m(0, 2),
            m(1, 0), m(1, 1), m(1, 2),
            m(2, 0), m(2, 1), m(2, 2)};
}

void writeTensorTarget(RawArrayWriter& raw,
                       RawArrayWriter& t0,
                       RawArrayWriter& t1,
                       RawArrayWriter& t2,
                       const model::Mat3& m,
                       bool present) {
    const std::array<double, 9> r = present ? matArray(m)
                                            : std::array<double, 9>{kNaN, kNaN, kNaN,
                                                                    kNaN, kNaN, kNaN,
                                                                    kNaN, kNaN, kNaN};
    raw.writeMany(r.data(), r.size());
    if (!present) {
        const double nan = kNaN;
        const std::array<double, 3> nan3{kNaN, kNaN, kNaN};
        const std::array<double, 5> nan5{kNaN, kNaN, kNaN, kNaN, kNaN};
        t0.writeOne(nan);
        t1.writeMany(nan3.data(), nan3.size());
        t2.writeMany(nan5.data(), nan5.size());
        return;
    }
    const model::SphericalTensor st = DecomposeLibrary(m);
    t0.writeOne(st.T0);
    t1.writeMany(st.T1.data(), st.T1.size());
    t2.writeMany(st.T2.data(), st.T2.size());
}

std::optional<std::pair<int32_t, int32_t>>
nativeEndpoints(const Body& body,
                const io::FieldSpec& spec,
                std::size_t nativeRow,
                std::size_t frame) {
    const StaticNpyArray* a = body.run.producerArray(spec.kind);
    if (!a || nativeRow >= a->rowsForFrame(frame)) return std::nullopt;
    const std::size_t row = a->rowFor(nativeRow, frame);
    auto at = [&](std::size_t col) -> double {
        if (col >= a->cols) return kNaN;
        if (!a->values.empty()) return a->values[row * a->cols + col];
        if (!a->floatValues.empty()) return a->floatValues[row * a->cols + col];
        return kNaN;
    };
    switch (spec.kind) {
    case io::FieldKind::MOPACBondOrders:
        return std::make_pair(static_cast<int32_t>(std::llround(at(0))),
                              static_cast<int32_t>(std::llround(at(1))));
    case io::FieldKind::MOPACTopologyBondOrdersFull:
        return std::make_pair(static_cast<int32_t>(std::llround(at(1))),
                              static_cast<int32_t>(std::llround(at(2))));
    case io::FieldKind::MOPACBondNeighbors:
        return std::make_pair(static_cast<int32_t>(std::llround(at(0))),
                              static_cast<int32_t>(std::llround(at(1))));
    case io::FieldKind::MOPACBondOrdersUnique:
        return std::make_pair(static_cast<int32_t>(std::llround(at(0))),
                              static_cast<int32_t>(std::llround(at(1))));
    default:
        break;
    }
    if (spec.axis == io::NativeAxis::Bond) {
        if (nativeRow >= body.run.protein->topology().bondCount()) return std::nullopt;
        const model::QtBond& b = body.run.protein->topology().bondAt(nativeRow);
        return std::make_pair(b.atomIndexA, b.atomIndexB);
    }
    return std::nullopt;
}

}  // namespace
}  // namespace h5reader::rediscover
