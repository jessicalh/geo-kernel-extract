#include "ReaderAlignmentExporter.h"

#include "CalcsetManifest.h"
#include "QtManifest.h"
#include "QtNpyReader.h"
#include "QtNpyWriter.h"
#include "QtProteinLoader.h"

#include "../model/QtProtein.h"
#include "../model/TrajectoryConformation.h"

#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonValue>
#include <QLoggingCategory>
#include <QSaveFile>
#include <QStringList>
#include <QUuid>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#ifndef H5READER_GIT_COMMIT
#define H5READER_GIT_COMMIT "unknown"
#endif

#ifndef H5READER_GIT_DIRTY
#define H5READER_GIT_DIRTY 1
#endif

namespace h5reader::io {

namespace {
Q_LOGGING_CATEGORY(cAlignmentExport, "h5reader.alignment.export")

constexpr const char* kReaderBaseCommit =
    "c51e9e382bac61a8557251dc6bbaf2dc82255be1";
constexpr const char* kProducerCatalogSha256 =
    "5238d7fff451afe6a816118369e12cbd322361491d7915d38344d558f65f64fc";
constexpr const char* kExtractorCommit =
    "91f8294d5b1c18030cec906e5ec8c08400fd8310";
constexpr const char* kRecipeId =
    "08de7358b2179445241f8e6ba18ecc04803309fe67bb98ad16575beb3453c117";
constexpr const char* kSourceProvenanceAuthority =
    "READER_ALIGNMENT_EXPORT_CONTRACT.md";

constexpr double kSourcePositionAbsoluteToleranceAngstrom = 1e-10;
constexpr double kSourcePositionRelativeTolerance = 1e-12;
constexpr std::size_t kExpectedFrameCount = 100;
constexpr std::uint64_t kFirstOriginalFrame = 15;
constexpr std::uint64_t kOriginalFrameStride = 15;
constexpr double kFirstTimePicoseconds = 150.0;
constexpr double kTimeStridePicoseconds = 150.0;
constexpr double kTimeAxisTolerancePicoseconds = 1e-9;

struct ExportedFile {
    QString name;
    QString dtype;
    std::vector<std::size_t> shape;
    bool byteCopy = false;
};

struct PositionSourceValidation {
    bool ok = false;
    QString error;
    double maxCoordinateDifferenceAngstrom = 0.0;
    double maxPointDifferenceAngstrom = 0.0;
    QStringList paths;
};

struct TransformValidation {
    bool ok = false;
    QString error;
    double maxOrthogonalityError = 0.0;
    double minDeterminant = std::numeric_limits<double>::infinity();
    double maxDeterminant = -std::numeric_limits<double>::infinity();
    double maxForwardReconstructionErrorAngstrom = 0.0;
    double maxInverseReconstructionErrorAngstrom = 0.0;
    double maxRoundTripReconstructionErrorAngstrom = 0.0;
};

class StagingDirectoryGuard final {
public:
    explicit StagingDirectoryGuard(QString path) : path_(std::move(path)) {}
    ~StagingDirectoryGuard() {
        if (active_ && !path_.isEmpty())
            QDir(path_).removeRecursively();
    }
    void release() { active_ = false; }

private:
    QString path_;
    bool active_ = true;
};

bool CancellationRequested(const ReaderAlignmentExportControl* control,
                           QString* error = nullptr) {
    if (!control || !control->cancelRequested())
        return false;
    if (error)
        *error = QStringLiteral("scientific alignment export cancelled");
    return true;
}

QString CleanAbsolutePath(const QString& path) {
    if (path.isEmpty())
        return {};

    QFileInfo info(path);
    if (info.exists()) {
        const QString canonical = info.canonicalFilePath();
        if (!canonical.isEmpty())
            return QDir::fromNativeSeparators(QDir::cleanPath(canonical));
    }

    QString unresolved = info.fileName();
    QFileInfo ancestor(info.absolutePath());
    while (!ancestor.exists()) {
        const QString name = ancestor.fileName();
        if (!name.isEmpty())
            unresolved = name + QLatin1Char('/') + unresolved;
        const QString parentPath = ancestor.absolutePath();
        if (parentPath == ancestor.filePath())
            break;
        ancestor = QFileInfo(parentPath);
    }
    const QString canonicalAncestor = ancestor.canonicalFilePath();
    const QString base = canonicalAncestor.isEmpty()
        ? ancestor.absoluteFilePath() : canonicalAncestor;
    return QDir::fromNativeSeparators(
        QDir::cleanPath(QDir(base).filePath(unresolved)));
}

Qt::CaseSensitivity PathCaseSensitivity() {
#ifdef Q_OS_WIN
    return Qt::CaseInsensitive;
#else
    return Qt::CaseSensitive;
#endif
}

bool SameOrChildPath(const QString& candidate, const QString& parent) {
    const QString cleanCandidate = CleanAbsolutePath(candidate);
    QString cleanParent = CleanAbsolutePath(parent);
    if (cleanCandidate.isEmpty() || cleanParent.isEmpty())
        return false;
    while (cleanParent.endsWith(QLatin1Char('/')) && cleanParent.size() > 1)
        cleanParent.chop(1);
    if (cleanCandidate.compare(cleanParent, PathCaseSensitivity()) == 0)
        return true;
    return cleanCandidate.startsWith(cleanParent + QLatin1Char('/'),
                                     PathCaseSensitivity());
}

bool PathsOverlap(const QString& left, const QString& right) {
    return SameOrChildPath(left, right) || SameOrChildPath(right, left);
}

bool NumericDtypeMatches(const std::string& descriptor,
                         const QString& expectedDtype) {
    const std::string expected = expectedDtype.toStdString();
    return descriptor == "'" + expected + "'"
        || descriptor == "\"" + expected + "\"";
}

bool ValidateRetainedFrameAxis(
    const std::vector<std::uint64_t>& originalFrames,
    const std::vector<double>& timesPicoseconds,
    QString* error) {
    if (originalFrames.size() != kExpectedFrameCount
        || timesPicoseconds.size() != kExpectedFrameCount) {
        *error = QStringLiteral(
            "retained frame axis must contain exactly %1 rows")
            .arg(static_cast<qulonglong>(kExpectedFrameCount));
        return false;
    }
    for (std::size_t row = 0; row < kExpectedFrameCount; ++row) {
        const std::uint64_t expectedFrame =
            kFirstOriginalFrame + row * kOriginalFrameStride;
        const double expectedTime =
            kFirstTimePicoseconds
            + static_cast<double>(row) * kTimeStridePicoseconds;
        if (originalFrames[row] != expectedFrame) {
            *error = QStringLiteral(
                "original frame axis mismatch at row %1: expected %2, got %3")
                .arg(static_cast<qulonglong>(row))
                .arg(static_cast<qulonglong>(expectedFrame))
                .arg(static_cast<qulonglong>(originalFrames[row]));
            return false;
        }
        if (!std::isfinite(timesPicoseconds[row])
            || std::abs(timesPicoseconds[row] - expectedTime)
                > kTimeAxisTolerancePicoseconds) {
            *error = QStringLiteral(
                "time axis mismatch at row %1: expected %2 ps, got %3 ps")
                .arg(static_cast<qulonglong>(row))
                .arg(expectedTime, 0, 'g', 17)
                .arg(timesPicoseconds[row], 0, 'g', 17);
            return false;
        }
    }
    return true;
}

bool ReadJsonObject(const QString& path, QJsonObject* object, QString* error) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        *error = QStringLiteral("could not open JSON %1: %2")
            .arg(path, file.errorString());
        return false;
    }
    QJsonParseError parseError{};
    const QJsonDocument document = QJsonDocument::fromJson(file.readAll(), &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject()) {
        *error = QStringLiteral("could not parse JSON %1: %2")
            .arg(path, parseError.errorString());
        return false;
    }
    *object = document.object();
    return true;
}

bool WriteJsonObject(const QString& path,
                     const QJsonObject& object,
                     QString* error) {
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
        *error = QStringLiteral("could not open JSON output %1: %2")
            .arg(path, file.errorString());
        return false;
    }
    const QByteArray bytes = QJsonDocument(object).toJson(QJsonDocument::Indented);
    if (file.write(bytes) != bytes.size()) {
        *error = QStringLiteral("could not write JSON output %1: %2")
            .arg(path, file.errorString());
        file.cancelWriting();
        return false;
    }
    if (!file.commit()) {
        *error = QStringLiteral("could not commit JSON output %1: %2")
            .arg(path, file.errorString());
        return false;
    }
    return true;
}

QString Sha256File(const QString& path,
                   QString* error,
                   const ReaderAlignmentExportControl* control = nullptr) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        *error = QStringLiteral("could not hash %1: %2")
            .arg(path, file.errorString());
        return {};
    }
    QCryptographicHash hash(QCryptographicHash::Sha256);
    constexpr qint64 blockSize = qint64{1024} * 1024;
    while (!file.atEnd()) {
        if (CancellationRequested(control, error))
            return {};
        const QByteArray block = file.read(blockSize);
        if (block.isEmpty() && file.error() != QFileDevice::NoError) {
            *error = QStringLiteral("could not read %1 while hashing: %2")
                .arg(path, file.errorString());
            return {};
        }
        hash.addData(block);
    }
    return QString::fromLatin1(hash.result().toHex());
}

QJsonArray ShapeToJson(const std::vector<std::size_t>& shape) {
    QJsonArray array;
    for (std::size_t dimension : shape)
        array.append(static_cast<qint64>(dimension));
    return array;
}

std::vector<std::size_t> ShapeFromJson(const QJsonValue& value,
                                       bool* ok) {
    std::vector<std::size_t> shape;
    *ok = false;
    if (!value.isArray())
        return shape;
    const QJsonArray dimensions = value.toArray();
    for (const auto dimension : dimensions) {
        if (!dimension.isDouble())
            return {};
        const qint64 raw = dimension.toInteger(-1);
        if (raw < 0)
            return {};
        shape.push_back(static_cast<std::size_t>(raw));
    }
    *ok = !shape.empty();
    return shape;
}

QJsonObject StatusMeanings() {
    QJsonObject meanings;
    for (int code = static_cast<int>(model::ScientificFitStatus::Valid);
         code <= static_cast<int>(model::ScientificFitStatus::NonFiniteOutput);
         ++code) {
        const auto status = static_cast<model::ScientificFitStatus>(code);
        meanings.insert(QString::number(code),
                        QString::fromLatin1(model::NameForScientificFitStatus(status)));
    }
    return meanings;
}

QJsonArray Vec3ToJson(const model::Vec3& value) {
    return QJsonArray{value.x(), value.y(), value.z()};
}

QJsonObject FitToJson(const model::ScientificAlignmentResult& alignment,
                      const QString& atomSelection) {
    return QJsonObject{
        {"atom_selection", atomSelection},
        {"fit_policy", QStringLiteral("proper Kabsch SVD; accept numerical rank >= 2")},
        {"seed_frame_row", static_cast<qint64>(alignment.policy.seedFrame)},
        {"centroid_anchor_A", Vec3ToJson(alignment.mean.centroidAnchor)},
        {"max_iterations", alignment.policy.maxIterations},
        {"convergence_tolerance_A", alignment.policy.convergenceToleranceAngstrom},
        {"converged", alignment.mean.converged},
        {"iterations", alignment.mean.iterations},
        {"final_delta_A", alignment.mean.finalDeltaAngstrom},
        {"reference_build_degeneracy_count",
         static_cast<qint64>(alignment.mean.referenceBuildDegeneracyCount)},
        {"temporal_smoothing", false},
        {"numerical_rank_rule",
         QStringLiteral("rank = count(sigma_i > max(abs_tol, rel_tol * sigma_max)); accept rank >= 2")},
        {"rank_relative_tolerance", alignment.policy.rankRelativeTolerance},
        {"rank_absolute_tolerance_A2", alignment.policy.rankAbsoluteTolerance},
        {"rotation_tolerance", alignment.policy.rotationTolerance},
    };
}

bool CopyFileExact(const QString& source,
                   const QString& destination,
                   QString* error) {
    if (!QFileInfo(source).isFile()) {
        *error = QStringLiteral("required source file is missing: %1").arg(source);
        return false;
    }
    if (!QFile::copy(source, destination)) {
        *error = QStringLiteral("could not copy %1 to %2")
            .arg(source, destination);
        return false;
    }
    return true;
}

bool AddInputHash(const QString& role,
                  const QString& path,
                  QJsonArray* inputs,
                  QString* error,
                  std::optional<std::size_t> frameRow = std::nullopt,
                  std::optional<std::uint64_t> originalFrame = std::nullopt,
                  const ReaderAlignmentExportControl* control = nullptr) {
    const QString hash = Sha256File(path, error, control);
    if (hash.isEmpty())
        return false;
    QJsonObject item{
        {"role", role},
        {"path", CleanAbsolutePath(path)},
        {"size", QFileInfo(path).size()},
        {"sha256", hash},
    };
    if (frameRow)
        item.insert(QStringLiteral("frame_row"), static_cast<qint64>(*frameRow));
    if (originalFrame)
        item.insert(QStringLiteral("original_frame_index"),
                    static_cast<qint64>(*originalFrame));
    inputs->append(item);
    return true;
}

bool NumericArrayMatches(const QString& path,
                         const QString& expectedDtype,
                         const std::vector<std::size_t>& expectedShape,
                         const std::vector<double>& expectedValues,
                         QString* error) {
    const QtNpyReader::NumericArray array =
        QtNpyReader::ReadNumericArrayWidened(path);
    if (!array.ok) {
        *error = array.error;
        return false;
    }
    if (!NumericDtypeMatches(array.descr, expectedDtype)) {
        *error = QStringLiteral("NPY dtype mismatch in %1: expected %2, got %3")
            .arg(path, expectedDtype, QString::fromStdString(array.descr));
        return false;
    }
    if (array.shape != expectedShape) {
        *error = QStringLiteral("NPY shape mismatch in %1").arg(path);
        return false;
    }
    if (array.data.size() != expectedValues.size()) {
        *error = QStringLiteral("NPY payload length mismatch in %1").arg(path);
        return false;
    }
    for (std::size_t i = 0; i < array.data.size(); ++i) {
        if (array.data[i] != expectedValues[i]) {
            *error = QStringLiteral("NPY payload mismatch in %1 at flattened index %2")
                .arg(path)
                .arg(static_cast<qulonglong>(i));
            return false;
        }
    }
    return true;
}

bool UIntArrayMatches(const QString& path,
                      const QString& expectedDtype,
                      const std::vector<std::size_t>& expectedShape,
                      const std::vector<std::uint64_t>& expectedValues,
                      QString* error) {
    std::vector<double> widened;
    widened.reserve(expectedValues.size());
    for (std::uint64_t value : expectedValues) {
        if (value > (std::uint64_t{1} << 53)) {
            *error = QStringLiteral("integer validation value exceeds exact double range");
            return false;
        }
        widened.push_back(static_cast<double>(value));
    }
    return NumericArrayMatches(path, expectedDtype, expectedShape, widened, error);
}

bool UInt8ArrayMatches(const QString& path,
                       const std::vector<std::size_t>& expectedShape,
                       const std::vector<std::uint8_t>& expectedValues,
                       QString* error) {
    std::vector<double> widened;
    widened.reserve(expectedValues.size());
    for (std::uint8_t value : expectedValues)
        widened.push_back(static_cast<double>(value));
    return NumericArrayMatches(path, QStringLiteral("|u1"), expectedShape,
                               widened, error);
}

std::vector<double> FlattenReference(
    const model::ScientificAlignmentResult& alignment) {
    std::vector<double> values;
    values.reserve(alignment.referencePositions.size() * 3);
    for (const model::Vec3& position : alignment.referencePositions) {
        values.push_back(position.x());
        values.push_back(position.y());
        values.push_back(position.z());
    }
    return values;
}

std::vector<std::uint64_t> WidenIndices(
    const std::vector<std::size_t>& indices) {
    std::vector<std::uint64_t> values;
    values.reserve(indices.size());
    for (std::size_t index : indices)
        values.push_back(static_cast<std::uint64_t>(index));
    return values;
}

PositionSourceValidation ValidatePositionSources(
    const QString& extractionDirectory,
    const model::ScientificPositionTable& positions,
    const std::vector<std::uint64_t>& originalFrames,
    const ReaderAlignmentExportControl* control = nullptr) {
    PositionSourceValidation result;
    if (originalFrames.size() != positions.frameCount) {
        result.error = QStringLiteral("H5 original frame axis is incomplete");
        return result;
    }

    result.paths.reserve(static_cast<qsizetype>(positions.frameCount));
    for (std::size_t frame = 0; frame < positions.frameCount; ++frame) {
        if (CancellationRequested(control, &result.error))
            return result;
        const QString path = QDir(extractionDirectory).filePath(
            QStringLiteral("npys/frame_%1/pos.npy")
                .arg(static_cast<qulonglong>(originalFrames[frame]),
                     6, 10, QLatin1Char('0')));
        const QtNpyReader::NumericArray source =
            QtNpyReader::ReadNumericArrayWidened(path);
        if (!source.ok) {
            result.error = source.error;
            return result;
        }
        if (!NumericDtypeMatches(source.descr, QStringLiteral("<f8"))
            || source.shape != std::vector<std::size_t>{positions.atomCount, 3}
            || source.data.size() != positions.atomCount * 3) {
            result.error = QStringLiteral(
                "frame position source has wrong dtype or shape: %1").arg(path);
            return result;
        }

        for (std::size_t atom = 0; atom < positions.atomCount; ++atom) {
            const model::Vec3& h5 = positions.at(frame, atom);
            const model::Vec3 npy(source.data[atom * 3],
                                  source.data[atom * 3 + 1],
                                  source.data[atom * 3 + 2]);
            if (!npy.allFinite()) {
                result.error = QStringLiteral("frame position source contains nonfinite data: %1")
                    .arg(path);
                return result;
            }
            const model::Vec3 difference = h5 - npy;
            result.maxPointDifferenceAngstrom =
                std::max(result.maxPointDifferenceAngstrom, difference.norm());
            for (int coordinate = 0; coordinate < 3; ++coordinate) {
                const double absoluteDifference = std::abs(difference(coordinate));
                result.maxCoordinateDifferenceAngstrom =
                    std::max(result.maxCoordinateDifferenceAngstrom,
                             absoluteDifference);
                const double scale = std::max(std::abs(h5(coordinate)),
                                              std::abs(npy(coordinate)));
                const double tolerance = kSourcePositionAbsoluteToleranceAngstrom
                    + kSourcePositionRelativeTolerance * scale;
                if (absoluteDifference > tolerance) {
                    result.error = QStringLiteral(
                        "H5 and pos.npy differ at frame %1 atom %2 coordinate %3: %4 A > %5 A")
                        .arg(static_cast<qulonglong>(frame))
                        .arg(static_cast<qulonglong>(atom))
                        .arg(coordinate)
                        .arg(absoluteDifference, 0, 'g', 17)
                        .arg(tolerance, 0, 'g', 17);
                    return result;
                }
            }
        }
        result.paths.push_back(path);
    }
    result.ok = true;
    return result;
}

TransformValidation ValidateTransforms(
    const model::ScientificPositionTable& raw,
    const model::ScientificAlignmentResult& alignment,
    const std::vector<double>* storedAligned) {
    TransformValidation result;
    if (!alignment.ok || alignment.frameFits.size() != raw.frameCount
        || (storedAligned
            && storedAligned->size() != raw.frameCount * raw.atomCount * 3)) {
        result.error = QStringLiteral("transform validation dimensions are inconsistent");
        return result;
    }

    for (std::size_t frame = 0; frame < raw.frameCount; ++frame) {
        const model::ScientificFrameFit& fit = alignment.frameFits[frame];
        if (!fit.valid()) {
            result.error = QStringLiteral("transform status is invalid at frame %1")
                .arg(static_cast<qulonglong>(frame));
            return result;
        }
        const double determinant = fit.rotation.determinant();
        result.minDeterminant = std::min(result.minDeterminant, determinant);
        result.maxDeterminant = std::max(result.maxDeterminant, determinant);
        result.maxOrthogonalityError = std::max(
            result.maxOrthogonalityError,
            (fit.rotation * fit.rotation.transpose()
                - model::Mat3::Identity()).norm());

        for (std::size_t atom = 0; atom < raw.atomCount; ++atom) {
            const std::size_t base = (frame * raw.atomCount + atom) * 3;
            const model::Vec3 forward =
                fit.rotation * raw.at(frame, atom) + fit.translation;
            const model::Vec3 roundTrip =
                fit.rotation.transpose() * (forward - fit.translation);
            result.maxRoundTripReconstructionErrorAngstrom = std::max(
                result.maxRoundTripReconstructionErrorAngstrom,
                (roundTrip - raw.at(frame, atom)).norm());
            if (storedAligned) {
                const model::Vec3 stored(
                    (*storedAligned)[base], (*storedAligned)[base + 1],
                    (*storedAligned)[base + 2]);
                const model::Vec3 inverse =
                    fit.rotation.transpose() * (stored - fit.translation);
                result.maxForwardReconstructionErrorAngstrom = std::max(
                    result.maxForwardReconstructionErrorAngstrom,
                    (stored - forward).norm());
                result.maxInverseReconstructionErrorAngstrom = std::max(
                    result.maxInverseReconstructionErrorAngstrom,
                    (inverse - raw.at(frame, atom)).norm());
            }
        }
    }
    if (std::abs(result.minDeterminant - 1.0) > alignment.policy.rotationTolerance
        || std::abs(result.maxDeterminant - 1.0) > alignment.policy.rotationTolerance
        || result.maxOrthogonalityError > alignment.policy.rotationTolerance
        || result.maxForwardReconstructionErrorAngstrom > 1e-12
        || result.maxInverseReconstructionErrorAngstrom > 1e-10
        || result.maxRoundTripReconstructionErrorAngstrom > 1e-10) {
        result.error = QStringLiteral("rotation or reconstruction validation exceeded tolerance");
        return result;
    }
    result.ok = true;
    return result;
}

bool WriteFloatFile(const QString& directory,
                    const QString& name,
                    const std::vector<std::size_t>& shape,
                    const std::vector<double>& values,
                    std::vector<ExportedFile>* files,
                    QString* error) {
    const QString path = QDir(directory).filePath(name);
    const QtNpyWriter::Result written =
        QtNpyWriter::WriteFloat64(path, shape, values);
    if (!written.ok) {
        *error = QStringLiteral("%1: %2").arg(name, written.error);
        return false;
    }
    if (!NumericArrayMatches(path, QStringLiteral("<f8"), shape, values, error))
        return false;
    files->push_back({name, QStringLiteral("<f8"), shape, false});
    return true;
}

bool WriteUIntFile(const QString& directory,
                   const QString& name,
                   const std::vector<std::size_t>& shape,
                   const std::vector<std::uint64_t>& values,
                   std::vector<ExportedFile>* files,
                   QString* error) {
    const QString path = QDir(directory).filePath(name);
    const QtNpyWriter::Result written =
        QtNpyWriter::WriteUInt64(path, shape, values);
    if (!written.ok) {
        *error = QStringLiteral("%1: %2").arg(name, written.error);
        return false;
    }
    if (!UIntArrayMatches(path, QStringLiteral("<u8"), shape, values, error))
        return false;
    files->push_back({name, QStringLiteral("<u8"), shape, false});
    return true;
}

bool WriteStatusFile(const QString& directory,
                     const QString& name,
                     const std::vector<std::size_t>& shape,
                     const std::vector<std::uint8_t>& values,
                     std::vector<ExportedFile>* files,
                     QString* error) {
    const QString path = QDir(directory).filePath(name);
    const QtNpyWriter::Result written =
        QtNpyWriter::WriteUInt8(path, shape, values);
    if (!written.ok) {
        *error = QStringLiteral("%1: %2").arg(name, written.error);
        return false;
    }
    if (!UInt8ArrayMatches(path, shape, values, error))
        return false;
    files->push_back({name, QStringLiteral("|u1"), shape, false});
    return true;
}

QJsonObject BuildFileManifest(const QString& directory,
                              const std::vector<ExportedFile>& files,
                              QString* error,
                              const ReaderAlignmentExportControl* control = nullptr) {
    QJsonObject manifest;
    for (const ExportedFile& file : files) {
        const QString path = QDir(directory).filePath(file.name);
        const QString hash = Sha256File(path, error, control);
        if (hash.isEmpty())
            return {};
        QJsonObject item{
            {"size", QFileInfo(path).size()},
            {"sha256", hash},
        };
        if (file.byteCopy) {
            item.insert(QStringLiteral("kind"), QStringLiteral("byte_copy"));
        } else {
            item.insert(QStringLiteral("kind"), QStringLiteral("npy"));
            item.insert(QStringLiteral("dtype"), file.dtype);
            item.insert(QStringLiteral("shape"), ShapeToJson(file.shape));
            item.insert(QStringLiteral("order"), QStringLiteral("C"));
        }
        manifest.insert(file.name, item);
    }
    return manifest;
}

bool VerifyRecipe(const QString& path, QString* error) {
    QJsonObject object;
    if (!ReadJsonObject(path, &object, error))
        return false;
    const QString recipe = object.value(QStringLiteral("recipe_id")).toString();
    if (recipe != QString::fromLatin1(kRecipeId)) {
        *error = QStringLiteral("native MOPAC recipe mismatch in %1: %2")
            .arg(path, recipe);
        return false;
    }
    return true;
}

bool ReadNumericContract(const QString& memberDirectory,
                         const QString& name,
                         const QString& dtype,
                         const std::vector<std::size_t>& shape,
                         QtNpyReader::NumericArray* array,
                         QString* error) {
    const QString path = QDir(memberDirectory).filePath(name);
    *array = QtNpyReader::ReadNumericArrayWidened(path);
    if (!array->ok) {
        *error = array->error;
        return false;
    }
    if (array->shape != shape
        || !NumericDtypeMatches(array->descr, dtype)) {
        *error = QStringLiteral("numeric sidecar contract mismatch: %1")
            .arg(path);
        return false;
    }
    return true;
}

bool JsonVec3(const QJsonValue& value, model::Vec3* result) {
    if (!value.isArray())
        return false;
    const QJsonArray array = value.toArray();
    if (array.size() != 3)
        return false;
    for (const auto coordinate : array) {
        if (!coordinate.isDouble() || !std::isfinite(coordinate.toDouble()))
            return false;
    }
    *result = model::Vec3(array[0].toDouble(), array[1].toDouble(),
                         array[2].toDouble());
    return true;
}

bool WidenedIndices(const QtNpyReader::NumericArray& array,
                    std::size_t upperBound,
                    std::vector<std::size_t>* indices,
                    QString* error) {
    std::vector<bool> seen(upperBound, false);
    indices->clear();
    indices->reserve(array.data.size());
    for (double raw : array.data) {
        if (!std::isfinite(raw) || std::floor(raw) != raw || raw < 0.0
            || raw >= static_cast<double>(upperBound)) {
            *error = QStringLiteral("fit atom index is not a valid atom-axis index");
            return false;
        }
        const std::size_t index = static_cast<std::size_t>(raw);
        if (seen[index]) {
            *error = QStringLiteral("fit atom index is duplicated");
            return false;
        }
        seen[index] = true;
        indices->push_back(index);
    }
    return true;
}

}  // namespace

CompletedAlignmentValidation ReaderAlignmentExporter::ValidateCompletedMember(
    const QString& memberDirectory,
    const ReaderAlignmentExportControl* control) {
    CompletedAlignmentValidation result;
    if (CancellationRequested(control, &result.error))
        return result;
    const QString exportPath = QDir(memberDirectory).filePath(
        QStringLiteral("export.json"));
    if (!QFileInfo(exportPath).isFile()) {
        result.error = QStringLiteral("completion record is missing: %1")
            .arg(exportPath);
        return result;
    }
    if (!ReadJsonObject(exportPath, &result.manifest, &result.error))
        return result;
    if (!result.manifest.value(QStringLiteral("complete")).toBool(false)
        || result.manifest.value(QStringLiteral("contains_targets")).toBool(true)
        || result.manifest.value(QStringLiteral("protected_targets_opened")).toInt(-1) != 0) {
        result.error = QStringLiteral("completion record has invalid safety flags: %1")
            .arg(exportPath);
        return result;
    }

    if (result.manifest.value(QStringLiteral("schema_version")).toInt(-1) != 1
        || result.manifest.value(QStringLiteral("transform_convention")).toString()
            != QStringLiteral(
                "aligned[t,a] = rotations[t] @ raw[t,a] + translations[t]")) {
        result.error = QStringLiteral("completion record has an unsupported schema or transform convention");
        return result;
    }
    const QJsonObject identity =
        result.manifest.value(QStringLiteral("identity")).toObject();
    const QString memberId = identity.value(QStringLiteral("member_id")).toString();
    const bool stagingDirectory =
        QFileInfo(memberDirectory).fileName().startsWith(QLatin1Char('.'));
    if (memberId.isEmpty()
        || (!stagingDirectory
            && memberId != QFileInfo(memberDirectory).fileName())) {
        result.error = QStringLiteral("completion record member identity does not match its directory");
        return result;
    }
    const QJsonObject source =
        result.manifest.value(QStringLiteral("source")).toObject();
    if (source.value(QStringLiteral("provenance_authority")).toString()
            != QString::fromLatin1(kSourceProvenanceAuthority)
        || source.value(QStringLiteral("producer_catalog_sha256")).toString()
            != QString::fromLatin1(kProducerCatalogSha256)
        || source.value(QStringLiteral("extractor_commit")).toString()
            != QString::fromLatin1(kExtractorCommit)
        || source.value(QStringLiteral("recipe_id")).toString()
            != QString::fromLatin1(kRecipeId)) {
        result.error = QStringLiteral("completion record source provenance is not authoritative");
        return result;
    }

    const QJsonObject dimensions =
        result.manifest.value(QStringLiteral("dimensions")).toObject();
    const qint64 frameCountRaw =
        dimensions.value(QStringLiteral("frames")).toInteger(-1);
    const qint64 atomCountRaw =
        dimensions.value(QStringLiteral("atoms")).toInteger(-1);
    const qint64 primaryCountRaw =
        dimensions.value(QStringLiteral("primary_fit_atoms")).toInteger(-1);
    const qint64 caCountRaw =
        dimensions.value(QStringLiteral("ca_fit_atoms")).toInteger(-1);
    const qint64 residueCountRaw =
        dimensions.value(QStringLiteral("residues")).toInteger(-1);
    const qint64 bondCountRaw =
        dimensions.value(QStringLiteral("bonds")).toInteger(-1);
    const qint64 aromaticRingCountRaw =
        dimensions.value(QStringLiteral("aromatic_rings")).toInteger(-1);
    const qint64 saturatedRingCountRaw =
        dimensions.value(QStringLiteral("saturated_rings")).toInteger(-1);
    const qint64 ringCountRaw =
        dimensions.value(QStringLiteral("rings")).toInteger(-1);
    const qint64 ringMembershipCountRaw =
        dimensions.value(QStringLiteral("ring_memberships")).toInteger(-1);
    if (frameCountRaw != static_cast<qint64>(kExpectedFrameCount)
        || atomCountRaw <= 0
        || primaryCountRaw < 3 || caCountRaw < 3
        || residueCountRaw <= 0 || bondCountRaw < 0
        || aromaticRingCountRaw < 0 || saturatedRingCountRaw < 0
        || ringCountRaw < 0 || ringMembershipCountRaw < 0
        || ringCountRaw != aromaticRingCountRaw + saturatedRingCountRaw) {
        result.error = QStringLiteral("completion record dimensions are invalid");
        return result;
    }
    const std::size_t frameCount = static_cast<std::size_t>(frameCountRaw);
    const std::size_t atomCount = static_cast<std::size_t>(atomCountRaw);
    const std::size_t primaryCount = static_cast<std::size_t>(primaryCountRaw);
    const std::size_t caCount = static_cast<std::size_t>(caCountRaw);

    const std::vector<ExportedFile> numericContracts{
        {QStringLiteral("aligned_positions.npy"), QStringLiteral("<f8"),
         {frameCount, atomCount, 3}, false},
        {QStringLiteral("rotations.npy"), QStringLiteral("<f8"),
         {frameCount, 3, 3}, false},
        {QStringLiteral("translations.npy"), QStringLiteral("<f8"),
         {frameCount, 3}, false},
        {QStringLiteral("original_frame_index.npy"), QStringLiteral("<u8"),
         {frameCount}, false},
        {QStringLiteral("time_ps.npy"), QStringLiteral("<f8"),
         {frameCount}, false},
        {QStringLiteral("fit_atom_index.npy"), QStringLiteral("<u8"),
         {primaryCount}, false},
        {QStringLiteral("fit_reference_positions.npy"), QStringLiteral("<f8"),
         {primaryCount, 3}, false},
        {QStringLiteral("fit_rmsd_A.npy"), QStringLiteral("<f8"),
         {frameCount}, false},
        {QStringLiteral("fit_singular_values.npy"), QStringLiteral("<f8"),
         {frameCount, 3}, false},
        {QStringLiteral("fit_status.npy"), QStringLiteral("|u1"),
         {frameCount}, false},
        {QStringLiteral("ca_rotations.npy"), QStringLiteral("<f8"),
         {frameCount, 3, 3}, false},
        {QStringLiteral("ca_translations.npy"), QStringLiteral("<f8"),
         {frameCount, 3}, false},
        {QStringLiteral("ca_fit_atom_index.npy"), QStringLiteral("<u8"),
         {caCount}, false},
        {QStringLiteral("ca_fit_reference_positions.npy"), QStringLiteral("<f8"),
         {caCount, 3}, false},
        {QStringLiteral("ca_fit_rmsd_A.npy"), QStringLiteral("<f8"),
         {frameCount}, false},
        {QStringLiteral("ca_fit_singular_values.npy"), QStringLiteral("<f8"),
         {frameCount, 3}, false},
        {QStringLiteral("ca_fit_status.npy"), QStringLiteral("|u1"),
         {frameCount}, false},
    };

    static const QStringList requiredFiles{
        QStringLiteral("aligned_positions.npy"),
        QStringLiteral("rotations.npy"),
        QStringLiteral("translations.npy"),
        QStringLiteral("original_frame_index.npy"),
        QStringLiteral("time_ps.npy"),
        QStringLiteral("fit_atom_index.npy"),
        QStringLiteral("fit_reference_positions.npy"),
        QStringLiteral("fit_rmsd_A.npy"),
        QStringLiteral("fit_singular_values.npy"),
        QStringLiteral("fit_status.npy"),
        QStringLiteral("ca_rotations.npy"),
        QStringLiteral("ca_translations.npy"),
        QStringLiteral("ca_fit_atom_index.npy"),
        QStringLiteral("ca_fit_reference_positions.npy"),
        QStringLiteral("ca_fit_rmsd_A.npy"),
        QStringLiteral("ca_fit_singular_values.npy"),
        QStringLiteral("ca_fit_status.npy"),
        QStringLiteral("atoms_category_info.npy"),
        QStringLiteral("residues.npy"),
        QStringLiteral("bonds.npy"),
        QStringLiteral("rings.npy"),
        QStringLiteral("ring_membership.npy"),
        QStringLiteral("extraction_manifest.json"),
        QStringLiteral("native_mopac_complete.json"),
    };

    const QJsonObject files = result.manifest.value(QStringLiteral("files")).toObject();
    for (const QString& name : requiredFiles) {
        if (!files.contains(name) || !files.value(name).isObject()) {
            result.error = QStringLiteral("completion record does not declare %1")
                .arg(name);
            return result;
        }
    }
    if (files.size() != requiredFiles.size()) {
        result.error = QStringLiteral("completion record declares an unexpected file set");
        return result;
    }

    for (auto it = files.begin(); it != files.end(); ++it) {
        if (CancellationRequested(control, &result.error))
            return result;
        const QString path = QDir(memberDirectory).filePath(it.key());
        const QJsonObject declared = it.value().toObject();
        const QFileInfo info(path);
        if (!info.isFile() || info.size() != declared.value(QStringLiteral("size")).toInteger(-1)) {
            result.error = QStringLiteral("declared file is missing or has wrong size: %1")
                .arg(path);
            return result;
        }
        QString hashError;
        const QString hash = Sha256File(path, &hashError, control);
        if (hash.isEmpty() || hash != declared.value(QStringLiteral("sha256")).toString()) {
            result.error = hashError.isEmpty()
                ? QStringLiteral("declared file hash mismatch: %1").arg(path)
                : hashError;
            return result;
        }
        const auto numeric = std::find_if(
            numericContracts.begin(), numericContracts.end(),
            [&it](const ExportedFile& contract) { return contract.name == it.key(); });
        if (numeric != numericContracts.end()) {
            bool declaredShapeOk = false;
            const std::vector<std::size_t> declaredShape = ShapeFromJson(
                declared.value(QStringLiteral("shape")), &declaredShapeOk);
            if (declared.value(QStringLiteral("kind")).toString()
                    != QStringLiteral("npy")
                || declared.value(QStringLiteral("dtype")).toString()
                    != numeric->dtype
                || declared.value(QStringLiteral("order")).toString()
                    != QStringLiteral("C")
                || !declaredShapeOk || declaredShape != numeric->shape) {
                result.error = QStringLiteral("declared NPY contract mismatch: %1")
                    .arg(path);
                return result;
            }
            const QtNpyReader::NumericArray array =
                QtNpyReader::ReadNumericArrayWidened(path);
            if (!array.ok || array.shape != numeric->shape
                || !NumericDtypeMatches(array.descr, numeric->dtype)
                || !std::all_of(array.data.begin(), array.data.end(),
                                [](double value) { return std::isfinite(value); })) {
                result.error = QStringLiteral("declared NPY contract mismatch: %1")
                    .arg(path);
                return result;
            }
        } else if (declared.value(QStringLiteral("kind")).toString()
                       != QStringLiteral("byte_copy")) {
            result.error = QStringLiteral("declared byte-copy contract mismatch: %1")
                .arg(path);
            return result;
        }
    }

    QtNpyReader::NumericArray primaryIndices;
    QtNpyReader::NumericArray caIndices;
    QtNpyReader::NumericArray originalFrames;
    QtNpyReader::NumericArray primaryStatus;
    QtNpyReader::NumericArray caStatus;
    QtNpyReader::NumericArray primaryRotations;
    QtNpyReader::NumericArray caRotations;
    QtNpyReader::NumericArray primarySingularValues;
    QtNpyReader::NumericArray caSingularValues;
    QtNpyReader::NumericArray times;
    QtNpyReader::NumericArray primaryRmsd;
    QtNpyReader::NumericArray caRmsd;
    QtNpyReader::NumericArray primaryReference;
    QtNpyReader::NumericArray caReference;
    QString contractError;
    if (!ReadNumericContract(memberDirectory, QStringLiteral("fit_atom_index.npy"),
                             QStringLiteral("<u8"), {primaryCount},
                             &primaryIndices, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("ca_fit_atom_index.npy"),
                                QStringLiteral("<u8"), {caCount},
                                &caIndices, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("original_frame_index.npy"),
                                QStringLiteral("<u8"), {frameCount},
                                &originalFrames, &contractError)
        || !ReadNumericContract(memberDirectory, QStringLiteral("time_ps.npy"),
                                QStringLiteral("<f8"), {frameCount},
                                &times, &contractError)
        || !ReadNumericContract(memberDirectory, QStringLiteral("fit_status.npy"),
                                QStringLiteral("|u1"), {frameCount},
                                &primaryStatus, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("ca_fit_status.npy"),
                                QStringLiteral("|u1"), {frameCount},
                                &caStatus, &contractError)
        || !ReadNumericContract(memberDirectory, QStringLiteral("rotations.npy"),
                                QStringLiteral("<f8"), {frameCount, 3, 3},
                                &primaryRotations, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("ca_rotations.npy"),
                                QStringLiteral("<f8"), {frameCount, 3, 3},
                                &caRotations, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("fit_singular_values.npy"),
                                QStringLiteral("<f8"), {frameCount, 3},
                                &primarySingularValues, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("ca_fit_singular_values.npy"),
                                QStringLiteral("<f8"), {frameCount, 3},
                                &caSingularValues, &contractError)
        || !ReadNumericContract(memberDirectory, QStringLiteral("fit_rmsd_A.npy"),
                                QStringLiteral("<f8"), {frameCount},
                                &primaryRmsd, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("ca_fit_rmsd_A.npy"),
                                QStringLiteral("<f8"), {frameCount},
                                &caRmsd, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("fit_reference_positions.npy"),
                                QStringLiteral("<f8"), {primaryCount, 3},
                                &primaryReference, &contractError)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("ca_fit_reference_positions.npy"),
                                QStringLiteral("<f8"), {caCount, 3},
                                &caReference, &contractError)) {
        result.error = contractError;
        return result;
    }
    std::vector<std::size_t> checkedIndices;
    if (!WidenedIndices(primaryIndices, atomCount, &checkedIndices,
                        &result.error)
        || !WidenedIndices(caIndices, atomCount, &checkedIndices,
                           &result.error)) {
        return result;
    }
    std::vector<std::uint64_t> frameIdentities;
    frameIdentities.reserve(frameCount);
    for (double raw : originalFrames.data) {
        if (std::floor(raw) != raw || raw < 0.0
            || raw > static_cast<double>(std::uint64_t{1} << 53)) {
            result.error = QStringLiteral("original frame identities are invalid");
            return result;
        }
        frameIdentities.push_back(static_cast<std::uint64_t>(raw));
    }
    if (!ValidateRetainedFrameAxis(frameIdentities, times.data,
                                   &result.error)) {
        return result;
    }
    if (!std::all_of(primaryRmsd.data.begin(), primaryRmsd.data.end(),
                     [](double value) { return value >= 0.0; })
        || !std::all_of(caRmsd.data.begin(), caRmsd.data.end(),
                       [](double value) { return value >= 0.0; })) {
        result.error = QStringLiteral("fit RMSD contains a negative value");
        return result;
    }
    if (result.manifest.value(QStringLiteral("fit_status_meanings")).toObject()
        != StatusMeanings()) {
        result.error = QStringLiteral("fit status meanings are incomplete or invalid");
        return result;
    }
    for (double status : primaryStatus.data) {
        if (status != static_cast<double>(model::ScientificFitStatus::Valid)) {
            result.error = QStringLiteral("primary fit contains a non-valid status");
            return result;
        }
    }
    for (double status : caStatus.data) {
        if (status != static_cast<double>(model::ScientificFitStatus::Valid)) {
            result.error = QStringLiteral("CA fit contains a non-valid status");
            return result;
        }
    }

    auto validateRotations = [frameCount, &result](
        const QtNpyReader::NumericArray& rotations,
        double tolerance,
        const QString& label,
        TransformValidation* summary) {
        if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
            result.error = QStringLiteral("%1 rotation tolerance is invalid")
                .arg(label);
            return false;
        }
        for (std::size_t frame = 0; frame < frameCount; ++frame) {
            model::Mat3 rotation;
            for (int row = 0; row < 3; ++row) {
                for (int column = 0; column < 3; ++column) {
                    rotation(row, column) = rotations.data[
                        frame * 9
                        + static_cast<std::size_t>(row * 3 + column)];
                }
            }
            const double orthogonalityError =
                (rotation * rotation.transpose() - model::Mat3::Identity()).norm();
            const double determinant = rotation.determinant();
            summary->maxOrthogonalityError = std::max(
                summary->maxOrthogonalityError, orthogonalityError);
            summary->minDeterminant = std::min(summary->minDeterminant,
                                               determinant);
            summary->maxDeterminant = std::max(summary->maxDeterminant,
                                               determinant);
            if (std::abs(determinant - 1.0) > tolerance
                || orthogonalityError > tolerance) {
                result.error = QStringLiteral("%1 rotation is invalid at frame %2")
                    .arg(label)
                    .arg(static_cast<qulonglong>(frame));
                return false;
            }
        }
        return true;
    };
    const QJsonObject primaryFit =
        result.manifest.value(QStringLiteral("primary_fit")).toObject();
    const QJsonObject caFit =
        result.manifest.value(QStringLiteral("ca_fit")).toObject();
    auto validateFitPolicy = [&result](
        const QJsonObject& fit,
        const QString& label,
        const QString& expectedAtomSelection,
        const QtNpyReader::NumericArray& reference) {
        model::Vec3 anchor;
        const double convergenceTolerance =
            fit.value(QStringLiteral("convergence_tolerance_A")).toDouble(-1.0);
        const double finalDelta =
            fit.value(QStringLiteral("final_delta_A")).toDouble(-1.0);
        const qint64 iterations =
            fit.value(QStringLiteral("iterations")).toInteger(-1);
        if (fit.value(QStringLiteral("atom_selection")).toString()
                != expectedAtomSelection
            || fit.value(QStringLiteral("fit_policy")).toString()
                != QStringLiteral(
                    "proper Kabsch SVD; accept numerical rank >= 2")
            || fit.value(QStringLiteral("seed_frame_row")).toInteger(-1) != 0
            || fit.value(QStringLiteral("max_iterations")).toInt(-1) != 20
            || convergenceTolerance != 1e-4
            || !fit.value(QStringLiteral("converged")).toBool(false)
            || iterations <= 0 || iterations > 20
            || !std::isfinite(finalDelta) || finalDelta < 0.0
            || finalDelta > convergenceTolerance
            || fit.value(QStringLiteral(
                    "reference_build_degeneracy_count")).toInteger(-1) != 0
            || fit.value(QStringLiteral("temporal_smoothing")).toBool(true)
            || fit.value(QStringLiteral("numerical_rank_rule")).toString()
                != QStringLiteral(
                    "rank = count(sigma_i > max(abs_tol, rel_tol * sigma_max)); accept rank >= 2")
            || fit.value(QStringLiteral("rank_relative_tolerance")).toDouble(-1.0)
                != 1e-12
            || fit.value(QStringLiteral("rank_absolute_tolerance_A2")).toDouble(-1.0)
                != 1e-12
            || fit.value(QStringLiteral("rotation_tolerance")).toDouble(-1.0)
                != 1e-8
            || !JsonVec3(fit.value(QStringLiteral("centroid_anchor_A")),
                         &anchor)) {
            result.error = QStringLiteral("%1 fit policy is invalid").arg(label);
            return false;
        }
        model::Vec3 referenceCentroid = model::Vec3::Zero();
        const std::size_t pointCount = reference.data.size() / 3;
        for (std::size_t point = 0; point < pointCount; ++point) {
            referenceCentroid += model::Vec3(
                reference.data[point * 3], reference.data[point * 3 + 1],
                reference.data[point * 3 + 2]);
        }
        referenceCentroid /= static_cast<double>(pointCount);
        if ((referenceCentroid - anchor).norm() > 1e-10) {
            result.error = QStringLiteral(
                "%1 reference centroid does not match its declared anchor")
                .arg(label);
            return false;
        }
        return true;
    };
    auto validateSingularValues = [frameCount, &result](
        const QtNpyReader::NumericArray& singularValues,
        const QString& label) {
        for (std::size_t frame = 0; frame < frameCount; ++frame) {
            const double first = singularValues.data[frame * 3];
            const double second = singularValues.data[frame * 3 + 1];
            const double third = singularValues.data[frame * 3 + 2];
            const double threshold = std::max(1e-12, 1e-12 * first);
            int rank = 0;
            for (double value : {first, second, third}) {
                if (value > threshold)
                    ++rank;
            }
            if (first < second || second < third || third < 0.0 || rank < 2) {
                result.error = QStringLiteral(
                    "%1 singular values are invalid at frame %2")
                    .arg(label)
                    .arg(static_cast<qulonglong>(frame));
                return false;
            }
        }
        return true;
    };
    TransformValidation primaryRotationSummary;
    TransformValidation caRotationSummary;
    if (!validateFitPolicy(
            primaryFit, QStringLiteral("primary"),
            QStringLiteral(
                "H5 /trajectory/rmsd_tracking/atom_indices; typed NCACO cross-check"),
            primaryReference)
        || !validateFitPolicy(
            caFit, QStringLiteral("CA"),
            QStringLiteral("typed CA subset of primary NCACO selection"),
            caReference)
        || !validateSingularValues(primarySingularValues,
                                   QStringLiteral("primary"))
        || !validateSingularValues(caSingularValues, QStringLiteral("CA"))
        || !validateRotations(
            primaryRotations,
            primaryFit.value(QStringLiteral("rotation_tolerance")).toDouble(-1.0),
            QStringLiteral("primary"), &primaryRotationSummary)
        || !validateRotations(
            caRotations,
            caFit.value(QStringLiteral("rotation_tolerance")).toDouble(-1.0),
            QStringLiteral("CA"), &caRotationSummary)) {
        if (result.error.isEmpty())
            result.error = QStringLiteral("fit completion policy is invalid");
        return result;
    }

    const QJsonObject validation =
        result.manifest.value(QStringLiteral("validation")).toObject();
    const double orthogonality = validation.value(
        QStringLiteral("max_rotation_orthogonality_frobenius_error"))
        .toDouble(-1.0);
    const double determinantMinimum =
        validation.value(QStringLiteral("determinant_min")).toDouble(-1.0);
    const double determinantMaximum =
        validation.value(QStringLiteral("determinant_max")).toDouble(-1.0);
    const double forwardError = validation.value(
        QStringLiteral("max_forward_reconstruction_error_A")).toDouble(-1.0);
    const double inverseError = validation.value(
        QStringLiteral("max_inverse_reconstruction_error_A")).toDouble(-1.0);
    const double roundTripError = validation.value(
        QStringLiteral("max_round_trip_reconstruction_error_A")).toDouble(-1.0);
    const double caOrthogonality = validation.value(
        QStringLiteral("ca_max_rotation_orthogonality_frobenius_error"))
        .toDouble(-1.0);
    const double caDeterminantMinimum =
        validation.value(QStringLiteral("ca_determinant_min")).toDouble(-1.0);
    const double caDeterminantMaximum =
        validation.value(QStringLiteral("ca_determinant_max")).toDouble(-1.0);
    const double caRoundTripError = validation.value(
        QStringLiteral("ca_max_round_trip_reconstruction_error_A"))
        .toDouble(-1.0);
    if (!std::isfinite(orthogonality) || orthogonality < 0.0
        || orthogonality > 1e-8
        || !std::isfinite(determinantMinimum)
        || !std::isfinite(determinantMaximum)
        || std::abs(determinantMinimum - 1.0) > 1e-8
        || std::abs(determinantMaximum - 1.0) > 1e-8
        || !std::isfinite(forwardError) || forwardError < 0.0
        || forwardError > 1e-12
        || !std::isfinite(inverseError) || inverseError < 0.0
        || inverseError > 1e-10
        || !std::isfinite(roundTripError) || roundTripError < 0.0
        || roundTripError > 1e-10
        || !std::isfinite(caOrthogonality) || caOrthogonality < 0.0
        || caOrthogonality > 1e-8
        || !std::isfinite(caDeterminantMinimum)
        || !std::isfinite(caDeterminantMaximum)
        || std::abs(caDeterminantMinimum - 1.0) > 1e-8
        || std::abs(caDeterminantMaximum - 1.0) > 1e-8
        || !std::isfinite(caRoundTripError) || caRoundTripError < 0.0
        || caRoundTripError > 1e-10
        || std::abs(orthogonality
                    - primaryRotationSummary.maxOrthogonalityError) > 1e-15
        || std::abs(determinantMinimum
                    - primaryRotationSummary.minDeterminant) > 1e-15
        || std::abs(determinantMaximum
                    - primaryRotationSummary.maxDeterminant) > 1e-15
        || std::abs(caOrthogonality
                    - caRotationSummary.maxOrthogonalityError) > 1e-15
        || std::abs(caDeterminantMinimum
                    - caRotationSummary.minDeterminant) > 1e-15
        || std::abs(caDeterminantMaximum
                    - caRotationSummary.maxDeterminant) > 1e-15
        || validation.value(QStringLiteral("matrix_serialization")).toString()
            != QStringLiteral(
                "explicit row-major loops; Eigen data() not used")) {
        result.error = QStringLiteral("recorded transform validation is invalid");
        return result;
    }

    result.ok = true;
    return result;
}

model::ScientificAlignmentResult
ReaderAlignmentExporter::LoadCompletedPrimaryAlignment(
    const QString& memberDirectory,
    const ReaderAlignmentExportControl* control) {
    model::ScientificAlignmentResult result;
    const CompletedAlignmentValidation completed =
        ValidateCompletedMember(memberDirectory, control);
    if (!completed.ok) {
        result.error = completed.error;
        return result;
    }

    const QJsonObject dimensions =
        completed.manifest.value(QStringLiteral("dimensions")).toObject();
    const qint64 frameCountRaw =
        dimensions.value(QStringLiteral("frames")).toInteger(-1);
    const qint64 atomCountRaw =
        dimensions.value(QStringLiteral("atoms")).toInteger(-1);
    const qint64 fitAtomCountRaw =
        dimensions.value(QStringLiteral("primary_fit_atoms")).toInteger(-1);
    if (frameCountRaw <= 0 || atomCountRaw <= 0 || fitAtomCountRaw < 3) {
        result.error = QStringLiteral("completed alignment dimensions are invalid");
        return result;
    }
    const std::size_t frameCount = static_cast<std::size_t>(frameCountRaw);
    const std::size_t atomCount = static_cast<std::size_t>(atomCountRaw);
    const std::size_t fitAtomCount = static_cast<std::size_t>(fitAtomCountRaw);

    QtNpyReader::NumericArray rotations;
    QtNpyReader::NumericArray translations;
    QtNpyReader::NumericArray fitIndices;
    QtNpyReader::NumericArray reference;
    QtNpyReader::NumericArray rmsd;
    QtNpyReader::NumericArray singularValues;
    QtNpyReader::NumericArray status;
    if (!ReadNumericContract(memberDirectory, QStringLiteral("rotations.npy"),
                             QStringLiteral("<f8"), {frameCount, 3, 3},
                             &rotations, &result.error)
        || !ReadNumericContract(memberDirectory, QStringLiteral("translations.npy"),
                                QStringLiteral("<f8"), {frameCount, 3},
                                &translations, &result.error)
        || !ReadNumericContract(memberDirectory, QStringLiteral("fit_atom_index.npy"),
                                QStringLiteral("<u8"), {fitAtomCount},
                                &fitIndices, &result.error)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("fit_reference_positions.npy"),
                                QStringLiteral("<f8"), {fitAtomCount, 3},
                                &reference, &result.error)
        || !ReadNumericContract(memberDirectory, QStringLiteral("fit_rmsd_A.npy"),
                                QStringLiteral("<f8"), {frameCount},
                                &rmsd, &result.error)
        || !ReadNumericContract(memberDirectory,
                                QStringLiteral("fit_singular_values.npy"),
                                QStringLiteral("<f8"), {frameCount, 3},
                                &singularValues, &result.error)
        || !ReadNumericContract(memberDirectory, QStringLiteral("fit_status.npy"),
                                QStringLiteral("|u1"), {frameCount},
                                &status, &result.error)) {
        return result;
    }
    if (!WidenedIndices(fitIndices, atomCount, &result.fitAtomIndices,
                        &result.error)) {
        return result;
    }

    const QJsonObject fit =
        completed.manifest.value(QStringLiteral("primary_fit")).toObject();
    result.policy.seedFrame = static_cast<std::size_t>(
        fit.value(QStringLiteral("seed_frame_row")).toInteger(-1));
    result.policy.maxIterations =
        fit.value(QStringLiteral("max_iterations")).toInt(0);
    result.policy.convergenceToleranceAngstrom =
        fit.value(QStringLiteral("convergence_tolerance_A")).toDouble(-1.0);
    result.policy.rankRelativeTolerance =
        fit.value(QStringLiteral("rank_relative_tolerance")).toDouble(-1.0);
    result.policy.rankAbsoluteTolerance =
        fit.value(QStringLiteral("rank_absolute_tolerance_A2")).toDouble(-1.0);
    result.policy.rotationTolerance =
        fit.value(QStringLiteral("rotation_tolerance")).toDouble(-1.0);
    result.mean.converged = fit.value(QStringLiteral("converged")).toBool(false);
    result.mean.iterations = fit.value(QStringLiteral("iterations")).toInt(0);
    result.mean.finalDeltaAngstrom =
        fit.value(QStringLiteral("final_delta_A")).toDouble(-1.0);
    const qint64 degeneracyCount = fit.value(
        QStringLiteral("reference_build_degeneracy_count")).toInteger(-1);
    if (!JsonVec3(fit.value(QStringLiteral("centroid_anchor_A")),
                  &result.mean.centroidAnchor)
        || result.policy.seedFrame >= frameCount
        || result.policy.maxIterations <= 0
        || !(result.policy.convergenceToleranceAngstrom > 0.0)
        || !(result.policy.rankRelativeTolerance >= 0.0)
        || !(result.policy.rankAbsoluteTolerance >= 0.0)
        || !(result.policy.rotationTolerance > 0.0)
        || !result.mean.converged || degeneracyCount < 0) {
        result.error = QStringLiteral("completed primary fit policy is invalid");
        return result;
    }
    result.mean.referenceBuildDegeneracyCount =
        static_cast<std::size_t>(degeneracyCount);

    result.referencePositions.reserve(fitAtomCount);
    for (std::size_t atom = 0; atom < fitAtomCount; ++atom) {
        result.referencePositions.emplace_back(reference.data[atom * 3],
                                               reference.data[atom * 3 + 1],
                                               reference.data[atom * 3 + 2]);
    }
    result.frameFits.reserve(frameCount);
    for (std::size_t frame = 0; frame < frameCount; ++frame) {
        model::ScientificFrameFit frameFit;
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                frameFit.rotation(row, column) =
                    rotations.data[frame * 9
                                   + static_cast<std::size_t>(row * 3 + column)];
            }
        }
        frameFit.translation = model::Vec3(
            translations.data[frame * 3], translations.data[frame * 3 + 1],
            translations.data[frame * 3 + 2]);
        frameFit.rmsdAngstrom = rmsd.data[frame];
        frameFit.singularValues = model::Vec3(
            singularValues.data[frame * 3],
            singularValues.data[frame * 3 + 1],
            singularValues.data[frame * 3 + 2]);
        const double rawStatus = status.data[frame];
        if (!std::isfinite(rawStatus) || std::floor(rawStatus) != rawStatus
            || rawStatus < 0.0
            || rawStatus > static_cast<double>(
                model::ScientificFitStatus::NonFiniteOutput)) {
            result.error = QStringLiteral("completed fit status is invalid");
            return result;
        }
        frameFit.status = static_cast<model::ScientificFitStatus>(
            static_cast<std::uint8_t>(rawStatus));
        const double largest = frameFit.singularValues(0);
        frameFit.rankThreshold = std::max(
            result.policy.rankAbsoluteTolerance,
            result.policy.rankRelativeTolerance * largest);
        for (int value = 0; value < 3; ++value) {
            if (frameFit.singularValues(value) > frameFit.rankThreshold)
                ++frameFit.numericalRank;
        }
        result.frameFits.push_back(std::move(frameFit));
    }

    result.ok = true;
    return result;
}

bool ReaderAlignmentExporter::RebuildMembersTable(const QString& outputRoot,
                                                  QString* error,
                                                  const ReaderAlignmentExportControl* control) {
    QString localError;
    QString* destinationError = error ? error : &localError;
    QDir root(outputRoot);
    if (!root.exists()) {
        *destinationError = QStringLiteral("alignment output root does not exist: %1")
            .arg(outputRoot);
        return false;
    }

    QStringList rows;
    rows.push_back(QStringLiteral(
        "member_id\trun_id\tframes\tatoms\tprimary_fit_atoms\tca_fit_atoms\tprimary_iterations\tprimary_final_delta_A\tgenerated_at_utc"));
    const QFileInfoList directories = root.entryInfoList(
        QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
    for (const QFileInfo& directory : directories) {
        if (CancellationRequested(control, destinationError))
            return false;
        if (directory.fileName().startsWith(QLatin1Char('.')))
            continue;
        const CompletedAlignmentValidation validation =
            ValidateCompletedMember(directory.absoluteFilePath(), control);
        if (!validation.ok) {
            *destinationError = QStringLiteral("invalid final member %1: %2")
                .arg(directory.absoluteFilePath(), validation.error);
            return false;
        }
        const QJsonObject identity =
            validation.manifest.value(QStringLiteral("identity")).toObject();
        const QJsonObject dimensions =
            validation.manifest.value(QStringLiteral("dimensions")).toObject();
        const QJsonObject primary =
            validation.manifest.value(QStringLiteral("primary_fit")).toObject();
        rows.push_back(QStringLiteral("%1\t%2\t%3\t%4\t%5\t%6\t%7\t%8\t%9")
            .arg(identity.value(QStringLiteral("member_id")).toString(),
                 identity.value(QStringLiteral("run_id")).toString())
            .arg(dimensions.value(QStringLiteral("frames")).toInteger())
            .arg(dimensions.value(QStringLiteral("atoms")).toInteger())
            .arg(dimensions.value(QStringLiteral("primary_fit_atoms")).toInteger())
            .arg(dimensions.value(QStringLiteral("ca_fit_atoms")).toInteger())
            .arg(primary.value(QStringLiteral("iterations")).toInteger())
            .arg(primary.value(QStringLiteral("final_delta_A")).toDouble(), 0, 'g', 17)
            .arg(validation.manifest.value(QStringLiteral("generated_at_utc")).toString()));
    }

    QSaveFile table(root.filePath(QStringLiteral("members.tsv")));
    if (!table.open(QIODevice::WriteOnly)) {
        *destinationError = QStringLiteral("could not open members.tsv: %1")
            .arg(table.errorString());
        return false;
    }
    const QByteArray bytes = (rows.join(QLatin1Char('\n')) + QLatin1Char('\n')).toUtf8();
    if (table.write(bytes) != bytes.size() || !table.commit()) {
        *destinationError = QStringLiteral("could not commit members.tsv: %1")
            .arg(table.errorString());
        return false;
    }
    return true;
}

ReaderAlignmentExportPreparation ReaderAlignmentExporter::Prepare(
    const QtLoadResult& loaded,
    const QString& outputRoot) {
    ReaderAlignmentExportPreparation prepared;
    if (!loaded.ok || !loaded.protein || !loaded.conformation) {
        prepared.error = QStringLiteral("no complete Reader load is available");
        return prepared;
    }
    if (loaded.manifest.kind != CalcsetManifest::Kind::Trajectory
        || !loaded.manifest.trajectory) {
        prepared.error = QStringLiteral(
            "scientific alignment export requires a trajectory LGS");
        return prepared;
    }
    const model::TrajectoryConformation* trajectory =
        loaded.conformation->asTrajectory();
    if (!trajectory || !trajectory->h5() || !trajectory->h5()->positions()) {
        prepared.error = QStringLiteral(
            "loaded trajectory does not expose H5 positions");
        return prepared;
    }
    if (outputRoot.trimmed().isEmpty()) {
        prepared.error = QStringLiteral("alignment output root is empty");
        return prepared;
    }

    const CalcsetManifest::Trajectory& source = *loaded.manifest.trajectory;
    const model::QtProtein& protein = *loaded.protein;
    const QtTrajectoryH5& h5 = *trajectory->h5();
    const std::size_t frameCount = h5.frameCount();
    const std::size_t atomCount = h5.atomCount();
    if (frameCount == 0 || atomCount == 0
        || protein.atomCount() != atomCount
        || loaded.conformation->frameCount() != frameCount) {
        prepared.error = QStringLiteral(
            "topology, conformation, and H5 axes do not agree");
        return prepared;
    }

    if (!ValidateRetainedFrameAxis(h5.frameIndices(), h5.frameTimes(),
                                   &prepared.error)) {
        return prepared;
    }
    ReaderAlignmentExportRequest& request = prepared.request;
    request.outputRoot = outputRoot;
    request.datasetId = loaded.manifest.dataset_id;
    request.runId = loaded.manifest.protein_id;
    request.humanName = loaded.manifest.human_name;
    request.lgsPath = loaded.manifest.lgs_path_abspath;
    request.calcsetRoot = loaded.manifest.calcset_root_abspath;
    request.extractionDirectory = source.extraction_dir_abspath;
    request.trajectoryH5 = source.trajectory_h5_abspath;
    request.extractionManifest = source.extraction_manifest_abspath;
    request.readerVersion = QCoreApplication::applicationVersion();
    request.processId = QCoreApplication::applicationPid();
    request.residueCount = protein.residueCount();
    request.bondCount = protein.bondCount();
    request.aromaticRingCount = protein.topology().aromaticRingCount();
    request.saturatedRingCount = protein.topology().saturatedRingCount();
    request.ringCount = protein.ringCount();
    request.ringMembershipCount = protein.ringMembershipCount();
    request.originalFrameIndices = h5.frameIndices();
    request.frameTimesPicoseconds = h5.frameTimes();
    request.positions.frameCount = frameCount;
    request.positions.atomCount = atomCount;
    request.positions.values.reserve(frameCount * atomCount);
    for (std::size_t frame = 0; frame < frameCount; ++frame) {
        for (std::size_t atom = 0; atom < atomCount; ++atom) {
            const model::Vec3 position = h5.positions()->at(atom, frame);
            if (!position.allFinite()) {
                prepared.error = QStringLiteral(
                    "H5 positions contain a nonfinite value at frame %1 atom %2")
                    .arg(static_cast<qulonglong>(frame))
                    .arg(static_cast<qulonglong>(atom));
                return prepared;
            }
            request.positions.values.push_back(position);
        }
    }

    const model::QtRmsdTracking* rmsdTracking = h5.rmsdTracking();
    if (!rmsdTracking || rmsdTracking->n_frames != frameCount
        || rmsdTracking->atom_indices.empty()) {
        prepared.error = QStringLiteral(
            "H5 rmsd_tracking selection is missing or has the wrong frame axis");
        return prepared;
    }

    request.primaryFitAtoms.reserve(rmsdTracking->atom_indices.size());
    std::vector<bool> seenAtom(atomCount, false);
    for (std::int32_t rawIndex : rmsdTracking->atom_indices) {
        if (rawIndex < 0 || static_cast<std::size_t>(rawIndex) >= atomCount) {
            prepared.error = QStringLiteral(
                "rmsd_tracking contains an out-of-range atom index");
            return prepared;
        }
        const std::size_t atom = static_cast<std::size_t>(rawIndex);
        if (seenAtom[atom]) {
            prepared.error = QStringLiteral(
                "rmsd_tracking contains a duplicate atom index");
            return prepared;
        }
        seenAtom[atom] = true;
        request.primaryFitAtoms.push_back(atom);
    }

    std::vector<std::size_t> expectedNcaco;
    expectedNcaco.reserve(request.primaryFitAtoms.size());
    for (std::size_t atom = 0; atom < atomCount; ++atom) {
        const model::QtAtom& identity = protein.atom(atom);
        if (identity.IsBackboneNitrogen()
            || identity.IsBackboneAlphaCarbon()
            || identity.IsBackboneCarbonylCarbon()
            || identity.IsBackboneCarbonylOxygen()) {
            expectedNcaco.push_back(atom);
        }
    }
    if (request.primaryFitAtoms != expectedNcaco) {
        prepared.error = QStringLiteral(
            "rmsd_tracking atom selection does not exactly match typed NCACO atoms");
        return prepared;
    }

    for (std::size_t atom : request.primaryFitAtoms) {
        if (protein.atom(atom).IsBackboneAlphaCarbon())
            request.caFitAtoms.push_back(atom);
    }
    if (request.caFitAtoms.size() < 3) {
        prepared.error = QStringLiteral(
            "CA control selection has fewer than three atoms");
        return prepared;
    }

    prepared.ok = true;
    return prepared;
}

ReaderAlignmentExportResult ReaderAlignmentExporter::Export(
    const ReaderAlignmentExportRequest& request,
    const ReaderAlignmentExportControl* control) {
    ReaderAlignmentExportResult result;
    const auto cancelled = [&]() {
        if (!CancellationRequested(control, &result.error))
            return false;
        result.cancelled = true;
        return true;
    };
    if (cancelled())
        return result;

    const QString sourceMember = CleanAbsolutePath(request.calcsetRoot);
    const QString cleanOutputRoot = CleanAbsolutePath(request.outputRoot);
    if (sourceMember.isEmpty() || cleanOutputRoot.isEmpty()) {
        result.error = QStringLiteral(
            "could not resolve source member or alignment output root");
        return result;
    }
    result.memberId = QFileInfo(sourceMember).fileName();
    if (result.memberId.isEmpty()) {
        result.error = QStringLiteral("could not derive member identity from %1")
            .arg(sourceMember);
        return result;
    }
    result.finalDirectory = QDir(cleanOutputRoot).filePath(result.memberId);
    if (PathsOverlap(sourceMember, result.finalDirectory)) {
        result.error = QStringLiteral(
            "source and final output paths overlap: %1 ; %2")
            .arg(sourceMember, result.finalDirectory);
        return result;
    }

    const std::size_t frameCount = request.positions.frameCount;
    const std::size_t atomCount = request.positions.atomCount;
    if (frameCount == 0 || atomCount == 0
        || request.positions.values.size() != frameCount * atomCount
        || request.originalFrameIndices.size() != frameCount
        || request.frameTimesPicoseconds.size() != frameCount) {
        result.error = QStringLiteral("prepared alignment axes do not agree");
        return result;
    }
    if (!ValidateRetainedFrameAxis(request.originalFrameIndices,
                                   request.frameTimesPicoseconds,
                                   &result.error)) {
        return result;
    }

    const QtManifest topologyManifest = QtManifest::Load(
        request.extractionManifest);
    if (!topologyManifest.ok) {
        result.error = QStringLiteral("could not validate topology axes: %1")
            .arg(topologyManifest.error);
        return result;
    }
    const QtManifestAxisSizes& axes = topologyManifest.axisSizes;
    QStringList axisMismatches;
    const auto compareAxis = [&axisMismatches](const char* name,
                                               std::size_t declared,
                                               std::size_t loadedCount) {
        if (declared != loadedCount) {
            axisMismatches.push_back(QStringLiteral("%1 %2 != %3")
                .arg(QString::fromLatin1(name))
                .arg(static_cast<qulonglong>(declared))
                .arg(static_cast<qulonglong>(loadedCount)));
        }
    };
    compareAxis("atom", axes.atom, atomCount);
    compareAxis("residue", axes.residue, request.residueCount);
    compareAxis("bond", axes.bond, request.bondCount);
    compareAxis("aromatic_ring", axes.aromaticRing,
                request.aromaticRingCount);
    compareAxis("saturated_ring", axes.saturatedRing,
                request.saturatedRingCount);
    compareAxis("ring", axes.ring, request.ringCount);
    compareAxis("ring_membership", axes.ringMembership,
                request.ringMembershipCount);
    if (!axisMismatches.isEmpty()) {
        result.error = QStringLiteral(
            "topology axis mismatch (manifest != loaded): %1")
            .arg(axisMismatches.join(QStringLiteral(", ")));
        return result;
    }

    if (!QDir().mkpath(cleanOutputRoot)) {
        result.error = QStringLiteral(
            "could not create alignment output root: %1")
            .arg(cleanOutputRoot);
        return result;
    }

    if (QFileInfo::exists(result.finalDirectory)) {
        const CompletedAlignmentValidation existing =
            ValidateCompletedMember(result.finalDirectory, control);
        if (!existing.ok) {
            if (cancelled())
                return result;
            result.error = QStringLiteral("final member exists but is invalid: %1")
                .arg(existing.error);
            return result;
        }
        const QJsonObject existingIdentity =
            existing.manifest.value(QStringLiteral("identity")).toObject();
        const QJsonObject existingSource =
            existing.manifest.value(QStringLiteral("source")).toObject();
        if (existingIdentity.value(QStringLiteral("member_id")).toString()
                != result.memberId
            || existingIdentity.value(QStringLiteral("dataset_id")).toString()
                != request.datasetId
            || existingIdentity.value(QStringLiteral("run_id")).toString()
                != request.runId
            || CleanAbsolutePath(existingIdentity.value(
                    QStringLiteral("lgs_path")).toString())
                .compare(CleanAbsolutePath(request.lgsPath),
                         PathCaseSensitivity()) != 0
            || CleanAbsolutePath(existingSource.value(
                    QStringLiteral("calcset_root")).toString())
                .compare(sourceMember, PathCaseSensitivity()) != 0
            || CleanAbsolutePath(existingSource.value(
                    QStringLiteral("trajectory_h5")).toString())
                .compare(CleanAbsolutePath(request.trajectoryH5),
                         PathCaseSensitivity()) != 0) {
            result.error = QStringLiteral(
                "completed member belongs to a different loaded source");
            return result;
        }
        QString tableError;
        if (!RebuildMembersTable(cleanOutputRoot, &tableError, control)) {
            if (cancelled())
                return result;
            result.error = tableError;
            return result;
        }
        result.ok = true;
        result.alreadyComplete = true;
        result.manifest = existing.manifest;
        return result;
    }

    const model::ScientificPositionTable& raw = request.positions;
    const std::vector<std::size_t>& primaryAtoms = request.primaryFitAtoms;
    const std::vector<std::size_t>& caAtoms = request.caFitAtoms;
    if (primaryAtoms.size() < 3 || caAtoms.size() < 3) {
        result.error = QStringLiteral("prepared fit selections are incomplete");
        return result;
    }

    const PositionSourceValidation sourcePositions =
        ValidatePositionSources(request.extractionDirectory, raw,
                                request.originalFrameIndices, control);
    if (!sourcePositions.ok) {
        if (cancelled())
            return result;
        result.error = sourcePositions.error;
        return result;
    }

    model::ScientificAlignmentPolicy policy;
    policy.seedFrame = 0;
    policy.maxIterations = 20;
    policy.convergenceToleranceAngstrom = 1e-4;
    result.primaryAlignment =
        model::BuildScientificAlignment(raw, primaryAtoms, policy);
    if (!result.primaryAlignment.ok) {
        result.error = QStringLiteral("primary scientific alignment failed: %1")
            .arg(result.primaryAlignment.error);
        return result;
    }
    result.caAlignment = model::BuildScientificAlignment(raw, caAtoms, policy);
    if (!result.caAlignment.ok) {
        result.error = QStringLiteral("CA scientific alignment failed: %1")
            .arg(result.caAlignment.error);
        return result;
    }
    if (cancelled())
        return result;

    const QString stageName = QStringLiteral(".%1.staging.%2.%3")
        .arg(result.memberId)
        .arg(request.processId)
        .arg(QUuid::createUuid().toString(QUuid::WithoutBraces));
    const QString stagingDirectory = QDir(cleanOutputRoot).filePath(stageName);
    if (!SameOrChildPath(stagingDirectory, cleanOutputRoot)
        || SameOrChildPath(stagingDirectory, sourceMember)
        || !QDir().mkpath(stagingDirectory)) {
        result.error = QStringLiteral("could not create a safe staging directory: %1")
            .arg(stagingDirectory);
        return result;
    }
    StagingDirectoryGuard stagingGuard(stagingDirectory);

    std::vector<double> alignedPositions;
    alignedPositions.reserve(frameCount * atomCount * 3);
    std::vector<double> rotations;
    rotations.reserve(frameCount * 9);
    std::vector<double> translations;
    translations.reserve(frameCount * 3);
    std::vector<double> fitRmsd;
    fitRmsd.reserve(frameCount);
    std::vector<double> fitSingularValues;
    fitSingularValues.reserve(frameCount * 3);
    std::vector<std::uint8_t> fitStatus;
    fitStatus.reserve(frameCount);

    for (std::size_t frame = 0; frame < frameCount; ++frame) {
        if (cancelled())
            return result;
        const model::ScientificFrameFit& fit =
            result.primaryAlignment.frameFits[frame];
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column)
                rotations.push_back(fit.rotation(row, column));
        }
        translations.push_back(fit.translation.x());
        translations.push_back(fit.translation.y());
        translations.push_back(fit.translation.z());
        fitRmsd.push_back(fit.rmsdAngstrom);
        fitSingularValues.push_back(fit.singularValues(0));
        fitSingularValues.push_back(fit.singularValues(1));
        fitSingularValues.push_back(fit.singularValues(2));
        fitStatus.push_back(static_cast<std::uint8_t>(fit.status));

        for (std::size_t atom = 0; atom < atomCount; ++atom) {
            const model::Vec3 aligned = fit.rotation * raw.at(frame, atom)
                + fit.translation;
            alignedPositions.push_back(aligned.x());
            alignedPositions.push_back(aligned.y());
            alignedPositions.push_back(aligned.z());
        }
    }

    std::vector<double> caRotations;
    caRotations.reserve(frameCount * 9);
    std::vector<double> caTranslations;
    caTranslations.reserve(frameCount * 3);
    std::vector<double> caRmsd;
    caRmsd.reserve(frameCount);
    std::vector<double> caSingularValues;
    caSingularValues.reserve(frameCount * 3);
    std::vector<std::uint8_t> caStatus;
    caStatus.reserve(frameCount);
    for (const model::ScientificFrameFit& fit : result.caAlignment.frameFits) {
        if (cancelled())
            return result;
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column)
                caRotations.push_back(fit.rotation(row, column));
        }
        caTranslations.push_back(fit.translation.x());
        caTranslations.push_back(fit.translation.y());
        caTranslations.push_back(fit.translation.z());
        caRmsd.push_back(fit.rmsdAngstrom);
        caSingularValues.push_back(fit.singularValues(0));
        caSingularValues.push_back(fit.singularValues(1));
        caSingularValues.push_back(fit.singularValues(2));
        caStatus.push_back(static_cast<std::uint8_t>(fit.status));
    }

    const TransformValidation transformValidation = ValidateTransforms(
        raw, result.primaryAlignment, &alignedPositions);
    if (!transformValidation.ok) {
        result.error = transformValidation.error;
        return result;
    }
    const TransformValidation caTransformValidation = ValidateTransforms(
        raw, result.caAlignment, nullptr);
    if (!caTransformValidation.ok) {
        result.error = caTransformValidation.error;
        return result;
    }

    std::vector<ExportedFile> files;
    QString writeError;
    const std::vector<std::uint64_t>& originalFrames =
        request.originalFrameIndices;
    const std::vector<double>& times = request.frameTimesPicoseconds;
    const std::vector<std::uint64_t> primaryIndices = WidenIndices(primaryAtoms);
    const std::vector<std::uint64_t> caIndices = WidenIndices(caAtoms);
    const std::vector<double> primaryReference = FlattenReference(result.primaryAlignment);
    const std::vector<double> caReference = FlattenReference(result.caAlignment);

    if (!WriteFloatFile(stagingDirectory, QStringLiteral("aligned_positions.npy"),
                        {frameCount, atomCount, 3}, alignedPositions, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("rotations.npy"),
                           {frameCount, 3, 3}, rotations, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("translations.npy"),
                           {frameCount, 3}, translations, &files, &writeError)
        || !WriteUIntFile(stagingDirectory, QStringLiteral("original_frame_index.npy"),
                          {frameCount}, originalFrames, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("time_ps.npy"),
                           {frameCount}, times, &files, &writeError)
        || !WriteUIntFile(stagingDirectory, QStringLiteral("fit_atom_index.npy"),
                          {primaryAtoms.size()}, primaryIndices, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("fit_reference_positions.npy"),
                           {primaryAtoms.size(), 3}, primaryReference, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("fit_rmsd_A.npy"),
                           {frameCount}, fitRmsd, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("fit_singular_values.npy"),
                           {frameCount, 3}, fitSingularValues, &files, &writeError)
        || !WriteStatusFile(stagingDirectory, QStringLiteral("fit_status.npy"),
                            {frameCount}, fitStatus, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("ca_rotations.npy"),
                           {frameCount, 3, 3}, caRotations, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("ca_translations.npy"),
                           {frameCount, 3}, caTranslations, &files, &writeError)
        || !WriteUIntFile(stagingDirectory, QStringLiteral("ca_fit_atom_index.npy"),
                          {caAtoms.size()}, caIndices, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("ca_fit_reference_positions.npy"),
                           {caAtoms.size(), 3}, caReference, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("ca_fit_rmsd_A.npy"),
                           {frameCount}, caRmsd, &files, &writeError)
        || !WriteFloatFile(stagingDirectory, QStringLiteral("ca_fit_singular_values.npy"),
                           {frameCount, 3}, caSingularValues, &files, &writeError)
        || !WriteStatusFile(stagingDirectory, QStringLiteral("ca_fit_status.npy"),
                            {frameCount}, caStatus, &files, &writeError)) {
        result.error = writeError;
        return result;
    }

    const std::vector<std::pair<QString, QString>> byteCopies{
        {QStringLiteral("atoms_category_info.npy"),
         QDir(request.extractionDirectory).filePath(QStringLiteral("atoms_category_info.npy"))},
        {QStringLiteral("residues.npy"),
         QDir(request.extractionDirectory).filePath(QStringLiteral("residues.npy"))},
        {QStringLiteral("bonds.npy"),
         QDir(request.extractionDirectory).filePath(QStringLiteral("bonds.npy"))},
        {QStringLiteral("rings.npy"),
         QDir(request.extractionDirectory).filePath(QStringLiteral("rings.npy"))},
        {QStringLiteral("ring_membership.npy"),
         QDir(request.extractionDirectory).filePath(QStringLiteral("ring_membership.npy"))},
        {QStringLiteral("extraction_manifest.json"), request.extractionManifest},
        {QStringLiteral("native_mopac_complete.json"),
         QDir(request.extractionDirectory).filePath(QStringLiteral("native_mopac_complete.json"))},
    };
    for (const auto& copy : byteCopies) {
        if (cancelled())
            return result;
        if (!CopyFileExact(copy.second,
                           QDir(stagingDirectory).filePath(copy.first),
                           &writeError)) {
            result.error = writeError;
            return result;
        }
        files.push_back({copy.first, {}, {}, true});
    }
    if (!VerifyRecipe(QDir(stagingDirectory).filePath(
            QStringLiteral("native_mopac_complete.json")), &writeError)) {
        result.error = writeError;
        return result;
    }

    QJsonArray inputFiles;
    if (!AddInputHash(QStringLiteral("lgs"), request.lgsPath,
                      &inputFiles, &writeError, std::nullopt, std::nullopt,
                      control)
        || !AddInputHash(QStringLiteral("trajectory_h5"), request.trajectoryH5,
                         &inputFiles, &writeError, std::nullopt, std::nullopt,
                         control)) {
        if (cancelled())
            return result;
        result.error = writeError;
        return result;
    }
    for (const auto& copy : byteCopies) {
        if (!AddInputHash(QStringLiteral("copied_source"), copy.second,
                          &inputFiles, &writeError, std::nullopt, std::nullopt,
                          control)) {
            if (cancelled())
                return result;
            result.error = writeError;
            return result;
        }
    }
    for (std::size_t frame = 0; frame < frameCount; ++frame) {
        if (!AddInputHash(QStringLiteral("frame_position"),
                          sourcePositions.paths[static_cast<qsizetype>(frame)],
                          &inputFiles, &writeError, frame, originalFrames[frame],
                          control)) {
            if (cancelled())
                return result;
            result.error = writeError;
            return result;
        }
    }

    const QJsonObject fileManifest =
        BuildFileManifest(stagingDirectory, files, &writeError, control);
    if (fileManifest.isEmpty()) {
        if (cancelled())
            return result;
        result.error = writeError;
        return result;
    }

    const QString generatedAt =
        QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs);
    QJsonObject exportManifest{
        {"schema_version", 1},
        {"complete", true},
        {"generated_at_utc", generatedAt},
        {"contains_targets", false},
        {"protected_targets_opened", 0},
        {"transform_convention",
         QStringLiteral("aligned[t,a] = rotations[t] @ raw[t,a] + translations[t]")},
        {"identity", QJsonObject{
            {"member_id", result.memberId},
            {"dataset_id", request.datasetId},
            {"run_id", request.runId},
            {"human_name", request.humanName},
            {"lgs_path", CleanAbsolutePath(request.lgsPath)},
        }},
        {"source", QJsonObject{
            {"calcset_root", sourceMember},
            {"extraction_directory", CleanAbsolutePath(request.extractionDirectory)},
            {"trajectory_h5", CleanAbsolutePath(request.trajectoryH5)},
            {"provenance_authority", QString::fromLatin1(kSourceProvenanceAuthority)},
            {"producer_catalog_sha256", QString::fromLatin1(kProducerCatalogSha256)},
            {"extractor_commit", QString::fromLatin1(kExtractorCommit)},
            {"recipe_id", QString::fromLatin1(kRecipeId)},
            {"input_files", inputFiles},
        }},
        {"reader", QJsonObject{
            {"version", request.readerVersion},
            {"base_commit", QString::fromLatin1(kReaderBaseCommit)},
            {"final_commit", QStringLiteral(H5READER_GIT_COMMIT)},
            {"dirty", static_cast<bool>(H5READER_GIT_DIRTY)},
        }},
        {"dimensions", QJsonObject{
            {"frames", static_cast<qint64>(frameCount)},
            {"atoms", static_cast<qint64>(atomCount)},
            {"residues", static_cast<qint64>(request.residueCount)},
            {"bonds", static_cast<qint64>(request.bondCount)},
            {"aromatic_rings", static_cast<qint64>(request.aromaticRingCount)},
            {"saturated_rings", static_cast<qint64>(request.saturatedRingCount)},
            {"rings", static_cast<qint64>(request.ringCount)},
            {"ring_memberships", static_cast<qint64>(request.ringMembershipCount)},
            {"primary_fit_atoms", static_cast<qint64>(primaryAtoms.size())},
            {"ca_fit_atoms", static_cast<qint64>(caAtoms.size())},
        }},
        {"primary_fit", FitToJson(result.primaryAlignment,
                                  QStringLiteral("H5 /trajectory/rmsd_tracking/atom_indices; typed NCACO cross-check"))},
        {"ca_fit", FitToJson(result.caAlignment,
                             QStringLiteral("typed CA subset of primary NCACO selection"))},
        {"fit_status_meanings", StatusMeanings()},
        {"validation", QJsonObject{
            {"source_position_comparison",
             QStringLiteral("all H5 positions compared with corresponding pos.npy using abs+relative coordinate tolerance")},
            {"source_position_abs_tolerance_A", kSourcePositionAbsoluteToleranceAngstrom},
            {"source_position_relative_tolerance", kSourcePositionRelativeTolerance},
            {"max_source_coordinate_difference_A",
             sourcePositions.maxCoordinateDifferenceAngstrom},
            {"max_source_point_difference_A",
             sourcePositions.maxPointDifferenceAngstrom},
            {"max_rotation_orthogonality_frobenius_error",
             transformValidation.maxOrthogonalityError},
            {"determinant_min", transformValidation.minDeterminant},
            {"determinant_max", transformValidation.maxDeterminant},
            {"max_forward_reconstruction_error_A",
             transformValidation.maxForwardReconstructionErrorAngstrom},
            {"max_inverse_reconstruction_error_A",
             transformValidation.maxInverseReconstructionErrorAngstrom},
            {"max_round_trip_reconstruction_error_A",
             transformValidation.maxRoundTripReconstructionErrorAngstrom},
            {"ca_max_rotation_orthogonality_frobenius_error",
             caTransformValidation.maxOrthogonalityError},
            {"ca_determinant_min", caTransformValidation.minDeterminant},
            {"ca_determinant_max", caTransformValidation.maxDeterminant},
            {"ca_max_round_trip_reconstruction_error_A",
             caTransformValidation.maxRoundTripReconstructionErrorAngstrom},
            {"matrix_serialization", QStringLiteral("explicit row-major loops; Eigen data() not used")},
        }},
        {"files", fileManifest},
    };

    const QString exportPath = QDir(stagingDirectory).filePath(
        QStringLiteral("export.json"));
    if (!WriteJsonObject(exportPath, exportManifest, &writeError)) {
        result.error = writeError;
        return result;
    }
    const CompletedAlignmentValidation staged =
        ValidateCompletedMember(stagingDirectory, control);
    if (!staged.ok) {
        if (cancelled())
            return result;
        result.error = QStringLiteral("staged member validation failed: %1")
            .arg(staged.error);
        return result;
    }

    QDir root(cleanOutputRoot);
    if (!root.rename(stageName, result.memberId)) {
        result.error = QStringLiteral("could not publish staged member %1 as %2")
            .arg(stagingDirectory, result.finalDirectory);
        return result;
    }
    stagingGuard.release();

    // Publication is the cancellation boundary. Once the complete staging
    // directory has its final name, finish validation and the members table so
    // the published collection cannot be left internally inconsistent.
    const CompletedAlignmentValidation published =
        ValidateCompletedMember(result.finalDirectory, nullptr);
    if (!published.ok) {
        result.error = QStringLiteral("published member validation failed: %1")
            .arg(published.error);
        return result;
    }
    if (!RebuildMembersTable(cleanOutputRoot, &writeError, nullptr)) {
        result.error = writeError;
        return result;
    }

    result.manifest = published.manifest;
    result.ok = true;
    qCInfo(cAlignmentExport).noquote()
        << "scientific alignment exported"
        << "| member=" << result.memberId
        << "| frames=" << static_cast<qulonglong>(frameCount)
        << "| atoms=" << static_cast<qulonglong>(atomCount)
        << "| primary_iterations=" << result.primaryAlignment.mean.iterations
        << "| ca_iterations=" << result.caAlignment.mean.iterations
        << "| output=" << result.finalDirectory;
    return result;
}

}  // namespace h5reader::io
