#include "ExperimentalShieldingMlStore.h"

#include "Conformation.h"
#include "QtConformationSnapshot.h"
#include "QtProtein.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ErrorBus.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/FrameFieldPolicy.h"

#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLoggingCategory>
#include <QMetaObject>
#include <QProcessEnvironment>
#include <QTextStream>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <utility>

namespace h5reader::model {

namespace {

Q_LOGGING_CATEGORY(cExperimentalMl, "h5reader.experimental_ml")

constexpr int kOutputColumns = 6;
constexpr qint64 kModelHashChunkBytes = 8LL * 1024LL * 1024LL;
constexpr double kSmoothFiniteNormalisation = 1.14136;
constexpr double kPi = 3.14159265358979323846;

bool checkedSizeToInt(qsizetype size,
                      const QString& context,
                      int& out,
                      QString& error) {
    if (size < 0
        || size > static_cast<qsizetype>(std::numeric_limits<int>::max())) {
        error = QStringLiteral("%1 exceeds the reader's integer channel limit")
                    .arg(context);
        return false;
    }
    out = static_cast<int>(size);
    return true;
}

bool readJsonObject(const QString& path, QJsonObject& out, QString& error) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        error = QStringLiteral("could not open JSON file: %1").arg(path);
        return false;
    }
    QJsonParseError parseError{};
    const QJsonDocument document = QJsonDocument::fromJson(file.readAll(), &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject()) {
        error = QStringLiteral("invalid JSON in %1: %2").arg(path, parseError.errorString());
        return false;
    }
    out = document.object();
    return true;
}

bool readStringList(const QJsonValue& value,
                    const QString& context,
                    QStringList& out,
                    QString& error) {
    if (!value.isArray()) {
        error = QStringLiteral("%1 is not an array").arg(context);
        return false;
    }
    for (const QJsonValue item : value.toArray()) {
        if (!item.isString()) {
            error = QStringLiteral("%1 contains a non-string value").arg(context);
            return false;
        }
        out.append(item.toString());
    }
    return true;
}

bool readIntVector(const QJsonValue& value,
                   const QString& context,
                   QVector<int>& out,
                   QString& error) {
    if (!value.isArray()) {
        error = QStringLiteral("%1 is not an array").arg(context);
        return false;
    }
    for (const QJsonValue item : value.toArray()) {
        if (!item.isDouble()) {
            error = QStringLiteral("%1 contains a non-numeric value").arg(context);
            return false;
        }
        const double raw = item.toDouble();
        if (!std::isfinite(raw)
            || raw < static_cast<double>(std::numeric_limits<int>::min())
            || raw > static_cast<double>(std::numeric_limits<int>::max())) {
            error = QStringLiteral("%1 contains an out-of-range integer value")
                        .arg(context);
            return false;
        }
        const int converted = static_cast<int>(raw);
        if (raw != static_cast<double>(converted)) {
            error = QStringLiteral("%1 contains a non-integer value").arg(context);
            return false;
        }
        out.push_back(converted);
    }
    return true;
}

bool readDoubleVector(const QJsonValue& value,
                      const QString& context,
                      QVector<double>& out,
                      QString& error) {
    if (!value.isArray()) {
        error = QStringLiteral("%1 is not an array").arg(context);
        return false;
    }
    for (const QJsonValue item : value.toArray()) {
        if (!item.isDouble() || !std::isfinite(item.toDouble())) {
            error = QStringLiteral("%1 contains a non-finite numeric value").arg(context);
            return false;
        }
        out.push_back(item.toDouble());
    }
    return true;
}

bool validateRuntimeManifest(const QJsonObject& manifest, QString& error) {
    const QJsonObject schema =
        manifest.value(QStringLiteral("inference_schema")).toObject();
    const QJsonObject runtime =
        manifest.value(QStringLiteral("runtime_contract")).toObject();
    if (schema.value(QStringLiteral("version")).toInt() != 2
        || schema.value(QStringLiteral("model_id")).toString().isEmpty()
        || !schema.value(QStringLiteral("numeric_features")).isArray()
        || !schema.value(QStringLiteral("label_keys")).isArray()
        || !schema.value(QStringLiteral("label_vocabs")).isObject()
        || !schema.value(QStringLiteral("expected_channels")).isObject()
        || !schema.value(QStringLiteral("expected_channel_names")).isObject()
        || runtime.value(QStringLiteral("helper_protocol")).toInt() != 2) {
        error = QStringLiteral("manifest inference_schema/runtime protocol 2 is missing");
        return false;
    }

    const QString modelFile =
        manifest.value(QStringLiteral("model_file")).toString();
    if (modelFile != QStringLiteral("model.ts")) {
        error = QStringLiteral("manifest model_file is %1; expected model.ts")
                    .arg(modelFile.isEmpty() ? QStringLiteral("<missing>") : modelFile);
        return false;
    }
    if (manifest.value(QStringLiteral("target")).toString()
        != QStringLiteral("ORCA total shielding: isotropic 0e plus traceless 2e")) {
        error = QStringLiteral("manifest target is not the F003 total shielding tensor");
        return false;
    }

    const QJsonObject output =
        manifest.value(QStringLiteral("output_contract")).toObject();
    if (output.value(QStringLiteral("dtype")).toString() != QStringLiteral("float32")
        || output.value(QStringLiteral("shape")).toString() != QStringLiteral("N x 6")
        || output.value(QStringLiteral("tensor_basis")).toString()
               != QStringLiteral("project 0e plus traceless 2e")
        || output.value(QStringLiteral("t1")).toString()
               != QStringLiteral("not predicted and not synthesized")) {
        error = QStringLiteral("manifest output_contract is not the F003 0e plus 2e contract");
        return false;
    }

    QStringList outputColumns;
    if (!readStringList(output.value(QStringLiteral("columns")),
                        QStringLiteral("output_contract.columns"),
                        outputColumns,
                        error)) {
        return false;
    }
    const QStringList expectedOutputColumns{
        QStringLiteral("sigma_iso_ppm"),
        QStringLiteral("T2_project_xy"),
        QStringLiteral("T2_project_yz"),
        QStringLiteral("T2_project_axial_z"),
        QStringLiteral("T2_project_xz"),
        QStringLiteral("T2_project_x2_minus_y2"),
    };
    if (outputColumns != expectedOutputColumns) {
        error = QStringLiteral("manifest output_contract.columns do not match the reader basis");
        return false;
    }

    const QJsonValue pythonRequired =
        runtime.value(QStringLiteral("python_required"));
    if (!pythonRequired.isBool() || pythonRequired.toBool()) {
        error = QStringLiteral("manifest runtime_contract must declare python_required=false");
        return false;
    }

    const QString modelId =
        schema.value(QStringLiteral("model_id")).toString();
    int matchingModels = 0;
    for (const QJsonValue value :
         manifest.value(QStringLiteral("models")).toArray()) {
        const QJsonObject model = value.toObject();
        if (model.value(QStringLiteral("id")).toString() != modelId)
            continue;
        ++matchingModels;
        if (model.value(QStringLiteral("model_file")).toString() != modelFile) {
            error = QStringLiteral("manifest model %1 does not use %2")
                        .arg(modelId, modelFile);
            return false;
        }
    }
    if (matchingModels != 1) {
        error = QStringLiteral("manifest model_id %1 has %2 matching model entries")
                    .arg(modelId)
                    .arg(matchingModels);
        return false;
    }
    return true;
}

bool validateModelArtifact(const QJsonObject& manifest,
                           const QString& modelPath,
                           QString& error) {
    const QString modelFile =
        manifest.value(QStringLiteral("model_file")).toString();
    if (QFileInfo(modelPath).fileName() != modelFile) {
        error = QStringLiteral("configured model is %1; manifest requires %2")
                    .arg(QFileInfo(modelPath).fileName(), modelFile);
        return false;
    }

    const QString expected =
        manifest.value(QStringLiteral("model_sha256")).toString().toLower();
    if (expected.size() != 64) {
        error = QStringLiteral("manifest model_sha256 is missing or malformed");
        return false;
    }

    QFile model(modelPath);
    if (!model.open(QIODevice::ReadOnly)) {
        error = QStringLiteral("could not open model for SHA-256 validation: %1")
                    .arg(modelPath);
        return false;
    }
    QCryptographicHash hash(QCryptographicHash::Sha256);
    while (!model.atEnd()) {
        const QByteArray chunk = model.read(kModelHashChunkBytes);
        if (chunk.isEmpty() && model.error() != QFileDevice::NoError) {
            error = QStringLiteral("could not read model for SHA-256 validation: %1")
                        .arg(modelPath);
            return false;
        }
        hash.addData(chunk);
    }
    const QString observed =
        QString::fromLatin1(hash.result().toHex());
    if (observed != expected) {
        error = QStringLiteral("model SHA-256 mismatch: manifest=%1 file=%2")
                    .arg(expected, observed);
        return false;
    }
    return true;
}

bool npyStem(const QString& fileName, QString& stem, QString& error) {
    if (fileName.isEmpty() || QFileInfo(fileName).fileName() != fileName
        || !fileName.endsWith(QStringLiteral(".npy"), Qt::CaseSensitive)) {
        error = QStringLiteral("model input is not a direct .npy filename: %1").arg(fileName);
        return false;
    }
    stem = fileName.left(fileName.size() - 4);
    return !stem.isEmpty();
}

bool writeRaw(const QString& path, const void* data, qsizetype bytes, QString& error) {
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
        error = QStringLiteral("could not create %1").arg(path);
        return false;
    }
    if (bytes > 0 && file.write(static_cast<const char*>(data), bytes) != bytes) {
        error = QStringLiteral("short write to %1").arg(path);
        return false;
    }
    return true;
}

double softUnitStep(double value) {
    return value > 0.0 ? std::exp(-1.0 / value) : 0.0;
}

float smoothFiniteRadial(double length, double radius, int index, int count) {
    const double step = radius / static_cast<double>(count + 1);
    const double center = step * static_cast<double>(index + 1);
    const double diff = (length - center) / step;
    const double basis = kSmoothFiniteNormalisation * std::exp(2.0)
                         * softUnitStep(diff + 1.0) * softUnitStep(1.0 - diff);
    const double ratio = std::min(length / radius, 1.0);
    const double cutoff = 0.5 * (std::cos(kPi * ratio) + 1.0);
    return static_cast<float>(basis * cutoff);
}

template <typename Enum>
QString enumNumber(Enum value) {
    return QString::number(static_cast<int>(value));
}

bool channelNamesMatch(const QString& kind,
                       const QStringList& observed,
                       const QStringList& expected,
                       QString& error) {
    if (observed == expected)
        return true;
    const qsizetype common = std::min(observed.size(), expected.size());
    qsizetype firstDifference = 0;
    while (firstDifference < common && observed[firstDifference] == expected[firstDifference])
        ++firstDifference;
    if (firstDifference < common) {
        error = QStringLiteral("%1 channel %2 is %3; model requires %4")
                    .arg(kind)
                    .arg(firstDifference)
                    .arg(observed[firstDifference], expected[firstDifference]);
    } else {
        error = QStringLiteral("%1 channel count is %2; model requires %3")
                    .arg(kind)
                    .arg(observed.size())
                    .arg(expected.size());
    }
    return false;
}

}  // namespace

ExperimentalShieldingMlStore::ExperimentalShieldingMlStore(
    const QtProtein* protein,
    Conformation* conformation,
    QString modelPath,
    QString runtimeManifestPath,
    QString extractionManifestPath,
    QString helperPath,
    QString device,
    QString fallbackHelperPath,
    QObject* parent)
    : QObject(parent)
    , protein_(protein)
    , conformation_(conformation)
    , modelPath_(std::move(modelPath))
    , runtimeManifestPath_(std::move(runtimeManifestPath))
    , extractionManifestPath_(std::move(extractionManifestPath))
    , helperPath_(std::move(helperPath))
    , device_(std::move(device))
    , fallbackHelperPath_(std::move(fallbackHelperPath))
    , process_(new QProcess(this)) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("ExperimentalShieldingMlStore"));

    if (!protein_ || !conformation_) {
        errorReason_ = QStringLiteral("protein or conformation is missing");
    } else if (!QFileInfo(modelPath_).isFile()) {
        errorReason_ = QStringLiteral("model is missing: %1").arg(modelPath_);
    } else if (!QFileInfo(helperPath_).isFile()) {
        errorReason_ = QStringLiteral("inference helper is missing: %1").arg(helperPath_);
    } else if (device_ != QStringLiteral("cpu") && device_ != QStringLiteral("rocm")) {
        errorReason_ = QStringLiteral("unsupported inference device: %1").arg(device_);
    } else if (!workRoot_.isValid()) {
        errorReason_ = QStringLiteral("could not create inference working directory");
    } else {
        ready_ = loadContract(runtimeManifestPath_)
                 && validateExtractionManifest(extractionManifestPath_)
                 && validateInitialFrameInputs();
    }

    ACONNECT(process_, &QProcess::finished, this, &ExperimentalShieldingMlStore::finishProcess);
    ACONNECT(process_, &QProcess::errorOccurred, this, [this](QProcess::ProcessError error) {
        if (error != QProcess::FailedToStart || !activeFrame_)
            return;
        const QString reason =
            QStringLiteral("inference helper failed to start: %1").arg(process_->errorString());
        if (!scheduleCpuFallback(reason))
            failActiveFrame(reason);
    });

    if (!ready_) {
        diagnostics::ErrorBus::Report(diagnostics::Severity::Error,
                                      QStringLiteral("ExperimentalShieldingMlStore"),
                                      errorReason_,
                                      runtimeManifestPath_);
    }
}

ExperimentalShieldingMlStore::~ExperimentalShieldingMlStore() {
    pendingFrames_.clear();
    activeFrame_.reset();
    QObject::disconnect(process_, nullptr, this, nullptr);
    process_->blockSignals(true);
    if (process_->state() != QProcess::NotRunning) {
        process_->kill();
        process_->waitForFinished();
    }
}

bool ExperimentalShieldingMlStore::ManifestHasInferenceSchema(
    const QString& manifestPath,
    QString* reason) {
    QJsonObject manifest;
    QString error;
    if (!readJsonObject(manifestPath, manifest, error)) {
        if (reason)
            *reason = error;
        return false;
    }
    const bool valid = validateRuntimeManifest(manifest, error);
    if (!valid && reason)
        *reason = error;
    return valid;
}

bool ExperimentalShieldingMlStore::isRunning() const {
    return process_->state() != QProcess::NotRunning;
}

bool ExperimentalShieldingMlStore::loadContract(const QString& manifestPath) {
    QJsonObject manifest;
    if (!readJsonObject(manifestPath, manifest, errorReason_))
        return false;
    if (!validateRuntimeManifest(manifest, errorReason_))
        return false;
    if (!validateModelArtifact(manifest, modelPath_, errorReason_))
        return false;

    const QJsonObject schema = manifest.value(QStringLiteral("inference_schema")).toObject();
    contract_.modelId = schema.value(QStringLiteral("model_id")).toString();

    const QJsonObject graph = schema.value(QStringLiteral("graph")).toObject();
    contract_.radius = graph.value(QStringLiteral("radius_angstrom")).toDouble();
    contract_.maxNeighbors = graph.value(QStringLiteral("max_neighbors")).toInt();
    contract_.radialDim = graph.value(QStringLiteral("radial_dim")).toInt();
    if (!(contract_.radius > 0.0) || !std::isfinite(contract_.radius)
        || contract_.maxNeighbors <= 0 || contract_.radialDim <= 1) {
        errorReason_ = QStringLiteral("invalid graph parameters in inference schema");
        return false;
    }

    if (!readStringList(schema.value(QStringLiteral("label_keys")),
                        QStringLiteral("label_keys"),
                        contract_.labelKeys,
                        errorReason_)) {
        return false;
    }
    if (contract_.labelKeys.size() != 15) {
        errorReason_ =
            QStringLiteral("F003 runtime requires 15 categorical labels; manifest has %1")
                .arg(contract_.labelKeys.size());
        return false;
    }
    const QJsonObject vocabularies = schema.value(QStringLiteral("label_vocabs")).toObject();
    for (const QString& key : contract_.labelKeys) {
        const QJsonObject source = vocabularies.value(key).toObject();
        if (source.isEmpty()) {
            errorReason_ = QStringLiteral("label vocabulary is absent: %1").arg(key);
            return false;
        }
        QHash<QString, std::int64_t> vocabulary;
        for (auto it = source.constBegin(); it != source.constEnd(); ++it) {
            if (!it.value().isDouble()) {
                errorReason_ = QStringLiteral("label vocabulary %1 has a non-numeric id")
                                   .arg(key);
                return false;
            }
            vocabulary.insert(it.key(), static_cast<std::int64_t>(it.value().toDouble()));
        }
        if (!vocabulary.contains(QStringLiteral("<UNK>"))) {
            errorReason_ = QStringLiteral("label vocabulary %1 has no <UNK>").arg(key);
            return false;
        }
        contract_.labelVocabs.insert(key, vocabulary);
    }

    if (!readStringList(schema.value(QStringLiteral("ring_type_order")),
                        QStringLiteral("ring_type_order"),
                        contract_.ringTypeOrder,
                        errorReason_)
        || !readIntVector(schema.value(QStringLiteral("ring_active")),
                          QStringLiteral("ring_active"),
                          contract_.ringActive,
                          errorReason_)
        || !readDoubleVector(schema.value(QStringLiteral("ring_intensity_nA_per_T")),
                             QStringLiteral("ring_intensity_nA_per_T"),
                             contract_.ringIntensity,
                             errorReason_)) {
        return false;
    }
    if (contract_.ringTypeOrder.isEmpty()
        || contract_.ringTypeOrder.size() != contract_.ringIntensity.size()
        || contract_.ringActive.isEmpty()) {
        errorReason_ = QStringLiteral("ring-current contract dimensions do not agree");
        return false;
    }
    for (const int index : contract_.ringActive) {
        if (index < 0 || index >= contract_.ringTypeOrder.size()) {
            errorReason_ = QStringLiteral("ring_active contains out-of-range index %1")
                               .arg(index);
            return false;
        }
    }

    const QJsonObject expectedCounts = schema.value(QStringLiteral("expected_channels")).toObject();
    contract_.expectedL1Channels = expectedCounts.value(QStringLiteral("l1")).toInt();
    contract_.expectedL2Channels = expectedCounts.value(QStringLiteral("l2")).toInt();
    contract_.expectedScalarChannels = expectedCounts.value(QStringLiteral("scalars")).toInt();
    contract_.expectedApplicabilityChannels =
        expectedCounts.value(QStringLiteral("applicability")).toInt();

    const QJsonObject expectedNames =
        schema.value(QStringLiteral("expected_channel_names")).toObject();
    if (!readStringList(expectedNames.value(QStringLiteral("l1")),
                        QStringLiteral("expected_channel_names.l1"),
                        contract_.expectedL1Names,
                        errorReason_)
        || !readStringList(expectedNames.value(QStringLiteral("l2")),
                           QStringLiteral("expected_channel_names.l2"),
                           contract_.expectedL2Names,
                           errorReason_)
        || !readStringList(expectedNames.value(QStringLiteral("scalar")),
                           QStringLiteral("expected_channel_names.scalar"),
                           contract_.expectedScalarNames,
                           errorReason_)
        || !readStringList(expectedNames.value(QStringLiteral("applicability")),
                           QStringLiteral("expected_channel_names.applicability"),
                           contract_.expectedApplicabilityNames,
                           errorReason_)) {
        return false;
    }
    if (contract_.expectedL1Channels != contract_.expectedL1Names.size()
        || contract_.expectedL2Channels != contract_.expectedL2Names.size()
        || contract_.expectedScalarChannels != contract_.expectedScalarNames.size()
        || contract_.expectedApplicabilityChannels
               != contract_.expectedApplicabilityNames.size()) {
        errorReason_ =
            QStringLiteral("expected channel counts do not match expected channel-name lists");
        return false;
    }

    const auto parseAxis = [](const QString& text,
                              FeatureAxis& out) {
        if (text == QStringLiteral("atom")) {
            out = FeatureAxis::Atom;
            return true;
        }
        if (text == QStringLiteral("residue")) {
            out = FeatureAxis::Residue;
            return true;
        }
        return false;
    };
    const auto parseKind = [](const QString& text,
                              FeatureKind& out) {
        if (text == QStringLiteral("scalar")) {
            out = FeatureKind::Scalar;
            return true;
        }
        if (text == QStringLiteral("l1")) {
            out = FeatureKind::L1;
            return true;
        }
        if (text == QStringLiteral("l2")) {
            out = FeatureKind::L2;
            return true;
        }
        return false;
    };
    const auto parseLayout = [](const QString& text,
                                FeatureLayout& out) {
        if (text == QStringLiteral("scalar")) out = FeatureLayout::Scalar;
        else if (text == QStringLiteral("scalar_cols")) out = FeatureLayout::ScalarColumns;
        else if (text == QStringLiteral("vector")) out = FeatureLayout::Vector;
        else if (text == QStringLiteral("t2")) out = FeatureLayout::T2;
        else if (text == QStringLiteral("full9_t0")) out = FeatureLayout::Full9T0;
        else if (text == QStringLiteral("full9_t2")) out = FeatureLayout::Full9T2;
        else if (text == QStringLiteral("ring_t0")) out = FeatureLayout::RingT0;
        else if (text == QStringLiteral("ring_t1")) out = FeatureLayout::RingT1;
        else if (text == QStringLiteral("ring_t2")) out = FeatureLayout::RingT2;
        else return false;
        return true;
    };

    const QJsonArray featureValues = schema.value(QStringLiteral("numeric_features")).toArray();
    if (featureValues.isEmpty()) {
        errorReason_ = QStringLiteral("numeric_features is empty");
        return false;
    }
    for (const QJsonValue featureValue : featureValues) {
        if (!featureValue.isObject()) {
            errorReason_ = QStringLiteral("numeric_features contains a non-object");
            return false;
        }
        const QJsonObject item = featureValue.toObject();
        FeatureSpec spec;
        spec.key = item.value(QStringLiteral("key")).toString();
        if (spec.key.isEmpty()
            || !parseAxis(item.value(QStringLiteral("axis")).toString(), spec.axis)
            || !parseKind(item.value(QStringLiteral("kind")).toString(), spec.kind)
            || !parseLayout(item.value(QStringLiteral("layout")).toString(), spec.layout)) {
            errorReason_ =
                QStringLiteral("invalid numeric feature declaration: %1").arg(spec.key);
            return false;
        }
        const QString validity = item.value(QStringLiteral("validity")).toString();
        if (validity == QStringLiteral("required"))
            spec.validity = FeatureValidity::Required;
        else if (validity == QStringLiteral("external_mask"))
            spec.validity = FeatureValidity::ExternalMask;
        else {
            errorReason_ = QStringLiteral("%1 uses unsupported validity policy %2")
                               .arg(spec.key, validity);
            return false;
        }
        const QString scale = item.value(QStringLiteral("scale")).toString();
        if (scale == QStringLiteral("none"))
            spec.scale = FeatureScale::None;
        else if (scale == QStringLiteral("manifest_bs_intensity"))
            spec.scale = FeatureScale::ManifestBsIntensity;
        else {
            errorReason_ =
                QStringLiteral("%1 uses unsupported scale %2").arg(spec.key, scale);
            return false;
        }
        spec.center = item.value(QStringLiteral("center")).toBool();
        spec.emitMask = item.value(QStringLiteral("emit_mask")).toBool();
        if (!readIntVector(item.value(QStringLiteral("columns")),
                           spec.key + QStringLiteral(".columns"),
                           spec.columns,
                           errorReason_)
            || !readIntVector(item.value(QStringLiteral("mask_columns")),
                              spec.key + QStringLiteral(".mask_columns"),
                              spec.maskColumns,
                              errorReason_)
            || !readIntVector(item.value(QStringLiteral("emit_mask_columns")),
                              spec.key + QStringLiteral(".emit_mask_columns"),
                              spec.emitMaskColumns,
                              errorReason_)) {
            return false;
        }

        QStringList sourceFiles;
        if (!readStringList(item.value(QStringLiteral("sources")),
                            spec.key + QStringLiteral(".sources"),
                            sourceFiles,
                            errorReason_)
            || sourceFiles.isEmpty()) {
            if (errorReason_.isEmpty())
                errorReason_ = QStringLiteral("%1 has no input sources").arg(spec.key);
            return false;
        }
        const io::NativeAxis expectedAxis =
            spec.axis == FeatureAxis::Atom ? io::NativeAxis::Atom : io::NativeAxis::Residue;
        for (const QString& fileName : sourceFiles) {
            FeatureSource source;
            source.fileName = fileName;
            if (!npyStem(fileName, source.stem, errorReason_))
                return false;
            const QByteArray stemBytes = source.stem.toUtf8();
            const std::optional<io::FieldKind> field =
                io::FindFieldByStem(std::string_view(stemBytes.constData(),
                                                     static_cast<std::size_t>(stemBytes.size())));
            if (!field) {
                errorReason_ = QStringLiteral("%1 names unknown catalog field %2")
                                   .arg(spec.key, source.stem);
                return false;
            }
            if (!io::ShouldLoadFrameField(*field)) {
                errorReason_ =
                    QStringLiteral("%1 input %2 is excluded from Reader snapshots")
                        .arg(spec.key, source.stem);
                return false;
            }
            if (io::FieldSpecFor(*field).axis != expectedAxis) {
                errorReason_ = QStringLiteral("%1 input %2 has the wrong catalog axis")
                                   .arg(spec.key, source.stem);
                return false;
            }
            source.field = *field;
            spec.sources.push_back(source);
        }

        spec.maskFileName = item.value(QStringLiteral("mask_source")).toString();
        if (spec.validity == FeatureValidity::ExternalMask) {
            if (!npyStem(spec.maskFileName, spec.maskStem, errorReason_))
                return false;
            const QByteArray stemBytes = spec.maskStem.toUtf8();
            spec.maskField =
                io::FindFieldByStem(std::string_view(stemBytes.constData(),
                                                     static_cast<std::size_t>(stemBytes.size())));
            if (!spec.maskField || !io::ShouldLoadFrameField(*spec.maskField)
                || io::FieldSpecFor(*spec.maskField).axis != expectedAxis) {
                errorReason_ = QStringLiteral("%1 mask %2 is not a loadable %3-axis field")
                                   .arg(spec.key,
                                        spec.maskStem,
                                        spec.axis == FeatureAxis::Atom
                                            ? QStringLiteral("atom")
                                            : QStringLiteral("residue"));
                return false;
            }
        } else if (!spec.maskFileName.isEmpty()) {
            errorReason_ =
                QStringLiteral("%1 declares a mask without external_mask validity")
                    .arg(spec.key);
            return false;
        }

        const bool kindLayoutMatches =
            (spec.kind == FeatureKind::Scalar
             && (spec.layout == FeatureLayout::Scalar
                 || spec.layout == FeatureLayout::ScalarColumns
                 || spec.layout == FeatureLayout::Full9T0
                 || spec.layout == FeatureLayout::RingT0))
            || (spec.kind == FeatureKind::L1
                && (spec.layout == FeatureLayout::Vector
                    || spec.layout == FeatureLayout::RingT1))
            || (spec.kind == FeatureKind::L2
                && (spec.layout == FeatureLayout::T2
                    || spec.layout == FeatureLayout::Full9T2
                    || spec.layout == FeatureLayout::RingT2));
        if (!kindLayoutMatches) {
            errorReason_ = QStringLiteral("%1 kind and projection layout disagree")
                               .arg(spec.key);
            return false;
        }
        if (spec.scale == FeatureScale::ManifestBsIntensity
            && spec.layout != FeatureLayout::RingT0
            && spec.layout != FeatureLayout::RingT1
            && spec.layout != FeatureLayout::RingT2) {
            errorReason_ = QStringLiteral("%1 applies ring intensity to a non-ring layout")
                               .arg(spec.key);
            return false;
        }
        contract_.features.push_back(std::move(spec));
    }
    return true;
}

bool ExperimentalShieldingMlStore::validateExtractionManifest(
    const QString& extractionManifestPath) {
    QJsonObject manifest;
    if (!readJsonObject(extractionManifestPath, manifest, errorReason_))
        return false;
    if (manifest.value(QStringLiteral("schema_version")).toString()
        != QStringLiteral("1.0")) {
        errorReason_ = QStringLiteral("unsupported extraction manifest schema in %1")
                           .arg(extractionManifestPath);
        return false;
    }
    const QJsonObject ring =
        manifest.value(QStringLiteral("feature_metadata"))
            .toObject()
            .value(QStringLiteral("ring_current"))
            .toObject();
    QStringList order;
    QVector<double> intensity;
    if (!readStringList(ring.value(QStringLiteral("ring_type_order")),
                        QStringLiteral("extraction ring_type_order"),
                        order,
                        errorReason_)
        || !readDoubleVector(ring.value(QStringLiteral("resolved_intensity_nA_per_T")),
                             QStringLiteral("extraction ring intensities"),
                             intensity,
                             errorReason_)) {
        return false;
    }
    if (order != contract_.ringTypeOrder || intensity != contract_.ringIntensity) {
        errorReason_ =
            QStringLiteral("extraction ring-current order or intensity differs from %1")
                .arg(contract_.modelId);
        return false;
    }
    const QString kernel = ring.value(QStringLiteral("kernel_contract")).toString();
    if (!kernel.contains(QStringLiteral("unit-current"), Qt::CaseSensitive)
        || !kernel.contains(QStringLiteral("multiply"), Qt::CaseSensitive)) {
        errorReason_ = QStringLiteral("unrecognised extraction ring-current kernel contract");
        return false;
    }
    return true;
}

bool ExperimentalShieldingMlStore::validateInitialFrameInputs() {
    if (conformation_->frameCount() == 0) {
        errorReason_ = QStringLiteral("conformation has no frames");
        return false;
    }

    conformation_->requestSnapshot(0);
    const std::shared_ptr<const QtConformationSnapshot> snapshot =
        conformation_->snapshot(0);
    if (!snapshot) {
        errorReason_ = QStringLiteral("frame 0 model input snapshot is unavailable");
        return false;
    }
    if (!snapshot->has(io::FieldKind::Pos)) {
        errorReason_ = QStringLiteral("frame 0 required model input pos.npy is absent");
        return false;
    }

    for (const FeatureSpec& spec : contract_.features) {
        for (const FeatureSource& source : spec.sources) {
            if (!snapshot->has(source.field)) {
                errorReason_ =
                    QStringLiteral("frame 0 required model input %1 is absent")
                        .arg(source.fileName);
                return false;
            }
        }
        if (spec.validity == FeatureValidity::ExternalMask
            && (!spec.maskField || !snapshot->has(*spec.maskField))) {
            errorReason_ =
                QStringLiteral("frame 0 required validity mask %1 is absent")
                    .arg(spec.maskFileName);
            return false;
        }
    }
    return true;
}

bool ExperimentalShieldingMlStore::buildInput(
    std::size_t frame,
    const QtConformationSnapshot& snapshot,
    const QString& inputDir,
    QString& error) const {
    struct ProjectedBlock {
        int channels = 0;
        int width = 0;
        std::vector<float> values;
        std::vector<float> valid;
        QStringList names;
        int emittedMaskChannels = 0;
        std::vector<float> emittedMasks;
        QStringList emittedMaskNames;
    };

    const std::size_t atomCount = protein_->atomCount();
    const std::size_t residueCount = protein_->residueCount();
    const QString frameLabel = QStringLiteral("frame %1").arg(frame);
    if (atomCount == 0 || residueCount == 0) {
        error = QStringLiteral("protein has no atoms or residues");
        return false;
    }
    QDir dir;
    if (!dir.mkpath(inputDir)) {
        error = QStringLiteral("could not create inference input directory");
        return false;
    }

    std::vector<std::size_t> atomResidues(atomCount);
    for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
        const int residueIndex = protein_->atom(atomIndex).residueIndex;
        if (residueIndex < 0
            || static_cast<std::size_t>(residueIndex) >= residueCount) {
            error = QStringLiteral("%1 atom %2 has out-of-range residue index %3")
                        .arg(frameLabel)
                        .arg(atomIndex)
                        .arg(residueIndex);
            return false;
        }
        atomResidues[atomIndex] = static_cast<std::size_t>(residueIndex);
    }

    if (!snapshot.has(io::FieldKind::Pos)) {
        error = QStringLiteral("%1 required model input pos.npy is absent").arg(frameLabel);
        return false;
    }
    const NpyColumn& positionColumn = snapshot.column(io::FieldKind::Pos);
    if (positionColumn.rows != static_cast<int>(atomCount)
        || positionColumn.cols != 3
        || positionColumn.data.size() != atomCount * 3) {
        error = QStringLiteral("%1 pos.npy shape is %2x%3; expected %4x3")
                    .arg(frameLabel)
                    .arg(positionColumn.rows)
                    .arg(positionColumn.cols)
                    .arg(atomCount);
        return false;
    }
    std::array<double, 3> centroid{0.0, 0.0, 0.0};
    for (std::size_t atom = 0; atom < atomCount; ++atom) {
        const double* row = positionColumn.row(atom);
        for (int axis = 0; axis < 3; ++axis) {
            if (!std::isfinite(row[axis])) {
                error = QStringLiteral("%1 pos.npy contains a non-finite value")
                            .arg(frameLabel);
                return false;
            }
            centroid[static_cast<std::size_t>(axis)] += row[axis];
        }
    }
    for (double& value : centroid)
        value /= static_cast<double>(atomCount);

    std::vector<float> positions(atomCount * 3);
    for (std::size_t atom = 0; atom < atomCount; ++atom) {
        const double* row = positionColumn.row(atom);
        for (int axis = 0; axis < 3; ++axis) {
            const float value = static_cast<float>(
                row[axis] - centroid[static_cast<std::size_t>(axis)]);
            if (!std::isfinite(value)) {
                error = QStringLiteral("%1 centered position overflows float32")
                            .arg(frameLabel);
                return false;
            }
            positions[atom * 3 + static_cast<std::size_t>(axis)] = value;
        }
    }

    std::vector<ProjectedBlock> scalarBlocks;
    std::vector<ProjectedBlock> l1Blocks;
    std::vector<ProjectedBlock> l2Blocks;
    std::vector<ProjectedBlock> applicabilityBlocks;

    for (const FeatureSpec& spec : contract_.features) {
        const std::size_t sourceRows =
            spec.axis == FeatureAxis::Atom ? atomCount : residueCount;
        const int width = spec.kind == FeatureKind::Scalar
                              ? 1
                              : (spec.kind == FeatureKind::L1 ? 3 : 5);

        struct SourceShape {
            const FeatureSource* source = nullptr;
            const NpyColumn* column = nullptr;
            int channels = 0;
        };
        std::vector<SourceShape> sources;
        int totalChannels = 0;
        for (const FeatureSource& source : spec.sources) {
            if (!snapshot.has(source.field)) {
                error = QStringLiteral("%1 required model input %2 is absent")
                            .arg(frameLabel, source.fileName);
                return false;
            }
            const NpyColumn& column = snapshot.column(source.field);
            if (column.rows != static_cast<int>(sourceRows)
                || column.cols <= 0
                || column.data.size()
                       != sourceRows * static_cast<std::size_t>(column.cols)) {
                error = QStringLiteral("%1/%2 has invalid shape %3x%4")
                            .arg(frameLabel, source.fileName)
                            .arg(column.rows)
                            .arg(column.cols);
                return false;
            }

            int channels = 0;
            switch (spec.layout) {
                case FeatureLayout::Scalar:
                    channels = column.cols;
                    break;
                case FeatureLayout::ScalarColumns:
                    if (spec.columns.isEmpty()) {
                        error = QStringLiteral("%1 has an empty scalar column selection")
                                    .arg(spec.key);
                        return false;
                    }
                    for (const int selected : spec.columns) {
                        if (selected < 0 || selected >= column.cols) {
                            error = QStringLiteral("%1 selects column %2 from %3 columns")
                                        .arg(spec.key)
                                        .arg(selected)
                                        .arg(column.cols);
                            return false;
                        }
                    }
                    if (!checkedSizeToInt(
                            spec.columns.size(),
                            QStringLiteral("%1 scalar columns").arg(spec.key),
                            channels,
                            error)) {
                        return false;
                    }
                    break;
                case FeatureLayout::Vector:
                    if (column.cols != 3) {
                        error = QStringLiteral("%1/%2 requires three vector columns")
                                    .arg(frameLabel, source.fileName);
                        return false;
                    }
                    channels = 1;
                    break;
                case FeatureLayout::T2:
                    if (column.cols != 5) {
                        error = QStringLiteral("%1/%2 requires five T2 columns")
                                    .arg(frameLabel, source.fileName);
                        return false;
                    }
                    channels = 1;
                    break;
                case FeatureLayout::Full9T0:
                case FeatureLayout::Full9T2:
                    if (column.cols != 9) {
                        error = QStringLiteral("%1/%2 requires nine spherical-tensor columns")
                                    .arg(frameLabel, source.fileName);
                        return false;
                    }
                    channels = 1;
                    break;
                case FeatureLayout::RingT0:
                    if (column.cols != contract_.ringTypeOrder.size()) {
                        error = QStringLiteral("%1/%2 requires %3 ring-type columns")
                                    .arg(frameLabel, source.fileName)
                                    .arg(contract_.ringTypeOrder.size());
                        return false;
                    }
                    if (!checkedSizeToInt(
                            contract_.ringActive.size(),
                            QStringLiteral("%1 active ring types").arg(spec.key),
                            channels,
                            error)) {
                        return false;
                    }
                    break;
                case FeatureLayout::RingT1:
                case FeatureLayout::RingT2: {
                    const int ringWidth =
                        spec.layout == FeatureLayout::RingT1 ? 3 : 5;
                    if (column.cols != contract_.ringTypeOrder.size() * ringWidth) {
                        error = QStringLiteral("%1/%2 has %3 columns; expected %4 ring tensor columns")
                                    .arg(frameLabel, source.fileName)
                                    .arg(column.cols)
                                    .arg(contract_.ringTypeOrder.size() * ringWidth);
                        return false;
                    }
                    if (!checkedSizeToInt(
                            contract_.ringActive.size(),
                            QStringLiteral("%1 active ring types").arg(spec.key),
                            channels,
                            error)) {
                        return false;
                    }
                    break;
                }
            }
            sources.push_back({&source, &column, channels});
            if (channels > std::numeric_limits<int>::max() - totalChannels) {
                error = QStringLiteral("%1 channel count exceeds the reader limit")
                            .arg(spec.key);
                return false;
            }
            totalChannels += channels;
        }
        if (totalChannels <= 0) {
            error = QStringLiteral("%1 projected no channels").arg(spec.key);
            return false;
        }

        std::vector<float> sourceValues(
            sourceRows * static_cast<std::size_t>(totalChannels)
                * static_cast<std::size_t>(width),
            0.0F);
        std::vector<std::uint8_t> finite(
            sourceRows * static_cast<std::size_t>(totalChannels),
            1);
        QStringList names;
        int channelOffset = 0;
        for (const SourceShape& sourceShape : sources) {
            const FeatureSource& source = *sourceShape.source;
            const NpyColumn& column = *sourceShape.column;
            for (int channel = 0; channel < sourceShape.channels; ++channel) {
                switch (spec.layout) {
                    case FeatureLayout::Scalar:
                        names.append(column.cols == 1
                                         ? QStringLiteral("%1:%2").arg(spec.key, source.stem)
                                         : QStringLiteral("%1:%2:%3")
                                               .arg(spec.key, source.stem)
                                               .arg(channel));
                        break;
                    case FeatureLayout::ScalarColumns:
                        names.append(QStringLiteral("%1:%2:col%3")
                                         .arg(spec.key, source.stem)
                                         .arg(spec.columns[channel]));
                        break;
                    case FeatureLayout::Vector:
                    case FeatureLayout::T2:
                        names.append(QStringLiteral("%1:%2").arg(spec.key, source.stem));
                        break;
                    case FeatureLayout::Full9T0:
                        names.append(QStringLiteral("%1:%2:T0").arg(spec.key, source.stem));
                        break;
                    case FeatureLayout::Full9T2:
                        names.append(QStringLiteral("%1:%2:T2").arg(spec.key, source.stem));
                        break;
                    case FeatureLayout::RingT0:
                    case FeatureLayout::RingT1:
                    case FeatureLayout::RingT2:
                        names.append(QStringLiteral("%1:%2")
                                         .arg(spec.key,
                                              contract_.ringTypeOrder[
                                                  contract_.ringActive[channel]]));
                        break;
                }
            }

            for (std::size_t rowIndex = 0; rowIndex < sourceRows; ++rowIndex) {
                const double* sourceRow = column.row(rowIndex);
                for (int channel = 0; channel < sourceShape.channels; ++channel) {
                    bool channelFinite = true;
                    for (int component = 0; component < width; ++component) {
                        int sourceColumn = 0;
                        switch (spec.layout) {
                            case FeatureLayout::Scalar:
                                sourceColumn = channel;
                                break;
                            case FeatureLayout::ScalarColumns:
                                sourceColumn = spec.columns[channel];
                                break;
                            case FeatureLayout::Vector:
                            case FeatureLayout::T2:
                                sourceColumn = component;
                                break;
                            case FeatureLayout::Full9T0:
                                sourceColumn = 0;
                                break;
                            case FeatureLayout::Full9T2:
                                sourceColumn = 4 + component;
                                break;
                            case FeatureLayout::RingT0:
                                sourceColumn = contract_.ringActive[channel];
                                break;
                            case FeatureLayout::RingT1:
                            case FeatureLayout::RingT2:
                                sourceColumn =
                                    contract_.ringActive[channel] * width + component;
                                break;
                        }
                        double value = sourceRow[sourceColumn];
                        if (spec.scale == FeatureScale::ManifestBsIntensity) {
                            value *= contract_.ringIntensity[
                                contract_.ringActive[channel]];
                        }
                        const float projected = static_cast<float>(value);
                        if (!std::isfinite(value) || !std::isfinite(projected))
                            channelFinite = false;
                        const std::size_t target =
                            ((rowIndex * static_cast<std::size_t>(totalChannels)
                              + static_cast<std::size_t>(channelOffset + channel))
                             * static_cast<std::size_t>(width))
                            + static_cast<std::size_t>(component);
                        sourceValues[target] = projected;
                    }
                    finite[rowIndex * static_cast<std::size_t>(totalChannels)
                           + static_cast<std::size_t>(channelOffset + channel)] =
                        channelFinite ? 1 : 0;
                }
            }
            channelOffset += sourceShape.channels;
        }

        std::vector<std::uint8_t> valid(
            sourceRows * static_cast<std::size_t>(totalChannels),
            1);
        std::vector<std::uint8_t> rawMask;
        int rawMaskChannels = 0;
        if (spec.validity == FeatureValidity::Required) {
            const std::size_t invalidCount =
                static_cast<std::size_t>(std::count(finite.cbegin(), finite.cend(), 0));
            if (invalidCount != 0) {
                error = QStringLiteral("%1 %2 has %3 undeclared non-finite rows")
                            .arg(frameLabel, spec.key)
                            .arg(invalidCount);
                return false;
            }
        } else {
            if (!spec.maskField || !snapshot.has(*spec.maskField)) {
                error = QStringLiteral("%1 required validity mask %2 is absent")
                            .arg(frameLabel, spec.maskFileName);
                return false;
            }
            const NpyColumn& mask = snapshot.column(*spec.maskField);
            if (mask.rows != static_cast<int>(sourceRows)
                || mask.cols <= 0
                || mask.data.size()
                       != sourceRows * static_cast<std::size_t>(mask.cols)) {
                error = QStringLiteral("%1/%2 has invalid shape %3x%4")
                            .arg(frameLabel, spec.maskFileName)
                            .arg(mask.rows)
                            .arg(mask.cols);
                return false;
            }
            QVector<int> maskColumns = spec.maskColumns;
            if (maskColumns.isEmpty()) {
                maskColumns.reserve(mask.cols);
                for (int column = 0; column < mask.cols; ++column)
                    maskColumns.push_back(column);
            }
            for (const int selected : maskColumns) {
                if (selected < 0 || selected >= mask.cols) {
                    error = QStringLiteral("%1 mask selects out-of-range column %2")
                                .arg(spec.key)
                                .arg(selected);
                    return false;
                }
            }
            if (!checkedSizeToInt(
                    maskColumns.size(),
                    QStringLiteral("%1 mask columns").arg(spec.key),
                    rawMaskChannels,
                    error)) {
                return false;
            }
            if (rawMaskChannels != 1 && rawMaskChannels != totalChannels) {
                error = QStringLiteral("%1 mask has %2 channels and cannot cover %3")
                            .arg(spec.key)
                            .arg(rawMaskChannels)
                            .arg(totalChannels);
                return false;
            }
            rawMask.resize(sourceRows * static_cast<std::size_t>(rawMaskChannels));
            for (std::size_t rowIndex = 0; rowIndex < sourceRows; ++rowIndex) {
                const double* maskRow = mask.row(rowIndex);
                for (int channel = 0; channel < rawMaskChannels; ++channel) {
                    const double value = maskRow[maskColumns[channel]];
                    if (!std::isfinite(value) || (value != 0.0 && value != 1.0)) {
                        error = QStringLiteral("%1/%2 is not a finite binary mask")
                                    .arg(frameLabel, spec.maskFileName);
                        return false;
                    }
                    rawMask[rowIndex * static_cast<std::size_t>(rawMaskChannels)
                            + static_cast<std::size_t>(channel)] =
                        value == 1.0 ? 1 : 0;
                }
            }
            for (std::size_t rowIndex = 0; rowIndex < sourceRows; ++rowIndex) {
                for (int channel = 0; channel < totalChannels; ++channel) {
                    const std::uint8_t isValid =
                        rawMask[rowIndex * static_cast<std::size_t>(rawMaskChannels)
                                + static_cast<std::size_t>(
                                    rawMaskChannels == 1 ? 0 : channel)];
                    valid[rowIndex * static_cast<std::size_t>(totalChannels)
                          + static_cast<std::size_t>(channel)] = isValid;
                    if (isValid
                        && !finite[rowIndex * static_cast<std::size_t>(totalChannels)
                                   + static_cast<std::size_t>(channel)]) {
                        error = QStringLiteral("%1 %2 is non-finite where %3 says valid")
                                    .arg(frameLabel, spec.key, spec.maskFileName);
                        return false;
                    }
                }
            }
        }

        for (std::size_t rowIndex = 0; rowIndex < sourceRows; ++rowIndex) {
            for (int channel = 0; channel < totalChannels; ++channel) {
                if (valid[rowIndex * static_cast<std::size_t>(totalChannels)
                          + static_cast<std::size_t>(channel)]) {
                    continue;
                }
                const std::size_t base =
                    (rowIndex * static_cast<std::size_t>(totalChannels)
                     + static_cast<std::size_t>(channel))
                    * static_cast<std::size_t>(width);
                std::fill_n(sourceValues.begin() + static_cast<std::ptrdiff_t>(base),
                            width,
                            0.0F);
            }
        }
        if (!std::all_of(sourceValues.cbegin(), sourceValues.cend(), [](float value) {
                return std::isfinite(value);
            })) {
            error = QStringLiteral("%1 %2 remains non-finite after masking")
                        .arg(frameLabel, spec.key);
            return false;
        }

        ProjectedBlock block;
        block.channels = totalChannels;
        block.width = width;
        block.names = names;
        block.values.resize(atomCount * static_cast<std::size_t>(totalChannels)
                            * static_cast<std::size_t>(width));
        block.valid.resize(atomCount * static_cast<std::size_t>(totalChannels));
        for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
            const std::size_t sourceRow =
                spec.axis == FeatureAxis::Atom ? atomIndex : atomResidues[atomIndex];
            for (int channel = 0; channel < totalChannels; ++channel) {
                const std::size_t sourceChannel =
                    sourceRow * static_cast<std::size_t>(totalChannels)
                    + static_cast<std::size_t>(channel);
                const std::size_t targetChannel =
                    atomIndex * static_cast<std::size_t>(totalChannels)
                    + static_cast<std::size_t>(channel);
                block.valid[targetChannel] = valid[sourceChannel] ? 1.0F : 0.0F;
                for (int component = 0; component < width; ++component) {
                    block.values[targetChannel * static_cast<std::size_t>(width)
                                 + static_cast<std::size_t>(component)] =
                        sourceValues[sourceChannel * static_cast<std::size_t>(width)
                                     + static_cast<std::size_t>(component)];
                }
            }
        }

        if (spec.emitMask) {
            QVector<int> emittedColumns = spec.emitMaskColumns;
            if (emittedColumns.isEmpty()) {
                emittedColumns.reserve(rawMaskChannels);
                for (int column = 0; column < rawMaskChannels; ++column)
                    emittedColumns.push_back(column);
            }
            if (rawMaskChannels == 0) {
                error = QStringLiteral("%1 requests an emitted mask without an external mask")
                            .arg(spec.key);
                return false;
            }
            for (const int selected : emittedColumns) {
                if (selected < 0 || selected >= rawMaskChannels) {
                    error = QStringLiteral("%1 emits out-of-range mask column %2")
                                .arg(spec.key)
                                .arg(selected);
                    return false;
                }
            }
            if (!checkedSizeToInt(
                    emittedColumns.size(),
                    QStringLiteral("%1 emitted mask columns").arg(spec.key),
                    block.emittedMaskChannels,
                    error)) {
                return false;
            }
            block.emittedMasks.resize(
                atomCount * static_cast<std::size_t>(block.emittedMaskChannels));
            for (const int selected : emittedColumns) {
                block.emittedMaskNames.append(
                    QStringLiteral("valid:%1:%2").arg(spec.key).arg(selected));
            }
            for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
                const std::size_t sourceRow =
                    spec.axis == FeatureAxis::Atom ? atomIndex : atomResidues[atomIndex];
                for (int channel = 0; channel < block.emittedMaskChannels; ++channel) {
                    block.emittedMasks[
                        atomIndex * static_cast<std::size_t>(block.emittedMaskChannels)
                        + static_cast<std::size_t>(channel)] =
                        rawMask[sourceRow * static_cast<std::size_t>(rawMaskChannels)
                                + static_cast<std::size_t>(emittedColumns[channel])]
                            ? 1.0F
                            : 0.0F;
                }
            }
            applicabilityBlocks.push_back(block);
        }

        if (spec.kind == FeatureKind::Scalar)
            scalarBlocks.push_back(std::move(block));
        else if (spec.kind == FeatureKind::L1)
            l1Blocks.push_back(std::move(block));
        else
            l2Blocks.push_back(std::move(block));
    }

    const auto concatenateBlocks =
        [atomCount, &error](const QString& kind,
                            const std::vector<ProjectedBlock>& blocks,
                            int width,
                            int expectedChannels,
                            const QStringList& expectedNames,
                            std::vector<float>& values,
                            std::vector<float>& valid) {
            int channels = 0;
            QStringList names;
            for (const ProjectedBlock& block : blocks) {
                if (block.width != width) {
                    error = QStringLiteral("%1 block width disagrees").arg(kind);
                    return false;
                }
                channels += block.channels;
                names.append(block.names);
            }
            if (channels != expectedChannels
                || !channelNamesMatch(kind, names, expectedNames, error)) {
                if (error.isEmpty()) {
                    error = QStringLiteral("%1 channel count is %2; expected %3")
                                .arg(kind)
                                .arg(channels)
                                .arg(expectedChannels);
                }
                return false;
            }
            values.resize(atomCount * static_cast<std::size_t>(channels)
                          * static_cast<std::size_t>(width));
            valid.resize(atomCount * static_cast<std::size_t>(channels));
            for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
                int channelOffset = 0;
                for (const ProjectedBlock& block : blocks) {
                    for (int channel = 0; channel < block.channels; ++channel) {
                        const std::size_t sourceChannel =
                            atomIndex * static_cast<std::size_t>(block.channels)
                            + static_cast<std::size_t>(channel);
                        const std::size_t targetChannel =
                            atomIndex * static_cast<std::size_t>(channels)
                            + static_cast<std::size_t>(channelOffset + channel);
                        valid[targetChannel] = block.valid[sourceChannel];
                        for (int component = 0; component < width; ++component) {
                            values[targetChannel * static_cast<std::size_t>(width)
                                   + static_cast<std::size_t>(component)] =
                                block.values[
                                    sourceChannel * static_cast<std::size_t>(width)
                                    + static_cast<std::size_t>(component)];
                        }
                    }
                    channelOffset += block.channels;
                }
            }
            return true;
        };

    std::vector<float> scalars;
    std::vector<float> scalarValid;
    std::vector<float> l1;
    std::vector<float> l1Valid;
    std::vector<float> l2;
    std::vector<float> l2Valid;
    if (!concatenateBlocks(QStringLiteral("scalar"),
                           scalarBlocks,
                           1,
                           contract_.expectedScalarChannels,
                           contract_.expectedScalarNames,
                           scalars,
                           scalarValid)
        || !concatenateBlocks(QStringLiteral("l1"),
                              l1Blocks,
                              3,
                              contract_.expectedL1Channels,
                              contract_.expectedL1Names,
                              l1,
                              l1Valid)
        || !concatenateBlocks(QStringLiteral("l2"),
                              l2Blocks,
                              5,
                              contract_.expectedL2Channels,
                              contract_.expectedL2Names,
                              l2,
                              l2Valid)) {
        return false;
    }

    int applicabilityChannels = 0;
    QStringList applicabilityNames;
    for (const ProjectedBlock& block : applicabilityBlocks) {
        applicabilityChannels += block.emittedMaskChannels;
        applicabilityNames.append(block.emittedMaskNames);
    }
    if (applicabilityChannels != contract_.expectedApplicabilityChannels
        || !channelNamesMatch(QStringLiteral("applicability"),
                              applicabilityNames,
                              contract_.expectedApplicabilityNames,
                              error)) {
        if (error.isEmpty()) {
            error = QStringLiteral("applicability channel count is %1; expected %2")
                        .arg(applicabilityChannels)
                        .arg(contract_.expectedApplicabilityChannels);
        }
        return false;
    }
    std::vector<float> applicability(
        atomCount * static_cast<std::size_t>(applicabilityChannels));
    for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
        int channelOffset = 0;
        for (const ProjectedBlock& block : applicabilityBlocks) {
            for (int channel = 0; channel < block.emittedMaskChannels; ++channel) {
                applicability[
                    atomIndex * static_cast<std::size_t>(applicabilityChannels)
                    + static_cast<std::size_t>(channelOffset + channel)] =
                    block.emittedMasks[
                        atomIndex * static_cast<std::size_t>(block.emittedMaskChannels)
                        + static_cast<std::size_t>(channel)];
            }
            channelOffset += block.emittedMaskChannels;
        }
    }

    std::vector<std::int64_t> labelIds(
        atomCount * static_cast<std::size_t>(contract_.labelKeys.size()));
    QHash<QString, int> unknownCounts;
    for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
        const QtAtom& atom = protein_->atom(atomIndex);
        const std::size_t residueIndex = atomResidues[atomIndex];
        const QtResidue& residue = protein_->residue(residueIndex);
        const QtResidueNames& residueNames = protein_->residueNames(residueIndex);
        const QtAtomNames& atomNames = protein_->atomNames(atomIndex);
        for (int labelIndex = 0; labelIndex < contract_.labelKeys.size(); ++labelIndex) {
            const QString& key = contract_.labelKeys[labelIndex];
            QString value;
            if (key == QStringLiteral("element"))
                value = QString::number(atom.AtomicNumber());
            else if (key == QStringLiteral("iupac_atom"))
                value = atomNames.iupac;
            else if (key == QStringLiteral("iupac_residue"))
                value = residueNames.iupac;
            else if (key == QStringLiteral("amber_residue"))
                value = residueNames.amber;
            else if (key == QStringLiteral("terminal_state"))
                value = enumNumber(residue.terminalState);
            else if (key == QStringLiteral("residue_variant"))
                value = QString::number(residue.protonationVariantIndex);
            else if (key == QStringLiteral("formal_charge"))
                value = QString::number(atom.formalCharge);
            else if (key == QStringLiteral("polar_h_kind"))
                value = enumNumber(atom.polarH);
            else if (key == QStringLiteral("backbone_role"))
                value = enumNumber(atom.backboneRole);
            else if (key == QStringLiteral("planar_group"))
                value = enumNumber(atom.planarGroup);
            else if (key == QStringLiteral("locant"))
                value = enumNumber(atom.locant);
            else if (key == QStringLiteral("ring_position"))
                value = enumNumber(atom.ringPositionPrimary);
            else if (key == QStringLiteral("aromatic"))
                value = atom.aromatic ? QStringLiteral("1") : QStringLiteral("0");
            else if (key == QStringLiteral("exchangeable"))
                value = atom.isExchangeable ? QStringLiteral("1") : QStringLiteral("0");
            else if (key == QStringLiteral("iupac_pair"))
                value = residueNames.iupac + QLatin1Char(':') + atomNames.iupac;
            else {
                error = QStringLiteral("unsupported categorical label key %1").arg(key);
                return false;
            }

            const auto vocabularyIt = contract_.labelVocabs.constFind(key);
            if (vocabularyIt == contract_.labelVocabs.constEnd()) {
                error = QStringLiteral("categorical vocabulary disappeared for %1").arg(key);
                return false;
            }
            const QHash<QString, std::int64_t>& vocabulary = vocabularyIt.value();
            const auto idIt = vocabulary.constFind(value);
            const std::int64_t id =
                idIt == vocabulary.constEnd()
                    ? vocabulary.value(QStringLiteral("<UNK>"))
                    : idIt.value();
            if (idIt == vocabulary.constEnd())
                unknownCounts[key] += 1;
            labelIds[atomIndex * static_cast<std::size_t>(contract_.labelKeys.size())
                     + static_cast<std::size_t>(labelIndex)] = id;
        }
    }
    for (auto it = unknownCounts.constBegin(); it != unknownCounts.constEnd(); ++it) {
        qCWarning(cExperimentalMl).noquote()
            << QStringLiteral("event=label_oov frame=%1 key=%2 rows=%3")
                   .arg(frame)
                   .arg(it.key())
                   .arg(it.value());
    }

    std::vector<std::int64_t> edgeSrc;
    std::vector<std::int64_t> edgeDst;
    std::vector<float> radial;
    struct Candidate {
        std::size_t source = 0;
        double distance = 0.0;
    };
    for (std::size_t destination = 0; destination < atomCount; ++destination) {
        std::vector<Candidate> candidates;
        const std::size_t destinationOffset = destination * 3;
        for (std::size_t source = 0; source < atomCount; ++source) {
            if (source == destination)
                continue;
            const std::size_t sourceOffset = source * 3;
            const double dx = static_cast<double>(positions[destinationOffset])
                              - static_cast<double>(positions[sourceOffset]);
            const double dy = static_cast<double>(positions[destinationOffset + 1])
                              - static_cast<double>(positions[sourceOffset + 1]);
            const double dz = static_cast<double>(positions[destinationOffset + 2])
                              - static_cast<double>(positions[sourceOffset + 2]);
            const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (distance < contract_.radius)
                candidates.push_back({source, distance});
        }
        if (candidates.size() > static_cast<std::size_t>(contract_.maxNeighbors)) {
            std::partial_sort(
                candidates.begin(),
                candidates.begin() + contract_.maxNeighbors,
                candidates.end(),
                [](const Candidate& a, const Candidate& b) {
                    if (a.distance != b.distance)
                        return a.distance < b.distance;
                    return a.source < b.source;
                });
            candidates.resize(static_cast<std::size_t>(contract_.maxNeighbors));
        }
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& a,
                                                           const Candidate& b) {
            return a.source < b.source;
        });
        for (const Candidate& candidate : candidates) {
            edgeSrc.push_back(static_cast<std::int64_t>(candidate.source));
            edgeDst.push_back(static_cast<std::int64_t>(destination));
            for (int radialIndex = 0; radialIndex < contract_.radialDim; ++radialIndex) {
                radial.push_back(smoothFiniteRadial(candidate.distance,
                                                    contract_.radius,
                                                    radialIndex,
                                                    contract_.radialDim));
            }
        }
    }

    const auto bytesFor = [&error](std::size_t count,
                                   std::size_t width,
                                   const QString& name) {
        if (count > static_cast<std::size_t>(std::numeric_limits<qsizetype>::max())
                        / width) {
            error = QStringLiteral("%1 input is too large").arg(name);
            return qsizetype{-1};
        }
        return static_cast<qsizetype>(count * width);
    };
    const qsizetype posBytes = bytesFor(positions.size(), sizeof(float), QStringLiteral("position"));
    const qsizetype l1Bytes = bytesFor(l1.size(), sizeof(float), QStringLiteral("l1"));
    const qsizetype l1ValidBytes = bytesFor(l1Valid.size(), sizeof(float), QStringLiteral("l1 validity"));
    const qsizetype l2Bytes = bytesFor(l2.size(), sizeof(float), QStringLiteral("l2"));
    const qsizetype l2ValidBytes = bytesFor(l2Valid.size(), sizeof(float), QStringLiteral("l2 validity"));
    const qsizetype scalarBytes = bytesFor(scalars.size(), sizeof(float), QStringLiteral("scalar"));
    const qsizetype scalarValidBytes = bytesFor(scalarValid.size(), sizeof(float), QStringLiteral("scalar validity"));
    const qsizetype applicabilityBytes =
        bytesFor(applicability.size(), sizeof(float), QStringLiteral("applicability"));
    const qsizetype labelBytes = bytesFor(labelIds.size(), sizeof(std::int64_t), QStringLiteral("label"));
    const qsizetype edgeBytes = bytesFor(edgeSrc.size(), sizeof(std::int64_t), QStringLiteral("edge"));
    const qsizetype radialBytes = bytesFor(radial.size(), sizeof(float), QStringLiteral("radial"));
    if (posBytes < 0 || l1Bytes < 0 || l1ValidBytes < 0
        || l2Bytes < 0 || l2ValidBytes < 0
        || scalarBytes < 0 || scalarValidBytes < 0
        || applicabilityBytes < 0 || labelBytes < 0
        || edgeBytes < 0 || radialBytes < 0) {
        return false;
    }

    const QDir output(inputDir);
    if (!writeRaw(output.filePath(QStringLiteral("pos.bin")), positions.data(), posBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("l1.bin")), l1.data(), l1Bytes, error)
        || !writeRaw(output.filePath(QStringLiteral("l1_valid.bin")), l1Valid.data(), l1ValidBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("l2.bin")), l2.data(), l2Bytes, error)
        || !writeRaw(output.filePath(QStringLiteral("l2_valid.bin")), l2Valid.data(), l2ValidBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("scalars.bin")), scalars.data(), scalarBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("scalar_valid.bin")), scalarValid.data(), scalarValidBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("applicability.bin")), applicability.data(), applicabilityBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("label_ids.bin")), labelIds.data(), labelBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("edge_src.bin")), edgeSrc.data(), edgeBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("edge_dst.bin")), edgeDst.data(), edgeBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("radial.bin")), radial.data(), radialBytes, error)) {
        return false;
    }

    QFile shapes(output.filePath(QStringLiteral("shapes.txt")));
    if (!shapes.open(QIODevice::WriteOnly | QIODevice::Text)) {
        error = QStringLiteral("could not create shapes.txt");
        return false;
    }
    QTextStream stream(&shapes);
    stream << "N " << atomCount << '\n'
           << "E " << edgeSrc.size() << '\n'
           << "C1 " << contract_.expectedL1Channels << '\n'
           << "C2 " << contract_.expectedL2Channels << '\n'
           << "C0 " << contract_.expectedScalarChannels << '\n'
           << "A " << contract_.expectedApplicabilityChannels << '\n'
           << "label_count " << contract_.labelKeys.size() << '\n'
           << "radial_dim " << contract_.radialDim << '\n';
    return true;
}

void ExperimentalShieldingMlStore::requestFrame(std::size_t frame) {
    ASSERT_THREAD(this);
    if (!ready_ || !conformation_ || frame >= conformation_->frameCount())
        return;
    if (residentFrame_ && *residentFrame_ == frame) {
        emit frameReady(frame);
        return;
    }
    if (std::find(failedFrames_.cbegin(), failedFrames_.cend(), frame)
        != failedFrames_.cend()) {
        emit frameReady(frame);
        return;
    }
    if (activeFrame_ && *activeFrame_ == frame)
        return;
    if (std::find(pendingFrames_.cbegin(), pendingFrames_.cend(), frame)
        != pendingFrames_.cend()) {
        return;
    }
    pendingFrames_.push_back(frame);
    startPendingFrame();
}

void ExperimentalShieldingMlStore::startFrame(std::size_t frame) {
    conformation_->requestSnapshot(frame);
    const std::shared_ptr<const QtConformationSnapshot> snapshot =
        conformation_->snapshot(frame);
    if (!snapshot) {
        failedFrames_.push_back(frame);
        emit frameReady(frame);
        startPendingFrame();
        return;
    }

    activeFrame_ = frame;
    activeDir_ =
        QDir(workRoot_.path()).filePath(QStringLiteral("frame_%1").arg(frame));
    activeOutput_ = QDir(activeDir_).filePath(QStringLiteral("output.bin"));
    QString error;
    if (!buildInput(frame,
                    *snapshot,
                    QDir(activeDir_).filePath(QStringLiteral("input")),
                    error)) {
        failActiveFrame(error);
        return;
    }
    launchActiveFrame();
}

void ExperimentalShieldingMlStore::launchActiveFrame() {
    if (!activeFrame_)
        return;
    QProcessEnvironment environment = QProcessEnvironment::systemEnvironment();
    const QString runtimeDir = QFileInfo(helperPath_).absolutePath();
    QStringList searchPath{runtimeDir};
    const QString applicationDir = QCoreApplication::applicationDirPath();
    if (!applicationDir.isEmpty() && applicationDir != runtimeDir)
        searchPath.append(applicationDir);
    const QString inheritedPath = environment.value(QStringLiteral("PATH"));
    if (!inheritedPath.isEmpty())
        searchPath.append(inheritedPath);
    environment.insert(QStringLiteral("PATH"),
                       searchPath.join(QDir::listSeparator()));
    process_->setProcessEnvironment(environment);
    process_->setWorkingDirectory(runtimeDir);
    process_->setProgram(helperPath_);
    process_->setArguments(
        {modelPath_,
         QDir(activeDir_).filePath(QStringLiteral("input")),
         activeOutput_,
         device_});
    qCInfo(cExperimentalMl).noquote()
        << QStringLiteral("event=inference_start frame=%1 model=%2 device=%3")
               .arg(*activeFrame_)
               .arg(contract_.modelId, device_);
    process_->start();
}

bool ExperimentalShieldingMlStore::scheduleCpuFallback(const QString& reason) {
    if (!activeFrame_ || fallbackAttempted_ || fallbackHelperPath_.isEmpty()
        || device_ != QStringLiteral("rocm")
        || !QFileInfo(fallbackHelperPath_).isFile()) {
        return false;
    }
    fallbackAttempted_ = true;
    helperPath_ = std::move(fallbackHelperPath_);
    fallbackHelperPath_.clear();
    device_ = QStringLiteral("cpu");
    qCWarning(cExperimentalMl).noquote()
        << QStringLiteral("event=inference_fallback frame=%1 from=rocm to=cpu reason=%2")
               .arg(*activeFrame_)
               .arg(reason);
    emit runtimeChanged();
    QMetaObject::invokeMethod(this,
                              [this]() { launchActiveFrame(); },
                              Qt::QueuedConnection);
    return true;
}

void ExperimentalShieldingMlStore::finishProcess(
    int exitCode,
    QProcess::ExitStatus exitStatus) {
    ASSERT_THREAD(this);
    if (!activeFrame_)
        return;
    const std::size_t frame = *activeFrame_;
    const QString stderrText =
        QString::fromUtf8(process_->readAllStandardError()).trimmed();
    if (!stderrText.isEmpty())
        qCInfo(cExperimentalMl).noquote() << stderrText;

    if (exitStatus != QProcess::NormalExit || exitCode != 0) {
        QString reason =
            QStringLiteral("inference helper exited with code %1").arg(exitCode);
        if (!stderrText.isEmpty())
            reason += QStringLiteral(": ") + stderrText;
        if (!scheduleCpuFallback(reason))
            failActiveFrame(reason);
        return;
    }
    QFile output(activeOutput_);
    if (!output.open(QIODevice::ReadOnly)) {
        failActiveFrame(QStringLiteral("inference helper did not create output.bin"));
        return;
    }
    const QByteArray bytes = output.readAll();
    const std::size_t expectedValues = protein_->atomCount() * kOutputColumns;
    const std::size_t expectedBytes = expectedValues * sizeof(float);
    if (static_cast<std::size_t>(bytes.size()) != expectedBytes) {
        failActiveFrame(QStringLiteral("inference output has %1 bytes; expected %2")
                            .arg(bytes.size())
                            .arg(expectedBytes));
        return;
    }
    residentValues_.resize(expectedValues);
    std::memcpy(residentValues_.data(), bytes.constData(), expectedBytes);
    if (!std::all_of(residentValues_.cbegin(),
                     residentValues_.cend(),
                     [](float value) { return std::isfinite(value); })) {
        failActiveFrame(QStringLiteral("inference output contains non-finite values"));
        return;
    }
    residentFrame_ = frame;
    activeFrame_.reset();
    QDir(activeDir_).removeRecursively();
    activeDir_.clear();
    activeOutput_.clear();
    qCInfo(cExperimentalMl).noquote()
        << QStringLiteral("event=inference_ready frame=%1 atoms=%2 model=%3 device=%4")
               .arg(frame)
               .arg(protein_->atomCount())
               .arg(contract_.modelId, device_);
    emit frameReady(frame);
    startPendingFrame();
}

void ExperimentalShieldingMlStore::failActiveFrame(const QString& reason) {
    if (!activeFrame_)
        return;
    const std::size_t frame = *activeFrame_;
    failedFrames_.push_back(frame);
    activeFrame_.reset();
    if (!activeDir_.isEmpty())
        QDir(activeDir_).removeRecursively();
    activeDir_.clear();
    activeOutput_.clear();
    diagnostics::ErrorBus::Report(diagnostics::Severity::Error,
                                  QStringLiteral("ExperimentalShieldingMlStore"),
                                  reason,
                                  modelPath_);
    emit frameReady(frame);
    startPendingFrame();
}

void ExperimentalShieldingMlStore::startPendingFrame() {
    if (activeFrame_ || pendingFrames_.empty())
        return;
    const std::size_t pending = pendingFrames_.front();
    pendingFrames_.pop_front();
    if (residentFrame_ && *residentFrame_ == pending) {
        emit frameReady(pending);
        startPendingFrame();
        return;
    }
    startFrame(pending);
}

std::optional<std::array<double, 6>>
ExperimentalShieldingMlStore::tensor(std::size_t frame, std::size_t atom) const {
    if (!residentFrame_ || *residentFrame_ != frame
        || atom >= protein_->atomCount()) {
        return std::nullopt;
    }
    const std::size_t offset = atom * kOutputColumns;
    std::array<double, 6> result{};
    for (int i = 0; i < kOutputColumns; ++i) {
        const double value = residentValues_[offset + static_cast<std::size_t>(i)];
        if (!std::isfinite(value))
            return std::nullopt;
        result[static_cast<std::size_t>(i)] = value;
    }
    return result;
}

std::optional<double> ExperimentalShieldingMlStore::sample(
    std::size_t frame,
    std::size_t atom,
    ExperimentalShieldingMlScalar scalar) const {
    const auto values = tensor(frame, atom);
    if (!values)
        return std::nullopt;
    if (scalar == ExperimentalShieldingMlScalar::Isotropic)
        return (*values)[0];
    if (scalar == ExperimentalShieldingMlScalar::T2Component)
        return (*values)[1];
    double squared = 0.0;
    for (std::size_t i = 1; i < values->size(); ++i)
        squared += (*values)[i] * (*values)[i];
    return std::sqrt(squared);
}

}  // namespace h5reader::model
