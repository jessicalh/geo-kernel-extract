#include "ExperimentalShieldingMlStore.h"

#include "Conformation.h"
#include "QtConformationSnapshot.h"
#include "QtProtein.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ErrorBus.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLoggingCategory>
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
constexpr double kSmoothFiniteNormalisation = 1.14136;
constexpr double kPi = 3.14159265358979323846;

bool readManifest(const QString& path, QJsonObject& out, QString& error) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        error = QStringLiteral("could not open manifest: %1").arg(path);
        return false;
    }
    QJsonParseError parseError{};
    const QJsonDocument document = QJsonDocument::fromJson(file.readAll(), &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject()) {
        error = QStringLiteral("invalid manifest JSON: %1").arg(parseError.errorString());
        return false;
    }
    out = document.object();
    return true;
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

double finiteOrZero(double value) {
    return std::isfinite(value) ? value : 0.0;
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

}  // namespace

ExperimentalShieldingMlStore::ExperimentalShieldingMlStore(
    const QtProtein* protein,
    Conformation* conformation,
    QString modelPath,
    QString manifestPath,
    QString helperPath,
    QString device,
    QString fallbackHelperPath,
    QString modelId,
    QObject* parent)
    : QObject(parent),
      protein_(protein),
      conformation_(conformation),
      modelPath_(std::move(modelPath)),
      manifestPath_(std::move(manifestPath)),
      helperPath_(std::move(helperPath)),
      device_(std::move(device)),
      fallbackHelperPath_(std::move(fallbackHelperPath)),
      modelId_(std::move(modelId)),
      process_(new QProcess(this)) {
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
        ready_ = loadContract(manifestPath_, modelId_);
    }

    ACONNECT(process_,
             &QProcess::finished,
             this,
             &ExperimentalShieldingMlStore::finishProcess);
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
                                      manifestPath_);
    }
}

ExperimentalShieldingMlStore::~ExperimentalShieldingMlStore() {
    if (process_->state() != QProcess::NotRunning) {
        process_->kill();
        process_->waitForFinished();
    }
}

bool ExperimentalShieldingMlStore::ManifestHasInferenceSchema(const QString& manifestPath,
                                                                QString* reason) {
    QJsonObject manifest;
    QString error;
    if (!readManifest(manifestPath, manifest, error)) {
        if (reason)
            *reason = error;
        return false;
    }
    const QJsonObject schema = manifest.value(QStringLiteral("inference_schema")).toObject();
    const QJsonObject models = schema.value(QStringLiteral("models")).toObject();
    const bool valid = schema.value(QStringLiteral("version")).toInt() == 1
                       && models.value(QStringLiteral("full")).isObject()
                       && models.value(QStringLiteral("no_mopac_no_tripeptide")).isObject()
                       && schema.value(QStringLiteral("label_vocabs")).isObject();
    if (!valid && reason)
        *reason = QStringLiteral("manifest inference_schema version 1 is missing");
    return valid;
}

bool ExperimentalShieldingMlStore::isRunning() const {
    return process_->state() != QProcess::NotRunning;
}

bool ExperimentalShieldingMlStore::loadContract(const QString& manifestPath,
                                                 const QString& modelId) {
    QJsonObject manifest;
    if (!readManifest(manifestPath, manifest, errorReason_))
        return false;

    const QJsonObject schema = manifest.value(QStringLiteral("inference_schema")).toObject();
    if (schema.value(QStringLiteral("version")).toInt() != 1) {
        errorReason_ = QStringLiteral("unsupported or missing inference_schema version");
        return false;
    }
    const QJsonObject graph = schema.value(QStringLiteral("graph")).toObject();
    contract_.radius = graph.value(QStringLiteral("radius_angstrom")).toDouble();
    contract_.maxNeighbors = graph.value(QStringLiteral("max_neighbors")).toInt();
    contract_.radialDim = graph.value(QStringLiteral("radial_dim")).toInt();
    if (!(contract_.radius > 0.0) || contract_.maxNeighbors <= 0 || contract_.radialDim <= 1) {
        errorReason_ = QStringLiteral("invalid graph parameters in inference_schema");
        return false;
    }

    for (const QJsonValue value : schema.value(QStringLiteral("label_keys")).toArray()) {
        if (value.isString())
            contract_.labelKeys.push_back(value.toString());
    }
    const QJsonObject vocabularies = schema.value(QStringLiteral("label_vocabs")).toObject();
    for (const QString& key : contract_.labelKeys) {
        QHash<QString, std::int64_t> vocabulary;
        const QJsonObject source = vocabularies.value(key).toObject();
        for (auto it = source.constBegin(); it != source.constEnd(); ++it) {
            if (it.value().isDouble())
                vocabulary.insert(it.key(), static_cast<std::int64_t>(it.value().toDouble()));
        }
        if (!vocabulary.contains(QStringLiteral("<UNK>"))) {
            errorReason_ = QStringLiteral("label vocabulary %1 has no <UNK>").arg(key);
            return false;
        }
        contract_.labelVocabs.insert(key, vocabulary);
    }
    if (contract_.labelKeys.size() != 10) {
        errorReason_ = QStringLiteral("inference schema requires exactly 10 label columns");
        return false;
    }

    const QJsonObject model = schema.value(QStringLiteral("models")).toObject()
                                  .value(modelId)
                                  .toObject();
    if (model.isEmpty()) {
        errorReason_ = QStringLiteral("inference schema has no model %1").arg(modelId);
        return false;
    }
    const auto parseSpecs = [this](const QJsonArray& values,
                                   QVector<FeatureSpec>& destination,
                                   QString& error) {
        for (const QJsonValue value : values) {
            const QJsonObject item = value.toObject();
            FeatureSpec spec;
            spec.stem = item.value(QStringLiteral("stem")).toString();
            spec.channels = item.value(QStringLiteral("channels")).toInt();
            const std::optional<io::FieldKind> field =
                io::FindFieldByStem(spec.stem.toStdString());
            if (spec.stem.isEmpty() || spec.channels <= 0 || !field) {
                error = QStringLiteral("invalid or unknown feature specification: %1").arg(spec.stem);
                return false;
            }
            spec.field = *field;
            destination.push_back(spec);
        }
        return true;
    };
    if (!parseSpecs(model.value(QStringLiteral("l1")).toArray(), contract_.l1, errorReason_)
        || !parseSpecs(model.value(QStringLiteral("l2")).toArray(), contract_.l2, errorReason_)
        || !parseSpecs(model.value(QStringLiteral("scalars")).toArray(), contract_.scalars, errorReason_)) {
        return false;
    }
    const QJsonObject widths = model.value(QStringLiteral("expected_channels")).toObject();
    contract_.expectedL1Channels = widths.value(QStringLiteral("l1")).toInt();
    contract_.expectedL2Channels = widths.value(QStringLiteral("l2")).toInt();
    contract_.expectedScalarChannels =
        widths.value(QStringLiteral("scalars_with_topology")).toInt();
    const auto channels = [](const QVector<FeatureSpec>& specs) {
        return std::accumulate(specs.cbegin(), specs.cend(), 0,
                               [](int sum, const FeatureSpec& spec) {
                                   return sum + spec.channels;
                               });
    };
    constexpr int kTopologyScalars = 5;
    if (contract_.expectedL1Channels <= 0 || contract_.expectedL2Channels <= 0
        || contract_.expectedScalarChannels <= kTopologyScalars
        || channels(contract_.l1) != contract_.expectedL1Channels
        || channels(contract_.l2) != contract_.expectedL2Channels
        || channels(contract_.scalars) + kTopologyScalars != contract_.expectedScalarChannels) {
        errorReason_ = QStringLiteral("inference schema feature widths do not match its feature lists");
        return false;
    }
    return true;
}

bool ExperimentalShieldingMlStore::buildInput(
    std::size_t frame,
    const QtConformationSnapshot& snapshot,
    const QString& inputDir,
    QString& error) const {
    const std::size_t atomCount = protein_->atomCount();
    if (atomCount == 0) {
        error = QStringLiteral("protein has no atoms");
        return false;
    }
    QDir dir;
    if (!dir.mkpath(inputDir)) {
        error = QStringLiteral("could not create inference input directory");
        return false;
    }

    std::vector<float> positions(atomCount * 3);
    Vec3 centroid = Vec3::Zero();
    for (std::size_t atom = 0; atom < atomCount; ++atom)
        centroid += conformation_->atomPosition(frame, atom);
    centroid /= static_cast<double>(atomCount);
    for (std::size_t atom = 0; atom < atomCount; ++atom) {
        const Vec3 position = conformation_->atomPosition(frame, atom) - centroid;
        for (int axis = 0; axis < 3; ++axis)
            positions[(atom * 3) + static_cast<std::size_t>(axis)] =
                static_cast<float>(position[axis]);
    }

    const auto featureChannels = [](const QVector<FeatureSpec>& specs) {
        return std::accumulate(specs.cbegin(), specs.cend(), 0,
                               [](int sum, const FeatureSpec& spec) { return sum + spec.channels; });
    };
    const int l1Channels = featureChannels(contract_.l1);
    const int l2Channels = featureChannels(contract_.l2);
    const int scalarFeatureChannels = featureChannels(contract_.scalars);
    constexpr int kTopologyScalars = 5;
    const int scalarChannels = scalarFeatureChannels + kTopologyScalars;
    std::vector<float> l1(atomCount * static_cast<std::size_t>(l1Channels) * 3, 0.0F);
    std::vector<float> l2(atomCount * static_cast<std::size_t>(l2Channels) * 5, 0.0F);
    std::vector<float> scalars(atomCount * static_cast<std::size_t>(scalarChannels), 0.0F);

    const auto copyFeatures = [&](const QVector<FeatureSpec>& specs,
                                  int width,
                                  std::vector<float>& destination,
                                  int destinationChannels,
                                  bool t2) {
        int channelOffset = 0;
        for (const FeatureSpec& spec : specs) {
            if (snapshot.has(spec.field)) {
                const NpyColumn& column = snapshot.column(spec.field);
                const bool packedT2 = t2 && spec.channels == 1 && column.cols == 9;
                const int expectedCols = spec.channels * width;
                const int sourceCols = packedT2 ? 9 : expectedCols;
                const std::size_t expectedValues =
                    atomCount * static_cast<std::size_t>(sourceCols);
                if (column.rows != static_cast<int>(atomCount)
                    || column.cols != sourceCols || column.data.size() != expectedValues) {
                    error = QStringLiteral("feature %1 shape %2x%3 does not match %4 atoms x %5 columns")
                                .arg(spec.stem)
                                .arg(column.rows)
                                .arg(column.cols)
                                .arg(atomCount)
                                .arg(expectedCols);
                    return false;
                }
                for (std::size_t atom = 0; atom < atomCount; ++atom) {
                    const double* row = column.row(atom);
                    for (int channel = 0; channel < spec.channels; ++channel) {
                        for (int component = 0; component < width; ++component) {
                            const int source = packedT2
                                                   ? 4 + component
                                                   : (channel * width) + component;
                            const std::size_t target =
                                ((atom * static_cast<std::size_t>(destinationChannels))
                                 + static_cast<std::size_t>(channelOffset + channel))
                                    * static_cast<std::size_t>(width)
                                + static_cast<std::size_t>(component);
                            destination[target] = static_cast<float>(finiteOrZero(row[source]));
                        }
                    }
                }
            }
            channelOffset += spec.channels;
        }
        return true;
    };
    if (!copyFeatures(contract_.l1, 3, l1, l1Channels, false)
        || !copyFeatures(contract_.l2, 5, l2, l2Channels, true)
        || !copyFeatures(contract_.scalars, 1, scalars, scalarChannels, false)) {
        return false;
    }

    std::vector<std::int64_t> labelIds(atomCount * static_cast<std::size_t>(contract_.labelKeys.size()));
    for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
        const QtAtom& atom = protein_->atom(atomIndex);
        const bool validResidue = atom.residueIndex >= 0
                                  && static_cast<std::size_t>(atom.residueIndex) < protein_->residueCount();
        const QtResidue* residue = validResidue
                                       ? &protein_->residue(static_cast<std::size_t>(atom.residueIndex))
                                       : nullptr;
        const QtResidueNames* residueNames = validResidue
                                                ? &protein_->residueNames(static_cast<std::size_t>(atom.residueIndex))
                                                : nullptr;

        const std::size_t scalarBase = (atomIndex * static_cast<std::size_t>(scalarChannels))
                                       + static_cast<std::size_t>(scalarFeatureChannels);
        scalars[scalarBase] = static_cast<float>(atom.formalCharge);
        scalars[scalarBase + 1] = residue ? static_cast<float>(residue->protonationVariantIndex) : -1.0F;
        scalars[scalarBase + 2] = static_cast<float>(atom.polarH);
        scalars[scalarBase + 3] = atom.aromatic ? 1.0F : 0.0F;
        scalars[scalarBase + 4] = atom.isExchangeable ? 1.0F : 0.0F;

        for (int labelIndex = 0; labelIndex < contract_.labelKeys.size(); ++labelIndex) {
            const QString& key = contract_.labelKeys[labelIndex];
            QString value;
            if (key == QStringLiteral("element"))
                value = QString::number(atom.AtomicNumber());
            else if (key == QStringLiteral("iupac"))
                value = protein_->atomNames(atomIndex).iupac;
            else if (key == QStringLiteral("residue"))
                value = residueNames ? residueNames->iupac : QStringLiteral("NA");
            else if (key == QStringLiteral("amber_residue"))
                value = residueNames ? residueNames->amber : QStringLiteral("NA");
            else if (key == QStringLiteral("terminal"))
                value = residue ? enumNumber(residue->terminalState) : QStringLiteral("NA");
            else if (key == QStringLiteral("variant"))
                value = residue ? QString::number(residue->protonationVariantIndex) : QStringLiteral("-1");
            else if (key == QStringLiteral("formal_charge"))
                value = QString::number(atom.formalCharge);
            else if (key == QStringLiteral("polar_h"))
                value = enumNumber(atom.polarH);
            else if (key == QStringLiteral("backbone_role"))
                value = enumNumber(atom.backboneRole);
            else if (key == QStringLiteral("planar_group"))
                value = enumNumber(atom.planarGroup);

            const QHash<QString, std::int64_t>& vocabulary = contract_.labelVocabs[key];
            labelIds[(atomIndex * static_cast<std::size_t>(contract_.labelKeys.size()))
                     + static_cast<std::size_t>(labelIndex)] =
                vocabulary.value(value, vocabulary.value(QStringLiteral("<UNK>")));
        }
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
        const Vec3 dst(positions[destinationOffset],
                       positions[destinationOffset + 1],
                       positions[destinationOffset + 2]);
        for (std::size_t source = 0; source < atomCount; ++source) {
            if (source == destination)
                continue;
            const std::size_t sourceOffset = source * 3;
            const Vec3 src(positions[sourceOffset],
                           positions[sourceOffset + 1],
                           positions[sourceOffset + 2]);
            const double distance = (dst - src).norm();
            if (distance < contract_.radius)
                candidates.push_back({source, distance});
        }
        if (candidates.size() > static_cast<std::size_t>(contract_.maxNeighbors)) {
            std::partial_sort(candidates.begin(),
                              candidates.begin() + contract_.maxNeighbors,
                              candidates.end(),
                              [](const Candidate& a, const Candidate& b) {
                                  if (a.distance != b.distance)
                                      return a.distance < b.distance;
                                  return a.source < b.source;
                              });
            candidates.resize(static_cast<std::size_t>(contract_.maxNeighbors));
        }
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
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

    const auto bytesFor = [&error](std::size_t count, std::size_t width, const QString& name) {
        if (count > static_cast<std::size_t>(std::numeric_limits<qsizetype>::max()) / width) {
            error = QStringLiteral("%1 input is too large").arg(name);
            return qsizetype{-1};
        }
        return static_cast<qsizetype>(count * width);
    };
    const qsizetype posBytes = bytesFor(positions.size(), sizeof(float), QStringLiteral("position"));
    const qsizetype l1Bytes = bytesFor(l1.size(), sizeof(float), QStringLiteral("l1"));
    const qsizetype l2Bytes = bytesFor(l2.size(), sizeof(float), QStringLiteral("l2"));
    const qsizetype scalarBytes = bytesFor(scalars.size(), sizeof(float), QStringLiteral("scalar"));
    const qsizetype labelBytes = bytesFor(labelIds.size(), sizeof(std::int64_t), QStringLiteral("label"));
    const qsizetype edgeBytes = bytesFor(edgeSrc.size(), sizeof(std::int64_t), QStringLiteral("edge"));
    const qsizetype radialBytes = bytesFor(radial.size(), sizeof(float), QStringLiteral("radial"));
    if (posBytes < 0 || l1Bytes < 0 || l2Bytes < 0 || scalarBytes < 0 || labelBytes < 0
        || edgeBytes < 0 || radialBytes < 0) {
        return false;
    }
    const QDir output(inputDir);
    if (!writeRaw(output.filePath(QStringLiteral("pos.bin")), positions.data(), posBytes, error)
        || !writeRaw(output.filePath(QStringLiteral("l1.bin")), l1.data(), l1Bytes, error)
        || !writeRaw(output.filePath(QStringLiteral("l2.bin")), l2.data(), l2Bytes, error)
        || !writeRaw(output.filePath(QStringLiteral("scalars.bin")), scalars.data(), scalarBytes, error)
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
           << "C1 " << l1Channels << '\n'
           << "C2 " << l2Channels << '\n'
           << "C0 " << scalarChannels << '\n'
           << "label_count " << contract_.labelKeys.size() << '\n'
           << "radial_dim " << contract_.radialDim << '\n';
    return true;
}

void ExperimentalShieldingMlStore::requestFrame(std::size_t frame) {
    ASSERT_THREAD(this);
    if (!ready_ || frame >= conformation_->frameCount())
        return;
    if (residentFrame_ && *residentFrame_ == frame) {
        emit frameReady(frame);
        return;
    }
    if (std::find(failedFrames_.cbegin(), failedFrames_.cend(), frame) != failedFrames_.cend()) {
        emit frameReady(frame);
        return;
    }
    if (activeFrame_) {
        pendingFrame_ = frame;
        return;
    }
    startFrame(frame);
}

void ExperimentalShieldingMlStore::startFrame(std::size_t frame) {
    conformation_->requestSnapshot(frame);
    const std::shared_ptr<const QtConformationSnapshot> snapshot = conformation_->snapshot(frame);
    if (!snapshot) {
        failedFrames_.push_back(frame);
        emit frameReady(frame);
        startPendingFrame();
        return;
    }

    activeFrame_ = frame;
    activeDir_ = QDir(workRoot_.path()).filePath(QStringLiteral("frame_%1").arg(frame));
    activeOutput_ = QDir(activeDir_).filePath(QStringLiteral("output.bin"));
    QString error;
    if (!buildInput(frame, *snapshot, QDir(activeDir_).filePath(QStringLiteral("input")), error)) {
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
    environment.insert(QStringLiteral("PATH"),
                       runtimeDir + QDir::listSeparator()
                           + environment.value(QStringLiteral("PATH")));
    process_->setProcessEnvironment(environment);
    process_->setWorkingDirectory(runtimeDir);
    process_->setProgram(helperPath_);
    process_->setArguments({modelPath_,
                            QDir(activeDir_).filePath(QStringLiteral("input")),
                            activeOutput_,
                            device_});
    qCInfo(cExperimentalMl).noquote()
        << QStringLiteral("event=inference_start frame=%1 model=%2 device=%3")
               .arg(*activeFrame_)
               .arg(modelId_)
               .arg(device_);
    process_->start();
}

bool ExperimentalShieldingMlStore::scheduleCpuFallback(const QString& reason) {
    if (!activeFrame_ || fallbackAttempted_ || fallbackHelperPath_.isEmpty()
        || device_ != QStringLiteral("rocm") || !QFileInfo(fallbackHelperPath_).isFile()) {
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
    QMetaObject::invokeMethod(this, [this]() { launchActiveFrame(); }, Qt::QueuedConnection);
    return true;
}

void ExperimentalShieldingMlStore::finishProcess(int exitCode, QProcess::ExitStatus exitStatus) {
    ASSERT_THREAD(this);
    if (!activeFrame_)
        return;
    const std::size_t frame = *activeFrame_;
    const QString stderrText = QString::fromUtf8(process_->readAllStandardError()).trimmed();
    if (!stderrText.isEmpty())
        qCInfo(cExperimentalMl).noquote() << stderrText;

    if (exitStatus != QProcess::NormalExit || exitCode != 0) {
        const QString reason = QStringLiteral("inference helper exited with code %1").arg(exitCode);
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
    residentFrame_ = frame;
    activeFrame_.reset();
    QDir(activeDir_).removeRecursively();
    activeDir_.clear();
    activeOutput_.clear();
    qCInfo(cExperimentalMl).noquote()
        << QStringLiteral("event=inference_ready frame=%1 atoms=%2 model=%3 device=%4")
               .arg(frame)
               .arg(protein_->atomCount())
               .arg(modelId_)
               .arg(device_);
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
    if (!pendingFrame_)
        return;
    const std::size_t pending = *pendingFrame_;
    pendingFrame_.reset();
    if (!residentFrame_ || *residentFrame_ != pending)
        startFrame(pending);
}

std::optional<std::array<double, 6>> ExperimentalShieldingMlStore::tensor(
    std::size_t frame,
    std::size_t atom) const {
    if (!residentFrame_ || *residentFrame_ != frame || atom >= protein_->atomCount())
        return std::nullopt;
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
