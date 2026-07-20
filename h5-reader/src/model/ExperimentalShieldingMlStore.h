#pragma once

#include "../io/QtFieldCatalog.gen.h"

#include <QHash>
#include <QObject>
#include <QPointer>
#include <QProcess>
#include <QString>
#include <QStringList>
#include <QTemporaryDir>
#include <QVector>

#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <optional>
#include <vector>

namespace h5reader::model {

class Conformation;
class QtConformationSnapshot;
class QtProtein;

enum class ExperimentalShieldingMlScalar : std::uint8_t {
    Isotropic,
    T2Magnitude,
    T2Component,
};

class ExperimentalShieldingMlStore final : public QObject {
    Q_OBJECT

public:
    ExperimentalShieldingMlStore(const QtProtein* protein,
                                 Conformation* conformation,
                                 QString modelPath,
                                 QString runtimeManifestPath,
                                 QString extractionManifestPath,
                                 QString helperPath,
                                 QString device,
                                 QString fallbackHelperPath,
                                 QObject* parent = nullptr);
    ~ExperimentalShieldingMlStore() override;

    static bool ManifestHasInferenceSchema(const QString& manifestPath,
                                           QString* reason = nullptr);

    bool isReady() const { return ready_; }
    QString errorReason() const { return errorReason_; }
    QString modelId() const { return contract_.modelId; }
    QString device() const { return device_; }
    bool usingFallback() const { return fallbackAttempted_; }
    bool isRunning() const;

    void requestFrame(std::size_t frame);
    std::optional<double> sample(std::size_t frame,
                                 std::size_t atom,
                                 ExperimentalShieldingMlScalar scalar) const;
    std::optional<std::array<double, 6>> tensor(std::size_t frame,
                                                std::size_t atom) const;

signals:
    void frameReady(std::size_t frame);
    void runtimeChanged();

private:
    enum class FeatureAxis : std::uint8_t {
        Atom,
        Residue,
    };

    enum class FeatureKind : std::uint8_t {
        Scalar,
        L1,
        L2,
    };

    enum class FeatureLayout : std::uint8_t {
        Scalar,
        ScalarColumns,
        Vector,
        T2,
        Full9T0,
        Full9T2,
        RingT0,
        RingT1,
        RingT2,
    };

    enum class FeatureValidity : std::uint8_t {
        Required,
        ExternalMask,
    };

    enum class FeatureScale : std::uint8_t {
        None,
        ManifestBsIntensity,
    };

    struct FeatureSource {
        QString fileName;
        QString stem;
        io::FieldKind field = io::FieldKind::Count;
    };

    struct FeatureSpec {
        QString key;
        FeatureAxis axis = FeatureAxis::Atom;
        FeatureKind kind = FeatureKind::Scalar;
        FeatureLayout layout = FeatureLayout::Scalar;
        FeatureValidity validity = FeatureValidity::Required;
        FeatureScale scale = FeatureScale::None;
        bool center = false;
        bool emitMask = false;
        QVector<int> columns;
        QVector<int> maskColumns;
        QVector<int> emitMaskColumns;
        QVector<FeatureSource> sources;
        QString maskFileName;
        QString maskStem;
        std::optional<io::FieldKind> maskField;
    };

    struct Contract {
        QString modelId;
        double radius = 0.0;
        int maxNeighbors = 0;
        int radialDim = 0;
        int expectedL1Channels = 0;
        int expectedL2Channels = 0;
        int expectedScalarChannels = 0;
        int expectedApplicabilityChannels = 0;
        QStringList expectedL1Names;
        QStringList expectedL2Names;
        QStringList expectedScalarNames;
        QStringList expectedApplicabilityNames;
        QStringList labelKeys;
        QHash<QString, QHash<QString, std::int64_t>> labelVocabs;
        QVector<FeatureSpec> features;
        QStringList ringTypeOrder;
        QVector<int> ringActive;
        QVector<double> ringIntensity;
    };

    bool loadContract(const QString& manifestPath);
    bool validateExtractionManifest(const QString& extractionManifestPath);
    bool validateInitialFrameInputs();
    bool buildInput(std::size_t frame,
                    const QtConformationSnapshot& snapshot,
                    const QString& inputDir,
                    QString& error) const;
    void startFrame(std::size_t frame);
    void launchActiveFrame();
    bool scheduleCpuFallback(const QString& reason);
    void finishProcess(int exitCode, QProcess::ExitStatus exitStatus);
    void failActiveFrame(const QString& reason);
    void startPendingFrame();

    const QtProtein* protein_ = nullptr;
    QPointer<Conformation> conformation_;
    QString modelPath_;
    QString runtimeManifestPath_;
    QString extractionManifestPath_;
    QString helperPath_;
    QString device_;
    QString fallbackHelperPath_;
    Contract contract_;
    bool ready_ = false;
    QString errorReason_;

    QProcess* process_ = nullptr;
    QTemporaryDir workRoot_;
    std::optional<std::size_t> activeFrame_;
    std::deque<std::size_t> pendingFrames_;
    QString activeDir_;
    QString activeOutput_;

    std::optional<std::size_t> residentFrame_;
    std::vector<float> residentValues_;
    std::vector<std::size_t> failedFrames_;
    bool fallbackAttempted_ = false;
};

}  // namespace h5reader::model
