#pragma once

#include "../io/QtFieldCatalog.gen.h"

#include <QObject>
#include <QHash>
#include <QPointer>
#include <QProcess>
#include <QString>
#include <QTemporaryDir>
#include <QVector>

#include <array>
#include <cstddef>
#include <cstdint>
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
                                 QString manifestPath,
                                 QString helperPath,
                                 QString device,
                                 QString fallbackHelperPath,
                                 QString modelId,
                                 QObject* parent = nullptr);
    ~ExperimentalShieldingMlStore() override;

    static bool ManifestHasInferenceSchema(const QString& manifestPath,
                                           QString* reason = nullptr);

    bool isReady() const { return ready_; }
    QString errorReason() const { return errorReason_; }
    QString modelId() const { return modelId_; }
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
    struct FeatureSpec {
        QString stem;
        io::FieldKind field = io::FieldKind::Count;
        int channels = 0;
    };

    struct Contract {
        double radius = 0.0;
        int maxNeighbors = 0;
        int radialDim = 0;
        int expectedL1Channels = 0;
        int expectedL2Channels = 0;
        int expectedScalarChannels = 0;
        QVector<QString> labelKeys;
        QHash<QString, QHash<QString, std::int64_t>> labelVocabs;
        QVector<FeatureSpec> l1;
        QVector<FeatureSpec> l2;
        QVector<FeatureSpec> scalars;
    };

    bool loadContract(const QString& manifestPath, const QString& modelId);
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
    QString manifestPath_;
    QString helperPath_;
    QString device_;
    QString fallbackHelperPath_;
    QString modelId_;
    Contract contract_;
    bool ready_ = false;
    QString errorReason_;

    QProcess* process_ = nullptr;
    QTemporaryDir workRoot_;
    std::optional<std::size_t> activeFrame_;
    std::optional<std::size_t> pendingFrame_;
    QString activeDir_;
    QString activeOutput_;

    std::optional<std::size_t> residentFrame_;
    std::vector<float> residentValues_;
    std::vector<std::size_t> failedFrames_;
    bool fallbackAttempted_ = false;
};

}  // namespace h5reader::model
