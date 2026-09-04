#pragma once

#include "../model/ScientificAlignment.h"

#include <QPointer>
#include <QString>
#include <QStringList>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace h5reader::model {
class Conformation;
class QtConformationSnapshot;
class QtProtein;
}  // namespace h5reader::model

namespace h5reader::io {

class ModelInputExporter final {
public:
    struct Result {
        bool ok = false;
        bool cancelled = false;
        QString error;
        std::size_t frameCount = 0;
        std::size_t atomCount = 0;
        std::size_t bondCount = 0;
        bool trajectory = false;
        model::ScientificAlignmentResult alignment;
    };

    ModelInputExporter(const model::QtProtein* protein, model::Conformation* conformation, QString outputDirectory);
    ~ModelInputExporter();

    bool prepare(QString* error);
    std::optional<Result> advance();
    void requestCancel() { cancelRequested_ = true; }

private:
    enum class Phase {
        Positions,
        Alignment,
        Features,
        Files,
        Finished,
    };

    bool collectPositionFrame(std::size_t frame, QString* error);
    bool buildAlignment(QString* error);
    bool projectFeatureFrame(std::size_t frame, QString* error);
    bool writeNextFile(QString* error);

    Result finishSuccess();
    Result finishFailure(const QString& error, bool cancelled = false);
    void removeWrittenFiles();

    const model::QtProtein* protein_ = nullptr;
    QPointer<model::Conformation> conformation_;
    QString outputDirectory_;

    Phase phase_ = Phase::Positions;
    bool prepared_ = false;
    bool cancelRequested_ = false;
    bool trajectory_ = false;
    std::size_t frameCount_ = 0;
    std::size_t atomCount_ = 0;
    std::size_t bondCount_ = 0;
    std::size_t nextFrame_ = 0;
    std::size_t nextFile_ = 0;
    float totalFormalCharge_ = 0.0f;

    model::ScientificPositionTable sourcePositions_;
    model::ScientificAlignmentResult alignment_;
    std::vector<std::size_t> fitAtomIndices_;
    std::vector<std::size_t> residueForAtom_;
    std::vector<std::uint64_t> originalFrameIndices_;
    std::vector<double> frameTimesPicoseconds_;

    std::vector<float> positions_;
    std::vector<float> scalars_;
    std::vector<std::uint8_t> scalarValid_;
    std::vector<float> applicability_;
    std::vector<float> l1_;
    std::vector<double> rotations_;
    std::vector<double> translations_;
    std::vector<std::uint8_t> fitStatus_;
    std::vector<double> fitRmsdAngstrom_;
    std::vector<int> atomRoles_;
    std::vector<int> hybridisations_;
    QStringList writtenFiles_;
};

}  // namespace h5reader::io
