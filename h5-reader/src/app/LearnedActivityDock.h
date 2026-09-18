#pragma once

#include "../model/Types.h"

#include <QDockWidget>
#include <QJsonArray>
#include <QJsonObject>
#include <QPointer>

#include <map>
#include <memory>
#include <vector>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;

namespace h5reader::model {
class Conformation;
class QtProtein;
class TransformedConformation;
}

namespace h5reader::app {
class MoleculeScene;
class QtPlaybackController;
class TensorGlyphActor;

// Optional, captured hidden 2e features. The model runs outside Reader; this
// dock binds the capture to the loaded atoms and draws it in the display frame.
class LearnedActivityDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit LearnedActivityDock(QWidget* parent = nullptr);
    ~LearnedActivityDock() override;
    void setContext(MoleculeScene* scene, const model::QtProtein* protein,
                    model::Conformation* raw, model::TransformedConformation* display,
                    QtPlaybackController* playback);
    bool load(const QString& path, QString* error);
    bool configure(const QJsonObject& options, QString* error);
    void clear();
    QJsonObject state() const;
    void openCapture();

private:
    struct Frame {
        std::vector<std::vector<model::Mat3>> channels;
    };
    void refresh();
    void updateControls();

    QPointer<MoleculeScene> scene_;
    const model::QtProtein* protein_ = nullptr;
    QPointer<model::Conformation> raw_;
    QPointer<model::TransformedConformation> display_;
    QPointer<QtPlaybackController> playback_;
    QMetaObject::Connection frameConnection_;
    QMetaObject::Connection transformConnection_;

    QString modelName_;
    QStringList channelNames_;
    std::vector<std::size_t> atoms_;
    std::map<std::size_t, Frame> frames_;
    std::vector<std::unique_ptr<TensorGlyphActor>> glyphs_;
    QJsonArray drawnSamples_;
    double referenceMagnitude_ = 1.0;

    QComboBox* channel_ = nullptr;
    QDoubleSpinBox* radius_ = nullptr;
    QDoubleSpinBox* opacity_ = nullptr;
    QCheckBox* visible_ = nullptr;
    QLabel* status_ = nullptr;
};
}  // namespace h5reader::app
