#pragma once

#include <QWidget>
#include <QJsonObject>

#include <cstddef>
#include <vector>

class QLabel;
class QAbstractButton;

namespace h5reader::model {
class AtomSelection;
class QtProtein;
class TrajectorySignalCatalog;
class TransformedConformation;
}

namespace h5reader::app {
class MoleculeScene;
class QtPlaybackController;

// Owned by ReaderMainWindow and destroyed before the loaded run it describes.
class SelectionContextWidget final : public QWidget {
    Q_OBJECT
public:
    SelectionContextWidget(const model::QtProtein& protein,
                           model::AtomSelection& selection,
                           model::TransformedConformation& conformation,
                           MoleculeScene& scene, QtPlaybackController& playback,
                           const model::TrajectorySignalCatalog& catalog,
                           QWidget* parent);

    void refresh();
    void setInteractionEnabled(bool enabled);
    bool triggerAction(const QString& action);
    QJsonObject measurementState() const;
    QJsonObject stateJson() const;

private:
    QString atomLabel(std::size_t atom) const;
    QString atomList(const std::vector<std::size_t>& atoms) const;

    const model::QtProtein& protein_;
    model::AtomSelection& selection_;
    model::TransformedConformation& conformation_;
    MoleculeScene& scene_;
    QtPlaybackController& playback_;
    const model::TrajectorySignalCatalog& catalog_;
    bool interactionEnabled_ = true;

    QLabel* atoms_ = nullptr;
    QLabel* kind_ = nullptr;
    QLabel* value_ = nullptr;
    QLabel* meaning_ = nullptr;
    QLabel* highlight_ = nullptr;
    QLabel* camera_ = nullptr;
    QAbstractButton* clearSelection_ = nullptr;
    QAbstractButton* clearHighlight_ = nullptr;
    QAbstractButton* releaseCamera_ = nullptr;
};
}  // namespace h5reader::app
