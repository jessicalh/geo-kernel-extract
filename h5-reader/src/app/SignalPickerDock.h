// SignalPickerDock -- discovery surface for dashboard strip bindings.
//
// Current atom selection remains the local focus. This dock builds a nearby
// atom/residue candidate set from that focus and lists signal/display choices
// that could be pinned to the dashboard later.

#pragma once

#include <QDockWidget>
#include <QPointer>

#include <cstddef>
#include <optional>

class QCheckBox;
class QDoubleSpinBox;
class QLabel;
class QListWidget;
class QPushButton;
class QTableView;

namespace h5reader::model {
class AtomSelection;
class Conformation;
class DftShieldingStore;
class QtProtein;
}

namespace h5reader::app {

class NearbySignalModel;

class SignalPickerDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit SignalPickerDock(QWidget* parent = nullptr);
    ~SignalPickerDock() override = default;

    void setContext(const model::QtProtein* protein, model::Conformation* conformation);
    void setSelection(model::AtomSelection* selection);
    void setDftStore(model::DftShieldingStore* store);

public slots:
    void setFrame(int frame);

private slots:
    void onFocusChanged(std::size_t atomIdx);
    void onSelectionCleared();
    void onLiveToggled(bool live);
    void onRadiusChanged(double radius);
    void onCandidateChanged();
    void onAddClicked();

private:
    void refreshFocusLabel();
    void refreshCandidateAnchor();
    void refreshStripChoices();
    QString atomDisplayLabel(std::size_t atomIdx) const;

    NearbySignalModel* nearbyModel_ = nullptr;
    QPointer<model::AtomSelection> selection_;
    const model::QtProtein* protein_ = nullptr;
    QPointer<model::Conformation> conformation_;
    QPointer<model::DftShieldingStore> dftStore_;

    QPointer<QLabel> focusLabel_;
    QPointer<QCheckBox> liveBox_;
    QPointer<QDoubleSpinBox> radiusSpin_;
    QPointer<QTableView> candidatesView_;
    QPointer<QListWidget> stripList_;
    QPointer<QPushButton> addButton_;

    std::optional<std::size_t> latestFocusAtom_;
    int frame_ = 0;
};

}  // namespace h5reader::app
