// SignalDisplayDialog -- scale-first signal/display management surface.
//
// The dialog is intentionally renderer-agnostic. It lets the user add catalog
// descriptors to the active dashboard signal model and toggle display mode ids;
// it does not create strip widgets, tables, overlays, or VTK scene objects.

#pragma once

#include <QDialog>

#include <cstddef>
#include <memory>

namespace h5reader::model {
class AtomSelection;
class Conformation;
class DashboardPanelModel;
class DashboardSignalModel;
class QtProtein;
class TrajectorySignalCatalog;
}

namespace h5reader::app {

class SignalDisplayDialog final : public QDialog {
    Q_OBJECT

public:
    explicit SignalDisplayDialog(QWidget* parent = nullptr);
    ~SignalDisplayDialog() override;

    void setTrajectorySignalCatalog(model::TrajectorySignalCatalog* catalog);
    void setDashboardSignalModel(model::DashboardSignalModel* model);
    void setDashboardPanelModel(model::DashboardPanelModel* panelModel);
    void setContext(const model::QtProtein* protein, model::Conformation* conformation);
    void setSelection(model::AtomSelection* selection);

    model::TrajectorySignalCatalog* trajectorySignalCatalog() const;
    model::DashboardSignalModel* dashboardSignalModel() const;

public slots:
    void refreshCatalog();
    void setFrame(int frame);

private slots:
    void onFocusChanged(std::size_t atomIdx);
    void onSelectionCleared();
    void onLiveToggled(bool live);
    void onRadiusChanged(double radius);
    void onAnchorSelectionChanged();
    void onCandidateSelectionChanged();
    void onCandidateModeChanged();
    void onAddSelected();
    void refreshPanelTargets();
    void onActiveSelectionChanged();
    void onActiveModeToggled(bool checked);
    void onRemoveActive();

private:
    struct Impl;
    std::unique_ptr<Impl> d_;
};

}  // namespace h5reader::app
