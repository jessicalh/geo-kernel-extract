#pragma once

#include <QDockWidget>
#include <QPointer>

class QCheckBox;
class QLabel;
class QPushButton;
class QSpinBox;
class QTabBar;
class QToolButton;

namespace h5reader::model {
class AtomSelection;
class Conformation;
class DashboardPanelModel;
class DashboardSignalModel;
class DftShieldingStore;
class QtProtein;
class SignalBinding;
class TrajectorySignalCatalog;
}

namespace h5reader::app {

class DashboardDisplayController;
struct DashboardSmokeSummary;
class StripStackWidget;
class TimeViewportController;

class DashboardStripDock final : public QDockWidget {
    Q_OBJECT

public:
    explicit DashboardStripDock(QWidget* parent = nullptr);
    ~DashboardStripDock() override = default;

    void setContext(const model::QtProtein* protein, model::Conformation* conformation);
    void setSignalModels(model::TrajectorySignalCatalog* catalog, model::DashboardSignalModel* activeModel);
    void setPanelModel(model::DashboardPanelModel* panelModel);
    void setSelection(model::AtomSelection* selection);
    void setDftStore(model::DftShieldingStore* store);
    void setTimeViewport(TimeViewportController* viewport);
    DashboardSmokeSummary smokeSummary() const;
    int stripDisplaySinkCount() const;
    int spectrumDisplaySinkCount() const;

signals:
    void revealRequested(const model::SignalBinding& binding);
    void metricPickerRequested();

public slots:
    void setFrame(int frame);

private slots:
    void refreshTracks();
    void updateViewportReadout(int first, int last);
    void syncPanelTabs();
    void onPanelTabChanged(int row);
    void onPanelTabCloseRequested(int row);
    void onAddPanelRequested();

private:
    DashboardDisplayController* controller_ = nullptr;
    QPointer<StripStackWidget> stackWidget_;
    QPointer<QTabBar> panelTabs_;
    QPointer<QToolButton> addPanelButton_;
    QPointer<QCheckBox> followBox_;
    QPointer<QPushButton> metricButton_;
    QPointer<QSpinBox> windowFramesSpin_;
    QPointer<QLabel> viewportReadout_;
    QPointer<QLabel> statusLabel_;
    QPointer<TimeViewportController> timeViewport_;
    QPointer<model::DashboardPanelModel> panelModel_;
    QPointer<model::AtomSelection> selection_;
    int frame_ = 0;
};

}  // namespace h5reader::app
