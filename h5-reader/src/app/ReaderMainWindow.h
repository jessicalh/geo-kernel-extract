// ReaderMainWindow — top-level QMainWindow for h5-reader.
//
// Holds the QVTKOpenGLNativeWidget as the central widget, owns the
// MoleculeScene and QtPlaybackController, wires frame-scrub controls
// in a toolbar and status bar. Designed to accept additional dock
// widgets (atom inspector, time-series tab) in later commits without
// restructuring the central layout.
//
// Shutdown protocol — see feedback_qt_discipline and the library
// viewer's MainWindow::shutdown(). When QApplication is about to
// quit, stop all timers, drop VTK references in order, and
// renderWindow_->Finalize() BEFORE Qt destroys the GL context.

#pragma once

#include <QMainWindow>
#include <QPointer>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkSmartPointer.h>

#include <memory>

class QLabel;
class QSlider;
class QSpinBox;
class QToolBar;
class QVTKOpenGLNativeWidget;

namespace h5reader::io {
struct QtLoadResult;
}

namespace h5reader::model {
class AtomSelection;
class DashboardPanelModel;
class DashboardSignalModel;
class DftShieldingStore;
class TrajectorySignalCatalog;
}

namespace h5reader::app {

class MoleculeScene;
class QtPlaybackController;
class DashboardDisplayController;
class TimeViewportController;

class ReaderMainWindow final : public QMainWindow {
    Q_OBJECT

public:
    // Takes the loader's result by rvalue — the window assumes ownership
    // of the protein and conformation (which owns the typed H5 reader).
    // The result must have ok=true; otherwise the caller should have
    // aborted already.
    explicit ReaderMainWindow(h5reader::io::QtLoadResult&& loaded, QWidget* parent = nullptr);
    ~ReaderMainWindow() override;

    // Start the embedded REST test surface bound to 127.0.0.1:<port>.
    // Port 0 asks the kernel to pick a free port. Returns the actually-bound
    // port, or 0 on failure. The bound port is also written to stderr as
    // `H5READER_REST_PORT=NNNNN\n` for the pytest fixture to scrape.
    quint16 startRestServer(quint16 port);

public slots:
    // Called from aboutToQuit. Stops timers, cancels any workers, and
    // finalises the VTK render window before Qt tears down the GL context.
    void shutdown();

protected:
    // Logs OpenGL vendor / renderer / version exactly once on the first
    // show. The GL context only exists after the widget has been mapped
    // and painted, so the log itself defers via QTimer::singleShot(0)
    // from inside this handler. Diagnostic-only — if Qt fell back to
    // ANGLE / software OpenGL, this is where it shows up.
    void showEvent(QShowEvent* event) override;

private slots:
    void onFrameChanged(int t);
    void onPlayPauseClicked();
    void onOpenDirectory();
    void onOpenSignalDisplays();
    void onPlaneLockToggled(bool checked);

private:
    void buildUi();
    void buildToolbar();
    void buildStatusBar();
    void updatePlaneLockAction();

    // The loaded model. Owned by the window for its lifetime.
    std::unique_ptr<h5reader::io::QtLoadResult> loaded_;

    // VTK viewport widget.
    QVTKOpenGLNativeWidget* vtkWidget_ = nullptr;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;

    // Scene + playback.
    MoleculeScene* scene_ = nullptr;
    QtPlaybackController* playback_ = nullptr;
    TimeViewportController* timeViewport_ = nullptr;

    // Atom picker + Atom Info dock. Picker is an event filter on the
    // VTK widget; Atom Info is tabified with the compact selection panel.
    class QtAtomPicker* picker_ = nullptr;
    class QtAtomInspectorDock* inspectorDock_ = nullptr;

    // Selection model — the QAbstractListModel for the ≤4-atom group — plus
    // its QListView panel. The picker feeds the model; the model fans focus
    // to Atom Info and the set to the measurement overlay/dashboard context.
    model::AtomSelection* selection_ = nullptr;
    class SelectionDock* selectionDock_ = nullptr;
    class SignalDisplayDialog* signalDisplayDialog_ = nullptr;

    // Scale-first dashboard signal state. AtomSelection supplies focus/context;
    // this model is the active signal/display set edited by SignalDisplayDialog.
    model::TrajectorySignalCatalog* signalCatalog_ = nullptr;
    model::DashboardSignalModel* dashboardSignals_ = nullptr;
    model::DashboardPanelModel* dashboardPanels_ = nullptr;

    // Unified strip dashboard. SignalDisplayDialog owns selection of active
    // signals/display modes; this dock renders strip-capable active signals.
    class DashboardStripDock* dashboardStripDock_ = nullptr;
    QPointer<DashboardDisplayController> dashboardController_;

    // DFT shielding provider for the strip chart's shielding panel — constructed
    // only when the run has a dft/ campaign (located by convention from the run
    // path). Window-owned (Qt parent); the dock holds a QPointer to it.
    model::DftShieldingStore* dftStore_ = nullptr;

    // Optional REST test surface — constructed by startRestServer(), only
    // when h5reader is launched with --rest <port>. Window-owned.
    class RestServer* restServer_ = nullptr;

    // Toolbar controls.
    QPointer<QSlider> frameSlider_;
    QPointer<QSpinBox> fpsSpinner_;
    QPointer<QAction> playAction_;
    QPointer<QAction> showRibbonAction_;
    QPointer<QAction> showRingsAction_;
    QPointer<QAction> showButterflyAction_;
    QPointer<QAction> showBFieldAction_;
    QPointer<QAction> signalDisplaysAction_;
    QPointer<QAction> planeLockAction_;

    // Status bar labels.
    QPointer<QLabel> proteinLabel_;
    QPointer<QLabel> frameLabel_;
    QPointer<QLabel> timeLabel_;

    bool shutdownDone_ = false;
    bool glInfoLogged_ = false;
};

}  // namespace h5reader::app
