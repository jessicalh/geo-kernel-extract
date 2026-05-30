// ReaderMainWindow — top-level QMainWindow for h5-reader.
//
// Holds the QVTKOpenGLNativeWidget as the central widget, owns the
// MoleculeScene and QtPlaybackController, wires frame-scrub controls
// in a toolbar and status bar. Designed to accept additional dock
// widgets (atom inspector, time-series tab) in later commits without
// restructuring the central layout.
//
// Shutdown protocol — see feedback_qt_discipline,
// spec/viewport_pipeline_2026-05-30.md §4.4, and the library viewer's
// MainWindow::shutdown(). When QApplication is about to quit, stop the
// REST server synchronously, stop all timers, then detach the render
// window from the widget. The explicit renderWindow_->Finalize() call
// is gone — setRenderWindow(nullptr) makes the GL context current and
// invokes Finalize through the QVTKRenderWindowAdapter's destructor in
// the right order (QVTKRenderWindowAdapter.cxx:150-166). Calling
// Finalize ourselves AFTER detach left the adapter holding a destroyed
// render window for the brief moment between the two calls.

#pragma once

#include <QMainWindow>
#include <QPointer>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkSmartPointer.h>

#include <memory>
#include <vector>

class QActionGroup;
class QDockWidget;
class QMenu;

class QLabel;
class QSlider;
class QSpinBox;
class QToolBar;
class QToolButton;
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
class TransformedConformation;
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

    // Hide or restore the docks (inspector, selection, dashboard strip)
    // wholesale. Hide preserves each dock's prior visibility so restore
    // returns each one to whatever it was — a dock that was already hidden
    // before hide() stays hidden after restore. Used by the viewport
    // harness (POST /docks/visible) to expand the central viewport so the
    // marker blob fits in more pixels and the drift detector finds it
    // reliably. No-op if the requested state matches the current state.
    void setDocksVisible(bool visible);

    // Read-only accessor for REST GET /docks/visible (if added) — true if
    // the docks are in the "visible" mode (some / all might be individually
    // hidden by the user, but we have not stashed and hidden them all).
    bool docksVisible() const { return docksHidden_ == false; }

    // Access to the wrapped TransformedConformation so REST handlers can
    // call setMode without re-walking the loader result. Non-null after
    // construction.
    h5reader::model::TransformedConformation* transformedConformation() const { return transformed_; }

public slots:
    // Called from aboutToQuit. Stops the REST server, stops timers, and
    // detaches the render window from the widget so the
    // QVTKRenderWindowAdapter's destructor calls Finalize in the right
    // GL context (per QVTKRenderWindowAdapter.cxx:150-166). The class
    // docstring above has the full reasoning.
    void shutdown();

protected:
    // Logs OpenGL vendor / renderer / version exactly once on the first
    // show. The GL context only exists after the widget has been mapped
    // and painted, so the log itself defers via QTimer::singleShot(0)
    // from inside this handler. Diagnostic-only — if Qt fell back to
    // ANGLE / software OpenGL, this is where it shows up.
    void showEvent(QShowEvent* event) override;

    // QSettings save runs here before the existing aboutToQuit → shutdown
    // chain. Tolerant: if save fails for any reason the user still gets
    // their window closed. event->accept() is unconditional.
    void closeEvent(QCloseEvent* event) override;

private slots:
    void onFrameChanged(int t);
    void onPlayPauseClicked();
    void onOpenDirectory();
    void onOpenSignalDisplays();
    void onPlaneLockTriggered();
    void onFocusCameraTriggered();
    void onNewmanProjectionTriggered();
    void onFreeCameraTriggered();
    void onInstrumentToggled(bool checked);

private:
    void buildUi();
    void buildToolbar();
    void buildStatusBar();
    // Update enabled state + checked state of the exclusive camera-mode
    // action group (Focus / Newman / Plane lock / Free). Gating: Focus
    // needs selection focus; Newman needs exactly 4 selected atoms; Plane
    // needs exactly 3; Free is always enabled. Checked state is sourced
    // from the composer's mode().kind so REST or programmatic changes
    // reflect in the toolbar too.
    void updateCameraModeActions();
    // Apply a TransformedConformation mode chosen from the Transform
    // popup menu. The mode int is the underlying enum value carried on
    // the QAction via setData().
    void applyTransformModeFromAction(QAction* action);

    // QSettings persistence — see kSettingsVersion in the .cpp for the
    // versioned QMainWindow state blob policy. Tolerant on restore (any
    // missing / mismatched key is silently skipped) so a fresh install
    // boots clean and an old install upgrades without losing usability.
    void saveAllSettings();
    void restoreAllSettings();

    // File ▸ Recent — prepend a path, dedupe, cap at 10, rebuild menu,
    // write to QSettings immediately. Called from the ctor with the
    // current runPath so the next session sees it at the top.
    void addToRecentFiles(const QString& path);
    void rebuildRecentFilesMenu(const QStringList& paths);
    void openRecentPath(const QString& path);

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

    // Camera input filter — Qt eventFilter on the VTK widget, intercepts
    // mouse + wheel before VTK's trackball. Routes gestures to the
    // CameraComposer (per spec/viewport_pipeline_2026-05-30.md §4).
    // Installed after the picker so Qt's filter chain runs THIS first;
    // double-clicks fall through to the picker.
    class CameraInputFilter* cameraInputFilter_ = nullptr;

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

    // Exclusive camera-mode action group. QActionGroup is the standard Qt
    // idiom for radio-style mutual exclusion across actions. Source of
    // truth is the composer's mode_; updateCameraModeActions() syncs the
    // checked state from composer->mode().kind whenever modeChanged fires.
    QActionGroup* cameraModeGroup_ = nullptr;
    QPointer<QAction> focusAction_;
    QPointer<QAction> newmanAction_;
    QPointer<QAction> planeLockAction_;
    QPointer<QAction> freeAction_;

    // Transform popup-menu button (Identity / Center COM / Fit backbone /
    // Fit selection). Drives TransformedConformation::setMode directly.
    QPointer<QToolButton> transformButton_;
    QPointer<QMenu> transformMenu_;
    QPointer<QActionGroup> transformGroup_;

    // Marker preset toggle for live demos. Same code path as
    // /selection/instrument REST endpoint.
    QPointer<QAction> instrumentAction_;

    // File ▸ Recent submenu — populated from QSettings on ctor restore.
    QPointer<QMenu> recentMenu_;

    // Playback toolbar — built by buildToolbar(); referenced from the
    // ctor after the docks exist to append their toggleViewAction()s.
    QPointer<QToolBar> playbackToolbar_;

    // Status bar labels.
    QPointer<QLabel> proteinLabel_;
    QPointer<QLabel> frameLabel_;
    QPointer<QLabel> timeLabel_;

    bool shutdownDone_ = false;
    bool glInfoLogged_ = false;

    // Wraps loaded_->conformation so consumers (scene, picker, overlays)
    // read positions through a runtime-switchable rigid-body transform.
    // Owned by the window. Built in the ctor immediately after the
    // loader returns; default mode is Identity so behaviour at startup
    // is identical to today.
    h5reader::model::TransformedConformation* transformed_ = nullptr;

    // Dock-hide state for setDocksVisible(). We stash each dock's
    // pre-hide visibility so restore puts a dock that was user-hidden
    // BEFORE setDocksVisible(false) back into the hidden state, not
    // a brittle "all visible" default. Empty when no hide is active.
    bool docksHidden_ = false;
    struct DockVis { QPointer<QDockWidget> dock; bool wasVisible; };
    std::vector<DockVis> stashedDockVisibility_;
};

}  // namespace h5reader::app
