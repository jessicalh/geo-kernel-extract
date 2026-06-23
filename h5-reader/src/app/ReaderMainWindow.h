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
#include <QString>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkSmartPointer.h>

#include "../model/CsaProbe.h"
#include "../model/VisualizationDefinition.h"

#include <memory>
#include <vector>

class QActionGroup;
class QDockWidget;
class QJsonObject;
class QJsonArray;
class QMenu;
class QString;

class QLabel;
class QSlider;
class QSpinBox;
class QToolBar;
class QToolButton;
class QWidget;
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
class TrajectoryFieldAvailability;
class TransformedConformation;
}

namespace h5reader::app {

class MoleculeScene;
class QtPlaybackController;
class DashboardDisplayController;
class DashboardSelectionController;
class TimeViewportController;

class ReaderMainWindow final : public QMainWindow {
    Q_OBJECT

public:
    // Constructs the chrome in a clean pre-load state. A calcset can be
    // loaded later through loadRunPath().
    explicit ReaderMainWindow(QWidget* parent = nullptr);
    ~ReaderMainWindow() override;

    // Compute one atom's CSA result (PAS shape + molecular frame) for the
    // current frame -- the same orchestration the focus-driven glyph uses.
    // Public so RestServer's GET /csa vets exactly what is drawn.
    model::AtomCsaResult probeAtomCsa(std::size_t atom);

    // Load or replace the current calcset in this window. The path is resolved
    // through QtProteinLoader::LoadRunPath. Returns false without changing the
    // current run when loading fails; lastLoadError() carries the loader error.
    bool loadRunPath(const QString& path);
    QString lastLoadError() const { return lastLoadError_; }

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

    // Toggle a named overlay by driving its toolbar action, so REST /
    // automation changes run the SAME path as a human click (per-frame
    // refresh for the kernel overlays) and keep the toolbar checkbox in
    // sync. name ∈ {ribbon, rings, butterfly, bfield, shadow} (+ aliases
    // fieldgrid/field/isosurface, streamlines/stream, occupancy/shells).
    // Returns false on an unknown name. Used by POST /overlay so the
    // headless snapshot harness can enable the field overlays, which
    // default off and are not persisted in QSettings.
    bool setOverlayVisible(const QString& name, bool on);

    // Set the butterfly (field-grid) isosurface |T0| contour threshold in ppm
    // and re-render. Runtime-tunable so the dominant-zone level can be swept
    // (POST /field/threshold) without a rebuild. Returns false if no field-grid
    // overlay is live (no scene loaded).
    bool setFieldThreshold(double ppm);

    // The two orthogonal butterfly knobs (POST /field/gaussian): the radial
    // Gaussian taper width σ in Å ("extent" — lobe reach) and the amplitude
    // gain ("peak"; 1.0 = true physics). Re-evaluate + re-render. Return false
    // if no field-grid overlay is live.
    bool setFieldExtent(double sigmaA);
    bool setFieldPeak(double amplitude);

    // Display-isolation: show only the given residues (their atoms) in the 3-D
    // view; an empty list restores the full structure. Maps residues → atoms and
    // drives MoleculeScene::setAtomFilter. Backs POST /filter and (soon) the
    // Filter toolbar toggle + nearby-residue checklist.
    void setResidueFilter(const std::vector<std::size_t>& residues);

    // Access to the wrapped TransformedConformation so REST handlers can
    // call setMode without re-walking the loader result. Null until a run
    // is loaded.
    h5reader::model::TransformedConformation* transformedConformation() const { return transformed_; }

    // Dashboard strip dock surface for the REST instrumentation harness.
    // Showing uses the same queued reveal path as in-product code, while
    // hiding is a plain dock visibility change.
    void setDashboardDockVisible(bool visible);
    bool dashboardDockVisible() const;
    int dashboardDockWidth() const;
    bool dashboardDockRaised() const;
    int dashboardOwnedPanelCount() const;
    int dashboardStripTrackCount() const;
    // Display manifest of the static panels (curve/spectrum/matrix/fixed-freq/
    // sequence-bar) -- the read-to-display hook that makes those panels, which
    // the strip-series path can't see, REST-visible for the metrics test.
    QJsonArray dashboardPanelManifest() const;
    bool openSignalDisplayPicker(QString* blockedReason = nullptr);
    QJsonObject signalDisplayPickerState() const;
    QJsonObject addSelectedSignalFromPicker();

    // Full operating-state + per-control {enabled,checked} snapshot. Backs
    // GET /ui/state and is the introspection used to verify the selectability
    // rules. Pure read of the live control states.
    QJsonObject uiStateJson() const;

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
    void onOpenFile();
    void onOpenDirectory();
    void onOpenSignalDisplays();
    void onFocusCameraTriggered();   // toggle: track focused atom / release to manual
    void onTransformFitClicked();

private:
    void buildUi();
    void buildToolbar();
    void buildStatusBar();
    void buildDocks();
    void installLoadedRun(h5reader::io::QtLoadResult&& loaded);
    void clearLoadedRun();
    // Recompute + redraw the focused atom's DFT shielding-tensor glyph (PAS
    // ellipsoid + molecular-frame axes) for the current frame; clears when
    // there is no focus or no DFT for the frame. Focus/frame driven.
    void updateCsaGlyph();
    void setEmptyState();
    // Single source of truth for control enable/checked state. Derives the
    // whole operating state (loaded / playable / playing / selection / data
    // capabilities) and sets every control accordingly. Called on every
    // state-changing signal — replaces the old scattered setLoadedControlsEnabled
    // + per-signal enable lambdas.
    void refreshControlStates();
    void applyOverlayActionState();
    void updateMutantAlternateAction(const QString& alternatePath);
    void syncRestServerContext();
    // Update enabled state + checked state of the exclusive camera-mode
    // action group (Focus / Newman / Plane lock / Free). Gating: Focus
    // needs selection focus; Newman needs exactly 4 selected atoms; Plane
    // needs exactly 3; Free is always enabled. Checked state is sourced
    // from the composer's mode().kind so REST or programmatic changes
    // reflect in the toolbar too.
    void updateCameraModeActions();
    void updateFitModeLabel();
    void revealDockQueued(QDockWidget* dock);
    void updateSelectionStatus();   // status-bar selection count + geometry kind
    void buildFilterMenu();         // (re)populate the Filter dropdown checklist
    void onFilterResidueToggled(std::size_t residue, bool on);
    void updateFilterButton();      // text/state of the Filter toolbar button
    void resetDashboardStateForRunLoad();
    // QSettings persistence — see kSettingsVersion in the .cpp for the
    // versioned QMainWindow state blob policy. Tolerant on restore (any
    // missing / mismatched key is silently skipped) so a fresh install
    // boots clean and an old install upgrades without losing usability.
    void saveAllSettings();
    void restoreAllSettings();

    // File ▸ Recent — prepend a path, dedupe, cap at 10, rebuild menu,
    // write to QSettings immediately. Called after a successful load so
    // the next session sees that calcset at the top.
    void addToRecentFiles(const QString& path);
    void rebuildRecentFilesMenu(const QStringList& paths);
    void openRecentPath(const QString& path);

    // The loaded model. Owned by the window for its lifetime.
    std::unique_ptr<h5reader::io::QtLoadResult> loaded_;

    // VTK viewport widget plus quiet empty-state placeholder.
    QPointer<QWidget> centralContainer_;
    QPointer<QLabel> emptyPlaceholder_;
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
    class NewmanDock* newmanDock_ = nullptr;

    // Camera input filter — Qt eventFilter on the VTK widget, intercepts
    // mouse + wheel before VTK's trackball. Routes gestures to the
    // CameraComposer (per spec/viewport_pipeline_2026-05-30.md §4).
    // Installed after the picker so Qt's filter chain runs THIS first;
    // double-clicks fall through to the picker.
    class CameraInputFilter* cameraInputFilter_ = nullptr;

    // Selection model — the QAbstractListModel for the ≤4-atom group. The picker
    // feeds the model; focus fans to Atom Info, the set to the measurement
    // overlay/dashboard context, and a count/kind summary to the status bar.
    // (The compact SelectionDock was retired — redundant with the in-scene
    // MeasurementOverlay.)
    model::AtomSelection* selection_ = nullptr;
    class SignalDisplayDialog* signalDisplayDialog_ = nullptr;

    // Scale-first dashboard signal state. AtomSelection supplies focus/context;
    // this model is the active signal/display set edited by SignalDisplayDialog.
    model::TrajectorySignalCatalog* signalCatalog_ = nullptr;
    model::DashboardSignalModel* dashboardSignals_ = nullptr;
    model::DashboardPanelModel* dashboardPanels_ = nullptr;
    std::shared_ptr<const model::TrajectoryFieldAvailability> fieldAvailability_;
    QPointer<DashboardSelectionController> dashboardSelectionController_;
    model::VisualizationContext visualizationContext_;

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
    QPointer<QAction> playBackAction_;      // continuous reverse
    QPointer<QAction> stepBackAction_;      // single frame back
    QPointer<QAction> stopAction_;          // stop / pause
    QPointer<QAction> stepForwardAction_;   // single frame forward
    QPointer<QAction> playForwardAction_;   // continuous forward
    QPointer<QAction> showRibbonAction_;
    QPointer<QAction> showRingsAction_;
    QPointer<QAction> showButterflyAction_;
    QPointer<QAction> showBFieldAction_;
    QPointer<QAction> showOccupancyAction_;
    QPointer<QAction> signalDisplaysAction_;

    // Display-isolation ("Filter"): a toolbar button whose dropdown is a live
    // checklist of residues near the focused atom (NearbySignalModel). The
    // checked set (filterResidues_) drives setResidueFilter.
    QPointer<QToolButton>    filterButton_;
    QPointer<QMenu>          filterMenu_;
    class NearbySignalModel* filterNearby_ = nullptr;
    std::vector<std::size_t> filterResidues_;

    // Focus is a self-contained checkable toggle (no action group): checked =
    // the composer tracks the focused atom, unchecked = manual camera.
    // updateCameraModeActions() syncs its checked/enabled state from
    // composer->mode().kind whenever modeChanged fires.
    QPointer<QAction> focusAction_;

    // One transform switch: text names the active stabilisation mode.
    QPointer<QAction> transformFitAction_;

    // File ▸ Recent submenu — populated from QSettings on ctor restore.
    QPointer<QMenu> fileMenu_;
    QPointer<QMenu> recentMenu_;
    QPointer<QAction> mutantAlternateAction_;

    // Toolbars — built by buildToolbar(); the ctor appends dock panel controls
    // once the docks and View -> Panels menu exist.
    QPointer<QToolBar> playbackToolbar_;
    QPointer<QToolBar> toolsToolbar_;

    // Status bar labels.
    QPointer<QLabel> selectionLabel_;
    QPointer<QLabel> proteinLabel_;
    QPointer<QLabel> frameLabel_;
    QPointer<QLabel> timeLabel_;

    bool shutdownDone_ = false;
    bool glInfoLogged_ = false;
    int lastDashboardSelectedCount_ = 0;
    QString lastLoadError_;

    // Wraps loaded_->conformation so consumers (scene, picker, overlays)
    // read positions through a runtime-switchable rigid-body transform.
    // Owned by the window and rebuilt on each successful load; startup
    // mode is FitSubset on the typed backbone.
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
