#include "ReaderMainWindow.h"

#include "CameraAnchorHelper.h"
#include "CameraComposer.h"
#include "CameraInputFilter.h"
#include "CameraMode.h"
#include "OrientationPolicy.h"
#include "MoleculeScene.h"
#include "QtAtomInspectorDock.h"
#include "QtAtomPicker.h"
#include "RestServer.h"
#include "QtBackboneRibbonOverlay.h"
#include "QtBFieldStreamOverlay.h"
#include "QtFieldGridOverlay.h"
#include "QtPlaybackController.h"
#include "TimeViewportController.h"
#include "MeasurementOverlay.h"
#include "QtRingPolygonOverlay.h"
#include "SelectionDock.h"
#include "DashboardStripDock.h"
#include "DashboardDisplayController.h"
#include "DashboardSelectionController.h"
#include "SignalDisplayDialog.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/DashboardLogging.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/StructuredLogger.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/QtProteinLoader.h"
#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/DashboardPanelModel.h"
#include "../model/DashboardSignalModel.h"
#include "../model/DftShieldingStore.h"
#include "../model/QtProtein.h"
#include "../model/TrajectoryConformation.h"
#include "../model/TrajectoryFieldAvailability.h"
#include "../model/TrajectorySignalCatalog.h"
#include "../model/TransformedConformation.h"
#include "../model/VisualizationRegistry.h"

#include <QDockWidget>

#include <QDir>
#include <QFileInfo>

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCloseEvent>
#include <QFileDialog>
#include <QFont>
#include <QJsonObject>
#include <QKeySequence>
#include <QLabel>
#include <QLoggingCategory>
#include <QMenu>
#include <QMenuBar>
#include <QProcess>
#include <QRegion>
#include <QSettings>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QStringList>
#include <QStatusBar>
#include <QStyle>
#include <QTimer>
#include <QToolBar>
#include <QToolButton>
#include <QUuid>
#include <QVariant>

#include <QVTKOpenGLNativeWidget.h>

#include <vtkRendererCollection.h>
#include <vtkCamera.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cWindow, "h5reader.window")

// QSettings — versioned state blob policy. Bump on dock-object
// additions or any layout-invalidating change so old blobs are
// silently discarded by QMainWindow::restoreState. Schema-evolution
// safe per ROBUSTNESS_BACKLOG_2026-05-30.md item 7.
constexpr int kSettingsVersion = 2;   // bumped: property docks now start hidden
constexpr int kMaxRecentFiles  = 10;

QString fitModeToolTip() {
    return QStringLiteral(
        "Stabilisation mode — click to switch.\n"
        "Locked backbone: Kabsch fit of the backbone (industry standard) — removes global tumbling; the backbone holds still while sidechains/loops move.\n"
        "Kabsch with give: all-atom fit — removes tumbling but lets real internal motion show.");
}

// Note: locateDftJobsDir was deleted as part of the 2026-05-31 SIMPLIFY
// pass; the DFT campaign now comes from the `.LGS` `dft.frames[]` array
// (see CalcsetManifest + DftShieldingStore).

}  // namespace

ReaderMainWindow::ReaderMainWindow(h5reader::io::QtLoadResult&& loaded,
                                   QWidget* parent)
    : QMainWindow(parent),
      loaded_(std::make_unique<h5reader::io::QtLoadResult>(std::move(loaded)))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("ReaderMainWindow"));

    qCInfo(cWindow).noquote() << "ctor entered";

    buildUi();
    buildToolbar();
    buildStatusBar();

    // Upstream data-transform layer (feedback_viewer_two_layers_transform_and_camera).
    // Wraps the loader's Conformation so consumers (scene, picker, overlays,
    // REST /positions) read positions through a runtime-switchable rigid-body
    // display transform. Startup mode is backbone fit with the iterative mean
    // seeded/anchored at frame 0 so the reader opens stationary.
    transformed_ = new h5reader::model::TransformedConformation(loaded_->conformation.get(), this);
    const auto backboneSubset =
        h5reader::model::TransformedConformation::BackboneSubset(*loaded_->protein);
    using TMode = h5reader::model::TransformedConformation::Mode;
    if (backboneSubset.size() >= 3) {
        transformed_->setMode(TMode::FitSubset, 0, backboneSubset);
    } else {
        qCWarning(cWindow).noquote()
            << "backbone fit unavailable at startup; falling back to all-atom fit";
        transformed_->setMode(TMode::FitReference, 0);
    }
    ACONNECT(transformed_, &h5reader::model::TransformedConformation::transformChanged,
             this, [this]() {
                 updateFitModeLabel();
                 if (scene_) scene_->refreshCurrentFrame();
             });
    updateFitModeLabel();

    // Scene binds to the VTK widget's render window. The scene reads
    // positions through the wrapped conformation so transform mode
    // changes are visible immediately. The widget is passed in so the
    // render scheduler (MoleculeScene::requestRender) can call VTK Render()
    // and let Qt blit the fresh FBO.
    scene_ = new MoleculeScene(vtkWidget_, renderWindow_, this);
    scene_->Build(*loaded_->protein, *transformed_);
    scene_->ResetCamera();
    // The plane-lock-specific signal is now a strict subset of the
    // composer's modeChanged below; updateCameraModeActions sources the
    // checked state from composer->mode() and gates the actions on the
    // current selection in one pass.
    ACONNECT(scene_, &MoleculeScene::cameraPlaneLockChanged,
             this, [this](bool) { updateCameraModeActions(); });
    if (scene_ && scene_->cameraComposer()) {
        ACONNECT(scene_->cameraComposer(), &CameraComposer::modeChanged,
                 this, [this]() { updateCameraModeActions(); });
    }

    // Playback controller — frameChanged drives the scene, which drives
    // the render. Toolbar controls drive the playback.
    const int T = static_cast<int>(loaded_->conformation->frameCount());
    playback_ = new QtPlaybackController(T, this);
    timeViewport_ = new TimeViewportController(T, this);

    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             scene_,    &MoleculeScene::setFrame);
    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             this,      &ReaderMainWindow::onFrameChanged);
    ACONNECT(playback_,     &QtPlaybackController::frameChanged,
             timeViewport_, &TimeViewportController::setCurrentFrame);
    ACONNECT(timeViewport_, &TimeViewportController::playbackFrameRequested,
             playback_,     &QtPlaybackController::setFrame);
    ACONNECT(playback_, &QtPlaybackController::playingChanged,
             this,      [this](bool playing) {
                 if (playAction_) {
                     playAction_->setIcon(style()->standardIcon(
                         playing ? QStyle::SP_MediaPause : QStyle::SP_MediaPlay));
                 }
             });

    // Slider ↔ playback. Slider emits while dragging; controller clamps
    // and re-emits frameChanged to everyone including us.
    if (frameSlider_) {
        frameSlider_->setRange(0, std::max(0, T - 1));
        ACONNECT(frameSlider_.data(), &QSlider::valueChanged,
                 playback_, &QtPlaybackController::setFrame);
    }
    if (fpsSpinner_) {
        fpsSpinner_->setRange(1, 60);
        fpsSpinner_->setValue(playback_->fps());
        ACONNECT(fpsSpinner_.data(), qOverload<int>(&QSpinBox::valueChanged),
                 playback_, &QtPlaybackController::setFps);
    }

    // Atom picker — event filter on the VTK widget. Emits
    // atomPicked(idx, modifiers) on double-click. It stays dumb: it does NOT
    // interpret the gesture. AtomSelection (below) is the sole consumer; it
    // applies the plain/Shift policy and fans typed changes to Atom Info,
    // the dashboard context, and the measurement overlay. Pull the renderer
    // smart-ptr out of the render window so the picker holds the same
    // instance MoleculeScene installed.
    auto* firstRenderer = renderWindow_->GetRenderers()->GetFirstRenderer();
    picker_ = new QtAtomPicker(vtkWidget_, firstRenderer,
                                loaded_->protein.get(),
                                transformed_,
                                playback_, this);

    // Camera input filter — installed AFTER the picker so Qt's filter
    // chain runs THIS first. Double-click events fall through to the
    // picker (which still owns the dbl-click → atomPicked path).
    cameraInputFilter_ = new CameraInputFilter(vtkWidget_, scene_,
                                                 scene_->cameraComposer(), this);

    // Atom Info dock — tabified on the LEFT alongside Selection + Strip.
    // Tracks the selection's FOCUS atom (one atom's full per-frame pile).
    // Starts empty; fills in on the first pick.
    inspectorDock_ = new QtAtomInspectorDock(this);
    inspectorDock_->setContext(loaded_->protein.get(),
                                transformed_);
    addDockWidget(Qt::LeftDockWidgetArea, inspectorDock_);

    ACONNECT(playback_,  &QtPlaybackController::frameChanged,
             inspectorDock_, &QtAtomInspectorDock::setFrame);

    // ---- Selection model — the single source of selection truth ----------
    //
    // The picker reports a pick + its keyboard modifiers; AtomSelection
    // interprets the gesture (plain = replace the focus; Shift = toggle the
    // atom in the ≤4 ordered set) and is itself the QAbstractListModel the
    // SelectionDock view binds to. It fans typed changes out:
    //   focusChanged → Atom Info retargets to the focus atom;
    //   cleared      → Atom Info clears;
    //   changed      → the measurement overlay rebuilds, and the scene
    //                  refreshes the current frame so the spheres reposition.
    selection_ = new model::AtomSelection(loaded_->protein.get(), this);

    signalCatalog_ = new model::TrajectorySignalCatalog(this);
    auto fieldAvailability = std::make_shared<model::TrajectoryFieldAvailability>(
        model::TrajectoryFieldAvailability::Build(loaded_->conformation.get(),
                                                  signalCatalog_->allDescriptorList()));
    signalCatalog_->setFieldAvailability(fieldAvailability);
    visualizationContext_ = {};
    visualizationContext_.availability = fieldAvailability.get();
    visualizationContext_.hasTrajectory = loaded_->conformation
        && loaded_->conformation->asTrajectory() != nullptr;
    visualizationContext_.tensorGlyphGestureEnabled = false;
    const QStringList unresolvedModes =
        model::VisualizationRegistry::instance().unresolvedStaticModes(*signalCatalog_);
    for (const QString& mode : unresolvedModes) {
        qCWarning(diagnostics::cDash).noquote()
            << QStringLiteral("event=viz_unresolved_static_mode mode=%1").arg(mode);
    }
    Q_ASSERT(unresolvedModes.isEmpty());
    inspectorDock_->setFieldAvailability(fieldAvailability);
    dashboardSignals_ = new model::DashboardSignalModel(this);
    dashboardSignals_->setFieldAvailability(fieldAvailability);
    dashboardPanels_ = new model::DashboardPanelModel(this);
    dashboardSelectionController_ =
        new DashboardSelectionController(signalCatalog_, dashboardSignals_, dashboardPanels_, this);
    lastDashboardSelectedCount_ = dashboardSelectionController_->selectedCount();
    signalDisplayDialog_ = new SignalDisplayDialog(this);
    signalDisplayDialog_->setTrajectorySignalCatalog(signalCatalog_);
    signalDisplayDialog_->setDashboardSignalModel(dashboardSignals_);
    signalDisplayDialog_->setDashboardPanelModel(dashboardPanels_);
    signalDisplayDialog_->setDashboardSelectionController(dashboardSelectionController_.data());
    signalDisplayDialog_->setContext(loaded_->protein.get(), transformed_);
    signalDisplayDialog_->setVisualizationContext(visualizationContext_);
    signalDisplayDialog_->setSelection(selection_);
    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             signalDisplayDialog_, &SignalDisplayDialog::setFrame);

    ACONNECT(picker_,    &QtAtomPicker::atomPicked,
             selection_, &model::AtomSelection::applyPick);
    ACONNECT(picker_, &QtAtomPicker::atomPicked,
             scene_,  &MoleculeScene::clearReveal);
    // Tag the render scheduler so the EndEvent observer logs source=picker
    // for the render that follows. selection_->applyPick triggers
    // refreshCurrentFrame which itself calls requestRender(Timer);
    // tagging Picker afterward overrides the source (requestRender is
    // coalescing — lastRenderSource_ is last-writer-wins, the queued
    // paint hasn't fired yet within this synchronous signal handling).
    ACONNECT(picker_, &QtAtomPicker::atomPicked,
             this,   [this](std::size_t, Qt::KeyboardModifiers) {
                 if (scene_) scene_->requestRender(
                     MoleculeScene::RenderSource::Picker);
             });

    ACONNECT(selection_, &model::AtomSelection::focusChanged,
             inspectorDock_, &QtAtomInspectorDock::setPickedAtom);
    ACONNECT(selection_, &model::AtomSelection::cleared,
             inspectorDock_, &QtAtomInspectorDock::clearSelection);
    const auto updateMetricAction = [this]() {
        if (signalDisplaysAction_)
            signalDisplaysAction_->setEnabled(selection_ && selection_->hasFocus());
    };
    ACONNECT(selection_, &model::AtomSelection::focusChanged, this, [this, updateMetricAction](std::size_t) {
        updateMetricAction();
        updateCameraModeActions();
    });
    ACONNECT(selection_, &model::AtomSelection::cleared, this, [this, updateMetricAction]() {
        updateMetricAction();
        updateCameraModeActions();
    });
    updateMetricAction();
    updateCameraModeActions();

    if (auto* meas = scene_->measurementOverlay()) {
        meas->setSelection(selection_);
        ACONNECT(selection_, &model::AtomSelection::changed,
                 meas,       &MeasurementOverlay::onSelectionChanged);
    }
    ACONNECT(selection_, &model::AtomSelection::changed,
             this, [this]() {
                 // Plane-lock release-on-selection-change is the documented
                 // behaviour (see CameraMode.h "lock release semantics" —
                 // Plane releases, Atom/Bond/Dihedral/Subset stay). The
                 // composer owns mode state, so explicitly drop the plane
                 // lock here only if it's currently active.
                 if (scene_ && scene_->cameraComposer()
                     && scene_->cameraComposer()->mode().kind
                            == CameraMode::Kind::Plane) {
                     const std::size_t t = playback_
                         ? static_cast<std::size_t>(playback_->currentFrame()) : 0u;
                     scene_->cameraComposer()->setMode(FreeMode(), FreePolicy(), t);
                 }
                 updateCameraModeActions();
                 if (scene_) scene_->refreshCurrentFrame();
             });

    // Selected-atoms panel — the QListView bound to the AtomSelection model
    // (slot colour swatch + residue:atom label + geometry kind). Tabified
    // with Atom Info in the left dock area.
    selectionDock_ = new SelectionDock(this);
    selectionDock_->setModel(selection_);
    addDockWidget(Qt::LeftDockWidgetArea, selectionDock_);
    tabifyDockWidget(inspectorDock_, selectionDock_);

    // Dashboard strips — active signals from SignalDisplayDialog rendered
    // through one shared strip surface and the shared TimeViewportController.
    // Tabified with Inspector + Selection on the left so the central
    // viewport gets the widest stable real estate; the user toggles each
    // panel via the toolbar buttons added below.
    dashboardStripDock_ = new DashboardStripDock(this);
    dashboardStripDock_->setContext(loaded_->protein.get(), transformed_);
    dashboardStripDock_->setSignalModels(signalCatalog_, dashboardSignals_);
    dashboardStripDock_->setPanelModel(dashboardPanels_);
    dashboardStripDock_->setSelectionController(dashboardSelectionController_.data());
    dashboardStripDock_->setSelection(selection_);
    dashboardStripDock_->setTimeViewport(timeViewport_);
    dashboardController_ = dashboardStripDock_->displayController();
    if (dashboardController_)
        dashboardController_->setVisualizationContext(visualizationContext_);
    addDockWidget(Qt::LeftDockWidgetArea, dashboardStripDock_);
    tabifyDockWidget(inspectorDock_, dashboardStripDock_);
    inspectorDock_->raise();
    resizeDocks({inspectorDock_}, {360}, Qt::Horizontal);

    // Start clean — no property docks open on launch. The user opens each via
    // the toolbar toggles (toggleViewAction) below; picking an atom updates the
    // Inspector's contents but does not force-show the dock. kSettingsVersion
    // was bumped so a stale "docks visible" layout is discarded by restoreState
    // rather than re-opening them; window geometry (restoreGeometry) persists.
    inspectorDock_->setVisible(false);
    selectionDock_->setVisible(false);
    dashboardStripDock_->setVisible(false);
    ACONNECT(dashboardStripDock_, &QDockWidget::visibilityChanged,
             this, [this](bool visible) {
                 qCInfo(diagnostics::cDash).noquote()
                     << QStringLiteral("event=dock_visibility_changed visible=%1 width=%2")
                            .arg(visible ? 1 : 0)
                            .arg(dashboardDockWidth());
             });
    ACONNECT(dashboardSelectionController_.data(),
             &DashboardSelectionController::selectedCountChanged,
             this,
             [this](int count) {
                 ASSERT_THREAD(this);
                 const bool added = count > lastDashboardSelectedCount_;
                 lastDashboardSelectedCount_ = count;
                 if (!added || !dashboardStripDock_ || dashboardStripDock_->isVisible())
                     return;
                 qCInfo(diagnostics::cDash).noquote()
                     << QStringLiteral("event=dock_reveal_on_add count=%1").arg(count);
                 revealDockQueued(dashboardStripDock_);
             });

    // Panel recovery — one QAction per dock from QDockWidget::toggleViewAction().
    // The same QAction instances live in View -> Panels and the toolbar menu,
    // so checked state stays correct no matter how the user toggles a panel.
    const struct { QDockWidget* dock; const char* label; const char* tip; } kPanels[] = {
        { inspectorDock_,      "Inspector",
          "Show or hide the Atom Info panel." },
        { selectionDock_,      "Selection",
          "Show or hide the Selected Atoms panel." },
        { dashboardStripDock_, "Strip",
          "Show or hide the time-series strip dock." },
    };
    if (panelsMenu_) {
        panelsMenu_->clear();
        for (const auto& p : kPanels) {
            if (!p.dock) continue;
            QAction* a = p.dock->toggleViewAction();
            a->setText(QString::fromUtf8(p.label));
            a->setToolTip(QString::fromUtf8(p.tip));
            ACONNECT(a, &QAction::triggered,
                     this, [this, dock = QPointer<QDockWidget>(p.dock)](bool checked) {
                         if (checked && dock)
                             revealDockQueued(dock.data());
                     });
            panelsMenu_->addAction(a);
        }
    }
    if (toolsToolbar_ && panelsMenu_) {
        toolsToolbar_->addSeparator();
        panelsButton_ = new QToolButton(toolsToolbar_);
        panelsButton_->setText(QStringLiteral("Panels"));
        panelsButton_->setToolTip(QStringLiteral("Show or hide panels."));
        panelsButton_->setPopupMode(QToolButton::InstantPopup);
        panelsButton_->setMenu(panelsMenu_);
        toolsToolbar_->addWidget(panelsButton_);
    }

    ACONNECT(dashboardStripDock_, &DashboardStripDock::revealRequested,
             scene_,              &MoleculeScene::revealBinding);
    // L-3a (2026-05-29): expose the scene's reveal overlay to the
    // dashboard controller so static.tensor mode on Reorient
    // orientation_tensor signals fires an ellipsoid glyph in the
    // 3-D view.
    if (scene_ && scene_->revealOverlay()) {
        dashboardStripDock_->setSceneOverlay(scene_->revealOverlay());
        visualizationContext_.hasSceneOverlay = true;
        if (dashboardController_)
            dashboardController_->setVisualizationContext(visualizationContext_);
        if (signalDisplayDialog_)
            signalDisplayDialog_->setVisualizationContext(visualizationContext_);
    }
    ACONNECT(dashboardStripDock_, &DashboardStripDock::metricPickerRequested,
             this,                &ReaderMainWindow::onOpenSignalDisplays);
    ACONNECT(playback_,           &QtPlaybackController::frameChanged,
             dashboardStripDock_, &DashboardStripDock::setFrame);
    if (frameSlider_) {
        ACONNECT(frameSlider_.data(), &QSlider::sliderPressed,
                 this, [this]() {
                     if (dashboardController_)
                         dashboardController_->setScrubActive(true);
                 });
        ACONNECT(frameSlider_.data(), &QSlider::sliderReleased,
                 this, [this]() {
                     if (dashboardController_)
                         dashboardController_->setScrubActive(false);
                 });
    }

    // DFT shielding campaign (optional): make the frame-local source
    // available to descriptor-family samplers. The `.LGS` carries the
    // typed `dft.frames[]` map — frame_index → meta.json — so the
    // store builds straight from it (no dir scanning, no name parsing).
    if (loaded_->manifest.dft.has_value()) {
        const auto& dft = *loaded_->manifest.dft;
        dftStore_ = new model::DftShieldingStore(loaded_->protein.get(), dft.frames, this);
        dashboardStripDock_->setDftStore(dftStore_);
        visualizationContext_.hasDftStore = true;
        if (dashboardController_)
            dashboardController_->setVisualizationContext(visualizationContext_);
        if (signalDisplayDialog_)
            signalDisplayDialog_->setVisualizationContext(visualizationContext_);
        qCInfo(cWindow).noquote() << "DFT shielding store wired from .LGS |"
                                  << "method=" << dft.method
                                  << "| frames=" << dftStore_->jobCount()
                                  << "| campaign_target=" << dft.campaign_target_frames;
    }

    // Mutant-pair alternate-pose action — when the manifest says we
    // auto-opened WT, expose a File menu action that launches a fresh
    // process on the ALA `.LGS`. Spawned the same way as Recent files.
    if (loaded_->manifest.kind == h5reader::io::CalcsetManifest::Kind::MutantPair
        && loaded_->manifest.mutant_pair.has_value()
        && !loaded_->manifest.mutant_pair->ala_lgs_abspath.isEmpty()) {
        const QString alt = loaded_->manifest.mutant_pair->ala_lgs_abspath;
        if (QMenuBar* mb = menuBar()) {
            for (QAction* a : mb->actions()) {
                if (a->menu() && a->text() == QStringLiteral("&File")) {
                    QAction* switchAct = a->menu()->addAction(
                        QStringLiteral("Open mutant alternate (ALA)…"));
                    switchAct->setToolTip(QStringLiteral(
                        "This run is a mutant pair; WT is opened in this window. "
                        "Click to launch a separate reader on the ALA pose: %1").arg(alt));
                    ACONNECT(switchAct, &QAction::triggered, this, [this, alt]() {
                        openRecentPath(alt);
                    });
                    break;
                }
            }
        }
    }

    // Initial status bar population.
    onFrameChanged(0);

    // Default size — wide enough for the playback + camera + transform +
    // metrics + overlays + panel controls to fit in one
    // toolbar without Qt's overflow chevron. QSettings restore overrides
    // this on later launches.
    resize(1600, 900);
    setWindowTitle(QStringLiteral("h5-reader — %1").arg(loaded_->proteinId));

    // QSettings restore — geometry, dock state, log mask, recent menu.
    // Tolerant: missing or version-mismatched blobs leave the ctor's
    // explicit defaults intact. Runs AFTER all docks/toolbars exist so
    // restoreState has named docks to bind to.
    restoreAllSettings();
    resetDashboardStateForRunLoad();
    if (!loaded_->runPath.isEmpty())
        addToRecentFiles(QDir(loaded_->runPath).absolutePath());

    qCInfo(cWindow).noquote() << "ctor done";
}

void ReaderMainWindow::showEvent(QShowEvent* event) {
    QMainWindow::showEvent(event);
    if (glInfoLogged_) return;
    glInfoLogged_ = true;

    // The VTK widget hasn't painted yet at showEvent time, so the GL
    // context isn't current. Defer the query to the next event-loop
    // tick — by then the first frame has rendered and ReportCapabilities
    // returns a populated string. Single-shot, never repeats.
    QPointer<ReaderMainWindow> self(this);
    QTimer::singleShot(0, this, [self]() {
        if (!self || !self->renderWindow_) return;
        const QString caps =
            QString::fromUtf8(self->renderWindow_->ReportCapabilities());
        const QStringList wanted = {
            QStringLiteral("OpenGL vendor"),
            QStringLiteral("OpenGL renderer"),
            QStringLiteral("OpenGL version"),
            QStringLiteral("OpenGL vendor-specific"),
        };
        for (const QString& line : caps.split(QChar('\n'))) {
            for (const QString& key : wanted) {
                if (line.contains(key, Qt::CaseInsensitive)) {
                    qCInfo(cWindow).noquote() << "GL:" << line.trimmed();
                    break;
                }
            }
        }
    });
}

ReaderMainWindow::~ReaderMainWindow() {
    // Most cleanup runs in shutdown(). The destructor only handles the
    // pathological case where shutdown() was never called (e.g. window
    // deleted outside the normal quit flow).
    if (!shutdownDone_) {
        qCWarning(cWindow).noquote()
            << "destructor called without prior shutdown(); running now";
        shutdown();
    }
}


quint16 ReaderMainWindow::startRestServer(quint16 port) {
    ASSERT_THREAD(this);
    if (!loaded_ || !loaded_->protein || !loaded_->conformation) {
        qCCritical(cWindow).noquote() << "REST start refused: loader result not wired";
        return 0;
    }
    if (restServer_) {
        qCWarning(cWindow).noquote() << "REST server already running; ignoring re-start";
        return 0;
    }
    restServer_ = new RestServer(this);
    restServer_->setContext(scene_,
                            selection_,
                            dashboardSignals_,
                            dashboardPanels_,
                            dashboardSelectionController_.data(),
                            dashboardController_.data(),
                            signalCatalog_,
                            playback_,
                            loaded_.get(),
                            this,
                            this,
                            transformed_);
    const quint16 bound = restServer_->listen(port);
    if (bound == 0) {
        qCCritical(cWindow).noquote() << "REST server failed to bind port" << port;
        restServer_->deleteLater();
        restServer_ = nullptr;
    }
    return bound;
}

void ReaderMainWindow::setDocksVisible(bool visible) {
    ASSERT_THREAD(this);

    // Hide path: stash each dock's pre-hide visibility so a later restore
    // can return individually-user-hidden docks to their hidden state.
    // No-op if already hidden (don't double-stash).
    if (!visible) {
        if (docksHidden_)
            return;
        stashedDockVisibility_.clear();
        const std::vector<QDockWidget*> docks = {
            inspectorDock_, selectionDock_, dashboardStripDock_
        };
        for (QDockWidget* d : docks) {
            if (!d) continue;
            stashedDockVisibility_.push_back({QPointer<QDockWidget>(d), d->isVisible()});
            d->setVisible(false);
        }
        docksHidden_ = true;
        qCInfo(cWindow).noquote()
            << "docks hidden | count=" << stashedDockVisibility_.size();
        return;
    }

    // Restore path: walk the stash; QPointer-safe iteration in case any
    // dock was destroyed since the hide. Each dock returns to its stashed
    // visibility, so a dock that was already hidden before the harness
    // requested hide stays hidden.
    if (!docksHidden_)
        return;
    for (const DockVis& dv : stashedDockVisibility_) {
        if (dv.dock)
            dv.dock->setVisible(dv.wasVisible);
    }
    qCInfo(cWindow).noquote()
        << "docks restored | count=" << stashedDockVisibility_.size();
    stashedDockVisibility_.clear();
    docksHidden_ = false;
}

void ReaderMainWindow::setDashboardDockVisible(bool visible) {
    ASSERT_THREAD(this);
    if (!dashboardStripDock_)
        return;
    if (visible) {
        revealDockQueued(dashboardStripDock_);
    } else {
        dashboardStripDock_->setVisible(false);
    }
}

bool ReaderMainWindow::dashboardDockVisible() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ && dashboardStripDock_->isVisible();
}

int ReaderMainWindow::dashboardDockWidth() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ ? dashboardStripDock_->width() : 0;
}

bool ReaderMainWindow::dashboardDockRaised() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ && dashboardStripDock_->isVisible()
        && !dashboardStripDock_->visibleRegion().isEmpty();
}

int ReaderMainWindow::dashboardOwnedPanelCount() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ ? dashboardStripDock_->ownedPanelCount() : 0;
}

int ReaderMainWindow::dashboardStripTrackCount() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ ? dashboardStripDock_->stripTrackCount() : 0;
}

bool ReaderMainWindow::openSignalDisplayPicker(QString* blockedReason) {
    ASSERT_THREAD(this);
    if (blockedReason)
        blockedReason->clear();
    if (!signalDisplayDialog_) {
        if (blockedReason)
            *blockedReason = QStringLiteral("signal display dialog not wired");
        return false;
    }
    if (!selection_ || !selection_->hasFocus()) {
        if (blockedReason)
            *blockedReason = QStringLiteral("Metrics action is disabled because AtomSelection has no focus.");
        return false;
    }
    if (signalDisplaysAction_ && !signalDisplaysAction_->isEnabled()) {
        if (blockedReason)
            *blockedReason = QStringLiteral("Metrics action is disabled.");
        return false;
    }
    if (playback_)
        signalDisplayDialog_->setFrame(playback_->currentFrame());
    signalDisplayDialog_->refreshCatalog();
    signalDisplayDialog_->show();
    signalDisplayDialog_->raise();
    signalDisplayDialog_->activateWindow();
    return true;
}

QJsonObject ReaderMainWindow::signalDisplayPickerState() const {
    ASSERT_THREAD(this);
    if (!signalDisplayDialog_)
        return QJsonObject{{"open", false}};
    return signalDisplayDialog_->pickerState();
}

QJsonObject ReaderMainWindow::addSelectedSignalFromPicker() {
    ASSERT_THREAD(this);
    if (signalDisplayDialog_)
        signalDisplayDialog_->onAddSelected();
    return signalDisplayPickerState();
}

void ReaderMainWindow::revealDockQueued(QDockWidget* dock) {
    ASSERT_THREAD(this);
    if (!dock)
        return;
    dock->setVisible(true);
    QTimer::singleShot(0, this, [this, dock = QPointer<QDockWidget>(dock)]() {
        if (!dock)
            return;
        resizeDocks({dock.data()}, {360}, Qt::Horizontal);
        dock->raise();
    });
}

void ReaderMainWindow::resetDashboardStateForRunLoad() {
    ASSERT_THREAD(this);
    if (dashboardSelectionController_) {
        dashboardSelectionController_->clearAllMetrics();
        lastDashboardSelectedCount_ = dashboardSelectionController_->selectedCount();
    }
    if (dashboardStripDock_)
        dashboardStripDock_->setVisible(false);
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=selection_reset_on_load count=%1 dock_visible=%2")
               .arg(dashboardSelectionController_ ? dashboardSelectionController_->selectedCount() : 0)
               .arg(dashboardStripDock_ && dashboardStripDock_->isVisible() ? 1 : 0);
}

void ReaderMainWindow::shutdown() {
    ASSERT_THREAD(this);
    if (shutdownDone_) return;
    shutdownDone_ = true;

    qCInfo(cWindow).noquote() << "shutdown entered";

    // Per spec/viewport_pipeline_2026-05-30.md §4.4:
    //
    // 1. Stop the REST server SYNCHRONOUSLY. The /shutdown endpoint
    //    fires from a request handler; the server needs to drain
    //    before timers stop so a follow-up request can't trigger a
    //    race with timer teardown.
    if (restServer_) {
        // RestServer doesn't expose stopListening(); the QHttpServer
        // owned by it tears down when the RestServer is deleted, but
        // deleteLater on shutdown is enough for this path because
        // aboutToQuit drains the event loop afterwards. We do hold a
        // direct pointer; do a synchronous delete here.
        delete restServer_;
        restServer_ = nullptr;
    }

    // 2. Stop every timer owned by us or our children. The generic
    //    findChildren sweep catches QtPlaybackController's timer too.
    const auto timers = findChildren<QTimer*>();
    for (auto* timer : timers) {
        if (timer->isActive()) timer->stop();
    }

    // 3. Detach the render window from the widget BEFORE dropping our
    //    smart pointer. setRenderWindow(nullptr) makes the context
    //    current and calls Finalize on the old render window via the
    //    adapter's destructor (QVTKRenderWindowAdapter.cxx:150-166).
    //    The explicit renderWindow_->Finalize() that used to live here
    //    is gone — doing it AFTER detaching the widget left the adapter
    //    holding a destroyed window for the brief moment between the
    //    two calls.
    if (vtkWidget_) {
        vtkWidget_->setRenderWindow(static_cast<vtkGenericOpenGLRenderWindow*>(nullptr));
    }

    qCInfo(cWindow).noquote() << "shutdown done";
}

void ReaderMainWindow::buildUi() {
    vtkWidget_    = new QVTKOpenGLNativeWidget(this);
    renderWindow_ = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkWidget_->setRenderWindow(renderWindow_);
    setCentralWidget(vtkWidget_);

    // File ▸ Open Directory… — point the reader at a run directory (a
    // trajectory or a single pose) or a trajectory.h5. Launches a fresh
    // reader process on the chosen path (multiple-instance safe).
    auto* fileMenu = menuBar()->addMenu(QStringLiteral("&File"));
    auto* openDirAct = fileMenu->addAction(QStringLiteral("Open Directory…"));
    openDirAct->setShortcut(QKeySequence::Open);
    ACONNECT(openDirAct, &QAction::triggered, this, &ReaderMainWindow::onOpenDirectory);

    // File ▸ Recent — populated from QSettings during restoreAllSettings.
    // Empty until then; each entry launches a fresh reader on click.
    recentMenu_ = fileMenu->addMenu(QStringLiteral("&Recent"));
    recentMenu_->setObjectName(QStringLiteral("RecentMenu"));

    auto* viewMenu = menuBar()->addMenu(QStringLiteral("&View"));
    panelsMenu_ = viewMenu->addMenu(QStringLiteral("&Panels"));
    panelsMenu_->setObjectName(QStringLiteral("PanelsMenu"));
}

void ReaderMainWindow::buildToolbar() {
    auto* tb = addToolBar(QStringLiteral("Playback"));
    tb->setObjectName(QStringLiteral("PlaybackToolbar"));
    tb->setMovable(false);
    playbackToolbar_ = tb;
    QFont toolbarFont = tb->font();
    if (toolbarFont.pointSize() > 8)
        toolbarFont.setPointSize(toolbarFont.pointSize() - 1);
    else if (toolbarFont.pixelSize() > 10)
        toolbarFont.setPixelSize(toolbarFont.pixelSize() - 1);
    tb->setFont(toolbarFont);

    playAction_ = tb->addAction(
        style()->standardIcon(QStyle::SP_MediaPlay),
        QStringLiteral("Play / Pause"));
    ACONNECT(playAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onPlayPauseClicked);

    auto* stepBack = tb->addAction(
        style()->standardIcon(QStyle::SP_MediaSeekBackward),
        QStringLiteral("Step back"));
    auto* stepFwd  = tb->addAction(
        style()->standardIcon(QStyle::SP_MediaSeekForward),
        QStringLiteral("Step forward"));

    tb->addSeparator();

    frameSlider_ = new QSlider(Qt::Horizontal, tb);
    frameSlider_->setMinimumWidth(400);
    tb->addWidget(frameSlider_);

    tb->addSeparator();
    tb->addWidget(new QLabel(QStringLiteral("fps"), tb));
    fpsSpinner_ = new QSpinBox(tb);
    fpsSpinner_->setSuffix(QStringLiteral(" /s"));
    tb->addWidget(fpsSpinner_);

    addToolBarBreak();
    tb = addToolBar(QStringLiteral("Tools"));
    tb->setObjectName(QStringLiteral("ToolsToolbar"));
    tb->setMovable(false);
    toolsToolbar_ = tb;
    tb->setFont(toolbarFont);

    tb->addSeparator();

    // Camera-mode action group — Focus / Newman / Plane lock / Free as
    // mutually-exclusive radio actions. QActionGroup with exclusive=true is
    // the standard Qt idiom; the visual checked state is the union of
    // user clicks and the composer's modeChanged signal (updateCameraModeActions
    // syncs to composer->mode().kind). Each action uses QAction::triggered
    // (user-only) not QAction::toggled (also fires on programmatic setChecked)
    // so the sync loop is closed.
    cameraModeGroup_ = new QActionGroup(this);
    cameraModeGroup_->setExclusive(true);

    focusAction_ = tb->addAction(QStringLiteral("Focus"));
    focusAction_->setCheckable(true);
    focusAction_->setEnabled(false);
    focusAction_->setToolTip(QStringLiteral(
        "Focus the camera on the selected atom's residue backbone plane (N, CA, C). "
        "Requires a focused atom."));
    cameraModeGroup_->addAction(focusAction_);
    ACONNECT(focusAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onFocusCameraTriggered);

    newmanAction_ = tb->addAction(QStringLiteral("Newman"));
    newmanAction_->setCheckable(true);
    newmanAction_->setEnabled(false);
    newmanAction_->setToolTip(QStringLiteral(
        "Newman projection — sight down the central bond of a 4-atom dihedral. "
        "Requires exactly 4 selected atoms."));
    cameraModeGroup_->addAction(newmanAction_);
    ACONNECT(newmanAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onNewmanProjectionTriggered);

    planeLockAction_ = tb->addAction(QStringLiteral("Plane lock"));
    planeLockAction_->setCheckable(true);
    planeLockAction_->setEnabled(false);
    planeLockAction_->setToolTip(QStringLiteral(
        "Keep the view centred and oriented to the plane defined by exactly three selected atoms."));
    cameraModeGroup_->addAction(planeLockAction_);
    ACONNECT(planeLockAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onPlaneLockTriggered);

    freeAction_ = tb->addAction(QStringLiteral("Free"));
    freeAction_->setCheckable(true);
    freeAction_->setChecked(true);   // composer's default mode is Free
    freeAction_->setToolTip(QStringLiteral(
        "Release any camera lock; mouse drag controls the view directly."));
    cameraModeGroup_->addAction(freeAction_);
    ACONNECT(freeAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onFreeCameraTriggered);

    tb->addSeparator();

    transformFitAction_ = tb->addAction(QStringLiteral("Mode: Locked backbone  ⇄"));
    transformFitAction_->setToolTip(fitModeToolTip());
    ACONNECT(transformFitAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onTransformFitClicked);

    // Harness marker preset. Kept as an action for the existing slot/REST
    // path, hidden from the normal toolbar.
    instrumentAction_ = tb->addAction(QStringLiteral("Harness marker"));
    instrumentAction_->setCheckable(true);
    instrumentAction_->setToolTip(QStringLiteral(
        "Enable the marker preset on the focus atom."));
    instrumentAction_->setVisible(false);
    ACONNECT(instrumentAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onInstrumentToggled);

    tb->addSeparator();

    signalDisplaysAction_ = tb->addAction(QStringLiteral("Metrics..."));
    signalDisplaysAction_->setEnabled(false);
    signalDisplaysAction_->setToolTip(QStringLiteral("Select a nearby atom or residue and add a metric display."));
    ACONNECT(signalDisplaysAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onOpenSignalDisplays);

    tb->addSeparator();

    // Overlay toggles — Ribbon, Rings, Butterfly. The scene + overlays
    // are constructed after the toolbar, so we defer connections until
    // after ctor via a zero-delay singleShot.
    showRibbonAction_ = tb->addAction(QStringLiteral("Ribbon"));
    showRibbonAction_->setCheckable(true);
    showRibbonAction_->setChecked(true);
    showRibbonAction_->setToolTip(QStringLiteral(
        "Backbone ribbon; secondary structure driven by per-frame DSSP."));

    showRingsAction_ = tb->addAction(QStringLiteral("Rings"));
    showRingsAction_->setCheckable(true);
    showRingsAction_->setChecked(true);
    showRingsAction_->setToolTip(QStringLiteral(
        "Aromatic ring polygons + normal arrows (per-frame ring_geometry)."));

    showButterflyAction_ = tb->addAction(QStringLiteral("Butterfly"));
    showButterflyAction_->setCheckable(true);
    showButterflyAction_->setChecked(false);   // off by default — expensive
    showButterflyAction_->setToolTip(QStringLiteral(
        "BS / HM volumetric isosurfaces around each aromatic ring. "
        "Re-evaluates closed-form kernel per frame on a 20³ grid."));

    showBFieldAction_ = tb->addAction(QStringLiteral("B-field"));
    showBFieldAction_->setCheckable(true);
    showBFieldAction_->setChecked(false);   // off by default — expensive
    showBFieldAction_->setToolTip(QStringLiteral(
        "Biot-Savart B-field streamlines around each aromatic ring, "
        "seeded on a circle at 1.5× ring radius, coloured by |B|."));

    // Wire step actions + overlay toggles once playback_ / scene_ exist.
    // Toolbar must be constructed before those objects because Qt needs
    // the action parent (the toolbar), but the connections need the
    // recipients; a zero-delay singleShot schedules us for the next
    // event-loop tick.
    QTimer::singleShot(0, this, [this, stepBack, stepFwd]() {
        if (!playback_) return;
        ACONNECT(stepBack, &QAction::triggered,
                 playback_, &QtPlaybackController::stepBackward);
        ACONNECT(stepFwd, &QAction::triggered,
                 playback_, &QtPlaybackController::stepForward);

        if (!scene_) return;

        // Ribbon / Rings — cheap visibility flips. No kernel eval to
        // rerun; just render the current scene with the new actor
        // visibility.
        ACONNECT(showRibbonAction_.data(), &QAction::toggled,
                 this, [this](bool on) {
                     if (!scene_ || !scene_->ribbonOverlay()) return;
                     scene_->ribbonOverlay()->setVisible(on);
                     scene_->requestRender();
                 });
        ACONNECT(showRingsAction_.data(), &QAction::toggled,
                 this, [this](bool on) {
                     if (!scene_ || !scene_->ringPolygonOverlay()) return;
                     scene_->ringPolygonOverlay()->setVisible(on);
                     scene_->requestRender();
                 });

        // Butterfly — the field-grid overlay skips kernel re-eval while
        // hidden. Turning on requires a refresh so its scalar arrays
        // populate for the current frame; turning off just needs a
        // render to flush the hidden actors.
        ACONNECT(showButterflyAction_.data(), &QAction::toggled,
                 this, [this](bool on) {
                     if (!scene_ || !scene_->fieldGridOverlay()) return;
                     scene_->fieldGridOverlay()->setVisible(on);
                     if (on) scene_->refreshCurrentFrame();
                     else    scene_->requestRender();
                 });

        // B-field streamlines — same expensive-when-visible pattern as
        // the butterfly. The overlay's vtkStreamTracer only runs once
        // the structured-grid input has non-zero vectors, which UpdateRing
        // fills in via the kernel eval during refreshCurrentFrame.
        ACONNECT(showBFieldAction_.data(), &QAction::toggled,
                 this, [this](bool on) {
                     if (!scene_ || !scene_->bfieldStreamOverlay()) return;
                     scene_->bfieldStreamOverlay()->setVisible(on);
                     if (on) scene_->refreshCurrentFrame();
                     else    scene_->requestRender();
                 });
    });
}

void ReaderMainWindow::buildStatusBar() {
    proteinLabel_ = new QLabel(loaded_->proteinId, this);
    frameLabel_   = new QLabel(QStringLiteral("frame —"), this);
    timeLabel_    = new QLabel(QStringLiteral("t=— ps"), this);

    statusBar()->addPermanentWidget(proteinLabel_);
    statusBar()->addPermanentWidget(frameLabel_);
    statusBar()->addPermanentWidget(timeLabel_);
}

void ReaderMainWindow::onFrameChanged(int t) {
    ASSERT_THREAD(this);
    const int T = static_cast<int>(loaded_->conformation->frameCount());
    const double t_ps = loaded_->conformation->timePicoseconds(
        static_cast<size_t>(std::clamp(t, 0, T - 1)));

    if (frameLabel_) {
        frameLabel_->setText(QStringLiteral("frame %1 / %2").arg(t + 1).arg(T));
    }
    if (timeLabel_) {
        timeLabel_->setText(QStringLiteral("t=%1 ps").arg(t_ps, 0, 'f', 1));
    }
    if (frameSlider_ && frameSlider_->value() != t) {
        const QSignalBlocker block(frameSlider_);
        frameSlider_->setValue(t);
    }
}

void ReaderMainWindow::updateCameraModeActions() {
    // Gating — what each action requires from the current selection.
    const bool hasFocus  = selection_ && selection_->hasFocus();
    const std::size_t n  = selection_ ? selection_->count() : 0;
    if (focusAction_)     focusAction_->setEnabled(scene_ && hasFocus);
    if (newmanAction_)    newmanAction_->setEnabled(scene_ && n == 4);
    if (planeLockAction_) planeLockAction_->setEnabled(scene_ && n == 3);
    if (freeAction_)      freeAction_->setEnabled(scene_ != nullptr);

    // Visual checked state — sourced from the composer. Programmatic
    // setChecked here would fire QAction::toggled but we connected via
    // QAction::triggered (user-only), so no loop. Use a signal blocker
    // anyway since QActionGroup itself emits triggered on exclusive change.
    if (!scene_ || !scene_->cameraComposer())
        return;
    const auto kind = scene_->cameraComposer()->mode().kind;
    const auto setOne = [](QAction* a, bool on) {
        if (!a) return;
        const QSignalBlocker block(a);
        a->setChecked(on);
    };
    setOne(focusAction_,     false);
    setOne(newmanAction_,    false);
    setOne(planeLockAction_, false);
    setOne(freeAction_,      false);
    switch (kind) {
        case CameraMode::Kind::Plane:
            setOne(planeLockAction_, true); break;
        case CameraMode::Kind::Dihedral:
            // Newman is the only dihedral path we expose in the toolbar; a
            // dihedral mode that came in via REST also shows here.
            setOne(newmanAction_, true); break;
        case CameraMode::Kind::Free:
            setOne(freeAction_, true); break;
        case CameraMode::Kind::Atom:
            break;
        case CameraMode::Kind::Bond:
            break;
        case CameraMode::Kind::Subset:
            break;
    }
}

void ReaderMainWindow::updateFitModeLabel() {
    if (!transformFitAction_ || !transformed_)
        return;

    using TMode = h5reader::model::TransformedConformation::Mode;
    transformFitAction_->setText(transformed_->mode() == TMode::FitSubset
        ? QStringLiteral("Mode: Locked backbone  ⇄")
        : QStringLiteral("Mode: Kabsch with give  ⇄"));
    transformFitAction_->setToolTip(fitModeToolTip());
}

void ReaderMainWindow::onPlaneLockTriggered() {
    ASSERT_THREAD(this);
    if (!scene_ || !selection_ || selection_->count() != 3
        || !scene_->lockCameraToSelectionPlane(selection_->atoms())) {
        // setMode failed (degenerate); make sure the toolbar reflects
        // composer truth — likely Free or whatever was active before.
        updateCameraModeActions();
        return;
    }
    // The composer's modeChanged signal will fire updateCameraModeActions
    // for us; nothing else to do.
}

void ReaderMainWindow::onFocusCameraTriggered() {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->cameraComposer() || !selection_ || !selection_->hasFocus()) {
        updateCameraModeActions();
        return;
    }
    auto result = h5reader::app::DeriveFocusAnchor(*loaded_->protein,
                                                    selection_->focus(),
                                                    FocusAnchorKind::Plane);
    if (result.outcome != FocusAnchorOutcome::Ok) {
        qCWarning(cWindow).noquote()
            << "Focus camera: derive failed | atom=" << selection_->focus()
            << "| outcome=" << static_cast<int>(result.outcome);
        updateCameraModeActions();
        return;
    }
    const std::size_t t = playback_ ? static_cast<std::size_t>(playback_->currentFrame()) : 0u;
    scene_->cameraComposer()->setMode(result.mode, result.policy, t);
    (void)scene_->cameraComposer()->write(t);
    scene_->syncCameraClippingRange();
    scene_->requestRender(MoleculeScene::RenderSource::External);
    updateCameraModeActions();
}

void ReaderMainWindow::onNewmanProjectionTriggered() {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->cameraComposer() || !selection_ || selection_->count() != 4) {
        updateCameraModeActions();
        return;
    }
    const auto& a = selection_->atoms();
    CameraMode m = DihedralMode(a[0], a[1], a[2], a[3]);
    OrientationPolicy p = DownAxisPolicy(a[1], a[2]);
    const std::size_t t = playback_ ? static_cast<std::size_t>(playback_->currentFrame()) : 0u;
    scene_->cameraComposer()->setMode(m, p, t);
    (void)scene_->cameraComposer()->write(t);
    scene_->syncCameraClippingRange();
    scene_->requestRender(MoleculeScene::RenderSource::External);
    updateCameraModeActions();
}

void ReaderMainWindow::onFreeCameraTriggered() {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->cameraComposer()) {
        updateCameraModeActions();
        return;
    }
    const std::size_t t = playback_ ? static_cast<std::size_t>(playback_->currentFrame()) : 0u;
    scene_->cameraComposer()->setMode(FreeMode(), FreePolicy(), t);
}

void ReaderMainWindow::onTransformFitClicked() {
    ASSERT_THREAD(this);
    using TMode = h5reader::model::TransformedConformation::Mode;
    if (!transformed_ || !loaded_ || !loaded_->protein)
        return;

    if (transformed_->mode() == TMode::FitReference) {
        auto subset = h5reader::model::TransformedConformation::BackboneSubset(*loaded_->protein);
        if (subset.size() < 3) {
            qCWarning(cWindow).noquote()
                << "backbone fit requested but subset has <3 atoms; keeping all-atom fit";
            transformed_->setMode(TMode::FitReference, 0);
            return;
        }
        transformed_->setMode(TMode::FitSubset, 0, std::move(subset));
    } else {
        transformed_->setMode(TMode::FitReference, 0);
    }
}

void ReaderMainWindow::onInstrumentToggled(bool checked) {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->measurementOverlay()) return;
    scene_->measurementOverlay()->setInstrumentMode(checked, /*focusOnly=*/true);
    scene_->requestRender(MoleculeScene::RenderSource::External);
}

void ReaderMainWindow::closeEvent(QCloseEvent* event) {
    ASSERT_THREAD(this);
    saveAllSettings();
    // Accept unconditionally — a failed save is logged but not allowed
    // to trap the user inside the application. aboutToQuit fires the
    // existing shutdown() chain after this returns.
    event->accept();
}

void ReaderMainWindow::saveAllSettings() {
    ASSERT_THREAD(this);
    QSettings s;   // org/app names set in main_reader.cpp
    s.setValue(QStringLiteral("viewer/window/geometry"), saveGeometry());
    s.setValue(QStringLiteral("viewer/window/state"),
               saveState(kSettingsVersion));
    s.setValue(QStringLiteral("viewer/log/mask"),
               static_cast<uint>(h5reader::diagnostics::StructuredLogger::CategoryMask()));
    // Recent files list is write-through (addToRecentFiles writes
    // immediately) so no batch write here.
    qCInfo(cWindow).noquote() << "settings saved | mask="
                              << h5reader::diagnostics::StructuredLogger::CategoryMask();
}

void ReaderMainWindow::restoreAllSettings() {
    ASSERT_THREAD(this);
    QSettings s;
    const QByteArray geom = s.value(QStringLiteral("viewer/window/geometry")).toByteArray();
    if (!geom.isEmpty())
        restoreGeometry(geom);
    const QByteArray state = s.value(QStringLiteral("viewer/window/state")).toByteArray();
    if (!state.isEmpty())
        restoreState(state, kSettingsVersion);
    const QVariant maskVar = s.value(QStringLiteral("viewer/log/mask"));
    if (maskVar.isValid()) {
        bool ok = false;
        const uint mask = maskVar.toUInt(&ok);
        if (ok)
            h5reader::diagnostics::StructuredLogger::SetCategoryMask(mask);
    }
    const QStringList recent = s.value(QStringLiteral("viewer/recent/files")).toStringList();
    rebuildRecentFilesMenu(recent);
}

void ReaderMainWindow::addToRecentFiles(const QString& path) {
    ASSERT_THREAD(this);
    if (path.isEmpty()) return;
    QSettings s;
    QStringList recent = s.value(QStringLiteral("viewer/recent/files")).toStringList();
    recent.removeAll(path);
    recent.prepend(path);
    while (recent.size() > kMaxRecentFiles)
        recent.removeLast();
    s.setValue(QStringLiteral("viewer/recent/files"), recent);
    rebuildRecentFilesMenu(recent);
}

void ReaderMainWindow::rebuildRecentFilesMenu(const QStringList& paths) {
    ASSERT_THREAD(this);
    if (!recentMenu_) return;
    recentMenu_->clear();
    if (paths.isEmpty()) {
        QAction* empty = recentMenu_->addAction(QStringLiteral("(none)"));
        empty->setEnabled(false);
        return;
    }
    for (const QString& path : paths) {
        QAction* a = recentMenu_->addAction(path);
        ACONNECT(a, &QAction::triggered, this, [this, path]() {
            openRecentPath(path);
        });
    }
}

void ReaderMainWindow::openRecentPath(const QString& path) {
    ASSERT_THREAD(this);
    // Same launch pattern as onOpenDirectory — multiple-instance safe.
    const bool ok = QProcess::startDetached(QCoreApplication::applicationFilePath(),
                                             QStringList{path});
    if (!ok)
        qCWarning(cWindow).noquote()
            << "failed to launch a reader for recent path" << path;
}

void ReaderMainWindow::onPlayPauseClicked() {
    ASSERT_THREAD(this);
    if (playback_) playback_->togglePlayPause();
}

void ReaderMainWindow::onOpenDirectory() {
    ASSERT_THREAD(this);
    const QString dir = QFileDialog::getExistingDirectory(
        this, QStringLiteral("Open a run directory (trajectory or single pose)"));
    if (dir.isEmpty())
        return;
    // Launch a fresh reader process on the chosen directory rather than
    // tearing down and rebuilding the scene in place — multiple-instance
    // safe (project_viewer_cli_needs). main() sniffs the path shape.
    const bool ok = QProcess::startDetached(QCoreApplication::applicationFilePath(),
                                            QStringList{dir});
    if (!ok)
        qCWarning(cWindow).noquote() << "failed to launch a reader for" << dir;
}

void ReaderMainWindow::onOpenSignalDisplays() {
    ASSERT_THREAD(this);
    openSignalDisplayPicker();
}

}  // namespace h5reader::app
