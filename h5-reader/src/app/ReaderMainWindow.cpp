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
#include "QtOccupancyShellsOverlay.h"
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
#include <QMessageBox>
#include <QRegion>
#include <QSettings>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QStackedLayout>
#include <QStringList>
#include <QStatusBar>
#include <QStyle>
#include <QTimer>
#include <QToolBar>
#include <QToolButton>
#include <QUuid>
#include <QVariant>
#include <QWidget>

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

ReaderMainWindow::ReaderMainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("ReaderMainWindow"));

    qCInfo(cWindow).noquote() << "ctor entered";

    buildUi();
    buildToolbar();
    buildStatusBar();
    buildDocks();

    // Default size — wide enough for the playback + camera + transform +
    // metrics + overlays + panel controls to fit in one
    // toolbar without Qt's overflow chevron. QSettings restore overrides
    // this on later launches.
    resize(1600, 900);
    setWindowTitle(QStringLiteral("h5-reader"));

    // QSettings restore — geometry, dock state, log mask, recent menu.
    // Tolerant: missing or version-mismatched blobs leave the ctor's
    // explicit defaults intact. Runs AFTER all docks/toolbars exist so
    // restoreState has named docks to bind to.
    restoreAllSettings();
    setEmptyState();

    qCInfo(cWindow).noquote() << "ctor done";
}

bool ReaderMainWindow::loadRunPath(const QString& path) {
    ASSERT_THREAD(this);
    lastLoadError_.clear();
    if (path.isEmpty()) {
        lastLoadError_ = QStringLiteral("No calcset path was provided.");
        return false;
    }

    qCInfo(cWindow).noquote() << "loading" << path;
    auto loaded = h5reader::io::QtProteinLoader::LoadRunPath(path);
    if (!loaded.ok) {
        lastLoadError_ = loaded.error.isEmpty()
            ? QStringLiteral("Load failed for %1").arg(path)
            : loaded.error;
        qCCritical(cWindow).noquote() << "Load failed:" << lastLoadError_;
        return false;
    }
    if (loaded.decodeWarnings > 0) {
        qCWarning(cWindow).noquote()
            << "Decode completed with" << loaded.decodeWarnings << "warnings";
    }

    installLoadedRun(std::move(loaded));
    lastLoadError_.clear();
    return true;
}

void ReaderMainWindow::installLoadedRun(h5reader::io::QtLoadResult&& loaded) {
    ASSERT_THREAD(this);
    clearLoadedRun();
    loaded_ = std::make_unique<h5reader::io::QtLoadResult>(std::move(loaded));

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
    // changes are visible immediately.
    scene_ = new MoleculeScene(vtkWidget_, renderWindow_, this);
    scene_->Build(*loaded_->protein, *transformed_);
    applyOverlayActionState();
    scene_->refreshCurrentFrame();
    scene_->ResetCamera();
    ACONNECT(scene_, &MoleculeScene::cameraPlaneLockChanged,
             this, [this](bool) { updateCameraModeActions(); });
    if (scene_->cameraComposer()) {
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

    if (frameSlider_) {
        const QSignalBlocker block(frameSlider_);
        frameSlider_->setRange(0, std::max(0, T - 1));
        frameSlider_->setValue(0);
    }
    if (fpsSpinner_) {
        const QSignalBlocker block(fpsSpinner_);
        fpsSpinner_->setRange(1, 60);
        fpsSpinner_->setValue(playback_->fps());
    }

    // Atom picker — event filter on the VTK widget. Emits atomPicked(idx,
    // modifiers) on double-click. It stays dumb; AtomSelection interprets
    // the gesture and fans typed changes out.
    auto* firstRenderer = renderWindow_->GetRenderers()->GetFirstRenderer();
    picker_ = new QtAtomPicker(vtkWidget_, firstRenderer,
                                loaded_->protein.get(),
                                transformed_,
                                playback_, this);

    // Camera input filter — installed AFTER the picker so Qt's filter
    // chain runs THIS first. Double-click events fall through to the picker.
    cameraInputFilter_ = new CameraInputFilter(vtkWidget_, scene_,
                                                 scene_->cameraComposer(), this);

    inspectorDock_->setContext(loaded_->protein.get(), transformed_);
    ACONNECT(playback_,  &QtPlaybackController::frameChanged,
             inspectorDock_, &QtAtomInspectorDock::setFrame);

    // ---- Selection model — the single source of selection truth ----------
    selection_ = new model::AtomSelection(loaded_->protein.get(), this);

    signalCatalog_ = new model::TrajectorySignalCatalog(this);
    fieldAvailability_ = std::make_shared<model::TrajectoryFieldAvailability>(
        model::TrajectoryFieldAvailability::Build(loaded_->conformation.get(),
                                                  signalCatalog_->allDescriptorList()));
    signalCatalog_->setFieldAvailability(fieldAvailability_);
    visualizationContext_ = {};
    visualizationContext_.availability = fieldAvailability_.get();
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
    inspectorDock_->setFieldAvailability(fieldAvailability_);
    dashboardSignals_ = new model::DashboardSignalModel(this);
    dashboardSignals_->setFieldAvailability(fieldAvailability_);
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
    // tagging Picker afterward overrides the source.
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

    if (auto* meas = scene_->measurementOverlay()) {
        meas->setSelection(selection_);
        ACONNECT(selection_, &model::AtomSelection::changed,
                 meas,       &MeasurementOverlay::onSelectionChanged);
    }

    // Occupancy-shells overlay: focus-driven (NOT changed() — focus-only, so a
    // plain pick doesn't rebuild twice) + transform-driven (a fit-mode change
    // moves every aligned position, so the whole-trajectory aggregate is
    // stale). The overlay self-guards single-pose / rigid atoms. A render is
    // requested after focus/clear rebuilds (coalesced); transformChanged is
    // already followed by refreshCurrentFrame above.
    if (auto* occ = scene_->occupancyShellsOverlay()) {
        occ->setSelection(selection_);
        ACONNECT(selection_, &model::AtomSelection::focusChanged,
                 occ,        &QtOccupancyShellsOverlay::onFocusChanged);
        ACONNECT(selection_, &model::AtomSelection::cleared,
                 occ,        &QtOccupancyShellsOverlay::onSelectionCleared);
        ACONNECT(transformed_, &model::TransformedConformation::transformChanged,
                 occ,          &QtOccupancyShellsOverlay::onTransformChanged);
        ACONNECT(selection_, &model::AtomSelection::focusChanged, this,
                 [this](std::size_t) {
                     if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
                 });
        ACONNECT(selection_, &model::AtomSelection::cleared, this,
                 [this]() {
                     if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
                 });
    }
    ACONNECT(selection_, &model::AtomSelection::changed,
             this, [this]() {
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

    selectionDock_->setModel(selection_);

    dashboardStripDock_->setContext(loaded_->protein.get(), transformed_);
    dashboardStripDock_->setSignalModels(signalCatalog_, dashboardSignals_);
    dashboardStripDock_->setPanelModel(dashboardPanels_);
    dashboardStripDock_->setSelectionController(dashboardSelectionController_.data());
    dashboardStripDock_->setSelection(selection_);
    dashboardStripDock_->setTimeViewport(timeViewport_);
    dashboardController_ = dashboardStripDock_->displayController();
    if (dashboardController_)
        dashboardController_->setVisualizationContext(visualizationContext_);

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
    ACONNECT(playback_,           &QtPlaybackController::frameChanged,
             dashboardStripDock_, &DashboardStripDock::setFrame);

    // L-3a (2026-05-29): expose the scene's reveal overlay to the dashboard
    // controller so static.tensor mode can fire an ellipsoid glyph in the 3-D view.
    if (scene_->revealOverlay()) {
        dashboardStripDock_->setSceneOverlay(scene_->revealOverlay());
        visualizationContext_.hasSceneOverlay = true;
        if (dashboardController_)
            dashboardController_->setVisualizationContext(visualizationContext_);
        if (signalDisplayDialog_)
            signalDisplayDialog_->setVisualizationContext(visualizationContext_);
    }

    // DFT shielding campaign (optional): make the frame-local source
    // available to descriptor-family samplers. The `.LGS` carries the
    // typed `dft.frames[]` map — frame_index → meta.json — so the store
    // builds straight from it.
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

    QString alternatePath;
    if (loaded_->manifest.kind == h5reader::io::CalcsetManifest::Kind::MutantPair
        && loaded_->manifest.mutant_pair.has_value()) {
        alternatePath = loaded_->manifest.mutant_pair->ala_lgs_abspath;
    }
    updateMutantAlternateAction(alternatePath);

    setLoadedControlsEnabled(true);
    updateMetricAction();
    updateCameraModeActions();
    onFrameChanged(0);
    resetDashboardStateForRunLoad();
    if (centralContainer_) {
        if (auto* stack = qobject_cast<QStackedLayout*>(centralContainer_->layout()))
            stack->setCurrentWidget(vtkWidget_);
    }
    if (emptyPlaceholder_)
        emptyPlaceholder_->setVisible(false);
    setWindowTitle(QStringLiteral("h5-reader — %1").arg(loaded_->proteinId));
    if (!loaded_->runPath.isEmpty())
        addToRecentFiles(QDir(loaded_->runPath).absolutePath());
    syncRestServerContext();

    qCInfo(cWindow).noquote()
        << "run loaded | protein=" << loaded_->proteinId
        << "| path=" << loaded_->runPath;
}

void ReaderMainWindow::clearLoadedRun() {
    ASSERT_THREAD(this);

    if (playback_)
        playback_->pause();

    if (signalDisplayDialog_) {
        signalDisplayDialog_->hide();
        delete signalDisplayDialog_;
        signalDisplayDialog_ = nullptr;
    }

    if (dashboardStripDock_) {
        dashboardStripDock_->setSceneOverlay(nullptr);
        dashboardStripDock_->setSelection(nullptr);
        dashboardStripDock_->setSelectionController(nullptr);
        dashboardStripDock_->setTimeViewport(nullptr);
        dashboardStripDock_->setDftStore(nullptr);
        dashboardStripDock_->setPanelModel(nullptr);
        dashboardStripDock_->setSignalModels(nullptr, nullptr);
        dashboardStripDock_->setContext(nullptr, nullptr);
    }
    if (dashboardController_)
        dashboardController_->setVisualizationContext({});
    if (selectionDock_)
        selectionDock_->setModel(nullptr);
    if (inspectorDock_) {
        inspectorDock_->setContext(nullptr, nullptr);
        inspectorDock_->setFieldAvailability({});
        inspectorDock_->clearSelection();
    }

    delete cameraInputFilter_;
    cameraInputFilter_ = nullptr;
    delete picker_;
    picker_ = nullptr;
    delete scene_;
    scene_ = nullptr;

    delete playback_;
    playback_ = nullptr;
    delete timeViewport_;
    timeViewport_ = nullptr;

    delete dftStore_;
    dftStore_ = nullptr;
    delete dashboardSelectionController_;
    dashboardSelectionController_ = nullptr;
    delete dashboardPanels_;
    dashboardPanels_ = nullptr;
    delete dashboardSignals_;
    dashboardSignals_ = nullptr;
    delete signalCatalog_;
    signalCatalog_ = nullptr;
    delete selection_;
    selection_ = nullptr;

    delete transformed_;
    transformed_ = nullptr;

    visualizationContext_ = {};
    fieldAvailability_.reset();
    lastDashboardSelectedCount_ = 0;
    loaded_.reset();
    updateMutantAlternateAction({});

    setEmptyState();
    syncRestServerContext();
}

void ReaderMainWindow::setEmptyState() {
    ASSERT_THREAD(this);
    setLoadedControlsEnabled(false);
    updateFitModeLabel();
    updateCameraModeActions();
    if (emptyPlaceholder_)
        emptyPlaceholder_->setVisible(true);
    if (centralContainer_) {
        if (auto* stack = qobject_cast<QStackedLayout*>(centralContainer_->layout()))
            stack->setCurrentWidget(emptyPlaceholder_);
    }
    if (proteinLabel_)
        proteinLabel_->setText(QStringLiteral("No calcset loaded"));
    if (frameLabel_)
        frameLabel_->setText(QStringLiteral("frame —"));
    if (timeLabel_)
        timeLabel_->setText(QStringLiteral("t=— ps"));
    if (frameSlider_) {
        const QSignalBlocker block(frameSlider_);
        frameSlider_->setRange(0, 0);
        frameSlider_->setValue(0);
    }
    if (fpsSpinner_) {
        const QSignalBlocker block(fpsSpinner_);
        fpsSpinner_->setRange(1, 60);
        fpsSpinner_->setValue(5);
    }
    if (playAction_) {
        playAction_->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    }
    setWindowTitle(QStringLiteral("h5-reader"));
}

void ReaderMainWindow::setLoadedControlsEnabled(bool enabled) {
    ASSERT_THREAD(this);
    if (playAction_) playAction_->setEnabled(enabled);
    if (stepBackAction_) stepBackAction_->setEnabled(enabled);
    if (stepForwardAction_) stepForwardAction_->setEnabled(enabled);
    if (frameSlider_) frameSlider_->setEnabled(enabled);
    if (fpsSpinner_) fpsSpinner_->setEnabled(enabled);

    if (transformFitAction_) transformFitAction_->setEnabled(enabled && transformed_);
    if (instrumentAction_) instrumentAction_->setEnabled(enabled && scene_);
    if (signalDisplaysAction_)
        signalDisplaysAction_->setEnabled(enabled && selection_ && selection_->hasFocus());

    if (showRibbonAction_) showRibbonAction_->setEnabled(enabled && scene_);
    if (showRingsAction_) showRingsAction_->setEnabled(enabled && scene_);
    if (showButterflyAction_) showButterflyAction_->setEnabled(enabled && scene_);
    if (showBFieldAction_) showBFieldAction_->setEnabled(enabled && scene_);
    if (showOccupancyAction_) showOccupancyAction_->setEnabled(enabled && scene_);
}

void ReaderMainWindow::applyOverlayActionState() {
    ASSERT_THREAD(this);
    if (!scene_)
        return;
    if (showRibbonAction_ && scene_->ribbonOverlay())
        scene_->ribbonOverlay()->setVisible(showRibbonAction_->isChecked());
    if (showRingsAction_ && scene_->ringPolygonOverlay())
        scene_->ringPolygonOverlay()->setVisible(showRingsAction_->isChecked());
    if (showButterflyAction_ && scene_->fieldGridOverlay())
        scene_->fieldGridOverlay()->setVisible(showButterflyAction_->isChecked());
    if (showBFieldAction_ && scene_->bfieldStreamOverlay())
        scene_->bfieldStreamOverlay()->setVisible(showBFieldAction_->isChecked());
    if (showOccupancyAction_ && scene_->occupancyShellsOverlay())
        scene_->occupancyShellsOverlay()->setVisible(showOccupancyAction_->isChecked());
}

void ReaderMainWindow::updateMutantAlternateAction(const QString& alternatePath) {
    ASSERT_THREAD(this);
    if (mutantAlternateAction_) {
        if (fileMenu_)
            fileMenu_->removeAction(mutantAlternateAction_);
        delete mutantAlternateAction_.data();
        mutantAlternateAction_ = nullptr;
    }
    if (alternatePath.isEmpty() || !fileMenu_)
        return;

    mutantAlternateAction_ = fileMenu_->addAction(
        QStringLiteral("Open mutant alternate (ALA)…"));
    mutantAlternateAction_->setToolTip(QStringLiteral(
        "This run is a mutant pair; WT is opened in this window. "
        "Click to load the ALA pose in this window: %1").arg(alternatePath));
    ACONNECT(mutantAlternateAction_.data(), &QAction::triggered, this, [this, alternatePath]() {
        if (!loadRunPath(alternatePath)) {
            QMessageBox::critical(this,
                                  QStringLiteral("Open calcset failed"),
                                  lastLoadError());
        }
    });
}

void ReaderMainWindow::syncRestServerContext() {
    ASSERT_THREAD(this);
    if (!restServer_)
        return;
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
    centralContainer_ = new QWidget(this);
    auto* stack = new QStackedLayout(centralContainer_);
    stack->setContentsMargins(0, 0, 0, 0);
    stack->setStackingMode(QStackedLayout::StackAll);

    vtkWidget_    = new QVTKOpenGLNativeWidget(centralContainer_);
    renderWindow_ = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkWidget_->setRenderWindow(renderWindow_);
    stack->addWidget(vtkWidget_);

    emptyPlaceholder_ = new QLabel(QStringLiteral("Open a calcset (.LGS) to begin."), centralContainer_);
    emptyPlaceholder_->setAlignment(Qt::AlignCenter);
    emptyPlaceholder_->setWordWrap(true);
    emptyPlaceholder_->setAttribute(Qt::WA_TransparentForMouseEvents);
    emptyPlaceholder_->setStyleSheet(QStringLiteral(
        "QLabel { color: #5f6872; background: #fafafa; font-size: 18px; }"));
    stack->addWidget(emptyPlaceholder_);
    stack->setCurrentWidget(emptyPlaceholder_);
    setCentralWidget(centralContainer_);

    // File ▸ Open… loads a calcset into this window.
    fileMenu_ = menuBar()->addMenu(QStringLiteral("&File"));
    auto* openFileAct = fileMenu_->addAction(QStringLiteral("Open…"));
    openFileAct->setShortcut(QKeySequence::Open);  // Ctrl+O — pick a .LGS file with the mouse
    ACONNECT(openFileAct, &QAction::triggered, this, &ReaderMainWindow::onOpenFile);

    auto* openDirAct = fileMenu_->addAction(QStringLiteral("Open Directory…"));
    ACONNECT(openDirAct, &QAction::triggered, this, &ReaderMainWindow::onOpenDirectory);

    // File ▸ Recent — populated from QSettings during restoreAllSettings.
    // Empty until then; each entry loads into this window on click.
    recentMenu_ = fileMenu_->addMenu(QStringLiteral("&Recent"));
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

    stepBackAction_ = tb->addAction(
        style()->standardIcon(QStyle::SP_MediaSeekBackward),
        QStringLiteral("Step back"));
    ACONNECT(stepBackAction_.data(), &QAction::triggered,
             this, [this]() {
                 if (playback_) playback_->stepBackward();
             });
    stepForwardAction_  = tb->addAction(
        style()->standardIcon(QStyle::SP_MediaSeekForward),
        QStringLiteral("Step forward"));
    ACONNECT(stepForwardAction_.data(), &QAction::triggered,
             this, [this]() {
                 if (playback_) playback_->stepForward();
             });

    tb->addSeparator();

    frameSlider_ = new QSlider(Qt::Horizontal, tb);
    frameSlider_->setMinimumWidth(400);
    tb->addWidget(frameSlider_);
    ACONNECT(frameSlider_.data(), &QSlider::valueChanged,
             this, [this](int frame) {
                 if (playback_) playback_->setFrame(frame);
             });

    tb->addSeparator();
    tb->addWidget(new QLabel(QStringLiteral("fps"), tb));
    fpsSpinner_ = new QSpinBox(tb);
    fpsSpinner_->setSuffix(QStringLiteral(" /s"));
    tb->addWidget(fpsSpinner_);
    ACONNECT(fpsSpinner_.data(), qOverload<int>(&QSpinBox::valueChanged),
             this, [this](int fps) {
                 if (playback_) playback_->setFps(fps);
             });

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

    // Overlay toggles — inert until a scene has been loaded.
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

    showOccupancyAction_ = tb->addAction(QStringLiteral("Shadow"));
    showOccupancyAction_->setCheckable(true);
    showOccupancyAction_->setChecked(false);   // off by default
    showOccupancyAction_->setToolTip(QStringLiteral(
        "Occupation-probability envelope shells for the FOCUSED atom: nested "
        "50% / 90% highest-density regions over the trajectory (backbone-aligned). "
        "Trajectory data only; rigid atoms are skipped."));

    ACONNECT(showRibbonAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->ribbonOverlay()) return;
                 scene_->ribbonOverlay()->setVisible(on);
                 scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showRingsAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->ringPolygonOverlay()) return;
                 scene_->ringPolygonOverlay()->setVisible(on);
                 scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showButterflyAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->fieldGridOverlay()) return;
                 scene_->fieldGridOverlay()->setVisible(on);
                 if (on) scene_->refreshCurrentFrame();
                 else    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showBFieldAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->bfieldStreamOverlay()) return;
                 scene_->bfieldStreamOverlay()->setVisible(on);
                 if (on) scene_->refreshCurrentFrame();
                 else    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showOccupancyAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->occupancyShellsOverlay()) return;
                 // setVisible(true) rebuilds for the current focus; shells are
                 // frame-invariant so a plain render (not refreshCurrentFrame)
                 // suffices.
                 scene_->occupancyShellsOverlay()->setVisible(on);
                 scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
}

void ReaderMainWindow::buildStatusBar() {
    proteinLabel_ = new QLabel(QStringLiteral("No calcset loaded"), this);
    frameLabel_   = new QLabel(QStringLiteral("frame —"), this);
    timeLabel_    = new QLabel(QStringLiteral("t=— ps"), this);

    statusBar()->addPermanentWidget(proteinLabel_);
    statusBar()->addPermanentWidget(frameLabel_);
    statusBar()->addPermanentWidget(timeLabel_);
}

void ReaderMainWindow::buildDocks() {
    // Atom Info dock — tabified on the LEFT alongside Selection + Strip.
    // It is constructed before load and starts with its own placeholder.
    inspectorDock_ = new QtAtomInspectorDock(this);
    addDockWidget(Qt::LeftDockWidgetArea, inspectorDock_);

    // Selected-atoms panel — the QListView binds to AtomSelection after load.
    selectionDock_ = new SelectionDock(this);
    selectionDock_->setModel(nullptr);
    addDockWidget(Qt::LeftDockWidgetArea, selectionDock_);
    tabifyDockWidget(inspectorDock_, selectionDock_);

    // Dashboard strips — the dock and controller are stable chrome; run-backed
    // models are swapped in by installLoadedRun().
    dashboardStripDock_ = new DashboardStripDock(this);
    dashboardStripDock_->setContext(nullptr, nullptr);
    dashboardStripDock_->setSignalModels(nullptr, nullptr);
    dashboardStripDock_->setPanelModel(nullptr);
    dashboardStripDock_->setSelectionController(nullptr);
    dashboardStripDock_->setSelection(nullptr);
    dashboardStripDock_->setTimeViewport(nullptr);
    dashboardController_ = dashboardStripDock_->displayController();
    if (dashboardController_)
        dashboardController_->setVisualizationContext({});
    addDockWidget(Qt::LeftDockWidgetArea, dashboardStripDock_);
    tabifyDockWidget(inspectorDock_, dashboardStripDock_);
    inspectorDock_->raise();
    resizeDocks({inspectorDock_}, {360}, Qt::Horizontal);

    // Start clean — no property docks open on launch. The user opens each via
    // View -> Panels or the toolbar menu. QSettings restore can override this
    // for users who intentionally left a dock visible.
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
    ACONNECT(dashboardStripDock_, &DashboardStripDock::revealRequested,
             this, [this](const model::SignalBinding& binding) {
                 if (scene_)
                     scene_->revealBinding(binding);
             });
    ACONNECT(dashboardStripDock_, &DashboardStripDock::metricPickerRequested,
             this, &ReaderMainWindow::onOpenSignalDisplays);

    // Panel recovery — one QAction per dock from QDockWidget::toggleViewAction().
    // The same QAction instances live in View -> Panels and the toolbar menu.
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
}

void ReaderMainWindow::onFrameChanged(int t) {
    ASSERT_THREAD(this);
    if (!loaded_ || !loaded_->conformation) {
        if (frameLabel_)
            frameLabel_->setText(QStringLiteral("frame —"));
        if (timeLabel_)
            timeLabel_->setText(QStringLiteral("t=— ps"));
        if (frameSlider_ && frameSlider_->value() != 0) {
            const QSignalBlocker block(frameSlider_);
            frameSlider_->setValue(0);
        }
        return;
    }

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
    const auto setOne = [](QAction* a, bool on) {
        if (!a) return;
        const QSignalBlocker block(a);
        a->setChecked(on);
    };
    if (!scene_ || !scene_->cameraComposer()) {
        setOne(focusAction_,     false);
        setOne(newmanAction_,    false);
        setOne(planeLockAction_, false);
        setOne(freeAction_,      true);
        return;
    }
    const auto kind = scene_->cameraComposer()->mode().kind;
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
    if (!transformFitAction_)
        return;
    if (!transformed_) {
        transformFitAction_->setText(QStringLiteral("Mode: no calcset"));
        transformFitAction_->setToolTip(fitModeToolTip());
        return;
    }

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
    if (!scene_ || !scene_->cameraComposer() || !selection_ || !selection_->hasFocus()
        || !loaded_ || !loaded_->protein) {
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
    if (!loadRunPath(path)) {
        QMessageBox::critical(this,
                              QStringLiteral("Open calcset failed"),
                              lastLoadError());
    }
}

void ReaderMainWindow::onPlayPauseClicked() {
    ASSERT_THREAD(this);
    if (playback_) playback_->togglePlayPause();
}

void ReaderMainWindow::onOpenFile() {
    ASSERT_THREAD(this);
    // Pick a .LGS calcset file directly with the mouse. CalcsetManifest::Load
    // (via ResolveLgsPath) accepts a .LGS file path, so load it into this
    // window.
    const QString file = QFileDialog::getOpenFileName(
        this, QStringLiteral("Open a calcset (.LGS) file"),
        qEnvironmentVariable("H5READER_OPEN_DIR"),
        QStringLiteral("Calcset manifest (*.lgs *.LGS);;All files (*)"));
    if (file.isEmpty())
        return;
    if (!loadRunPath(file)) {
        QMessageBox::critical(this,
                              QStringLiteral("Open calcset failed"),
                              lastLoadError());
    }
}

void ReaderMainWindow::onOpenDirectory() {
    ASSERT_THREAD(this);
    const QString dir = QFileDialog::getExistingDirectory(
        this, QStringLiteral("Open a run directory (trajectory or single pose)"));
    if (dir.isEmpty())
        return;
    if (!loadRunPath(dir)) {
        QMessageBox::critical(this,
                              QStringLiteral("Open calcset failed"),
                              lastLoadError());
    }
}

void ReaderMainWindow::onOpenSignalDisplays() {
    ASSERT_THREAD(this);
    openSignalDisplayPicker();
}

}  // namespace h5reader::app
