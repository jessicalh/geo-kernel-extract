#include "ReaderMainWindow.h"

#include "CameraAnchorHelper.h"
#include "CameraComposer.h"
#include "CameraInputFilter.h"
#include "CameraMode.h"
#include "OrientationPolicy.h"
#include "MoleculeScene.h"
#include "NearbySignalModel.h"
#include "NewmanDock.h"
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
#include "../model/CsaShape.h"
#include "../model/MolecularFrame.h"
#include "../model/MolecularFrameSelect.h"
#include "../model/ConformationGeometry.h"
#include "CsaTensorOverlay.h"

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
#include <QMouseEvent>
#include <QMessageBox>
#include <QRegion>
#include <QSettings>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QStackedLayout>
#include <QStringList>
#include <QStatusBar>
#include <QColor>
#include <QIcon>
#include <QPainter>
#include <QPalette>
#include <QPixmap>
#include <QPolygonF>
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
             this,      [this](bool) { refreshControlStates(); });

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

    // Click on empty space (no atom hit, no drag) stops/restarts the animation.
    ACONNECT(cameraInputFilter_, &CameraInputFilter::viewportClicked,
             this, [this](QPointF pos) {
                 if (!picker_ || !playback_) return;
                 const auto hit = picker_->atomAt(
                     static_cast<int>(pos.x()), static_cast<int>(pos.y()));
                 if (!hit) playback_->togglePlayPause();
             });

    inspectorDock_->setContext(loaded_->protein.get(), transformed_);
    ACONNECT(playback_,  &QtPlaybackController::frameChanged,
             inspectorDock_, &QtAtomInspectorDock::setFrame);

    if (newmanDock_) {
        newmanDock_->setContext(loaded_->protein.get(), loaded_->conformation.get());
        ACONNECT(playback_, &QtPlaybackController::frameChanged,
                 newmanDock_, &NewmanDock::setFrame);
    }

    // ---- Selection model — the single source of selection truth ----------
    selection_ = new model::AtomSelection(loaded_->protein.get(), this);

    // Nearby-residue source for the Filter checklist. Prefer the transformed
    // (displayed) conformation; fall back to the base conformation if no
    // backbone fit is active so the checklist still populates. Lazily created
    // once and re-pointed on each load (parented to this; never per-load leak).
    if (!filterNearby_)
        filterNearby_ = new NearbySignalModel(this);
    filterNearby_->setContext(
        loaded_->protein.get(),
        transformed_ ? static_cast<model::Conformation*>(transformed_)
                     : loaded_->conformation.get());
    filterNearby_->setRadiusAngstrom(5.0);
    filterResidues_.clear();

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

    // Reveal-on-pick: picking an atom brings up its Inspector (it starts hidden;
    // this is the dock's reveal path now that the Panels menu is gone).
    ACONNECT(picker_, &QtAtomPicker::atomPicked,
             this,   [this](std::size_t, Qt::KeyboardModifiers) {
                 if (inspectorDock_ && !inspectorDock_->isVisible())
                     revealDockQueued(inspectorDock_);
             });

    ACONNECT(selection_, &model::AtomSelection::focusChanged,
             inspectorDock_, &QtAtomInspectorDock::setPickedAtom);
    ACONNECT(selection_, &model::AtomSelection::cleared,
             inspectorDock_, &QtAtomInspectorDock::clearSelection);

    if (newmanDock_) {
        ACONNECT(selection_, &model::AtomSelection::focusChanged,
                 newmanDock_, &NewmanDock::setFocusAtom);
        ACONNECT(selection_, &model::AtomSelection::cleared,
                 newmanDock_, &NewmanDock::clear);
        // Reveal on first focus (3D pick or REST), like the Inspector.
        ACONNECT(selection_, &model::AtomSelection::focusChanged, this,
                 [this](std::size_t) {
                     if (newmanDock_ && !newmanDock_->isVisible())
                         revealDockQueued(newmanDock_);
                 });
    }
    ACONNECT(selection_, &model::AtomSelection::focusChanged, this,
             [this](std::size_t) { refreshControlStates(); });
    ACONNECT(selection_, &model::AtomSelection::cleared, this,
             [this]() { refreshControlStates(); });

    // CSA tensor glyph (mode-2): focus + frame driven; honest gap on a missing
    // DFT frame; raw->display alignment via the molecular frame.
    ACONNECT(selection_, &model::AtomSelection::focusChanged, this,
             [this](std::size_t) { updateCsaGlyph(); });
    ACONNECT(selection_, &model::AtomSelection::cleared, this, [this]() {
        if (scene_ && scene_->csaOverlay()) {
            scene_->csaOverlay()->clear();
            scene_->requestRender(MoleculeScene::RenderSource::Overlay);
        }
    });
    ACONNECT(playback_, &QtPlaybackController::frameChanged, this,
             [this](int) { updateCsaGlyph(); });

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
                 refreshControlStates();
                 if (scene_) scene_->refreshCurrentFrame();
             });

    // Selection summary in the status bar (count + measurement kind).
    ACONNECT(selection_, &model::AtomSelection::changed, this, [this]() { updateSelectionStatus(); });
    ACONNECT(selection_, &model::AtomSelection::cleared, this, [this]() { updateSelectionStatus(); });
    updateSelectionStatus();

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

    refreshControlStates();
    onFrameChanged(0);
    resetDashboardStateForRunLoad();
    if (centralContainer_) {
        if (auto* stack = qobject_cast<QStackedLayout*>(centralContainer_->layout()))
            stack->setCurrentWidget(vtkWidget_);
    }
    if (emptyPlaceholder_)
        emptyPlaceholder_->setVisible(false);
    setWindowTitle(QStringLiteral("h5-reader — %1").arg(loaded_->proteinId));
    if (proteinLabel_)
        proteinLabel_->setText(loaded_->proteinId);
    if (!loaded_->runPath.isEmpty())
        addToRecentFiles(QDir(loaded_->runPath).absolutePath());
    syncRestServerContext();

    qCInfo(cWindow).noquote()
        << "run loaded | protein=" << loaded_->proteinId
        << "| path=" << loaded_->runPath;
}

void ReaderMainWindow::updateCsaGlyph() {
    ASSERT_THREAD(this);
    CsaTensorOverlay* overlay = scene_ ? scene_->csaOverlay() : nullptr;
    if (!overlay)
        return;
    auto redraw = [this] {
        if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    };
    auto hide = [&] { overlay->clear(); redraw(); };

    if (!loaded_ || !dftStore_ || !transformed_ || !selection_->hasFocus()) {
        hide();
        return;
    }
    const model::QtProtein* protein = loaded_->protein.get();
    model::Conformation* rawConf = loaded_->conformation.get();
    if (!protein || !rawConf) {
        hide();
        return;
    }
    const std::size_t atom = selection_->focus();
    if (atom >= protein->atomCount()) {
        hide();
        return;
    }
    const int frameI = playback_ ? playback_->currentFrame() : 0;
    const std::size_t frame = static_cast<std::size_t>(frameI < 0 ? 0 : frameI);

    // DFT shielding tensor for this atom + frame; the campaign is partial, so a
    // missing frame is an honest gap (clear), never a faked value.
    const std::size_t original = rawConf->originalFrameIndex(frame);
    if (!dftStore_->hasJob(original)) {
        hide();
        return;
    }
    dftStore_->requestFrame(original);
    const model::DftShieldingFrame* dft = dftStore_->frame(original);
    if (!dft || !dft->valid || atom >= dft->atoms.size()) {
        hide();
        return;
    }
    const model::Mat3 sigmaRaw = dft->atoms[atom].total_raw;

    // Align the raw-frame tensor onto the displayed (Kabsch-stabilized)
    // molecule via the molecular frame: sigma_mol = M_raw^T sigma M_raw is
    // rotation-invariant (and IS the molecular-frame component basis the stats
    // report), so lift it back out through the molecular frame built from
    // DISPLAY positions. Framed atoms only; unframed fall back to the raw PAS.
    model::CsaShape shape;
    std::optional<model::Mat3> molDisp;
    if (const auto spec = model::SelectMolecularFrameSpec(*protein, atom)) {
        const auto posRaw = [&](std::int32_t i) -> std::optional<model::Vec3> {
            if (i < 0 || static_cast<std::size_t>(i) >= protein->atomCount())
                return std::nullopt;
            return rawConf->atomPosition(frame, static_cast<std::size_t>(i));
        };
        const auto ringRaw =
            [&](std::int32_t r) -> std::optional<std::pair<model::Vec3, model::Vec3>> {
            if (r < 0) return std::nullopt;
            const model::RingGeometry g =
                model::RingGeometryAt(*rawConf, static_cast<std::size_t>(r), frame);
            return std::make_pair(g.center, g.normal);
        };
        const auto posDisp = [&](std::int32_t i) -> std::optional<model::Vec3> {
            if (i < 0 || static_cast<std::size_t>(i) >= protein->atomCount())
                return std::nullopt;
            return transformed_->atomPosition(frame, static_cast<std::size_t>(i));
        };
        const auto ringDisp =
            [&](std::int32_t r) -> std::optional<std::pair<model::Vec3, model::Vec3>> {
            if (r < 0) return std::nullopt;
            const model::RingGeometry g =
                model::RingGeometryAt(*transformed_, static_cast<std::size_t>(r), frame);
            return std::make_pair(g.center, g.normal);
        };
        model::MolFrameContinuity contRaw;
        model::MolFrameContinuity contDisp;
        const auto mRaw = model::BuildMolecularFrameAxes(*spec, posRaw, ringRaw, contRaw);
        const auto mDisp = model::BuildMolecularFrameAxes(*spec, posDisp, ringDisp, contDisp);
        if (mRaw && mDisp) {
            const model::Mat3 sigmaMol = mRaw->transpose() * sigmaRaw * (*mRaw);
            const model::Mat3 sigmaDisp = (*mDisp) * sigmaMol * mDisp->transpose();
            shape = model::ComputeCsaShape(sigmaDisp);
            if (shape.valid)
                molDisp = *mDisp;
        }
    }
    if (!shape.valid)
        shape = model::ComputeCsaShape(sigmaRaw);  // unframed: raw-frame PAS
    if (!shape.valid) {
        hide();
        return;
    }

    const model::Vec3 atomPos = transformed_->atomPosition(frame, atom);
    qCDebug(cWindow).noquote()
        << "CSA glyph | atom=" << atom << "| frame=" << frame
        << "| framed=" << molDisp.has_value() << "| iso=" << shape.sigma_iso
        << "| s11/s22/s33=" << shape.principal_values[0] << shape.principal_values[1]
        << shape.principal_values[2] << "| eta=" << shape.eta << "| span=" << shape.span;
    overlay->show(atomPos, shape, molDisp);
    redraw();
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
    if (selectionLabel_)
        selectionLabel_->setText(QStringLiteral("no selection"));
    if (inspectorDock_) {
        inspectorDock_->setContext(nullptr, nullptr);
        inspectorDock_->setFieldAvailability({});
        inspectorDock_->clearSelection();
    }
    if (newmanDock_)
        newmanDock_->setContext(nullptr, nullptr);

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
    refreshControlStates();
    updateFitModeLabel();
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
    setWindowTitle(QStringLiteral("h5-reader"));
}

void ReaderMainWindow::refreshControlStates() {
    ASSERT_THREAD(this);
    // Derive the whole operating state once — single source of truth.
    const bool loaded   = scene_ != nullptr;
    const int  frames   = (loaded && loaded_ && loaded_->conformation)
                            ? static_cast<int>(loaded_->conformation->frameCount()) : 0;
    const bool playable = loaded && frames > 1;            // single pose: nothing to play
    const bool playing  = playback_ && playback_->isPlaying();
    const bool traj     = loaded && visualizationContext_.hasTrajectory;
    const bool hasFocus = selection_ && selection_->hasFocus();
    const bool hasRings = loaded && loaded_ && loaded_->protein
                            && loaded_->protein->ringCount() > 0;
    // Filter (isolation) mode hides the whole-molecule overlays, so their
    // toggles can't do anything while it's on — keep them on the toolbar
    // (consistent layout) but disabled.
    const bool filtered = loaded && scene_->atomFilterActive();

    const auto en = [](QAction* a, bool e) { if (a) a->setEnabled(e); };

    // Transport — only meaningful for a multi-frame trajectory.
    en(playBackAction_,    playable);
    en(stepBackAction_,    playable);
    en(stepForwardAction_, playable);
    en(playForwardAction_, playable);
    en(stopAction_,        playable && playing);   // stop is inert when stopped
    if (frameSlider_) frameSlider_->setEnabled(playable);
    if (fpsSpinner_)  fpsSpinner_->setEnabled(playable);

    // Analysis controls.
    en(transformFitAction_,   loaded && transformed_ != nullptr);
    en(signalDisplaysAction_, loaded && hasFocus);
    en(instrumentAction_,     loaded);

    // Overlays — gated on the data that makes each one mean something.
    en(showRibbonAction_,    loaded   && !filtered);
    en(showRingsAction_,     hasRings && !filtered);
    en(showButterflyAction_, hasRings && !filtered);
    en(showBFieldAction_,    hasRings && !filtered);
    en(showOccupancyAction_, traj     && !filtered);   // per-atom; needs trajectory

    // Camera modes (enable gating + checked-state sync).
    updateCameraModeActions();

    // Filter button: usable whenever a run is loaded (so an active filter can
    // always be cleared) — its dropdown content adapts to the current focus.
    if (filterButton_) filterButton_->setEnabled(loaded);
    updateFilterButton();
}

QJsonObject ReaderMainWindow::uiStateJson() const {
    const bool loaded  = scene_ != nullptr;
    const int  frames  = (loaded && loaded_ && loaded_->conformation)
                           ? static_cast<int>(loaded_->conformation->frameCount()) : 0;
    const bool playing = playback_ && playback_->isPlaying();

    const auto ctl = [](const QAction* a) {
        QJsonObject o;
        o[QStringLiteral("present")] = (a != nullptr);
        o[QStringLiteral("enabled")] = (a != nullptr && a->isEnabled());
        if (a && a->isCheckable())
            o[QStringLiteral("checked")] = a->isChecked();
        return o;
    };

    QJsonObject controls;
    controls[QStringLiteral("playBack")]     = ctl(playBackAction_);
    controls[QStringLiteral("stepBack")]     = ctl(stepBackAction_);
    controls[QStringLiteral("stop")]         = ctl(stopAction_);
    controls[QStringLiteral("stepForward")]  = ctl(stepForwardAction_);
    controls[QStringLiteral("playForward")]  = ctl(playForwardAction_);
    controls[QStringLiteral("focus")]        = ctl(focusAction_);
    controls[QStringLiteral("transformFit")] = ctl(transformFitAction_);
    controls[QStringLiteral("metrics")]      = ctl(signalDisplaysAction_);
    controls[QStringLiteral("ribbon")]       = ctl(showRibbonAction_);
    controls[QStringLiteral("rings")]        = ctl(showRingsAction_);
    controls[QStringLiteral("butterfly")]    = ctl(showButterflyAction_);
    controls[QStringLiteral("bfield")]       = ctl(showBFieldAction_);
    controls[QStringLiteral("occupancy")]    = ctl(showOccupancyAction_);

    // Filter is a QToolButton (not a QAction), so report it explicitly:
    // present/enabled like the others, plus whether isolation is active and
    // how many residues are pinned. Keeps /ui/state an honest mirror of the UI.
    QJsonObject filter;
    filter[QStringLiteral("present")] = (filterButton_ != nullptr);
    filter[QStringLiteral("enabled")] = (filterButton_ && filterButton_->isEnabled());
    filter[QStringLiteral("active")]  = (scene_ && scene_->atomFilterActive());
    filter[QStringLiteral("count")]   = static_cast<int>(filterResidues_.size());
    controls[QStringLiteral("filter")] = filter;

    QJsonObject sel;
    sel[QStringLiteral("count")] = static_cast<int>(selection_ ? selection_->count() : 0);
    sel[QStringLiteral("focus")] = (selection_ && selection_->hasFocus());

    QJsonObject out;
    out[QStringLiteral("loaded")]        = loaded;
    out[QStringLiteral("protein")]       = (loaded && loaded_) ? loaded_->proteinId : QString();
    out[QStringLiteral("frames")]        = frames;
    out[QStringLiteral("currentFrame")]  = playback_ ? playback_->currentFrame() : 0;
    out[QStringLiteral("playing")]       = playing;
    out[QStringLiteral("playDirection")] = playback_ ? playback_->direction() : 1;
    out[QStringLiteral("selection")]     = sel;
    out[QStringLiteral("cameraMode")]    =
        (loaded && scene_ && scene_->cameraComposer())
            ? QString::fromLatin1(NameFor(scene_->cameraComposer()->mode().kind))
            : QStringLiteral("none");
    out[QStringLiteral("controls")]      = controls;
    return out;
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

bool ReaderMainWindow::setOverlayVisible(const QString& name, bool on) {
    ASSERT_THREAD(this);
    // Map a stable automation key → the toolbar QAction, then setChecked()
    // so the already-connected QAction::toggled handler runs the real overlay
    // logic (setVisible + refreshCurrentFrame for the per-frame kernel
    // overlays). This keeps REST control and the toolbar UI on one code path.
    const QString key = name.toLower();
    QPointer<QAction> a;
    if (key == QStringLiteral("ribbon"))
        a = showRibbonAction_;
    else if (key == QStringLiteral("rings"))
        a = showRingsAction_;
    else if (key == QStringLiteral("butterfly") || key == QStringLiteral("fieldgrid")
             || key == QStringLiteral("field") || key == QStringLiteral("isosurface"))
        a = showButterflyAction_;
    else if (key == QStringLiteral("bfield") || key == QStringLiteral("streamlines")
             || key == QStringLiteral("stream"))
        a = showBFieldAction_;
    else if (key == QStringLiteral("shadow") || key == QStringLiteral("occupancy")
             || key == QStringLiteral("shells"))
        a = showOccupancyAction_;
    if (!a)
        return false;
    if (a->isChecked() != on)
        a->setChecked(on);   // fires QAction::toggled → overlay setVisible + refresh
    return true;
}

bool ReaderMainWindow::setFieldThreshold(double ppm) {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->fieldGridOverlay())
        return false;
    scene_->fieldGridOverlay()->setThresholdPpm(ppm);
    // Only the isovalue changed (scalar field unchanged), so a plain render
    // re-runs the contour filters — no need to re-evaluate the kernel grid.
    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    return true;
}

bool ReaderMainWindow::setFieldExtent(double sigmaA) {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->fieldGridOverlay())
        return false;
    // setGaussianExtent re-evaluates the scalar grid (the taper reshapes the
    // field); a render then re-contours it.
    scene_->fieldGridOverlay()->setGaussianExtent(sigmaA);
    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    return true;
}

bool ReaderMainWindow::setFieldPeak(double amplitude) {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->fieldGridOverlay())
        return false;
    scene_->fieldGridOverlay()->setGaussianPeak(amplitude);
    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    return true;
}

void ReaderMainWindow::setResidueFilter(const std::vector<std::size_t>& residues) {
    ASSERT_THREAD(this);
    if (!scene_ || !loaded_ || !loaded_->protein)
        return;

    // The whole-molecule overlays (ribbon, rings, field, occupancy) would keep
    // drawing the entire structure and defeat the isolation, so hide them while
    // filtered and restore them (per the toolbar toggles) when the filter clears.
    const auto restoreOverlays = [this]() { applyOverlayActionState(); };
    const auto hideOverlays = [this]() {
        if (scene_->ribbonOverlay())          scene_->ribbonOverlay()->setVisible(false);
        if (scene_->ringPolygonOverlay())     scene_->ringPolygonOverlay()->setVisible(false);
        if (scene_->fieldGridOverlay())       scene_->fieldGridOverlay()->setVisible(false);
        if (scene_->bfieldStreamOverlay())    scene_->bfieldStreamOverlay()->setVisible(false);
        if (scene_->occupancyShellsOverlay()) scene_->occupancyShellsOverlay()->setVisible(false);
    };

    const auto* protein = loaded_->protein.get();
    std::vector<std::size_t> atoms;
    if (!residues.empty()) {
        std::vector<char> keep(protein->residueCount(), 0);
        for (std::size_t r : residues)
            if (r < keep.size()) keep[r] = 1;
        for (std::size_t i = 0; i < protein->atomCount(); ++i) {
            const auto& atom = protein->atom(i);
            if (atom.residueIndex >= 0
                && static_cast<std::size_t>(atom.residueIndex) < keep.size()
                && keep[static_cast<std::size_t>(atom.residueIndex)])
                atoms.push_back(i);
        }
    }

    if (atoms.empty()) {
        filterResidues_.clear();
        restoreOverlays();
        scene_->clearAtomFilter();
    } else {
        filterResidues_ = residues;   // single source of truth: keeps the button
        hideOverlays();               // label + checklist consistent whether the
        scene_->setAtomFilter(atoms); // filter was driven by the menu or REST
    }
    refreshControlStates();   // grey overlay toggles while filtered; restore on clear
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
            inspectorDock_, dashboardStripDock_
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

void ReaderMainWindow::updateSelectionStatus() {
    ASSERT_THREAD(this);
    if (!selectionLabel_)
        return;
    if (!selection_ || selection_->empty()) {
        selectionLabel_->setText(QStringLiteral("no selection"));
        return;
    }
    const std::size_t        n = selection_->count();
    const model::GeometryKind k = selection_->geometryKind();
    if (k == model::GeometryKind::None)
        selectionLabel_->setText(QStringLiteral("%1 atom selected").arg(n));
    else
        selectionLabel_->setText(QStringLiteral("%1 atoms · %2")
            .arg(n).arg(QString::fromLatin1(model::NameForGeometryKind(k))));
}

void ReaderMainWindow::buildFilterMenu() {
    ASSERT_THREAD(this);
    if (!filterMenu_)
        return;
    filterMenu_->clear();

    // Leaving filter mode is always offered, enabled only while a filter is on.
    QAction* showAll = filterMenu_->addAction(QStringLiteral("Show whole structure"));
    showAll->setEnabled(scene_ && scene_->atomFilterActive());
    ACONNECT(showAll, &QAction::triggered, this, [this]() {
        setResidueFilter({});      // empty restores the full structure + overlays
    });
    filterMenu_->addSeparator();

    // The checklist is the residues near the focused atom. No focus → nothing
    // to list; say so rather than present an empty menu.
    if (!selection_ || !selection_->hasFocus() || !filterNearby_) {
        QAction* hint = filterMenu_->addAction(
            QStringLiteral("Select an atom to list nearby residues"));
        hint->setEnabled(false);
        return;
    }

    // Rebuild the neighbourhood for the current focus + frame, then offer one
    // checkable row per nearby residue, nearest first.
    filterNearby_->setAnchor(selection_->focus(),
                             playback_ ? playback_->currentFrame() : 0);

    struct Row { std::size_t residue; double dist; QString label; };
    std::vector<Row> rows;
    const int rowN = filterNearby_->rowCount();
    for (int r = 0; r < rowN; ++r) {
        const NearbySignalModel::Candidate* c =
            filterNearby_->candidateAt(filterNearby_->index(r, 0));
        if (!c || c->kind != NearbySignalModel::CandidateKind::Residue
              || !c->residueContext)
            continue;
        rows.push_back({*c->residueContext, c->distanceAngstrom, c->label});
    }
    std::sort(rows.begin(), rows.end(),
              [](const Row& a, const Row& b) { return a.dist < b.dist; });

    if (rows.empty()) {
        QAction* none = filterMenu_->addAction(
            QStringLiteral("No residues within %1 Å")
                .arg(filterNearby_->radiusAngstrom(), 0, 'f', 1));
        none->setEnabled(false);
        return;
    }

    for (const Row& row : rows) {
        QAction* a = filterMenu_->addAction(
            QStringLiteral("%1 · %2 Å").arg(row.label).arg(row.dist, 0, 'f', 1));
        a->setCheckable(true);
        a->setChecked(std::find(filterResidues_.begin(), filterResidues_.end(),
                                row.residue) != filterResidues_.end());
        const std::size_t residue = row.residue;
        ACONNECT(a, &QAction::toggled, this, [this, residue](bool on) {
            onFilterResidueToggled(residue, on);
        });
    }
}

void ReaderMainWindow::onFilterResidueToggled(std::size_t residue, bool on) {
    ASSERT_THREAD(this);
    // Build the desired set locally; setResidueFilter is the sole writer of
    // filterResidues_, so it (not this handler) keeps the button label synced.
    std::vector<std::size_t> next = filterResidues_;
    const auto it = std::find(next.begin(), next.end(), residue);
    if (on) {
        if (it == next.end())
            next.push_back(residue);
    } else if (it != next.end()) {
        next.erase(it);
    }
    setResidueFilter(next);   // empty set restores the whole structure
}

void ReaderMainWindow::updateFilterButton() {
    ASSERT_THREAD(this);
    if (!filterButton_)
        return;
    const bool active = scene_ && scene_->atomFilterActive();
    filterButton_->setText(active
        ? QStringLiteral("Filter (%1)").arg(filterResidues_.size())
        : QStringLiteral("Filter"));
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
}

namespace {
// Transport-control glyphs drawn in a given colour so they stay legible on any
// palette — Qt's SP_Media* standard icons render near-black on the dark Fusion
// toolbar. Centred in a 32 px square.
enum class TransportGlyph { PlayForward, PlayBackward, StepForward, StepBackward, Stop };

QIcon makeTransportIcon(TransportGlyph kind, const QColor& color) {
    constexpr int s = 32;
    QPixmap pm(s, s);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setPen(Qt::NoPen);
    p.setBrush(color);
    const double m   = s * 0.30;             // margin
    const double x0  = m, y0 = m, x1 = s - m, y1 = s - m, yc = s * 0.5;
    const double bar = (x1 - x0) * 0.32;     // bar thickness for the step glyphs
    switch (kind) {
        case TransportGlyph::PlayForward: {
            QPolygonF t; t << QPointF(x0, y0) << QPointF(x1, yc) << QPointF(x0, y1);
            p.drawPolygon(t);
            break;
        }
        case TransportGlyph::PlayBackward: {
            QPolygonF t; t << QPointF(x1, y0) << QPointF(x0, yc) << QPointF(x1, y1);
            p.drawPolygon(t);
            break;
        }
        case TransportGlyph::StepForward: {
            QPolygonF t; t << QPointF(x0, y0) << QPointF(x1 - bar, yc) << QPointF(x0, y1);
            p.drawPolygon(t);
            p.drawRect(QRectF(x1 - bar, y0, bar, y1 - y0));
            break;
        }
        case TransportGlyph::StepBackward: {
            QPolygonF t; t << QPointF(x1, y0) << QPointF(x0 + bar, yc) << QPointF(x1, y1);
            p.drawPolygon(t);
            p.drawRect(QRectF(x0, y0, bar, y1 - y0));
            break;
        }
        case TransportGlyph::Stop: {
            p.drawRect(QRectF(x0, y0, x1 - x0, y1 - y0));
            break;
        }
    }
    p.end();
    return QIcon(pm);
}

// A QMenu that does NOT close when a checkable item is clicked, so the Filter
// checklist lets you tick several residues in one go.
class FilterMenu final : public QMenu {
public:
    using QMenu::QMenu;
protected:
    void mouseReleaseEvent(QMouseEvent* e) override {
        QAction* a = activeAction();
        if (a && a->isEnabled() && a->isCheckable()) {
            a->trigger();   // toggle + fire; keep the menu open
            return;
        }
        QMenu::mouseReleaseEvent(e);
    }
};
}  // namespace

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

    // Transport: ⏪ play-back · ◀ step-back · ■ stop · ▶ step-fwd · ⏩ play-fwd.
    // Custom glyphs in the palette text colour (legible on the dark toolbar).
    const QColor glyph = palette().color(QPalette::ButtonText);

    playBackAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::PlayBackward, glyph),
        QStringLiteral("Play backward"));
    playBackAction_->setToolTip(QStringLiteral("Play continuously, backward in time."));
    ACONNECT(playBackAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->playBackward(); });

    stepBackAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::StepBackward, glyph),
        QStringLiteral("Step back"));
    stepBackAction_->setToolTip(QStringLiteral("Step one frame back."));
    ACONNECT(stepBackAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->stepBackward(); });

    stopAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::Stop, glyph),
        QStringLiteral("Stop"));
    stopAction_->setToolTip(QStringLiteral("Stop the animation."));
    ACONNECT(stopAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->pause(); });

    stepForwardAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::StepForward, glyph),
        QStringLiteral("Step forward"));
    stepForwardAction_->setToolTip(QStringLiteral("Step one frame forward."));
    ACONNECT(stepForwardAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->stepForward(); });

    playForwardAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::PlayForward, glyph),
        QStringLiteral("Play forward"));
    playForwardAction_->setToolTip(QStringLiteral("Play continuously, forward in time."));
    ACONNECT(playForwardAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->playForward(); });

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

    // Camera-mode cluster removed: Focus is now a single de-emphasized toggle
    // at the toolbar tail (built below). The Newman / Plane-lock / Free buttons
    // are gone; the composer keeps every mode for the dashboard-reveal + REST
    // paths, which are their real consumers.
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

    // Display isolation ("Filter"): a dropdown checklist of residues near the
    // focused atom. Ticking residues enters filter mode (only those render);
    // "Show whole structure" leaves it. Disabled until a run is loaded; the
    // dropdown is (re)built lazily on aboutToShow so it tracks the current
    // focus + frame. First QToolButton-with-menu on this toolbar.
    filterMenu_ = new FilterMenu(this);
    filterButton_ = new QToolButton(this);
    filterButton_->setText(QStringLiteral("Filter"));
    filterButton_->setToolButtonStyle(Qt::ToolButtonTextOnly);
    filterButton_->setPopupMode(QToolButton::InstantPopup);
    filterButton_->setMenu(filterMenu_);
    filterButton_->setEnabled(false);
    filterButton_->setToolTip(QStringLiteral(
        "Show only chosen residues near the selected atom and step through "
        "frames isolated. Select an atom first."));
    tb->addWidget(filterButton_);
    ACONNECT(filterMenu_.data(), &QMenu::aboutToShow,
             this, &ReaderMainWindow::buildFilterMenu);

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

    showOccupancyAction_ = tb->addAction(QStringLiteral("Occupancy"));
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

    // Focus — a self-contained toggle at the toolbar tail (deliberately
    // de-emphasized; the Filter feature now covers the common "see one atom +
    // neighbours" need). Checked = track the focused atom, keeping it centred
    // as frames play; unchecked = manual mouse camera. No separate Free/lock
    // mode is surfaced — the button's own state is the only indicator, and it
    // also releases a lock a dashboard reveal engaged.
    tb->addSeparator();
    focusAction_ = tb->addAction(QStringLiteral("Focus"));
    focusAction_->setCheckable(true);
    focusAction_->setEnabled(false);
    focusAction_->setToolTip(QStringLiteral(
        "Track the focused atom — keep it centred as frames play. "
        "Toggle off for free mouse control. Requires a focused atom."));
    ACONNECT(focusAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onFocusCameraTriggered);
}

void ReaderMainWindow::buildStatusBar() {
    selectionLabel_ = new QLabel(QStringLiteral("no selection"), this);
    proteinLabel_   = new QLabel(QStringLiteral("No calcset loaded"), this);
    frameLabel_     = new QLabel(QStringLiteral("frame —"), this);
    timeLabel_      = new QLabel(QStringLiteral("t=— ps"), this);

    // Selection summary on the LEFT; identity/frame/time pinned on the right.
    statusBar()->addWidget(selectionLabel_);
    statusBar()->addPermanentWidget(proteinLabel_);
    statusBar()->addPermanentWidget(frameLabel_);
    statusBar()->addPermanentWidget(timeLabel_);
}

void ReaderMainWindow::buildDocks() {
    // Atom Info dock — tabified on the LEFT alongside Selection + Strip.
    // It is constructed before load and starts with its own placeholder.
    inspectorDock_ = new QtAtomInspectorDock(this);
    addDockWidget(Qt::LeftDockWidgetArea, inspectorDock_);

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

    // Newman dock — tabified alongside; reveals on focus like the Inspector,
    // and redraws the focused residue's backbone phi/psi live as frames play.
    newmanDock_ = new NewmanDock(this);
    addDockWidget(Qt::LeftDockWidgetArea, newmanDock_);
    tabifyDockWidget(inspectorDock_, newmanDock_);

    inspectorDock_->raise();
    resizeDocks({inspectorDock_}, {360}, Qt::Horizontal);

    // Start clean — docks hidden on launch. The Inspector reveals on atom pick,
    // the Strip dock when a metric is added. QSettings restore can override this
    // for users who intentionally left a dock visible.
    inspectorDock_->setVisible(false);
    dashboardStripDock_->setVisible(false);
    newmanDock_->setVisible(false);

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

    // The "Panels" menu/toolbar button was removed: it exposed dock toggles that
    // greyed out with no working route. Docks reveal themselves where it makes
    // sense — the Strip dock when a metric is added, the Inspector on atom pick;
    // the Selection dock was retired (redundant with the in-scene measurements).

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
    // Focus is the only camera control now. Enabled when there's a focused atom
    // to track OR the composer is already tracking (so a lock can always be
    // released); checked iff the composer is tracking anything — any non-Free
    // mode, including one a dashboard reveal engaged. Signal-blocked because we
    // connect via QAction::triggered (user-only), not toggled.
    if (!focusAction_)
        return;
    const bool hasFocus = selection_ && selection_->hasFocus();
    const bool tracking = scene_ && scene_->cameraComposer()
        && scene_->cameraComposer()->mode().kind != CameraMode::Kind::Free;
    focusAction_->setEnabled(scene_ && (hasFocus || tracking));
    const QSignalBlocker block(focusAction_);
    focusAction_->setChecked(tracking);
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

void ReaderMainWindow::onFocusCameraTriggered() {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->cameraComposer()) {
        updateCameraModeActions();
        return;
    }
    const std::size_t t = playback_ ? static_cast<std::size_t>(playback_->currentFrame()) : 0u;
    // Toggle: unchecked → release to manual (Free); modeChanged re-syncs the
    // button. This also releases a lock a dashboard reveal engaged, so the one
    // toggle is the universal "stop tracking" control.
    if (focusAction_ && !focusAction_->isChecked()) {
        scene_->cameraComposer()->setMode(FreeMode(), FreePolicy(), t);
        return;
    }
    if (!selection_ || !selection_->hasFocus() || !loaded_ || !loaded_->protein) {
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
    scene_->cameraComposer()->setMode(result.mode, result.policy, t);
    (void)scene_->cameraComposer()->write(t);
    scene_->syncCameraClippingRange();
    scene_->requestRender(MoleculeScene::RenderSource::External);
    updateCameraModeActions();
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
