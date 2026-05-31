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
#include "QtSelectionOverlay.h"
#include "SelectionDock.h"
#include "DashboardStripDock.h"
#include "DashboardDisplayController.h"
#include "SignalDisplayDialog.h"

#include "../diagnostics/ConnectionAuditor.h"
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
#include "../model/TrajectorySignalCatalog.h"
#include "../model/TransformedConformation.h"

#include <QDockWidget>

#include <QDir>
#include <QFileInfo>

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCloseEvent>
#include <QFileDialog>
#include <QFont>
#include <QKeySequence>
#include <QLabel>
#include <QLoggingCategory>
#include <QMenu>
#include <QMenuBar>
#include <QProcess>
#include <QSet>
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
constexpr int kSettingsVersion = 1;
constexpr int kMaxRecentFiles  = 10;

// Note: locateDftJobsDir was deleted as part of the 2026-05-31 SIMPLIFY
// pass; the DFT campaign now comes from the `.LGS` `dft.frames[]` array
// (see CalcsetManifest + DftShieldingStore).

QUuid addInitialGenericDashboardSignal(model::TrajectorySignalCatalog* catalog,
                                       model::DashboardSignalModel* activeModel,
                                       model::DashboardPanelModel* panelModel,
                                       const QString& descriptorId,
                                       const model::SignalAnchor& anchor,
                                       const QStringList& displayModes,
                                       bool followsFocus,
                                       const QString& label = QString())
{
    if (!catalog || !activeModel)
        return {};
    const model::SignalDescriptor* descriptor = catalog->findDescriptor(descriptorId);
    if (!descriptor)
        return {};
    const QUuid id = activeModel->addSignal(*descriptor, anchor, QString(), displayModes, followsFocus, label);
    if (panelModel && !id.isNull()) {
        panelModel->addDisplayRefs(panelModel->activePanelId(),
                                   model::DisplayRefsForSignal(id, *descriptor, displayModes));
    }
    return id;
}

// Owns the one intentional signal/panel cleanup loop for dashboard models:
// removing a signal removes its display refs, and removing the last display
// ref prunes the now-unreferenced signal.
class DashboardSignalPanelCoordinator final : public QObject {
public:
    DashboardSignalPanelCoordinator(model::DashboardSignalModel* signalModel,
                                    model::DashboardPanelModel* panelModel,
                                    QObject* parent)
        : QObject(parent)
        , signals_(signalModel)
        , panels_(panelModel)
    {
        if (signals_) {
            ACONNECT(signals_.data(), &model::DashboardSignalModel::signalRemoved,
                     this, [this](const QUuid& id) { onSignalRemoved(id); });
        }
        if (panels_) {
            ACONNECT(panels_.data(), &model::DashboardPanelModel::displayRefRemoved,
                     this, [this](const QUuid&, const model::DashboardDisplayRef& ref) {
                         onDisplayRefRemoved(ref);
                     });
        }
    }

private:
    void onSignalRemoved(const QUuid& id)
    {
        if (!panels_ || id.isNull())
            return;
        signalsBeingRemoved_.insert(id);
        panels_->removeDisplayRefsForSignal(id);
        signalsBeingRemoved_.remove(id);
    }

    void onDisplayRefRemoved(const model::DashboardDisplayRef& ref)
    {
        if (!signals_ || !panels_ || ref.signalId.isNull())
            return;
        if (signalsBeingRemoved_.contains(ref.signalId) || !signals_->signalById(ref.signalId))
            return;
        if (panels_->signalReferenceCount(ref.signalId) == 0)
            signals_->removeSignal(ref.signalId);
    }

    QPointer<model::DashboardSignalModel> signals_;
    QPointer<model::DashboardPanelModel> panels_;
    QSet<QUuid> signalsBeingRemoved_;
};
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
    // transform (identity / center_com / fit_reference / fit_subset). Default
    // mode is Identity, so behaviour at startup is identical to today; flipping
    // the mode via POST /transform produces stabilised display positions
    // without re-extracting the trajectory.
    transformed_ = new h5reader::model::TransformedConformation(loaded_->conformation.get(), this);
    ACONNECT(transformed_, &h5reader::model::TransformedConformation::transformChanged,
             this, [this]() {
                 if (scene_) scene_->refreshCurrentFrame();
             });

    // Scene binds to the VTK widget's render window. The scene reads
    // positions through the wrapped conformation so transform mode
    // changes are visible immediately. The widget is passed in so the
    // render scheduler (MoleculeScene::requestRender) can call
    // vtkWidget_->update() — the only render verb in app code per
    // spec/viewport_pipeline_2026-05-30.md §2.5.
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
    dashboardSignals_ = new model::DashboardSignalModel(this);
    dashboardPanels_ = new model::DashboardPanelModel(this);
    signalDisplayDialog_ = new SignalDisplayDialog(this);
    signalDisplayDialog_->setTrajectorySignalCatalog(signalCatalog_);
    signalDisplayDialog_->setDashboardSignalModel(dashboardSignals_);
    signalDisplayDialog_->setDashboardPanelModel(dashboardPanels_);
    signalDisplayDialog_->setContext(loaded_->protein.get(), transformed_);
    signalDisplayDialog_->setSelection(selection_);
    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             signalDisplayDialog_, &SignalDisplayDialog::setFrame);
    addInitialGenericDashboardSignal(signalCatalog_, dashboardSignals_,
                                     dashboardPanels_,
                                     QStringLiteral("npy:dssp_chi"),
                                     model::ResidueAnchor{},
                                     {QStringLiteral("strip.per-class")},
                                     true,
                                     QStringLiteral("Generic NPY DSSP chi"));
    new DashboardSignalPanelCoordinator(dashboardSignals_, dashboardPanels_, this);

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
    dashboardStripDock_->setSelection(selection_);
    dashboardStripDock_->setTimeViewport(timeViewport_);
    dashboardController_ = dashboardStripDock_->displayController();
    addDockWidget(Qt::LeftDockWidgetArea, dashboardStripDock_);
    tabifyDockWidget(inspectorDock_, dashboardStripDock_);
    inspectorDock_->raise();
    resizeDocks({inspectorDock_}, {360}, Qt::Horizontal);

    // Allow the tabbed dock group to shrink to a thin tab-only reminder
    // strip if the user drags the splitter in. By default each dock's
    // contained widget reports a minimumSizeHint from its layout, which
    // QDockWidget honours and prevents below. Override to (0, 0) on the
    // dock AND its inner widget so the splitter's only constraint is the
    // tab bar's own width.
    for (QDockWidget* d : std::vector<QDockWidget*>{
             inspectorDock_, selectionDock_, dashboardStripDock_}) {
        if (!d) continue;
        d->setMinimumSize(0, 0);
        if (QWidget* inner = d->widget())
            inner->setMinimumSize(0, 0);
    }

    // Panel-toggle buttons — appended to the toolbar now that all three
    // docks exist. QDockWidget::toggleViewAction() is the standard Qt
    // primitive for two-way visibility binding; relabel the actions with
    // short text (the dock's title is verbose) and add them as a group
    // at the right end of the toolbar.
    if (playbackToolbar_) {
        playbackToolbar_->addSeparator();
        const struct { QDockWidget* dock; const char* label; const char* tip; } kPanels[] = {
            { inspectorDock_,      "Inspector",
              "Show / hide the Atom Info panel (focus atom's per-frame state)." },
            { selectionDock_,      "Selection",
              "Show / hide the Selected Atoms panel (≤4 ordered atoms with slot colours)." },
            { dashboardStripDock_, "Strip",
              "Show / hide the time-series strip dock." },
        };
        for (const auto& p : kPanels) {
            if (!p.dock) continue;
            QAction* a = p.dock->toggleViewAction();
            a->setText(QString::fromUtf8(p.label));
            a->setToolTip(QString::fromUtf8(p.tip));
            playbackToolbar_->addAction(a);
        }
    }

    ACONNECT(dashboardStripDock_, &DashboardStripDock::revealRequested,
             scene_,              &MoleculeScene::revealBinding);
    // L-3a (2026-05-29): expose the scene's reveal overlay to the
    // dashboard controller so static.tensor mode on Reorient
    // orientation_tensor signals fires an ellipsoid glyph in the
    // 3-D view.
    if (scene_ && scene_->revealOverlay())
        dashboardStripDock_->setSceneOverlay(scene_->revealOverlay());
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
    // instrument + metrics + overlays + panel-toggle row to fit in one
    // toolbar without Qt's overflow chevron. QSettings restore overrides
    // this on later launches.
    resize(1600, 900);
    setWindowTitle(QStringLiteral("h5-reader — %1").arg(loaded_->proteinId));

    // QSettings restore — geometry, dock state, log mask, recent menu.
    // Tolerant: missing or version-mismatched blobs leave the ctor's
    // explicit defaults intact. Runs AFTER all docks/toolbars exist so
    // restoreState has named docks to bind to.
    restoreAllSettings();
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

    // Transform popup-menu button. Drives TransformedConformation::setMode
    // directly. Radio items inside the menu via a second exclusive group;
    // each action's data() carries the underlying Mode enum value.
    transformMenu_ = new QMenu(QStringLiteral("Transform"), this);
    transformGroup_ = new QActionGroup(this);
    transformGroup_->setExclusive(true);

    using TMode = h5reader::model::TransformedConformation::Mode;
    const struct { TMode mode; const char* label; const char* tip; } kTransforms[] = {
        { TMode::Identity,     "Identity",
          "No transform applied; positions come straight from the trajectory." },
        { TMode::CenterCom,    "Center COM",
          "Translate every frame so the centre of mass sits at the origin." },
        { TMode::FitReference, "Fit reference",
          "Kabsch-fit every frame against the reference frame on all atoms." },
        { TMode::FitSubset,    "Fit backbone",
          "Kabsch-fit every frame against the reference frame on backbone atoms only — "
          "removes rigid-body drift while keeping sidechain motion." },
    };
    for (const auto& entry : kTransforms) {
        QAction* act = transformMenu_->addAction(QString::fromUtf8(entry.label));
        act->setCheckable(true);
        act->setData(static_cast<int>(entry.mode));
        act->setToolTip(QString::fromUtf8(entry.tip));
        transformGroup_->addAction(act);
        if (entry.mode == TMode::Identity)
            act->setChecked(true);   // matches TransformedConformation default
        ACONNECT(act, &QAction::triggered, this, [this, act]() {
            applyTransformModeFromAction(act);
        });
    }

    transformButton_ = new QToolButton(tb);
    transformButton_->setText(QStringLiteral("Transform"));
    transformButton_->setToolTip(QStringLiteral(
        "Choose a rigid-body transform applied to positions before display."));
    transformButton_->setPopupMode(QToolButton::InstantPopup);
    transformButton_->setMenu(transformMenu_);
    tb->addWidget(transformButton_);

    // Instrument-mode toggle — toggles MeasurementOverlay focus-only marker.
    // Same code path as POST /selection/instrument; useful for live demos.
    instrumentAction_ = tb->addAction(QStringLiteral("Instrument"));
    instrumentAction_->setCheckable(true);
    instrumentAction_->setToolTip(QStringLiteral(
        "Enable the marker preset on the focus atom: magenta, high opacity, "
        "fixed radius. Used by the harness; useful for live demos."));
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
        case CameraMode::Kind::Bond:
        case CameraMode::Kind::Subset:
            // No dedicated toolbar action; leave the group with nothing
            // checked. Atom/Bond/Subset arrive only via REST or reveal
            // bindings today.
            break;
    }
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

void ReaderMainWindow::applyTransformModeFromAction(QAction* action) {
    ASSERT_THREAD(this);
    if (!action || !transformed_) return;
    bool ok = false;
    const int raw = action->data().toInt(&ok);
    if (!ok) return;
    using TMode = h5reader::model::TransformedConformation::Mode;
    const TMode mode = static_cast<TMode>(raw);
    if (mode == TMode::FitSubset) {
        transformed_->setMode(mode, 0,
            h5reader::model::TransformedConformation::BackboneSubset(*loaded_->protein));
    } else {
        transformed_->setMode(mode);
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
    if (!signalDisplayDialog_)
        return;
    if (!selection_ || !selection_->hasFocus())
        return;
    if (playback_)
        signalDisplayDialog_->setFrame(playback_->currentFrame());
    signalDisplayDialog_->refreshCatalog();
    signalDisplayDialog_->show();
    signalDisplayDialog_->raise();
    signalDisplayDialog_->activateWindow();
}

}  // namespace h5reader::app
