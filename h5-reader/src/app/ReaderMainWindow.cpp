#include "ReaderMainWindow.h"

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

#include <QDir>
#include <QFileInfo>

#include <QAction>
#include <QApplication>
#include <QFileDialog>
#include <QFont>
#include <QKeySequence>
#include <QLabel>
#include <QLoggingCategory>
#include <QMenuBar>
#include <QProcess>
#include <QSet>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QStringList>
#include <QStatusBar>
#include <QStyle>
#include <QTimer>
#include <QToolBar>
#include <QUuid>

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

// Locate the DFT campaign's jobs dir relative to the opened run path, by a
// BOUNDED documented-convention check (not file discovery): the dataset root
// holds extract/ and dft/ as siblings, so dft/jobs is either <runDir>/dft/jobs
// or <runDir>/../dft/jobs (when the run path points inside extract/). Empty if
// neither exists — then the run has no DFT and the shielding panel stays hidden.
// A head-of-directory TOML descriptor (#25) will state this path explicitly.
QString locateDftJobsDir(const QString& runPath) {
    if (runPath.isEmpty()) return {};
    const QFileInfo fi(runPath);
    const QDir runDir = fi.isDir() ? QDir(runPath) : fi.absoluteDir();
    for (const QString& rel : {QStringLiteral("dft/jobs"), QStringLiteral("../dft/jobs")}) {
        const QString cand = runDir.absoluteFilePath(rel);
        if (QFileInfo(cand).isDir()) return QDir(cand).absolutePath();
    }
    return {};
}

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

    // Scene binds to the VTK widget's render window.
    scene_ = new MoleculeScene(renderWindow_, this);
    scene_->Build(*loaded_->protein, *loaded_->conformation);
    scene_->ResetCamera();
    ACONNECT(scene_, &MoleculeScene::cameraPlaneLockChanged,
             this, [this](bool active) {
                 if (planeLockAction_ && planeLockAction_->isChecked() != active) {
                     const QSignalBlocker blocker(planeLockAction_);
                     planeLockAction_->setChecked(active);
                 }
                 updatePlaneLockAction();
             });

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
                                loaded_->conformation.get(),
                                playback_, this);

    // Atom Info dock — tabified on the right. Tracks the selection's
    // FOCUS atom (one atom's full per-frame pile). Starts empty; fills in on
    // the first pick.
    inspectorDock_ = new QtAtomInspectorDock(this);
    inspectorDock_->setContext(loaded_->protein.get(),
                                loaded_->conformation.get());
    addDockWidget(Qt::RightDockWidgetArea, inspectorDock_);

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
    signalDisplayDialog_->setContext(loaded_->protein.get(), loaded_->conformation.get());
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

    ACONNECT(selection_, &model::AtomSelection::focusChanged,
             inspectorDock_, &QtAtomInspectorDock::setPickedAtom);
    ACONNECT(selection_, &model::AtomSelection::cleared,
             inspectorDock_, &QtAtomInspectorDock::clearSelection);
    const auto updateMetricAction = [this]() {
        if (signalDisplaysAction_)
            signalDisplaysAction_->setEnabled(selection_ && selection_->hasFocus());
    };
    ACONNECT(selection_, &model::AtomSelection::focusChanged, this, [updateMetricAction](std::size_t) {
        updateMetricAction();
    });
    ACONNECT(selection_, &model::AtomSelection::cleared, this, updateMetricAction);
    updateMetricAction();
    updatePlaneLockAction();

    if (auto* meas = scene_->measurementOverlay()) {
        meas->setSelection(selection_);
        ACONNECT(selection_, &model::AtomSelection::changed,
                 meas,       &MeasurementOverlay::onSelectionChanged);
    }
    ACONNECT(selection_, &model::AtomSelection::changed,
             this, [this]() {
                 if (planeLockAction_ && planeLockAction_->isChecked())
                     planeLockAction_->setChecked(false);
                 updatePlaneLockAction();
                 if (scene_) scene_->refreshCurrentFrame();
             });

    // Selected-atoms panel — the QListView bound to the AtomSelection model
    // (slot colour swatch + residue:atom label + geometry kind). The model's
    // first view; tabified with Atom Info.
    selectionDock_ = new SelectionDock(this);
    selectionDock_->setModel(selection_);
    addDockWidget(Qt::RightDockWidgetArea, selectionDock_);
    tabifyDockWidget(inspectorDock_, selectionDock_);
    inspectorDock_->raise();

    // Dashboard strips — active signals from SignalDisplayDialog rendered
    // through one shared strip surface and the shared TimeViewportController.
    dashboardStripDock_ = new DashboardStripDock(this);
    dashboardStripDock_->setContext(loaded_->protein.get(), loaded_->conformation.get());
    dashboardStripDock_->setSignalModels(signalCatalog_, dashboardSignals_);
    dashboardStripDock_->setPanelModel(dashboardPanels_);
    dashboardStripDock_->setSelection(selection_);
    dashboardStripDock_->setTimeViewport(timeViewport_);
    addDockWidget(Qt::BottomDockWidgetArea, dashboardStripDock_);
    resizeDocks({dashboardStripDock_}, {300}, Qt::Vertical);
    resizeDocks({inspectorDock_}, {320}, Qt::Horizontal);

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

    // DFT shielding campaign (optional): make the frame-local source available
    // to descriptor-family samplers. Its presence does not auto-pin a DFT strip.
    if (const QString dftJobs = locateDftJobsDir(loaded_->runPath); !dftJobs.isEmpty()) {
        dftStore_ = new model::DftShieldingStore(loaded_->protein.get(), dftJobs, this);
        dashboardStripDock_->setDftStore(dftStore_);
        qCInfo(cWindow).noquote() << "DFT shielding store wired |" << dftJobs
                                  << "| jobs=" << dftStore_->jobCount();
    }

    // Initial status bar population.
    onFrameChanged(0);

    resize(1200, 800);
    setWindowTitle(QStringLiteral("h5-reader — %1").arg(loaded_->proteinId));

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
                            this);
    const quint16 bound = restServer_->listen(port);
    if (bound == 0) {
        qCCritical(cWindow).noquote() << "REST server failed to bind port" << port;
        restServer_->deleteLater();
        restServer_ = nullptr;
    }
    return bound;
}

void ReaderMainWindow::shutdown() {
    ASSERT_THREAD(this);
    if (shutdownDone_) return;
    shutdownDone_ = true;

    qCInfo(cWindow).noquote() << "shutdown entered";

    // 1. Stop every timer owned by us or our children. The generic
    //    findChildren sweep catches QtPlaybackController's timer too.
    const auto timers = findChildren<QTimer*>();
    for (auto* timer : timers) {
        if (timer->isActive()) timer->stop();
    }

    // 2. Finalise VTK before Qt destroys the GL context. Without this,
    //    vtkSmartPointer destructors downstream touch dead OpenGL
    //    resources and crash.
    if (renderWindow_) {
        renderWindow_->Finalize();
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
}

void ReaderMainWindow::buildToolbar() {
    auto* tb = addToolBar(QStringLiteral("Playback"));
    tb->setObjectName(QStringLiteral("PlaybackToolbar"));
    tb->setMovable(false);
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

    signalDisplaysAction_ = tb->addAction(QStringLiteral("Metrics..."));
    signalDisplaysAction_->setEnabled(false);
    signalDisplaysAction_->setToolTip(QStringLiteral("Select a nearby atom or residue and add a metric display."));
    ACONNECT(signalDisplaysAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onOpenSignalDisplays);

    planeLockAction_ = tb->addAction(QStringLiteral("Plane lock"));
    planeLockAction_->setCheckable(true);
    planeLockAction_->setEnabled(false);
    planeLockAction_->setToolTip(QStringLiteral(
        "Keep the view centred and oriented to the plane defined by exactly three selected atoms."));
    ACONNECT(planeLockAction_.data(), &QAction::toggled,
             this, &ReaderMainWindow::onPlaneLockToggled);

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

void ReaderMainWindow::updatePlaneLockAction() {
    if (!planeLockAction_)
        return;
    const bool canLock = scene_ && selection_ && selection_->count() == 3;
    planeLockAction_->setEnabled(canLock);
}

void ReaderMainWindow::onPlaneLockToggled(bool checked) {
    ASSERT_THREAD(this);
    if (!scene_)
        return;

    if (!checked) {
        scene_->clearCameraPlaneLock();
        updatePlaneLockAction();
        return;
    }

    if (!selection_ || selection_->count() != 3
        || !scene_->lockCameraToSelectionPlane(selection_->atoms())) {
        const QSignalBlocker blocker(planeLockAction_);
        planeLockAction_->setChecked(false);
        scene_->clearCameraPlaneLock();
    }
    updatePlaneLockAction();
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
