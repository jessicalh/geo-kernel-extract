#include "ReaderMainWindow.h"

#include "MoleculeScene.h"
#include "QtAtomInspectorDock.h"
#include "QtAtomPicker.h"
#include "QtAtomTimeSeriesDock.h"
#include "QtBackboneRibbonOverlay.h"
#include "QtBFieldStreamOverlay.h"
#include "QtFieldGridOverlay.h"
#include "QtPlaybackController.h"
#include "QtRingPolygonOverlay.h"
#include "QtSelectionOverlay.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/QtProteinLoader.h"
#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QAction>
#include <QApplication>
#include <QLabel>
#include <QLoggingCategory>
#include <QSlider>
#include <QSpinBox>
#include <QStatusBar>
#include <QStyle>
#include <QTimer>
#include <QToolBar>

#include <QVTKOpenGLNativeWidget.h>

#include <vtkRendererCollection.h>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cWindow, "h5reader.window")
}

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

    // Playback controller — frameChanged drives the scene, which drives
    // the render. Toolbar controls drive the playback.
    const int T = static_cast<int>(loaded_->conformation->frameCount());
    playback_ = new QtPlaybackController(T, this);

    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             scene_,    &MoleculeScene::setFrame);
    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             this,      &ReaderMainWindow::onFrameChanged);
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

    // Atom picker — event filter on the VTK widget. Emits atomPicked
    // on double-click. Two listeners:
    //   - MoleculeScene's selection overlay (yellow sphere on the atom).
    //   - QtAtomInspectorDock (rebuild tree with per-frame values).
    // Pull the renderer smart-ptr out of the render window so the
    // picker holds the same instance MoleculeScene installed.
    auto* firstRenderer = renderWindow_->GetRenderers()->GetFirstRenderer();
    picker_ = new QtAtomPicker(vtkWidget_, firstRenderer,
                                loaded_->protein.get(),
                                loaded_->conformation.get(),
                                playback_, this);
    if (auto* sel = scene_->selectionOverlay()) {
        ACONNECT(picker_, &QtAtomPicker::atomPicked,
                 sel,     &QtSelectionOverlay::setPickedAtom);
    }

    // Atom inspector dock — tabified on the right. Starts with a
    // placeholder; fills in when the first pick arrives.
    inspectorDock_ = new QtAtomInspectorDock(this);
    inspectorDock_->setContext(loaded_->protein.get(),
                                loaded_->conformation.get());
    addDockWidget(Qt::RightDockWidgetArea, inspectorDock_);

    ACONNECT(picker_,    &QtAtomPicker::atomPicked,
             inspectorDock_, &QtAtomInspectorDock::setPickedAtom);
    ACONNECT(playback_,  &QtPlaybackController::frameChanged,
             inspectorDock_, &QtAtomInspectorDock::setFrame);

    // Time-series dock — per-atom scalar-vs-frame line chart via
    // Qt6 Charts. Tabified with the inspector so they share the
    // right dock area; user clicks the tab to switch between
    // "current values" and "trace across the trajectory".
    timeSeriesDock_ = new QtAtomTimeSeriesDock(this);
    timeSeriesDock_->setContext(loaded_->protein.get(),
                                 loaded_->conformation.get());
    addDockWidget(Qt::RightDockWidgetArea, timeSeriesDock_);
    tabifyDockWidget(inspectorDock_, timeSeriesDock_);
    inspectorDock_->raise();   // inspector on top by default

    ACONNECT(picker_,    &QtAtomPicker::atomPicked,
             timeSeriesDock_, &QtAtomTimeSeriesDock::setPickedAtom);
    ACONNECT(playback_,  &QtPlaybackController::frameChanged,
             timeSeriesDock_, &QtAtomTimeSeriesDock::setFrame);

    // One more wire: when the selection overlay receives a pick, the
    // scene needs to reposition the sphere for the CURRENT frame.
    // Hook it through the scene's render-request path — on pick,
    // re-run setFrame so the selection overlay sees the pick AND
    // the position update in one pass.
    ACONNECT(picker_, &QtAtomPicker::atomPicked,
             this,   [this](std::size_t /*a*/) {
                 if (scene_ && playback_)
                     scene_->refreshCurrentFrame();
             });

    // Initial status bar population.
    onFrameChanged(0);

    resize(1200, 800);
    setWindowTitle(QStringLiteral("h5-reader — %1").arg(loaded_->proteinId));

    qCInfo(cWindow).noquote() << "ctor done";
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
}

void ReaderMainWindow::buildToolbar() {
    auto* tb = addToolBar(QStringLiteral("Playback"));
    tb->setObjectName(QStringLiteral("PlaybackToolbar"));
    tb->setMovable(false);

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
    const double t_ps = loaded_->conformation->frame(
        static_cast<size_t>(std::clamp(t, 0, T - 1))).timePicoseconds();

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

void ReaderMainWindow::onPlayPauseClicked() {
    ASSERT_THREAD(this);
    if (playback_) playback_->togglePlayPause();
}

void ReaderMainWindow::onFpsChanged(int /*fps*/) {
    // Reserved for future display (e.g., a "current fps" readout that
    // differs from the requested fps when frame rendering is slower
    // than the interval).
}

}  // namespace h5reader::app
