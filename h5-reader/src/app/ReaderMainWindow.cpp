#include "ReaderMainWindow.h"

#include "MoleculeScene.h"
#include "QtAtomInspectorDock.h"
#include "QtAtomPicker.h"
#include "QtAtomTimeSeriesDock.h"
#include "QtBackboneRibbonOverlay.h"
#include "QtBFieldStreamOverlay.h"
#include "QtFieldGridOverlay.h"
#include "QtPlaybackController.h"
#include "MeasurementOverlay.h"
#include "QtRingPolygonOverlay.h"
#include "QtSelectionOverlay.h"
#include "SelectionDock.h"
#include "StripChartDock.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/QtProteinLoader.h"
#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/DftShieldingStore.h"
#include "../model/QtProtein.h"

#include <QDir>
#include <QFileInfo>

#include <QAction>
#include <QApplication>
#include <QFileDialog>
#include <QKeySequence>
#include <QLabel>
#include <QLoggingCategory>
#include <QMenuBar>
#include <QProcess>
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

    // Atom picker — event filter on the VTK widget. Emits
    // atomPicked(idx, modifiers) on double-click. It stays dumb: it does NOT
    // interpret the gesture. AtomSelection (below) is the sole consumer; it
    // applies the plain/Shift policy and fans typed changes to the inspector,
    // the time-series dock, and the measurement overlay. Pull the renderer
    // smart-ptr out of the render window so the picker holds the same
    // instance MoleculeScene installed.
    auto* firstRenderer = renderWindow_->GetRenderers()->GetFirstRenderer();
    picker_ = new QtAtomPicker(vtkWidget_, firstRenderer,
                                loaded_->protein.get(),
                                loaded_->conformation.get(),
                                playback_, this);

    // Atom inspector dock — tabified on the right. Tracks the selection's
    // FOCUS atom (one atom's full per-frame pile). Starts empty; fills in on
    // the first pick.
    inspectorDock_ = new QtAtomInspectorDock(this);
    inspectorDock_->setContext(loaded_->protein.get(),
                                loaded_->conformation.get());
    addDockWidget(Qt::RightDockWidgetArea, inspectorDock_);

    ACONNECT(playback_,  &QtPlaybackController::frameChanged,
             inspectorDock_, &QtAtomInspectorDock::setFrame);

    // Time-series dock — per-atom scalar-vs-frame line chart via Qt6 Charts.
    // Also tracks the FOCUS atom. Tabified with the inspector.
    timeSeriesDock_ = new QtAtomTimeSeriesDock(this);
    timeSeriesDock_->setContext(loaded_->protein.get(),
                                 loaded_->conformation.get());
    addDockWidget(Qt::RightDockWidgetArea, timeSeriesDock_);
    tabifyDockWidget(inspectorDock_, timeSeriesDock_);
    inspectorDock_->raise();   // inspector on top by default

    ACONNECT(playback_,  &QtPlaybackController::frameChanged,
             timeSeriesDock_, &QtAtomTimeSeriesDock::setFrame);

    // ---- Selection model — the single source of selection truth ----------
    //
    // The picker reports a pick + its keyboard modifiers; AtomSelection
    // interprets the gesture (plain = replace the focus; Shift = toggle the
    // atom in the ≤4 ordered set) and is itself the QAbstractListModel the
    // SelectionDock view binds to. It fans typed changes out:
    //   focusChanged → inspector + time-series retarget to the focus atom;
    //   cleared      → those two clear;
    //   changed      → the measurement overlay rebuilds, and the scene
    //                  refreshes the current frame so the spheres reposition.
    selection_ = new model::AtomSelection(loaded_->protein.get(), this);

    ACONNECT(picker_,    &QtAtomPicker::atomPicked,
             selection_, &model::AtomSelection::applyPick);

    ACONNECT(selection_, &model::AtomSelection::focusChanged,
             inspectorDock_, &QtAtomInspectorDock::setPickedAtom);
    ACONNECT(selection_, &model::AtomSelection::focusChanged,
             timeSeriesDock_, &QtAtomTimeSeriesDock::setPickedAtom);
    ACONNECT(selection_, &model::AtomSelection::cleared,
             inspectorDock_, &QtAtomInspectorDock::clearSelection);
    ACONNECT(selection_, &model::AtomSelection::cleared,
             timeSeriesDock_, &QtAtomTimeSeriesDock::clearSelection);

    if (auto* meas = scene_->measurementOverlay()) {
        meas->setSelection(selection_);
        ACONNECT(selection_, &model::AtomSelection::changed,
                 meas,       &MeasurementOverlay::onSelectionChanged);
    }
    ACONNECT(selection_, &model::AtomSelection::changed,
             this, [this]() {
                 if (scene_) scene_->refreshCurrentFrame();
             });

    // Selected-atoms panel — the QListView bound to the AtomSelection model
    // (slot colour swatch + residue:atom label + geometry kind). The model's
    // first view; tabified with the inspector + time-series.
    selectionDock_ = new SelectionDock(this);
    selectionDock_->setModel(selection_);
    addDockWidget(Qt::RightDockWidgetArea, selectionDock_);
    tabifyDockWidget(inspectorDock_, selectionDock_);
    inspectorDock_->raise();

    // Strip-chart dock — the trajectory geometry instrument (QtCharts, NOT a
    // second VTK surface; decision 2026-05-27). Charts the SELECTION's derived
    // geometry (distance / angle / dihedral) over frames with a scrolling
    // window, science grids, and a digital readout. setSelection ACONNECTs the
    // selection's changed()/cleared() itself; playback drives the playhead.
    stripChartDock_ = new StripChartDock(this);
    stripChartDock_->setContext(loaded_->protein.get(), loaded_->conformation.get());
    stripChartDock_->setSelection(selection_);
    addDockWidget(Qt::RightDockWidgetArea, stripChartDock_);
    tabifyDockWidget(inspectorDock_, stripChartDock_);
    inspectorDock_->raise();

    ACONNECT(playback_,       &QtPlaybackController::frameChanged,
             stripChartDock_, &StripChartDock::setFrame);

    // DFT shielding campaign (optional): if the run sits in a dataset with a
    // sibling dft/jobs directory, wire the per-frame ORCA shielding into the
    // strip chart's shielding panel. Absent dft/ leaves that panel hidden.
    if (const QString dftJobs = locateDftJobsDir(loaded_->runPath); !dftJobs.isEmpty()) {
        dftStore_ = new model::DftShieldingStore(loaded_->protein.get(), dftJobs, this);
        stripChartDock_->setDftStore(dftStore_);
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

void ReaderMainWindow::onPlayPauseClicked() {
    ASSERT_THREAD(this);
    if (playback_) playback_->togglePlayPause();
}

void ReaderMainWindow::onFpsChanged(int /*fps*/) {
    // Reserved for future display (e.g., a "current fps" readout that
    // differs from the requested fps when frame rendering is slower
    // than the interval).
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

}  // namespace h5reader::app
