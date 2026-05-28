#include "ReaderMainWindow.h"

#include "MoleculeScene.h"
#include "QtAtomInspectorDock.h"
#include "QtAtomPicker.h"
#include "QtAtomTimeSeriesDock.h"
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
#include <QKeySequence>
#include <QLabel>
#include <QLoggingCategory>
#include <QMenuBar>
#include <QProcess>
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

#include <algorithm>
#include <utility>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cWindow, "h5reader.window")
Q_LOGGING_CATEGORY(cDashboardSmoke, "h5reader.dashboard.smoke")

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

bool isStripMode(const QString& mode) {
    return mode.startsWith(QStringLiteral("strip."));
}

bool isSmokeStripMode(const QString& mode) {
    return isStripMode(mode) && mode != QStringLiteral("strip.spectrum");
}

QStringList stripModesForDescriptor(const model::SignalDescriptor& descriptor) {
    QStringList modes;
    for (const QString& mode : model::AllDisplayModes(descriptor)) {
        if (isStripMode(mode) && !modes.contains(mode))
            modes.push_back(mode);
    }
    return modes;
}

QStringList smokeStripModesForDescriptor(const model::SignalDescriptor& descriptor) {
    QStringList modes;
    for (const QString& mode : stripModesForDescriptor(descriptor)) {
        if (isSmokeStripMode(mode))
            modes.push_back(mode);
    }
    return modes;
}

model::SignalAnchor defaultSmokeAnchor(model::SignalAxis axis, const model::QtProtein* protein) {
    switch (axis) {
    case model::SignalAxis::Atom:
        return model::AtomAnchor{0};
    case model::SignalAxis::Residue:
        if (protein && protein->atomCount() > 0) {
            const int residue = protein->atom(0).residueIndex;
            if (residue >= 0)
                return model::ResidueAnchor{static_cast<std::size_t>(residue)};
        }
        return model::ResidueAnchor{0};
    case model::SignalAxis::AtomTuple: {
        std::vector<std::size_t> atoms;
        const std::size_t n = protein ? std::min<std::size_t>(protein->atomCount(), model::AtomSelection::kMaxAtoms)
                                      : 0;
        atoms.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            atoms.push_back(i);
        return model::AtomTupleAnchor{std::move(atoms)};
    }
    case model::SignalAxis::Bond:
        return model::BondAnchor{0};
    case model::SignalAxis::Ring:
        return model::RingAnchor{0};
    case model::SignalAxis::AromaticRing:
        return model::AromaticRingAnchor{0};
    case model::SignalAxis::SaturatedRing:
        return model::SaturatedRingAnchor{0};
    case model::SignalAxis::RingContributionPair:
        return model::RingContributionPairAnchor{0};
    case model::SignalAxis::RingMembership:
        return model::RingMembershipAnchor{0};
    case model::SignalAxis::MutationMatchPair:
        return model::MutationMatchPairAnchor{0};
    case model::SignalAxis::Protein:
        return model::ProteinAnchor{};
    case model::SignalAxis::System:
        return model::SystemAnchor{};
    case model::SignalAxis::Event:
        return model::EventAnchor{};
    case model::SignalAxis::None:
        return model::NoneAnchor{};
    }
    return model::NoneAnchor{};
}

std::vector<int> dashboardSmokeFrames(const model::Conformation* conformation,
                                      int firstFrame,
                                      int frameCount) {
    std::vector<int> frames;
    if (!conformation)
        return frames;

    const int totalFrames = static_cast<int>(conformation->frameCount());
    if (totalFrames <= 0 || firstFrame < 0 || firstFrame >= totalFrames)
        return frames;

    const int availableFrames = totalFrames - firstFrame;
    const int requestedFrames = std::clamp(frameCount > 0 ? frameCount : 10, 1, availableFrames);
    frames.reserve(requestedFrames);
    for (int frame = firstFrame; frame < firstFrame + requestedFrames; ++frame)
        frames.push_back(frame);
    return frames;
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
    ACONNECT(dashboardPanels_, &model::DashboardPanelModel::displayRefRemoved,
             this, [this](const QUuid&, const model::DashboardDisplayRef& ref) {
                 if (dashboardPanels_ && dashboardSignals_
                     && dashboardPanels_->signalReferenceCount(ref.signalId) == 0) {
                     dashboardSignals_->removeSignal(ref.signalId);
                 }
             });

    ACONNECT(picker_,    &QtAtomPicker::atomPicked,
             selection_, &model::AtomSelection::applyPick);
    ACONNECT(picker_, &QtAtomPicker::atomPicked,
             scene_,  &MoleculeScene::clearReveal);

    ACONNECT(selection_, &model::AtomSelection::focusChanged,
             inspectorDock_, &QtAtomInspectorDock::setPickedAtom);
    ACONNECT(selection_, &model::AtomSelection::focusChanged,
             timeSeriesDock_, &QtAtomTimeSeriesDock::setPickedAtom);
    ACONNECT(selection_, &model::AtomSelection::cleared,
             inspectorDock_, &QtAtomInspectorDock::clearSelection);
    ACONNECT(selection_, &model::AtomSelection::cleared,
             timeSeriesDock_, &QtAtomTimeSeriesDock::clearSelection);
    const auto updateMetricAction = [this]() {
        if (signalDisplaysAction_)
            signalDisplaysAction_->setEnabled(selection_ && selection_->hasFocus());
    };
    ACONNECT(selection_, &model::AtomSelection::focusChanged, this, [updateMetricAction](std::size_t) {
        updateMetricAction();
    });
    ACONNECT(selection_, &model::AtomSelection::cleared, this, updateMetricAction);
    updateMetricAction();

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

    // Dashboard strips — active signals from SignalDisplayDialog rendered
    // through one shared strip surface and the shared TimeViewportController.
    dashboardStripDock_ = new DashboardStripDock(this);
    dashboardStripDock_->setContext(loaded_->protein.get(), loaded_->conformation.get());
    dashboardStripDock_->setSignalModels(signalCatalog_, dashboardSignals_);
    dashboardStripDock_->setPanelModel(dashboardPanels_);
    dashboardStripDock_->setSelection(selection_);
    dashboardStripDock_->setTimeViewport(timeViewport_);
    addDockWidget(Qt::RightDockWidgetArea, dashboardStripDock_);
    tabifyDockWidget(inspectorDock_, dashboardStripDock_);
    inspectorDock_->raise();

    ACONNECT(dashboardStripDock_, &DashboardStripDock::revealRequested,
             scene_,              &MoleculeScene::revealBinding);
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

bool ReaderMainWindow::runDashboardPathSmoke(int firstFrame,
                                             int frameCount,
                                             bool requireFrameSnapshots) {
    ASSERT_THREAD(this);

    if (!loaded_ || !loaded_->protein || !loaded_->conformation || !playback_
        || !signalCatalog_ || !dashboardSignals_ || !dashboardPanels_ || !dashboardStripDock_) {
        qCCritical(cDashboardSmoke).noquote()
            << "dashboard path smoke cannot run: main window is not fully wired";
        return false;
    }

    playback_->pause();

    const auto atomCount = loaded_->protein->atomCount();
    if (selection_) {
        selection_->clear();
        if (atomCount > 0) {
            selection_->applyPick(0, Qt::NoModifier);
            const std::size_t selected = std::min<std::size_t>(atomCount, model::AtomSelection::kMaxAtoms);
            for (std::size_t atom = 1; atom < selected; ++atom)
                selection_->applyPick(atom, Qt::ShiftModifier);
        }
    }

    playback_->setFrame(0);
    dashboardStripDock_->setFrame(0);

    dashboardSignals_->clear();
    dashboardPanels_->clear();

    int stripDescriptorCount = 0;
    int activatedDescriptorCount = 0;
    int activatedStripDisplayModeCount = 0;
    int advertisedStripDisplayModeCount = 0;
    int advertisedSpectrumModeCount = 0;
    int skippedSpectrumOnlyCount = 0;
    int bindFailureCount = 0;
    for (const model::SignalDescriptor& descriptor : signalCatalog_->descriptorList()) {
        const QStringList stripModes = stripModesForDescriptor(descriptor);
        if (stripModes.isEmpty())
            continue;

        ++stripDescriptorCount;
        advertisedStripDisplayModeCount += stripModes.size();
        for (const QString& mode : stripModes) {
            if (mode == QStringLiteral("strip.spectrum"))
                ++advertisedSpectrumModeCount;
        }

        const QStringList smokeModes = smokeStripModesForDescriptor(descriptor);
        if (smokeModes.isEmpty()) {
            ++skippedSpectrumOnlyCount;
            continue;
        }

        const model::SignalAxis anchorAxis = descriptor.requiredAnchor != model::SignalAxis::None
                                                 ? descriptor.requiredAnchor
                                                 : descriptor.nativeAxis;
        const QString label = descriptor.label.isEmpty() ? descriptor.id : descriptor.label;
        const model::SignalAnchor anchor = defaultSmokeAnchor(anchorAxis, loaded_->protein.get());
        QStringList acceptedModes;
        for (const QString& mode : smokeModes) {
            model::DisplaySignalBinding binding;
            binding.sourceKind = descriptor.sourceKind;
            binding.descriptorId = descriptor.id;
            binding.conceptKey = descriptor.conceptKey;
            binding.displayModeId = mode;
            binding.anchor = anchor;
            binding.followsFocus = false;
            if (!signalCatalog_->canBind(binding)) {
                ++bindFailureCount;
                qCWarning(cDashboardSmoke).noquote()
                    << "dashboard path smoke descriptor failed binding validation"
                    << "| descriptor=" << descriptor.id
                    << "| mode=" << mode
                    << "| axis=" << model::ToString(anchorAxis);
                continue;
            }
            acceptedModes.push_back(mode);
        }
        if (acceptedModes.isEmpty())
            continue;

        const QUuid signalId =
            dashboardSignals_->addSignal(descriptor,
                                         anchor,
                                         QString(),
                                         acceptedModes,
                                         false,
                                         QStringLiteral("[smoke %1] %2 | %3 | %4 | %5 mode%6")
                                             .arg(activatedDescriptorCount + 1, 3, 10, QChar('0'))
                                             .arg(label,
                                                  model::ToString(descriptor.sourceKind),
                                                  model::ToString(anchorAxis),
                                                  QString::number(acceptedModes.size()),
                                                  acceptedModes.size() == 1 ? QString() : QStringLiteral("s")));
        dashboardPanels_->addDisplayRefs(dashboardPanels_->activePanelId(),
                                         model::DisplayRefsForSignal(signalId, descriptor, acceptedModes));
        ++activatedDescriptorCount;
        activatedStripDisplayModeCount += acceptedModes.size();
    }

    const std::vector<int> framesToRun = dashboardSmokeFrames(loaded_->conformation.get(),
                                                              firstFrame,
                                                              frameCount);
    const int smokeFirstFrame = framesToRun.empty() ? -1 : framesToRun.front();
    const int smokeLastFrame = framesToRun.empty() ? -1 : framesToRun.back();
    const long long expectedWindowFrames = static_cast<long long>(framesToRun.size());
    const long long expectedBufferFrames = framesToRun.empty()
                                               ? 0LL
                                               : static_cast<long long>(smokeLastFrame) + 1;

    qCInfo(cDashboardSmoke).noquote()
        << "dashboard path smoke started"
        << "| catalog_descriptors=" << signalCatalog_->descriptorList().size()
        << "| strip_descriptors=" << stripDescriptorCount
        << "| active_descriptors=" << activatedDescriptorCount
        << "| active_strip_display_modes=" << activatedStripDisplayModeCount
        << "| advertised_strip_display_modes=" << advertisedStripDisplayModeCount
        << "| advertised_spectrum_modes=" << advertisedSpectrumModeCount
        << "| skipped_spectrum_only=" << skippedSpectrumOnlyCount
        << "| bind_failures=" << bindFailureCount
        << "| active_signals=" << dashboardSignals_->rowCount()
        << "| visited_frames=" << framesToRun.size()
        << "| require_frame_snapshots=" << requireFrameSnapshots
        << "| requested_first_frame=" << firstFrame
        << "| requested_frames=" << frameCount
        << "| window_first_frame=" << smokeFirstFrame
        << "| window_last_frame=" << smokeLastFrame
        << "| window_frames=" << expectedWindowFrames
        << "| prefix_backfill_frames=" << std::max(0, smokeFirstFrame)
        << "| buffer_frames=" << expectedBufferFrames;

    int playbackFramesObserved = 0;
    int playbackFramesExpected = 0;
    const QMetaObject::Connection frameObserver =
        QObject::connect(playback_, &QtPlaybackController::frameChanged, this,
                         [&playbackFramesObserved](int) { ++playbackFramesObserved; });

    int snapshotsResident = 0;
    int snapshotsAbsent = 0;
    for (int i = 0; i < static_cast<int>(framesToRun.size()); ++i) {
        const int frame = framesToRun[static_cast<std::size_t>(i)];
        loaded_->conformation->requestSnapshot(static_cast<std::size_t>(frame));
        if (loaded_->conformation->snapshot(static_cast<std::size_t>(frame)))
            ++snapshotsResident;
        else
            ++snapshotsAbsent;

        if (playback_->currentFrame() == frame) {
            dashboardStripDock_->setFrame(frame);
        } else {
            ++playbackFramesExpected;
            playback_->setFrame(frame);
        }
    }
    QObject::disconnect(frameObserver);

    const DashboardSmokeSummary summary = dashboardStripDock_->smokeSummary(smokeFirstFrame, smokeLastFrame);
    const int stripDisplaySinkCount = dashboardStripDock_->stripDisplaySinkCount();
    const int spectrumDisplaySinkCount = dashboardStripDock_->spectrumDisplaySinkCount();
    const long long expectedSamples = static_cast<long long>(summary.seriesCount) * expectedWindowFrames;

    qCInfo(cDashboardSmoke).noquote()
        << "dashboard path smoke summary"
        << "| series=" << summary.seriesCount
        << "| strip_display_sinks=" << stripDisplaySinkCount
        << "| spectrum_display_sinks=" << spectrumDisplaySinkCount
        << "| with_samples=" << summary.seriesWithSamples
        << "| with_valid=" << summary.seriesWithValidSamples
        << "| pending_only=" << summary.seriesPendingOnly
        << "| dense_series=" << summary.denseSeries
        << "| sparse_series=" << summary.sparseSeries
        << "| all_gap_series=" << summary.allGapSeries
        << "| frame_source_absent_series=" << summary.seriesWithFrameSourceAbsentGaps
        << "| frame_npy_frame_source_absent_series=" << summary.frameNpySeriesWithFrameSourceAbsentGaps
        << "| orca_dft_frame_source_absent_series=" << summary.orcaDftSeriesWithFrameSourceAbsentGaps
        << "| source_absent_series=" << summary.seriesWithSourceAbsentGaps
        << "| anchor_unavailable_series=" << summary.seriesWithAnchorUnavailableGaps
        << "| max_gap_run=" << summary.maxLongestGapRun
        << "| mismatched_buffers=" << summary.seriesWithMismatchedBuffers
        << "| samples=" << summary.samples
        << "| channel_values=" << summary.channelValues
        << "| channel_validity=" << summary.channelValidity
        << "| valid=" << summary.validSamples
        << "| gaps=" << summary.gapSamples
        << "| pending_gaps=" << summary.pendingGapSamples
        << "| source_absent_gaps=" << summary.sourceAbsentGapSamples
        << "| frame_source_absent_gaps=" << summary.frameSourceAbsentGapSamples
        << "| frame_npy_frame_source_absent_gaps=" << summary.frameNpyFrameSourceAbsentGapSamples
        << "| orca_dft_frame_source_absent_gaps=" << summary.orcaDftFrameSourceAbsentGapSamples
        << "| anchor_unavailable_gaps=" << summary.anchorUnavailableGapSamples
        << "| invalid=" << summary.invalidSamples
        << "| snapshots_resident=" << snapshotsResident
        << "| snapshots_absent=" << snapshotsAbsent
        << "| playback_frame_changed=" << playbackFramesObserved
        << "| expected_playback_frame_changed=" << playbackFramesExpected;

    for (const DashboardSmokeSummary::SeriesSparseness& s : summary.seriesSparseness) {
        if (s.samples == 0 || (s.gapSamples == 0 && s.invalidSamples == 0))
            continue;
        const QString channel = s.channelLabel.isEmpty() ? s.channelId : s.channelLabel;
        const double validPct = s.samples > 0
                                    ? 100.0 * static_cast<double>(s.validSamples) / static_cast<double>(s.samples)
                                    : 0.0;
        qCInfo(cDashboardSmoke).noquote()
            << "dashboard signal coverage"
            << "| label=" << s.signalLabel
            << "| descriptor=" << s.descriptorId
            << "| concept=" << s.conceptKey
            << "| source=" << s.sourceKind
            << "| storage=" << s.storagePath
            << "| display=" << s.displayModeId
            << "| channel=" << channel
            << "| samples=" << s.samples
            << "| valid=" << s.validSamples
            << "| valid_pct=" << QString::number(validPct, 'f', 1)
            << "| gaps=" << s.gapSamples
            << "| invalid=" << s.invalidSamples
            << "| first_valid=" << s.firstValidFrame
            << "| last_valid=" << s.lastValidFrame
            << "| longest_valid_run=" << s.longestValidRun
            << "| longest_gap_run=" << s.longestGapRun
            << "| pending=" << s.pendingGapSamples
            << "| source_absent=" << s.sourceAbsentGapSamples
            << "| frame_source_absent=" << s.frameSourceAbsentGapSamples
            << "| source_mask_off=" << s.sourceMaskOffGapSamples
            << "| anchor_unavailable=" << s.anchorUnavailableGapSamples
            << "| not_applicable=" << s.notApplicableGapSamples
            << "| nan=" << s.nanSentinelGapSamples
            << "| malformed=" << s.malformedSourceGapSamples;
    }

    if (summary.validSamples == 0) {
        qCWarning(cDashboardSmoke).noquote()
            << "dashboard path smoke produced no valid samples; pending gaps are expected until source samplers are wired";
    }

    const bool frameSnapshotRequirementOk =
        !requireFrameSnapshots || (snapshotsAbsent == 0 && summary.frameNpyFrameSourceAbsentGapSamples == 0);

    const bool ok = stripDescriptorCount > 0
                    && !framesToRun.empty()
                    && activatedDescriptorCount > 0
                    && bindFailureCount == 0
                    && dashboardSignals_->rowCount() > 0
                    && summary.seriesCount > 0
                    && stripDisplaySinkCount == summary.seriesCount
                    && spectrumDisplaySinkCount == 0
                    && summary.seriesWithMismatchedBuffers == 0
                    && summary.invalidSamples == 0
                    && summary.samples == expectedSamples
                    && summary.channelValues == expectedSamples
                    && summary.channelValidity == expectedSamples
                    && playbackFramesObserved == playbackFramesExpected
                    && frameSnapshotRequirementOk;
    if (!ok) {
        qCCritical(cDashboardSmoke).noquote()
            << "dashboard path smoke failed"
            << "| expected_samples=" << expectedSamples
            << "| frame_snapshot_requirement_ok=" << frameSnapshotRequirementOk;
    }
    return ok;
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
