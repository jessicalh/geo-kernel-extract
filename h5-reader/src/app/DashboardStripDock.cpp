#include "DashboardStripDock.h"

#include "DashboardDisplayController.h"
#include "StripStackWidget.h"
#include "TimeViewportController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"

#include <QCheckBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QVBoxLayout>
#include <QWidget>

namespace h5reader::app {

DashboardStripDock::DashboardStripDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Dashboard Strips"), parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("DashboardStripDock"));
    setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);

    controller_ = new DashboardDisplayController(this);

    auto* container = new QWidget(this);
    auto* layout = new QVBoxLayout(container);
    layout->setContentsMargins(4, 4, 4, 4);
    layout->setSpacing(4);

    auto* top = new QHBoxLayout;
    followBox_ = new QCheckBox(QStringLiteral("Follow"), container);
    followBox_->setChecked(true);
    followBox_->setToolTip(QStringLiteral("Keep the visible strip window attached to the playback frame."));
    top->addWidget(followBox_);

    top->addWidget(new QLabel(QStringLiteral("Window"), container));
    windowFramesSpin_ = new QSpinBox(container);
    windowFramesSpin_->setRange(8, 1000000);
    windowFramesSpin_->setSingleStep(10);
    windowFramesSpin_->setAccelerated(true);
    windowFramesSpin_->setSuffix(QStringLiteral(" frames"));
    windowFramesSpin_->setToolTip(QStringLiteral("Number of frames visible in follow mode."));
    top->addWidget(windowFramesSpin_);

    viewportReadout_ = new QLabel(QStringLiteral("f1-f1"), container);
    viewportReadout_->setToolTip(QStringLiteral("Visible frame range on the dashboard strips."));
    top->addWidget(viewportReadout_);
    top->addStretch(1);
    metricButton_ = new QPushButton(QStringLiteral("Metrics..."), container);
    metricButton_->setEnabled(false);
    metricButton_->setToolTip(QStringLiteral("Select a nearby atom or residue and add a metric strip."));
    top->addWidget(metricButton_);
    layout->addLayout(top);

    stackWidget_ = new StripStackWidget(container);
    layout->addWidget(stackWidget_, 1);

    statusLabel_ = new QLabel(QStringLiteral("No active strip signals."), container);
    statusLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(statusLabel_);

    setWidget(container);

    ACONNECT(controller_, &DashboardDisplayController::stripTracksChanged,
             this, &DashboardStripDock::refreshTracks);
    ACONNECT(stackWidget_.data(), &StripStackWidget::revealRequested,
             this, &DashboardStripDock::revealRequested);
    ACONNECT(metricButton_.data(), &QPushButton::clicked,
             this, &DashboardStripDock::metricPickerRequested);
    ACONNECT(followBox_.data(), &QCheckBox::toggled, this, [this](bool on) {
        if (windowFramesSpin_)
            windowFramesSpin_->setEnabled(on);
        if (timeViewport_)
            timeViewport_->setFollowPlayhead(on);
    });
    ACONNECT(windowFramesSpin_.data(), qOverload<int>(&QSpinBox::valueChanged), this, [this](int frames) {
        if (timeViewport_)
            timeViewport_->setWindowFrames(frames);
    });
}

void DashboardStripDock::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    ASSERT_THREAD(this);
    controller_->setContext(protein, conformation);
}

void DashboardStripDock::setSignalModels(model::TrajectorySignalCatalog* catalog,
                                         model::DashboardSignalModel* activeModel) {
    ASSERT_THREAD(this);
    controller_->setSignalModels(catalog, activeModel);
}

void DashboardStripDock::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    if (selection_)
        disconnect(selection_, nullptr, this, nullptr);
    selection_ = selection;
    controller_->setSelection(selection);
    const auto updateButton = [this]() {
        if (metricButton_)
            metricButton_->setEnabled(selection_ && selection_->hasFocus());
    };
    if (selection_) {
        ACONNECT(selection_.data(), &model::AtomSelection::focusChanged, this, [updateButton](std::size_t) {
            updateButton();
        });
        ACONNECT(selection_.data(), &model::AtomSelection::cleared, this, updateButton);
    }
    updateButton();
}

void DashboardStripDock::setDftStore(model::DftShieldingStore* store) {
    ASSERT_THREAD(this);
    controller_->setDftStore(store);
}

void DashboardStripDock::setTimeViewport(TimeViewportController* viewport) {
    ASSERT_THREAD(this);
    if (timeViewport_)
        disconnect(timeViewport_, nullptr, this, nullptr);
    timeViewport_ = viewport;
    if (timeViewport_) {
        ACONNECT(timeViewport_.data(), &TimeViewportController::visibleRangeChanged,
                 this, &DashboardStripDock::updateViewportReadout);
        ACONNECT(timeViewport_.data(), &TimeViewportController::windowFramesChanged,
                 this, [this](int frames) {
                     if (!windowFramesSpin_)
                         return;
                     const QSignalBlocker block(windowFramesSpin_);
                     windowFramesSpin_->setValue(frames);
                 });
        ACONNECT(timeViewport_.data(), &TimeViewportController::followPlayheadChanged,
                 this, [this](bool follow) {
                     if (!followBox_)
                         return;
                     const QSignalBlocker block(followBox_);
                     followBox_->setChecked(follow);
                     if (windowFramesSpin_)
                         windowFramesSpin_->setEnabled(follow);
                 });
        if (windowFramesSpin_) {
            windowFramesSpin_->setRange(1, timeViewport_->frameCount());
            const QSignalBlocker block(windowFramesSpin_);
            windowFramesSpin_->setValue(timeViewport_->windowFrames());
        }
        if (followBox_) {
            const QSignalBlocker block(followBox_);
            followBox_->setChecked(timeViewport_->followsPlayhead());
        }
        updateViewportReadout(timeViewport_->visibleStart(), timeViewport_->visibleEnd());
    }
    if (stackWidget_)
        stackWidget_->setTimeViewport(timeViewport_);
}

DashboardSmokeSummary DashboardStripDock::smokeSummary() const {
    if (!controller_)
        return {};
    return controller_->smokeSummary();
}

int DashboardStripDock::stripDisplaySinkCount() const {
    return stackWidget_ ? stackWidget_->trackCount() : 0;
}

int DashboardStripDock::spectrumDisplaySinkCount() const {
    return stackWidget_ ? stackWidget_->spectrumTrackCount() : 0;
}

void DashboardStripDock::setFrame(int frame) {
    ASSERT_THREAD(this);
    frame_ = frame;
    controller_->setFrame(frame);
    if (stackWidget_)
        stackWidget_->setCurrentFrame(frame);
}

void DashboardStripDock::refreshTracks() {
    ASSERT_THREAD(this);
    if (!stackWidget_)
        return;

    QVector<StripStackWidget::Track> tracks;
    const QVector<DashboardDisplayController::StripTrack> controllerTracks = controller_->stripTracks();
    tracks.reserve(controllerTracks.size());
    for (const DashboardDisplayController::StripTrack& source : controllerTracks) {
        StripStackWidget::Track track;
        track.buffer = source.buffer;
        track.color = source.color;
        track.hasBinding = source.hasBinding;
        track.binding = source.binding;
        tracks.push_back(track);
    }
    stackWidget_->setTracks(std::move(tracks));
    stackWidget_->setSpectrumTracks({});
    stackWidget_->setCurrentFrame(frame_);
    if (statusLabel_)
        statusLabel_->setText(controller_->statusText());
}

void DashboardStripDock::updateViewportReadout(int first, int last) {
    if (viewportReadout_)
        viewportReadout_->setText(QStringLiteral("f%1-f%2").arg(first + 1).arg(last + 1));
}

}  // namespace h5reader::app
