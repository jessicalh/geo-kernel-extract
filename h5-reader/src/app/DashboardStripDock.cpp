#include "DashboardStripDock.h"

#include "DashboardDisplayController.h"
#include "StripStackWidget.h"
#include "TimeViewportController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"
#include "../model/DashboardPanelModel.h"

#include <QAbstractItemModel>
#include <QCheckBox>
#include <QFont>
#include <QFrame>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QSignalBlocker>
#include <QSizePolicy>
#include <QSpinBox>
#include <QTabBar>
#include <QToolButton>
#include <QVBoxLayout>
#include <QWidget>

namespace h5reader::app {

DashboardStripDock::DashboardStripDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Dashboard Strips"), parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("DashboardStripDock"));
    setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    setMinimumHeight(64);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Ignored);
    QFont compactFont = font();
    if (compactFont.pointSize() > 8)
        compactFont.setPointSize(compactFont.pointSize() - 1);
    else if (compactFont.pixelSize() > 10)
        compactFont.setPixelSize(compactFont.pixelSize() - 1);
    setFont(compactFont);

    controller_ = new DashboardDisplayController(this);

    auto* container = new QWidget(this);
    container->setMinimumHeight(0);
    container->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Ignored);
    auto* layout = new QVBoxLayout(container);
    layout->setContentsMargins(4, 4, 4, 4);
    layout->setSpacing(4);

    auto* tabRow = new QHBoxLayout;
    tabRow->setSpacing(4);
    panelTabs_ = new QTabBar(container);
    panelTabs_->setDocumentMode(true);
    panelTabs_->setExpanding(false);
    panelTabs_->setMovable(true);
    panelTabs_->setTabsClosable(true);
    tabRow->addWidget(panelTabs_, 1);
    addPanelButton_ = new QToolButton(container);
    addPanelButton_->setText(QStringLiteral("+"));
    addPanelButton_->setToolTip(QStringLiteral("Create a dashboard panel."));
    tabRow->addWidget(addPanelButton_);
    layout->addLayout(tabRow);

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

    auto* scroll = new QScrollArea(container);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setWidgetResizable(true);
    scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scroll->setMinimumHeight(48);
    scroll->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    stackWidget_ = new StripStackWidget(scroll);
    scroll->setWidget(stackWidget_);
    layout->addWidget(scroll, 1);

    statusLabel_ = new QLabel(QStringLiteral("No active strip signals."), container);
    statusLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(statusLabel_);

    setWidget(container);

    ACONNECT(controller_, &DashboardDisplayController::stripTracksChanged,
             this, &DashboardStripDock::refreshTracks);
    // Owned panels (SequenceBarPanel etc.) survive playhead ticks; only
    // a controller rebuild produces a new set. The dock pulls them out
    // (move-out) only on this signal, NOT on every frame change —
    // otherwise the first setFrame() drains them and the panel
    // disappears.
    ACONNECT(controller_, &DashboardDisplayController::ownedPanelsChanged,
             this, [this]() {
        if (stackWidget_)
            stackWidget_->setOwnedPanels(controller_->takeOwnedPanels());
    });
    ACONNECT(stackWidget_.data(), &StripStackWidget::revealRequested,
             this, &DashboardStripDock::revealRequested);
    ACONNECT(metricButton_.data(), &QPushButton::clicked,
             this, &DashboardStripDock::metricPickerRequested);
    ACONNECT(panelTabs_.data(), &QTabBar::currentChanged,
             this, &DashboardStripDock::onPanelTabChanged);
    ACONNECT(panelTabs_.data(), &QTabBar::tabCloseRequested,
             this, &DashboardStripDock::onPanelTabCloseRequested);
    ACONNECT(addPanelButton_.data(), &QToolButton::clicked,
             this, &DashboardStripDock::onAddPanelRequested);
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

DashboardDisplayController* DashboardStripDock::displayController() const {
    return controller_;
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

void DashboardStripDock::setPanelModel(model::DashboardPanelModel* panelModel) {
    ASSERT_THREAD(this);
    if (panelModel_)
        disconnect(panelModel_, nullptr, this, nullptr);
    panelModel_ = panelModel;
    controller_->setPanelModel(panelModel);
    if (panelModel_) {
        ACONNECT(panelModel_.data(), &QAbstractItemModel::rowsInserted,
                 this, &DashboardStripDock::syncPanelTabs);
        ACONNECT(panelModel_.data(), &QAbstractItemModel::rowsRemoved,
                 this, &DashboardStripDock::syncPanelTabs);
        ACONNECT(panelModel_.data(), &QAbstractItemModel::modelReset,
                 this, &DashboardStripDock::syncPanelTabs);
        ACONNECT(panelModel_.data(), &QAbstractItemModel::dataChanged,
                 this, [this](const QModelIndex&, const QModelIndex&, const QList<int>&) {
                     syncPanelTabs();
                 });
        ACONNECT(panelModel_.data(), &model::DashboardPanelModel::activePanelChanged,
                 this, [this](const QUuid&) { syncPanelTabs(); });
    }
    syncPanelTabs();
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

void DashboardStripDock::setSceneOverlay(SceneRevealOverlay* overlay) {
    ASSERT_THREAD(this);
    controller_->setSceneOverlay(overlay);
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

DashboardSmokeSummary DashboardStripDock::smokeSummary(int firstFrame, int lastFrame) const {
    if (!controller_)
        return {};
    return controller_->smokeSummary(firstFrame, lastFrame);
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
    // Owned panels are NOT touched here. They flow through the
    // separate ownedPanelsChanged → setOwnedPanels path so they
    // survive setFrame()-only ticks.
    stackWidget_->setCurrentFrame(frame_);
    if (statusLabel_)
        statusLabel_->setText(controller_->statusText());
}

void DashboardStripDock::updateViewportReadout(int first, int last) {
    if (viewportReadout_)
        viewportReadout_->setText(QStringLiteral("f%1-f%2").arg(first + 1).arg(last + 1));
}

void DashboardStripDock::syncPanelTabs() {
    ASSERT_THREAD(this);
    if (!panelTabs_)
        return;

    const QSignalBlocker block(panelTabs_);
    while (panelTabs_->count() > 0)
        panelTabs_->removeTab(0);
    if (!panelModel_) {
        panelTabs_->setVisible(false);
        if (addPanelButton_)
            addPanelButton_->setVisible(false);
        return;
    }

    panelTabs_->setVisible(true);
    if (addPanelButton_)
        addPanelButton_->setVisible(true);
    for (int row = 0; row < panelModel_->rowCount(); ++row) {
        const QModelIndex index = panelModel_->index(row, 0);
        const QString name = panelModel_->data(index, model::DashboardPanelModel::NameRole).toString();
        const QUuid id = panelModel_->data(index, model::DashboardPanelModel::UuidRole).toUuid();
        panelTabs_->addTab(name);
        panelTabs_->setTabData(row, id);
    }
    panelTabs_->setTabsClosable(panelModel_->rowCount() > 1);
    const int activeRow = panelModel_->activePanelRow();
    if (activeRow >= 0)
        panelTabs_->setCurrentIndex(activeRow);
}

void DashboardStripDock::onPanelTabChanged(int row) {
    ASSERT_THREAD(this);
    if (!panelModel_ || row < 0)
        return;
    const QUuid id = panelTabs_ ? panelTabs_->tabData(row).toUuid() : QUuid{};
    if (!id.isNull())
        panelModel_->setActivePanel(id);
}

void DashboardStripDock::onPanelTabCloseRequested(int row) {
    ASSERT_THREAD(this);
    if (!panelModel_ || row < 0)
        return;
    panelModel_->removePanelAt(row);
}

void DashboardStripDock::onAddPanelRequested() {
    ASSERT_THREAD(this);
    if (!panelModel_)
        return;
    bool ok = false;
    const QString name = QInputDialog::getText(this,
                                               QStringLiteral("New dashboard panel"),
                                               QStringLiteral("Name"),
                                               QLineEdit::Normal,
                                               QStringLiteral("Panel %1").arg(panelModel_->rowCount() + 1),
                                               &ok);
    if (!ok)
        return;
    const QUuid id = panelModel_->addPanel(name);
    panelModel_->setActivePanel(id);
}

}  // namespace h5reader::app
