#include "StripChartDock.h"

#include "StripStackWidget.h"
#include "TimeViewportController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"
#include "../model/DftShieldingStore.h"
#include "../model/QtProtein.h"
#include "../model/QtResidueNames.h"
#include "../model/SpectralAnalysis.h"
#include "../model/TrajectoryConformation.h"

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

#include <QCheckBox>
#include <QColor>
#include <QHBoxLayout>
#include <QLabel>
#include <QList>
#include <QLoggingCategory>
#include <QMargins>
#include <QPainter>
#include <QPen>
#include <QPointF>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QStringList>
#include <QVBoxLayout>
#include <QWidget>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cStrip, "h5reader.stripchart")

// Scrolling time window (frames). When more than this is collected, the x-axis
// shows a fixed-width window following the playhead instead of scrunching the
// whole trace into the panel. "Fit all" overrides.
constexpr int kWindowFrames = 300;

// Display cap for ONE rendered window. The authoritative ChannelBuffer is never
// decimated (the user's directive); but pushing a million QPointF into a
// software QLineSeries would choke Qt Charts, so when a single visible window
// exceeds this (only "fit all" on a microsecond run), we stride the DISPLAY.
// At demo scale (≤ a few thousand) stride is 1 and every point is drawn.
constexpr int kMaxDisplayPoints = 4000;

// FFT recompute throttle: redo the spectrum every this-many newly collected
// frames during play (plus always on selection change). Cheap at demo scale; a
// microsecond run should move continuous-play FFT to an on-pause recompute.
constexpr int kFftRecomputeStride = 64;
constexpr double kRadToDeg = 180.0 / 3.141592653589793238462643383279502884;

const QColor kGeomColor(86, 166, 244);    // structural trace: high-contrast blue
const QColor kDftColor(255, 179, 87);     // shielding trace: high-contrast amber
const QColor kCursorColor(217, 26, 26);   // playhead: red dashed
const QColor kFftColor(0, 128, 96);       // spectrum: teal, distinct from above
const QColor kPhiColor(229, 99, 99);
const QColor kPsiColor(120, 184, 92);
const QColor kOmegaColor(180, 131, 230);
const QColor kChiColor(94, 170, 220);

int decimalsForUnit(const QString& unit) {
    if (unit == QString::fromUtf8("Å"))  return 3;
    if (unit == QStringLiteral("ppm"))   return 2;
    return 1;  // degrees
}

double angleToDegrees(double value, const QString& units) {
    const QString u = units.trimmed().toLower();
    if (u.isEmpty() || u == QStringLiteral("rad") || u == QStringLiteral("radian")
        || u == QStringLiteral("radians"))
        return value * kRadToDeg;
    return value;
}
}  // namespace

StripChartDock::StripChartDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Strip Chart"), parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("StripChartDock"));
    setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);

    auto* container = new QWidget(this);
    auto* vbox      = new QVBoxLayout(container);
    vbox->setContentsMargins(4, 4, 4, 4);
    vbox->setSpacing(4);

    // Top row: collected-history overview vs trailing strip window.
    auto* top  = new QHBoxLayout;
    followBox_ = new QCheckBox(QStringLiteral("Follow"), container);
    followBox_->setChecked(true);
    followBox_->setToolTip(QStringLiteral("Keep the visible frame window attached to the playback frame."));
    top->addWidget(followBox_);
    fitAllBox_ = new QCheckBox(QStringLiteral("Fit collected"), container);
    fitAllBox_->setToolTip(QStringLiteral(
        "Off: a trailing time window scrolls with playback. "
        "On: all collected past frames are shown at once."));
    top->addWidget(fitAllBox_);
    residueModeBox_ = new QCheckBox(QStringLiteral("Residue"), container);
    residueModeBox_->setToolTip(QStringLiteral(
        "Use the focused atom's residue and plot phi, psi, omega, and chi dihedrals."));
    top->addWidget(residueModeBox_);
    auto* windowLabel = new QLabel(QStringLiteral("Window"), container);
    top->addWidget(windowLabel);
    windowFramesSpin_ = new QSpinBox(container);
    windowFramesSpin_->setRange(8, 1000000);
    windowFramesSpin_->setSingleStep(10);
    windowFramesSpin_->setAccelerated(true);
    windowFramesSpin_->setSuffix(QStringLiteral(" frames"));
    windowFramesSpin_->setToolTip(QStringLiteral("Number of frames visible in strip-chart follow mode."));
    top->addWidget(windowFramesSpin_);
    viewportReadout_ = new QLabel(QStringLiteral("f1-f1"), container);
    viewportReadout_->setToolTip(QStringLiteral("Visible frame range on the strip-chart x-axis."));
    top->addWidget(viewportReadout_);
    top->addStretch(1);
    vbox->addLayout(top);

    // ---- Custom strip stack: structural + DFT tracks share one time viewport.
    stackWidget_ = new StripStackWidget(container);
    vbox->addWidget(stackWidget_, 4);

    setWidget(container);

    ACONNECT(fitAllBox_.data(), &QCheckBox::toggled, this, &StripChartDock::onFitToggled);
    ACONNECT(stackWidget_.data(), &StripStackWidget::revealRequested,
             this, &StripChartDock::revealRequested);
    ACONNECT(residueModeBox_.data(), &QCheckBox::toggled, this, &StripChartDock::onResidueModeToggled);
    ACONNECT(followBox_.data(), &QCheckBox::toggled, this, [this](bool on) {
        if (on && fitAllBox_ && fitAllBox_->isChecked()) {
            const QSignalBlocker block(fitAllBox_);
            fitAllBox_->setChecked(false);
        }
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

void StripChartDock::buildTrack(Track& tr, QVBoxLayout* into, const QColor& traceColor, int stretch) {
    tr.chart = new QChart();
    tr.chart->legend()->hide();
    tr.chart->setMargins(QMargins(4, 4, 4, 4));

    tr.series = new QLineSeries(tr.chart);
    {
        QPen pen(traceColor);
        pen.setWidthF(2.0);
        tr.series->setPen(pen);
    }
    tr.chart->addSeries(tr.series);

    tr.cursor = new QLineSeries(tr.chart);
    {
        QPen pen(kCursorColor);
        pen.setWidthF(1.5);
        pen.setStyle(Qt::DashLine);
        tr.cursor->setPen(pen);
    }
    tr.chart->addSeries(tr.cursor);

    tr.xAxis = new QValueAxis(tr.chart);
    tr.xAxis->setTitleText(QStringLiteral("frame"));
    tr.xAxis->setLabelFormat(QStringLiteral("%d"));
    tr.xAxis->setGridLineVisible(true);
    tr.xAxis->setMinorTickCount(1);
    tr.xAxis->setMinorGridLineVisible(true);
    tr.chart->addAxis(tr.xAxis, Qt::AlignBottom);

    tr.yAxis = new QValueAxis(tr.chart);
    tr.yAxis->setTitleText(QStringLiteral("value"));
    tr.yAxis->setGridLineVisible(true);
    tr.yAxis->setMinorTickCount(1);
    tr.yAxis->setMinorGridLineVisible(true);
    tr.chart->addAxis(tr.yAxis, Qt::AlignLeft);

    tr.series->attachAxis(tr.xAxis);
    tr.series->attachAxis(tr.yAxis);
    tr.cursor->attachAxis(tr.xAxis);
    tr.cursor->attachAxis(tr.yAxis);

    tr.view = new QChartView(tr.chart);  // takes ownership of tr.chart
    tr.view->setRenderHint(QPainter::Antialiasing);
    into->addWidget(tr.view, stretch);

    tr.readout = new QLabel(QStringLiteral("—"));
    tr.readout->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    into->addWidget(tr.readout);
}

void StripChartDock::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    ASSERT_THREAD(this);
    protein_      = protein;
    conformation_ = conformation;
}

void StripChartDock::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    // Re-bind guard: ACONNECT cannot pass Qt::UniqueConnection, so drop any
    // prior selection's connections to this dock before rebinding (architecture
    // skill §4 — re-init paths must not stack duplicate connections).
    if (selection_)
        disconnect(selection_, nullptr, this, nullptr);
    selection_ = selection;
    if (selection) {
        // `changed` fires on every membership change, and every focus change
        // accompanies a membership change under the current gesture policy
        // (plain pick replaces, Shift+pick toggles) — so one rebuild on `changed`
        // refreshes BOTH the structural channel (from atoms) and the DFT channel
        // (from focus). No separate focusChanged wire needed.
        ACONNECT(selection, &model::AtomSelection::changed, this, &StripChartDock::onSelectionChanged);
        ACONNECT(selection, &model::AtomSelection::cleared, this, &StripChartDock::clearSelection);
    }
    rebuildChannels();
}

void StripChartDock::setDftStore(model::DftShieldingStore* store) {
    ASSERT_THREAD(this);
    dftStore_ = store;
    const bool have = (store != nullptr);
    if (dft_.view) dft_.view->setVisible(have);
    if (dft_.readout) dft_.readout->setVisible(have);
    rebuildChannels();
    refreshStackTracks();
}

void StripChartDock::setTimeViewport(TimeViewportController* viewport) {
    ASSERT_THREAD(this);
    if (timeViewport_)
        disconnect(timeViewport_, nullptr, this, nullptr);
    timeViewport_ = viewport;
    if (timeViewport_) {
        ACONNECT(timeViewport_.data(),
                 &TimeViewportController::visibleRangeChanged,
                 this,
                 &StripChartDock::onVisibleRangeChanged);
        if (windowFramesSpin_) {
            windowFramesSpin_->setRange(1, timeViewport_->frameCount());
            {
                const QSignalBlocker block(windowFramesSpin_);
                windowFramesSpin_->setValue(timeViewport_->windowFrames());
            }
        }
        ACONNECT(timeViewport_.data(), &TimeViewportController::windowFramesChanged, this, [this](int frames) {
            if (!windowFramesSpin_)
                return;
            const QSignalBlocker block(windowFramesSpin_);
            windowFramesSpin_->setValue(frames);
        });
        ACONNECT(timeViewport_.data(), &TimeViewportController::followPlayheadChanged, this, [this](bool on) {
            if (followBox_) {
                const QSignalBlocker block(followBox_);
                followBox_->setChecked(on);
            }
            if (on && fitAllBox_) {
                const QSignalBlocker block(fitAllBox_);
                fitAllBox_->setChecked(false);
            }
            if (windowFramesSpin_)
                windowFramesSpin_->setEnabled(on);
        });
        if (followBox_) {
            const QSignalBlocker block(followBox_);
            followBox_->setChecked(timeViewport_->followsPlayhead());
        }
        updateViewportReadout(timeViewport_->visibleStart(), timeViewport_->visibleEnd());
    }
    if (stackWidget_)
        stackWidget_->setTimeViewport(timeViewport_);
    renderTrack(structural_, frame_);
    renderTrack(dft_, frame_);
    refreshStackTracks();
}

void StripChartDock::setFrame(int t) {
    ASSERT_THREAD(this);
    frame_ = t;
    extendTo(t);                  // forward motion appends; backward is a no-op
    renderTrack(structural_, t);  // window + cursor + readout (cheap)
    renderTrack(dft_, t);
    if (dftFftStrip_)
        refreshStackTracks();
    else if (stackWidget_)
        stackWidget_->setCurrentFrame(t);
}

void StripChartDock::onSelectionChanged() {
    ASSERT_THREAD(this);
    rebuildChannels();
}

void StripChartDock::clearSelection() {
    ASSERT_THREAD(this);
    for (Track* tr : {&structural_, &dft_}) {
        clearTrack(*tr);
    }
    clearDftSignalStrips();
    clearResidueDihedralChannels();
    if (fftSeries_) fftSeries_->clear();
    if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
    fftPoints_.clear();
    fftReadoutText_ = QStringLiteral("—");
    dftFftPoints_.clear();
    dftFftReadoutText_ = QStringLiteral("—");
    fftRecomputeAtSize_ = 0;
    refreshStackTracks();
}

void StripChartDock::onFitToggled(bool on) {
    ASSERT_THREAD(this);
    if (windowFramesSpin_)
        windowFramesSpin_->setEnabled(!on);
    if (followBox_) {
        const QSignalBlocker block(followBox_);
        followBox_->setChecked(!on);
    }
    if (timeViewport_) {
        if (on)
            timeViewport_->fitCollectedRange();
        else
            timeViewport_->setFollowPlayhead(true);
    } else {
        fitAll_ = on;
    }
    renderTrack(structural_, frame_);
    renderTrack(dft_, frame_);
    if (stackWidget_)
        stackWidget_->update();
}

void StripChartDock::onResidueModeToggled(bool on) {
    ASSERT_THREAD(this);
    if (residueMode_ == on)
        return;
    residueMode_ = on;
    rebuildChannels();
}

void StripChartDock::onVisibleRangeChanged(int, int) {
    ASSERT_THREAD(this);
    if (timeViewport_)
        updateViewportReadout(timeViewport_->visibleStart(), timeViewport_->visibleEnd());
    renderTrack(structural_, frame_);
    renderTrack(dft_, frame_);
    if (stackWidget_)
        stackWidget_->update();
}

void StripChartDock::updateViewportReadout(int first, int last) {
    if (!viewportReadout_)
        return;
    viewportReadout_->setText(QStringLiteral("f%1-f%2").arg(first + 1).arg(last + 1));
}

std::optional<model::SignalBinding> StripChartDock::currentDftBinding() const {
    if (!dftStore_ || !selection_ || !selection_->hasFocus())
        return std::nullopt;
    model::SignalBinding binding;
    binding.key = model::SignalKey{
        model::SignalFamily::DftShielding,
        QStringLiteral("sigma"),
        QStringLiteral("total.T0"),
    };
    binding.anchorKind = model::SignalAnchorKind::Atom;
    binding.atom = selection_->focus();
    binding.followsFocus = true;
    return binding;
}

std::optional<model::SignalBinding> StripChartDock::sourceBindingForFft() const {
    if (structural_.active && structural_.binding)
        return structural_.binding;
    for (const Track& tr : residueDihedrals_) {
        if (tr.active && tr.binding)
            return tr.binding;
    }
    return std::nullopt;
}

void StripChartDock::clearTrack(Track& tr, const QString& readout) {
    tr.active = false;
    tr.binding.reset();
    tr.buffer.clear();
    tr.source.sample = nullptr;
    if (tr.series) tr.series->clear();
    if (tr.cursor) tr.cursor->clear();
    if (tr.readout) tr.readout->setText(readout);
}

void StripChartDock::clearResidueDihedralChannels() {
    for (Track& tr : residueDihedrals_)
        clearTrack(tr);
    residueDihedrals_.clear();
}

void StripChartDock::clearDftSignalStrips() {
    dftTimeStrip_.reset();
    dftFftStrip_.reset();
    dftFftPoints_.clear();
    dftFftReadoutText_ = QStringLiteral("—");
}

QString StripChartDock::residueDisplayLabel(std::size_t residueIdx) const {
    if (!protein_ || residueIdx >= protein_->residueCount())
        return QStringLiteral("residue %1").arg(residueIdx + 1);
    const auto& res = protein_->residue(residueIdx);
    const QString res3 = protein_->residueLabel(residueIdx,
                                                model::NamingConvention::Amber,
                                                model::NamingSource::Verbatim);
    const QString chain = res.address.chainId.isEmpty()
                              ? QString()
                              : QStringLiteral("%1:").arg(res.address.chainId);
    return QStringLiteral("%1%2%3").arg(chain, res3).arg(res.address.residueNumber);
}

QString StripChartDock::atomDisplayLabel(std::size_t atomIdx) const {
    if (!protein_ || atomIdx >= protein_->atomCount())
        return QStringLiteral("#%1").arg(atomIdx);
    const auto& atom = protein_->atom(atomIdx);
    const QString atomName = protein_->atomNames(atomIdx).amber;
    if (atom.residueIndex >= 0
        && static_cast<std::size_t>(atom.residueIndex) < protein_->residueCount()) {
        return QStringLiteral("%1:%2").arg(
            residueDisplayLabel(static_cast<std::size_t>(atom.residueIndex)),
            atomName);
    }
    return QStringLiteral("#%1:%2").arg(atomIdx).arg(atomName);
}

QString StripChartDock::selectionTupleLabel(const std::vector<std::size_t>& atoms) const {
    QStringList labels;
    labels.reserve(static_cast<int>(atoms.size()));
    for (std::size_t atom : atoms)
        labels.push_back(atomDisplayLabel(atom));
    return labels.join(QStringLiteral(" - "));
}

void StripChartDock::rebuildChannels() {
    ASSERT_THREAD(this);

    // Reset the active chart source. Normal mode shows an ordered atom-tuple
    // measurement; residue mode shows all available dihedrals for the focused
    // atom's residue. DFT stays with normal atom mode for now.
    clearTrack(structural_);
    if (residueMode_)
        clearTrack(dft_);
    clearDftSignalStrips();
    clearResidueDihedralChannels();
    fftRecomputeAtSize_       = 0;
    fftPoints_.clear();
    fftReadoutText_ = QStringLiteral("—");
    dftFftPoints_.clear();
    dftFftReadoutText_ = QStringLiteral("—");

    model::Conformation* conf     = conformation_;  // QPointer -> raw
    const bool           haveTraj = protein_ && conf && conf->frameCount() > 1;

    if (residueMode_) {
        bindResidueDihedralChannels();
        extendTo(frame_);
        rebuildFft();
        refreshStackTracks();
        return;
    }

    // ----- Structural channel: a >=2-atom geometry over a trajectory. -----
    QString structMsg;
    if (!haveTraj) {
        // Multi-frame only (gate on frameCount, NOT asTrajectory, so any future
        // multi-frame backing engages): a single pose has no frame axis.
        structMsg = (conf && conf->frameCount() <= 1) ? QStringLiteral("single pose — no trajectory")
                                                      : QStringLiteral("—");
    } else if (!selection_ || selection_->geometryKind() == model::GeometryKind::None) {
        structMsg = (selection_ && !selection_->empty()) ? QStringLiteral("select ≥ 2 atoms")
                                                         : QStringLiteral("—");
    } else {
        // Build the structural channel for THIS selection tuple. The source
        // captures a copy of the ordered atoms and reads positions through the
        // conformation's atomPosition seam via model::Measure — works for any
        // conformation backing, not just the trajectory H5.
        const model::GeometryKind kind   = selection_->geometryKind();
        const bool                isDist = (kind == model::GeometryKind::Distance);
        structural_.buffer.id    = QStringLiteral("geometry");
        structural_.buffer.label = QStringLiteral("%1: %2")
                                       .arg(QString::fromLatin1(model::NameForGeometryKind(kind)),
                                            selectionTupleLabel(selection_->atoms()));
        structural_.buffer.unit  = isDist ? QString::fromUtf8("Å") : QString::fromUtf8("°");

        const std::vector<std::size_t> atoms = selection_->atoms();
        structural_.source.id    = structural_.buffer.id;
        structural_.source.label = structural_.buffer.label;
        structural_.source.unit  = structural_.buffer.unit;
        structural_.source.sample = [this, atoms](std::size_t frame) -> std::optional<double> {
            model::Conformation* c = conformation_;  // QPointer null-checked each call
            if (!c) return std::nullopt;
            const auto m = model::Measure(*c, frame, atoms);
            return m.valid ? std::optional<double>(m.value) : std::nullopt;
        };
        model::SignalBinding binding;
        binding.key.family = model::SignalFamily::Geometry;
        switch (kind) {
            case model::GeometryKind::Distance:
                binding.key.name = QStringLiteral("distance");
                break;
            case model::GeometryKind::Angle:
                binding.key.name = QStringLiteral("angle");
                break;
            case model::GeometryKind::Dihedral:
                binding.key.name = QStringLiteral("dihedral");
                break;
            case model::GeometryKind::None:
                break;
        }
        binding.anchorKind = model::SignalAnchorKind::AtomTuple;
        binding.atomTuple = atoms;
        binding.followsFocus = true;
        structural_.binding = std::move(binding);
        structural_.active = true;
        if (structural_.yAxis) {
            structural_.yAxis->setTitleText(structural_.buffer.label.toLower() +
                                            QStringLiteral(" (") + structural_.buffer.unit + QStringLiteral(")"));
        }
        qCDebug(cStrip).noquote() << "structural channel |" << structural_.buffer.label
                                  << "| frames=" << conf->frameCount() << "| backfill to" << frame_;
    }
    if (!structural_.active) {
        if (structural_.series) structural_.series->clear();
        if (structural_.cursor) structural_.cursor->clear();
        if (structural_.readout) structural_.readout->setText(structMsg);
        if (fftSeries_) fftSeries_->clear();
        if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
        fftPoints_.clear();
        fftReadoutText_ = QStringLiteral("—");
    }

    // ----- DFT channel: the focus atom's σ over the trajectory (needs a store).
    bindDftChannel();

    // ----- Fill 0..currentFrame (real past, never the future) + render. -----
    extendTo(frame_);
    renderTrack(structural_, frame_);
    renderTrack(dft_, frame_);
    refreshStackTracks();
}

void StripChartDock::bindDftChannel() {
    // Reset, then (re)bind to the current FOCUS atom. The DFT channel is
    // INDEPENDENT of the structural one: a single picked atom defines no geometry
    // (structural blank) but still has a focus whose shielding we can chart.
    clearTrack(dft_);
    clearDftSignalStrips();

    model::Conformation* conf     = conformation_;
    const bool           haveTraj = protein_ && conf && conf->frameCount() > 1;
    if (!dftStore_ || !haveTraj || !selection_ || !selection_->hasFocus()) {
        if (dft_.series) dft_.series->clear();
        if (dft_.cursor) dft_.cursor->clear();
        if (dft_.readout)
            dft_.readout->setText(dftStore_ ? QStringLiteral("select an atom") : QStringLiteral("—"));
        return;
    }

    const std::size_t atom = selection_->focus();
    const auto maybeBinding = currentDftBinding();
    if (!maybeBinding)
        return;
    const model::SignalBinding binding = *maybeBinding;

    dftTimeStrip_ = std::make_unique<model::DftSigmaAtomTimeStrip>(conf, dftStore_);
    if (dftTimeStrip_->canBind(binding)) {
        dftTimeStrip_->bind(binding);
    } else {
        dftTimeStrip_.reset();
    }

    dftFftStrip_ = std::make_unique<model::DftSigmaAtomFftStrip>(conf, dftStore_);
    if (dftFftStrip_->canBind(binding)) {
        dftFftStrip_->bind(binding);
    } else {
        dftFftStrip_.reset();
    }

    qCDebug(cStrip).noquote() << "dft channel | focus atom" << atom << "| jobs=" << dftStore_->jobCount();
}

void StripChartDock::bindResidueDihedralChannels() {
    if (!protein_ || !conformation_ || !selection_ || !selection_->hasFocus())
        return;

    const auto* traj = conformation_->asTrajectory();
    const auto* ts = (traj && traj->h5()) ? traj->h5()->dihedrals() : nullptr;
    if (!ts || ts->n_residues == 0 || ts->n_frames == 0)
        return;

    const std::size_t focusAtom = selection_->focus();
    if (focusAtom >= protein_->atomCount())
        return;
    const int resIdx = protein_->atom(focusAtom).residueIndex;
    if (resIdx < 0 || static_cast<std::size_t>(resIdx) >= protein_->residueCount()
        || static_cast<std::size_t>(resIdx) >= ts->n_residues)
        return;

    const std::size_t r = static_cast<std::size_t>(resIdx);
    const QString residue = residueDisplayLabel(r);
    const QString units = ts->units;

    auto addAngle = [this, ts, r, residue, units](QString id,
                                                  QString shortLabel,
                                                  auto sampleValue) {
        Track tr;
        tr.active = true;
        tr.buffer.id = std::move(id);
        tr.buffer.label = QStringLiteral("%1 %2").arg(residue, shortLabel);
        tr.buffer.unit = QString::fromUtf8("°");
        model::SignalBinding binding;
        binding.key = model::SignalKey{
            model::SignalFamily::ResidueDihedral,
            shortLabel,
            QString(),
        };
        binding.anchorKind = model::SignalAnchorKind::Residue;
        binding.residue = r;
        binding.followsFocus = true;
        tr.binding = std::move(binding);
        tr.source.id = tr.buffer.id;
        tr.source.label = tr.buffer.label;
        tr.source.unit = tr.buffer.unit;
        tr.source.sample = [ts, r, units, sampleValue](std::size_t frame) -> std::optional<double> {
            if (frame >= ts->n_frames || !ts->sourceAttachedAt(frame))
                return std::nullopt;
            const double value = sampleValue(*ts, r, frame);
            if (!std::isfinite(value))
                return std::nullopt;
            return angleToDegrees(value, units);
        };
        residueDihedrals_.push_back(std::move(tr));
    };

    addAngle(QStringLiteral("residue.phi"), QStringLiteral("phi"),
             [](const model::QtDihedralTimeSeries& series, std::size_t rr, std::size_t frame) {
                 return series.phiAt(rr, frame);
             });
    addAngle(QStringLiteral("residue.psi"), QStringLiteral("psi"),
             [](const model::QtDihedralTimeSeries& series, std::size_t rr, std::size_t frame) {
                 return series.psiAt(rr, frame);
             });
    addAngle(QStringLiteral("residue.omega"), QStringLiteral("omega"),
             [](const model::QtDihedralTimeSeries& series, std::size_t rr, std::size_t frame) {
                 return series.omegaAt(rr, frame);
             });
    for (int chi = 0; chi < 4; ++chi) {
        const std::size_t existsIdx = r * 4 + static_cast<std::size_t>(chi);
        const bool haveMask = existsIdx < ts->chi_exists.size();
        const bool haveValues = !ts->chi.empty();
        if (!haveValues || (haveMask && ts->chi_exists[existsIdx] == 0))
            continue;
        addAngle(QStringLiteral("residue.chi%1").arg(chi + 1),
                 QStringLiteral("chi%1").arg(chi + 1),
                 [chi](const model::QtDihedralTimeSeries& series, std::size_t rr, std::size_t frame) {
                     return series.chiAt(rr, frame, chi);
                 });
    }
}

bool StripChartDock::extendDenseTrackTo(Track& tr, int t) {
    if (t < 0 || !tr.active || !tr.source.sample)
        return false;
    bool extended = false;
    for (long long f = tr.buffer.lastFrame() + 1; f <= t; ++f) {
        tr.buffer.append(tr.source.sample(static_cast<std::size_t>(f)));
        extended = true;
    }
    return extended;
}

void StripChartDock::extendTo(int t) {
    if (t < 0) return;

    // Structural/residue channels are dense and cheap: every frame in
    // (lastFrame, t] is sampled from resident H5 positions or per-residue TRs.
    // Forward play / a forward scrub appends; backward motion only moves the
    // cursor. No frame > t is ever sampled ("you cannot see the future").
    bool structExt = extendDenseTrackTo(structural_, t);
    bool residueExt = false;
    for (Track& tr : residueDihedrals_)
        residueExt = extendDenseTrackTo(tr, t) || residueExt;

    // DFT strips preserve the sparse-load policy: on forward jumps they append
    // gaps for skipped frames and parse only the landing frame; continuous play
    // fills one frame at a time. The time strip and FFT strip both bind to the
    // same DFT signal and sample the store independently.
    if (dftTimeStrip_)
        dftTimeStrip_->extendToFrame(static_cast<std::size_t>(t));
    if (dftFftStrip_)
        dftFftStrip_->extendToFrame(static_cast<std::size_t>(t));

    const Track* fftTrack = structural_.active ? &structural_
                                               : (!residueDihedrals_.empty() ? &residueDihedrals_.front() : nullptr);
    if ((structExt || residueExt) && fftTrack
        && static_cast<int>(fftTrack->buffer.size()) >= fftRecomputeAtSize_) {
        rebuildFft();
        fftRecomputeAtSize_ = static_cast<int>(fftTrack->buffer.size()) + kFftRecomputeStride;
    }
}

void StripChartDock::renderTrack(Track& tr, int t) {
    ASSERT_THREAD(this);
    if (!tr.active || !tr.xAxis || !tr.series) return;

    const long long last = tr.buffer.lastFrame();  // max frame reached (the past we hold)
    if (last < 0) {
        if (tr.series) tr.series->clear();
        if (tr.cursor) tr.cursor->clear();
        if (tr.readout) tr.readout->setText(QStringLiteral("—"));
        return;
    }

    // X window is shared application state when TimeViewportController is
    // bound. Fallback preserves the pre-controller local behavior.
    int       axisLo;
    int       axisHi;
    const int n = static_cast<int>(last) + 1;
    if (timeViewport_) {
        axisLo = timeViewport_->visibleStart();
        axisHi = timeViewport_->visibleEnd();
    } else if (fitAll_ || n <= kWindowFrames) {
        axisLo = 0;
        axisHi = static_cast<int>(last);
    } else {
        const int loMax = std::max(0, static_cast<int>(last) - kWindowFrames + 1);
        axisLo = std::clamp(t - kWindowFrames / 2, 0, loMax);
        axisHi = std::min(static_cast<int>(last), axisLo + kWindowFrames - 1);
    }
    if (axisLo > axisHi)
        std::swap(axisLo, axisHi);
    tr.xAxis->setRange(axisLo, axisHi);

    // Slice the visible window from the authoritative buffer. The BUFFER is
    // never decimated; only a huge DISPLAYED window strides (µs fit-all). A gap
    // (NaN) is kept so Qt Charts breaks the line rather than bridging a missing
    // value.
    const int      dataLo = std::max(0, axisLo);
    const int      dataHi = std::min(static_cast<int>(last), axisHi);
    const int      span   = std::max(1, dataHi - dataLo + 1);
    const int      stride = std::max(1, span / kMaxDisplayPoints);
    QList<QPointF> pts;
    pts.reserve(span / stride + 1);
    for (int f = dataLo; f <= dataHi; f += stride)
        pts.append(QPointF(f, tr.buffer.values[static_cast<std::size_t>(f)]));  // NaN at gaps
    tr.series->replace(pts);

    // Y range from the collected valid data (padded); neutral if all-gaps.
    double yMin = tr.buffer.yMin;
    double yMax = tr.buffer.yMax;
    if (!tr.buffer.hasRange) {
        yMin = 0.0;
        yMax = 1.0;
    } else if (yMax - yMin < 1e-9) {
        yMin -= 1.0;
        yMax += 1.0;
    } else {
        const double pad = 0.05 * (yMax - yMin);
        yMin -= pad;
        yMax += pad;
    }
    tr.yAxis->setRange(yMin, yMax);

    // Playhead at the current frame (clamped into the data we hold).
    const int tc = std::clamp(t, 0, static_cast<int>(last));
    if (tr.cursor)
        tr.cursor->replace(QList<QPointF>{QPointF(tc, yMin), QPointF(tc, yMax)});

    // Digital readout — the value at the current frame (the number we keep OUT
    // of the moving 3-D scene).
    if (tr.readout) {
        const bool haveVal = tc <= last && tr.buffer.valid[static_cast<std::size_t>(tc)];
        if (haveVal) {
            tr.readout->setText(QStringLiteral("%1: %2 %3   (frame %4 / %5)")
                                    .arg(tr.buffer.label)
                                    .arg(tr.buffer.values[static_cast<std::size_t>(tc)], 0, 'f',
                                         decimalsForUnit(tr.buffer.unit))
                                    .arg(tr.buffer.unit)
                                    .arg(tc + 1)
                                    .arg(n));
        } else {
            tr.readout->setText(QStringLiteral("%1: —   (frame %2 / %3)")
                                    .arg(tr.buffer.label)
                                    .arg(tc + 1)
                                    .arg(n));
        }
    }
}

void StripChartDock::refreshStackTracks() {
    if (!stackWidget_)
        return;

    QVector<StripStackWidget::Track> tracks;
    auto addTimeTrack = [&tracks](const model::ChannelBuffer* buffer,
                                  const QColor& color,
                                  const std::optional<model::SignalBinding>& binding) {
        StripStackWidget::Track item;
        item.buffer = buffer;
        item.color = color;
        if (binding) {
            item.hasBinding = true;
            item.binding = *binding;
        }
        tracks.push_back(std::move(item));
    };
    if (structural_.active) {
        addTimeTrack(&structural_.buffer, kGeomColor, structural_.binding);
    }
    const auto dftBinding = currentDftBinding();
    if (dftTimeStrip_) {
        addTimeTrack(&dftTimeStrip_->buffer(), kDftColor, dftBinding);
    } else if (dft_.active) {
        addTimeTrack(&dft_.buffer, kDftColor, dft_.binding);
    }
    const QColor residueColors[] = {kPhiColor, kPsiColor, kOmegaColor, kChiColor,
                                    kChiColor.lighter(125), kChiColor.darker(120), kChiColor.lighter(150)};
    for (std::size_t i = 0; i < residueDihedrals_.size(); ++i) {
        Track& tr = residueDihedrals_[i];
        if (!tr.active)
            continue;
        addTimeTrack(&tr.buffer,
                     residueColors[i % (sizeof(residueColors) / sizeof(residueColors[0]))],
                     tr.binding);
    }
    stackWidget_->setTracks(std::move(tracks));
    QVector<StripStackWidget::SpectrumTrack> spectrumTracks;
    auto addSpectrumTrack = [&spectrumTracks](const QVector<QPointF>* points,
                                              QString label,
                                              QString xUnit,
                                              QString yUnit,
                                              QString readout,
                                              const QColor& color,
                                              const std::optional<model::SignalBinding>& binding) {
        StripStackWidget::SpectrumTrack item;
        item.points = points;
        item.label = std::move(label);
        item.xUnit = std::move(xUnit);
        item.yUnit = std::move(yUnit);
        item.readout = std::move(readout);
        item.color = color;
        if (binding) {
            item.hasBinding = true;
            item.binding = *binding;
        }
        spectrumTracks.push_back(std::move(item));
    };
    if (!fftPoints_.empty()) {
        addSpectrumTrack(&fftPoints_,
                         QStringLiteral("Power spectrum"),
                         QStringLiteral("1/ns"),
                         QStringLiteral("power"),
                         fftReadoutText_,
                         kFftColor,
                         sourceBindingForFft());
    }
    if (dftFftStrip_) {
        const model::StripRenderData dftFft = dftFftStrip_->renderData();
        dftFftPoints_ = dftFft.points;
        dftFftReadoutText_ = dftFft.readout;
        if (!dftFftPoints_.empty()) {
            addSpectrumTrack(&dftFftPoints_,
                             dftFft.label,
                             dftFft.xUnit,
                             dftFft.yUnit,
                             dftFftReadoutText_,
                             kDftColor.darker(120),
                             dftBinding);
        }
    }
    stackWidget_->setSpectrumTracks(std::move(spectrumTracks));
    stackWidget_->setCurrentFrame(frame_);
}

void StripChartDock::rebuildFft() {
    ASSERT_THREAD(this);

    auto blankFft = [&]() {
        fftPoints_.clear();
        fftReadoutText_ = QStringLiteral("—");
        if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
        if (stackWidget_)
            stackWidget_->setSpectrumTracks(QVector<StripStackWidget::SpectrumTrack>{});
    };

    const Track* sourceTrack = structural_.active ? &structural_
                                                  : (!residueDihedrals_.empty() ? &residueDihedrals_.front() : nullptr);
    if (!sourceTrack || !sourceTrack->active) {
        blankFft();
        return;
    }

    // FFT over the structural track's collected VALID values (the geometry
    // channel is gapless in practice; gather defensively). A dihedral straddling
    // ±180° must read as one smooth signal, so phase-unwrap a copy before the
    // transform (the time-domain panel keeps the physical wrapped values).
    const auto&         buf = sourceTrack->buffer;
    std::vector<double> input;
    input.reserve(buf.values.size());
    for (std::size_t i = 0; i < buf.values.size(); ++i)
        if (buf.valid[i])
            input.push_back(buf.values[i]);

    const bool unwrapAngles = sourceTrack != &dft_ && sourceTrack->buffer.unit == QString::fromUtf8("°");
    if (unwrapAngles) {
        for (std::size_t k = 1; k < input.size(); ++k) {
            double d = input[k] - input[k - 1];
            while (d > 180.0)  { input[k] -= 360.0; d -= 360.0; }
            while (d < -180.0) { input[k] += 360.0; d += 360.0; }
        }
    }

    model::Conformation* conf  = conformation_;
    const double         dtPs  = (conf && conf->frameCount() >= 2)
                                     ? (conf->timePicoseconds(1) - conf->timePicoseconds(0))
                                     : 0.0;
    const model::PowerSpectrum spec = model::ComputePowerSpectrum(input, dtPs);

    if (!spec.valid || spec.power.size() < 2) {
        blankFft();
        return;
    }

    fftPoints_.clear();
    fftPoints_.reserve(static_cast<int>(spec.power.size()));
    for (std::size_t k = 1; k < spec.power.size(); ++k) {  // skip k=0 (DC, ≈0 after mean-subtract)
        fftPoints_.append(QPointF(spec.frequencyPerNs[k], spec.power[k]));
    }

    // dominantPeriodPs == 0 means ComputePowerSpectrum found no non-DC peak
    // (a flat / featureless series) — say so rather than print a bogus
    // record-length period (the old false-period bug).
    if (spec.dominantPeriodPs <= 0.0) {
        fftReadoutText_ = QStringLiteral("no dominant period");
    } else {
        const double p = spec.dominantPeriodPs;
        fftReadoutText_ = p >= 1000.0
                              ? QStringLiteral("%1 ns").arg(p / 1000.0, 0, 'f', 2)
                              : QStringLiteral("%1 ps").arg(p, 0, 'f', 0);
    }
    if (stackWidget_)
        refreshStackTracks();
}

}  // namespace h5reader::app
