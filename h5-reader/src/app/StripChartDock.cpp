#include "StripChartDock.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"
#include "../model/DftShieldingStore.h"
#include "../model/QtProtein.h"
#include "../model/SpectralAnalysis.h"

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
#include <QVBoxLayout>
#include <QWidget>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
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

const QColor kGeomColor(33, 102, 171);    // structural trace: saturated blue
const QColor kDftColor(204, 102, 0);      // shielding trace: burnt orange (Okabe-Ito)
const QColor kCursorColor(217, 26, 26);   // playhead: red dashed
const QColor kFftColor(0, 128, 96);       // spectrum: teal, distinct from above

int decimalsForUnit(const QString& unit) {
    if (unit == QString::fromUtf8("Å"))  return 3;
    if (unit == QStringLiteral("ppm"))   return 2;
    return 1;  // degrees
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

    // Top row: Fit-all toggle (overview vs scrolling window).
    auto* top  = new QHBoxLayout;
    fitAllBox_ = new QCheckBox(QStringLiteral("Fit all frames"), container);
    fitAllBox_->setToolTip(QStringLiteral(
        "Off: a fixed-width window scrolls with the playhead. "
        "On: the whole collected trace is shown at once (overview)."));
    top->addWidget(fitAllBox_);
    top->addStretch(1);
    vbox->addLayout(top);

    // ---- Structural panel: the selection's geometry observable vs frame. ----
    buildTrack(structural_, vbox, kGeomColor, /*stretch=*/2);

    // ---- Shielding panel: the focus atom's DFT σ vs frame. Hidden until a DFT
    //      store is bound (many runs have no dft/ campaign) — see setDftStore.
    buildTrack(dft_, vbox, kDftColor, /*stretch=*/2);
    dft_.view->setVisible(false);
    dft_.readout->setVisible(false);

    // ---- Frequency-domain panel: the power spectrum of the structural series.
    fftChart_ = new QChart();
    fftChart_->legend()->hide();
    fftChart_->setMargins(QMargins(4, 4, 4, 4));

    fftSeries_ = new QLineSeries(fftChart_);
    {
        QPen pen(kFftColor);
        pen.setWidthF(1.5);
        fftSeries_->setPen(pen);
    }
    fftChart_->addSeries(fftSeries_);

    fftXAxis_ = new QValueAxis(fftChart_);
    fftXAxis_->setTitleText(QStringLiteral("frequency (1/ns)"));
    fftXAxis_->setGridLineVisible(true);
    fftXAxis_->setMinorTickCount(1);
    fftXAxis_->setMinorGridLineVisible(true);
    fftChart_->addAxis(fftXAxis_, Qt::AlignBottom);

    fftYAxis_ = new QValueAxis(fftChart_);
    fftYAxis_->setTitleText(QStringLiteral("power"));
    fftYAxis_->setGridLineVisible(true);
    fftYAxis_->setLabelFormat(QStringLiteral("%.2g"));
    fftChart_->addAxis(fftYAxis_, Qt::AlignLeft);

    fftSeries_->attachAxis(fftXAxis_);
    fftSeries_->attachAxis(fftYAxis_);

    fftView_ = new QChartView(fftChart_);  // takes ownership of fftChart_
    fftView_->setRenderHint(QPainter::Antialiasing);
    vbox->addWidget(fftView_, 1);

    fftReadout_ = new QLabel(QStringLiteral("—"), container);
    fftReadout_->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    vbox->addWidget(fftReadout_);

    setWidget(container);

    ACONNECT(fitAllBox_.data(), &QCheckBox::toggled, this, &StripChartDock::onFitToggled);
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
}

void StripChartDock::setFrame(int t) {
    ASSERT_THREAD(this);
    frame_ = t;
    extendTo(t);                  // forward motion appends; backward is a no-op
    renderTrack(structural_, t);  // window + cursor + readout (cheap)
    renderTrack(dft_, t);
}

void StripChartDock::onSelectionChanged() {
    ASSERT_THREAD(this);
    rebuildChannels();
}

void StripChartDock::clearSelection() {
    ASSERT_THREAD(this);
    for (Track* tr : {&structural_, &dft_}) {
        tr->active        = false;
        tr->buffer.clear();
        tr->source.sample = nullptr;
        if (tr->series) tr->series->clear();
        if (tr->cursor) tr->cursor->clear();
        if (tr->readout) tr->readout->setText(QStringLiteral("—"));
    }
    if (fftSeries_) fftSeries_->clear();
    if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
    fftRecomputeAtSize_ = 0;
}

void StripChartDock::onFitToggled(bool on) {
    ASSERT_THREAD(this);
    fitAll_ = on;
    renderTrack(structural_, frame_);
    renderTrack(dft_, frame_);
}

void StripChartDock::rebuildChannels() {
    ASSERT_THREAD(this);

    // Reset the structural track (the DFT track is reset inside bindDftChannel).
    structural_.active        = false;
    structural_.buffer.clear();
    structural_.source.sample = nullptr;
    fftRecomputeAtSize_       = 0;

    model::Conformation* conf     = conformation_;  // QPointer -> raw
    const bool           haveTraj = protein_ && conf && conf->frameCount() > 1;

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
        structural_.buffer.label = QString::fromLatin1(model::NameForGeometryKind(kind));
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
        structural_.active = true;
        structural_.yAxis->setTitleText(structural_.buffer.label.toLower() +
                                        QStringLiteral(" (") + structural_.buffer.unit + QStringLiteral(")"));
        qCDebug(cStrip).noquote() << "structural channel |" << structural_.buffer.label
                                  << "| frames=" << conf->frameCount() << "| backfill to" << frame_;
    }
    if (!structural_.active) {
        if (structural_.series) structural_.series->clear();
        if (structural_.cursor) structural_.cursor->clear();
        if (structural_.readout) structural_.readout->setText(structMsg);
        if (fftSeries_) fftSeries_->clear();
        if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
    }

    // ----- DFT channel: the focus atom's σ over the trajectory (needs a store).
    bindDftChannel();

    // ----- Fill 0..currentFrame (real past, never the future) + render. -----
    extendTo(frame_);
    renderTrack(structural_, frame_);
    renderTrack(dft_, frame_);
}

void StripChartDock::bindDftChannel() {
    // Reset, then (re)bind to the current FOCUS atom. The DFT channel is
    // INDEPENDENT of the structural one: a single picked atom defines no geometry
    // (structural blank) but still has a focus whose shielding we can chart.
    dft_.active        = false;
    dft_.buffer.clear();
    dft_.source.sample = nullptr;

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
    dft_.buffer.id    = QStringLiteral("dft.total.T0");
    dft_.buffer.label = QStringLiteral("DFT σ (total)");
    dft_.buffer.unit  = QStringLiteral("ppm");
    dft_.source.id    = dft_.buffer.id;
    dft_.source.label = dft_.buffer.label;
    dft_.source.unit  = dft_.buffer.unit;
    dft_.source.sample = [this, atom](std::size_t frame) -> std::optional<double> {
        if (!dftStore_ || !conformation_) return std::nullopt;
        const std::size_t orig = conformation_->originalFrameIndex(frame);
        // Cached-or-null: extendTo() does the requestFrame() that actually parses.
        return dftStore_->sample(orig, atom, model::DftPart::Total, model::DftScalar::IsotropicT0);
    };
    dft_.active = true;
    dft_.yAxis->setTitleText(QStringLiteral("σ total (ppm)"));
    qCDebug(cStrip).noquote() << "dft channel | focus atom" << atom << "| jobs=" << dftStore_->jobCount();
}

void StripChartDock::extendTo(int t) {
    if (t < 0) return;

    // Structural: every frame in (lastFrame, t] is cheap (positions are in the
    // eager H5), so backfill them all. Forward play / a forward scrub appends;
    // a backward scrub leaves lastFrame ahead of t so the loop is a no-op and the
    // past remains. No frame > t is ever sampled ("you cannot see the future").
    bool structExt = false;
    if (structural_.active && structural_.source.sample) {
        for (long long f = structural_.buffer.lastFrame() + 1; f <= t; ++f) {
            structural_.buffer.append(structural_.source.sample(static_cast<std::size_t>(f)));
            structExt = true;
        }
    }

    // DFT: parsing a frame's ORCA .out is expensive, so do NOT parse the skipped
    // past on a forward jump — append GAPS for (lastFrame, t-1] and request+sample
    // only the LANDING frame t. Continuous play (one-frame steps) thus fills every
    // frame; a scrub leaves gaps that fill if revisited. v1 parses synchronously
    // here; frameReady() is the seam for a future prefetch worker. The two buffers
    // stay frame-aligned (each grows to t+1) but are read independently.
    if (dft_.active && dft_.source.sample && dftStore_) {
        const long long from = dft_.buffer.lastFrame() + 1;
        for (long long f = from; f < t; ++f)
            dft_.buffer.append(std::nullopt);  // skipped past -> gap (left unparsed)
        if (from <= t) {
            const std::size_t orig =
                conformation_ ? conformation_->originalFrameIndex(static_cast<std::size_t>(t))
                              : static_cast<std::size_t>(t);
            dftStore_->requestFrame(orig);  // parse the landing frame (synchronous v1)
            dft_.buffer.append(dft_.source.sample(static_cast<std::size_t>(t)));
        }
    }

    if (structExt && static_cast<int>(structural_.buffer.size()) >= fftRecomputeAtSize_) {
        rebuildFft();
        fftRecomputeAtSize_ = static_cast<int>(structural_.buffer.size()) + kFftRecomputeStride;
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

    // X window over 0..maxReached (NOT 0..t): once played, the trace stays and
    // the cursor marks the present. Full range when fitting or short, else a
    // fixed-width window scrolling with the playhead (the anti-scrunch).
    int       lo;
    int       hi;
    const int n = static_cast<int>(last) + 1;
    if (fitAll_ || n <= kWindowFrames) {
        lo = 0;
        hi = static_cast<int>(last);
    } else {
        const int loMax = std::max(0, static_cast<int>(last) - kWindowFrames + 1);
        lo = std::clamp(t - kWindowFrames / 2, 0, loMax);
        hi = std::min(static_cast<int>(last), lo + kWindowFrames - 1);
    }
    tr.xAxis->setRange(lo, hi);

    // Slice the visible window from the authoritative buffer. The BUFFER is
    // never decimated; only a huge DISPLAYED window strides (µs fit-all). A gap
    // (NaN) is kept so Qt Charts breaks the line rather than bridging a missing
    // value.
    const int      span   = hi - lo + 1;
    const int      stride = std::max(1, span / kMaxDisplayPoints);
    QList<QPointF> pts;
    pts.reserve(span / stride + 1);
    for (int f = lo; f <= hi; f += stride)
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

void StripChartDock::rebuildFft() {
    ASSERT_THREAD(this);
    if (!fftSeries_) return;

    auto blankFft = [&]() {
        fftSeries_->clear();
        if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
    };

    if (!structural_.active || !selection_) {
        blankFft();
        return;
    }

    // FFT over the structural track's collected VALID values (the geometry
    // channel is gapless in practice; gather defensively). A dihedral straddling
    // ±180° must read as one smooth signal, so phase-unwrap a copy before the
    // transform (the time-domain panel keeps the physical wrapped values).
    const auto&         buf = structural_.buffer;
    std::vector<double> input;
    input.reserve(buf.values.size());
    for (std::size_t i = 0; i < buf.values.size(); ++i)
        if (buf.valid[i])
            input.push_back(buf.values[i]);

    if (selection_->geometryKind() == model::GeometryKind::Dihedral) {
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

    QList<QPointF> fftPts;
    fftPts.reserve(static_cast<int>(spec.power.size()));
    double maxP = 0.0;
    for (std::size_t k = 1; k < spec.power.size(); ++k) {  // skip k=0 (DC, ≈0 after mean-subtract)
        fftPts.append(QPointF(spec.frequencyPerNs[k], spec.power[k]));
        maxP = std::max(maxP, spec.power[k]);
    }
    fftSeries_->replace(fftPts);

    const double fMax = spec.frequencyPerNs.back();
    fftXAxis_->setRange(0.0, fMax > 0.0 ? fMax : 1.0);
    fftYAxis_->setRange(0.0, maxP > 0.0 ? maxP * 1.05 : 1.0);

    if (fftReadout_) {
        // dominantPeriodPs == 0 means ComputePowerSpectrum found no non-DC peak
        // (a flat / featureless series) — say so rather than print a bogus
        // record-length period (the old false-period bug).
        if (spec.dominantPeriodPs <= 0.0) {
            fftReadout_->setText(QStringLiteral("no dominant period"));
        } else {
            const double p = spec.dominantPeriodPs;
            fftReadout_->setText(p >= 1000.0
                                     ? QStringLiteral("dominant period: %1 ns").arg(p / 1000.0, 0, 'f', 2)
                                     : QStringLiteral("dominant period: %1 ps").arg(p, 0, 'f', 0));
        }
    }
}

}  // namespace h5reader::app
