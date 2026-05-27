#include "StripChartDock.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"
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
#include <limits>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cStrip, "h5reader.stripchart")

// Scrolling time window (frames). When the trajectory is longer than this, the
// x-axis shows a fixed-width window that follows the playhead instead of
// cramming all frames into the panel ("scrunch as you go"). "Fit all" overrides.
constexpr int kWindowFrames = 300;

const QColor kGeomColor(33, 102, 171);    // time-domain trace: saturated blue
const QColor kCursorColor(217, 26, 26);   // playhead: red dashed
const QColor kFftColor(0, 128, 96);       // spectrum: teal, distinct from above
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
        "On: the whole trajectory is shown at once (overview)."));
    top->addWidget(fitAllBox_);
    top->addStretch(1);
    vbox->addLayout(top);

    // ---- Time-domain panel: the geometry observable vs frame. -------------
    chart_ = new QChart();
    chart_->legend()->hide();
    chart_->setMargins(QMargins(4, 4, 4, 4));

    geomSeries_ = new QLineSeries(chart_);
    {
        QPen pen(kGeomColor);
        pen.setWidthF(2.0);
        geomSeries_->setPen(pen);
    }
    chart_->addSeries(geomSeries_);

    cursorSeries_ = new QLineSeries(chart_);
    {
        QPen pen(kCursorColor);
        pen.setWidthF(1.5);
        pen.setStyle(Qt::DashLine);
        cursorSeries_->setPen(pen);
    }
    chart_->addSeries(cursorSeries_);

    // Science grids: major + minor on both axes.
    xAxis_ = new QValueAxis(chart_);
    xAxis_->setTitleText(QStringLiteral("frame"));
    xAxis_->setLabelFormat(QStringLiteral("%d"));
    xAxis_->setGridLineVisible(true);
    xAxis_->setMinorTickCount(1);
    xAxis_->setMinorGridLineVisible(true);
    chart_->addAxis(xAxis_, Qt::AlignBottom);

    yAxis_ = new QValueAxis(chart_);
    yAxis_->setTitleText(QStringLiteral("value"));
    yAxis_->setGridLineVisible(true);
    yAxis_->setMinorTickCount(1);
    yAxis_->setMinorGridLineVisible(true);
    chart_->addAxis(yAxis_, Qt::AlignLeft);

    geomSeries_->attachAxis(xAxis_);
    geomSeries_->attachAxis(yAxis_);
    cursorSeries_->attachAxis(xAxis_);
    cursorSeries_->attachAxis(yAxis_);

    chartView_ = new QChartView(chart_, container);  // takes ownership of chart_
    chartView_->setRenderHint(QPainter::Antialiasing);
    vbox->addWidget(chartView_, 2);

    // Time-domain digital readout — the measured value at the current frame,
    // the number we deliberately keep OUT of the (moving) 3-D scene.
    readout_ = new QLabel(QStringLiteral("—"), container);
    readout_->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    vbox->addWidget(readout_);

    // ---- Frequency-domain panel: the power spectrum of the series. --------
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

    fftView_ = new QChartView(fftChart_, container);  // takes ownership of fftChart_
    fftView_->setRenderHint(QPainter::Antialiasing);
    vbox->addWidget(fftView_, 1);

    fftReadout_ = new QLabel(QStringLiteral("—"), container);
    fftReadout_->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    vbox->addWidget(fftReadout_);

    setWidget(container);

    ACONNECT(fitAllBox_.data(), &QCheckBox::toggled, this, &StripChartDock::onFitToggled);
}

void StripChartDock::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    ASSERT_THREAD(this);
    protein_      = protein;
    conformation_ = conformation;
}

void StripChartDock::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    selection_ = selection;
    if (selection) {
        ACONNECT(selection, &model::AtomSelection::changed, this, &StripChartDock::onSelectionChanged);
        ACONNECT(selection, &model::AtomSelection::cleared, this, &StripChartDock::clearSelection);
    }
    rebuildGeometrySeries();
}

void StripChartDock::setFrame(int t) {
    ASSERT_THREAD(this);
    frame_ = t;
    updateCursorAndWindow();
}

void StripChartDock::onSelectionChanged() {
    ASSERT_THREAD(this);
    rebuildGeometrySeries();
}

void StripChartDock::clearSelection() {
    ASSERT_THREAD(this);
    haveSeries_ = false;
    if (geomSeries_) geomSeries_->clear();
    if (cursorSeries_) cursorSeries_->clear();
    if (fftSeries_) fftSeries_->clear();
    if (readout_) readout_->setText(QStringLiteral("—"));
    if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
}

void StripChartDock::onFitToggled(bool on) {
    ASSERT_THREAD(this);
    fitAll_ = on;
    updateCursorAndWindow();
}

void StripChartDock::rebuildGeometrySeries() {
    ASSERT_THREAD(this);
    haveSeries_ = false;

    model::Conformation* conf = conformation_;  // QPointer -> raw via operator T*
    if (!protein_ || !conf || !selection_ || !geomSeries_) {
        if (readout_) readout_->setText(QStringLiteral("—"));
        return;
    }

    auto blank = [&](const QString& msg) {
        geomSeries_->clear();
        cursorSeries_->clear();
        if (fftSeries_) fftSeries_->clear();
        if (readout_) readout_->setText(msg);
        if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
    };

    // Trajectory-only: a single pose has no frame axis worth a strip chart
    // (same gate as the time-series dock).
    if (!conf->asTrajectory()) {
        blank(QStringLiteral("single pose — no trajectory"));
        return;
    }

    const model::GeometryKind kind = selection_->geometryKind();
    if (kind == model::GeometryKind::None) {
        blank(selection_->empty() ? QStringLiteral("—") : QStringLiteral("select ≥ 2 atoms"));
        return;
    }

    const int           T = static_cast<int>(conf->frameCount());
    QList<QPointF>      pts;
    std::vector<double> values;  // the same series, for the FFT
    pts.reserve(T);
    values.reserve(T);

    double yMin      = std::numeric_limits<double>::infinity();
    double yMax      = -std::numeric_limits<double>::infinity();
    double lastValid = 0.0;
    for (int t = 0; t < T; ++t) {
        const auto m = model::Measure(*conf, static_cast<std::size_t>(t), selection_->atoms());
        // Degenerate geometry (coincident/collinear atoms) is undefined; it
        // does not occur in real structures, but if it did we carry the last
        // valid value forward rather than draw a spurious spike to zero.
        if (m.valid) {
            lastValid = m.value;
            yMin      = std::min(yMin, m.value);
            yMax      = std::max(yMax, m.value);
        }
        pts.append(QPointF(t, lastValid));
        values.push_back(lastValid);
    }
    geomSeries_->replace(pts);

    if (!std::isfinite(yMin) || !std::isfinite(yMax)) {
        yMin = 0.0;
        yMax = 1.0;
    }
    if (yMax - yMin < 1e-9) {
        yMin -= 1.0;
        yMax += 1.0;
    } else {
        const double pad = 0.05 * (yMax - yMin);
        yMin -= pad;
        yMax += pad;
    }
    yMin_ = yMin;
    yMax_ = yMax;

    const bool    isDist = (kind == model::GeometryKind::Distance);
    const QString unit   = isDist ? QString::fromUtf8("Å") : QString::fromUtf8("°");
    yAxis_->setTitleText(QString::fromLatin1(model::NameForGeometryKind(kind)).toLower() +
                         QStringLiteral(" (") + unit + QStringLiteral(")"));
    yAxis_->setRange(yMin_, yMax_);

    haveSeries_ = true;
    qCDebug(cStrip).noquote() << "rebuild |" << model::NameForGeometryKind(kind)
                              << "| frames=" << T << "| y=[" << yMin_ << "," << yMax_ << "]";

    // ---- FFT panel: power spectrum of the geometry series. ----------------
    // A dihedral oscillating across ±180° must read as a smooth signal, not a
    // sawtooth — phase-unwrap (a copy) before the transform so the spectrum is
    // not flooded with wrap harmonics. The time-domain panel keeps the physical
    // wrapped values.
    std::vector<double> fftInput = values;
    if (kind == model::GeometryKind::Dihedral) {
        for (std::size_t k = 1; k < fftInput.size(); ++k) {
            double d = fftInput[k] - fftInput[k - 1];
            while (d > 180.0)  { fftInput[k] -= 360.0; d -= 360.0; }
            while (d < -180.0) { fftInput[k] += 360.0; d += 360.0; }
        }
    }

    const double dtPs = (T >= 2)
        ? (conf->timePicoseconds(1) - conf->timePicoseconds(0))
        : 0.0;
    const model::PowerSpectrum spec = model::ComputePowerSpectrum(fftInput, dtPs);

    if (fftSeries_ && spec.valid && spec.power.size() >= 2) {
        QList<QPointF> fftPts;
        fftPts.reserve(static_cast<int>(spec.power.size()));
        double maxP = 0.0;
        // Skip k=0 (the DC bin, ≈ 0 after mean subtraction).
        for (std::size_t k = 1; k < spec.power.size(); ++k) {
            fftPts.append(QPointF(spec.frequencyPerNs[k], spec.power[k]));
            maxP = std::max(maxP, spec.power[k]);
        }
        fftSeries_->replace(fftPts);

        const double fMax = spec.frequencyPerNs.back();
        fftXAxis_->setRange(0.0, fMax > 0.0 ? fMax : 1.0);
        fftYAxis_->setRange(0.0, maxP > 0.0 ? maxP * 1.05 : 1.0);

        if (fftReadout_) {
            const double p = spec.dominantPeriodPs;
            fftReadout_->setText(p >= 1000.0
                ? QStringLiteral("dominant period: %1 ns").arg(p / 1000.0, 0, 'f', 2)
                : QStringLiteral("dominant period: %1 ps").arg(p, 0, 'f', 0));
        }
    } else {
        if (fftSeries_) fftSeries_->clear();
        if (fftReadout_) fftReadout_->setText(QStringLiteral("—"));
    }

    updateCursorAndWindow();
}

void StripChartDock::updateCursorAndWindow() {
    ASSERT_THREAD(this);
    model::Conformation* conf = conformation_;
    if (!xAxis_ || !conf) return;

    const int T = static_cast<int>(conf->frameCount());
    const int t = std::clamp(frame_, 0, std::max(0, T - 1));

    // X window: full range when fitting or short; otherwise a fixed-width
    // window that scrolls to keep the playhead in view (no "scrunch").
    if (fitAll_ || T <= kWindowFrames) {
        xAxis_->setRange(0, std::max(0, T - 1));
    } else {
        const int lo = std::clamp(t - kWindowFrames / 2, 0, T - 1 - kWindowFrames);
        xAxis_->setRange(lo, lo + kWindowFrames);
    }

    // Playhead at x = t spanning the data range.
    if (haveSeries_ && cursorSeries_) {
        cursorSeries_->replace(QList<QPointF>{QPointF(t, yMin_), QPointF(t, yMax_)});
    }

    // Time-domain readout — the measured value at the current frame.
    if (readout_ && haveSeries_ && selection_) {
        const auto m = model::Measure(*conf, static_cast<std::size_t>(t), selection_->atoms());
        if (m.valid) {
            const bool    isDist = (selection_->geometryKind() == model::GeometryKind::Distance);
            const QString unit   = isDist ? QString::fromUtf8("Å") : QString::fromUtf8("°");
            readout_->setText(QStringLiteral("%1: %2 %3   (frame %4 / %5)")
                                  .arg(QString::fromLatin1(model::NameForGeometryKind(selection_->geometryKind())))
                                  .arg(m.value, 0, 'f', isDist ? 3 : 1)
                                  .arg(unit)
                                  .arg(t + 1)
                                  .arg(T));
        }
    }
}

}  // namespace h5reader::app
