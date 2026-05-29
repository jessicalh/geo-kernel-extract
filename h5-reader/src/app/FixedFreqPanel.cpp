#include "FixedFreqPanel.h"

#include <QPainter>
#include <QPen>
#include <QPolygonF>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace h5reader::app {

namespace {

constexpr double kMarkerRadiusPx = 4.0;
const QColor kTraceColor = QColor(74, 184, 220);  // soft cyan
const QColor kMarkerEdge = QColor(255, 255, 255, 200);

// Map ω (rad/s) to a log-x position within the plot rect. The J(ω)
// frequencies span ~0 (the constant term) to ω_H+ω_N (~10⁹ rad/s); log
// scale collapses that span into a readable axis. The ω=0 frequency
// is special-cased to the left margin since log(0) is undefined.
double xForFreq(double omega, double logMin, double logMax, const QRectF& plot)
{
    if (!std::isfinite(omega) || omega <= 0.0)
        return plot.left();
    const double l = std::log10(omega);
    const double t = (l - logMin) / std::max(1e-9, logMax - logMin);
    return plot.left() + std::clamp(t, 0.0, 1.0) * plot.width();
}

double yForValue(double v, double yMin, double yMax, const QRectF& plot)
{
    if (!std::isfinite(v))
        return plot.bottom();
    const double t = (v - yMin) / std::max(1e-9, yMax - yMin);
    return plot.bottom() - std::clamp(t, 0.0, 1.0) * plot.height();
}

}  // namespace

FixedFreqPanel::FixedFreqPanel(QString label,
                               std::unique_ptr<model::QtPerBondVectorFixedFreqBlock> data,
                               std::size_t row,
                               model::SignalBinding revealBinding)
    : label_(std::move(label)),
      data_(std::move(data)),
      row_(row),
      reveal_(std::move(revealBinding)),
      haveBinding_(!reveal_.descriptorId.isEmpty())
{}

void FixedFreqPanel::paint(QPainter& p,
                           const PanelGeometry& geometry,
                           const PaintContext& context) const
{
    paintPanelBackground(p, geometry);
    const QString title = data_ && !data_->units.isEmpty()
                              ? QStringLiteral("%1 (%2)").arg(label_, data_->units)
                              : label_;
    paintHeader(p,
                geometry,
                HeaderText{title, QString()},
                12,
                kStripTextMuted,
                haveBinding_,
                context.hasHover && geometry.reveal.contains(context.hoverPos));
    paintGrid(p, geometry.plot);

    if (!data_ || data_->n_freqs == 0 || row_ >= data_->n_vectors) {
        paintPlotBorder(p, geometry.plot);
        return;
    }

    // Collect finite (ω, J) pairs. NH-only J(ω) means non-NH rows are
    // NaN — we filter those at paint so the panel just stays empty for
    // a Cα-Hα / C=O pick rather than drawing a flat zero line.
    std::vector<std::pair<double, double>> samples;
    samples.reserve(data_->n_freqs);
    for (std::size_t i = 0; i < data_->n_freqs; ++i) {
        const double omega = (i < data_->freq_values.size()) ? data_->freq_values[i] : 0.0;
        const double v = data_->at(row_, i);
        if (!std::isfinite(omega) || !std::isfinite(v))
            continue;
        samples.emplace_back(omega, v);
    }
    if (samples.empty()) {
        paintPlotBorder(p, geometry.plot);
        return;
    }

    // Log-x axis bounds. Skip ω=0 in the bounds (it's pinned to the
    // left margin); positive frequencies set the log range.
    double logMin = std::numeric_limits<double>::infinity();
    double logMax = -std::numeric_limits<double>::infinity();
    for (const auto& [omega, v] : samples) {
        if (omega <= 0.0) continue;
        const double l = std::log10(omega);
        logMin = std::min(logMin, l);
        logMax = std::max(logMax, l);
    }
    if (!std::isfinite(logMin) || !std::isfinite(logMax)) {
        // All ω=0 (shouldn't happen for KTB J(ω) but defensive). Just
        // pin everything to the left margin.
        logMin = logMax = 0.0;
    }
    // Pad the log range so the leftmost/rightmost markers aren't on
    // the axis ticks.
    const double pad = 0.15 * std::max(1.0, logMax - logMin);
    logMin -= pad;
    logMax += pad;

    // Y bounds: include 0 so the J=0 floor is meaningful; pad a bit.
    double yMin = 0.0;
    double yMax = 0.0;
    for (const auto& [_, v] : samples) {
        yMin = std::min(yMin, v);
        yMax = std::max(yMax, v);
    }
    if (yMax - yMin < 1e-12) yMax = yMin + 1.0;
    const double ypad = (yMax - yMin) * 0.1;
    yMin -= ypad;
    yMax += ypad;

    ValueRange range;
    range.min = yMin;
    range.max = yMax;
    range.valid = true;
    paintYAxisLabels(p, geometry, range);

    // Polyline through the samples (left → right by frequency).
    std::vector<std::pair<double, double>> sorted = samples;
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    QPolygonF poly;
    poly.reserve(static_cast<int>(sorted.size()));
    for (const auto& [omega, v] : sorted) {
        const double x = xForFreq(omega, logMin, logMax, geometry.plot);
        const double y = yForValue(v, yMin, yMax, geometry.plot);
        poly.append(QPointF(x, y));
    }
    p.setRenderHint(QPainter::Antialiasing, true);
    if (poly.size() >= 2) {
        p.setPen(QPen(kTraceColor, 1.8, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        p.drawPolyline(poly);
    }

    // Markers (5 discrete points — the KTB Larmor combinations).
    p.setPen(QPen(kMarkerEdge, 1.0));
    p.setBrush(kTraceColor);
    for (const auto& [omega, v] : sorted) {
        const QPointF c(xForFreq(omega, logMin, logMax, geometry.plot),
                        yForValue(v, yMin, yMax, geometry.plot));
        p.drawEllipse(c, kMarkerRadiusPx, kMarkerRadiusPx);
    }

    // Codex LATER-3 (2026-05-29): label each marker with its KTB
    // Larmor combination so the user knows which discrete ω the dot
    // represents. The 5 KTB frequencies (per producer convention) are
    // [0, ω_N, ω_H−ω_N, ω_H, ω_H+ω_N] in order. Label the FIRST 5
    // markers with that order; any additional markers (defensive
    // case) get a generic ω_n label.
    static const char* kKtbLabels[] = {"0", "ωN", "ωH-ωN", "ωH", "ωH+ωN"};
    p.setFont(monoFont(9));
    p.setPen(kStripTextMuted);
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        const auto& [omega, v] = sorted[i];
        const double x = xForFreq(omega, logMin, logMax, geometry.plot);
        const QString label = (i < 5)
                                  ? QString::fromUtf8(kKtbLabels[i])
                                  : QStringLiteral("ω%1").arg(i);
        const QRectF labelRect(x - 30.0, geometry.plot.bottom() + 4.0, 60.0, 14.0);
        p.drawText(labelRect, Qt::AlignHCenter | Qt::AlignVCenter, label);
    }

    paintPlotBorder(p, geometry.plot);
}

}  // namespace h5reader::app
