#include "PowerSpectrumPanel.h"

#include <QFontMetricsF>
#include <QMouseEvent>
#include <QPainter>
#include <QPen>
#include <QPolygonF>
#include <QRectF>
#include <QString>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

namespace h5reader::app {

namespace {

// Distinct colours for up to 13 channels; chosen for legibility on the
// dark canvas. Wraps if there are more channels (none today).
const std::array<QColor, 13> kChannelPalette{
    QColor(115, 229, 214), QColor(255, 175, 76),  QColor(135, 211, 124),
    QColor(241, 91, 181),  QColor(118, 188, 247), QColor(255, 232, 105),
    QColor(170, 116, 247), QColor(245, 121, 121), QColor(76, 200, 209),
    QColor(255, 145, 145), QColor(199, 217, 92),  QColor(159, 199, 255),
    QColor(255, 215, 158),
};

}  // namespace

PowerSpectrumPanel::PowerSpectrumPanel(QString label,
                                       const model::QtPerAtomChannelCurve* data,
                                       std::size_t atomRow,
                                       model::SignalBinding revealBinding)
    : label_(std::move(label)),
      data_(data),
      atomRow_(atomRow),
      reveal_(std::move(revealBinding)),
      haveBinding_(!reveal_.descriptorId.isEmpty())
{}

QColor PowerSpectrumPanel::channelColor(std::size_t channelIndex) const
{
    return kChannelPalette[channelIndex % kChannelPalette.size()];
}

void PowerSpectrumPanel::paint(QPainter& p,
                               const PanelGeometry& geometry,
                               const PaintContext& context) const
{
    paintPanelBackground(p, geometry);

    const QString title = label_.isEmpty()
                              ? QStringLiteral("Power spectrum")
                              : label_;
    const QString readout = data_
        ? QStringLiteral("atom %1 · %2 channels")
              .arg(atomRow_).arg(data_->n_channels)
        : QString();
    paintHeader(p, geometry,
                HeaderText{title, readout},
                160,
                kStripTextMuted,
                /*hasBinding=*/haveBinding_,
                context.hasHover && geometry.reveal.contains(context.hoverPos));
    paintGrid(p, geometry.plot);

    if (!data_ || data_->n_samples < 2 || data_->n_channels == 0
        || atomRow_ >= data_->n_atoms) {
        paintPlotBorder(p, geometry.plot);
        return;
    }

    // Frequency axis range — always starts at 0 for a one-sided PSD.
    double xMin = data_->axis_values.empty() ? 0.0 : data_->axis_values.front();
    double xMax = data_->axis_values.empty() ? 1.0 : data_->axis_values.back();
    if (xMin > 0.0) xMin = 0.0;
    if (xMax <= xMin) xMax = xMin + 1.0;

    // Value-range walk across all channels for this atom.
    ValueRange range;
    range.valid = false;
    for (std::size_t c = 0; c < data_->n_channels; ++c) {
        for (std::size_t s = 0; s < data_->n_samples; ++s) {
            const double v = data_->at(atomRow_, c, s);
            if (!std::isfinite(v)) continue;
            if (!range.valid) { range.min = range.max = v; range.valid = true; }
            else {
                range.min = std::min(range.min, v);
                range.max = std::max(range.max, v);
            }
        }
    }
    if (range.valid && range.min > 0.0) range.min = 0.0;
    if (range.valid) range = padRange(range);

    paintSpectrumTicks(p, geometry.plot, xMin, xMax);
    paintYAxisLabels(p, geometry, range);

    if (range.valid) {
        p.setRenderHint(QPainter::Antialiasing, true);
        for (std::size_t c = 0; c < data_->n_channels; ++c) {
            QPolygonF poly;
            poly.reserve(static_cast<int>(data_->n_samples));
            for (std::size_t s = 0; s < data_->n_samples; ++s) {
                const double freq = (s < data_->axis_values.size())
                                        ? data_->axis_values[s] : 0.0;
                const double v = data_->at(atomRow_, c, s);
                if (!std::isfinite(v) || !std::isfinite(freq)) continue;
                const double x = geometry.plot.left()
                                 + geometry.plot.width() * (freq - xMin) / (xMax - xMin);
                const double y = yForValue(v, range, geometry.plot);
                poly.append(QPointF(x, y));
            }
            if (poly.size() < 2) continue;
            p.setPen(QPen(channelColor(c), 1.4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
            p.drawPolyline(poly);
        }
        // Legend strip along the top of the plot.
        p.setFont(monoFont(9));
        double legendX = geometry.plot.left() + 4.0;
        const double legendY = geometry.plot.top() + 4.0;
        for (std::size_t c = 0; c < data_->n_channels; ++c) {
            const QString name = (c < static_cast<std::size_t>(data_->channel_names.size()))
                                     ? data_->channel_names[c]
                                     : QStringLiteral("ch%1").arg(c);
            p.setPen(channelColor(c));
            p.drawText(QRectF(legendX, legendY, 80.0, 12.0),
                       Qt::AlignLeft | Qt::AlignVCenter, name);
            legendX += QFontMetricsF(p.font()).horizontalAdvance(name) + 12.0;
            if (legendX > geometry.plot.right() - 80.0) break;
        }
    }

    paintPlotBorder(p, geometry.plot);
    (void)context;
}

}  // namespace h5reader::app
