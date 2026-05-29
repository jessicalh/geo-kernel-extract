#include "LagDecayPanel.h"

#include <QFontMetricsF>
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

const std::array<QColor, 13> kChannelPalette{
    QColor(115, 229, 214), QColor(255, 175, 76),  QColor(135, 211, 124),
    QColor(241, 91, 181),  QColor(118, 188, 247), QColor(255, 232, 105),
    QColor(170, 116, 247), QColor(245, 121, 121), QColor(76, 200, 209),
    QColor(255, 145, 145), QColor(199, 217, 92),  QColor(159, 199, 255),
    QColor(255, 215, 158),
};

}  // namespace

LagDecayPanel::LagDecayPanel(QString label,
                             const model::QtPerAtomChannelCurve* data,
                             std::size_t atomRow,
                             model::SignalBinding revealBinding)
    : label_(std::move(label)),
      data_(data),
      atomRow_(atomRow),
      reveal_(std::move(revealBinding)),
      haveBinding_(!reveal_.descriptorId.isEmpty())
{}

LagDecayPanel::LagDecayPanel(QString label,
                             std::unique_ptr<model::QtPerAtomChannelCurve> ownedData,
                             std::size_t atomRow,
                             model::SignalBinding revealBinding)
    : label_(std::move(label)),
      owned_(std::move(ownedData)),
      atomRow_(atomRow),
      reveal_(std::move(revealBinding)),
      haveBinding_(!reveal_.descriptorId.isEmpty())
{
    data_ = owned_.get();
}

void LagDecayPanel::paint(QPainter& p,
                          const PanelGeometry& geometry,
                          const PaintContext& context) const
{
    paintPanelBackground(p, geometry);

    const QString title = label_.isEmpty()
                              ? QStringLiteral("Autocorrelation (lag)")
                              : label_;
    const QString readout = data_
        ? QStringLiteral("atom %1 · %2 channels").arg(atomRow_).arg(data_->n_channels)
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

    const double xMin = 0.0;
    double xMax = data_->axis_values.empty()
                      ? static_cast<double>(data_->n_samples)
                      : data_->axis_values.back();
    if (xMax <= xMin) xMax = xMin + 1.0;

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
    if (range.valid) range = padRange(range);

    paintSpectrumTicks(p, geometry.plot, xMin, xMax);
    paintYAxisLabels(p, geometry, range);

    if (range.valid) {
        p.setRenderHint(QPainter::Antialiasing, true);
        for (std::size_t c = 0; c < data_->n_channels; ++c) {
            QPolygonF poly;
            poly.reserve(static_cast<int>(data_->n_samples));
            for (std::size_t s = 0; s < data_->n_samples; ++s) {
                const double lag = (s < data_->axis_values.size())
                                       ? data_->axis_values[s]
                                       : static_cast<double>(s);
                const double v = data_->at(atomRow_, c, s);
                if (!std::isfinite(v) || !std::isfinite(lag)) continue;
                const double x = geometry.plot.left()
                                 + geometry.plot.width() * (lag - xMin) / (xMax - xMin);
                const double y = yForValue(v, range, geometry.plot);
                poly.append(QPointF(x, y));
            }
            if (poly.size() < 2) continue;
            const QColor& col = kChannelPalette[c % kChannelPalette.size()];
            p.setPen(QPen(col, 1.4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
            p.drawPolyline(poly);
        }
    }

    // Animated lag cursor: map currentFrame → lag axis position.
    // Currentframe runs over the trajectory; we plot it against the
    // lag axis (also in ps) so the cursor walks left-to-right with
    // playback. Clamp at xMax so the cursor parks at the right edge
    // for frames past the lag window.
    if (context.viewport && range.valid) {
        // sample interval = (xMax - xMin) / (n_samples - 1) — approximate
        // mapping from frame index to lag time. Honest fallback when the
        // axis_values aren't dense.
        const double dt = (data_->n_samples > 1)
                              ? (data_->axis_values.back() - data_->axis_values.front())
                                    / static_cast<double>(data_->n_samples - 1)
                              : 1.0;
        const double lagAtFrame = std::clamp(
            static_cast<double>(context.currentFrame) * dt, xMin, xMax);
        const double cursorX = geometry.plot.left()
            + geometry.plot.width() * (lagAtFrame - xMin) / (xMax - xMin);
        p.setPen(QPen(kStripCursor, 1.2, Qt::DashLine));
        p.drawLine(QPointF(cursorX + 0.5, geometry.plot.top()),
                   QPointF(cursorX + 0.5, geometry.plot.bottom()));
    }

    paintPlotBorder(p, geometry.plot);
}

}  // namespace h5reader::app
