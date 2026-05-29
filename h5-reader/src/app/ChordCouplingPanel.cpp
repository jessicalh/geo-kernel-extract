#include "ChordCouplingPanel.h"

#include <QFontMetricsF>
#include <QPainter>
#include <QPainterPath>
#include <QPen>
#include <QPointF>
#include <QRectF>
#include <QString>
#include <QStringList>

#include <algorithm>
#include <cmath>
#include <utility>

namespace h5reader::app {

ChordCouplingPanel::ChordCouplingPanel(QString label,
                                       const model::QtPerAtomMatrix* matrix,
                                       std::size_t atomRow,
                                       double threshold,
                                       model::SignalBinding revealBinding)
    : label_(std::move(label)),
      matrix_(matrix),
      atomRow_(atomRow),
      threshold_(threshold),
      reveal_(std::move(revealBinding)),
      haveBinding_(!reveal_.descriptorId.isEmpty())
{}

void ChordCouplingPanel::paint(QPainter& p,
                                const PanelGeometry& geometry,
                                const PaintContext& context) const
{
    paintPanelBackground(p, geometry);

    const QString title = label_.isEmpty()
                              ? QStringLiteral("Channel coherence")
                              : label_;
    const QString readout = matrix_
        ? QStringLiteral("atom %1 · |r| ≥ %2")
              .arg(atomRow_).arg(threshold_, 0, 'f', 2)
        : QString();
    paintHeader(p, geometry,
                HeaderText{title, readout},
                160,
                kStripTextMuted,
                /*hasBinding=*/haveBinding_,
                context.hasHover && geometry.reveal.contains(context.hoverPos));

    if (!matrix_ || matrix_->n_channels == 0 || atomRow_ >= matrix_->n_atoms) {
        paintPlotBorder(p, geometry.plot);
        return;
    }

    const QRectF& plot = geometry.plot;
    const QPointF center = plot.center();
    const double radius = std::min(plot.width(), plot.height()) * 0.5 - 36.0;
    if (radius < 20.0) {
        paintPlotBorder(p, geometry.plot);
        return;
    }

    const std::size_t C = matrix_->n_channels;
    p.setRenderHint(QPainter::Antialiasing, true);

    // Compute channel anchor positions around the circle (top = 12 o'clock).
    std::vector<QPointF> nodePositions(C);
    for (std::size_t c = 0; c < C; ++c) {
        const double theta = -M_PI / 2.0 + 2.0 * M_PI * c / static_cast<double>(C);
        nodePositions[c] = QPointF(center.x() + radius * std::cos(theta),
                                    center.y() + radius * std::sin(theta));
    }

    // Draw arcs for |r| ≥ threshold. Hue: positive=teal, negative=magenta.
    // Thickness scales with |r|; we draw an inner cubic Bezier toward the
    // center to give the classic chord curvature.
    for (std::size_t a = 0; a < C; ++a) {
        for (std::size_t b = a + 1; b < C; ++b) {
            const double r = matrix_->at(atomRow_, a, b);
            if (!std::isfinite(r) || std::abs(r) < threshold_) continue;
            const double mag = std::abs(r);
            const double width = 1.0 + 5.0 * mag;
            QColor col = (r >= 0.0) ? QColor(115, 229, 214) : QColor(241, 91, 181);
            col.setAlphaF(std::clamp(0.45 + 0.55 * mag, 0.45, 1.0));
            QPainterPath path(nodePositions[a]);
            path.cubicTo(center, center, nodePositions[b]);
            p.setPen(QPen(col, width, Qt::SolidLine, Qt::RoundCap));
            p.drawPath(path);
        }
    }

    // Draw channel nodes (small circles + labels).
    p.setFont(monoFont(9));
    for (std::size_t c = 0; c < C; ++c) {
        const QPointF& pt = nodePositions[c];
        p.setPen(Qt::NoPen);
        p.setBrush(kStripText);
        p.drawEllipse(pt, 3.0, 3.0);

        const QString name = (c < static_cast<std::size_t>(matrix_->channel_names.size()))
                                 ? matrix_->channel_names[c]
                                 : QStringLiteral("ch%1").arg(c);
        const double theta = -M_PI / 2.0 + 2.0 * M_PI * c / static_cast<double>(C);
        const QPointF labelPt(center.x() + (radius + 14.0) * std::cos(theta),
                               center.y() + (radius + 14.0) * std::sin(theta));
        const QFontMetricsF fm(p.font());
        const double w = fm.horizontalAdvance(name);
        const QRectF labelRect(labelPt.x() - w * 0.5, labelPt.y() - 8.0, w + 2.0, 16.0);
        p.setPen(kStripTextMuted);
        p.drawText(labelRect, Qt::AlignCenter, name);
    }

    paintPlotBorder(p, geometry.plot);
    (void)context;
}

}  // namespace h5reader::app
