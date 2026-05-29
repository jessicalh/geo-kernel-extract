// AbstractStripPanel.cpp — paint + geometry helper implementations.
// Direct lift of the helpers from StripStackWidget.cpp's previous
// anonymous namespace (PATTERNS.md 17: code moved verbatim, no
// behaviour change), now visible to every panel subclass.

#include "AbstractStripPanel.h"

#include "TimeViewportController.h"

#include <QFont>
#include <QFontMetricsF>
#include <QPainter>
#include <QStringList>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

QFont uiFont(int px, QFont::Weight weight)
{
    QFont f;
    f.setFamilies(QStringList{QStringLiteral("Inter"), QStringLiteral("Helvetica"),
                              QStringLiteral("Arial"), QStringLiteral("sans-serif")});
    f.setPixelSize(px);
    f.setWeight(weight);
    f.setStyleStrategy(QFont::PreferAntialias);
    return f;
}

QFont monoFont(int px)
{
    QFont f;
    f.setFamilies(QStringList{QStringLiteral("Roboto Mono"), QStringLiteral("DejaVu Sans Mono"),
                              QStringLiteral("monospace")});
    f.setPixelSize(px);
    f.setStyleStrategy(QFont::PreferAntialias);
    return f;
}

QString fmtValue(double value, int decimals)
{
    return QString::number(value, 'f', decimals);
}

double TimeScale::xForFrame(int frame, const QRectF& plot) const
{
    if (last <= first)
        return plot.left() + plot.width() * 0.5;
    return plot.left() + plot.width() * (frame - first) / std::max(1, last - first);
}

int TimeScale::frameAt(const QPoint& pos, const QRectF& plot) const
{
    const double xNorm = std::clamp(
        (pos.x() - plot.left()) / std::max(1.0, plot.width()), 0.0, 1.0);
    return first + static_cast<int>(std::round(xNorm * std::max(0, last - first)));
}

QRectF panelRectForIndex(const StackGeometry& stack, int panelIndex)
{
    const int n = std::max(1, stack.panelCount);
    const double h = (stack.viewportSize.height() - (n + 1) * kStripPanelGap)
                     / static_cast<double>(n);
    return QRectF(kStripPanelGap,
                  kStripPanelGap + panelIndex * (h + kStripPanelGap),
                  stack.viewportSize.width() - 2 * kStripPanelGap,
                  h);
}

QRectF plotRectForPanel(const QRectF& r)
{
    if (r.height() < 96.0)
        return r.adjusted(48, 28, -10, -20);
    return r.adjusted(52, 30, -12, -24);
}

QRectF revealRectForPanel(const QRectF& r)
{
    return QRectF(r.left() + 8.0, r.top() + 6.0, 18.0, 18.0);
}

QRectF letterboxedPanelRect(const QRectF& panel, std::optional<double> aspect)
{
    if (!aspect || !(*aspect > 0.0))
        return panel;
    const double a = *aspect;                       // width / height
    const double w = panel.width();
    const double h = panel.height();
    if (w / h > a) {                                // too wide → narrow horizontally
        const double newW = h * a;
        return QRectF(panel.left() + (w - newW) * 0.5, panel.top(), newW, h);
    }
    const double newH = w / a;                      // too tall → narrow vertically
    return QRectF(panel.left(), panel.top() + (h - newH) * 0.5, w, newH);
}

PanelGeometry panelGeometryForIndex(const StackGeometry& stack,
                                    int panelIndex,
                                    std::optional<double> preferredAspect)
{
    PanelGeometry geometry;
    const QRectF rawPanel = panelRectForIndex(stack, panelIndex);
    geometry.panel = letterboxedPanelRect(rawPanel, preferredAspect);
    geometry.plot = plotRectForPanel(geometry.panel);
    geometry.reveal = revealRectForPanel(geometry.panel);
    return geometry;
}

ValueRange padRange(ValueRange range)
{
    if (!range.valid)
        return range;
    if (std::abs(range.max - range.min) < 1e-12) {
        range.min -= 1.0;
        range.max += 1.0;
    } else {
        const double pad = (range.max - range.min) * 0.08;
        range.min -= pad;
        range.max += pad;
    }
    return range;
}

double yForValue(double value, const ValueRange& range, const QRectF& plot)
{
    return plot.bottom() - plot.height() * (value - range.min)
                              / std::max(1e-12, range.max - range.min);
}

void drawRevealButton(QPainter& p, const QRectF& r, bool hover)
{
    p.save();
    p.setRenderHint(QPainter::Antialiasing, true);
    QColor fill = kStripReveal;
    fill.setAlpha(hover ? 58 : 30);
    p.setBrush(fill);
    p.setPen(QPen(kStripReveal.lighter(120), 1.1));
    p.drawEllipse(r.adjusted(1.0, 1.0, -1.0, -1.0));

    const QPointF c = r.center();
    p.setBrush(Qt::NoBrush);
    p.setPen(QPen(kStripReveal.lighter(140), 1.1));
    p.drawEllipse(c, 3.6, 3.6);
    p.drawLine(QPointF(c.x() - 6.0, c.y()), QPointF(c.x() - 3.5, c.y()));
    p.drawLine(QPointF(c.x() + 3.5, c.y()), QPointF(c.x() + 6.0, c.y()));
    p.drawLine(QPointF(c.x(), c.y() - 6.0), QPointF(c.x(), c.y() - 3.5));
    p.drawLine(QPointF(c.x(), c.y() + 3.5), QPointF(c.x(), c.y() + 6.0));
    p.restore();
}

void paintPanelBackground(QPainter& p, const PanelGeometry& geometry)
{
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setPen(QPen(kStripPanelBorder, 1.0));
    p.setBrush(kStripPanel);
    p.drawRoundedRect(geometry.panel.adjusted(0.5, 0.5, -0.5, -0.5), 5, 5);
}

void paintHeader(QPainter& p,
                 const PanelGeometry& geometry,
                 const HeaderText& text,
                 int readoutWidth,
                 const QColor& readoutColor,
                 bool hasBinding,
                 bool revealHover)
{
    if (hasBinding)
        drawRevealButton(p, geometry.reveal, revealHover);

    p.setFont(uiFont(13, QFont::Medium));
    p.setPen(kStripText);
    const double titleLeft = geometry.panel.left() + (hasBinding ? 34.0 : 10.0);
    const QRectF titleRect(titleLeft,
                           geometry.panel.top() + 6.0,
                           std::max(20.0, geometry.panel.right() - titleLeft - readoutWidth - 10.0),
                           22.0);
    p.drawText(titleRect,
               Qt::AlignLeft | Qt::AlignVCenter,
               QFontMetricsF(p.font()).elidedText(text.title,
                                                  Qt::ElideMiddle,
                                                  static_cast<int>(titleRect.width())));

    if (!text.readout.isEmpty()) {
        p.setFont(monoFont(12));
        p.setPen(readoutColor);
        p.drawText(geometry.panel.adjusted(10, 6, -10, -geometry.panel.height() + 28),
                   Qt::AlignRight | Qt::AlignVCenter,
                   text.readout);
    }
}

void paintGrid(QPainter& p, const QRectF& plot)
{
    p.setRenderHint(QPainter::Antialiasing, false);
    p.setPen(QPen(kStripGrid, 1.0));
    for (int g = 0; g <= 3; ++g) {
        const double y = plot.top() + plot.height() * g / 3.0;
        p.drawLine(QPointF(plot.left(), std::round(y) + 0.5),
                   QPointF(plot.right(), std::round(y) + 0.5));
    }
    for (int g = 0; g <= 4; ++g) {
        const double x = plot.left() + plot.width() * g / 4.0;
        p.drawLine(QPointF(std::round(x) + 0.5, plot.top()),
                   QPointF(std::round(x) + 0.5, plot.bottom()));
    }
}

void paintYAxisLabels(QPainter& p, const PanelGeometry& geometry, const ValueRange& range)
{
    if (!range.valid)
        return;
    p.setFont(monoFont(10));
    p.setPen(kStripTextMuted);
    p.drawText(QRectF(geometry.panel.left() + 8, geometry.plot.top() - 8, 38, 16),
               Qt::AlignRight | Qt::AlignVCenter, fmtValue(range.max, 2));
    p.drawText(QRectF(geometry.panel.left() + 8, geometry.plot.bottom() - 8, 38, 16),
               Qt::AlignRight | Qt::AlignVCenter, fmtValue(range.min, 2));
}

void paintPlotBorder(QPainter& p, const QRectF& plot)
{
    p.setRenderHint(QPainter::Antialiasing, false);
    p.setPen(QPen(kStripPanelBorder, 1.0));
    p.setBrush(Qt::NoBrush);
    p.drawRect(plot.adjusted(0.5, 0.5, -0.5, -0.5));
}

void paintTimeTicks(QPainter& p, const QRectF& plot, const TimeScale& time)
{
    p.setFont(monoFont(10));
    p.setPen(kStripTextMuted);
    for (int g = 0; g <= 4; ++g) {
        const double frac = g / 4.0;
        const double x = plot.left() + plot.width() * frac;
        const int tickFrame = time.first + static_cast<int>(std::round((time.last - time.first) * frac));
        const QRectF labelRect(x - 28.0, plot.bottom() + 4.0, 56.0, 16.0);
        p.drawText(labelRect, Qt::AlignHCenter | Qt::AlignVCenter,
                   QStringLiteral("f%1").arg(tickFrame + 1));
    }
}

void paintSpectrumTicks(QPainter& p, const QRectF& plot, double xMin, double xMax)
{
    p.setFont(monoFont(10));
    p.setPen(kStripTextMuted);
    for (int g = 0; g <= 4; ++g) {
        const double frac = g / 4.0;
        const double x = plot.left() + plot.width() * frac;
        const double tick = xMin + (xMax - xMin) * frac;
        const QRectF labelRect(x - 32.0, plot.bottom() + 4.0, 64.0, 16.0);
        p.drawText(labelRect, Qt::AlignHCenter | Qt::AlignVCenter,
                   QStringLiteral("%1").arg(tick, 0, 'f', tick < 10.0 ? 2 : 1));
    }
}

void paintSelectedTimeRange(QPainter& p, const PanelGeometry& geometry, const PaintContext& context)
{
    if (!context.viewport || !context.viewport->hasSelectedRange())
        return;
    const double a = context.time.xForFrame(context.viewport->selectedStart(), geometry.plot);
    const double b = context.time.xForFrame(context.viewport->selectedEnd(), geometry.plot);
    p.fillRect(QRectF(QPointF(std::min(a, b), geometry.plot.top()),
                      QPointF(std::max(a, b), geometry.plot.bottom())),
               kStripSelection);
}

void paintTemporalCursor(QPainter& p, const PanelGeometry& geometry, const PaintContext& context)
{
    const double cursorX = context.time.xForFrame(context.currentFrame, geometry.plot);
    p.setPen(QPen(kStripCursor, 1.2, Qt::DashLine));
    p.drawLine(QPointF(cursorX + 0.5, geometry.plot.top()),
               QPointF(cursorX + 0.5, geometry.plot.bottom()));

    if (context.hasHover && geometry.plot.contains(context.hoverPos)) {
        p.setPen(QPen(kStripHover, 1.0, Qt::DotLine));
        p.drawLine(QPointF(context.hoverPos.x() + 0.5, geometry.plot.top()),
                   QPointF(context.hoverPos.x() + 0.5, geometry.plot.bottom()));
    }
}

}  // namespace h5reader::app
