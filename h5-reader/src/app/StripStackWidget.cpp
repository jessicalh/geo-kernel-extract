#include "StripStackWidget.h"

#include "TimeViewportController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/StripChartChannel.h"

#include <QApplication>
#include <QEvent>
#include <QFont>
#include <QFontMetricsF>
#include <QMouseEvent>
#include <QPainter>
#include <QPolygonF>
#include <QSize>
#include <QStringList>
#include <QToolTip>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace h5reader::app {

namespace {
constexpr int kMinTrackHeight = 72;
constexpr int kGap = 6;

const QColor kCanvas(17, 18, 23);
const QColor kPanel(24, 27, 31);
const QColor kPanelBorder(43, 48, 56);
const QColor kGrid(48, 54, 66);
const QColor kText(216, 217, 218);
const QColor kTextMuted(154, 160, 166);
const QColor kCursor(217, 26, 26);
const QColor kSelection(87, 148, 242, 55);
const QColor kHover(199, 208, 217);
const QColor kReveal(115, 229, 214);

QFont uiFont(int px, QFont::Weight weight = QFont::Normal)
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

QString fmt(double value, int decimals = 2)
{
    return QString::number(value, 'f', decimals);
}

void drawRevealButton(QPainter& p, const QRectF& r, bool hover)
{
    p.save();
    p.setRenderHint(QPainter::Antialiasing, true);
    QColor fill = kReveal;
    fill.setAlpha(hover ? 58 : 30);
    p.setBrush(fill);
    p.setPen(QPen(kReveal.lighter(120), 1.1));
    p.drawEllipse(r.adjusted(1.0, 1.0, -1.0, -1.0));

    const QPointF c = r.center();
    p.setBrush(Qt::NoBrush);
    p.setPen(QPen(kReveal.lighter(140), 1.1));
    p.drawEllipse(c, 3.6, 3.6);
    p.drawLine(QPointF(c.x() - 6.0, c.y()), QPointF(c.x() - 3.5, c.y()));
    p.drawLine(QPointF(c.x() + 3.5, c.y()), QPointF(c.x() + 6.0, c.y()));
    p.drawLine(QPointF(c.x(), c.y() - 6.0), QPointF(c.x(), c.y() - 3.5));
    p.drawLine(QPointF(c.x(), c.y() + 3.5), QPointF(c.x(), c.y() + 6.0));
    p.restore();
}

struct ValueRange {
    double min = 0.0;
    double max = 1.0;
    bool valid = false;
};

struct StackGeometry {
    QSize viewportSize;
    int panelCount = 0;
};

struct PanelGeometry {
    QRectF panel;
    QRectF plot;
    QRectF reveal;
};

struct HeaderText {
    QString title;
    QString readout;
};

struct FrameWindow {
    int first = 0;
    int last = 0;

    bool valid() const { return first <= last; }
};

struct TimeScale {
    int first = 0;
    int last = 0;

    double xForFrame(int frame, const QRectF& plot) const
    {
        if (last <= first)
            return plot.left() + plot.width() * 0.5;
        return plot.left() + plot.width() * (frame - first) / std::max(1, last - first);
    }

    int frameAt(const QPoint& pos, const QRectF& plot) const
    {
        const double xNorm = std::clamp((pos.x() - plot.left()) / std::max(1.0, plot.width()), 0.0, 1.0);
        return first + static_cast<int>(std::round(xNorm * std::max(0, last - first)));
    }
};

struct PaintContext {
    const TimeViewportController* viewport = nullptr;
    int currentFrame = 0;
    bool hasHover = false;
    QPoint hoverPos;
    TimeScale time;
};

QRectF panelRectForIndex(const StackGeometry& stack, int panelIndex)
{
    const int n = std::max(1, stack.panelCount);
    const double h = (stack.viewportSize.height() - (n + 1) * kGap) / static_cast<double>(n);
    return QRectF(kGap, kGap + panelIndex * (h + kGap), stack.viewportSize.width() - 2 * kGap, h);
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

PanelGeometry panelGeometryForIndex(const StackGeometry& stack, int panelIndex)
{
    PanelGeometry geometry;
    geometry.panel = panelRectForIndex(stack, panelIndex);
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

ValueRange visibleTemporalRange(const model::ChannelBuffer& buffer, const TimeScale& time)
{
    ValueRange range;
    if (buffer.empty())
        return range;

    const int lo = std::max(0, time.first);
    const int hi = std::min(static_cast<int>(buffer.lastFrame()), time.last);
    for (int frame = lo; frame <= hi; ++frame) {
        const std::size_t idx = static_cast<std::size_t>(frame);
        if (idx >= buffer.valid.size() || !buffer.valid[idx])
            continue;
        const double v = buffer.values[idx];
        if (!std::isfinite(v))
            continue;
        if (!range.valid) {
            range.min = range.max = v;
            range.valid = true;
        } else {
            range.min = std::min(range.min, v);
            range.max = std::max(range.max, v);
        }
    }
    return padRange(range);
}

double yForValue(double value, const ValueRange& range, const QRectF& plot)
{
    return plot.bottom() - plot.height() * (value - range.min) / std::max(1e-12, range.max - range.min);
}

void paintPanelBackground(QPainter& p, const PanelGeometry& geometry)
{
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setPen(QPen(kPanelBorder, 1.0));
    p.setBrush(kPanel);
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
    p.setPen(kText);
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
    p.setPen(QPen(kGrid, 1.0));
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
    p.setPen(kTextMuted);
    p.drawText(QRectF(geometry.panel.left() + 8, geometry.plot.top() - 8, 38, 16),
               Qt::AlignRight | Qt::AlignVCenter, fmt(range.max, 2));
    p.drawText(QRectF(geometry.panel.left() + 8, geometry.plot.bottom() - 8, 38, 16),
               Qt::AlignRight | Qt::AlignVCenter, fmt(range.min, 2));
}

void paintPlotBorder(QPainter& p, const QRectF& plot)
{
    p.setRenderHint(QPainter::Antialiasing, false);
    p.setPen(QPen(kPanelBorder, 1.0));
    p.setBrush(Qt::NoBrush);
    p.drawRect(plot.adjusted(0.5, 0.5, -0.5, -0.5));
}

void paintTimeTicks(QPainter& p, const QRectF& plot, const TimeScale& time)
{
    p.setFont(monoFont(10));
    p.setPen(kTextMuted);
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
    p.setPen(kTextMuted);
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
               kSelection);
}

void paintTemporalCursor(QPainter& p, const PanelGeometry& geometry, const PaintContext& context)
{
    const double cursorX = context.time.xForFrame(context.currentFrame, geometry.plot);
    p.setPen(QPen(kCursor, 1.2, Qt::DashLine));
    p.drawLine(QPointF(cursorX + 0.5, geometry.plot.top()),
               QPointF(cursorX + 0.5, geometry.plot.bottom()));

    if (context.hasHover && geometry.plot.contains(context.hoverPos)) {
        p.setPen(QPen(kHover, 1.0, Qt::DotLine));
        p.drawLine(QPointF(context.hoverPos.x() + 0.5, geometry.plot.top()),
                   QPointF(context.hoverPos.x() + 0.5, geometry.plot.bottom()));
    }
}

class AbstractStripPanel {
public:
    virtual ~AbstractStripPanel() = default;
    virtual bool hasRevealBinding() const = 0;
    virtual model::SignalBinding revealBinding() const = 0;
    virtual QString tooltipLine(int) const { return {}; }
    virtual void paint(QPainter& p, const PanelGeometry& geometry, const PaintContext& context) const = 0;

    bool revealContains(const PanelGeometry& geometry, const QPoint& pos) const
    {
        return hasRevealBinding() && geometry.reveal.contains(pos);
    }

    virtual bool plotContains(const PanelGeometry&, const QPoint&) const
    {
        return false;
    }
};

class TemporalStripPanel final : public AbstractStripPanel {
public:
    explicit TemporalStripPanel(const StripStackWidget::Track& track)
        : track_(track)
    {}

    bool hasRevealBinding() const override { return track_.hasBinding; }
    model::SignalBinding revealBinding() const override { return track_.binding; }

    bool plotContains(const PanelGeometry& geometry, const QPoint& pos) const override
    {
        return geometry.plot.contains(pos);
    }

    QString tooltipLine(int frame) const override
    {
        const model::ChannelBuffer* buffer = track_.buffer;
        if (!buffer || frame < 0 || frame >= static_cast<int>(buffer->size()))
            return {};
        const std::size_t idx = static_cast<std::size_t>(frame);
        if (idx >= buffer->valid.size() || !buffer->valid[idx])
            return {};
        return QStringLiteral("<span style='color:%1'>■</span> %2: <b>%3 %4</b><br/>")
            .arg(track_.color.name(), buffer->label, fmt(buffer->values[idx], 3), buffer->unit);
    }

    void paint(QPainter& p, const PanelGeometry& geometry, const PaintContext& context) const override
    {
        paintPanelBackground(p, geometry);
        if (!track_.buffer)
            return;

        const auto& buffer = *track_.buffer;
        const ValueRange range = visibleTemporalRange(buffer, context.time);
        const bool hoverInPlot = context.hasHover && geometry.plot.contains(context.hoverPos);
        const int readoutFrame = hoverInPlot
                                     ? context.time.frameAt(context.hoverPos, geometry.plot)
                                     : context.currentFrame;
        QString readout;
        if (readoutFrame >= 0 && readoutFrame < static_cast<int>(buffer.size())) {
            const std::size_t idx = static_cast<std::size_t>(readoutFrame);
            if (idx < buffer.valid.size() && buffer.valid[idx])
                readout = QStringLiteral("f%1  %2").arg(readoutFrame + 1).arg(fmt(buffer.values[idx], 3));
        }

        paintHeader(p,
                    geometry,
                    HeaderText{QStringLiteral("%1 (%2)").arg(buffer.label, buffer.unit), readout},
                    170,
                    hoverInPlot ? kHover : kTextMuted,
                    track_.hasBinding,
                    context.hasHover && geometry.reveal.contains(context.hoverPos));
        paintGrid(p, geometry.plot);
        paintTimeTicks(p, geometry.plot, context.time);
        paintSelectedTimeRange(p, geometry, context);
        paintYAxisLabels(p, geometry, range);
        paintTrace(p, geometry, context, buffer, range);
        paintPlotBorder(p, geometry.plot);
        paintTemporalCursor(p, geometry, context);
    }

private:
    static bool isValidAt(const model::ChannelBuffer& buffer, int frame)
    {
        const std::size_t idx = static_cast<std::size_t>(frame);
        return idx < buffer.valid.size() && buffer.valid[idx] && std::isfinite(buffer.values[idx]);
    }

    void paintDirectTrace(QPainter& p,
                          const PanelGeometry& geometry,
                          const PaintContext& context,
                          const model::ChannelBuffer& buffer,
                          const ValueRange& range,
                          FrameWindow dataWindow,
                          int frameSpan) const
    {
        p.setRenderHint(QPainter::Antialiasing, true);
        QPointF previous;
        int previousFrame = -1;
        bool havePrevious = false;
        QPen tracePen(track_.color, 2.2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
        QColor markerFill = track_.color.lighter(115);

        for (int f = dataWindow.first; f <= dataWindow.last; ++f) {
            if (!isValidAt(buffer, f)) {
                havePrevious = false;
                previousFrame = -1;
                continue;
            }
            const std::size_t idx = static_cast<std::size_t>(f);
            const QPointF pt(context.time.xForFrame(f, geometry.plot),
                             yForValue(buffer.values[idx], range, geometry.plot));
            if (havePrevious && previousFrame + 1 == f) {
                p.setPen(tracePen);
                p.drawLine(previous, pt);
            }
            previous = pt;
            previousFrame = f;
            havePrevious = true;
        }

        const int markerStride = std::max(1, frameSpan / 220);
        p.setPen(QPen(kPanel, 1.0));
        p.setBrush(markerFill);
        for (int f = dataWindow.first; f <= dataWindow.last; f += markerStride) {
            if (!isValidAt(buffer, f))
                continue;
            const std::size_t idx = static_cast<std::size_t>(f);
            const QPointF pt(context.time.xForFrame(f, geometry.plot),
                             yForValue(buffer.values[idx], range, geometry.plot));
            p.drawEllipse(pt, 2.8, 2.8);
        }
    }

    void paintDecimatedTrace(QPainter& p,
                             const PanelGeometry& geometry,
                             const PaintContext& context,
                             const model::ChannelBuffer& buffer,
                             const ValueRange& range,
                             FrameWindow dataWindow,
                             int pxCount) const
    {
        QPolygonF segment;
        segment.reserve(pxCount + 1);
        auto flushSegment = [&]() {
            if (segment.size() >= 2) {
                p.setPen(QPen(track_.color, 1.8, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
                p.drawPolyline(segment);
            }
            segment.clear();
        };

        QColor spikeColor = track_.color;
        spikeColor.setAlpha(125);
        p.setPen(QPen(spikeColor, 1.0));

        for (int px = 0; px <= pxCount; ++px) {
            const double f0d = context.time.first
                               + (context.time.last - context.time.first) * px
                                     / static_cast<double>(pxCount);
            const double f1d = context.time.first
                               + (context.time.last - context.time.first) * (px + 1)
                                     / static_cast<double>(pxCount);
            int f0 = std::max(dataWindow.first, static_cast<int>(std::floor(f0d)));
            int f1 = std::min(dataWindow.last, std::max(f0, static_cast<int>(std::ceil(f1d))));
            if (f0 > dataWindow.last || f1 < dataWindow.first) {
                flushSegment();
                continue;
            }

            bool valid = false;
            double vMin = std::numeric_limits<double>::infinity();
            double vMax = -std::numeric_limits<double>::infinity();
            for (int f = f0; f <= f1; ++f) {
                if (!isValidAt(buffer, f))
                    continue;
                const double v = buffer.values[static_cast<std::size_t>(f)];
                vMin = std::min(vMin, v);
                vMax = std::max(vMax, v);
                valid = true;
            }
            if (!valid) {
                flushSegment();
                continue;
            }

            const double x = geometry.plot.left() + geometry.plot.width() * px / static_cast<double>(pxCount);
            const double yMin = yForValue(vMin, range, geometry.plot);
            const double yMax = yForValue(vMax, range, geometry.plot);
            if (std::abs(yMax - yMin) > 1.0) {
                p.setPen(QPen(spikeColor, 1.0));
                p.drawLine(QPointF(x, yMin), QPointF(x, yMax));
            }

            segment.append(QPointF(x, (yMin + yMax) * 0.5));
        }
        flushSegment();
    }

    void paintTrace(QPainter& p,
                    const PanelGeometry& geometry,
                    const PaintContext& context,
                    const model::ChannelBuffer& buffer,
                    const ValueRange& range) const
    {
        p.setRenderHint(QPainter::Antialiasing, true);
        if (!range.valid || buffer.empty())
            return;

        const FrameWindow dataWindow{
            std::max(0, context.time.first),
            std::min(static_cast<int>(buffer.lastFrame()), context.time.last),
        };
        if (!dataWindow.valid())
            return;

        p.save();

        const int frameSpan = std::max(1, context.time.last - context.time.first + 1);
        const int pxCount = std::max(1, static_cast<int>(std::ceil(geometry.plot.width())));
        const bool directSamples = frameSpan <= pxCount * 2;

        if (directSamples) {
            paintDirectTrace(p, geometry, context, buffer, range, dataWindow, frameSpan);
        } else {
            paintDecimatedTrace(p, geometry, context, buffer, range, dataWindow, pxCount);
        }

        if (context.currentFrame >= dataWindow.first
            && context.currentFrame <= dataWindow.last
            && isValidAt(buffer, context.currentFrame)) {
            const std::size_t idx = static_cast<std::size_t>(context.currentFrame);
            const QPointF pt(context.time.xForFrame(context.currentFrame, geometry.plot),
                             yForValue(buffer.values[idx], range, geometry.plot));
            p.setPen(QPen(kCanvas, 2.0));
            p.setBrush(track_.color);
            p.drawEllipse(pt, 4.0, 4.0);
        }

        p.restore();
    }

    const StripStackWidget::Track& track_;
};

class SpectrumStripPanel final : public AbstractStripPanel {
public:
    explicit SpectrumStripPanel(const StripStackWidget::SpectrumTrack& track)
        : track_(track)
    {}

    bool hasRevealBinding() const override { return track_.hasBinding; }
    model::SignalBinding revealBinding() const override { return track_.binding; }

    void paint(QPainter& p, const PanelGeometry& geometry, const PaintContext& context) const override
    {
        paintPanelBackground(p, geometry);
        const QString title = track_.yUnit.isEmpty()
                                  ? track_.label
                                  : QStringLiteral("%1 (%2)").arg(track_.label, track_.yUnit);
        paintHeader(p,
                    geometry,
                    HeaderText{title, track_.readout},
                    track_.readout.isEmpty() ? 12 : 170,
                    kTextMuted,
                    track_.hasBinding,
                    context.hasHover && geometry.reveal.contains(context.hoverPos));
        paintGrid(p, geometry.plot);

        double xMin = 0.0;
        double xMax = 1.0;
        ValueRange range = spectrumRange(&xMin, &xMax);
        paintSpectrumTicks(p, geometry.plot, xMin, xMax);
        paintYAxisLabels(p, geometry, range);
        paintTrace(p, geometry, range, xMin, xMax);
        paintPlotBorder(p, geometry.plot);
    }

private:
    ValueRange spectrumRange(double* xMin, double* xMax) const
    {
        ValueRange range;
        *xMin = 0.0;
        *xMax = 1.0;
        if (track_.points && !track_.points->empty()) {
            bool haveX = false;
            for (const QPointF& pt : *track_.points) {
                if (!std::isfinite(pt.x()) || !std::isfinite(pt.y()))
                    continue;
                if (!haveX) {
                    *xMin = *xMax = pt.x();
                    range.min = range.max = pt.y();
                    haveX = true;
                    range.valid = true;
                } else {
                    *xMin = std::min(*xMin, pt.x());
                    *xMax = std::max(*xMax, pt.x());
                    range.min = std::min(range.min, pt.y());
                    range.max = std::max(range.max, pt.y());
                }
            }
            *xMin = std::min(0.0, *xMin);
        }
        if (*xMax <= *xMin)
            *xMax = *xMin + 1.0;
        if (range.valid && std::abs(range.max - range.min) < 1e-12) {
            range.max += 1.0;
            range.min = std::min(0.0, range.min);
        } else if (range.valid) {
            const double pad = (range.max - range.min) * 0.08;
            range.max += pad;
            range.min = std::min(0.0, range.min - pad);
        }
        return range;
    }

    void paintTrace(QPainter& p,
                    const PanelGeometry& geometry,
                    const ValueRange& range,
                    double xMin,
                    double xMax) const
    {
        if (!range.valid || !track_.points || track_.points->size() < 2)
            return;

        QPolygonF polyline;
        polyline.reserve(track_.points->size());
        for (const QPointF& pt : *track_.points) {
            if (!std::isfinite(pt.x()) || !std::isfinite(pt.y()))
                continue;
            const double x = geometry.plot.left()
                             + geometry.plot.width() * (pt.x() - xMin) / (xMax - xMin);
            const double y = yForValue(pt.y(), range, geometry.plot);
            polyline.append(QPointF(x, y));
        }
        if (polyline.size() >= 2) {
            p.setRenderHint(QPainter::Antialiasing, true);
            p.setPen(QPen(track_.color, 1.8, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
            p.drawPolyline(polyline);
        }
    }

    const StripStackWidget::SpectrumTrack& track_;
};
}  // namespace

StripStackWidget::StripStackWidget(QWidget* parent)
    : QWidget(parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("StripStackWidget"));
    setMouseTracking(true);
    setMinimumHeight(kMinTrackHeight * 2 + kGap);
    setFocusPolicy(Qt::StrongFocus);
}

void StripStackWidget::setTracks(QVector<Track> tracks)
{
    ASSERT_THREAD(this);
    tracks_ = std::move(tracks);
    updateMinimumHeight();
    updateGeometry();
    update();
}

void StripStackWidget::setSpectrumTracks(QVector<SpectrumTrack> tracks)
{
    ASSERT_THREAD(this);
    spectrumTracks_ = std::move(tracks);
    updateMinimumHeight();
    updateGeometry();
    update();
}

void StripStackWidget::setTimeViewport(TimeViewportController* viewport)
{
    ASSERT_THREAD(this);
    if (viewport_)
        disconnect(viewport_, nullptr, this, nullptr);
    viewport_ = viewport;
    if (viewport_) {
        ACONNECT(viewport_.data(), &TimeViewportController::currentFrameChanged, this, [this](int frame) {
            currentFrame_ = frame;
            update();
        });
        ACONNECT(viewport_.data(), &TimeViewportController::visibleRangeChanged, this, [this](int, int) {
            update();
        });
        ACONNECT(viewport_.data(), &TimeViewportController::selectedRangeChanged, this, [this](int, int, bool) {
            update();
        });
    }
    update();
}

void StripStackWidget::setCurrentFrame(int frame)
{
    ASSERT_THREAD(this);
    currentFrame_ = frame;
    update();
}

QRectF StripStackWidget::trackRect(int index) const
{
    return panelRectForIndex(StackGeometry{size(), panelCount()}, index);
}

QRectF StripStackWidget::plotRect(const QRectF& r) const
{
    return plotRectForPanel(r);
}

QRectF StripStackWidget::revealRect(const QRectF& r) const
{
    return revealRectForPanel(r);
}

int StripStackWidget::panelCount() const
{
    return static_cast<int>(tracks_.size() + spectrumTracks_.size());
}

void StripStackWidget::updateMinimumHeight()
{
    const int n = std::max(2, panelCount());
    setMinimumHeight(kMinTrackHeight * n + kGap * (n + 1));
}

bool StripStackWidget::timePlotContains(const QPoint& pos) const
{
    const StackGeometry stack{size(), panelCount()};
    for (int i = 0; i < tracks_.size(); ++i) {
        const TemporalStripPanel panel(tracks_[i]);
        if (panel.plotContains(panelGeometryForIndex(stack, i), pos))
            return true;
    }
    return false;
}

bool StripStackWidget::revealAt(const QPoint& pos, model::SignalBinding* binding) const
{
    const StackGeometry stack{size(), panelCount()};
    for (int i = 0; i < tracks_.size(); ++i) {
        const TemporalStripPanel panel(tracks_[i]);
        if (panel.revealContains(panelGeometryForIndex(stack, i), pos)) {
            if (binding)
                *binding = panel.revealBinding();
            return true;
        }
    }
    for (int i = 0; i < spectrumTracks_.size(); ++i) {
        const SpectrumStripPanel panel(spectrumTracks_[i]);
        const int panelIndex = static_cast<int>(tracks_.size()) + i;
        if (panel.revealContains(panelGeometryForIndex(stack, panelIndex), pos)) {
            if (binding)
                *binding = panel.revealBinding();
            return true;
        }
    }
    return false;
}

int StripStackWidget::frameAt(const QPoint& pos) const
{
    const TimeScale time{viewport_ ? viewport_->visibleStart() : 0,
                         viewport_ ? viewport_->visibleEnd() : std::max(0, currentFrame_)};
    const QRectF plot = plotRect(trackRect(0));
    return time.frameAt(pos, plot);
}

QString StripStackWidget::tooltipText(int frame) const
{
    QString text = QStringLiteral("<b>frame %1</b><br/>").arg(frame + 1);
    for (const auto& tr : tracks_) {
        const TemporalStripPanel panel(tr);
        text += panel.tooltipLine(frame);
    }
    return text;
}

void StripStackWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing, true);
    p.fillRect(rect(), kCanvas);

    if (tracks_.empty() && spectrumTracks_.empty()) {
        p.setFont(uiFont(13));
        p.setPen(kTextMuted);
        p.drawText(rect(), Qt::AlignCenter, QStringLiteral("Select atoms to start a trajectory strip"));
        return;
    }

    const PaintContext context{
        viewport_.data(),
        currentFrame_,
        hasHover_,
        hoverPos_,
        TimeScale{viewport_ ? viewport_->visibleStart() : 0,
                  viewport_ ? viewport_->visibleEnd() : std::max(0, currentFrame_)},
    };
    const StackGeometry stack{size(), panelCount()};

    for (int i = 0; i < tracks_.size(); ++i) {
        const TemporalStripPanel panel(tracks_[i]);
        panel.paint(p, panelGeometryForIndex(stack, i), context);
    }

    for (int i = 0; i < spectrumTracks_.size(); ++i) {
        const SpectrumStripPanel panel(spectrumTracks_[i]);
        const int panelIndex = static_cast<int>(tracks_.size()) + i;
        panel.paint(p, panelGeometryForIndex(stack, panelIndex), context);
    }
}

void StripStackWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        model::SignalBinding binding;
        if (revealAt(event->pos(), &binding)) {
            emit revealRequested(binding);
            event->accept();
            return;
        }
    }

    if (event->button() != Qt::LeftButton || !viewport_ || !timePlotContains(event->pos())) {
        QWidget::mousePressEvent(event);
        return;
    }
    selecting_ = true;
    dragSelecting_ = false;
    pressWasFollowing_ = viewport_->followsPlayhead();
    pressPos_ = event->pos();
    selectionAnchor_ = frameAt(event->pos());
    event->accept();
}

void StripStackWidget::mouseMoveEvent(QMouseEvent* event)
{
    hasHover_ = true;
    hoverPos_ = event->pos();
    const int frame = frameAt(event->pos());
    model::SignalBinding revealBinding;
    const bool overReveal = revealAt(event->pos(), &revealBinding);
    if (overReveal) {
        setCursor(Qt::PointingHandCursor);
        QToolTip::showText(event->globalPosition().toPoint(), QStringLiteral("Reveal in 3-D scene"), this);
    } else {
        unsetCursor();
    }
    if (selecting_ && viewport_) {
        if (!dragSelecting_ &&
            (event->pos() - pressPos_).manhattanLength() >= QApplication::startDragDistance()) {
            dragSelecting_ = true;
            viewport_->setFollowPlayhead(false);
            viewport_->setSelectedRange(selectionAnchor_, frame);
        } else if (dragSelecting_) {
            viewport_->setSelectedRange(selectionAnchor_, frame);
        }
    } else if (!overReveal && !tracks_.empty() && timePlotContains(event->pos())) {
        QToolTip::showText(event->globalPosition().toPoint(), tooltipText(frame), this);
    } else if (!overReveal) {
        QToolTip::hideText();
    }
    update();
}

void StripStackWidget::mouseReleaseEvent(QMouseEvent* event)
{
    if (!selecting_ || event->button() != Qt::LeftButton || !viewport_) {
        QWidget::mouseReleaseEvent(event);
        return;
    }
    selecting_ = false;
    const int frame = frameAt(event->pos());
    if (!dragSelecting_ || std::abs(frame - selectionAnchor_) <= 1) {
        viewport_->clearSelectedRange();
        viewport_->requestPlaybackFrame(frame);
        if (pressWasFollowing_)
            viewport_->setFollowPlayhead(true);
    } else {
        viewport_->setSelectedRange(selectionAnchor_, frame);
        viewport_->requestPlaybackFrame(std::min(selectionAnchor_, frame));
    }
    dragSelecting_ = false;
    event->accept();
}

void StripStackWidget::mouseDoubleClickEvent(QMouseEvent* event)
{
    if (viewport_ && timePlotContains(event->pos())) {
        viewport_->clearSelectedRange();
        viewport_->setFollowPlayhead(true);
        viewport_->requestPlaybackFrame(frameAt(event->pos()));
        event->accept();
        return;
    }
    QWidget::mouseDoubleClickEvent(event);
}

void StripStackWidget::leaveEvent(QEvent*)
{
    hasHover_ = false;
    unsetCursor();
    QToolTip::hideText();
    update();
}

void StripStackWidget::wheelEvent(QWheelEvent* event)
{
    if (!viewport_ || !timePlotContains(event->position().toPoint())) {
        QWidget::wheelEvent(event);
        return;
    }
    const int steps = event->angleDelta().y() / 120;
    if (steps == 0)
        return;

    if (event->modifiers() & Qt::ShiftModifier) {
        const int span = viewport_->visibleEnd() - viewport_->visibleStart() + 1;
        viewport_->panFrames(-steps * std::max(1, span / 8));
    } else {
        const double factor = steps > 0 ? 0.75 : 1.333333333333;
        viewport_->zoomAround(frameAt(event->position().toPoint()), factor);
    }
    event->accept();
}

}  // namespace h5reader::app
