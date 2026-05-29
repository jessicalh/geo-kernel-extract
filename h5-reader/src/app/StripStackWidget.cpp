#include "StripStackWidget.h"

#include "AbstractStripPanel.h"
#include "TimeViewportController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/StripChartChannel.h"

#include <QApplication>
#include <QEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPolygonF>
#include <QPointF>
#include <QToolTip>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace h5reader::app {

namespace {

// Theme + fonts + geometry types + paint helpers all live in
// AbstractStripPanel.h now. What remains here is panel-specific:
// the two concrete panels that ride Track / SpectrumTrack data, plus
// the temporal-axis ValueRange helper that walks a ChannelBuffer.

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

class TemporalStripPanel final : public AbstractStripPanel {
public:
    explicit TemporalStripPanel(StripStackWidget::Track track)
        : track_(std::move(track))
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
            .arg(track_.color.name(), buffer->label, fmtValue(buffer->values[idx], 3), buffer->unit);
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
                readout = QStringLiteral("f%1  %2").arg(readoutFrame + 1).arg(fmtValue(buffer.values[idx], 3));
        }

        paintHeader(p,
                    geometry,
                    HeaderText{QStringLiteral("%1 (%2)").arg(buffer.label, buffer.unit), readout},
                    170,
                    hoverInPlot ? kStripHover : kStripTextMuted,
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
        p.setPen(QPen(kStripPanel, 1.0));
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
            p.setPen(QPen(kStripCanvas, 2.0));
            p.setBrush(track_.color);
            p.drawEllipse(pt, 4.0, 4.0);
        }

        p.restore();
    }

    // L-2b (2026-05-29): owned by value; was const-ref. The panel
    // outlives the QVector<Track> the controller built it from, so a
    // setTracks() reallocation can't dangle the reference.
    StripStackWidget::Track track_;
};

class SpectrumStripPanel final : public AbstractStripPanel {
public:
    explicit SpectrumStripPanel(StripStackWidget::SpectrumTrack track)
        : track_(std::move(track))
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
                    kStripTextMuted,
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

    // L-2b: same lifetime fix as TemporalStripPanel — owned by value.
    StripStackWidget::SpectrumTrack track_;
};

}  // namespace

StripStackWidget::StripStackWidget(QWidget* parent)
    : QWidget(parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("StripStackWidget"));
    setMouseTracking(true);
    setMinimumHeight(kStripMinPanelHeight * 2 + kStripPanelGap);
    setFocusPolicy(Qt::StrongFocus);
}

void StripStackWidget::setTracks(QVector<Track> tracks)
{
    ASSERT_THREAD(this);
    // Replace the [0..n_temporal_) section only; keep the spectrum and
    // owned tails intact in their existing order.
    std::vector<std::unique_ptr<AbstractStripPanel>> next;
    const std::size_t newCount = static_cast<std::size_t>(tracks.size());
    next.reserve(newCount + (panels_.size() - n_temporal_));
    for (Track& t : tracks)
        next.push_back(std::make_unique<TemporalStripPanel>(std::move(t)));
    for (std::size_t i = n_temporal_; i < panels_.size(); ++i)
        next.push_back(std::move(panels_[i]));
    panels_ = std::move(next);
    n_temporal_ = newCount;
    updateMinimumHeight();
    updateGeometry();
    update();
}

void StripStackWidget::setSpectrumTracks(QVector<SpectrumTrack> tracks)
{
    ASSERT_THREAD(this);
    // Replace the [n_temporal_..n_temporal_+n_spectrum_) section only.
    std::vector<std::unique_ptr<AbstractStripPanel>> next;
    const std::size_t newCount = static_cast<std::size_t>(tracks.size());
    const std::size_t ownedStart = n_temporal_ + n_spectrum_;
    next.reserve(n_temporal_ + newCount + (panels_.size() - ownedStart));
    for (std::size_t i = 0; i < n_temporal_; ++i)
        next.push_back(std::move(panels_[i]));
    for (SpectrumTrack& t : tracks)
        next.push_back(std::make_unique<SpectrumStripPanel>(std::move(t)));
    for (std::size_t i = ownedStart; i < panels_.size(); ++i)
        next.push_back(std::move(panels_[i]));
    panels_ = std::move(next);
    n_spectrum_ = newCount;
    updateMinimumHeight();
    updateGeometry();
    update();
}

void StripStackWidget::setOwnedPanels(std::vector<std::unique_ptr<AbstractStripPanel>> panels)
{
    ASSERT_THREAD(this);
    // Replace the trailing section; preserve temporal + spectrum heads.
    std::vector<std::unique_ptr<AbstractStripPanel>> next;
    const std::size_t headSize = n_temporal_ + n_spectrum_;
    next.reserve(headSize + panels.size());
    for (std::size_t i = 0; i < headSize && i < panels_.size(); ++i)
        next.push_back(std::move(panels_[i]));
    for (auto& p : panels)
        next.push_back(std::move(p));
    panels_ = std::move(next);
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
    return static_cast<int>(panels_.size());
}

void StripStackWidget::updateMinimumHeight()
{
    const int n = std::max(2, panelCount());
    setMinimumHeight(kStripMinPanelHeight * n + kStripPanelGap * (n + 1));
}

bool StripStackWidget::timePlotContains(const QPoint& pos) const
{
    // Time-axis drag/select gestures only apply to temporal panels.
    // n_temporal_ tracks the front section of panels_ that hold
    // TemporalStripPanel instances.
    const StackGeometry stack{size(), panelCount()};
    for (std::size_t i = 0; i < n_temporal_ && i < panels_.size(); ++i) {
        const PanelGeometry geom =
            panelGeometryForIndex(stack, static_cast<int>(i), panels_[i]->preferredAspect());
        if (panels_[i]->plotContains(geom, pos))
            return true;
    }
    return false;
}

bool StripStackWidget::revealAt(const QPoint& pos, model::SignalBinding* binding) const
{
    const StackGeometry stack{size(), panelCount()};
    for (std::size_t i = 0; i < panels_.size(); ++i) {
        const PanelGeometry geom =
            panelGeometryForIndex(stack, static_cast<int>(i), panels_[i]->preferredAspect());
        if (panels_[i]->revealContains(geom, pos)) {
            if (binding)
                *binding = panels_[i]->revealBinding();
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
    // Temporal tooltip composer — only temporal panels contribute a
    // per-frame value line. Other panel types return an empty string
    // from tooltipLine() so they're filtered naturally; restricting the
    // walk to [0..n_temporal_) makes the intent explicit.
    QString text = QStringLiteral("<b>frame %1</b><br/>").arg(frame + 1);
    for (std::size_t i = 0; i < n_temporal_ && i < panels_.size(); ++i)
        text += panels_[i]->tooltipLine(frame);
    return text;
}

void StripStackWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing, true);
    p.fillRect(rect(), kStripCanvas);

    if (panels_.empty()) {
        p.setFont(uiFont(13));
        p.setPen(kStripTextMuted);
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

    // Unified single loop — panels_ holds temporal, spectrum, and
    // owned in that order. Each panel's preferredAspect() lets the
    // geometry helper letterbox aspect-locked panels (chord) within
    // their assigned rect.
    for (std::size_t i = 0; i < panels_.size(); ++i) {
        const PanelGeometry geom =
            panelGeometryForIndex(stack, static_cast<int>(i), panels_[i]->preferredAspect());
        panels_[i]->paint(p, geom, context);
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
        // Forward to non-temporal panels (spectrum, sequence-bar,
        // chord, etc.) for in-plot click handling. Temporal panels'
        // default mousePressInPlot is a no-op and they own the
        // drag-select-range gesture below via timePlotContains, so
        // restricting the forward loop to the spectrum+owned section
        // keeps the two interaction models from racing.
        const StackGeometry stack{size(), panelCount()};
        for (std::size_t i = n_temporal_; i < panels_.size(); ++i) {
            auto& panel = *panels_[i];
            const PanelGeometry geom =
                panelGeometryForIndex(stack, static_cast<int>(i), panel.preferredAspect());
            if (geom.plot.contains(event->pos())) {
                if (auto pressBinding = panel.mousePressInPlot(event, geom); pressBinding) {
                    emit revealRequested(*pressBinding);
                    event->accept();
                    return;
                }
                if (event->isAccepted())
                    return;
            }
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
        // Forward to non-temporal panels for in-plot hover handling.
        const StackGeometry stack{size(), panelCount()};
        for (std::size_t i = n_temporal_; i < panels_.size(); ++i) {
            auto& panel = *panels_[i];
            const PanelGeometry geom =
                panelGeometryForIndex(stack, static_cast<int>(i), panel.preferredAspect());
            if (geom.plot.contains(event->pos()))
                panel.mouseMoveInPlot(event, geom);
        }
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
    } else if (!overReveal && n_temporal_ > 0 && timePlotContains(event->pos())) {
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
