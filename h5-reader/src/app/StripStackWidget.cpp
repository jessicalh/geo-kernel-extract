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
    const int n = std::max(1, panelCount());
    const double h = (height() - (n + 1) * kGap) / static_cast<double>(n);
    return QRectF(kGap, kGap + index * (h + kGap), width() - 2 * kGap, h);
}

QRectF StripStackWidget::plotRect(const QRectF& r) const
{
    if (r.height() < 96.0)
        return r.adjusted(48, 28, -10, -20);
    return r.adjusted(52, 30, -12, -24);
}

QRectF StripStackWidget::revealRect(const QRectF& r) const
{
    return QRectF(r.left() + 8.0, r.top() + 6.0, 18.0, 18.0);
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
    for (int i = 0; i < tracks_.size(); ++i) {
        if (plotRect(trackRect(i)).contains(pos))
            return true;
    }
    return false;
}

bool StripStackWidget::revealAt(const QPoint& pos, model::SignalBinding* binding) const
{
    for (int i = 0; i < tracks_.size(); ++i) {
        const auto& tr = tracks_[i];
        if (tr.hasBinding && revealRect(trackRect(i)).contains(pos)) {
            if (binding)
                *binding = tr.binding;
            return true;
        }
    }
    for (int i = 0; i < spectrumTracks_.size(); ++i) {
        const auto& tr = spectrumTracks_[i];
        const int panelIndex = static_cast<int>(tracks_.size()) + i;
        if (tr.hasBinding && revealRect(trackRect(panelIndex)).contains(pos)) {
            if (binding)
                *binding = tr.binding;
            return true;
        }
    }
    return false;
}

StripStackWidget::Range StripStackWidget::visibleRange(const model::ChannelBuffer& buffer) const
{
    Range range;
    if (buffer.empty())
        return range;

    const int first = viewport_ ? viewport_->visibleStart() : 0;
    const int last = viewport_ ? viewport_->visibleEnd() : static_cast<int>(buffer.lastFrame());
    const int lo = std::max(0, first);
    const int hi = std::min(static_cast<int>(buffer.lastFrame()), last);
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

int StripStackWidget::frameAt(const QPoint& pos) const
{
    const int first = viewport_ ? viewport_->visibleStart() : 0;
    const int last = viewport_ ? viewport_->visibleEnd() : std::max(0, currentFrame_);
    const QRectF plot = plotRect(trackRect(0));
    const double xNorm = std::clamp((pos.x() - plot.left()) / std::max(1.0, plot.width()), 0.0, 1.0);
    return first + static_cast<int>(std::round(xNorm * std::max(0, last - first)));
}

double StripStackWidget::xForFrame(int frame, const QRectF& plot) const
{
    const int first = viewport_ ? viewport_->visibleStart() : 0;
    const int last = viewport_ ? viewport_->visibleEnd() : std::max(0, currentFrame_);
    if (last <= first)
        return plot.left() + plot.width() * 0.5;
    return plot.left() + plot.width() * (frame - first) / std::max(1, last - first);
}

double StripStackWidget::yForValue(double value, const Range& range, const QRectF& plot) const
{
    return plot.bottom() - plot.height() * (value - range.min) / std::max(1e-12, range.max - range.min);
}

QString StripStackWidget::tooltipText(int frame) const
{
    QString text = QStringLiteral("<b>frame %1</b><br/>").arg(frame + 1);
    for (const auto& tr : tracks_) {
        if (!tr.buffer || frame < 0 || frame >= static_cast<int>(tr.buffer->size()))
            continue;
        const std::size_t idx = static_cast<std::size_t>(frame);
        if (idx >= tr.buffer->valid.size() || !tr.buffer->valid[idx])
            continue;
        text += QStringLiteral("<span style='color:%1'>■</span> %2: <b>%3 %4</b><br/>")
                    .arg(tr.color.name(), tr.buffer->label, fmt(tr.buffer->values[idx], 3), tr.buffer->unit);
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

    for (int i = 0; i < tracks_.size(); ++i) {
        const auto& tr = tracks_[i];
        const QRectF panel = trackRect(i);
        const QRectF plot = plotRect(panel);

        p.setPen(QPen(kPanelBorder, 1.0));
        p.setBrush(kPanel);
        p.drawRoundedRect(panel.adjusted(0.5, 0.5, -0.5, -0.5), 5, 5);

        if (!tr.buffer) {
            continue;
        }

        const auto& buffer = *tr.buffer;
        const Range range = visibleRange(buffer);

        const QRectF targetRect = revealRect(panel);
        if (tr.hasBinding)
            drawRevealButton(p, targetRect, hasHover_ && targetRect.contains(hoverPos_));

        p.setFont(uiFont(13, QFont::Medium));
        p.setPen(kText);
        const QString title = QStringLiteral("%1 (%2)").arg(buffer.label, buffer.unit);
        const int readoutWidth = 170;
        const double titleLeft = panel.left() + (tr.hasBinding ? 34.0 : 10.0);
        const QRectF titleRect(titleLeft, panel.top() + 6.0,
                               std::max(20.0, panel.right() - titleLeft - readoutWidth - 10.0),
                               22.0);
        p.drawText(titleRect,
                   Qt::AlignLeft | Qt::AlignVCenter,
                   QFontMetricsF(p.font()).elidedText(title, Qt::ElideMiddle,
                                                      static_cast<int>(titleRect.width())));

        const bool hoverInTrack = hasHover_ && plot.contains(hoverPos_);
        const int readoutFrame = hoverInTrack ? frameAt(hoverPos_) : currentFrame_;
        if (readoutFrame >= 0 && readoutFrame < static_cast<int>(buffer.size())) {
            const std::size_t idx = static_cast<std::size_t>(readoutFrame);
            if (idx < buffer.valid.size() && buffer.valid[idx]) {
                p.setFont(monoFont(12));
                p.setPen(hoverInTrack ? kHover : kTextMuted);
                p.drawText(panel.adjusted(10, 6, -10, -panel.height() + 28),
                           Qt::AlignRight | Qt::AlignVCenter,
                           QStringLiteral("f%1  %2").arg(readoutFrame + 1).arg(fmt(buffer.values[idx], 3)));
            }
        }

        p.setRenderHint(QPainter::Antialiasing, false);
        p.setPen(QPen(kGrid, 1.0));
        const int first = viewport_ ? viewport_->visibleStart() : 0;
        const int last = viewport_ ? viewport_->visibleEnd() : std::max(0, currentFrame_);
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

        p.setFont(monoFont(10));
        p.setPen(kTextMuted);
        for (int g = 0; g <= 4; ++g) {
            const double frac = g / 4.0;
            const double x = plot.left() + plot.width() * frac;
            const int tickFrame = first + static_cast<int>(std::round((last - first) * frac));
            const QRectF labelRect(x - 28.0, plot.bottom() + 4.0, 56.0, 16.0);
            p.drawText(labelRect, Qt::AlignHCenter | Qt::AlignVCenter,
                       QStringLiteral("f%1").arg(tickFrame + 1));
        }

        if (viewport_ && viewport_->hasSelectedRange()) {
            const double a = xForFrame(viewport_->selectedStart(), plot);
            const double b = xForFrame(viewport_->selectedEnd(), plot);
            p.fillRect(QRectF(QPointF(std::min(a, b), plot.top()),
                              QPointF(std::max(a, b), plot.bottom())),
                       kSelection);
        }

        p.setFont(monoFont(10));
        p.setPen(kTextMuted);
        if (range.valid) {
            p.drawText(QRectF(panel.left() + 8, plot.top() - 8, 38, 16),
                       Qt::AlignRight | Qt::AlignVCenter, fmt(range.max, 2));
            p.drawText(QRectF(panel.left() + 8, plot.bottom() - 8, 38, 16),
                       Qt::AlignRight | Qt::AlignVCenter, fmt(range.min, 2));
        }

        p.setRenderHint(QPainter::Antialiasing, true);
        if (range.valid && !buffer.empty()) {
            const int dataLo = std::max(0, first);
            const int dataHi = std::min(static_cast<int>(buffer.lastFrame()), last);

            auto isValidAt = [&buffer](int frame) {
                const std::size_t idx = static_cast<std::size_t>(frame);
                return idx < buffer.valid.size() && buffer.valid[idx] && std::isfinite(buffer.values[idx]);
            };

            if (dataLo <= dataHi) {
                p.save();

                const int frameSpan = std::max(1, last - first + 1);
                const int pxCount = std::max(1, static_cast<int>(std::ceil(plot.width())));
                const bool directSamples = frameSpan <= pxCount * 2;

                if (directSamples) {
                    p.setRenderHint(QPainter::Antialiasing, true);
                    QPointF previous;
                    int previousFrame = -1;
                    bool havePrevious = false;
                    QPen tracePen(tr.color, 2.2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
                    QColor markerFill = tr.color.lighter(115);

                    for (int f = dataLo; f <= dataHi; ++f) {
                        if (!isValidAt(f)) {
                            havePrevious = false;
                            previousFrame = -1;
                            continue;
                        }
                        const std::size_t idx = static_cast<std::size_t>(f);
                        const QPointF pt(xForFrame(f, plot), yForValue(buffer.values[idx], range, plot));
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
                    for (int f = dataLo; f <= dataHi; f += markerStride) {
                        if (!isValidAt(f))
                            continue;
                        const std::size_t idx = static_cast<std::size_t>(f);
                        const QPointF pt(xForFrame(f, plot), yForValue(buffer.values[idx], range, plot));
                        p.drawEllipse(pt, 2.8, 2.8);
                    }
                } else {
                    QPolygonF segment;
                    segment.reserve(pxCount + 1);
                    auto flushSegment = [&]() {
                        if (segment.size() >= 2) {
                            p.setPen(QPen(tr.color, 1.8, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
                            p.drawPolyline(segment);
                        }
                        segment.clear();
                    };

                    QColor spikeColor = tr.color;
                    spikeColor.setAlpha(125);
                    p.setPen(QPen(spikeColor, 1.0));

                    for (int px = 0; px <= pxCount; ++px) {
                        const double f0d = first + (last - first) * px / static_cast<double>(pxCount);
                        const double f1d = first + (last - first) * (px + 1) / static_cast<double>(pxCount);
                        int f0 = std::max(dataLo, static_cast<int>(std::floor(f0d)));
                        int f1 = std::min(dataHi, std::max(f0, static_cast<int>(std::ceil(f1d))));
                        if (f0 > dataHi || f1 < dataLo) {
                            flushSegment();
                            continue;
                        }

                        bool valid = false;
                        double vMin = std::numeric_limits<double>::infinity();
                        double vMax = -std::numeric_limits<double>::infinity();
                        for (int f = f0; f <= f1; ++f) {
                            if (!isValidAt(f))
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

                        const double x = plot.left() + plot.width() * px / static_cast<double>(pxCount);
                        const double yMin = yForValue(vMin, range, plot);
                        const double yMax = yForValue(vMax, range, plot);
                        if (std::abs(yMax - yMin) > 1.0) {
                            p.setPen(QPen(spikeColor, 1.0));
                            p.drawLine(QPointF(x, yMin), QPointF(x, yMax));
                        }

                        const QPointF mid(x, (yMin + yMax) * 0.5);
                        segment.append(mid);
                    }
                    flushSegment();
                }

                if (currentFrame_ >= dataLo && currentFrame_ <= dataHi && isValidAt(currentFrame_)) {
                    const std::size_t idx = static_cast<std::size_t>(currentFrame_);
                    const QPointF pt(xForFrame(currentFrame_, plot), yForValue(buffer.values[idx], range, plot));
                    p.setPen(QPen(kCanvas, 2.0));
                    p.setBrush(tr.color);
                    p.drawEllipse(pt, 4.0, 4.0);
                }

                p.restore();
            }
        }

        p.setRenderHint(QPainter::Antialiasing, false);
        p.setPen(QPen(kPanelBorder, 1.0));
        p.setBrush(Qt::NoBrush);
        p.drawRect(plot.adjusted(0.5, 0.5, -0.5, -0.5));

        const double cursorX = xForFrame(currentFrame_, plot);
        p.setPen(QPen(kCursor, 1.2, Qt::DashLine));
        p.drawLine(QPointF(cursorX + 0.5, plot.top()), QPointF(cursorX + 0.5, plot.bottom()));

        if (hasHover_ && plot.contains(hoverPos_)) {
            p.setPen(QPen(kHover, 1.0, Qt::DotLine));
            p.drawLine(QPointF(hoverPos_.x() + 0.5, plot.top()),
                       QPointF(hoverPos_.x() + 0.5, plot.bottom()));
        }
    }

    for (int i = 0; i < spectrumTracks_.size(); ++i) {
        const auto& tr = spectrumTracks_[i];
        const QRectF panel = trackRect(static_cast<int>(tracks_.size()) + i);
        const QRectF plot = plotRect(panel);

        p.setRenderHint(QPainter::Antialiasing, true);
        p.setPen(QPen(kPanelBorder, 1.0));
        p.setBrush(kPanel);
        p.drawRoundedRect(panel.adjusted(0.5, 0.5, -0.5, -0.5), 5, 5);

        p.setFont(uiFont(13, QFont::Medium));
        p.setPen(kText);
        const QString title = tr.yUnit.isEmpty()
                                  ? tr.label
                                  : QStringLiteral("%1 (%2)").arg(tr.label, tr.yUnit);
        const int readoutWidth = tr.readout.isEmpty() ? 12 : 170;
        const QRectF targetRect = revealRect(panel);
        if (tr.hasBinding)
            drawRevealButton(p, targetRect, hasHover_ && targetRect.contains(hoverPos_));
        const double titleLeft = panel.left() + (tr.hasBinding ? 34.0 : 10.0);
        const QRectF titleRect(titleLeft, panel.top() + 6,
                               std::max(20.0, panel.right() - titleLeft - readoutWidth - 10.0), 22);
        p.drawText(titleRect,
                   Qt::AlignLeft | Qt::AlignVCenter,
                   QFontMetricsF(p.font()).elidedText(title, Qt::ElideMiddle,
                                                      static_cast<int>(titleRect.width())));
        if (!tr.readout.isEmpty()) {
            p.setFont(monoFont(12));
            p.setPen(kTextMuted);
            p.drawText(panel.adjusted(10, 6, -10, -panel.height() + 28),
                       Qt::AlignRight | Qt::AlignVCenter,
                       tr.readout);
        }

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

        double xMin = 0.0;
        double xMax = 1.0;
        Range range;
        if (tr.points && !tr.points->empty()) {
            bool haveX = false;
            for (const QPointF& pt : *tr.points) {
                if (!std::isfinite(pt.x()) || !std::isfinite(pt.y()))
                    continue;
                if (!haveX) {
                    xMin = xMax = pt.x();
                    range.min = range.max = pt.y();
                    haveX = true;
                    range.valid = true;
                } else {
                    xMin = std::min(xMin, pt.x());
                    xMax = std::max(xMax, pt.x());
                    range.min = std::min(range.min, pt.y());
                    range.max = std::max(range.max, pt.y());
                }
            }
            xMin = std::min(0.0, xMin);
        }
        if (xMax <= xMin)
            xMax = xMin + 1.0;
        if (range.valid && std::abs(range.max - range.min) < 1e-12) {
            range.max += 1.0;
            range.min = std::min(0.0, range.min);
        } else if (range.valid) {
            const double pad = (range.max - range.min) * 0.08;
            range.max += pad;
            range.min = std::min(0.0, range.min - pad);
        }

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
        if (range.valid) {
            p.drawText(QRectF(panel.left() + 8, plot.top() - 8, 38, 16),
                       Qt::AlignRight | Qt::AlignVCenter, fmt(range.max, 2));
            p.drawText(QRectF(panel.left() + 8, plot.bottom() - 8, 38, 16),
                       Qt::AlignRight | Qt::AlignVCenter, fmt(range.min, 2));
        }

        if (range.valid && tr.points && tr.points->size() >= 2) {
            QPolygonF polyline;
            polyline.reserve(tr.points->size());
            for (const QPointF& pt : *tr.points) {
                if (!std::isfinite(pt.x()) || !std::isfinite(pt.y()))
                    continue;
                const double x = plot.left() + plot.width() * (pt.x() - xMin) / (xMax - xMin);
                const double y = yForValue(pt.y(), range, plot);
                polyline.append(QPointF(x, y));
            }
            if (polyline.size() >= 2) {
                p.setRenderHint(QPainter::Antialiasing, true);
                p.setPen(QPen(tr.color, 1.8, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
                p.drawPolyline(polyline);
            }
        }

        p.setRenderHint(QPainter::Antialiasing, false);
        p.setPen(QPen(kPanelBorder, 1.0));
        p.setBrush(Qt::NoBrush);
        p.drawRect(plot.adjusted(0.5, 0.5, -0.5, -0.5));
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
