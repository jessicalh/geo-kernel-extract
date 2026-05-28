#include "TimeViewportController.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QLoggingCategory>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cTimeViewport, "h5reader.timeviewport")
}

TimeViewportController::TimeViewportController(int frameCount, QObject* parent)
    : QObject(parent),
      frameCount_(std::max(1, frameCount))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("TimeViewportController"));

    windowFrames_ = clampSpan(windowFrames_);
    updateFollowWindow();

    qCInfo(cTimeViewport).noquote()
        << "created | frames=" << frameCount_
        << "| window=" << windowFrames_
        << "| visible=" << visibleStart_ << ".." << visibleEnd_;
}

int TimeViewportController::clampFrame(int frame) const
{
    return std::clamp(frame, 0, frameCount_ - 1);
}

int TimeViewportController::clampSpan(int frames) const
{
    return std::clamp(frames, 1, frameCount_);
}

void TimeViewportController::setCurrentFrame(int frame)
{
    ASSERT_THREAD(this);
    const int clamped = clampFrame(frame);
    if (clamped != currentFrame_) {
        currentFrame_ = clamped;
        emit currentFrameChanged(currentFrame_);
    }
    if (fitsCollected_)
        updateCollectedWindow();
    else if (followsPlayhead_)
        updateFollowWindow();
}

void TimeViewportController::setWindowFrames(int frames)
{
    ASSERT_THREAD(this);
    const int clamped = clampSpan(frames);
    if (clamped == windowFrames_)
        return;
    windowFrames_ = clamped;
    emit windowFramesChanged(windowFrames_);
    if (followsPlayhead_)
        updateFollowWindow();
}

void TimeViewportController::setFollowPlayhead(bool follow)
{
    ASSERT_THREAD(this);
    if (follow)
        fitsCollected_ = false;
    if (follow == followsPlayhead_)
        return;
    followsPlayhead_ = follow;
    emit followPlayheadChanged(followsPlayhead_);
    if (followsPlayhead_)
        updateFollowWindow();
}

void TimeViewportController::followPlayhead()
{
    setFollowPlayhead(true);
}

void TimeViewportController::fitCollectedRange()
{
    ASSERT_THREAD(this);
    const bool followChanged = followsPlayhead_;
    followsPlayhead_ = false;
    fitsCollected_ = true;
    if (followChanged)
        emit followPlayheadChanged(false);
    updateCollectedWindow();
}

void TimeViewportController::fitFullRange()
{
    ASSERT_THREAD(this);
    fitsCollected_ = false;
    setFollowPlayhead(false);
    setVisibleRangeInternal(0, frameCount_ - 1);
}

void TimeViewportController::setVisibleRange(int first, int last)
{
    ASSERT_THREAD(this);
    fitsCollected_ = false;
    setFollowPlayhead(false);
    setVisibleRangeInternal(first, last);
}

void TimeViewportController::panFrames(int deltaFrames)
{
    ASSERT_THREAD(this);
    if (deltaFrames == 0)
        return;
    setVisibleRange(visibleStart_ + deltaFrames, visibleEnd_ + deltaFrames);
}

void TimeViewportController::zoomAround(int anchorFrame, double factor)
{
    ASSERT_THREAD(this);
    if (!std::isfinite(factor) || factor <= 0.0)
        return;

    const int oldSpan = visibleEnd_ - visibleStart_ + 1;
    const int newSpan = clampSpan(static_cast<int>(std::round(oldSpan * factor)));
    const int anchor = clampFrame(anchorFrame);
    const double anchorFrac = oldSpan > 1
                                  ? static_cast<double>(anchor - visibleStart_) / static_cast<double>(oldSpan - 1)
                                  : 0.5;
    const int newStart = anchor - static_cast<int>(std::round(anchorFrac * (newSpan - 1)));
    setVisibleRange(newStart, newStart + newSpan - 1);
}

void TimeViewportController::setSelectedRange(int first, int last)
{
    ASSERT_THREAD(this);
    int a = clampFrame(first);
    int b = clampFrame(last);
    if (a > b)
        std::swap(a, b);

    if (hasSelectedRange_ && a == selectedStart_ && b == selectedEnd_)
        return;

    selectedStart_ = a;
    selectedEnd_ = b;
    hasSelectedRange_ = true;
    emit selectedRangeChanged(selectedStart_, selectedEnd_, true);
}

void TimeViewportController::clearSelectedRange()
{
    ASSERT_THREAD(this);
    if (!hasSelectedRange_)
        return;
    hasSelectedRange_ = false;
    emit selectedRangeChanged(selectedStart_, selectedEnd_, false);
}

void TimeViewportController::requestPlaybackFrame(int frame)
{
    ASSERT_THREAD(this);
    emit playbackFrameRequested(clampFrame(frame));
}

void TimeViewportController::updateFollowWindow()
{
    const int span = clampSpan(windowFrames_);
    if (frameCount_ <= span) {
        setVisibleRangeInternal(0, frameCount_ - 1);
        return;
    }

    int first = 0;
    int last = span - 1;
    if (currentFrame_ >= span) {
        last = currentFrame_;
        first = last - span + 1;
    }

    setVisibleRangeInternal(first, last);
}

void TimeViewportController::updateCollectedWindow()
{
    setVisibleRangeInternal(0, currentFrame_);
}

void TimeViewportController::setVisibleRangeInternal(int first, int last)
{
    int a = clampFrame(first);
    int b = clampFrame(last);
    if (a > b)
        std::swap(a, b);

    if (a == visibleStart_ && b == visibleEnd_)
        return;

    visibleStart_ = a;
    visibleEnd_ = b;
    emit visibleRangeChanged(visibleStart_, visibleEnd_);
}

}  // namespace h5reader::app
