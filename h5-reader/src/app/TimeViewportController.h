// TimeViewportController -- shared trajectory-time viewport state.
//
// This is application state, not diagnostics. ErrorBus remains the sink for
// degraded state and loader failures; playback, brush selection, and chart
// ranges move through typed Qt signals from this controller.

#pragma once

#include <QObject>

namespace h5reader::app {

class TimeViewportController final : public QObject {
    Q_OBJECT

public:
    explicit TimeViewportController(int frameCount, QObject* parent = nullptr);
    ~TimeViewportController() override = default;

    int frameCount() const { return frameCount_; }
    int currentFrame() const { return currentFrame_; }
    int visibleStart() const { return visibleStart_; }
    int visibleEnd() const { return visibleEnd_; }
    int windowFrames() const { return windowFrames_; }
    bool followsPlayhead() const { return followsPlayhead_; }

    bool hasSelectedRange() const { return hasSelectedRange_; }
    int selectedStart() const { return selectedStart_; }
    int selectedEnd() const { return selectedEnd_; }

public slots:
    void setCurrentFrame(int frame);
    void setWindowFrames(int frames);
    void setFollowPlayhead(bool follow);
    void followPlayhead();
    void fitCollectedRange();
    void fitFullRange();
    void setVisibleRange(int first, int last);
    void panFrames(int deltaFrames);
    void zoomAround(int anchorFrame, double factor);
    void setSelectedRange(int first, int last);
    void clearSelectedRange();
    void requestPlaybackFrame(int frame);

signals:
    void currentFrameChanged(int frame);
    void visibleRangeChanged(int first, int last);
    void windowFramesChanged(int frames);
    void followPlayheadChanged(bool follow);
    void selectedRangeChanged(int first, int last, bool active);
    void playbackFrameRequested(int frame);

private:
    int clampFrame(int frame) const;
    int clampSpan(int frames) const;
    void updateFollowWindow();
    void updateCollectedWindow();
    void setVisibleRangeInternal(int first, int last);

    int frameCount_ = 1;
    int currentFrame_ = 0;
    int visibleStart_ = 0;
    int visibleEnd_ = 0;
    int windowFrames_ = 120;
    bool followsPlayhead_ = true;
    bool fitsCollected_ = false;

    bool hasSelectedRange_ = false;
    int selectedStart_ = 0;
    int selectedEnd_ = 0;
};

}  // namespace h5reader::app
