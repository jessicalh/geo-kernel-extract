// QtPlaybackController — frame-scrub + play/pause for trajectory playback.
//
// The ONE class whose QTimer is legitimate (see feedback_qt_discipline).
// Timing is its responsibility — everyone else connects to frameChanged
// and reacts. No interactive control anywhere else in the reader owns a
// timer.
//
// Signals:
//   frameChanged(int t)     — emitted on every frame advance, scrub, or
//                              explicit setFrame call. Sole driver of
//                              MoleculeScene, overlays, atom inspector.
//   playingChanged(bool)    — emitted when play/pause state changes.
//   fpsChanged(int)         — emitted when playback speed changes.
//
// Thread affinity: GUI thread only. All public slots and signals
// marshal through the standard Qt event loop.

#pragma once

#include <QObject>
#include <QTimer>

namespace h5reader::app {

class QtPlaybackController final : public QObject {
    Q_OBJECT

public:
    // frameCount is the trajectory's QtConformation::frameCount(). fps
    // defaults to 10 frames/sec — user can override via setFps.
    explicit QtPlaybackController(int frameCount, QObject* parent = nullptr);
    ~QtPlaybackController() override = default;

    // ----- Current state -----
    int  currentFrame() const { return currentFrame_; }
    int  frameCount()   const { return frameCount_; }
    bool isPlaying()    const { return timer_.isActive(); }
    int  fps()          const { return fps_; }
    bool isLooping()    const { return looping_; }

public slots:
    void play();
    void pause();
    void togglePlayPause();
    void setFrame(int t);       // clamps to [0, frameCount-1]
    void stepForward();         // one frame, irrespective of play state
    void stepBackward();
    void setFps(int fps);       // clamped to [1, 60]; adjusts timer interval
    void setLooping(bool loop);

signals:
    void frameChanged(int t);
    void playingChanged(bool playing);
    void fpsChanged(int fps);

private slots:
    void advance();   // QTimer::timeout slot

private:
    QTimer timer_;
    int    frameCount_;
    int    currentFrame_ = 0;
    // Default 5 fps: at typical extraction stride of ~40 ps/frame, 5 fps
    // gives 200 ms of screen time per 40 ps of simulation — enough to
    // read each configuration without missing the dynamics. Users who
    // want a fast survey spin up the fpsSpinner; 10+ fps works but
    // individual configurations blur into each other.
    int    fps_          = 5;
    bool   looping_      = true;
};

}  // namespace h5reader::app
