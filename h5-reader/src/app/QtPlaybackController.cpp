#include "QtPlaybackController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QLoggingCategory>

#include <algorithm>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cPlayback, "h5reader.playback")
}

QtPlaybackController::QtPlaybackController(int frameCount, QObject* parent)
    : QObject(parent),
      frameCount_(std::max(1, frameCount))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtPlaybackController"));

    timer_.setInterval(1000 / fps_);
    timer_.setTimerType(Qt::PreciseTimer);

    ACONNECT(&timer_, &QTimer::timeout, this, &QtPlaybackController::advance);

    qCInfo(cPlayback).noquote()
        << "created | frames=" << frameCount_ << "| fps=" << fps_;
}

void QtPlaybackController::play() {
    ASSERT_THREAD(this);
    if (timer_.isActive()) return;
    timer_.start();
    qCInfo(cPlayback).noquote() << "play";
    emit playingChanged(true);
}

void QtPlaybackController::pause() {
    ASSERT_THREAD(this);
    if (!timer_.isActive()) return;
    timer_.stop();
    qCInfo(cPlayback).noquote() << "pause @ frame" << currentFrame_;
    emit playingChanged(false);
}

void QtPlaybackController::togglePlayPause() {
    if (isPlaying()) pause();
    else             play();
}

void QtPlaybackController::setFrame(int t) {
    ASSERT_THREAD(this);
    const int clamped = std::clamp(t, 0, frameCount_ - 1);
    if (clamped == currentFrame_) return;
    currentFrame_ = clamped;
    emit frameChanged(currentFrame_);
}

void QtPlaybackController::stepForward() {
    setFrame(currentFrame_ + 1);
}

void QtPlaybackController::stepBackward() {
    setFrame(currentFrame_ - 1);
}

void QtPlaybackController::setFps(int fps) {
    ASSERT_THREAD(this);
    const int clamped = std::clamp(fps, 1, 60);
    if (clamped == fps_) return;
    fps_ = clamped;
    timer_.setInterval(1000 / fps_);
    qCInfo(cPlayback).noquote() << "fps →" << fps_;
    emit fpsChanged(fps_);
}

void QtPlaybackController::setLooping(bool loop) {
    ASSERT_THREAD(this);
    looping_ = loop;
}

void QtPlaybackController::advance() {
    ASSERT_THREAD(this);
    int next = currentFrame_ + 1;
    if (next >= frameCount_) {
        if (looping_) {
            next = 0;
        } else {
            pause();
            return;
        }
    }
    currentFrame_ = next;
    emit frameChanged(currentFrame_);
}

}  // namespace h5reader::app
