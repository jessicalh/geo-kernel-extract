#pragma once

#include <QImage>
#include <QMediaRecorder>
#include <QObject>
#include <QPointer>
#include <QSize>
#include <QString>
#include <QVideoFrame>

class QMediaCaptureSession;
class QVideoFrameInput;

namespace h5reader::app {

class MoleculeScene;
class QtPlaybackController;

struct SceneVideoExportRequest {
    QString outputPath;
    int startFrame = 0;
    int endFrame = 0;
    int frameStep = 1;
    int framesPerSecond = 10;
};

enum class SceneVideoExportState {
    Idle,
    Starting,
    Recording,
    Stopping,
    Completed,
    Stopped,
    Failed,
};

struct SceneVideoExportStatus {
    SceneVideoExportState state = SceneVideoExportState::Idle;
    QString outputPath;
    QString codec;
    QString error;
    QSize resolution;
    int startFrame = 0;
    int endFrame = 0;
    int frameStep = 1;
    int framesPerSecond = 10;
    int framesTotal = 0;
    int framesWritten = 0;
    int lastFrame = -1;
    qint64 fileSizeBytes = 0;
};

QString SceneVideoExportStateName(SceneVideoExportState state);

// Records the exact current Reader scene over a trajectory frame range.
// Rendering, frame changes, framebuffer reads, and Qt Multimedia objects all
// stay on the GUI thread. Back-pressure is driven by QVideoFrameInput's real
// readiness signal; each capture is driven by MoleculeScene's VTK EndEvent.
class SceneVideoExporter final : public QObject {
    Q_OBJECT

public:
    explicit SceneVideoExporter(QObject* parent = nullptr);
    ~SceneVideoExporter() override;

    void setContext(MoleculeScene* scene, QtPlaybackController* playback);

    bool start(const SceneVideoExportRequest& request, QString* error = nullptr);
    bool requestStop(bool restoreDisplay = true);
    void stopForShutdown();

    bool isActive() const { return active_; }
    const SceneVideoExportStatus& status() const { return status_; }

signals:
    void finished();

private:
    void onInputReady();
    void onRenderCompleted(int frame);
    void onRecorderStateChanged(QMediaRecorder::RecorderState state);
    void onRecorderError(QMediaRecorder::Error error, const QString& message);

    void pump();
    void requestNextRender();
    void finishInput();
    void fail(const QString& message);
    void complete();
    void restoreDisplay();
    void destroyMediaObjects();
    QImage captureCurrentScene() const;

    QPointer<MoleculeScene> scene_;
    QPointer<QtPlaybackController> playback_;
    QPointer<MoleculeScene> exportScene_;
    QPointer<QtPlaybackController> exportPlayback_;

    QMediaCaptureSession* captureSession_ = nullptr;
    QMediaRecorder* recorder_ = nullptr;
    QVideoFrameInput* videoInput_ = nullptr;
    QMetaObject::Connection renderCompletedConnection_;
    QMetaObject::Connection sceneDestroyedConnection_;
    QMetaObject::Connection playbackDestroyedConnection_;

    SceneVideoExportStatus status_;
    QVideoFrame pendingFrame_;
    int nextFrame_ = 0;
    int requestedFrame_ = -1;
    int originalFrame_ = 0;
    int originalDirection_ = 1;
    bool originalPlaying_ = false;
    bool active_ = false;
    bool inputReady_ = false;
    bool waitingForRender_ = false;
    bool ending_ = false;
    bool endOfStreamSent_ = false;
    bool stoppedByRequest_ = false;
    bool restoreDisplayOnFinish_ = true;
};

}  // namespace h5reader::app
