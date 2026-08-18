#include "SceneVideoExporter.h"

#include "MoleculeScene.h"
#include "QtPlaybackController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QLoggingCategory>
#include <QMediaCaptureSession>
#include <QMediaFormat>
#include <QUrl>
#include <QVideoFrameFormat>
#include <QVideoFrameInput>

#include <vtkImageData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkWindowToImageFilter.h>

#include <cstring>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cSceneVideo, "h5reader.video")

QMediaFormat SelectMediaFormat(QString* codecName) {
    QMediaFormat format(QMediaFormat::MPEG4);
    format.setVideoCodec(QMediaFormat::VideoCodec::H264);
    if (format.isSupported(QMediaFormat::Encode)) {
        if (codecName) *codecName = QStringLiteral("H.264 / MPEG-4");
        return format;
    }

    format.setVideoCodec(QMediaFormat::VideoCodec::MPEG4);
    if (format.isSupported(QMediaFormat::Encode)) {
        if (codecName) *codecName = QStringLiteral("MPEG-4 Part 2 / MPEG-4");
        return format;
    }

    return {};
}
}  // namespace

QString SceneVideoExportStateName(SceneVideoExportState state) {
    switch (state) {
    case SceneVideoExportState::Idle:      return QStringLiteral("idle");
    case SceneVideoExportState::Starting:  return QStringLiteral("starting");
    case SceneVideoExportState::Recording: return QStringLiteral("recording");
    case SceneVideoExportState::Stopping:  return QStringLiteral("stopping");
    case SceneVideoExportState::Completed: return QStringLiteral("completed");
    case SceneVideoExportState::Stopped:   return QStringLiteral("stopped");
    case SceneVideoExportState::Failed:    return QStringLiteral("failed");
    }
    return QStringLiteral("unknown");
}

SceneVideoExporter::SceneVideoExporter(QObject* parent)
    : QObject(parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("SceneVideoExporter"));
}

SceneVideoExporter::~SceneVideoExporter() {
    ASSERT_THREAD(this);
    if (active_) {
        qCWarning(cSceneVideo).noquote()
            << "video exporter destroyed while active; stopping recorder";
        if (recorder_)
            recorder_->stop();
    }
    destroyMediaObjects();
}

void SceneVideoExporter::setContext(MoleculeScene* scene,
                                    QtPlaybackController* playback) {
    ASSERT_THREAD(this);
    if (active_ && (scene != exportScene_ || playback != exportPlayback_))
        fail(QStringLiteral("the loaded scene changed during video export"));
    scene_ = scene;
    playback_ = playback;
}

bool SceneVideoExporter::start(const SceneVideoExportRequest& request,
                               QString* error) {
    ASSERT_THREAD(this);
    const auto reject = [error](const QString& message) {
        if (error) *error = message;
        return false;
    };

    if (active_)
        return reject(QStringLiteral("a scene video export is already running"));
    if (!scene_ || !playback_)
        return reject(QStringLiteral("no loaded trajectory scene is available"));
    if (request.frameStep < 1)
        return reject(QStringLiteral("frame_step must be at least 1"));
    if (request.framesPerSecond < 1 || request.framesPerSecond > 60)
        return reject(QStringLiteral("frames_per_second must be between 1 and 60"));
    if (request.startFrame < 0 || request.endFrame < request.startFrame
        || request.endFrame >= playback_->frameCount()) {
        return reject(QStringLiteral(
            "start_frame and end_frame must name an ascending range within the loaded trajectory"));
    }
    if (request.outputPath.trimmed().isEmpty())
        return reject(QStringLiteral("output_path must not be empty"));

    QFileInfo outputInfo(QDir::cleanPath(QFileInfo(request.outputPath).absoluteFilePath()));
    if (outputInfo.suffix().compare(QStringLiteral("mp4"), Qt::CaseInsensitive) != 0)
        return reject(QStringLiteral("output_path must end in .mp4"));
    if (outputInfo.exists())
        return reject(QStringLiteral("output_path already exists"));
    if (!QDir().mkpath(outputInfo.absolutePath()))
        return reject(QStringLiteral("could not create the output directory"));

    vtkRenderWindow* renderWindow = scene_->Renderer()
        ? scene_->Renderer()->GetRenderWindow() : nullptr;
    if (!renderWindow)
        return reject(QStringLiteral("the scene has no render window"));
    int* renderSize = renderWindow->GetSize();
    if (!renderSize || renderSize[0] < 2 || renderSize[1] < 2)
        return reject(QStringLiteral("the scene render window has no usable size"));

    const QSize resolution(renderSize[0] & ~1, renderSize[1] & ~1);
    if (resolution.width() < 2 || resolution.height() < 2)
        return reject(QStringLiteral("the scene render size is too small for video"));

    QString codecName;
    const QMediaFormat mediaFormat = SelectMediaFormat(&codecName);
    if (!mediaFormat.isSupported(QMediaFormat::Encode)) {
        return reject(QStringLiteral(
            "Qt Multimedia has no supported H.264 or MPEG-4 video encoder"));
    }

    destroyMediaObjects();

    status_ = {};
    status_.state = SceneVideoExportState::Starting;
    status_.outputPath = outputInfo.absoluteFilePath();
    status_.codec = codecName;
    status_.resolution = resolution;
    status_.startFrame = request.startFrame;
    status_.endFrame = request.endFrame;
    status_.frameStep = request.frameStep;
    status_.framesPerSecond = request.framesPerSecond;
    status_.framesTotal =
        (request.endFrame - request.startFrame) / request.frameStep + 1;

    exportScene_ = scene_;
    exportPlayback_ = playback_;
    originalFrame_ = exportPlayback_->currentFrame();
    originalPlaying_ = exportPlayback_->isPlaying();
    originalDirection_ = exportPlayback_->direction();
    exportPlayback_->pause();

    nextFrame_ = request.startFrame;
    requestedFrame_ = -1;
    pendingFrame_ = {};
    active_ = true;
    inputReady_ = false;
    waitingForRender_ = false;
    ending_ = false;
    endOfStreamSent_ = false;
    stoppedByRequest_ = false;
    restoreDisplayOnFinish_ = true;

    QVideoFrameFormat frameFormat(
        resolution, QVideoFrameFormat::Format_RGBA8888);
    frameFormat.setStreamFrameRate(request.framesPerSecond);

    captureSession_ = new QMediaCaptureSession(this);
    recorder_ = new QMediaRecorder(this);
    videoInput_ = new QVideoFrameInput(frameFormat, this);
    captureSession_->setRecorder(recorder_);
    captureSession_->setVideoFrameInput(videoInput_);

    recorder_->setMediaFormat(mediaFormat);
    recorder_->setEncodingMode(QMediaRecorder::ConstantQualityEncoding);
    recorder_->setQuality(QMediaRecorder::HighQuality);
    recorder_->setVideoResolution(resolution);
    recorder_->setVideoFrameRate(request.framesPerSecond);
    recorder_->setAutoStop(true);
    recorder_->setOutputLocation(QUrl::fromLocalFile(status_.outputPath));

    ACONNECT(videoInput_, &QVideoFrameInput::readyToSendVideoFrame,
             this, &SceneVideoExporter::onInputReady);
    ACONNECT(recorder_, &QMediaRecorder::recorderStateChanged,
             this, &SceneVideoExporter::onRecorderStateChanged);
    ACONNECT(recorder_, &QMediaRecorder::errorOccurred,
             this, &SceneVideoExporter::onRecorderError);
    renderCompletedConnection_ =
        ACONNECT(exportScene_, &MoleculeScene::renderCompleted,
                 this, &SceneVideoExporter::onRenderCompleted);
    sceneDestroyedConnection_ =
        ACONNECT(exportScene_, &QObject::destroyed, this, [this]() {
            fail(QStringLiteral("the scene was closed during video export"));
        });
    playbackDestroyedConnection_ =
        ACONNECT(exportPlayback_, &QObject::destroyed, this, [this]() {
            fail(QStringLiteral("the trajectory was closed during video export"));
        });

    qCInfo(cSceneVideo).noquote()
        << "scene video export started"
        << "| output=" << status_.outputPath
        << "| frames=" << status_.startFrame << ".." << status_.endFrame
        << "| step=" << status_.frameStep
        << "| fps=" << status_.framesPerSecond
        << "| resolution=" << status_.resolution
        << "| codec=" << status_.codec;

    recorder_->record();
    if (recorder_->error() != QMediaRecorder::NoError) {
        const QString message = recorder_->errorString().isEmpty()
            ? QStringLiteral("Qt Multimedia could not start the recorder")
            : recorder_->errorString();
        fail(message);
        if (error) *error = message;
        return false;
    }
    return true;
}

bool SceneVideoExporter::requestStop(bool restoreDisplay) {
    ASSERT_THREAD(this);
    if (!active_)
        return false;

    stoppedByRequest_ = true;
    restoreDisplayOnFinish_ = restoreDisplay;
    ending_ = true;
    status_.state = SceneVideoExportState::Stopping;
    pendingFrame_ = {};
    waitingForRender_ = false;

    qCInfo(cSceneVideo).noquote()
        << "scene video export stop requested"
        << "| frames_written=" << status_.framesWritten;

    if (!recorder_ || recorder_->recorderState() != QMediaRecorder::RecordingState) {
        if (recorder_)
            recorder_->stop();
        if (!recorder_ || recorder_->recorderState() == QMediaRecorder::StoppedState)
            complete();
        return true;
    }

    pump();
    return true;
}

void SceneVideoExporter::stopForShutdown() {
    requestStop(false);
}

void SceneVideoExporter::onInputReady() {
    ASSERT_THREAD(this);
    if (!active_)
        return;
    inputReady_ = true;
    pump();
}

void SceneVideoExporter::onRenderCompleted(int frame) {
    ASSERT_THREAD(this);
    if (!active_ || !waitingForRender_)
        return;

    if (frame != requestedFrame_) {
        qCInfo(cSceneVideo).noquote()
            << "scene video export render superseded; requesting frame again"
            << "| requested=" << requestedFrame_
            << "| rendered=" << frame;
        waitingForRender_ = false;
        requestNextRender();
        return;
    }

    waitingForRender_ = false;
    if (ending_) {
        pump();
        return;
    }

    const QImage image = captureCurrentScene();
    if (image.isNull()) {
        fail(QStringLiteral("the rendered scene could not be read back"));
        return;
    }

    pendingFrame_ = QVideoFrame(image);
    if (!pendingFrame_.isValid()) {
        fail(QStringLiteral("Qt Multimedia rejected the rendered image format"));
        return;
    }

    const qint64 startMicroseconds =
        static_cast<qint64>(status_.framesWritten) * 1'000'000
        / status_.framesPerSecond;
    const qint64 endMicroseconds =
        static_cast<qint64>(status_.framesWritten + 1) * 1'000'000
        / status_.framesPerSecond;
    pendingFrame_.setStartTime(startMicroseconds);
    pendingFrame_.setEndTime(endMicroseconds);
    pendingFrame_.setStreamFrameRate(status_.framesPerSecond);
    pump();
}

void SceneVideoExporter::onRecorderStateChanged(
    QMediaRecorder::RecorderState state) {
    ASSERT_THREAD(this);
    if (!active_)
        return;

    if (state == QMediaRecorder::RecordingState) {
        if (!ending_)
            status_.state = SceneVideoExportState::Recording;
        pump();
        return;
    }
    if (state == QMediaRecorder::StoppedState)
        complete();
}

void SceneVideoExporter::onRecorderError(QMediaRecorder::Error,
                                         const QString& message) {
    ASSERT_THREAD(this);
    fail(message.isEmpty()
             ? QStringLiteral("Qt Multimedia reported a recording error")
             : message);
}

void SceneVideoExporter::pump() {
    ASSERT_THREAD(this);
    if (!active_ || !recorder_ || !videoInput_)
        return;
    if (recorder_->recorderState() != QMediaRecorder::RecordingState)
        return;

    if (ending_) {
        finishInput();
        return;
    }
    if (!inputReady_ || waitingForRender_)
        return;

    if (pendingFrame_.isValid()) {
        if (!videoInput_->sendVideoFrame(pendingFrame_)) {
            inputReady_ = false;
            return;
        }

        pendingFrame_ = {};
        ++status_.framesWritten;
        status_.lastFrame = requestedFrame_;
        requestedFrame_ = -1;
        nextFrame_ += status_.frameStep;
        if (nextFrame_ > status_.endFrame) {
            ending_ = true;
            status_.state = SceneVideoExportState::Stopping;
            finishInput();
            return;
        }
    }

    requestNextRender();
}

void SceneVideoExporter::requestNextRender() {
    ASSERT_THREAD(this);
    if (!active_ || ending_ || waitingForRender_ || !inputReady_)
        return;
    if (!exportScene_ || !exportPlayback_) {
        fail(QStringLiteral("the loaded scene became unavailable"));
        return;
    }

    requestedFrame_ = nextFrame_;
    waitingForRender_ = true;
    if (exportPlayback_->currentFrame() == requestedFrame_)
        exportScene_->refreshCurrentFrame();
    else
        exportPlayback_->setFrame(requestedFrame_);
}

void SceneVideoExporter::finishInput() {
    ASSERT_THREAD(this);
    if (!active_ || endOfStreamSent_ || !inputReady_ || !videoInput_)
        return;

    if (!videoInput_->sendVideoFrame(QVideoFrame{})) {
        inputReady_ = false;
        return;
    }
    endOfStreamSent_ = true;
}

void SceneVideoExporter::fail(const QString& message) {
    ASSERT_THREAD(this);
    if (!active_)
        return;
    if (status_.error.isEmpty())
        status_.error = message;
    status_.state = SceneVideoExportState::Stopping;
    ending_ = true;
    pendingFrame_ = {};
    waitingForRender_ = false;

    qCWarning(cSceneVideo).noquote()
        << "scene video export failed"
        << "| output=" << status_.outputPath
        << "| error=" << status_.error;

    if (recorder_ && recorder_->recorderState() != QMediaRecorder::StoppedState)
        recorder_->stop();
    else
        complete();
}

void SceneVideoExporter::complete() {
    ASSERT_THREAD(this);
    if (!active_)
        return;

    active_ = false;
    QObject::disconnect(renderCompletedConnection_);
    QObject::disconnect(sceneDestroyedConnection_);
    QObject::disconnect(playbackDestroyedConnection_);

    if (recorder_ && recorder_->actualLocation().isLocalFile())
        status_.outputPath = recorder_->actualLocation().toLocalFile();
    if (stoppedByRequest_ && status_.framesWritten == 0
        && status_.error.isEmpty()) {
        status_.error = QStringLiteral(
            "video export stopped before the first frame was written");
    }
    if (!status_.error.isEmpty()) {
        status_.state = SceneVideoExportState::Failed;
        if (QFileInfo::exists(status_.outputPath)
            && !QFile::remove(status_.outputPath)) {
            qCWarning(cSceneVideo).noquote()
                << "could not remove failed scene video output"
                << "| output=" << status_.outputPath;
        }
    } else if (stoppedByRequest_) {
        status_.state = SceneVideoExportState::Stopped;
    } else {
        status_.state = SceneVideoExportState::Completed;
    }

    const QFileInfo outputInfo(status_.outputPath);
    status_.fileSizeBytes = outputInfo.exists() ? outputInfo.size() : 0;

    if (restoreDisplayOnFinish_)
        restoreDisplay();

    qCInfo(cSceneVideo).noquote()
        << "scene video export finished"
        << "| state=" << SceneVideoExportStateName(status_.state)
        << "| frames=" << status_.framesWritten << "/" << status_.framesTotal
        << "| bytes=" << status_.fileSizeBytes
        << "| output=" << status_.outputPath;

    exportScene_.clear();
    exportPlayback_.clear();
    emit finished();
}

void SceneVideoExporter::restoreDisplay() {
    ASSERT_THREAD(this);
    if (!exportScene_ || !exportPlayback_)
        return;

    if (exportPlayback_->currentFrame() == originalFrame_)
        exportScene_->refreshCurrentFrame();
    else
        exportPlayback_->setFrame(originalFrame_);

    if (originalPlaying_) {
        if (originalDirection_ < 0)
            exportPlayback_->playBackward();
        else
            exportPlayback_->playForward();
    }
}

void SceneVideoExporter::destroyMediaObjects() {
    ASSERT_THREAD(this);
    QObject::disconnect(renderCompletedConnection_);
    QObject::disconnect(sceneDestroyedConnection_);
    QObject::disconnect(playbackDestroyedConnection_);

    if (captureSession_) {
        captureSession_->setVideoFrameInput(nullptr);
        captureSession_->setRecorder(nullptr);
    }
    delete videoInput_;
    videoInput_ = nullptr;
    delete recorder_;
    recorder_ = nullptr;
    delete captureSession_;
    captureSession_ = nullptr;
}

QImage SceneVideoExporter::captureCurrentScene() const {
    ASSERT_THREAD(this);
    if (!exportScene_ || !exportScene_->Renderer()
        || !exportScene_->Renderer()->GetRenderWindow()) {
        return {};
    }

    auto capture = vtkSmartPointer<vtkWindowToImageFilter>::New();
    capture->SetInput(exportScene_->Renderer()->GetRenderWindow());
    capture->SetInputBufferTypeToRGBA();
    capture->ReadFrontBufferOff();
    capture->ShouldRerenderOff();
    capture->Update();

    vtkImageData* pixels = capture->GetOutput();
    if (!pixels || pixels->GetNumberOfScalarComponents() != 4)
        return {};

    int dimensions[3] = {};
    pixels->GetDimensions(dimensions);
    if (dimensions[0] < 1 || dimensions[1] < 1)
        return {};

    const auto* source = static_cast<const unsigned char*>(
        pixels->GetScalarPointer());
    if (!source)
        return {};

    QImage image(dimensions[0], dimensions[1], QImage::Format_RGBA8888);
    if (image.isNull())
        return {};

    const qsizetype sourceStride = static_cast<qsizetype>(dimensions[0]) * 4;
    for (int y = 0; y < dimensions[1]; ++y) {
        const unsigned char* sourceRow =
            source + static_cast<qsizetype>(dimensions[1] - 1 - y) * sourceStride;
        std::memcpy(image.scanLine(y), sourceRow,
                    static_cast<std::size_t>(sourceStride));
    }

    if (image.size() == status_.resolution)
        return image;
    if (image.width() >= status_.resolution.width()
        && image.height() >= status_.resolution.height()) {
        return image.copy(QRect(QPoint(0, 0), status_.resolution));
    }
    return image.scaled(status_.resolution, Qt::IgnoreAspectRatio,
                        Qt::SmoothTransformation);
}

}  // namespace h5reader::app
