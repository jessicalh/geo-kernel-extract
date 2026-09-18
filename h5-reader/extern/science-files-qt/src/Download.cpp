#include <sciencefiles/Download.h>

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QLoggingCategory>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QSaveFile>
#include <QSslError>
#include <array>
#include <memory>

namespace sciencefiles {

Q_LOGGING_CATEGORY(downloadLog, "science.download")

namespace {
constexpr qint64 ReplyBufferBytes = 256 * 1024;
constexpr size_t CopyBufferBytes = 64 * 1024;
} // namespace

class DownloadWorker final : public QObject {
    Q_OBJECT
public:
    using Result = Download::Result;

    DownloadWorker(const QSslConfiguration &tls, std::chrono::milliseconds inactivity)
        : network_(this), tls_(tls), inactivity_(inactivity) {}
    ~DownloadWorker() override;

    void start(const QUrl &url, const QString &destination);
    void cancel();
    void shutdown();

signals:
    void progress(qint64 received, qint64 total);
    void finished(sciencefiles::Download::Result result, const QString &detail);

private:
    void readAvailable();
    void responseFinished();
    void complete(Result result, const QString &detail);

    QNetworkAccessManager network_;
    QSslConfiguration tls_;
    std::chrono::milliseconds inactivity_;
    QNetworkReply *reply_ = nullptr;
    std::unique_ptr<QSaveFile> output_;
    QUrl url_;
    QString destination_;
    QString writeError_;
    qint64 received_ = 0;
    bool cancelled_ = false;
    bool readScheduled_ = false;
    bool stopping_ = false;
    QElapsedTimer elapsed_;
};

DownloadWorker::~DownloadWorker() {
    if (reply_) {
        disconnect(reply_, nullptr, this, nullptr);
        reply_->abort();
        qCInfo(downloadLog) << "Download discarded during destruction:" << destination_;
    }
}

void DownloadWorker::start(const QUrl &url, const QString &destination) {
    Q_ASSERT(!reply_);
    url_ = url;
    destination_ = destination;
    received_ = 0;
    cancelled_ = false;
    readScheduled_ = false;
    writeError_.clear();
    elapsed_.start();
    if (!url.isValid() || url.scheme() != "https" || url.host().isEmpty() || !url.userInfo().isEmpty() ||
        destination.isEmpty() || inactivity_.count() <= 0) {
        complete(Result::InvalidRequest, tr("An HTTPS URL and a destination file are required."));
        return;
    }
    output_ = std::make_unique<QSaveFile>(destination);
    // Never fall back to overwriting an existing file while the download is incomplete.
    output_->setDirectWriteFallback(false);
    if (!output_->open(QIODevice::WriteOnly)) {
        complete(Result::FileError, output_->errorString());
        return;
    }
    QNetworkRequest request(url);
    request.setSslConfiguration(tls_);
    request.setTransferTimeout(inactivity_);
    request.setRawHeader("Accept-Encoding", "identity");
    request.setAttribute(QNetworkRequest::RedirectPolicyAttribute, QNetworkRequest::ManualRedirectPolicy);
    reply_ = network_.get(request);
    reply_->setReadBufferSize(ReplyBufferBytes);
    connect(reply_, &QNetworkReply::readyRead, this, &DownloadWorker::readAvailable);
    connect(reply_, &QNetworkReply::downloadProgress, this, &DownloadWorker::progress);
    connect(reply_, &QNetworkReply::finished, this, &DownloadWorker::responseFinished);
    connect(reply_, &QNetworkReply::metaDataChanged, this, [this] {
        const int status = reply_->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
        // Reject from the headers rather than leave an unwanted body filling the read buffer.
        if (status > 200)
            reply_->abort();
    });
    connect(reply_, &QNetworkReply::sslErrors, this, [](const QList<QSslError> &errors) {
        for (const auto &error : errors)
            qCWarning(downloadLog) << "TLS:" << error.errorString();
    });
    qCInfo(downloadLog) << "Started:" << url_ << "destination=" << destination_;
}

void DownloadWorker::cancel() {
    if (!reply_)
        return;
    cancelled_ = true;
    if (reply_->isFinished())
        responseFinished();
    else
        reply_->abort();
}

void DownloadWorker::shutdown() {
    stopping_ = true;
    if (reply_)
        cancel();
    else
        deleteLater();
}

void DownloadWorker::readAvailable() {
    if (!reply_ || cancelled_ || !writeError_.isEmpty() || reply_->error() != QNetworkReply::NoError)
        return;
    if (reply_->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt() != 200)
        return;
    std::array<char, CopyBufferBytes> buffer;
    qint64 remaining = ReplyBufferBytes;
    while (reply_->bytesAvailable() > 0 && remaining > 0) {
        const qint64 count = reply_->read(buffer.data(), qMin(remaining, static_cast<qint64>(buffer.size())));
        if (count < 0) {
            writeError_ = tr("Could not read the response: %1").arg(reply_->errorString());
        } else if (count == 0) {
            break;
        } else if (output_->write(buffer.data(), count) != count) {
            writeError_ = tr("Could not write the download: %1").arg(output_->errorString());
        } else {
            received_ += count;
            remaining -= count;
        }
        if (!writeError_.isEmpty()) {
            if (!reply_->isFinished())
                reply_->abort();
            return;
        }
    }
    if (reply_->bytesAvailable() > 0 && !readScheduled_) {
        // A reply can expose more than its requested buffer size. Yield between disk-copy batches.
        readScheduled_ = true;
        QMetaObject::invokeMethod(
            reply_,
            [this, reply = reply_] {
                if (reply_ != reply)
                    return;
                readScheduled_ = false;
                if (reply->isFinished())
                    responseFinished();
                else
                    readAvailable();
            },
            Qt::QueuedConnection);
    }
}

void DownloadWorker::responseFinished() {
    if (cancelled_) {
        complete(Result::Cancelled, tr("Cancelled. No incomplete file was kept."));
        return;
    }
    if (!writeError_.isEmpty()) {
        complete(Result::FileError, writeError_);
        return;
    }
    const int status = reply_->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
    if (status != 0 && status != 200) {
        const QByteArray code = reply_->rawHeader("X-Error-Code");
        complete(Result::HttpError,
                 tr("HTTP %1%2")
                     .arg(status)
                     .arg(code.isEmpty() ? QString() : tr(" (server error %1)").arg(QString::fromLatin1(code))));
        return;
    }
    if (reply_->error() != QNetworkReply::NoError) {
        complete(Result::NetworkError, reply_->errorString());
        return;
    }
    // The finished signal can precede consumption of the last buffered bytes.
    readAvailable();
    if (!writeError_.isEmpty()) {
        complete(Result::FileError, writeError_);
        return;
    }
    if (reply_->bytesAvailable() > 0)
        return;
    bool lengthValid = false;
    const qint64 expected = reply_->header(QNetworkRequest::ContentLengthHeader).toLongLong(&lengthValid);
    const QByteArray encoding = reply_->rawHeader("Content-Encoding");
    if (status != 200 || !lengthValid || expected != received_ || (!encoding.isEmpty() && encoding != "identity")) {
        complete(Result::HttpError, tr("Incomplete or unexpected response (%1 bytes received).").arg(received_));
        return;
    }
    if (!output_->commit()) {
        complete(Result::FileError, output_->errorString());
        return;
    }
    complete(Result::Saved, tr("Saved %1 bytes to %2").arg(received_).arg(destination_));
}

void DownloadWorker::complete(Result result, const QString &detail) {
    // Copy before releasing output_: detail can refer to one of our error strings.
    const QString message = detail;
    output_.reset();
    if (reply_) {
        disconnect(reply_, nullptr, this, nullptr);
        reply_->deleteLater();
        reply_ = nullptr;
    }
    qCInfo(downloadLog) << "Finished: code=" << static_cast<int>(result) << "url=" << url_ << "bytes=" << received_
                        << "elapsed_ms=" << elapsed_.elapsed() << message;
    emit finished(result, message);
    if (stopping_)
        deleteLater();
}

Download::Download(QObject *parent, const QSslConfiguration &tls, std::chrono::milliseconds inactivity)
    : QObject(parent), workerThread_(this), worker_(new DownloadWorker(tls, inactivity)) {
    worker_->moveToThread(&workerThread_);
    connect(&workerThread_, &QThread::finished, worker_, &QObject::deleteLater);
    // Destroy network resources on their own thread before ending its event loop.
    connect(worker_, &QObject::destroyed, &workerThread_, &QThread::quit, Qt::DirectConnection);
    connect(&workerThread_, &QThread::finished, this, [this] {
        shutdownComplete_ = true;
        qCInfo(downloadLog) << "Download worker stopped.";
        emit shutdownFinished();
    });
    connect(worker_, &DownloadWorker::progress, this, &Download::progress, Qt::QueuedConnection);
    connect(
        worker_, &DownloadWorker::finished, this,
        [this](Result result, const QString &detail) {
            active_ = false;
            emit finished(result, detail);
        },
        Qt::QueuedConnection);
    connect(qApp, &QCoreApplication::aboutToQuit, this, &Download::shutdown);
    workerThread_.start();
}

Download::~Download() {
    Q_ASSERT(thread() == QThread::currentThread());
    // Fallback for destruction without an asynchronous shutdown; no GUI callback is needed.
    shutdown();
    workerThread_.wait();
}

bool Download::start(const QUrl &url, const QString &destination) {
    Q_ASSERT(thread() == QThread::currentThread());
    if (active_ || shutdownRequested_) {
        qCWarning(downloadLog) << "Download active or shutting down; request rejected:" << url;
        return false;
    }
    active_ = true;
    QMetaObject::invokeMethod(
        worker_, [worker = worker_, url, destination] { worker->start(url, destination); }, Qt::QueuedConnection);
    return true;
}

void Download::cancel() {
    Q_ASSERT(thread() == QThread::currentThread());
    if (active_ && !shutdownRequested_)
        QMetaObject::invokeMethod(worker_, &DownloadWorker::cancel, Qt::QueuedConnection);
}

void Download::shutdown() {
    Q_ASSERT(thread() == QThread::currentThread());
    if (shutdownRequested_)
        return;
    shutdownRequested_ = true;
    qCInfo(downloadLog) << "Download shutdown requested.";
    QMetaObject::invokeMethod(worker_, &DownloadWorker::shutdown, Qt::QueuedConnection);
}

} // namespace sciencefiles

#include "Download.moc"
