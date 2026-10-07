#include <sciencefiles/Download.h>
#include "ResumeFile.h"

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QFileInfo>
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

    void start(const QUrl &url, const QString &destination, bool resume);
    void cancel();
    void shutdown();

signals:
    void progress(qint64 received, qint64 total);
    void finished(sciencefiles::Download::Result result, const QString &detail);

private:
    bool prepareResponse();
    void fail(Result result, const QString &detail);
    void readAvailable();
    void responseFinished();
    void complete(Result result, const QString &detail);

    QNetworkAccessManager network_;
    QSslConfiguration tls_;
    std::chrono::milliseconds inactivity_;
    QNetworkReply *reply_ = nullptr;
    std::unique_ptr<QSaveFile> output_;
    std::unique_ptr<ResumeFile> resume_;
    QUrl url_;
    QString destination_;
    QString failure_;
    Result failureResult_ = Result::FileError;
    qint64 received_ = 0;
    bool cancelled_ = false;
    bool readScheduled_ = false;
    bool headersAccepted_ = false;
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

void DownloadWorker::start(const QUrl &url, const QString &destination, bool resume) {
    Q_ASSERT(!reply_);
    url_ = url;
    destination_ = destination;
    received_ = 0;
    cancelled_ = false;
    readScheduled_ = false;
    headersAccepted_ = false;
    failure_.clear();
    elapsed_.start();
    const bool localFile = url.isLocalFile();
    const bool https = url.scheme() == "https" && !url.host().isEmpty();
    if (!url.isValid() || (!https && !localFile) || (localFile && url.toLocalFile().isEmpty()) ||
        !url.userInfo().isEmpty() || destination.isEmpty() || inactivity_.count() <= 0 || (resume && !https)) {
        complete(Result::InvalidRequest, tr("An HTTPS or file URL and a destination file are required."));
        return;
    }
    if (localFile) {
        const QString sourcePath = QFileInfo(url.toLocalFile()).canonicalFilePath();
        const QString destinationPath = QFileInfo(destination).canonicalFilePath();
#ifdef Q_OS_WIN
        constexpr auto pathCase = Qt::CaseInsensitive;
#else
        constexpr auto pathCase = Qt::CaseSensitive;
#endif
        if (!sourcePath.isEmpty() && sourcePath.compare(destinationPath, pathCase) == 0) {
            complete(Result::InvalidRequest, tr("Source and destination must be different files."));
            return;
        }
    }
    if (resume) {
        resume_ = std::make_unique<ResumeFile>(url, destination);
        if (!resume_->prepare()) {
            complete(resume_->errorResult(), resume_->errorString());
            return;
        }
        received_ = resume_->offset();
    } else {
        output_ = std::make_unique<QSaveFile>(destination);
        // Never fall back to overwriting an existing file while the download is incomplete.
        output_->setDirectWriteFallback(false);
        if (!output_->open(QIODevice::WriteOnly)) {
            complete(Result::FileError, output_->errorString());
            return;
        }
    }
    QNetworkRequest request(url);
    request.setSslConfiguration(tls_);
    request.setTransferTimeout(inactivity_);
    request.setRawHeader("Accept-Encoding", "identity");
    request.setAttribute(QNetworkRequest::RedirectPolicyAttribute, QNetworkRequest::ManualRedirectPolicy);
    if (resume_)
        resume_->setRequestHeaders(request);
    reply_ = network_.get(request);
    reply_->setReadBufferSize(ReplyBufferBytes);
    connect(reply_, &QNetworkReply::readyRead, this, &DownloadWorker::readAvailable);
    if (!localFile && !resume_)
        connect(reply_, &QNetworkReply::downloadProgress, this, &DownloadWorker::progress);
    // abort() may emit finished synchronously. Finish only after the current read/header handler returns.
    connect(reply_, &QNetworkReply::finished, this, [this, reply = reply_] {
        if (reply_ == reply)
            responseFinished();
    }, Qt::QueuedConnection);
    connect(reply_, &QNetworkReply::metaDataChanged, this, [this] { prepareResponse(); });
    connect(reply_, &QNetworkReply::sslErrors, this, [](const QList<QSslError> &errors) {
        for (const auto &error : errors)
            qCWarning(downloadLog) << "TLS:" << error.errorString();
    });
    qCInfo(downloadLog) << "Started:" << url_ << "destination=" << destination_ << "offset=" << received_;
}

void DownloadWorker::fail(Result result, const QString &detail) {
    if (!failure_.isEmpty())
        return;
    failureResult_ = result;
    failure_ = detail;
    if (reply_ && !reply_->isFinished())
        reply_->abort();
}

bool DownloadWorker::prepareResponse() {
    if (!reply_ || cancelled_ || !failure_.isEmpty())
        return false;
    if (headersAccepted_)
        return true;
    if (!url_.isLocalFile()) {
        const int status = reply_->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
        if (status == 0)
            return false;
        if (resume_) {
            if (!resume_->openResponse(*reply_)) {
                fail(resume_->errorResult(), resume_->errorString());
                return false;
            }
            received_ = resume_->offset();
            emit progress(received_, resume_->length());
        } else if (status != 200) {
            const QByteArray code = reply_->rawHeader("X-Error-Code");
            fail(Result::HttpError, tr("HTTP %1%2").arg(status).arg(
                code.isEmpty() ? QString() : tr(" (server error %1)").arg(QString::fromLatin1(code))));
            return false;
        }
    }
    headersAccepted_ = true;
    return true;
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
    if (!reply_ || cancelled_ || !failure_.isEmpty() || reply_->error() != QNetworkReply::NoError || !prepareResponse())
        return;
    std::array<char, CopyBufferBytes> buffer;
    qint64 remaining = ReplyBufferBytes;
    const qint64 previous = received_;
    while (reply_->bytesAvailable() > 0 && remaining > 0) {
        const qint64 count = reply_->read(buffer.data(), qMin(remaining, static_cast<qint64>(buffer.size())));
        if (count < 0) {
            fail(Result::FileError, tr("Could not read the response: %1").arg(reply_->errorString()));
            return;
        } else if (count == 0) {
            break;
        }
        if (resume_) {
            if (!resume_->write(buffer.data(), count)) {
                fail(resume_->errorResult(), resume_->errorString());
                return;
            }
        } else if (output_->write(buffer.data(), count) != count) {
            fail(Result::FileError, tr("Could not write the download: %1").arg(output_->errorString()));
            return;
        }
        received_ += count;
        remaining -= count;
    }
    if (resume_ && received_ != previous) {
        emit progress(received_, resume_->length());
    } else if (url_.isLocalFile() && received_ != previous) {
        bool lengthValid = false;
        const qint64 total = reply_->header(QNetworkRequest::ContentLengthHeader).toLongLong(&lengthValid);
        emit progress(received_, lengthValid ? total : -1);
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
    if (!reply_)
        return;
    if (cancelled_) {
        complete(Result::Cancelled, resume_ ? tr("Cancelled.") : tr("Cancelled. No incomplete file was kept."));
        return;
    }
    prepareResponse();
    if (!failure_.isEmpty()) {
        complete(failureResult_, failure_);
        return;
    }
    const int status = reply_->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
    if (reply_->error() != QNetworkReply::NoError) {
        complete(Result::NetworkError, reply_->errorString());
        return;
    }
    if (!headersAccepted_) {
        complete(Result::HttpError, tr("Missing HTTP response headers."));
        return;
    }
    // The finished signal can precede consumption of the last buffered bytes.
    readAvailable();
    if (!failure_.isEmpty()) {
        complete(failureResult_, failure_);
        return;
    }
    if (reply_->bytesAvailable() > 0)
        return;
    if (resume_) {
        if (!resume_->commit()) {
            complete(resume_->errorResult(), resume_->errorString());
            return;
        }
        complete(Result::Saved, tr("Saved %1 bytes to %2").arg(received_).arg(destination_));
        return;
    }
    bool lengthValid = false;
    const qint64 expected = reply_->header(QNetworkRequest::ContentLengthHeader).toLongLong(&lengthValid);
    const QByteArray encoding = reply_->rawHeader("Content-Encoding");
    if ((!url_.isLocalFile() && status != 200) || !lengthValid || expected != received_ ||
        (!encoding.isEmpty() && encoding != "identity")) {
        complete(url_.isLocalFile() ? Result::FileError : Result::HttpError,
                 tr("Incomplete or unexpected response (%1 bytes received).").arg(received_));
        return;
    }
    if (!output_->commit()) {
        complete(Result::FileError, output_->errorString());
        return;
    }
    if (url_.isLocalFile() && received_ == 0)
        emit progress(0, 0);
    complete(Result::Saved, tr("Saved %1 bytes to %2").arg(received_).arg(destination_));
}

void DownloadWorker::complete(Result result, const QString &detail) {
    // Copy before releasing output_: detail can refer to one of our error strings.
    QString message = detail;
    output_.reset();
    if (resume_) {
        if (!resume_->finish()) {
            result = Result::FileError;
            message = resume_->errorString();
        }
        resume_.reset();
    }
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

bool Download::start(const QUrl &url, const QString &destination, bool resume) {
    Q_ASSERT(thread() == QThread::currentThread());
    if (active_ || shutdownRequested_) {
        qCWarning(downloadLog) << "Download active or shutting down; request rejected:" << url;
        return false;
    }
    active_ = true;
    QMetaObject::invokeMethod(
        worker_, [worker = worker_, url, destination, resume] { worker->start(url, destination, resume); }, Qt::QueuedConnection);
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
