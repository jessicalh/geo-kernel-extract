#include <sciencefiles/BundleCache.h>
#include "ArchiveTask.h"
#include "BundlePaths.h"

#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QLoggingCategory>
#include <QPointer>

namespace sciencefiles {
Q_LOGGING_CATEGORY(cacheLog, "science.cache")

BundleCache::BundleCache(const QString &root, QObject *parent, const QSslConfiguration &tls)
    : QObject(parent), root_(QDir(root).absolutePath()), download_(std::make_unique<Download>(this, tls)) {
    connect(download_.get(), &Download::progress, this, [this](qint64 bytes, qint64 total) {
        if (state_ == State::Downloading)
            emit progress(bytes, total);
    });
    connect(download_.get(), &Download::finished, this, &BundleCache::downloaded);
    connect(download_.get(), &Download::shutdownFinished, this, [this] {
        shutdownComplete_ = true;
        qCInfo(cacheLog) << "Cache shutdown complete.";
        emit shutdownFinished();
    });
    connect(qApp, &QCoreApplication::aboutToQuit, this, &BundleCache::shutdown);
}

BundleCache::~BundleCache() {
    disconnect(download_.get(), nullptr, this, nullptr);
    // Join the transfer before removing staging: its worker can still have the archive open.
    download_.reset();
    if (disk_) {
        disconnect(disk_.get(), nullptr, this, nullptr);
        disk_->requestInterruption();
        disk_->wait();
    }
    if (staging_ && !staging_->remove())
        qCCritical(cacheLog) << "Could not remove interrupted download:" << staging_->path();
}

QString BundleCache::cachedPath(const Bundle &bundle) const {
    if (!safeBundleName(bundle.key) || !safeRelativePath(bundle.entryPoint))
        return {};
    const QString path = QDir(root_).filePath(bundle.key);
    const QFileInfo directory(path);
    const QFileInfo entry(QDir(path).filePath(bundle.entryPoint));
    return directory.isDir() && !directory.isSymLink() && entry.isFile() && !entry.isSymLink() ? path : QString();
}

bool BundleCache::fetch(const Bundle &bundle) {
    Q_ASSERT(thread() == QThread::currentThread());
    if (busy() || shutdownRequested_) {
        qCWarning(cacheLog) << "Cache busy or shutting down; fetch rejected:" << bundle.key;
        return false;
    }
    bundle_ = bundle;
    cancelled_ = false;
    // Queue entry so even validation and cache hits have the same asynchronous API.
    QMetaObject::invokeMethod(this, &BundleCache::beginDownload, Qt::QueuedConnection);
    setState(State::Downloading);
    return true;
}

void BundleCache::beginDownload() {
    if (cancelled_) {
        complete(Result::Cancelled, tr("Cancelled."));
        return;
    }
    if (!safeBundleName(bundle_.key) || !safeRelativePath(bundle_.entryPoint) || bundle_.expandedBytes < 0) {
        complete(Result::InvalidRequest, tr("Invalid bundle key, entry path or expanded byte count."));
        return;
    }
    const QString cached = cachedPath(bundle_);
    if (!cached.isEmpty()) {
        complete(Result::Ready, tr("Using the cached bundle."), cached);
        return;
    }
    if (QFileInfo::exists(QDir(root_).filePath(bundle_.key))) {
        complete(Result::FileError,
                 tr("An incomplete or unexpected cache entry exists. Clear it explicitly before retrying."));
        return;
    }
    if (makeStaging())
        download_->start(bundle_.url, staging_->filePath("download.archive"));
}

bool BundleCache::makeStaging() {
    if (!QDir().mkpath(root_)) {
        complete(Result::FileError, tr("Could not create cache directory: %1").arg(root_));
        return false;
    }
    staging_ = std::make_unique<QTemporaryDir>(QDir(root_).filePath(".partial-XXXXXX"));
    if (!staging_->isValid()) {
        const QString error = staging_->errorString();
        staging_.reset();
        complete(Result::FileError, error);
        return false;
    }
    // Cleanup is explicit and checked, normally on the disk worker.
    staging_->setAutoRemove(false);
    return true;
}

void BundleCache::downloaded(Download::Result result, const QString &detail) {
    if (result == Download::Result::Saved) {
        if (cancelled_) {
            startDiskTask(false, Result::Cancelled, tr("Cancelled. No incomplete bundle was kept."));
            return;
        }
        startDiskTask(true, Result::Ready, {});
        return;
    }
    Result outcome = Result::TransferError;
    if (result == Download::Result::Cancelled)
        outcome = Result::Cancelled;
    else if (result == Download::Result::InvalidRequest)
        outcome = Result::InvalidRequest;
    else if (result == Download::Result::FileError)
        outcome = Result::FileError;
    startDiskTask(false, outcome, detail);
}

void BundleCache::startDiskTask(bool extract, Result result, const QString &detail) {
    disk_ = std::make_unique<ArchiveTask>(std::move(staging_), bundle_, QDir(root_).filePath(bundle_.key), extract,
                                          result, detail);
    connect(disk_.get(), &ArchiveTask::progress, this, [this](qint64 bytes, qint64 total) {
        if (state_ == State::Extracting)
            emit progress(bytes, total);
    });
    connect(disk_.get(), &ArchiveTask::cleaning, this, [this] { setState(State::Cleaning); });
    connect(disk_.get(), &QThread::finished, this, &BundleCache::diskFinished);
    disk_->start();
    if (cancelled_ && extract)
        disk_->requestInterruption();
    setState(extract ? State::Extracting : result == Result::Cleared ? State::Clearing : State::Cleaning);
}

void BundleCache::diskFinished() {
    disk_->wait();
    const auto result = disk_->result();
    const QString detail = disk_->detail();
    const QString path = disk_->publishedPath();
    disk_.reset();
    complete(result, detail, path);
}

bool BundleCache::clear(const QString &key) {
    Q_ASSERT(thread() == QThread::currentThread());
    if (busy() || shutdownRequested_) {
        qCWarning(cacheLog) << "Cache busy or shutting down; clear rejected:" << key;
        return false;
    }
    bundle_ = Bundle{key, {}, {}, 0};
    cancelled_ = false;
    QMetaObject::invokeMethod(
        this,
        [this] {
            if (!safeBundleName(bundle_.key)) {
                complete(Result::InvalidRequest, tr("Invalid cache key."));
                return;
            }
            const QString path = QDir(root_).filePath(bundle_.key);
            if (!QFileInfo::exists(path)) {
                complete(Result::Cleared, tr("No local bundle to clear."));
                return;
            }
            if (!makeStaging())
                return;
            // Remove the entry from view before deleting its files. Never expose half a bundle.
            if (!QDir().rename(path, staging_->filePath("removed"))) {
                startDiskTask(false, Result::FileError, tr("Could not move the cache entry for removal: %1").arg(path));
                return;
            }
            startDiskTask(false, Result::Cleared, tr("Local bundle cleared. Opening it again requires a download."));
        },
        Qt::QueuedConnection);
    setState(State::Clearing);
    return true;
}

void BundleCache::cancel() {
    Q_ASSERT(thread() == QThread::currentThread());
    if (cancelled_ || (state_ != State::Downloading && state_ != State::Extracting))
        return;
    cancelled_ = true;
    qCInfo(cacheLog) << "Cancellation requested:" << bundle_.key;
    if (disk_)
        disk_->requestInterruption();
    download_->cancel();
    setState(State::Cleaning);
}

void BundleCache::setState(State state) {
    if (state_ == state)
        return;
    state_ = state;
    qCInfo(cacheLog) << "State:" << bundle_.key << state;
    emit stateChanged(state);
}

void BundleCache::shutdown() {
    Q_ASSERT(thread() == QThread::currentThread());
    if (shutdownRequested_)
        return;
    shutdownRequested_ = true;
    qCInfo(cacheLog) << "Cache shutdown requested.";
    cancel();
    if (!busy())
        download_->shutdown();
}

void BundleCache::complete(Result result, const QString &detail, const QString &path) {
    qCInfo(cacheLog) << "Completed:" << bundle_.key << "code=" << static_cast<int>(result) << detail << path;
    // finished is the operation boundary; clients may start their next request there.
    const QPointer<BundleCache> alive(this);
    setState(State::Idle);
    if (!alive)
        return;
    if (shutdownRequested_)
        download_->shutdown();
    emit finished(result, detail, path);
}

} // namespace sciencefiles
