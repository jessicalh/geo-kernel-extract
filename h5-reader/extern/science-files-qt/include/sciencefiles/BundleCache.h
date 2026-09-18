#pragma once

#include <sciencefiles/Download.h>
#include <QTemporaryDir>
#include <QUrl>
#include <memory>

namespace sciencefiles {

struct Bundle {
    QString key;
    QUrl url;
    QString entryPoint;
    qint64 expandedBytes = 0;
};

class ArchiveTask;

// One caller owns a dedicated cache root. Use this object only in its Qt thread.
class BundleCache final : public QObject {
    Q_OBJECT
public:
    enum class State { Idle, Downloading, Extracting, Clearing, Cleaning };
    Q_ENUM(State)
    enum class Result {
        Ready = 0,
        Cleared = 1,
        Cancelled = 2101,
        InvalidRequest = 2102,
        TransferError = 2103,
        ArchiveError = 2104,
        FileError = 2105
    };
    Q_ENUM(Result)

    explicit BundleCache(const QString &root, QObject *parent = nullptr,
                         const QSslConfiguration &tls = QSslConfiguration::defaultConfiguration());
    ~BundleCache() override;

    // False means busy or shutting down; accepted operations finish asynchronously.
    bool fetch(const Bundle &bundle);
    bool clear(const QString &key);
    State state() const { return state_; }
    bool busy() const { return state_ != State::Idle; }
    bool isShutdown() const { return shutdownComplete_; }
    QString cachedPath(const Bundle &bundle) const;

public slots:
    void cancel();
    void shutdown();

signals:
    void stateChanged(sciencefiles::BundleCache::State state);
    void progress(qint64 completedBytes, qint64 totalBytes);
    void finished(sciencefiles::BundleCache::Result result, const QString &detail, const QString &localDirectory);
    void shutdownFinished();

private:
    void beginDownload();
    void downloaded(Download::Result result, const QString &detail);
    void startDiskTask(bool extract, Result result, const QString &detail);
    void diskFinished();
    void setState(State state);
    void complete(Result result, const QString &detail, const QString &path = {});
    bool makeStaging();

    QString root_;
    std::unique_ptr<Download> download_;
    Bundle bundle_;
    State state_ = State::Idle;
    std::unique_ptr<QTemporaryDir> staging_;
    std::unique_ptr<ArchiveTask> disk_;
    bool cancelled_ = false;
    bool shutdownRequested_ = false;
    bool shutdownComplete_ = false;
};

} // namespace sciencefiles
