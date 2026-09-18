#pragma once

#include <QObject>
#include <QSslConfiguration>
#include <QThread>
#include <QUrl>
#include <chrono>

namespace sciencefiles {

class DownloadWorker;

// One transfer at a time. Keep this object in its creating Qt thread.
class Download final : public QObject {
    Q_OBJECT
public:
    enum class Result {
        Saved = 0,
        Cancelled = 2001,
        NetworkError = 2002,
        HttpError = 2003,
        FileError = 2004,
        InvalidRequest = 2005
    };
    Q_ENUM(Result)

    explicit Download(QObject *parent = nullptr,
                      const QSslConfiguration &tls = QSslConfiguration::defaultConfiguration(),
                      std::chrono::milliseconds inactivity = std::chrono::seconds(60));
    ~Download() override;

    // False means active or shutting down; accepted starts finish asynchronously.
    bool start(const QUrl &url, const QString &destination);
    bool active() const { return active_; }
    bool isShutdown() const { return shutdownComplete_; }

public slots:
    void cancel();
    void shutdown();

signals:
    void progress(qint64 received, qint64 total);
    void finished(Download::Result result, const QString &detail);
    void shutdownFinished();

private:
    QThread workerThread_;
    DownloadWorker *worker_;
    bool active_ = false;
    bool shutdownRequested_ = false;
    bool shutdownComplete_ = false;
};

} // namespace sciencefiles
