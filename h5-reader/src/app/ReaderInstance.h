#pragma once

#include <QLocalServer>
#include <QLockFile>
#include <QObject>

class QLocalSocket;

namespace h5reader::app {

class ReaderInstance final : public QObject {
    Q_OBJECT
public:
    enum class Result { Primary, Forwarded, Failed };
    Q_ENUM(Result)

    explicit ReaderInstance(QObject* parent = nullptr, const QString& root = {});
    ~ReaderInstance() override;

    // Call once, after setting the application identity and before creating the window.
    // Forwarded acknowledges receipt, not successful loading. Empty path means activate.
    Result start(const QString& path, bool allowForward = true);
    QString errorString() const { return error_; }

public slots:
    // Stop handoffs during window cleanup, retaining ownership until destruction.
    void shutdown();

signals:
    void openRequested(const QString& path);

private:
    Result fail(const QString& message);
    Result forward(const QString& path);
    void readRequest(QLocalSocket* socket);

    QString root_;
    QString serverName_;
    QLockFile lock_;
    QLocalServer server_;
    QString error_;
    bool stopping_ = false;
};

} // namespace h5reader::app
