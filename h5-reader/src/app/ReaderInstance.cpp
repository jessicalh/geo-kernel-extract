#include "ReaderInstance.h"

#include <QDeadlineTimer>
#include <QDir>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLocalSocket>
#include <QStandardPaths>

namespace h5reader::app {
namespace {
constexpr qint64 MaxRequestBytes = 64 * 1024;
constexpr int ForwardTimeoutMs = 20000;
}

ReaderInstance::ReaderInstance(QObject* parent, const QString& root)
    : QObject(parent),
      root_(QDir(root.isEmpty() ? QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) : root)
                .absolutePath()),
      serverName_(QDir(root_).filePath(QStringLiteral("reader-instance"))),
      lock_(QDir(root_).filePath(QStringLiteral("reader-instance.lock"))), server_(this) {
    lock_.setStaleLockTime(0);
    server_.setSocketOptions(QLocalServer::UserAccessOption);
    connect(&server_, &QLocalServer::newConnection, this, [this] {
        while (auto* socket = server_.nextPendingConnection()) {
            socket->setReadBufferSize(MaxRequestBytes + 1);
            connect(socket, &QLocalSocket::disconnected, socket, &QObject::deleteLater);
            connect(socket, &QLocalSocket::readyRead, this, [this, socket] { readRequest(socket); });
            readRequest(socket);
        }
    });
}

ReaderInstance::~ReaderInstance() {
    shutdown();
}

ReaderInstance::Result ReaderInstance::fail(const QString& message) {
    error_ = message;
    return Result::Failed;
}

ReaderInstance::Result ReaderInstance::start(const QString& path, bool allowForward) {
    if (!QDir().mkpath(root_))
        return fail(tr("Cannot create Reader's instance directory: %1").arg(root_));
    if (!lock_.tryLock(0)) {
        if (lock_.error() != QLockFile::LockFailedError)
            return fail(tr("Cannot lock Reader's instance directory: %1").arg(root_));
        if (!allowForward)
            return fail(tr("Reader is already running. A second --rest launch cannot change its REST server."));
        return forward(path);
    }

    // Only the legitimate owner may remove a socket left by a crashed primary.
    QLocalServer::removeServer(serverName_);
    if (!server_.listen(serverName_))
        return fail(tr("Cannot listen for Reader file requests: %1").arg(server_.errorString()));
    return Result::Primary;
}

ReaderInstance::Result ReaderInstance::forward(const QString& path) {
    const QString absolutePath = path.isEmpty() ? QString() : QFileInfo(path).absoluteFilePath();
    const QByteArray request = QJsonDocument(QJsonObject{{"path", absolutePath}}).toJson(QJsonDocument::Compact) + '\n';
    if (request.size() > MaxRequestBytes)
        return fail(tr("The file path is too long to send to Reader."));

    QDeadlineTimer deadline(ForwardTimeoutMs);
    QLocalSocket socket;
    socket.connectToServer(serverName_);
    if (!socket.waitForConnected(static_cast<int>(deadline.remainingTime())))
        return fail(tr("Reader is busy or closing. Open the file from the running window."));
    if (socket.write(request) != request.size())
        return fail(tr("Cannot send the file request to Reader: %1").arg(socket.errorString()));
    while (socket.bytesToWrite() > 0) {
        if (!socket.waitForBytesWritten(static_cast<int>(deadline.remainingTime())))
            return fail(tr("Cannot send the file request to Reader: %1").arg(socket.errorString()));
    }
    if (socket.bytesAvailable() == 0 && !socket.waitForReadyRead(static_cast<int>(deadline.remainingTime())))
        return fail(tr("Reader is busy or closing. Open the file from the running window."));
    if (socket.read(1) != QByteArrayLiteral("1"))
        return fail(tr("The running Reader rejected the file request."));
    return Result::Forwarded;
}

void ReaderInstance::readRequest(QLocalSocket* socket) {
    if (stopping_ || socket->state() != QLocalSocket::ConnectedState || socket->bytesAvailable() > MaxRequestBytes) {
        socket->abort();
        return;
    }
    if (!socket->canReadLine())
        return;
    const auto document = QJsonDocument::fromJson(socket->readLine(MaxRequestBytes + 1));
    const auto object = document.object();
    if (!document.isObject() || object.size() != 1 || !object.value("path").isString()) {
        socket->abort();
        return;
    }
    const QString path = object.value("path").toString();
    disconnect(socket, &QLocalSocket::readyRead, this, nullptr);
    // Open only after writing the receipt ACK, not after a disconnected sender gives up.
    connect(socket, &QLocalSocket::bytesWritten, socket, [this, socket, path] {
        socket->disconnectFromServer();
        if (!stopping_)
            emit openRequested(path);
    }, Qt::SingleShotConnection);
    if (socket->write("1", 1) != 1) {
        socket->abort();
        return;
    }
    socket->flush();
}

void ReaderInstance::shutdown() {
    stopping_ = true;
    server_.close();
    for (auto* socket : server_.findChildren<QLocalSocket*>()) {
        socket->disconnect(this);
        socket->abort();
    }
}

} // namespace h5reader::app
