#pragma once
//
// REST-like TCP command server for programmatic viewer control.
//
// Listens on localhost:9147. Accepts newline-delimited JSON commands,
// dispatches to MainWindow slots, returns JSON responses.
//
// All VTK operations happen on the main thread via queued invocations.
//

#include <QObject>
#include <QTcpServer>
#include <QTcpSocket>
#include <QJsonObject>

class MainWindow;

class RestServer : public QObject {
    Q_OBJECT
public:
    explicit RestServer(MainWindow* mainWindow, quint16 port = 9147,
                        QObject* parent = nullptr);

    quint16 actualPort() const { return actualPort_; }

private slots:
    void onNewConnection();
    void onReadyRead();
    void onDisconnected();

private:
    QJsonObject dispatch(const QJsonObject& cmd);

    // Command handlers — each returns a result JSON object
    QJsonObject cmdStatus();
    QJsonObject cmdLoadPdb(const QJsonObject& cmd);
    QJsonObject cmdLoadProteinDir(const QJsonObject& cmd);
    QJsonObject cmdSetOverlay(const QJsonObject& cmd);
    QJsonObject cmdSetRenderMode(const QJsonObject& cmd);
    QJsonObject cmdScreenshot(const QJsonObject& cmd);
    QJsonObject cmdOrbit(const QJsonObject& cmd);
    QJsonObject cmdResetView(const QJsonObject& cmd);
    QJsonObject cmdSetCalculators(const QJsonObject& cmd);
    QJsonObject cmdShowRings(const QJsonObject& cmd);
    QJsonObject cmdShowBonds(const QJsonObject& cmd);
    QJsonObject cmdShowButterfly(const QJsonObject& cmd);
    QJsonObject cmdSetGlyphScale(const QJsonObject& cmd);
    QJsonObject cmdSetOpacity(const QJsonObject& cmd);
    QJsonObject cmdSetIsoThreshold(const QJsonObject& cmd);
    QJsonObject cmdShowFieldGrid(const QJsonObject& cmd);
    QJsonObject cmdGetCamera(const QJsonObject& cmd);
    QJsonObject cmdSetCamera(const QJsonObject& cmd);
    QJsonObject cmdLookAtRing(const QJsonObject& cmd);
    QJsonObject cmdLookAtAtom(const QJsonObject& cmd);
    QJsonObject cmdListRings(const QJsonObject& cmd);

    void sendResponse(QTcpSocket* socket, const QJsonObject& response);

    QTcpServer* server_;
    MainWindow* mainWindow_;
    QList<QTcpSocket*> clients_;
    quint16 actualPort_ = 0;
};
