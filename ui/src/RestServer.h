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
    QJsonObject cmdGetLog(const QJsonObject& cmd);
    QJsonObject cmdExportFeatures(const QJsonObject& cmd);

    // ─── Programmatic introspection ──────────────────────────────────
    //
    // General-purpose endpoints that return per-atom and per-residue
    // typed data as JSON. Designed for scripting and test clients
    // (ui/tests/), not for any single test — they expose the same
    // typed object surface the Atom Inspector shows, in machine-readable
    // form. No test-specific shortcuts here.
    //
    //   atom_dump   {"atom":N}              → full inspector tree
    //   list_atoms  {filter?: {…}}          → concise per-atom records
    //
    QJsonObject cmdAtomDump(const QJsonObject& cmd);
    QJsonObject cmdListAtoms(const QJsonObject& cmd);

    // Polite shutdown. cmdQuit() flips `shutdown_after_reply_` and
    // returns the reply JSON; the dispatch caller (onReadyRead) sends
    // the reply, then drives the socket through `disconnectFromHost`
    // — Qt drains the write buffer before emitting `disconnected`, and
    // the same socket's `disconnected` signal is connected to
    // QCoreApplication::quit (UniqueConnection, so first disconnect
    // wins). Pure signal/slot chain: no QTimer::singleShot, no sync
    // waitForBytesWritten, no raw socket flush gymnastics.
    QJsonObject cmdQuit();

    void sendResponse(QTcpSocket* socket, const QJsonObject& response);

    QTcpServer* server_;
    MainWindow* mainWindow_;
    QList<QTcpSocket*> clients_;
    quint16 actualPort_ = 0;

    // Set by cmdQuit; observed by onReadyRead after sendResponse so the
    // graceful shutdown sequence runs only after the reply is queued
    // on the socket's write buffer.
    bool shutdown_after_reply_ = false;
};
