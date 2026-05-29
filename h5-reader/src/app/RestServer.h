// RestServer — read-only-by-default HTTP control surface for h5-reader,
// used by the pytest-driven REST smoke suite (tests/rest/) in place of the
// in-binary --dashboard-path-smoke / --camera-plane-lock-smoke runners.
//
// Lifecycle: created and owned by ReaderMainWindow when --rest <port> is
// passed on the CLI. All routes execute on the GUI thread (QHttpServer's
// default handler thread is the thread the server was constructed on);
// handlers may directly read/mutate model + scene state without marshalling.
//
// Binding: loopback only (QHostAddress::LocalHost). No authentication. This
// is a test fixture, not a public API.
//
// JSON: QJsonDocument / QJsonObject / QJsonArray only. No third-party JSON.
//
// Endpoints (initial inventory; expand only when a test needs more):
//
//   GET    /health                       → {"ok": true, "version": "..."}
//   GET    /protein/atoms                → {"count": N}
//   GET    /frame/current                → {"frame": int, "time_ps": float, "count": int}
//   POST   /frame/set                    → 204 (body: {"frame": int})
//   GET    /selection                    → {"atoms": [...], "focus": int|null, "count": int}
//   POST   /selection/pick               → 204 (body: {"atom": int, "modifiers": "none"|"shift"})
//   POST   /selection/clear              → 204
//   GET    /plane-lock                   → {"active": bool, "atoms": [...]|null}
//   POST   /plane-lock/enable            → 204 or 409 (body: {"atoms": [a,b,c]})
//   POST   /plane-lock/disable           → 204
//   GET    /scene/camera                 → {"focal":[x,y,z], "position":[x,y,z], "view_up":[x,y,z], "direction":[x,y,z]}
//   GET    /dashboard/signals            → [{"id": uuid, "descriptor_id": ..., "modes": [...], "label": ...}, ...]
//   POST   /screenshot                   → image/png (body: {"target":"scene"|"window"})

#pragma once

#include "../diagnostics/ObjectCensus.h"

#include <QHostAddress>
#include <QObject>
#include <QPointer>

#include <memory>

class QHttpServer;
class QWidget;

namespace h5reader {
namespace io { struct QtLoadResult; }
namespace model {
class AtomSelection;
class DashboardPanelModel;
class DashboardSignalModel;
class TrajectorySignalCatalog;
}
}

namespace h5reader::app {

class MoleculeScene;
class QtPlaybackController;

class RestServer final : public QObject {
    Q_OBJECT

public:
    explicit RestServer(QObject* parent = nullptr);
    ~RestServer() override;

    // Wire the dependencies this server reads / mutates. Must be called
    // before listen(); all pointers are stored as QPointer / raw and the
    // server never outlives them (it is owned by ReaderMainWindow).
    void setContext(MoleculeScene* scene,
                    model::AtomSelection* selection,
                    model::DashboardSignalModel* signalModel,
                    model::DashboardPanelModel* panelModel,
                    const model::TrajectorySignalCatalog* catalog,
                    QtPlaybackController* playback,
                    io::QtLoadResult* loaded,
                    QWidget* mainWindow);

    // Bind to QHostAddress::LocalHost on the requested port. Port 0 asks
    // the kernel to pick. Returns the actually-bound port, or 0 on failure.
    // Emits an info log line `H5READER_REST_PORT=<port>` to stderr on
    // success so the pytest fixture can scrape it.
    quint16 listen(quint16 port);

private:
    void registerRoutes();

    std::unique_ptr<QHttpServer>                server_;
    QPointer<MoleculeScene>                     scene_;
    QPointer<model::AtomSelection>              selection_;
    QPointer<model::DashboardSignalModel>       signalModel_;
    QPointer<model::DashboardPanelModel>        panelModel_;
    const model::TrajectorySignalCatalog*       catalog_ = nullptr;
    QPointer<QtPlaybackController>              playback_;
    io::QtLoadResult*                           loaded_ = nullptr;
    QPointer<QWidget>                           mainWindow_;
    bool                                        contextSet_ = false;
};

}  // namespace h5reader::app
