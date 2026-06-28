// RestServer — read-only-by-default HTTP control surface for h5-reader,
// used by the pytest-driven REST smoke suite (tests/rest/) in place of the
// in-binary --dashboard-path-smoke / --camera-plane-lock-smoke runners.
//
// Lifecycle: created and owned by ReaderMainWindow when --rest <port> is
// passed on the CLI. All routes execute on the GUI thread (QHttpServer's
// default handler thread is the thread the server was constructed on);
// handlers may directly read/mutate model + scene state without marshalling.
//
// Binding: loopback only (QHostAddress::LocalHost). No authentication; this is
// a local operator/API harness, not a remote public service.
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
//   POST   /selection/atoms              → 204 (body: {"atoms": [int, ...]})  bulk replace
//   POST   /selection/clear              → 204
//   POST   /selection/instrument         → 204 (body: {"enabled": bool, "focus_only": bool})  marker preset
//   POST   /docks/visible                → 204 (body: {"visible": bool})  hide/restore docks
//   GET    /transform                    → {"kind": "...", "reference_frame": int, "subset_atoms": [...],
//                                            "subset_size": int, "window": int}
    //   POST   /transform                    → 204 (body: {"kind": "all_atom_fit"|"backbone_fit",
    //                                                       "reference_frame": int seed/anchor, "subset_atoms": [int, ...],
    //                                                       "backbone_only": bool})
    //   POST   /transform/smoothing          → 204 (body: {"window": int})  rotation-only
//   GET    /plane-lock                   → {"active": bool, "atoms": [...]|null}
//   POST   /plane-lock/enable            → 204 or 409 (body: {"atoms": [a,b,c]})
//   POST   /plane-lock/disable           → 204
//   GET    /scene/camera                 → {"focal":[x,y,z], "position":[x,y,z], "view_up":[x,y,z], "direction":[x,y,z]}
//   GET    /dashboard/signals            → [{"id": uuid, "descriptor_id": ..., "modes": [...], "label": ...}, ...]
//   GET    /dashboard/state              → selected metrics + dock/render state
//   GET    /dashboard/strip/series       → active strip ChannelBuffer values/valid arrays
//   GET    /dashboard/picker             → live Metric Picker selector state
//   POST   /dashboard/picker/open        → live selector state (body: {"atom": int} optional)
//   POST   /dashboard/picker/add         → user-equivalent Add Signal from current picker state
//   POST   /dashboard/metric             → {"id": uuid, "added_refs": int}
//                                            body: {"descriptor_id": "...", "anchor": {...}, "modes": [...]}
//   POST   /dashboard/metric/remove      → {"removed": true}  body: {"id": uuid}
//   POST   /dashboard/metric/mode        → {"modes": [...]}  body: {"id": uuid, "mode": "...", "enabled": bool}
//   POST   /dashboard/dock               → {"visible": bool, "width": int}  body: {"visible": bool}
//   GET    /api/interface                -> namespace + route contract map
//   POST   /api/screenshot               -> image/png (body: {"target":"scene"|"window",
//                                                            "scale": int})
//   POST   /diagnostics/screenshot       -> image/png (adds target:"widget" and force_render)

#pragma once

#include "MoleculeScene.h"

#include "../diagnostics/ObjectCensus.h"

#include <QHostAddress>
#include <QObject>
#include <QPointer>

#include <cstddef>
#include <memory>
#include <optional>

class QHttpServer;
class QWidget;
class QTcpSocket;

namespace h5reader {
namespace io { struct QtLoadResult; }
namespace model {
class AtomSelection;
class DashboardPanelModel;
class DashboardSignalModel;
class TrajectorySignalCatalog;
class TransformedConformation;
}
}

namespace h5reader::app {

class DashboardSelectionController;
class DashboardDisplayController;
class QtPlaybackController;
class ReaderMainWindow;
class AngleCollarActor;
class AtomTrackOverlay;
class HeroshotButterflyOverlay;
class TensorGhostTrail;

class RestServer final : public QObject {
    Q_OBJECT

public:
    explicit RestServer(QObject* parent = nullptr);
    ~RestServer() override;

    // Wire the dependencies this server reads / mutates. Must be called
    // before listen(); all pointers are stored as QPointer / raw and the
    // server never outlives them (it is owned by ReaderMainWindow).
    //
    // `readerWindow` is the typed ReaderMainWindow* (also stored as
    // `mainWindow` as a QWidget* for the screenshot path). It's separate
    // because the dock-visible endpoint needs the typed ReaderMainWindow
    // surface; the screenshot path only needs QWidget::grab().
    //
    // `transformed` is the TransformedConformation wrapping the loader's
    // Conformation — handed to POST /transform so the harness can flip
    // the rigid-body transform mode at runtime.
    void setContext(MoleculeScene* scene,
                    model::AtomSelection* selection,
                    model::DashboardSignalModel* signalModel,
                    model::DashboardPanelModel* panelModel,
                    DashboardSelectionController* selectionController,
                    DashboardDisplayController* dashboardController,
                    const model::TrajectorySignalCatalog* catalog,
                    QtPlaybackController* playback,
                    io::QtLoadResult* loaded,
                    QWidget* mainWindow,
                    ReaderMainWindow* readerWindow = nullptr,
                    model::TransformedConformation* transformed = nullptr);

    // Bind to QHostAddress::LocalHost on the requested port. Port 0 asks
    // the kernel to pick. Returns the actually-bound port, or 0 on failure.
    // Emits an info log line `H5READER_REST_PORT=<port>` to stderr on
    // success so the pytest fixture can scrape it.
    quint16 listen(quint16 port);

private:
    void registerRoutes();

    std::unique_ptr<QHttpServer>                server_;
    // Most recently accepted REST socket (captured at connection time). POST
    // /shutdown waits on its flush event to exit cleanly — see RestServer.cpp.
    QPointer<QTcpSocket>                        activeSocket_;
    QPointer<MoleculeScene>                     scene_;
    QPointer<model::AtomSelection>              selection_;
    QPointer<model::DashboardSignalModel>       signalModel_;
    QPointer<model::DashboardPanelModel>        panelModel_;
    QPointer<DashboardSelectionController>      selectionController_;
    QPointer<DashboardDisplayController>        dashboardController_;
    const model::TrajectorySignalCatalog*       catalog_ = nullptr;
    QPointer<QtPlaybackController>              playback_;
    io::QtLoadResult*                           loaded_ = nullptr;
    QPointer<QWidget>                           mainWindow_;
    QPointer<ReaderMainWindow>                  readerWindow_;
    QPointer<model::TransformedConformation>    transformed_;
    // Resthero layer (transient figure FX, never part of the reader UI): the
    // tensor ghost trail built on demand by POST /resthero/ghost_trail.
    // Rebuilt against the live scene renderer each call; cleared by
    // /resthero/clear.
    std::unique_ptr<AtomTrackOverlay>            heroshotAtomTrack_;
    std::unique_ptr<HeroshotButterflyOverlay>    heroshotButterfly_;
    std::unique_ptr<TensorGhostTrail>           heroshotTrail_;
    std::unique_ptr<AngleCollarActor>           heroshotAngleCollar_;
    std::optional<bool>                         heroshotMeasurementVisibleBefore_;
    std::optional<MoleculeScene::MoleculeStyle> heroshotMoleculeStyleBefore_;
    std::optional<std::size_t>                  heroshotFieldRingBefore_;
    bool                                        heroshotFieldRingWasSet_ = false;
    bool                                        contextSet_ = false;
};

}  // namespace h5reader::app
