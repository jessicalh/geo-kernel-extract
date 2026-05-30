#include "RestServer.h"

#include "CameraAnchorHelper.h"
#include "CameraComposer.h"
#include "MeasurementOverlay.h"
#include "MoleculeScene.h"
#include "QtPlaybackController.h"
#include "ReaderMainWindow.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/StructuredLogger.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/QtProteinLoader.h"
#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/DashboardPanelModel.h"
#include "../model/DashboardSignal.h"
#include "../model/DashboardSignalModel.h"
#include "../model/QtProtein.h"
#include "../model/TrajectorySignalCatalog.h"
#include "../model/TransformedConformation.h"

#include <QBuffer>
#include <QByteArray>
#include <QCoreApplication>
#include <QHttpServer>
#include <QHttpServerRequest>
#include <QHttpServerResponse>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QLoggingCategory>
#include <QPixmap>
#include <QString>
#include <QTimer>
#include <QUuid>
#include <QVariant>
#include <QWidget>

#include <vtkCamera.h>
#include <vtkPNGWriter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkWindowToImageFilter.h>

#include <cstddef>
#include <cstdio>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cRest, "h5reader.rest")

constexpr const char* kMimeJson = "application/json";
constexpr const char* kMimePng = "image/png";

QHttpServerResponse jsonResponse(const QJsonObject& obj,
                                 QHttpServerResponse::StatusCode code = QHttpServerResponse::StatusCode::Ok) {
    return QHttpServerResponse(obj, code);
}

QHttpServerResponse jsonResponse(const QJsonArray& arr,
                                 QHttpServerResponse::StatusCode code = QHttpServerResponse::StatusCode::Ok) {
    return QHttpServerResponse(arr, code);
}

QHttpServerResponse errorResponse(const QString& message,
                                  QHttpServerResponse::StatusCode code) {
    return jsonResponse(QJsonObject{{"error", message}}, code);
}

QJsonObject parseJsonBody(const QHttpServerRequest& request, bool* ok) {
    *ok = false;
    QJsonParseError err{};
    const QJsonDocument doc = QJsonDocument::fromJson(request.body(), &err);
    if (err.error != QJsonParseError::NoError || !doc.isObject())
        return {};
    *ok = true;
    return doc.object();
}

QJsonArray vec3FromRaw(const double raw[3]) {
    return QJsonArray{raw[0], raw[1], raw[2]};
}

QJsonObject anchorToJson(const model::SignalAnchor& anchor) {
    using namespace model;
    QJsonObject out;
    if (std::holds_alternative<NoneAnchor>(anchor)) {
        out["kind"] = "none";
    } else if (const auto* a = std::get_if<AtomAnchor>(&anchor)) {
        out["kind"] = "atom"; out["atom"] = static_cast<qint64>(a->atom);
    } else if (const auto* r = std::get_if<ResidueAnchor>(&anchor)) {
        out["kind"] = "residue"; out["residue"] = static_cast<qint64>(r->residue);
    } else if (const auto* t = std::get_if<AtomTupleAnchor>(&anchor)) {
        QJsonArray atoms;
        for (auto a : t->atoms) atoms.append(static_cast<qint64>(a));
        out["kind"] = "atom_tuple"; out["atoms"] = atoms;
    } else if (const auto* b = std::get_if<BondAnchor>(&anchor)) {
        out["kind"] = "bond"; out["bond"] = static_cast<qint64>(b->bond);
    } else if (const auto* v = std::get_if<BondVectorAnchor>(&anchor)) {
        out["kind"] = "bond_vector";
        out["residue"] = static_cast<qint64>(v->residue);
        out["kind_id"] = static_cast<qint64>(v->kind);
    } else if (const auto* r = std::get_if<RingAnchor>(&anchor)) {
        out["kind"] = "ring"; out["ring"] = static_cast<qint64>(r->ring);
    } else if (const auto* r = std::get_if<AromaticRingAnchor>(&anchor)) {
        out["kind"] = "aromatic_ring"; out["ring"] = static_cast<qint64>(r->ring);
    } else if (const auto* r = std::get_if<SaturatedRingAnchor>(&anchor)) {
        out["kind"] = "saturated_ring"; out["ring"] = static_cast<qint64>(r->ring);
    } else if (const auto* p = std::get_if<RingContributionPairAnchor>(&anchor)) {
        out["kind"] = "ring_contribution_pair"; out["pair"] = static_cast<qint64>(p->pair);
    } else if (const auto* m = std::get_if<RingMembershipAnchor>(&anchor)) {
        out["kind"] = "ring_membership"; out["membership"] = static_cast<qint64>(m->membership);
    } else if (const auto* p = std::get_if<MutationMatchPairAnchor>(&anchor)) {
        out["kind"] = "mutation_match_pair"; out["pair"] = static_cast<qint64>(p->pair);
    } else if (std::holds_alternative<ProteinAnchor>(anchor)) {
        out["kind"] = "protein";
    } else if (std::holds_alternative<SystemAnchor>(anchor)) {
        out["kind"] = "system";
    } else if (std::holds_alternative<EventAnchor>(anchor)) {
        out["kind"] = "event";
    }
    return out;
}

// Capture the current VTK render window into a PNG byte buffer.
//
// forceRender (default true, back-compat with prior calls):
//   true  → leave vtkWindowToImageFilter's default ShouldRerenderOn — forces
//           a fresh Render() before reading pixels, so the snapshot reflects
//           the live scene state.
//   false → call ShouldRerenderOff() before Update() — read whatever pixels
//           are currently in the framebuffer. The right mode for the
//           paint-cycle-inversion experiment (VIEWPORT_OBSERVATIONS §5b);
//           lets the harness distinguish "the synchronous Render reached
//           the back buffer" from "we read the post-render FBO".
//
// Thread: VTK render/read must happen on the GUI thread. ASSERT_THREAD against
// the scene's affinity catches a future regression where a route handler
// might be routed off the GUI thread by QHttpServer.
QByteArray captureScenePng(MoleculeScene* scene, bool forceRender = true) {
    if (!scene || !scene->Renderer() || !scene->Renderer()->GetRenderWindow())
        return {};
    ASSERT_THREAD(scene);

    auto w2i = vtkSmartPointer<vtkWindowToImageFilter>::New();
    w2i->SetInput(scene->Renderer()->GetRenderWindow());
    w2i->SetInputBufferTypeToRGB();
    w2i->ReadFrontBufferOff();
    if (!forceRender)
        w2i->ShouldRerenderOff();
    w2i->Update();

    auto writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetWriteToMemory(1);
    writer->SetInputConnection(w2i->GetOutputPort());
    writer->Write();

    vtkUnsignedCharArray* memArr = writer->GetResult();
    if (!memArr || memArr->GetNumberOfTuples() == 0)
        return {};

    return QByteArray(reinterpret_cast<const char*>(memArr->GetPointer(0)),
                      static_cast<int>(memArr->GetNumberOfTuples()));
}

// Capture the top-level window via Qt (works for the whole UI, includes
// dock widgets and dashboard strip — what a human screenshots).
QByteArray captureWindowPng(QWidget* window) {
    if (!window)
        return {};
    const QPixmap pix = window->grab();
    if (pix.isNull())
        return {};
    QByteArray buf;
    QBuffer buffer(&buf);
    buffer.open(QIODevice::WriteOnly);
    if (!pix.save(&buffer, "PNG"))
        return {};
    return buf;
}
}  // namespace

RestServer::RestServer(QObject* parent)
    : QObject(parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("RestServer"));
}

RestServer::~RestServer() = default;

void RestServer::setContext(MoleculeScene* scene,
                            model::AtomSelection* selection,
                            model::DashboardSignalModel* signalModel,
                            model::DashboardPanelModel* panelModel,
                            const model::TrajectorySignalCatalog* catalog,
                            QtPlaybackController* playback,
                            io::QtLoadResult* loaded,
                            QWidget* mainWindow,
                            ReaderMainWindow* readerWindow,
                            model::TransformedConformation* transformed) {
    ASSERT_THREAD(this);
    scene_ = scene;
    selection_ = selection;
    signalModel_ = signalModel;
    panelModel_ = panelModel;
    catalog_ = catalog;
    playback_ = playback;
    loaded_ = loaded;
    mainWindow_ = mainWindow;
    readerWindow_ = readerWindow;
    transformed_ = transformed;
    contextSet_ = true;
}

quint16 RestServer::listen(quint16 port) {
    ASSERT_THREAD(this);
    if (!contextSet_) {
        qCCritical(cRest).noquote() << "listen() called before setContext()";
        return 0;
    }

    server_ = std::make_unique<QHttpServer>(this);
    registerRoutes();

    const quint16 bound = server_->listen(QHostAddress::LocalHost, port);
    if (bound == 0) {
        qCCritical(cRest).noquote()
            << "REST server failed to bind 127.0.0.1 port" << port;
        server_.reset();
        return 0;
    }

    qCInfo(cRest).noquote() << "REST server listening on 127.0.0.1:" << bound;
    // Handshake line for the pytest fixture to scrape.
    std::fprintf(stderr, "H5READER_REST_PORT=%u\n", static_cast<unsigned>(bound));
    std::fflush(stderr);
    return bound;
}

void RestServer::registerRoutes() {
    using SC = QHttpServerResponse::StatusCode;
    using Method = QHttpServerRequest::Method;

    // ---- health ---------------------------------------------------------

    server_->route(QStringLiteral("/health"), [this]() {
        ASSERT_THREAD(this);
        return jsonResponse(QJsonObject{
            {"ok", true},
            {"version", QStringLiteral(H5READER_VERSION)},
        });
    });

    // ---- protein / atoms inventory --------------------------------------

    server_->route(QStringLiteral("/protein/atoms"), [this]() {
        ASSERT_THREAD(this);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        if (!protein)
            return errorResponse(QStringLiteral("no protein loaded"), SC::ServiceUnavailable);
        return jsonResponse(QJsonObject{
            {"count", static_cast<qint64>(protein->atomCount())},
        });
    });

    // ---- frame ----------------------------------------------------------

    server_->route(QStringLiteral("/frame/current"), [this]() {
        ASSERT_THREAD(this);
        const auto* conf = loaded_ ? loaded_->conformation.get() : nullptr;
        if (!playback_ || !conf)
            return errorResponse(QStringLiteral("playback not wired"), SC::ServiceUnavailable);
        const int frame = playback_->currentFrame();
        const double time_ps = (frame >= 0 && static_cast<std::size_t>(frame) < conf->frameCount())
                                   ? conf->timePicoseconds(static_cast<std::size_t>(frame))
                                   : 0.0;
        return jsonResponse(QJsonObject{
            {"frame", frame},
            {"time_ps", time_ps},
            {"count", playback_->frameCount()},
        });
    });

    server_->route(QStringLiteral("/frame/set"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!playback_)
            return errorResponse(QStringLiteral("playback not wired"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("frame"))
            return errorResponse(QStringLiteral("body must be {\"frame\": int}"), SC::BadRequest);
        const int frame = body.value("frame").toInt();
        playback_->pause();
        playback_->setFrame(frame);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- selection ------------------------------------------------------

    server_->route(QStringLiteral("/selection"), [this]() {
        ASSERT_THREAD(this);
        if (!selection_)
            return errorResponse(QStringLiteral("selection not wired"), SC::ServiceUnavailable);
        QJsonArray atoms;
        for (std::size_t a : selection_->atoms())
            atoms.append(static_cast<qint64>(a));
        QJsonObject out{
            {"atoms", atoms},
            {"count", static_cast<qint64>(selection_->count())},
        };
        out["focus"] = selection_->hasFocus()
                           ? QJsonValue(static_cast<qint64>(selection_->focus()))
                           : QJsonValue(QJsonValue::Null);
        return jsonResponse(out);
    });

    server_->route(QStringLiteral("/selection/pick"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!selection_)
            return errorResponse(QStringLiteral("selection not wired"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("atom"))
            return errorResponse(QStringLiteral("body must be {\"atom\": int, \"modifiers\": \"none\"|\"shift\"}"),
                                 SC::BadRequest);
        const auto atom = static_cast<std::size_t>(body.value("atom").toInteger());
        const QString modStr = body.value("modifiers").toString(QStringLiteral("none"));
        const Qt::KeyboardModifiers mods = (modStr == QStringLiteral("shift"))
                                               ? Qt::ShiftModifier
                                               : Qt::NoModifier;
        selection_->applyPick(atom, mods);
        return QHttpServerResponse(SC::NoContent);
    });

    server_->route(QStringLiteral("/selection/atoms"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!selection_)
            return errorResponse(QStringLiteral("selection not wired"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.value("atoms").isArray())
            return errorResponse(QStringLiteral("body must be {\"atoms\": [int, ...]}"),
                                 SC::BadRequest);
        const QJsonArray arr = body.value("atoms").toArray();
        if (arr.size() > static_cast<int>(model::AtomSelection::kMaxAtoms))
            return errorResponse(
                QStringLiteral("atoms array exceeds kMaxAtoms=%1")
                    .arg(model::AtomSelection::kMaxAtoms),
                SC::BadRequest);
        // Primary validation: bounds-check every index against the loaded
        // protein up front so a partial bulkSet doesn't silently drop
        // entries. Mirrors the existing /plane-lock/enable pattern.
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        if (!protein)
            return errorResponse(QStringLiteral("no protein loaded"), SC::ServiceUnavailable);
        std::vector<std::size_t> atoms;
        atoms.reserve(static_cast<std::size_t>(arr.size()));
        for (const QJsonValue& v : arr) {
            const qint64 raw = v.toInteger(-1);
            if (raw < 0 || static_cast<std::size_t>(raw) >= protein->atomCount())
                return errorResponse(
                    QStringLiteral("atom index out of range: %1 (atomCount=%2)")
                        .arg(raw).arg(protein->atomCount()),
                    SC::BadRequest);
            atoms.push_back(static_cast<std::size_t>(raw));
        }
        selection_->bulkSet(atoms);
        return QHttpServerResponse(SC::NoContent);
    });

    server_->route(QStringLiteral("/selection/clear"), Method::Post,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        if (!selection_)
            return errorResponse(QStringLiteral("selection not wired"), SC::ServiceUnavailable);
        selection_->clear();
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- selection / instrument preset (harness marker mode) -----------
    //
    // Switches the MeasurementOverlay's 4 sphere colours to a CPK-distinct
    // table (magenta / spring green / deep pink / vivid violet) at opacity
    // 1.0 and radius 1.5 Å. Designed so a Python harness can locate the
    // marker via connected-component blob analysis on a snapshot PNG: the
    // colours are outside every CPK element colour so a hue threshold isolates
    // the marker against any rendered scene. Reversible: {"enabled": false}
    // restores the Okabe-Ito palette + default opacity/radius.
    //
    // `focus_only` (default false, back-compat): when true AND enabled is
    // true, all four sphere actors get the magenta colour and only the
    // focus-slot sphere renders — eliminates the slot-1-eclipses-slot-0
    // problem the no-lock baseline run hit (VIEWPORT_OBSERVATIONS §5b).
    server_->route(QStringLiteral("/selection/instrument"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->measurementOverlay())
            return errorResponse(QStringLiteral("measurement overlay not wired"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("enabled") || !body.value("enabled").isBool())
            return errorResponse(QStringLiteral("body must be {\"enabled\": bool, \"focus_only\": bool}"),
                                 SC::BadRequest);
        const bool on = body.value("enabled").toBool();
        const bool focusOnly = (body.contains("focus_only") && body.value("focus_only").isBool())
                                   ? body.value("focus_only").toBool()
                                   : false;
        scene_->measurementOverlay()->setInstrumentMode(on, focusOnly);
        // The overlay does not Render itself (overlay contract,
        // MoleculeScene.h §1-5); we flush via the scene the same way the
        // ribbon/rings visibility toggles do at ReaderMainWindow.cpp:591.
        scene_->requestRender();
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- docks visibility (harness viewport-maximise) ------------------
    //
    // POST /docks/visible {"visible": bool} — hides or restores the three
    // ReaderMainWindow dock widgets (inspector, selection, dashboard strip)
    // wholesale. Hide stashes each dock's pre-hide visibility so restore
    // returns each dock to its prior state. Used by the harness to expand
    // the central VTK viewport so the marker blob fits in more pixels.
    server_->route(QStringLiteral("/docks/visible"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("visible") || !body.value("visible").isBool())
            return errorResponse(QStringLiteral("body must be {\"visible\": bool}"),
                                 SC::BadRequest);
        const bool visible = body.value("visible").toBool();
        readerWindow_->setDocksVisible(visible);
        // The viewport widget's geometry doesn't change without a paint
        // refresh; ask the scene to render so the new frame fills the
        // expanded central widget on the same tick the docks vanish.
        if (scene_) scene_->requestRender();
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- transform (upstream data-layer transform — TransformedConformation) -
    //
    // POST /transform {"kind": "identity"|"center_com"|"fit_reference"|"fit_subset",
    //                   "reference_frame": int (default 0),
    //                   "subset_atoms": [int, ...] (FitSubset only),
    //                   "backbone_only": bool (FitSubset shorthand) }
    // Switches the wrapped Conformation's transform mode. Fire-and-forget;
    // the wrapper emits transformChanged() which is connected (in
    // ReaderMainWindow) to scene_->refreshCurrentFrame so the molecule
    // re-renders in the new frame without further client involvement.
    //
    // GET /transform → returns the current mode + parameters for the
    // harness's reproducibility manifest.
    server_->route(QStringLiteral("/transform"), Method::Get,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        if (!transformed_)
            return errorResponse(QStringLiteral("transformed conformation not wired"),
                                 SC::ServiceUnavailable);
        QString kind;
        switch (transformed_->mode()) {
            case model::TransformedConformation::Mode::Identity:     kind = QStringLiteral("identity"); break;
            case model::TransformedConformation::Mode::CenterCom:    kind = QStringLiteral("center_com"); break;
            case model::TransformedConformation::Mode::FitReference: kind = QStringLiteral("fit_reference"); break;
            case model::TransformedConformation::Mode::FitSubset:    kind = QStringLiteral("fit_subset"); break;
        }
        QJsonArray subsetArr;
        for (std::size_t a : transformed_->subsetAtoms())
            subsetArr.append(static_cast<qint64>(a));
        return jsonResponse(QJsonObject{
            {"kind", kind},
            {"reference_frame", static_cast<qint64>(transformed_->referenceFrame())},
            {"subset_atoms", subsetArr},
            {"subset_size", subsetArr.size()},
        });
    });

    server_->route(QStringLiteral("/transform"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!transformed_)
            return errorResponse(QStringLiteral("transformed conformation not wired"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("kind"))
            return errorResponse(QStringLiteral("body must be {\"kind\": ..., ...}"),
                                 SC::BadRequest);
        const QString kindStr = body.value("kind").toString();
        const std::size_t referenceFrame = body.contains("reference_frame")
            ? static_cast<std::size_t>(body.value("reference_frame").toInteger(0))
            : 0;

        using Mode = model::TransformedConformation::Mode;
        Mode mode = Mode::Identity;
        std::vector<std::size_t> subset;

        if (kindStr == QStringLiteral("identity")) {
            mode = Mode::Identity;
        } else if (kindStr == QStringLiteral("center_com")) {
            mode = Mode::CenterCom;
        } else if (kindStr == QStringLiteral("fit_reference")) {
            mode = Mode::FitReference;
        } else if (kindStr == QStringLiteral("fit_subset")) {
            mode = Mode::FitSubset;
            // backbone_only shorthand: compute the subset from the typed
            // QtAtom::IsBackbone() flag (no string parsing of atom names;
            // chemistry identity comes from the typed substrate).
            const bool backboneOnly = body.contains("backbone_only") && body.value("backbone_only").toBool();
            if (backboneOnly) {
                const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
                if (!protein)
                    return errorResponse(QStringLiteral("no protein loaded for backbone_only"),
                                         SC::ServiceUnavailable);
                subset = model::TransformedConformation::BackboneSubset(*protein);
                if (subset.size() < 3)
                    return errorResponse(QStringLiteral("backbone subset has <3 atoms — "
                                                        "fit underdetermined"),
                                         SC::Conflict);
            } else if (body.value("subset_atoms").isArray()) {
                const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
                if (!protein)
                    return errorResponse(QStringLiteral("no protein loaded"),
                                         SC::ServiceUnavailable);
                const QJsonArray arr = body.value("subset_atoms").toArray();
                subset.reserve(static_cast<std::size_t>(arr.size()));
                for (const QJsonValue& v : arr) {
                    const qint64 raw = v.toInteger(-1);
                    if (raw < 0 || static_cast<std::size_t>(raw) >= protein->atomCount())
                        return errorResponse(QStringLiteral("subset_atoms index out of range: %1").arg(raw),
                                             SC::BadRequest);
                    subset.push_back(static_cast<std::size_t>(raw));
                }
                if (subset.size() < 3)
                    return errorResponse(QStringLiteral("subset must have >= 3 atoms (Kabsch underdetermined)"),
                                         SC::BadRequest);
            } else {
                return errorResponse(QStringLiteral("fit_subset needs subset_atoms or backbone_only"),
                                     SC::BadRequest);
            }
        } else {
            return errorResponse(QStringLiteral("unknown transform kind: %1").arg(kindStr),
                                 SC::BadRequest);
        }

        transformed_->setMode(mode, referenceFrame, std::move(subset));
        // setMode emits transformChanged → ReaderMainWindow connects this
        // to scene_->refreshCurrentFrame. No explicit render here; the
        // connected slot handles it.
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- plane lock -----------------------------------------------------

    server_->route(QStringLiteral("/plane-lock"), [this]() {
        ASSERT_THREAD(this);
        if (!scene_)
            return errorResponse(QStringLiteral("scene not wired"), SC::ServiceUnavailable);
        QJsonObject out;
        out["active"] = scene_->isCameraPlaneLocked();
        QJsonArray atoms;
        for (std::size_t a : scene_->cameraPlaneLockAtoms())
            atoms.append(static_cast<qint64>(a));
        out["atoms"] = atoms;
        return jsonResponse(out);
    });

    server_->route(QStringLiteral("/plane-lock/enable"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_)
            return errorResponse(QStringLiteral("scene not wired"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.value("atoms").isArray())
            return errorResponse(QStringLiteral("body must be {\"atoms\": [a, b, c]}"), SC::BadRequest);
        const QJsonArray arr = body.value("atoms").toArray();
        if (arr.size() != 3)
            return errorResponse(QStringLiteral("plane lock requires exactly three atom indices"), SC::BadRequest);
        std::vector<std::size_t> atoms;
        atoms.reserve(3);
        for (const QJsonValue& v : arr)
            atoms.push_back(static_cast<std::size_t>(v.toInteger()));
        if (!scene_->lockCameraToSelectionPlane(atoms))
            return errorResponse(QStringLiteral("scene rejected the lock (degenerate plane or invalid atoms)"),
                                 SC::Conflict);
        return QHttpServerResponse(SC::NoContent);
    });

    server_->route(QStringLiteral("/plane-lock/disable"), Method::Post,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        if (!scene_)
            return errorResponse(QStringLiteral("scene not wired"), SC::ServiceUnavailable);
        scene_->clearCameraPlaneLock();
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- camera mode (typed CameraMode + OrientationPolicy) ------------
    //
    // GET /camera/mode → {"mode": "free|atom|bond|dihedral|plane|subset",
    //                     "atoms": [...], "policy": "..."}
    //
    // POST /camera/mode {"mode": "...", "atoms": [...], "orientation":
    //                    {"kind": "default|free|perp_bond|down_axis|perp_plane",
    //                     "axis_atoms": [a,b]}}
    // POST /camera/clear  — equivalent to setMode(Free, Default).
    //
    // Per spec/viewport_pipeline_2026-05-30.md §I (REST surface). The
    // typed CameraMode replaces ad-hoc camera-lock endpoints; the
    // existing /plane-lock/* endpoints continue to work as shims.
    server_->route(QStringLiteral("/camera/mode"), Method::Get,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->cameraComposer())
            return errorResponse(QStringLiteral("camera composer not wired"), SC::ServiceUnavailable);
        const auto* composer = scene_->cameraComposer();
        QJsonArray atoms;
        for (std::size_t a : composer->mode().atoms)
            atoms.append(static_cast<qint64>(a));
        QJsonObject policy{
            {"kind", QString::fromLatin1(NameFor(composer->policy().kind))},
            {"axis_atoms", QJsonArray{
                static_cast<qint64>(composer->policy().axisAtoms[0]),
                static_cast<qint64>(composer->policy().axisAtoms[1]),
            }},
        };
        return jsonResponse(QJsonObject{
            {"mode", QString::fromLatin1(NameFor(composer->mode().kind))},
            {"atoms", atoms},
            {"orientation", policy},
        });
    });

    server_->route(QStringLiteral("/camera/mode"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->cameraComposer())
            return errorResponse(QStringLiteral("camera composer not wired"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("mode"))
            return errorResponse(
                QStringLiteral("body must be {\"mode\": str, ...}"),
                SC::BadRequest);
        const QString modeStr = body.value("mode").toString().toLower();
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        if (!protein)
            return errorResponse(QStringLiteral("no protein loaded"), SC::ServiceUnavailable);
        const std::size_t atomCount = protein->atomCount();

        auto readAtomIdx = [&](const QString& key) -> std::optional<std::size_t> {
            if (!body.contains(key)) return std::nullopt;
            const qint64 raw = body.value(key).toInteger(-1);
            if (raw < 0 || static_cast<std::size_t>(raw) >= atomCount)
                return std::nullopt;
            return static_cast<std::size_t>(raw);
        };
        auto readAtomArray = [&]() -> std::optional<std::vector<std::size_t>> {
            if (!body.value("atoms").isArray()) return std::nullopt;
            const QJsonArray arr = body.value("atoms").toArray();
            std::vector<std::size_t> atoms;
            atoms.reserve(static_cast<std::size_t>(arr.size()));
            for (const QJsonValue& v : arr) {
                const qint64 raw = v.toInteger(-1);
                if (raw < 0 || static_cast<std::size_t>(raw) >= atomCount)
                    return std::nullopt;
                atoms.push_back(static_cast<std::size_t>(raw));
            }
            return atoms;
        };

        // Round-trip note: GET /camera/mode returns {"mode","atoms",...}
        // with atoms as an array. POST accepts the same shape for every
        // mode that takes atoms, AND also accepts the human-friendly
        // per-mode keys (atom / a / b / c / d) as alternatives. Either
        // form is valid input; the array form lets a client GET → POST
        // round-trip without rewriting the payload.
        CameraMode mode;
        const auto atomsArray = readAtomArray();   // nullopt if "atoms" missing or invalid
        if (modeStr == QStringLiteral("free")) {
            mode = FreeMode();
        } else if (modeStr == QStringLiteral("atom")) {
            std::optional<std::size_t> a;
            if (atomsArray && atomsArray->size() == 1) {
                a = atomsArray->at(0);
            } else {
                a = readAtomIdx(QStringLiteral("atom"));
            }
            if (!a)
                return errorResponse(QStringLiteral(
                    "atom mode needs {\"atom\": int} or {\"atoms\": [int]}"), SC::BadRequest);
            mode = AtomMode(*a);
        } else if (modeStr == QStringLiteral("bond")) {
            std::optional<std::size_t> a, b;
            if (atomsArray && atomsArray->size() == 2) {
                a = atomsArray->at(0);
                b = atomsArray->at(1);
            } else {
                a = readAtomIdx(QStringLiteral("a"));
                b = readAtomIdx(QStringLiteral("b"));
            }
            if (!a || !b)
                return errorResponse(QStringLiteral(
                    "bond mode needs {\"a\": int, \"b\": int} or {\"atoms\": [a, b]}"),
                    SC::BadRequest);
            mode = BondMode(*a, *b);
        } else if (modeStr == QStringLiteral("dihedral")) {
            std::optional<std::size_t> a, b, c, d;
            if (atomsArray && atomsArray->size() == 4) {
                a = atomsArray->at(0);
                b = atomsArray->at(1);
                c = atomsArray->at(2);
                d = atomsArray->at(3);
            } else {
                a = readAtomIdx(QStringLiteral("a"));
                b = readAtomIdx(QStringLiteral("b"));
                c = readAtomIdx(QStringLiteral("c"));
                d = readAtomIdx(QStringLiteral("d"));
            }
            if (!a || !b || !c || !d)
                return errorResponse(QStringLiteral(
                    "dihedral needs {\"a\"..\"d\": int} or {\"atoms\": [a, b, c, d]}"),
                    SC::BadRequest);
            mode = DihedralMode(*a, *b, *c, *d);
        } else if (modeStr == QStringLiteral("plane")) {
            std::optional<std::size_t> a, b, c;
            if (atomsArray && atomsArray->size() == 3) {
                a = atomsArray->at(0);
                b = atomsArray->at(1);
                c = atomsArray->at(2);
            } else {
                a = readAtomIdx(QStringLiteral("a"));
                b = readAtomIdx(QStringLiteral("b"));
                c = readAtomIdx(QStringLiteral("c"));
            }
            if (!a || !b || !c)
                return errorResponse(QStringLiteral(
                    "plane needs {\"a\", \"b\", \"c\": int} or {\"atoms\": [a, b, c]}"),
                    SC::BadRequest);
            mode = PlaneMode(*a, *b, *c);
        } else if (modeStr == QStringLiteral("subset")) {
            // Backbone-only shortcut mirrors the /transform endpoint.
            const bool backboneOnly = body.contains("backbone_only")
                                      && body.value("backbone_only").toBool();
            std::vector<std::size_t> atoms;
            if (backboneOnly) {
                atoms = model::TransformedConformation::BackboneSubset(*protein);
                if (atoms.size() < 3)
                    return errorResponse(QStringLiteral("backbone subset has <3 atoms"),
                                         SC::Conflict);
            } else {
                if (!atomsArray)
                    return errorResponse(QStringLiteral("subset needs {\"atoms\": [...]} or backbone_only=true"),
                                         SC::BadRequest);
                atoms = *atomsArray;
                if (atoms.size() < 3)
                    return errorResponse(QStringLiteral("subset needs >=3 atoms"),
                                         SC::BadRequest);
            }
            mode = SubsetMode(std::move(atoms));
        } else {
            return errorResponse(QStringLiteral("unknown mode: %1").arg(modeStr), SC::BadRequest);
        }

        OrientationPolicy policy = DefaultPolicy();
        if (body.contains("orientation") && body.value("orientation").isObject()) {
            const QJsonObject po = body.value("orientation").toObject();
            const QString kindStr = po.value("kind").toString().toLower();
            if (kindStr == QStringLiteral("default")) {
                policy = DefaultPolicy();
            } else if (kindStr == QStringLiteral("free")) {
                policy = FreePolicy();
            } else if (kindStr == QStringLiteral("perp_bond")
                       || kindStr == QStringLiteral("perpendiculartobond")) {
                policy = PerpToBondPolicy();
            } else if (kindStr == QStringLiteral("perp_plane")
                       || kindStr == QStringLiteral("perpendiculartoplane")) {
                policy = PerpToPlanePolicy();
            } else if (kindStr == QStringLiteral("down_axis")
                       || kindStr == QStringLiteral("downaxis")) {
                if (!po.value("axis_atoms").isArray() || po.value("axis_atoms").toArray().size() != 2)
                    return errorResponse(QStringLiteral("down_axis needs axis_atoms: [a, b]"), SC::BadRequest);
                const QJsonArray arr = po.value("axis_atoms").toArray();
                const qint64 a0 = arr.at(0).toInteger(-1);
                const qint64 a1 = arr.at(1).toInteger(-1);
                if (a0 < 0 || a1 < 0 || static_cast<std::size_t>(a0) >= atomCount
                    || static_cast<std::size_t>(a1) >= atomCount)
                    return errorResponse(QStringLiteral("down_axis axis_atoms out of range"), SC::BadRequest);
                policy = DownAxisPolicy(static_cast<std::size_t>(a0), static_cast<std::size_t>(a1));
            } else {
                return errorResponse(QStringLiteral("unknown orientation kind: %1").arg(kindStr), SC::BadRequest);
            }
        }

        const std::size_t currentFrame = playback_ ? static_cast<std::size_t>(playback_->currentFrame()) : 0;
        scene_->cameraComposer()->setMode(std::move(mode), policy, currentFrame);
        // One render to surface the new camera state on this tick.
        scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return QHttpServerResponse(SC::NoContent);
    });

    server_->route(QStringLiteral("/camera/clear"), Method::Post,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->cameraComposer())
            return errorResponse(QStringLiteral("camera composer not wired"), SC::ServiceUnavailable);
        const std::size_t currentFrame = playback_ ? static_cast<std::size_t>(playback_->currentFrame()) : 0;
        scene_->cameraComposer()->setMode(FreeMode(), DefaultPolicy(), currentFrame);
        scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- camera focus atom convenience -----------------------------------
    //
    // POST /camera/focus_atom {"atom": int, "kind": "plane"|"dihedral_phi"
    //                          |"dihedral_psi"}
    //
    // Derives a typed CameraMode + OrientationPolicy from a focus atom by
    // reaching into its residue's backbone-atom-index cache (built at
    // load time, no string scan). Plane = N/CA/C plane lock (the
    // canonical "focus atom + local neighborhood coherent" recipe per
    // HARNESS_BASELINE_PIPELINE doc); dihedral_phi/psi sight down the
    // residue's phi or psi torsion (Newman-projection view).
    //
    // Failure shapes:
    //   atom out of range            -> 400
    //   atom has no residue          -> 422 (data shape issue)
    //   residue missing N / CA / C   -> 422 (e.g. unusual terminal)
    //   dihedral missing flanking    -> 422 (terminal residue)
    server_->route(QStringLiteral("/camera/focus_atom"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->cameraComposer())
            return errorResponse(QStringLiteral("camera composer not wired"),
                                 SC::ServiceUnavailable);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        if (!protein)
            return errorResponse(QStringLiteral("no protein loaded"),
                                 SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("atom") || !body.contains("kind"))
            return errorResponse(QStringLiteral(
                "body must be {\"atom\": int, \"kind\": str}"), SC::BadRequest);

        const qint64 rawAtom = body.value("atom").toInteger(-1);
        if (rawAtom < 0 || static_cast<std::size_t>(rawAtom) >= protein->atomCount())
            return errorResponse(QStringLiteral("atom out of range"), SC::BadRequest);
        const std::size_t atomIdx = static_cast<std::size_t>(rawAtom);

        const QString kindStr = body.value("kind").toString().toLower();
        FocusAnchorKind kind;
        if (kindStr == QStringLiteral("plane")) {
            kind = FocusAnchorKind::Plane;
        } else if (kindStr == QStringLiteral("dihedral_phi")
                   || kindStr == QStringLiteral("dihedral")) {
            // "dihedral" alias defaults to phi — the more common torsion
            // for inspecting local backbone behaviour.
            kind = FocusAnchorKind::DihedralPhi;
        } else if (kindStr == QStringLiteral("dihedral_psi")) {
            kind = FocusAnchorKind::DihedralPsi;
        } else {
            return errorResponse(QStringLiteral(
                "kind must be plane | dihedral | dihedral_phi | dihedral_psi"),
                SC::BadRequest);
        }

        const auto result = DeriveFocusAnchor(*protein, atomIdx, kind);
        switch (result.outcome) {
            case FocusAnchorOutcome::Ok:
                break;
            case FocusAnchorOutcome::AtomIndexOutOfRange:
                return errorResponse(QStringLiteral("atom out of range"),
                                     SC::BadRequest);
            case FocusAnchorOutcome::AtomHasNoResidue:
                return errorResponse(QStringLiteral(
                    "focus atom has no residue (residue_index < 0)"),
                    SC::UnprocessableEntity);
            case FocusAnchorOutcome::MissingBackboneAtoms:
                return errorResponse(QStringLiteral(
                    "focus atom's residue lacks N / CA / C backbone atoms"),
                    SC::UnprocessableEntity);
            case FocusAnchorOutcome::MissingDihedralNeighbor:
                return errorResponse(QStringLiteral(
                    "focus atom's residue has no flanking residue for the "
                    "requested dihedral (terminal residue)"),
                    SC::UnprocessableEntity);
        }

        const std::size_t currentFrame = playback_
            ? static_cast<std::size_t>(playback_->currentFrame()) : 0;
        scene_->cameraComposer()->setMode(std::move(result.mode),
                                            result.policy,
                                            currentFrame);
        scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- log mask (bitmask gate for StructuredLogger) ------------------
    //
    // GET /log/mask → {"mask": int, "categories": ["FRAME", "CAMERA", ...]}
    // POST /log/mask {"mask": int} OR {"categories": [...]}
    //
    // Per spec/viewport_pipeline_2026-05-30.md §H + implementation prompt
    // §3 (bitmask logging instead of UDP throttling). RENDER (0x01) is
    // off by default; flip it on when debugging the render scheduler.
    server_->route(QStringLiteral("/log/mask"), Method::Get,
                   [](const QHttpServerRequest&) {
        const std::uint32_t mask = diagnostics::StructuredLogger::CategoryMask();
        QJsonArray cats;
        for (const QString& n : diagnostics::StructuredLogger::SymbolicNamesFromMask(mask))
            cats.append(n);
        return jsonResponse(QJsonObject{
            {"mask", static_cast<qint64>(mask)},
            {"categories", cats},
        });
    });

    server_->route(QStringLiteral("/log/mask"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(
                QStringLiteral("body must be {\"mask\": int} or {\"categories\": [...]}"),
                SC::BadRequest);
        std::uint32_t mask = diagnostics::StructuredLogger::CategoryMask();
        if (body.contains("mask")) {
            const qint64 raw = body.value("mask").toInteger(-1);
            if (raw < 0)
                return errorResponse(QStringLiteral("mask must be a non-negative integer"),
                                     SC::BadRequest);
            mask = static_cast<std::uint32_t>(raw);
        } else if (body.value("categories").isArray()) {
            QStringList names;
            for (const QJsonValue& v : body.value("categories").toArray())
                names.append(v.toString());
            mask = diagnostics::StructuredLogger::MaskFromSymbolicNames(names);
        } else {
            return errorResponse(
                QStringLiteral("body must be {\"mask\": int} or {\"categories\": [...]}"),
                SC::BadRequest);
        }
        diagnostics::StructuredLogger::SetCategoryMask(mask);
        QJsonArray cats;
        for (const QString& n : diagnostics::StructuredLogger::SymbolicNamesFromMask(mask))
            cats.append(n);
        return jsonResponse(QJsonObject{
            {"mask", static_cast<qint64>(mask)},
            {"categories", cats},
        });
    });

    // ---- atom positions (per-frame, for tests that need plane math) -----
    //
    // Reads through the TransformedConformation wrapper so the positions
    // returned are in the active display frame (matches what the renderer
    // shows). Tests that need raw H5 positions would need a different
    // endpoint (none today; add when a use case appears).
    server_->route(QStringLiteral("/positions"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        const model::Conformation* conf = transformed_
            ? static_cast<const model::Conformation*>(transformed_.data())
            : (loaded_ ? loaded_->conformation.get() : nullptr);
        if (!protein || !conf)
            return errorResponse(QStringLiteral("protein/conformation not wired"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("atoms") || !body.contains("frame"))
            return errorResponse(QStringLiteral("body must be {\"atoms\": [...], \"frame\": int}"),
                                 SC::BadRequest);
        const QJsonArray atomsArr = body.value("atoms").toArray();
        const int frame = body.value("frame").toInt();
        if (frame < 0 || static_cast<std::size_t>(frame) >= conf->frameCount())
            return errorResponse(QStringLiteral("frame out of range"), SC::BadRequest);
        QJsonArray out;
        for (const QJsonValue& v : atomsArr) {
            const auto atom = static_cast<std::size_t>(v.toInteger());
            if (atom >= protein->atomCount())
                return errorResponse(QStringLiteral("atom out of range"), SC::BadRequest);
            const model::Vec3 p = conf->atomPosition(static_cast<std::size_t>(frame), atom);
            out.append(QJsonObject{
                {"atom", static_cast<qint64>(atom)},
                {"position", QJsonArray{p.x(), p.y(), p.z()}},
            });
        }
        return jsonResponse(QJsonObject{
            {"frame", frame},
            {"positions", out},
        });
    });

    // ---- scene camera readback ------------------------------------------

    server_->route(QStringLiteral("/scene/camera"), [this]() {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->Renderer())
            return errorResponse(QStringLiteral("scene not wired"), SC::ServiceUnavailable);
        auto* camera = scene_->Renderer()->GetActiveCamera();
        if (!camera)
            return errorResponse(QStringLiteral("no active camera"), SC::ServiceUnavailable);
        double focal[3]{}, position[3]{}, viewUp[3]{}, direction[3]{};
        camera->GetFocalPoint(focal);
        camera->GetPosition(position);
        camera->GetViewUp(viewUp);
        camera->GetDirectionOfProjection(direction);
        return jsonResponse(QJsonObject{
            {"focal", vec3FromRaw(focal)},
            {"position", vec3FromRaw(position)},
            {"view_up", vec3FromRaw(viewUp)},
            {"direction", vec3FromRaw(direction)},
        });
    });

    // ---- dashboard signals (read-only listing) --------------------------

    server_->route(QStringLiteral("/dashboard/signals"), [this]() {
        ASSERT_THREAD(this);
        if (!signalModel_)
            return errorResponse(QStringLiteral("signal model not wired"), SC::ServiceUnavailable);
        QJsonArray out;
        for (const model::DashboardSignal& s : signalModel_->activeSignals()) {
            QJsonArray modes;
            for (const QString& m : s.displayModeIds) modes.append(m);
            out.append(QJsonObject{
                {"id", s.id.toString(QUuid::WithoutBraces)},
                {"label", s.label},
                {"descriptor_id", s.binding.descriptorId},
                {"concept_key", s.binding.conceptKey},
                {"display_mode_id", s.binding.displayModeId},
                {"display_modes", modes},
                {"anchor", anchorToJson(s.binding.anchor)},
                {"enabled", s.enabled},
                {"follows_focus", s.binding.followsFocus},
            });
        }
        return jsonResponse(out);
    });

    // ---- shutdown -------------------------------------------------------

    // POST /shutdown — graceful exit so the operator (or a test
    // harness) doesn't need pkill. Returns 204 immediately, then
    // posts a deferred quit() to the event loop so the response can
    // flush first. Async-fire vs synchronous quit() matters because
    // QHttpServer holds the request handler frame; calling quit()
    // synchronously would tear down before the response bytes ship.
    server_->route(QStringLiteral("/shutdown"), Method::Post,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        qCInfo(cRest).noquote() << "REST /shutdown — quitting on next event loop tick";
        QTimer::singleShot(0, qApp, &QCoreApplication::quit);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- screenshot -----------------------------------------------------

    server_->route(QStringLiteral("/screenshot"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        const QString target = (ok && body.contains("target"))
                                   ? body.value("target").toString()
                                   : QStringLiteral("scene");
        // force_render: default true (back-compat). false skips the
        // vtkWindowToImageFilter::ShouldRerender step so the snapshot reads
        // whatever pixels are currently in the framebuffer — the harness
        // mode for the paint-cycle-inversion experiment (VIEWPORT
        // OBSERVATIONS §5b). Only meaningful for target="scene".
        const bool forceRender = (ok && body.contains("force_render") && body.value("force_render").isBool())
                                     ? body.value("force_render").toBool()
                                     : true;
        QByteArray png;
        if (target == QStringLiteral("window")) {
            png = captureWindowPng(mainWindow_.data());
        } else if (target == QStringLiteral("scene")) {
            png = captureScenePng(scene_.data(), forceRender);
        } else {
            return errorResponse(QStringLiteral("target must be \"scene\" or \"window\""), SC::BadRequest);
        }
        if (png.isEmpty())
            return errorResponse(QStringLiteral("screenshot capture failed"), SC::InternalServerError);
        return QHttpServerResponse(QByteArray(kMimePng), png);
    });
}

}  // namespace h5reader::app
