#include "RestServer.h"
#include "MainWindow.h"

// Library headers — REST server reads the protein model directly
#include "AminoAcidType.h"
#include "Atom.h"
#include "Bond.h"
#include "ConformationAtom.h"
#include "ConformationResult.h"
#include "DsspResult.h"
#include "FieldGridOverlay.h"
#include "JobSpec.h"
#include "LegacyAmberTopology.h"
#include "MopacResult.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "Ring.h"
#include "SemanticEnums.h"
#include "Types.h"

#include <filesystem>

#include <QJsonDocument>
#include <QJsonArray>
#include <QJsonValue>
#include <QApplication>
#include <QComboBox>
#include <QCheckBox>
#include <QSlider>
#include <QDoubleSpinBox>
#include <QPlainTextEdit>
#include <QThread>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkCamera.h>
#include <iostream>
#include <cmath>

using namespace nmr;

RestServer::RestServer(MainWindow* mainWindow, quint16 port, QObject* parent)
    : QObject(parent)
    , server_(new QTcpServer(this))
    , mainWindow_(mainWindow)
{
    connect(server_, &QTcpServer::newConnection,
            this, &RestServer::onNewConnection);

    // Try requested port, then up to 9 more (9147, 9148, ..., 9156)
    for (quint16 p = port; p < port + 10; ++p) {
        if (server_->listen(QHostAddress::LocalHost, p)) {
            actualPort_ = p;
            std::cout << "REST server listening on localhost:" << p;
            if (p != port) std::cout << " (requested " << port << " was in use)";
            std::cout << "\n";
            return;
        }
    }
    std::cerr << "REST server: could not bind ports " << port
              << "-" << (port + 9) << "\n";
}

void RestServer::onNewConnection() {
    while (server_->hasPendingConnections()) {
        auto* socket = server_->nextPendingConnection();
        clients_.append(socket);
        connect(socket, &QTcpSocket::readyRead,
                this, &RestServer::onReadyRead);
        connect(socket, &QTcpSocket::disconnected,
                this, &RestServer::onDisconnected);
    }
}

void RestServer::onReadyRead() {
    auto* socket = qobject_cast<QTcpSocket*>(sender());
    if (!socket) return;

    while (socket->canReadLine()) {
        QByteArray const line = socket->readLine().trimmed();
        if (line.isEmpty()) continue;

        QJsonParseError err;
        QJsonDocument const doc = QJsonDocument::fromJson(line, &err);
        if (err.error != QJsonParseError::NoError || !doc.isObject()) {
            QJsonObject resp;
            resp["ok"] = false;
            resp["error"] = QString("JSON parse error: %1").arg(err.errorString());
            sendResponse(socket, resp);
            continue;
        }

        QJsonObject const result = dispatch(doc.object());
        sendResponse(socket, result);

        // cmdQuit sets shutdown_after_reply_. After sendResponse queues
        // the reply bytes on the socket, kick disconnectFromHost which
        // Qt documents as "If there is pending data waiting to be
        // written, QTcpSocket will enter ClosingState and wait until
        // all data has been written. Eventually, it will enter
        // UnconnectedState and emit the disconnected() signal."
        //
        // The existing onDisconnected() slot (connected at socket
        // accept) sees the shutdown flag and queues QCoreApplication
        // ::quit via QueuedConnection — by then the bytes are out, the
        // socket is closed, and quit can fire cleanly without racing
        // the write. No new connection needed; no timer; no sync wait.
        if (shutdown_after_reply_) {
            socket->disconnectFromHost();
            return;
        }
    }
}

void RestServer::onDisconnected() {
    auto* socket = qobject_cast<QTcpSocket*>(sender());
    if (socket) {
        clients_.removeAll(socket);
        socket->deleteLater();
    }
    // If cmdQuit set the shutdown flag, the disconnect we just
    // processed was triggered by RestServer::onReadyRead's
    // socket->disconnectFromHost() AFTER the reply bytes were queued.
    // The write buffer has drained (per Qt's ClosingState →
    // UnconnectedState contract) so the reply is on the wire and quit
    // can fire cleanly. Queued invocation lets the current slot finish
    // (including the deleteLater above) before quit walks the event
    // loop.
    if (shutdown_after_reply_) {
        QMetaObject::invokeMethod(qApp, &QCoreApplication::quit, Qt::QueuedConnection);
    }
}

void RestServer::sendResponse(QTcpSocket* socket, const QJsonObject& response) {
    QJsonDocument const doc(response);
    socket->write(doc.toJson(QJsonDocument::Compact));
    socket->write("\n");
    socket->flush();
}

QJsonObject RestServer::dispatch(const QJsonObject& cmd) {
    QString const action = cmd["cmd"].toString();

    if (action == "status")           return cmdStatus();
    if (action == "load_pdb")         return cmdLoadPdb(cmd);
    if (action == "load_protein_dir") return cmdLoadProteinDir(cmd);
    if (action == "set_overlay")      return cmdSetOverlay(cmd);
    if (action == "set_render_mode")  return cmdSetRenderMode(cmd);
    if (action == "screenshot")       return cmdScreenshot(cmd);
    if (action == "orbit")            return cmdOrbit(cmd);
    if (action == "reset_view")       return cmdResetView(cmd);
    if (action == "set_calculators")  return cmdSetCalculators(cmd);
    if (action == "show_rings")       return cmdShowRings(cmd);
    if (action == "show_bonds")       return cmdShowBonds(cmd);
    if (action == "show_butterfly")   return cmdShowButterfly(cmd);
    if (action == "set_glyph_scale")  return cmdSetGlyphScale(cmd);
    if (action == "set_opacity")      return cmdSetOpacity(cmd);
    if (action == "set_iso_threshold") return cmdSetIsoThreshold(cmd);
    if (action == "show_field_grid")  return cmdShowFieldGrid(cmd);
    if (action == "get_camera")       return cmdGetCamera(cmd);
    if (action == "set_camera")       return cmdSetCamera(cmd);
    if (action == "look_at_ring")     return cmdLookAtRing(cmd);
    if (action == "look_at_atom")     return cmdLookAtAtom(cmd);
    if (action == "list_rings")       return cmdListRings(cmd);
    if (action == "get_log")          return cmdGetLog(cmd);
    if (action == "export_features")  return cmdExportFeatures(cmd);
    if (action == "atom_dump") {
        return cmdAtomDump(cmd);
    }
    if (action == "list_atoms") {
        return cmdListAtoms(cmd);
    }
    if (action == "quit") {
        return cmdQuit();
    }

    QJsonObject resp;
    resp["ok"] = false;
    resp["error"] = QString("Unknown command: %1").arg(action);
    return resp;
}

// ---- Command implementations ----

// Polite shutdown — flip the flag, return the reply. onReadyRead()
// observes the flag after sendResponse() has queued the bytes on the
// socket's write buffer and drives the disconnected → quit signal
// chain (no timers, no sync waits). main_viewer.cpp's `aboutToQuit`
// then triggers MainWindow::shutdown to stop timers, workers, VTK.
QJsonObject RestServer::cmdQuit() {
    shutdown_after_reply_ = true;
    QJsonObject resp;
    resp["ok"] = true;
    QJsonObject result;
    result["message"] = "shutting down";
    resp["result"] = result;
    return resp;
}

QJsonObject RestServer::cmdStatus() {
    QJsonObject resp;
    resp["ok"] = true;
    QJsonObject result;

    // Report computation state
    result["computing"] = (mainWindow_->workerThread_ != nullptr &&
                           mainWindow_->workerThread_->isRunning());

    // Read directly from the library protein model. HeuristicTier counts
    // (n_report / n_pass / n_silent) were removed from this payload —
    // pre-kernel-catalogue prediction framing per UI_ROADMAP Known Issues
    // #1; no current ConformationResult writes those fields. Removing
    // before downstream REST clients pin a contract on dead data.
    auto& protein = mainWindow_->protein_;
    if (protein) {
        result["protein"] = QString::fromStdString(mainWindow_->currentProteinId_);
        result["n_atoms"] = static_cast<int>(protein->AtomCount());
        result["n_bonds"] = static_cast<int>(protein->BondCount());
        result["n_rings"] = static_cast<int>(protein->RingCount());
        result["n_residues"] = static_cast<int>(protein->ResidueCount());
    } else {
        result["protein"] = "";
        result["n_atoms"] = 0;
        result["n_bonds"] = 0;
        result["n_rings"] = 0;
        result["n_residues"] = 0;
    }

    result["overlay_mode"] = 0;  // overlay modes removed
    result["n_field_grids"] = static_cast<int>(mainWindow_->fieldGrids_.size());
    result["n_butterfly_grids"] = static_cast<int>(mainWindow_->butterflyFields_.size());
    result["field_grid_overlay"] = (mainWindow_->fieldGridOverlay_ != nullptr);
    result["butterfly_overlay"] = (mainWindow_->butterflyOverlay_ != nullptr);
    if (!mainWindow_->fieldGrids_.empty()) {
        double gmin = 1e30;
        double gmax = -1e30;
        int nz = 0;
        for (const auto& g : mainWindow_->fieldGrids_) {
            for (double const v : g.T0) {
                if (v < gmin) gmin = v;
                if (v > gmax) gmax = v;
                if (std::abs(v) > 1e-10) nz++;
            }
        }
        result["grid_T0_min"] = gmin;
        result["grid_T0_max"] = gmax;
        result["grid_T0_nonzero"] = nz;
    }
    resp["result"] = result;
    return resp;
}

QJsonObject RestServer::cmdLoadPdb(const QJsonObject& cmd) {
    QString const path = cmd["path"].toString();
    if (path.isEmpty()) {
        return QJsonObject{{"ok", false}, {"error", "missing 'path'"}};
    }
    nmr::JobSpec spec;
    spec.mode = nmr::JobMode::Pdb;
    spec.pdb_path = path.toStdString();
    mainWindow_->loadFromJobSpec(spec);
    return QJsonObject{{"ok", true}, {"result", QJsonObject{{"status", "loading"}}}};
}

QJsonObject RestServer::cmdLoadProteinDir(const QJsonObject& cmd) {
    // Deprecated: use --orca --root or --mutant --wt/--ala instead.
    // For REST backwards compatibility, treat as a PDB load of the first .pdb found.
    QString const path = cmd["path"].toString();
    if (path.isEmpty()) {
        return QJsonObject{{"ok", false}, {"error", "missing 'path'"}};
    }
    nmr::JobSpec spec;
    spec.mode = nmr::JobMode::Pdb;
    spec.pdb_path = path.toStdString();
    mainWindow_->loadFromJobSpec(spec);
    return QJsonObject{{"ok", true}, {"result", QJsonObject{{"status", "loading"},
        {"note", "loadProteinDir deprecated — use --orca --root or --mutant from CLI"}}}};
}

QJsonObject RestServer::cmdSetOverlay(const QJsonObject& cmd) {
    QString const mode = cmd["mode"].toString();
    // Overlay modes removed — per-calculator visualizations replace them
    return QJsonObject{{"ok", true}, {"result", QJsonObject{{"mode", "none"},
        {"note", "overlay modes removed; use per-calculator toggles"}}}};
}

QJsonObject RestServer::cmdSetRenderMode(const QJsonObject& cmd) {
    QString const mode = cmd["mode"].toString();
    static const QMap<QString, int> modes = {
        {"ball_stick", 0}, {"stick", 1}, {"liquorice", 1}
    };
    auto it = modes.find(mode);
    if (it == modes.end()) {
        return QJsonObject{{"ok", false}, {"error", "unknown render mode"}};
    }
    mainWindow_->renderModeCombo_->setCurrentIndex(it.value());
    mainWindow_->onRenderModeChanged(it.value());
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdScreenshot(const QJsonObject& cmd) {
    QString const path = cmd["path"].toString();
    if (path.isEmpty()) {
        return QJsonObject{{"ok", false}, {"error", "missing 'path'"}};
    }

    auto filter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    filter->SetInput(mainWindow_->renderWindow_);
    filter->SetScale(1);
    filter->SetInputBufferTypeToRGBA();
    filter->ReadFrontBufferOn();
    filter->ShouldRerenderOff();
    filter->Update();

    auto writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(path.toStdString().c_str());
    writer->SetInputConnection(filter->GetOutputPort());
    writer->Write();

    int* sz = mainWindow_->renderWindow_->GetSize();
    QJsonObject result;
    result["path"] = path;
    result["width"] = sz[0];
    result["height"] = sz[1];
    return QJsonObject{{"ok", true}, {"result", result}};
}

QJsonObject RestServer::cmdOrbit(const QJsonObject& cmd) {
    double const azimuth = cmd.value("azimuth").toDouble(0);
    double const elevation = cmd.value("elevation").toDouble(0);

    auto* camera = mainWindow_->renderer_->GetActiveCamera();
    camera->Azimuth(azimuth);
    camera->Elevation(elevation);
    camera->OrthogonalizeViewUp();
    mainWindow_->renderWindow_->Render();

    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdResetView(const QJsonObject& cmd) {
    Q_UNUSED(cmd);
    mainWindow_->renderer_->ResetCamera();
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetCalculators(const QJsonObject& cmd) {
    Q_UNUSED(cmd);
    // Per-calculator toggles will replace this when visualizations are added
    return QJsonObject{{"ok", true}, {"note", "calculator toggles pending per-calculator viz"}};
}

QJsonObject RestServer::cmdShowRings(const QJsonObject& cmd) {
    bool const visible = cmd["visible"].toBool(true);
    mainWindow_->showRingsCheck_->setChecked(visible);
    mainWindow_->onShowRingsToggled(visible);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdShowBonds(const QJsonObject& cmd) {
    bool const visible = cmd["visible"].toBool(true);
    mainWindow_->showPeptideBondsCheck_->setChecked(visible);
    mainWindow_->onShowPeptideBondsToggled(visible);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdShowButterfly(const QJsonObject& cmd) {
    bool const visible = cmd["visible"].toBool(true);
    mainWindow_->showButterflyCheck_->setChecked(visible);
    mainWindow_->onShowButterflyToggled(visible);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetGlyphScale(const QJsonObject& cmd) {
    double const scale = cmd["scale"].toDouble(0.5);
    int const sliderVal = static_cast<int>(scale * 100);
    mainWindow_->glyphScaleSlider_->setValue(sliderVal);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetOpacity(const QJsonObject& cmd) {
    double const opacity = cmd["value"].toDouble(0.7);
    int const sliderVal = static_cast<int>(opacity * 100);
    mainWindow_->opacitySlider_->setValue(sliderVal);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetIsoThreshold(const QJsonObject& cmd) {
    double const threshold = cmd["value"].toDouble(0.1);
    int const sliderVal = static_cast<int>(threshold * 100);
    mainWindow_->isoThresholdSlider_->setValue(sliderVal);
    mainWindow_->onIsoThresholdChanged();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdShowFieldGrid(const QJsonObject& cmd) {
    if (cmd.contains("shielded")) {
        bool const vis = cmd["shielded"].toBool(true);
        mainWindow_->showFieldGridCheck_->setChecked(vis);
        if (mainWindow_->fieldGridOverlay_) {
            mainWindow_->fieldGridOverlay_->setShieldedVisible(vis);
        }
    }
    if (cmd.contains("deshielded")) {
        bool const vis = cmd["deshielded"].toBool(true);
        mainWindow_->showDeshieldedCheck_->setChecked(vis);
        if (mainWindow_->fieldGridOverlay_) {
            mainWindow_->fieldGridOverlay_->setDeshieldedVisible(vis);
        }
    }
    if (!cmd.contains("shielded") && !cmd.contains("deshielded")) {
        bool const vis = cmd["visible"].toBool(true);
        mainWindow_->showFieldGridCheck_->setChecked(vis);
        mainWindow_->showDeshieldedCheck_->setChecked(vis);
        if (mainWindow_->fieldGridOverlay_) {
            mainWindow_->fieldGridOverlay_->setVisible(vis);
        }
    }
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdGetCamera(const QJsonObject&) {
    auto* camera = mainWindow_->renderer_->GetActiveCamera();
    double pos[3];
    double foc[3];
    double up[3];
    camera->GetPosition(pos);
    camera->GetFocalPoint(foc);
    camera->GetViewUp(up);

    QJsonObject result;
    result["position"] = QJsonArray{pos[0], pos[1], pos[2]};
    result["focal_point"] = QJsonArray{foc[0], foc[1], foc[2]};
    result["view_up"] = QJsonArray{up[0], up[1], up[2]};
    result["distance"] = camera->GetDistance();
    result["view_angle"] = camera->GetViewAngle();
    return QJsonObject{{"ok", true}, {"result", result}};
}

QJsonObject RestServer::cmdSetCamera(const QJsonObject& cmd) {
    auto* camera = mainWindow_->renderer_->GetActiveCamera();

    if (cmd.contains("position")) {
        QJsonArray p = cmd["position"].toArray();
        camera->SetPosition(p[0].toDouble(), p[1].toDouble(), p[2].toDouble());
    }
    if (cmd.contains("focal_point")) {
        QJsonArray f = cmd["focal_point"].toArray();
        camera->SetFocalPoint(f[0].toDouble(), f[1].toDouble(), f[2].toDouble());
    }
    if (cmd.contains("view_up")) {
        QJsonArray u = cmd["view_up"].toArray();
        camera->SetViewUp(u[0].toDouble(), u[1].toDouble(), u[2].toDouble());
    }
    camera->OrthogonalizeViewUp();
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdLookAtRing(const QJsonObject& cmd) {
    auto& protein = mainWindow_->protein_;
    if (!protein) {
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};
    }

    int const ringIdx = cmd["ring"].toInt(-1);
    if (ringIdx < 0 || ringIdx >= static_cast<int>(protein->RingCount())) {
        return QJsonObject{{"ok", false}, {"error", "invalid ring index"}};
    }

    const auto& conf = protein->Conformation();
    const auto& geo = conf.ring_geometries[ringIdx];
    const auto& ring = protein->RingAt(ringIdx);

    Vec3 center = geo.center;
    Vec3 normal = geo.normal.normalized();

    // Build orthonormal basis in ring plane
    Vec3 const arbitrary = (std::abs(normal.x()) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    Vec3 const u = normal.cross(arbitrary).normalized();
    Vec3 const v = normal.cross(u);

    // View direction: default "side" (perpendicular to normal, sees butterfly lobes)
    // "top" looks down the normal, "edge" looks along the other in-plane axis
    QString const view = cmd.value("view").toString("side");
    double const distance = cmd.value("distance").toDouble(15.0);

    Vec3 camDir;
    Vec3 upVec;
    if (view == "top") {
        camDir = normal;
        upVec = u;
    } else if (view == "edge") {
        camDir = v;
        upVec = normal;
    } else {  // "side"
        camDir = u;
        upVec = normal;
    }

    Vec3 camPos = center + distance * camDir;

    auto* camera = mainWindow_->renderer_->GetActiveCamera();
    camera->SetPosition(camPos.x(), camPos.y(), camPos.z());
    camera->SetFocalPoint(center.x(), center.y(), center.z());
    camera->SetViewUp(upVec.x(), upVec.y(), upVec.z());
    camera->OrthogonalizeViewUp();
    mainWindow_->renderer_->ResetCameraClippingRange();
    mainWindow_->renderWindow_->Render();

    QJsonObject result;
    result["ring_index"] = ringIdx;
    result["ring_type"] = QString::fromStdString(ring.TypeName());
    result["center"] = QJsonArray{center.x(), center.y(), center.z()};
    result["normal"] = QJsonArray{normal.x(), normal.y(), normal.z()};
    result["view"] = view;
    return QJsonObject{{"ok", true}, {"result", result}};
}

QJsonObject RestServer::cmdLookAtAtom(const QJsonObject& cmd) {
    auto& protein = mainWindow_->protein_;
    if (!protein) {
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};
    }

    int const atomIdx = cmd["atom"].toInt(-1);
    if (atomIdx < 0 || atomIdx >= static_cast<int>(protein->AtomCount())) {
        return QJsonObject{{"ok", false}, {"error", "invalid atom index"}};
    }

    const auto& conf = protein->Conformation();
    Vec3 pos = conf.AtomAt(atomIdx).Position();
    double const distance = cmd.value("distance").toDouble(15.0);

    // Keep current camera direction, just re-center on this atom
    auto* camera = mainWindow_->renderer_->GetActiveCamera();
    double camPos[3];
    double foc[3];
    camera->GetPosition(camPos);
    camera->GetFocalPoint(foc);
    Vec3 dir(camPos[0] - foc[0], camPos[1] - foc[1], camPos[2] - foc[2]);
    dir.normalize();

    Vec3 newCam = pos + distance * dir;
    camera->SetFocalPoint(pos.x(), pos.y(), pos.z());
    camera->SetPosition(newCam.x(), newCam.y(), newCam.z());
    camera->OrthogonalizeViewUp();
    mainWindow_->renderer_->ResetCameraClippingRange();
    mainWindow_->renderWindow_->Render();

    const auto& id = protein->AtomAt(atomIdx);
    const auto& res = protein->ResidueAt(id.residue_index);
    QJsonObject result;
    result["atom_index"] = atomIdx;
    result["element"] = QString::fromStdString(SymbolForElement(id.element));
    result["pdb_name"] = QString::fromStdString(id.pdb_atom_name);
    result["residue"] = QString("%1-%2").arg(
        QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type)))
        .arg(res.sequence_number);
    return QJsonObject{{"ok", true}, {"result", result}};
}

QJsonObject RestServer::cmdListRings(const QJsonObject&) {
    auto& protein = mainWindow_->protein_;
    if (!protein) {
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};
    }

    const auto& conf = protein->Conformation();
    QJsonArray rings;
    for (size_t i = 0; i < protein->RingCount(); ++i) {
        const auto& ring = protein->RingAt(i);
        const auto& geo = conf.ring_geometries[i];
        const auto& res = protein->ResidueAt(ring.parent_residue_index);
        QJsonObject r;
        r["index"] = static_cast<int>(i);
        r["type"] = QString::fromStdString(ring.TypeName());
        r["residue"] = QString("%1-%2").arg(
            QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type)))
            .arg(res.sequence_number);
        r["center"] = QJsonArray{geo.center.x(), geo.center.y(), geo.center.z()};
        r["normal"] = QJsonArray{geo.normal.x(), geo.normal.y(), geo.normal.z()};
        r["radius"] = geo.radius;
        r["intensity"] = ring.Intensity();
        rings.append(r);
    }
    return QJsonObject{{"ok", true}, {"result", rings}};
}

QJsonObject RestServer::cmdGetLog(const QJsonObject& cmd) {
    QString const text = mainWindow_->logText_->toPlainText();
    QStringList all = text.split('\n');
    int const total = all.size();

    // "lines":N — last N lines (convenience for tail)
    // "first":F, "last":L — specific range [F, L] inclusive, 0-based
    // No args — everything
    int first = 0;
    int last = total - 1;
    if (cmd.contains("lines")) {
        int const n = cmd["lines"].toInt(50);
        first = std::max(0, total - n);
    } else if (cmd.contains("first") || cmd.contains("last")) {
        first = cmd.value("first").toInt(0);
        last = cmd.value("last").toInt(total - 1);
    }
    first = std::max(0, std::min(first, total - 1));
    last = std::max(first, std::min(last, total - 1));

    QJsonArray lines;
    for (int i = first; i <= last; ++i) {
        lines.append(all[i]);
    }

    QJsonObject result;
    result["total_lines"] = total;
    result["first"] = first;
    result["last"] = last;
    result["returned"] = lines.size();
    result["lines"] = lines;
    return QJsonObject{{"ok", true}, {"result", result}};
}

QJsonObject RestServer::cmdExportFeatures(const QJsonObject& cmd) {
    auto& protein = mainWindow_->protein_;
    if (!protein) {
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};
    }

    QString const path = cmd["path"].toString();
    if (path.isEmpty()) {
        return QJsonObject{{"ok", false}, {"error", "missing 'path'"}};
    }

    std::string const outDir = path.toStdString();
    int totalArrays = 0;

    std::filesystem::create_directories(outDir);
    auto& conf = protein->Conformation();
    totalArrays = nmr::ConformationResult::WriteAllFeatures(conf, outDir);

    QJsonObject result;
    result["path"] = path;
    result["arrays"] = totalArrays;
    return QJsonObject{{"ok", true}, {"result", result}};
}

// ============================================================================
// Introspection endpoints: atom_dump + list_atoms
//
// These are general-purpose API surface, not test-specific helpers. Any
// scripted consumer (test harness, batch dump tool, ad-hoc debugger)
// uses them; tests are just one client. The shape mirrors what the
// Atom Inspector tree displays, in machine-readable form.
// ============================================================================

namespace {

QJsonArray Vec3Json(const Vec3& v) {
    return QJsonArray{v.x(), v.y(), v.z()};
}

QJsonObject SphericalTensorJson(const SphericalTensor& s) {
    QJsonObject obj;
    obj["T0"] = s.T0;
    obj["T1"] = QJsonArray{s.T1[0], s.T1[1], s.T1[2]};
    obj["T2"] = QJsonArray{s.T2[0], s.T2[1], s.T2[2], s.T2[3], s.T2[4]};
    return obj;
}

// Local NameFor switches mirror MainWindow.cpp's. Specifically:
//   - NameForPlanarGroupKindJ, NameForPolarHKindJ, NameForRingSystemKindJ,
//     NameForRingPositionLabelJ intentionally diverge from MainWindow's
//     versions: this side emits machine-friendly names (`"PyrroleAlpha"`,
//     `"Imidazole_His"`, `"NotInRing"`), MainWindow emits display names
//     (`"pyrrole α"`, `"His imidazole"`, `"—"`). The orthogonality is
//     the point — JSON consumers parse on stable identifiers, humans
//     read prose.
//   - NameForAtomRoleJ and NameForBondCategoryJ (the latter inherited
//     from MainWindow) are BYTE-IDENTICAL to MainWindow's. Drift risk
//     lives here. If these ever diverge that's a bug, not a feature.
// The library-side fix (NameFor* in Types.h next to SymbolForElement)
// is the clean answer; until that lands, duplication is preferred over
// a ui/-only shared header per PATTERNS §17.

const char* NameForAtomRoleJ(AtomRole r) {
    switch (r) {
    case AtomRole::BackboneN:
        return "BackboneN";
    case AtomRole::BackboneCA:
        return "BackboneCA";
    case AtomRole::BackboneC:
        return "BackboneC";
    case AtomRole::BackboneO:
        return "BackboneO";
    case AtomRole::SidechainC:
        return "SidechainC";
    case AtomRole::SidechainN:
        return "SidechainN";
    case AtomRole::SidechainO:
        return "SidechainO";
    case AtomRole::SidechainS:
        return "SidechainS";
    case AtomRole::AromaticC:
        return "AromaticC";
    case AtomRole::AromaticN:
        return "AromaticN";
    case AtomRole::AmideH:
        return "AmideH";
    case AtomRole::AlphaH:
        return "AlphaH";
    case AtomRole::MethylH:
        return "MethylH";
    case AtomRole::AromaticH:
        return "AromaticH";
    case AtomRole::HydroxylH:
        return "HydroxylH";
    case AtomRole::OtherH:
        return "OtherH";
    case AtomRole::Unknown:
        return "Unknown";
    }
    return "?";
}

const char* NameForBondCategoryJ(BondCategory c) {
    switch (c) {
    case BondCategory::PeptideCO:
        return "PeptideCO";
    case BondCategory::PeptideCN:
        return "PeptideCN";
    case BondCategory::BackboneOther:
        return "BackboneOther";
    case BondCategory::SidechainCO:
        return "SidechainCO";
    case BondCategory::Aromatic:
        return "Aromatic";
    case BondCategory::Disulfide:
        return "Disulfide";
    case BondCategory::SidechainOther:
        return "SidechainOther";
    case BondCategory::Unknown:
        return "Unknown";
    }
    return "?";
}

const char* NameForPlanarGroupKindJ(PlanarGroupKind k) {
    switch (k) {
    case PlanarGroupKind::None:
        return "None";
    case PlanarGroupKind::PeptideAmide:
        return "PeptideAmide";
    case PlanarGroupKind::SidechainAmide:
        return "SidechainAmide";
    case PlanarGroupKind::Guanidinium:
        return "Guanidinium";
    case PlanarGroupKind::Imidazole:
        return "Imidazole";
    case PlanarGroupKind::Aromatic6Ring:
        return "Aromatic6Ring";
    case PlanarGroupKind::Aromatic5Ring:
        return "Aromatic5Ring";
    case PlanarGroupKind::Carboxylate:
        return "Carboxylate";
    case PlanarGroupKind::AromaticHydroxyl:
        return "AromaticHydroxyl";
    case PlanarGroupKind::AromaticOxide:
        return "AromaticOxide";
    }
    return "?";
}

const char* NameForPolarHKindJ(PolarHKind k) {
    switch (k) {
    case PolarHKind::NotPolar:
        return "NotPolar";
    case PolarHKind::BackboneAmide:
        return "BackboneAmide";
    case PolarHKind::SidechainPrimaryAmide:
        return "SidechainPrimaryAmide";
    case PolarHKind::IndoleNH:
        return "IndoleNH";
    case PolarHKind::AmmoniumNH:
        return "AmmoniumNH";
    case PolarHKind::GuanidiniumNH:
        return "GuanidiniumNH";
    case PolarHKind::ImidazoleNH:
        return "ImidazoleNH";
    case PolarHKind::CarboxylOH:
        return "CarboxylOH";
    case PolarHKind::HydroxylOH_Aliphatic:
        return "HydroxylOH_Aliphatic";
    case PolarHKind::HydroxylOH_Aromatic:
        return "HydroxylOH_Aromatic";
    case PolarHKind::ThiolSH:
        return "ThiolSH";
    case PolarHKind::AmineNH:
        return "AmineNH";
    case PolarHKind::OtherPolarH:
        return "OtherPolarH";
    }
    return "?";
}

const char* NameForProchiralStereoJ(ProchiralStereo s) {
    switch (s) {
    case ProchiralStereo::NotProchiral:
        return "NotProchiral";
    case ProchiralStereo::ProR:
        return "ProR";
    case ProchiralStereo::ProS:
        return "ProS";
    case ProchiralStereo::Unassigned:
        return "Unassigned";
    }
    return "?";
}

const char* NameForPseudoatomKindJ(PseudoatomKind k) {
    switch (k) {
    case PseudoatomKind::None:
        return "None";
    case PseudoatomKind::M:
        return "M";
    case PseudoatomKind::Q:
        return "Q";
    case PseudoatomKind::R:
        return "R";
    }
    return "?";
}

const char* NameForRingSystemKindJ(RingSystemKind k) {
    switch (k) {
    case RingSystemKind::NotInRing:
        return "NotInRing";
    case RingSystemKind::Benzene_Phe:
        return "Benzene_Phe";
    case RingSystemKind::Benzene_Tyr:
        return "Benzene_Tyr";
    case RingSystemKind::Imidazole_His:
        return "Imidazole_His";
    case RingSystemKind::Indole_Trp_5:
        return "Indole_Trp_5";
    case RingSystemKind::Indole_Trp_6:
        return "Indole_Trp_6";
    case RingSystemKind::Pyrrolidine_Pro:
        return "Pyrrolidine_Pro";
    case RingSystemKind::Indole_Trp_9:
        return "Indole_Trp_9";
    }
    return "?";
}

const char* NameForRingPositionLabelJ(RingPositionLabel p) {
    switch (p) {
    case RingPositionLabel::NotInRing:
        return "NotInRing";
    case RingPositionLabel::Ipso:
        return "Ipso";
    case RingPositionLabel::Ortho1:
        return "Ortho1";
    case RingPositionLabel::Ortho2:
        return "Ortho2";
    case RingPositionLabel::Meta1:
        return "Meta1";
    case RingPositionLabel::Meta2:
        return "Meta2";
    case RingPositionLabel::Para:
        return "Para";
    case RingPositionLabel::PyrroleAlpha:
        return "PyrroleAlpha";
    case RingPositionLabel::PyrroleBeta:
        return "PyrroleBeta";
    case RingPositionLabel::BridgeFusion:
        return "BridgeFusion";
    case RingPositionLabel::Heteroatom_NH:
        return "Heteroatom_NH";
    case RingPositionLabel::Heteroatom_NoH:
        return "Heteroatom_NoH";
    case RingPositionLabel::Heteroatom_OH:
        return "Heteroatom_OH";
    case RingPositionLabel::Saturated:
        return "Saturated";
    case RingPositionLabel::ProRingNitrogen:
        return "ProRingNitrogen";
    case RingPositionLabel::ProRingAlphaCarbon:
        return "ProRingAlphaCarbon";
    case RingPositionLabel::ProRingBeta:
        return "ProRingBeta";
    case RingPositionLabel::ProRingPuckerPivot:
        return "ProRingPuckerPivot";
    case RingPositionLabel::ProRingDelta:
        return "ProRingDelta";
    case RingPositionLabel::PerimeterMember:
        return "PerimeterMember";
    }
    return "?";
}

QJsonObject RingMembershipJson(const RingMembership& m) {
    QJsonObject obj;
    obj["ring_system"] = NameForRingSystemKindJ(m.ring);
    obj["position"] = NameForRingPositionLabelJ(m.position);
    obj["ring_size"] = static_cast<int>(m.ring_size);
    obj["aromatic"] = m.aromatic;
    obj["planar"] = m.planar;
    obj["n_heteroatoms"] = static_cast<int>(m.n_heteroatoms);
    return obj;
}

}  // namespace

QJsonObject RestServer::cmdAtomDump(const QJsonObject& cmd) {
    auto& protein = mainWindow_->protein_;
    if (!protein) {
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};
    }

    int const atomIdx = cmd.value("atom").toInt(-1);
    if (atomIdx < 0 || atomIdx >= static_cast<int>(protein->AtomCount())) {
        return QJsonObject{
            {"ok", false},
            {"error", QString("atom index out of range: %1 (n_atoms=%2)").arg(atomIdx).arg(protein->AtomCount())}};
    }

    const auto& conf = protein->Conformation();
    const auto& id = protein->AtomAt(atomIdx);
    const auto& ca = conf.AtomAt(atomIdx);
    const auto& res = protein->ResidueAt(id.residue_index);

    QJsonObject result;
    result["index"] = atomIdx;

    // Identity
    {
        QJsonObject identity;
        identity["element"] = QString::fromStdString(SymbolForElement(id.element));
        identity["atomic_number"] = AtomicNumberForElement(id.element);
        identity["pdb_atom_name"] = QString::fromStdString(id.pdb_atom_name);
        identity["residue_index"] = static_cast<int>(id.residue_index);
        identity["residue_type"] = QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type));
        identity["residue_sequence"] = res.sequence_number;
        identity["chain_id"] = QString::fromStdString(res.chain_id);
        identity["terminal_state"] = ResidueTerminalStateName(res.terminal_state);
        identity["protonation_variant_index"] = res.protonation_variant_index;
        if (res.protonation_variant_index >= 0) {
            const auto& aaType = GetAminoAcidType(res.type);
            if (res.protonation_variant_index < static_cast<int>(aaType.variants.size())) {
                const auto& v = aaType.variants[res.protonation_variant_index];
                identity["protonation_variant_name"] = QString::fromUtf8(v.name);
                identity["protonation_variant_description"] = QString::fromUtf8(v.description);
            }
        }
        identity["role"] = NameForAtomRoleJ(ca.role);
        identity["is_backbone"] = ca.is_backbone;
        identity["is_amide_H"] = ca.is_amide_H;
        identity["is_alpha_H"] = ca.is_alpha_H;
        identity["is_aromatic_H"] = ca.is_aromatic_H;
        identity["is_methyl"] = ca.is_methyl;
        identity["is_on_aromatic_residue"] = ca.is_on_aromatic_residue;
        identity["is_hbond_donor"] = ca.is_hbond_donor;
        identity["is_hbond_acceptor"] = ca.is_hbond_acceptor;
        identity["parent_atom_index"] = (id.parent_atom_index == SIZE_MAX) ? -1 : static_cast<int>(id.parent_atom_index);
        identity["position"] = Vec3Json(ca.Position());
        result["identity"] = identity;
    }

    // AMBER substrate
    {
        const auto& topo = protein->LegacyAmber();
        if (topo.HasAtomSemantic()) {
            const auto& sem = topo.SemanticAt(atomIdx);
            QJsonObject sub;
            sub["planar_group"] = NameForPlanarGroupKindJ(sem.planar_group);
            sub["polar_h"] = NameForPolarHKindJ(sem.polar_h);
            sub["prochiral"] = NameForProchiralStereoJ(sem.prochiral);
            sub["formal_charge"] = sem.formal_charge;
            sub["is_exchangeable"] = sem.is_exchangeable;
            sub["aromatic_flag"] = sem.aromatic;
            sub["equivalence_class"] = sem.equivalence_class;

            QJsonObject rp;
            rp["in_any_ring"] = sem.ring_position.IsInAnyRing();
            rp["membership_count"] = sem.ring_position.MembershipCount();
            if (sem.ring_position.HasPrimaryRing()) {
                rp["primary"] = RingMembershipJson(sem.ring_position.primary);
            }
            if (sem.ring_position.HasSecondaryRing()) {
                rp["secondary"] = RingMembershipJson(sem.ring_position.secondary);
            }
            if (sem.ring_position.HasTertiaryRing()) {
                rp["tertiary"] = RingMembershipJson(sem.ring_position.tertiary);
            }
            sub["ring_position"] = rp;

            QJsonObject ps;
            ps["kind"] = NameForPseudoatomKindJ(sem.pseudoatom.kind);
            ps["locant"] = static_cast<int>(sem.pseudoatom.locant);
            ps["branch"] = static_cast<int>(sem.pseudoatom.branch);
            ps["in_super_group"] = sem.pseudoatom.in_super_group;
            sub["pseudoatom"] = ps;

            result["amber_substrate"] = sub;
        }
    }

    // Charges
    {
        QJsonObject c;
        c["ff14sb"] = ca.partial_charge;
        c["aimnet2_hirshfeld"] = ca.aimnet2_charge;
        c["eeq_charge"] = ca.eeq_charge;
        c["eeq_cn"] = ca.eeq_cn;
        c["mopac_pm7"] = ca.mopac_charge;
        c["mopac_s_pop"] = ca.mopac_s_pop;
        c["mopac_p_pop"] = ca.mopac_p_pop;
        c["mopac_valency"] = ca.mopac_valency;
        c["pb_radius_A"] = ca.pb_radius;
        result["charges"] = c;
    }

    // Shielding contributions (all SphericalTensor)
    {
        QJsonObject sh;
        sh["bs"] = SphericalTensorJson(ca.bs_shielding_contribution);
        sh["hm"] = SphericalTensorJson(ca.hm_shielding_contribution);
        sh["mc"] = SphericalTensorJson(ca.mc_shielding_contribution);
        sh["dispersion"] = SphericalTensorJson(ca.disp_shielding_contribution);
        sh["coulomb"] = SphericalTensorJson(ca.coulomb_shielding_contribution);
        sh["piquad"] = SphericalTensorJson(ca.piquad_shielding_contribution);
        sh["ringchi"] = SphericalTensorJson(ca.ringchi_shielding_contribution);
        sh["hbond_kernel"] = SphericalTensorJson(ca.hbond_shielding_contribution);
        sh["hbond_larsen"] = SphericalTensorJson(ca.larsen_hbond_shielding_spherical);
        sh["tripeptide_bb"] = SphericalTensorJson(ca.tripeptide_bb_shielding_spherical);
        sh["tripeptide_neighbor"] = SphericalTensorJson(ca.tripeptide_neighbor_shielding_spherical);
        sh["aimnet2_efg"] = SphericalTensorJson(ca.aimnet2_shielding_contribution);
        sh["mopac_coulomb"] = SphericalTensorJson(ca.mopac_coulomb_shielding_contribution);
        sh["mopac_mc"] = SphericalTensorJson(ca.mopac_mc_shielding_contribution);
        result["shielding_contributions"] = sh;
    }

    // Vector fields + EFG tensors
    {
        QJsonObject vf;
        vf["B_BS"] = Vec3Json(ca.total_B_field);
        vf["E_ff14sb"] = Vec3Json(ca.coulomb_E_total);
        vf["E_ff14sb_magnitude"] = ca.coulomb_E_magnitude;
        vf["E_ff14sb_backbone_frac"] = ca.coulomb_E_backbone_frac;
        vf["E_mopac"] = Vec3Json(ca.mopac_coulomb_E_total);
        vf["E_apbs"] = Vec3Json(ca.apbs_efield);
        vf["EFG_apbs"] = SphericalTensorJson(ca.apbs_efg_spherical);
        vf["EFG_aimnet2_total"] = SphericalTensorJson(ca.aimnet2_EFG_total_spherical);
        vf["EFG_aimnet2_backbone"] = SphericalTensorJson(ca.aimnet2_EFG_backbone_spherical);
        vf["EFG_aimnet2_aromatic"] = SphericalTensorJson(ca.aimnet2_EFG_aromatic_spherical);
        result["vector_fields"] = vf;
    }

    // Local geometry
    {
        QJsonObject g;
        g["pyramidalization_A"] = ca.pyramidalization;
        g["atom_sasa_A2"] = ca.atom_sasa;
        g["sasa_normal"] = Vec3Json(ca.sasa_normal);
        result["local_geometry"] = g;
    }

    // AIMNet2 polarisability (charge-response gradient)
    {
        QJsonObject p;
        p["gradient_vector"] = Vec3Json(ca.aimnet2_charge_response_gradient_vector);
        p["gradient_magnitude"] = ca.aimnet2_charge_response_gradient_scalar;
        result["aimnet2_polarisability"] = p;
    }

    // AIMNet2 embedding L2² + first 4 components (full 256-dim too noisy)
    {
        double n2 = 0.0;
        for (size_t k = 0; k < AIMNET2_AIM_DIMS; ++k) {
            const double v = static_cast<double>(ca.aimnet2_aim[k]);
            n2 += v * v;
        }
        QJsonObject e;
        e["dims"] = static_cast<int>(AIMNET2_AIM_DIMS);
        e["l2_squared"] = n2;
        QJsonArray first4;
        for (size_t k = 0; k < 4; ++k) {
            first4.append(static_cast<double>(ca.aimnet2_aim[k]));
        }
        e["first_4"] = first4;
        result["aimnet2_embedding"] = e;
    }

    // Graph topology (MolecularGraphResult)
    {
        QJsonObject g;
        g["graph_dist_ring"] = ca.graph_dist_ring;
        g["graph_dist_N"] = ca.graph_dist_N;
        g["graph_dist_O"] = ca.graph_dist_O;
        g["bfs_decay"] = ca.bfs_decay;
        g["eneg_sum_1"] = ca.eneg_sum_1;
        g["eneg_sum_2"] = ca.eneg_sum_2;
        g["n_pi_bonds_3"] = ca.n_pi_bonds_3;
        g["is_conjugated"] = ca.is_conjugated;
        result["graph_topology"] = g;
    }

    // Tripeptide DFT (BB + Neighbor)
    {
        QJsonObject t;
        QJsonObject bb;
        bb["has_match"] = ca.tripeptide_bb_has_match;
        bb["shielding"] = SphericalTensorJson(ca.tripeptide_bb_shielding_spherical);
        bb["match_distance_A"] = ca.tripeptide_bb_match_distance;
        bb["residual_vec"] = Vec3Json(ca.tripeptide_bb_residual_vec);
        bb["method_tag"] = static_cast<int>(ca.tripeptide_bb_method_tag);
        bb["method_name"] = ca.tripeptide_bb_method_tag == 1   ? "gaussian_standard_orientation"
                            : ca.tripeptide_bb_method_tag == 2 ? "orca_input_orientation"
                                                               : "none";
        t["backbone"] = bb;

        QJsonObject nb;
        nb["has_match"] = ca.tripeptide_neighbor_has_match;
        nb["shielding_sum"] = SphericalTensorJson(ca.tripeptide_neighbor_shielding_spherical);
        nb["residual_vec_prev"] = Vec3Json(ca.tripeptide_neighbor_residual_vec_prev);
        nb["residual_vec_next"] = Vec3Json(ca.tripeptide_neighbor_residual_vec_next);
        t["neighbor"] = nb;
        result["tripeptide_dft"] = t;
    }

    // Larsen H-bond (4 contribution classes + water term + diagnostic)
    {
        QJsonObject lhb;
        lhb["n_pairs"] = ca.larsen_hbond_n_pairs;
        lhb["water_term_ppm"] = ca.larsen_hbond_water_term;
        lhb["any_corner_imputed"] = ca.larsen_hbond_any_corner_imputed;
        lhb["sum"] = SphericalTensorJson(ca.larsen_hbond_shielding_spherical);
        lhb["delta_1pHB"] = SphericalTensorJson(ca.larsen_hbond_1pHB_spherical);
        lhb["delta_2pHB"] = SphericalTensorJson(ca.larsen_hbond_2pHB_spherical);
        lhb["delta_1pHaB"] = SphericalTensorJson(ca.larsen_hbond_1pHaB_spherical);
        lhb["delta_2pHaB"] = SphericalTensorJson(ca.larsen_hbond_2pHaB_spherical);
        lhb["diagnostic_CB"] = SphericalTensorJson(ca.larsen_hbond_diagnostic_CB_spherical);
        result["larsen_hbond"] = lhb;
    }

    // HBond kernel info
    {
        QJsonObject h;
        h["nearest_dist_A"] = ca.hbond_nearest_dist;
        h["nearest_dir"] = Vec3Json(ca.hbond_nearest_dir);
        h["count_within_3_5A"] = ca.hbond_count_within_3_5A;
        h["is_donor"] = ca.hbond_is_donor;
        h["is_acceptor"] = ca.hbond_is_acceptor;
        h["is_backbone"] = ca.hbond_is_backbone;
        result["hbond_kernel"] = h;
    }

    // ORCA DFT (if loaded — only in --orca / mutant modes)
    if (ca.has_orca_shielding) {
        QJsonObject o;
        o["total"] = SphericalTensorJson(ca.orca_shielding_total_spherical);
        o["diamagnetic"] = SphericalTensorJson(ca.orca_shielding_diamagnetic_spherical);
        o["paramagnetic"] = SphericalTensorJson(ca.orca_shielding_paramagnetic_spherical);
        result["orca_dft"] = o;
    }

    // Ring neighbours (per-atom-per-ring full structured array)
    {
        QJsonArray arr;
        for (const auto& rn : ca.ring_neighbours) {
            QJsonObject r;
            r["ring_index"] = static_cast<int>(rn.ring_index);
            r["ring_type"] = QString::fromUtf8(protein->RingAt(rn.ring_index).TypeName());
            r["distance_to_center_A"] = rn.distance_to_center;
            r["direction_to_center"] = Vec3Json(rn.direction_to_center);
            r["rho_A"] = rn.rho;
            r["z_A"] = rn.z;
            r["theta_rad"] = rn.theta;
            r["G_spherical"] = SphericalTensorJson(rn.G_spherical);
            r["hm_G_spherical"] = SphericalTensorJson(rn.hm_G_spherical);
            r["B_field"] = Vec3Json(rn.B_field);
            r["B_cylindrical"] = Vec3Json(rn.B_cylindrical);
            r["chi_scalar"] = rn.chi_scalar;
            r["quad_scalar"] = rn.quad_scalar;
            r["disp_scalar"] = rn.disp_scalar;
            r["disp_contacts"] = rn.disp_contacts;
            arr.append(r);
        }
        result["ring_neighbours"] = arr;
    }

    // Bonds — three parallel arrays mirroring the Atom Inspector's
    // bond-side panels:
    //   covalent:  direct covalent topology (from id.bond_indices),
    //              with category / length / direction / order / midpoint.
    //   mcconnell: McConnell dipolar bond neighbours within cutoff
    //              (from ca.bond_neighbours), one per nearby bond.
    //   mopac:     MOPAC Wiberg bond orders to specific other atoms
    //              (from ca.mopac_bond_neighbours).
    // The MopacResult-derived Wiberg bond-order lookup is reused across
    // covalent + mcconnell when MopacResult is attached (it isn't in
    // viewer's default skip_mopac=true path; the bond_order field stays
    // 0.0).
    {
        QJsonObject bonds;
        const bool has_mopac = conf.HasResult<MopacResult>();

        // covalent: every bond this atom participates in
        {
            QJsonArray arr;
            for (size_t const bi : id.bond_indices) {
                if (bi >= protein->BondCount()) {
                    continue;  // defensive
                }
                const Bond& bond = protein->BondAt(bi);
                const size_t other_idx = (bond.atom_index_a == static_cast<size_t>(atomIdx)) ? bond.atom_index_b
                                                                                             : bond.atom_index_a;
                const auto& other = protein->AtomAt(other_idx);
                QJsonObject b;
                b["bond_index"] = static_cast<int>(bi);
                b["other_atom"] = static_cast<int>(other_idx);
                b["other_element"] = QString::fromStdString(SymbolForElement(other.element));
                b["other_pdb_name"] = QString::fromStdString(other.pdb_atom_name);
                b["category"] = NameForBondCategoryJ(bond.category);
                b["length_A"] = conf.bond_lengths[bi];
                b["direction"] = Vec3Json(conf.bond_directions[bi]);
                b["midpoint"] = Vec3Json(conf.bond_midpoints[bi]);
                b["is_rotatable"] = bond.is_rotatable;
                if (has_mopac) {
                    b["wiberg_order"] = conf.Result<MopacResult>().TopologyBondOrder(bi);
                }
                arr.append(b);
            }
            bonds["covalent"] = arr;
        }

        // mcconnell: bond_neighbours has every bond contributing
        // dipolar shielding at this atom (within cutoff)
        {
            QJsonArray arr;
            for (const auto& bn : ca.bond_neighbours) {
                if (bn.bond_index >= protein->BondCount()) {
                    continue;
                }
                const Bond& bond = protein->BondAt(bn.bond_index);
                QJsonObject b;
                b["bond_index"] = static_cast<int>(bn.bond_index);
                b["category"] = NameForBondCategoryJ(bn.bond_category);
                b["distance_to_midpoint_A"] = bn.distance_to_midpoint;
                b["direction_to_midpoint"] = Vec3Json(bn.direction_to_midpoint);
                b["mcconnell_scalar"] = bn.mcconnell_scalar;
                b["dipolar_spherical"] = SphericalTensorJson(bn.dipolar_spherical);
                // Topology-bond category cross-check: bn.bond_category
                // should equal protein->BondAt(bn.bond_index).category.
                // Emit topology side too so consumers can verify.
                b["topology_category"] = NameForBondCategoryJ(bond.category);
                if (has_mopac) {
                    b["wiberg_order"] = conf.Result<MopacResult>().TopologyBondOrder(bn.bond_index);
                }
                arr.append(b);
            }
            bonds["mcconnell"] = arr;
        }

        // mopac: per-other-atom Wiberg bond orders (sorted descending in
        // the library). Empty when MOPAC is skipped.
        {
            QJsonArray arr;
            for (const auto& mb : ca.mopac_bond_neighbours) {
                QJsonObject b;
                b["other_atom"] = static_cast<int>(mb.other_atom);
                b["wiberg_order"] = mb.wiberg_order;
                b["topology_bond_index"] = (mb.topology_bond_index == SIZE_MAX) ? -1 : static_cast<int>(mb.topology_bond_index);
                if (mb.topology_bond_index != SIZE_MAX && mb.topology_bond_index < protein->BondCount()) {
                    b["topology_category"] = NameForBondCategoryJ(protein->BondAt(mb.topology_bond_index).category);
                }
                arr.append(b);
            }
            bonds["mopac"] = arr;
        }

        result["bonds"] = bonds;
    }

    // DSSP per-residue snapshot for this atom's residue
    if (conf.HasResult<DsspResult>()) {
        const auto& dssp = conf.Result<DsspResult>();
        QJsonObject d;
        d["secondary_structure"] = QString(QChar(dssp.SecondaryStructure(id.residue_index)));
        d["phi_deg"] = dssp.Phi(id.residue_index);
        d["psi_deg"] = dssp.Psi(id.residue_index);
        d["sasa_A2"] = dssp.SASA(id.residue_index);
        result["dssp"] = d;
    }

    // Trajectory H5 companion — present iff --analysis-h5 supplied AND
    // the identity check passed at ComputeWorker Phase 2b. Section
    // absent (not null) otherwise; consumers gate on .get("h5"). One
    // sub-key per TR group in the file (sparse-tolerant: absent groups
    // simply don't appear in the payload).
    {
        const auto& binding = mainWindow_->analysisBinding_;
        if (binding.Valid()) {
            const auto& h5 = *binding.h5;
            const size_t h5idx = binding.H5IndexFor(static_cast<size_t>(atomIdx));

            QJsonObject h5j;
            h5j["protein_id"] = QString::fromStdString(h5.ProteinId());
            h5j["n_atoms"] = static_cast<int>(h5.AtomCount());
            h5j["n_frames"] = static_cast<int>(h5.FrameCount());
            h5j["frame_time_ps_0"] = h5.FrameTimePs(0);
            h5j["h5_atom_index"] = static_cast<int>(h5idx);
            h5j["h5_atom_name"] = QString::fromStdString(h5.AtomNameAt(h5idx));
            h5j["h5_element"] = h5.ElementAt(h5idx);

            QJsonArray groups;
            for (const auto& g : h5.GroupsPresent()) {
                groups.append(QString::fromStdString(g));
            }
            h5j["groups_present"] = groups;

            auto addShieldingWelford = [&](const QString& key, std::optional<TrajectoryH5::ShieldingWelfordRow> w) {
                if (!w)
                    return;
                QJsonObject s;
                s["t0_mean"] = w->t0.mean;
                s["t0_std"] = w->t0.std;
                s["t2magnitude_mean"] = w->t2magnitude.mean;
                s["t2magnitude_std"] = w->t2magnitude.std;
                h5j[key] = s;
            };
            auto addShieldingFrame0 = [&](const QString& key, std::optional<TrajectoryH5::ShieldingFrame0Row> f) {
                if (!f)
                    return;
                QJsonObject s;
                s["T0"] = f->T0;
                s["T2_magnitude"] = f->T2_magnitude;
                h5j[key] = s;
            };

            addShieldingWelford("bs_welford", h5.BsWelford(h5idx));
            addShieldingWelford("hm_welford", h5.HmWelford(h5idx));
            addShieldingWelford("mc_welford", h5.McWelford(h5idx));

            addShieldingFrame0("bs_shielding_frame_0", h5.BsShieldingFrame0(h5idx));
            addShieldingFrame0("hm_shielding_frame_0", h5.HmShieldingFrame0(h5idx));
            addShieldingFrame0("mc_shielding_frame_0", h5.McShieldingFrame0(h5idx));
            addShieldingFrame0("piquad_shielding_frame_0", h5.PiQuadShieldingFrame0(h5idx));
            addShieldingFrame0("ringchi_shielding_frame_0", h5.RingChiShieldingFrame0(h5idx));
            addShieldingFrame0("disp_shielding_frame_0", h5.DispShieldingFrame0(h5idx));
            addShieldingFrame0("hbond_shielding_frame_0", h5.HBondShieldingFrame0(h5idx));

            if (auto sw = h5.SasaWelford(h5idx)) {
                h5j["sasa_welford"] = QJsonObject{{"mean", sw->sasa.mean}, {"std", sw->sasa.std}};
            }
            if (auto sf = h5.SasaFrame0(h5idx)) {
                h5j["sasa_frame_0"] = *sf;
            }
            if (auto ew = h5.EeqWelford(h5idx)) {
                h5j["eeq_welford"] = QJsonObject{{"mean", ew->charge.mean}, {"std", ew->charge.std}};
            }
            if (auto ac = h5.Aimnet2ChargeFrame0(h5idx)) {
                h5j["aimnet2_charge_frame_0"] = *ac;
            }
            if (auto hc = h5.HBondCountWelford(h5idx)) {
                h5j["hbond_count_welford"] = QJsonObject{{"mean", hc->count.mean}, {"std", hc->count.std}};
            }
            if (auto ef = h5.ApbsEfieldFrame0(h5idx)) {
                h5j["apbs_efield_frame_0"] = QJsonObject{{"x", ef->x}, {"y", ef->y}, {"z", ef->z}};
            }

            result["h5"] = h5j;
        }
    }

    return QJsonObject{{"ok", true}, {"result", result}};
}

QJsonObject RestServer::cmdListAtoms(const QJsonObject& cmd) {
    auto& protein = mainWindow_->protein_;
    if (!protein) {
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};
    }

    const QJsonObject filter = cmd.value("filter").toObject();
    const QString fRole = filter.value("role").toString();
    const QString fElement = filter.value("element").toString();
    const int fResidueIdx = filter.contains("residue_index") ? filter.value("residue_index").toInt() : -1;
    const QString fResidueType = filter.value("residue_type").toString();
    const int fResidueSeq = filter.contains("residue_sequence") ? filter.value("residue_sequence").toInt() : -1;
    const bool hasInRing = filter.contains("in_ring");
    const bool fInRing = hasInRing ? filter.value("in_ring").toBool() : false;
    const QString fPlanarGroup = filter.value("planar_group").toString();
    const QString fPolarH = filter.value("polar_h").toString();
    const bool hasIsAmideH = filter.contains("is_amide_H");
    const bool fIsAmideH = hasIsAmideH ? filter.value("is_amide_H").toBool() : false;
    const bool hasIsMethyl = filter.contains("is_methyl");
    const bool fIsMethyl = hasIsMethyl ? filter.value("is_methyl").toBool() : false;

    const auto& conf = protein->Conformation();
    const auto& topo = protein->LegacyAmber();
    const bool hasSubstrate = topo.HasAtomSemantic();

    QJsonArray arr;
    for (size_t i = 0; i < protein->AtomCount(); ++i) {
        const auto& id = protein->AtomAt(i);
        const auto& ca = conf.AtomAt(i);
        const auto& res = protein->ResidueAt(id.residue_index);

        if (!fRole.isEmpty() && fRole != NameForAtomRoleJ(ca.role)) {
            continue;
        }
        if (!fElement.isEmpty() && fElement != QString::fromStdString(SymbolForElement(id.element))) {
            continue;
        }
        if (fResidueIdx >= 0 && static_cast<int>(id.residue_index) != fResidueIdx) {
            continue;
        }
        if (!fResidueType.isEmpty() && fResidueType != QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type))) {
            continue;
        }
        if (fResidueSeq >= 0 && res.sequence_number != fResidueSeq) {
            continue;
        }
        if (hasIsAmideH && fIsAmideH != ca.is_amide_H) {
            continue;
        }
        if (hasIsMethyl && fIsMethyl != ca.is_methyl) {
            continue;
        }

        if (hasSubstrate) {
            const auto& sem = topo.SemanticAt(i);
            if (hasInRing && fInRing != sem.ring_position.IsInAnyRing()) {
                continue;
            }
            if (!fPlanarGroup.isEmpty() && fPlanarGroup != NameForPlanarGroupKindJ(sem.planar_group)) {
                continue;
            }
            if (!fPolarH.isEmpty() && fPolarH != NameForPolarHKindJ(sem.polar_h)) {
                continue;
            }
        } else if (hasInRing || !fPlanarGroup.isEmpty() || !fPolarH.isEmpty()) {
            // Substrate filter requested but substrate not populated → no match
            continue;
        }

        QJsonObject a;
        a["index"] = static_cast<int>(i);
        a["element"] = QString::fromStdString(SymbolForElement(id.element));
        a["pdb_atom_name"] = QString::fromStdString(id.pdb_atom_name);
        a["residue_index"] = static_cast<int>(id.residue_index);
        a["residue_type"] = QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type));
        a["residue_sequence"] = res.sequence_number;
        a["role"] = NameForAtomRoleJ(ca.role);
        if (hasSubstrate) {
            const auto& sem = topo.SemanticAt(i);
            a["planar_group"] = NameForPlanarGroupKindJ(sem.planar_group);
            a["polar_h"] = NameForPolarHKindJ(sem.polar_h);
            a["in_ring"] = sem.ring_position.IsInAnyRing();
        }
        arr.append(a);
    }

    QJsonObject result;
    result["count"] = arr.size();
    result["atoms"] = arr;
    return QJsonObject{{"ok", true}, {"result", result}};
}
