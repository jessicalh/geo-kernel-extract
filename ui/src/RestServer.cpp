#include "RestServer.h"
#include "MainWindow.h"

// Library headers — REST server reads the protein model directly
#include "Protein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "Ring.h"
#include "Atom.h"
#include "Residue.h"
#include "Types.h"
#include "ConformationResult.h"
#include "JobSpec.h"
#include "FieldGridOverlay.h"

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
        QByteArray line = socket->readLine().trimmed();
        if (line.isEmpty()) continue;

        QJsonParseError err;
        QJsonDocument doc = QJsonDocument::fromJson(line, &err);
        if (err.error != QJsonParseError::NoError || !doc.isObject()) {
            QJsonObject resp;
            resp["ok"] = false;
            resp["error"] = QString("JSON parse error: %1").arg(err.errorString());
            sendResponse(socket, resp);
            continue;
        }

        QJsonObject result = dispatch(doc.object());
        sendResponse(socket, result);
    }
}

void RestServer::onDisconnected() {
    auto* socket = qobject_cast<QTcpSocket*>(sender());
    if (socket) {
        clients_.removeAll(socket);
        socket->deleteLater();
    }
}

void RestServer::sendResponse(QTcpSocket* socket, const QJsonObject& response) {
    QJsonDocument doc(response);
    socket->write(doc.toJson(QJsonDocument::Compact));
    socket->write("\n");
    socket->flush();
}

QJsonObject RestServer::dispatch(const QJsonObject& cmd) {
    QString action = cmd["cmd"].toString();

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

    QJsonObject resp;
    resp["ok"] = false;
    resp["error"] = QString("Unknown command: %1").arg(action);
    return resp;
}

// ---- Command implementations ----

QJsonObject RestServer::cmdStatus() {
    QJsonObject resp;
    resp["ok"] = true;
    QJsonObject result;

    // Report computation state
    result["computing"] = (mainWindow_->workerThread_ != nullptr &&
                           mainWindow_->workerThread_->isRunning());

    // Read directly from the library protein model
    auto& protein = mainWindow_->protein_;
    if (protein) {
        const auto& conf = protein->Conformation();
        result["protein"] = QString::fromStdString(mainWindow_->currentProteinId_);
        result["n_atoms"] = (int)protein->AtomCount();
        result["n_rings"] = (int)protein->RingCount();
        result["n_residues"] = (int)protein->ResidueCount();

        // Tier counts: read from library (set by PredictionResult)
        int nReport = 0, nPass = 0, nSilent = 0;
        for (size_t i = 0; i < conf.AtomCount(); ++i) {
            switch (conf.AtomAt(i).tier) {
                case HeuristicTier::REPORT: nReport++; break;
                case HeuristicTier::PASS:   nPass++;   break;
                case HeuristicTier::SILENT: nSilent++; break;
            }
        }
        result["n_report"] = nReport;
        result["n_pass"] = nPass;
        result["n_silent"] = nSilent;
    } else {
        result["protein"] = "";
        result["n_atoms"] = 0;
        result["n_rings"] = 0;
        result["n_residues"] = 0;
        result["n_report"] = 0;
        result["n_pass"] = 0;
        result["n_silent"] = 0;
    }

    result["overlay_mode"] = 0;  // overlay modes removed
    result["n_field_grids"] = (int)mainWindow_->fieldGrids_.size();
    result["n_butterfly_grids"] = (int)mainWindow_->butterflyFields_.size();
    result["field_grid_overlay"] = (mainWindow_->fieldGridOverlay_ != nullptr);
    result["butterfly_overlay"] = (mainWindow_->butterflyOverlay_ != nullptr);
    if (!mainWindow_->fieldGrids_.empty()) {
        double gmin = 1e30, gmax = -1e30;
        int nz = 0;
        for (const auto& g : mainWindow_->fieldGrids_) {
            for (double v : g.T0) {
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
    QString path = cmd["path"].toString();
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
    QString path = cmd["path"].toString();
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
    QString mode = cmd["mode"].toString();
    // Overlay modes removed — per-calculator visualizations replace them
    return QJsonObject{{"ok", true}, {"result", QJsonObject{{"mode", "none"},
        {"note", "overlay modes removed; use per-calculator toggles"}}}};
}

QJsonObject RestServer::cmdSetRenderMode(const QJsonObject& cmd) {
    QString mode = cmd["mode"].toString();
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
    QString path = cmd["path"].toString();
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
    double azimuth = cmd.value("azimuth").toDouble(0);
    double elevation = cmd.value("elevation").toDouble(0);

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
    bool visible = cmd["visible"].toBool(true);
    mainWindow_->showRingsCheck_->setChecked(visible);
    mainWindow_->onShowRingsToggled(visible);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdShowBonds(const QJsonObject& cmd) {
    bool visible = cmd["visible"].toBool(true);
    mainWindow_->showPeptideBondsCheck_->setChecked(visible);
    mainWindow_->onShowPeptideBondsToggled(visible);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdShowButterfly(const QJsonObject& cmd) {
    bool visible = cmd["visible"].toBool(true);
    mainWindow_->showButterflyCheck_->setChecked(visible);
    mainWindow_->onShowButterflyToggled(visible);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetGlyphScale(const QJsonObject& cmd) {
    double scale = cmd["scale"].toDouble(0.5);
    int sliderVal = static_cast<int>(scale * 100);
    mainWindow_->glyphScaleSlider_->setValue(sliderVal);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetOpacity(const QJsonObject& cmd) {
    double opacity = cmd["value"].toDouble(0.7);
    int sliderVal = static_cast<int>(opacity * 100);
    mainWindow_->opacitySlider_->setValue(sliderVal);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetIsoThreshold(const QJsonObject& cmd) {
    double threshold = cmd["value"].toDouble(0.1);
    int sliderVal = static_cast<int>(threshold * 100);
    mainWindow_->isoThresholdSlider_->setValue(sliderVal);
    mainWindow_->onIsoThresholdChanged();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdShowFieldGrid(const QJsonObject& cmd) {
    if (cmd.contains("shielded")) {
        bool vis = cmd["shielded"].toBool(true);
        mainWindow_->showFieldGridCheck_->setChecked(vis);
        if (mainWindow_->fieldGridOverlay_)
            mainWindow_->fieldGridOverlay_->setShieldedVisible(vis);
    }
    if (cmd.contains("deshielded")) {
        bool vis = cmd["deshielded"].toBool(true);
        mainWindow_->showDeshieldedCheck_->setChecked(vis);
        if (mainWindow_->fieldGridOverlay_)
            mainWindow_->fieldGridOverlay_->setDeshieldedVisible(vis);
    }
    if (!cmd.contains("shielded") && !cmd.contains("deshielded")) {
        bool vis = cmd["visible"].toBool(true);
        mainWindow_->showFieldGridCheck_->setChecked(vis);
        mainWindow_->showDeshieldedCheck_->setChecked(vis);
        if (mainWindow_->fieldGridOverlay_)
            mainWindow_->fieldGridOverlay_->setVisible(vis);
    }
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdGetCamera(const QJsonObject&) {
    auto* camera = mainWindow_->renderer_->GetActiveCamera();
    double pos[3], foc[3], up[3];
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
    if (!protein)
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};

    int ringIdx = cmd["ring"].toInt(-1);
    if (ringIdx < 0 || ringIdx >= (int)protein->RingCount())
        return QJsonObject{{"ok", false}, {"error", "invalid ring index"}};

    const auto& conf = protein->Conformation();
    const auto& geo = conf.ring_geometries[ringIdx];
    const auto& ring = protein->RingAt(ringIdx);

    Vec3 center = geo.center;
    Vec3 normal = geo.normal.normalized();

    // Build orthonormal basis in ring plane
    Vec3 arbitrary = (std::abs(normal.x()) < 0.9) ? Vec3(1,0,0) : Vec3(0,1,0);
    Vec3 u = normal.cross(arbitrary).normalized();
    Vec3 v = normal.cross(u);

    // View direction: default "side" (perpendicular to normal, sees butterfly lobes)
    // "top" looks down the normal, "edge" looks along the other in-plane axis
    QString view = cmd.value("view").toString("side");
    double distance = cmd.value("distance").toDouble(15.0);

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
    if (!protein)
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};

    int atomIdx = cmd["atom"].toInt(-1);
    if (atomIdx < 0 || atomIdx >= (int)protein->AtomCount())
        return QJsonObject{{"ok", false}, {"error", "invalid atom index"}};

    const auto& conf = protein->Conformation();
    Vec3 pos = conf.AtomAt(atomIdx).Position();
    double distance = cmd.value("distance").toDouble(15.0);

    // Keep current camera direction, just re-center on this atom
    auto* camera = mainWindow_->renderer_->GetActiveCamera();
    double camPos[3], foc[3];
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
    if (!protein)
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};

    const auto& conf = protein->Conformation();
    QJsonArray rings;
    for (size_t i = 0; i < protein->RingCount(); ++i) {
        const auto& ring = protein->RingAt(i);
        const auto& geo = conf.ring_geometries[i];
        const auto& res = protein->ResidueAt(ring.parent_residue_index);
        QJsonObject r;
        r["index"] = (int)i;
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
    QString text = mainWindow_->logText_->toPlainText();
    QStringList all = text.split('\n');
    int total = all.size();

    // "lines":N — last N lines (convenience for tail)
    // "first":F, "last":L — specific range [F, L] inclusive, 0-based
    // No args — everything
    int first = 0, last = total - 1;
    if (cmd.contains("lines")) {
        int n = cmd["lines"].toInt(50);
        first = std::max(0, total - n);
    } else if (cmd.contains("first") || cmd.contains("last")) {
        first = cmd.value("first").toInt(0);
        last = cmd.value("last").toInt(total - 1);
    }
    first = std::max(0, std::min(first, total - 1));
    last = std::max(first, std::min(last, total - 1));

    QJsonArray lines;
    for (int i = first; i <= last; ++i)
        lines.append(all[i]);

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
    if (!protein)
        return QJsonObject{{"ok", false}, {"error", "no protein loaded"}};

    QString path = cmd["path"].toString();
    if (path.isEmpty())
        return QJsonObject{{"ok", false}, {"error", "missing 'path'"}};

    std::string outDir = path.toStdString();
    int totalArrays = 0;

    if (mainWindow_->pendingSpec_.mode == nmr::JobMode::Fleet) {
        for (size_t i = 0; i < protein->ConformationCount(); ++i) {
            auto& conf = protein->ConformationAt(i);
            std::string frameDir = outDir + "/frame_" +
                std::to_string(i + 1);
            std::filesystem::create_directories(frameDir);
            totalArrays += nmr::ConformationResult::WriteAllFeatures(
                conf, frameDir);
        }
    } else {
        std::filesystem::create_directories(outDir);
        auto& conf = protein->Conformation();
        totalArrays = nmr::ConformationResult::WriteAllFeatures(
            conf, outDir);
    }

    QJsonObject result;
    result["path"] = path;
    result["arrays"] = totalArrays;
    if (mainWindow_->pendingSpec_.mode == nmr::JobMode::Fleet)
        result["frames"] = (int)protein->ConformationCount();
    return QJsonObject{{"ok", true}, {"result", result}};
}
