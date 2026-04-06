#include "RestServer.h"
#include "MainWindow.h"

// Library headers — REST server reads the protein model directly
#include "Protein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "Types.h"

#include <QJsonDocument>
#include <QJsonArray>
#include <QApplication>
#include <QComboBox>
#include <QCheckBox>
#include <QSlider>
#include <QDoubleSpinBox>
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

    result["overlay_mode"] = mainWindow_->overlayModeCombo_->currentIndex();
    resp["result"] = result;
    return resp;
}

QJsonObject RestServer::cmdLoadPdb(const QJsonObject& cmd) {
    QString path = cmd["path"].toString();
    if (path.isEmpty()) {
        return QJsonObject{{"ok", false}, {"error", "missing 'path'"}};
    }
    mainWindow_->loadPdb(path.toStdString());
    return QJsonObject{{"ok", true}, {"result", QJsonObject{{"status", "loading"}}}};
}

QJsonObject RestServer::cmdLoadProteinDir(const QJsonObject& cmd) {
    QString path = cmd["path"].toString();
    if (path.isEmpty()) {
        return QJsonObject{{"ok", false}, {"error", "missing 'path'"}};
    }
    mainWindow_->loadProteinDir(path.toStdString());
    return QJsonObject{{"ok", true}, {"result", QJsonObject{{"status", "loading"}}}};
}

QJsonObject RestServer::cmdSetOverlay(const QJsonObject& cmd) {
    QString mode = cmd["mode"].toString();
    static const QMap<QString, int> modes = {
        {"none", 0}, {"heuristic", 1}, {"classical", 2},
        {"dft_delta", 3}, {"residual", 4}
    };
    auto it = modes.find(mode);
    if (it == modes.end()) {
        return QJsonObject{{"ok", false},
            {"error", QString("unknown mode: %1. Options: none, heuristic, classical, dft_delta, residual").arg(mode)}};
    }
    mainWindow_->overlayModeCombo_->setCurrentIndex(it.value());
    mainWindow_->onOverlayModeChanged(it.value());
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}, {"result", QJsonObject{{"mode", mode}}}};
}

QJsonObject RestServer::cmdSetRenderMode(const QJsonObject& cmd) {
    QString mode = cmd["mode"].toString();
    static const QMap<QString, int> modes = {
        {"ball_stick", 0}, {"stick", 1}, {"vdw", 2}, {"backbone", 3}
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

    int scale = cmd.value("scale").toInt(2);

    mainWindow_->renderWindow_->Render();

    auto filter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    filter->SetInput(mainWindow_->renderWindow_);
    filter->SetScale(scale);
    filter->SetInputBufferTypeToRGBA();
    filter->ReadFrontBufferOff();
    filter->Update();

    auto writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(path.toStdString().c_str());
    writer->SetInputConnection(filter->GetOutputPort());
    writer->Write();

    int* sz = mainWindow_->renderWindow_->GetSize();
    QJsonObject result;
    result["path"] = path;
    result["width"] = sz[0] * scale;
    result["height"] = sz[1] * scale;
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
    const char* names[] = {"bs", "hm", "mc", "ld", "ce", "pq", "rsa", "hb"};
    for (int i = 0; i < 8; i++) {
        if (cmd.contains(names[i])) {
            mainWindow_->physicsChecks_[i]->setChecked(cmd[names[i]].toBool());
        }
    }
    mainWindow_->onPhysicsCheckChanged();
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
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
    mainWindow_->onGlyphScaleChanged(sliderVal);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetOpacity(const QJsonObject& cmd) {
    double opacity = cmd["value"].toDouble(0.7);
    int sliderVal = static_cast<int>(opacity * 100);
    mainWindow_->opacitySlider_->setValue(sliderVal);
    mainWindow_->onOpacityChanged(sliderVal);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}

QJsonObject RestServer::cmdSetIsoThreshold(const QJsonObject& cmd) {
    double threshold = cmd["value"].toDouble(0.5);
    mainWindow_->isoThreshold_->setValue(threshold);
    mainWindow_->onIsoThresholdChanged(threshold);
    mainWindow_->renderWindow_->Render();
    return QJsonObject{{"ok", true}};
}
