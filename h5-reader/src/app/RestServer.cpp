#include "RestServer.h"

#include "AngleCollarActor.h"
#include "CameraAnchorHelper.h"
#include "CameraComposer.h"
#include "DashboardDisplayController.h"
#include "DashboardSelectionController.h"
#include "MeasurementOverlay.h"
#include "MoleculeScene.h"
#include "NewmanProjection.h"
#include "QtPlaybackController.h"
#include "ReaderMainWindow.h"
#include "SignalDisplayDialog.h"
#include "TensorGhostTrail.h"

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
#include "../model/DisplayPolicy.h"
#include "../model/MetricTaxonomy.h"
#include "../model/QtProtein.h"
#include "../model/TrajectoryConformation.h"
#include "../model/VisualizationRegistry.h"
#include "../model/TrajectorySignalCatalog.h"
#include "../model/TransformedConformation.h"

#include <QBuffer>
#include <QByteArray>
#include <QCoreApplication>
#include <QHttpServer>
#include <QHttpServerRequest>
#include <QHttpServerResponse>
#include <QJsonArray>
#include <QUrlQuery>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QLoggingCategory>
#include <QPixmap>
#include <QString>
#include <QTcpServer>
#include <QTcpSocket>
#include <QUuid>
#include <QVariant>
#include <QApplication>
#include <QWidget>

#include <vtkCamera.h>
#include <vtkMoleculeMapper.h>
#include <vtkPNGWriter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkWindowToImageFilter.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cRest, "h5reader.rest")

// A QTcpServer that reports each accepted socket as QAbstractHttpServer pulls it,
// so /shutdown can hook the live socket's flush event and exit on the real
// "response has left the wire" moment instead of a timer.
class TrackingTcpServer : public QTcpServer {
public:
    explicit TrackingTcpServer(QObject* parent = nullptr) : QTcpServer(parent) {}
    std::function<void(QTcpSocket*)> onConnection;
    QTcpSocket* nextPendingConnection() override {
        QTcpSocket* socket = QTcpServer::nextPendingConnection();
        if (socket && onConnection)
            onConnection(socket);
        return socket;
    }
};

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

QJsonArray vec3FromEigen(const model::Vec3& v) {
    return QJsonArray{v.x(), v.y(), v.z()};
}

std::optional<std::size_t> residueIndexFromRequest(const model::QtProtein& protein,
                                                   const QJsonObject& body,
                                                   const model::AtomSelection* selection) {
    if (body.contains(QStringLiteral("residue"))) {
        const int r = body.value(QStringLiteral("residue")).toInt(-1);
        if (r >= 0 && static_cast<std::size_t>(r) < protein.residueCount())
            return static_cast<std::size_t>(r);
        return std::nullopt;
    }
    if (body.contains(QStringLiteral("residue_number"))) {
        const int residueNumber = body.value(QStringLiteral("residue_number")).toInt(-1);
        for (std::size_t i = 0; i < protein.residueCount(); ++i) {
            if (protein.residue(i).address.residueNumber == residueNumber)
                return i;
        }
        return std::nullopt;
    }
    if (selection && selection->hasFocus()) {
        const std::size_t atom = selection->focus();
        if (atom < protein.atomCount()) {
            const int r = protein.atom(atom).residueIndex;
            if (r >= 0 && static_cast<std::size_t>(r) < protein.residueCount())
                return static_cast<std::size_t>(r);
        }
    }
    return std::nullopt;
}

std::optional<std::size_t> findResidueAtomByName(const model::QtProtein& protein,
                                                 std::size_t residue,
                                                 const QString& name) {
    if (residue >= protein.residueCount())
        return std::nullopt;
    const model::QtResidue& r = protein.residue(residue);
    for (int32_t atomIndex : r.atomIndices) {
        if (atomIndex < 0)
            continue;
        const std::size_t a = static_cast<std::size_t>(atomIndex);
        if (a >= protein.atomCount())
            continue;
        const auto& names = protein.atomNames(a);
        if (names.amber == name || names.iupac == name || names.bmrb == name)
            return a;
    }
    return std::nullopt;
}

QString rotamerState(double radians) {
    constexpr double kPiLocal = 3.141592653589793238462643383279502884;
    if (!std::isfinite(radians))
        return QStringLiteral("unknown");
    const double deg = (std::remainder(radians, 2.0 * kPiLocal)) * 180.0 / kPiLocal;
    if (deg < -120.0 || deg >= 120.0)
        return QStringLiteral("trans");
    if (deg < 0.0)
        return QStringLiteral("gminus");
    return QStringLiteral("gplus");
}

double radiansToDegrees(double radians) {
    constexpr double kRadToDeg = 57.295779513082320876798;
    return radians * kRadToDeg;
}

unsigned char byteColorComponent(const QJsonValue& value, unsigned char fallback) {
    if (!value.isDouble())
        return fallback;
    const double raw = value.toDouble();
    const double scaled = raw <= 1.0 ? raw * 255.0 : raw;
    return static_cast<unsigned char>(std::clamp(std::lround(scaled), 0L, 255L));
}

std::array<unsigned char, 3> colorFromJson(const QJsonValue& value,
                                           std::array<unsigned char, 3> fallback) {
    const QJsonArray arr = value.toArray();
    if (arr.size() != 3)
        return fallback;
    return {byteColorComponent(arr.at(0), fallback[0]),
            byteColorComponent(arr.at(1), fallback[1]),
            byteColorComponent(arr.at(2), fallback[2])};
}

QJsonArray colorToJson(const std::array<unsigned char, 3>& color) {
    return QJsonArray{static_cast<int>(color[0]), static_cast<int>(color[1]),
                      static_cast<int>(color[2])};
}

QJsonObject moleculeStyleToJson(const MoleculeScene::MoleculeStyle& style) {
    return QJsonObject{
        {"render_atoms", style.renderAtoms},
        {"render_bonds", style.renderBonds},
        {"use_multi_bonds", style.useMultiCylindersForBonds},
        {"atomic_radius_type", style.atomicRadiusType},
        {"atom_color_mode", style.atomColorMode},
        {"bond_color_mode", style.bondColorMode},
        {"atom_radius_scale", style.atomicRadiusScaleFactor},
        {"bond_radius", style.bondRadius},
        {"atom_color", colorToJson(style.atomColor)},
        {"bond_color", colorToJson(style.bondColor)},
    };
}

void applyMoleculeStylePreset(MoleculeScene::MoleculeStyle& style, const QString& preset) {
    if (preset == QStringLiteral("sticks")) {
        style.renderAtoms = false;
        style.renderBonds = true;
        style.useMultiCylindersForBonds = false;
        style.bondColorMode = vtkMoleculeMapper::SingleColor;
        style.bondRadius = 0.035f;
        style.bondColor = {176, 180, 184};
        return;
    }
    if (preset == QStringLiteral("scaffold") || preset.isEmpty()) {
        style.renderAtoms = true;
        style.renderBonds = true;
        style.useMultiCylindersForBonds = false;
        style.atomicRadiusType = vtkMoleculeMapper::UnitRadius;
        style.atomColorMode = vtkMoleculeMapper::SingleColor;
        style.bondColorMode = vtkMoleculeMapper::SingleColor;
        style.atomicRadiusScaleFactor = 0.11f;
        style.bondRadius = 0.030f;
        style.atomColor = {222, 224, 226};
        style.bondColor = {170, 174, 178};
    }
}

int radiusTypeFromString(const QString& value, int fallback) {
    if (value == QStringLiteral("covalent"))
        return vtkMoleculeMapper::CovalentRadius;
    if (value == QStringLiteral("vdw"))
        return vtkMoleculeMapper::VDWRadius;
    if (value == QStringLiteral("unit"))
        return vtkMoleculeMapper::UnitRadius;
    return fallback;
}

int colorModeFromString(const QString& value, int fallback) {
    if (value == QStringLiteral("single"))
        return vtkMoleculeMapper::SingleColor;
    if (value == QStringLiteral("element") || value == QStringLiteral("discrete"))
        return vtkMoleculeMapper::DiscreteByAtom;
    return fallback;
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

QJsonArray stringListToJson(const QStringList& values) {
    QJsonArray out;
    for (const QString& value : values)
        out.append(value);
    return out;
}

QString uuidToString(const QUuid& id) {
    return id.toString(QUuid::WithoutBraces);
}

QJsonObject stripTrackToJson(const DashboardDisplayController::StripTrack& track) {
    QJsonArray values;
    QJsonArray valid;
    const model::ChannelBuffer* buffer = track.buffer;

    if (buffer) {
        for (const double value : buffer->values) {
            values.append(std::isfinite(value) ? QJsonValue(value) : QJsonValue(QJsonValue::Null));
        }
        for (const std::uint8_t flag : buffer->valid) {
            valid.append(flag ? 1 : 0);
        }
    }

    return QJsonObject{
        {"signal_id", uuidToString(track.signalId)},
        {"descriptor_id", track.descriptorId},
        {"mode", track.displayModeId},
        {"channel", track.channelId},
        {"label", track.label},
        {"last_frame", buffer ? static_cast<qint64>(buffer->lastFrame()) : static_cast<qint64>(-1)},
        {"count", buffer ? static_cast<qint64>(buffer->size()) : static_cast<qint64>(0)},
        {"values", values},
        {"valid", valid},
    };
}

std::optional<QUuid> uuidFromBody(const QJsonObject& body, QString* error) {
    const QString raw = body.value("id").toString().trimmed();
    if (raw.isEmpty()) {
        *error = QStringLiteral("body must include non-empty string field \"id\"");
        return std::nullopt;
    }
    QUuid id(raw);
    if (id.isNull() && !raw.startsWith(QLatin1Char('{')))
        id = QUuid(QStringLiteral("{%1}").arg(raw));
    if (id.isNull()) {
        *error = QStringLiteral("invalid uuid: %1").arg(raw);
        return std::nullopt;
    }
    return id;
}

bool readNonNegativeInteger(const QJsonObject& obj,
                            const QString& key,
                            qint64* out,
                            QString* error,
                            bool required = true) {
    if (!obj.contains(key)) {
        if (required)
            *error = QStringLiteral("missing integer field \"%1\"").arg(key);
        return !required;
    }
    const QJsonValue value = obj.value(key);
    if (!value.isDouble()) {
        *error = QStringLiteral("field \"%1\" must be an integer").arg(key);
        return false;
    }
    const qint64 raw = value.toInteger(-1);
    if (raw < 0) {
        *error = QStringLiteral("field \"%1\" must be >= 0").arg(key);
        return false;
    }
    *out = raw;
    return true;
}

bool validateOptionalFrame(const QJsonObject& anchor,
                           const io::QtLoadResult* loaded,
                           QString* error) {
    if (!anchor.contains(QStringLiteral("frame")))
        return true;
    qint64 frame = -1;
    if (!readNonNegativeInteger(anchor, QStringLiteral("frame"), &frame, error))
        return false;
    const auto* conformation = loaded ? loaded->conformation.get() : nullptr;
    if (conformation && static_cast<std::size_t>(frame) >= conformation->frameCount()) {
        *error = QStringLiteral("anchor frame out of range: %1").arg(frame);
        return false;
    }
    return true;
}

bool indexInRange(qint64 index, std::size_t count, const QString& label, QString* error) {
    if (index < 0 || static_cast<std::size_t>(index) >= count) {
        *error = QStringLiteral("%1 index out of range: %2").arg(label).arg(index);
        return false;
    }
    return true;
}

struct ParsedDashboardAnchor {
    model::SignalAnchor anchor = model::NoneAnchor{};
    bool followsFocus = false;
    QString label = QStringLiteral("No anchor");
};

bool parseDashboardAnchor(const QJsonObject& body,
                          const io::QtLoadResult* loaded,
                          ParsedDashboardAnchor* out,
                          QString* error) {
    const QJsonValue anchorValue = body.value(QStringLiteral("anchor"));
    if (!anchorValue.isObject()) {
        *error = QStringLiteral("body must include object field \"anchor\"");
        return false;
    }
    const QJsonObject anchor = anchorValue.toObject();
    if (anchor.value(QStringLiteral("follows_focus")).toBool(false)) {
        out->anchor = model::NoneAnchor{};
        out->followsFocus = true;
        out->label = QStringLiteral("focus");
        return true;
    }
    if (!validateOptionalFrame(anchor, loaded, error))
        return false;

    const auto* protein = loaded ? loaded->protein.get() : nullptr;
    const QString kind = anchor.value(QStringLiteral("kind")).toString();
    auto readIndex = [&](const QString& key, qint64* value) {
        return readNonNegativeInteger(anchor, key, value, error);
    };

    qint64 raw = -1;
    if (kind == QStringLiteral("none")) {
        out->anchor = model::NoneAnchor{};
    } else if (kind == QStringLiteral("protein")) {
        out->anchor = model::ProteinAnchor{};
    } else if (kind == QStringLiteral("system")) {
        out->anchor = model::SystemAnchor{};
    } else if (kind == QStringLiteral("event")) {
        out->anchor = model::EventAnchor{};
    } else if (kind == QStringLiteral("atom") || (kind.isEmpty() && anchor.contains(QStringLiteral("atom")))) {
        if (!readIndex(QStringLiteral("atom"), &raw))
            return false;
        if (protein && !indexInRange(raw, protein->atomCount(), QStringLiteral("atom"), error))
            return false;
        out->anchor = model::AtomAnchor{static_cast<std::size_t>(raw)};
    } else if (kind == QStringLiteral("residue") || (kind.isEmpty() && anchor.contains(QStringLiteral("residue")))) {
        if (!readIndex(QStringLiteral("residue"), &raw))
            return false;
        if (protein && !indexInRange(raw, protein->residueCount(), QStringLiteral("residue"), error))
            return false;
        out->anchor = model::ResidueAnchor{static_cast<std::size_t>(raw)};
    } else if (kind == QStringLiteral("bond") || (kind.isEmpty() && anchor.contains(QStringLiteral("bond")))) {
        if (!readIndex(QStringLiteral("bond"), &raw))
            return false;
        if (protein && !indexInRange(raw, protein->bondCount(), QStringLiteral("bond"), error))
            return false;
        out->anchor = model::BondAnchor{static_cast<std::size_t>(raw)};
    } else if (kind == QStringLiteral("bond_vector")
               || (kind.isEmpty() && anchor.contains(QStringLiteral("kind_id")))) {
        qint64 residue = -1;
        qint64 kindId = 0;
        if (!readIndex(QStringLiteral("residue"), &residue))
            return false;
        if (anchor.contains(QStringLiteral("kind_id"))) {
            if (!readNonNegativeInteger(anchor, QStringLiteral("kind_id"), &kindId, error))
                return false;
        }
        if (protein && !indexInRange(residue, protein->residueCount(), QStringLiteral("residue"), error))
            return false;
        if (kindId > 255) {
            *error = QStringLiteral("kind_id must be <= 255");
            return false;
        }
        out->anchor = model::BondVectorAnchor{static_cast<std::size_t>(residue),
                                              static_cast<std::uint8_t>(kindId)};
    } else if (kind == QStringLiteral("ring") || (kind.isEmpty() && anchor.contains(QStringLiteral("ring")))) {
        if (!readIndex(QStringLiteral("ring"), &raw))
            return false;
        if (protein && !indexInRange(raw, protein->ringCount(), QStringLiteral("ring"), error))
            return false;
        out->anchor = model::RingAnchor{static_cast<std::size_t>(raw)};
    } else if (kind == QStringLiteral("aromatic_ring")) {
        if (!readIndex(QStringLiteral("ring"), &raw))
            return false;
        if (protein && !indexInRange(raw, protein->ringCount(), QStringLiteral("ring"), error))
            return false;
        out->anchor = model::AromaticRingAnchor{static_cast<std::size_t>(raw)};
    } else if (kind == QStringLiteral("saturated_ring")) {
        if (!readIndex(QStringLiteral("ring"), &raw))
            return false;
        if (protein && !indexInRange(raw, protein->ringCount(), QStringLiteral("ring"), error))
            return false;
        out->anchor = model::SaturatedRingAnchor{static_cast<std::size_t>(raw)};
    } else if (kind == QStringLiteral("ring_contribution_pair")
               || (kind.isEmpty() && anchor.contains(QStringLiteral("pair")))) {
        if (!readIndex(QStringLiteral("pair"), &raw))
            return false;
        out->anchor = model::RingContributionPairAnchor{static_cast<std::size_t>(raw)};
    } else if (kind == QStringLiteral("ring_membership")
               || (kind.isEmpty() && anchor.contains(QStringLiteral("membership")))) {
        if (!readIndex(QStringLiteral("membership"), &raw))
            return false;
        out->anchor = model::RingMembershipAnchor{static_cast<std::size_t>(raw)};
    } else if (kind == QStringLiteral("mutation_match_pair")) {
        if (!readIndex(QStringLiteral("pair"), &raw))
            return false;
        out->anchor = model::MutationMatchPairAnchor{static_cast<std::size_t>(raw)};
    } else {
        *error = QStringLiteral("unsupported anchor object");
        return false;
    }

    out->followsFocus = false;
    out->label = model::AnchorLabel(out->anchor);
    return true;
}

bool parseModeArray(const QJsonObject& body, QStringList* modes, QString* error) {
    const QJsonValue modesValue = body.value(QStringLiteral("modes"));
    if (!modesValue.isArray()) {
        *error = QStringLiteral("body must include array field \"modes\"");
        return false;
    }
    for (const QJsonValue& value : modesValue.toArray()) {
        const QString mode = value.toString().trimmed();
        if (mode.isEmpty()) {
            *error = QStringLiteral("modes must contain non-empty strings");
            return false;
        }
        if (!modes->contains(mode))
            modes->push_back(mode);
    }
    if (modes->isEmpty()) {
        *error = QStringLiteral("modes must contain at least one mode");
        return false;
    }
    return true;
}

const model::SignalDescriptor* findAnyDescriptor(const model::TrajectorySignalCatalog* catalog,
                                                 const QString& descriptorId) {
    if (!catalog)
        return nullptr;
    if (const auto* descriptor = catalog->findDescriptor(descriptorId))
        return descriptor;
    for (const model::SignalDescriptor& descriptor : catalog->allDescriptorList()) {
        if (descriptor.id == descriptorId)
            return &descriptor;
    }
    return nullptr;
}

QString availabilityString(const model::TrajectorySignalCatalog* catalog,
                           const QString& descriptorId) {
    if (!catalog || !catalog->fieldAvailability())
        return QStringLiteral("unknown");
    return QString::fromLatin1(model::ToString(catalog->fieldAvailability()->stateForDescriptor(descriptorId)));
}

QJsonArray modeStateArray(const model::DashboardSignal& signal,
                          const model::SignalDescriptor* descriptor) {
    QStringList modes = descriptor ? model::AllDisplayModes(*descriptor) : QStringList{};
    for (const QString& enabledMode : signal.displayModeIds) {
        if (!modes.contains(enabledMode))
            modes.push_back(enabledMode);
    }

    QJsonArray out;
    for (const QString& mode : modes) {
        const model::DashboardSignalModel::ModeRenderability renderability =
            model::DashboardSignalModel::ModeRenderabilityFor(mode);
        out.append(QJsonObject{
            {"mode", mode},
            {"enabled", signal.displayModeIds.contains(mode)},
            {"is_panel_ref", renderability.emitsPanelRef},
            {"renderable_panel", renderability.buildsPanelWidget},
            {"has_visible_surface", renderability.hasVisibleSurface},
        });
    }
    return out;
}

QJsonArray selectedStateArray(const model::DashboardSignalModel* signalModel,
                              const model::TrajectorySignalCatalog* catalog) {
    QJsonArray selected;
    if (!signalModel)
        return selected;
    const QVector<model::DashboardSignal>& selectedSignals = signalModel->activeSignals();
    for (int row = 0; row < selectedSignals.size(); ++row) {
        const model::DashboardSignal& signal = selectedSignals.at(row);
        const model::SignalDescriptor* descriptor = findAnyDescriptor(catalog, signal.binding.descriptorId);
        selected.append(QJsonObject{
            {"id", uuidToString(signal.id)},
            {"descriptor_id", signal.binding.descriptorId},
            {"concept_key", signal.binding.conceptKey},
            {"label", signal.label},
            {"anchor", anchorToJson(signal.binding.anchor)},
            {"enabled", signal.enabled},
            {"availability", signalModel->availabilityName(row)},
            {"modes", modeStateArray(signal, descriptor)},
        });
    }
    return selected;
}

QJsonArray expectedEmptyArray(const DashboardSmokeSummary& summary) {
    QJsonArray out;
    for (const DashboardSmokeSummary::ExpectedButEmpty& record : summary.expectedButEmpty) {
        out.append(QJsonObject{
            {"descriptor_id", record.descriptorId},
            {"storage_path", record.storagePath},
            {"visualization_type", record.visualizationType},
            {"canonical_state", record.canonicalState},
            {"storage_path_state", record.storagePathState},
        });
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
QByteArray captureScenePng(MoleculeScene* scene, bool forceRender = true, int scale = 1) {
    if (!scene || !scene->Renderer() || !scene->Renderer()->GetRenderWindow())
        return {};
    ASSERT_THREAD(scene);

    auto w2i = vtkSmartPointer<vtkWindowToImageFilter>::New();
    w2i->SetInput(scene->Renderer()->GetRenderWindow());
    w2i->SetInputBufferTypeToRGB();
    w2i->ReadFrontBufferOff();
    if (scale > 1) {
        // Poster / print export: render the scene at scale x scale tiles, then
        // stitch. The extra pixels are effective supersampling (SSAA) -- edges
        // smoother than the live FXAA, and a high-DPI image fit for print.
        // FixBoundary removes the seams between tiles. Scale > 1 always forces
        // a re-render (the tiles must be drawn), so ShouldRerender is moot here.
        w2i->SetScale(scale);
        w2i->FixBoundaryOn();
    } else if (!forceRender) {
        w2i->ShouldRerenderOff();
    }
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

// Find a live widget by objectName for targeted snapshots — top-level
// widgets (e.g. the modeless Metric Picker dialog) first, then descendants
// of the main window (docks, strip panels). Returns nullptr if none match.
QWidget* findWidgetByObjectName(QWidget* mainWindow, const QString& name) {
    if (name.isEmpty())
        return nullptr;
    const QWidgetList topLevels = QApplication::topLevelWidgets();
    for (QWidget* w : topLevels) {
        if (w && w->objectName() == name)
            return w;
    }
    if (mainWindow) {
        if (QWidget* child = mainWindow->findChild<QWidget*>(name))
            return child;
    }
    return nullptr;
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
                            DashboardSelectionController* selectionController,
                            DashboardDisplayController* dashboardController,
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
    selectionController_ = selectionController;
    dashboardController_ = dashboardController;
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

#if QT_VERSION >= QT_VERSION_CHECK(6, 10, 0)
    // Qt 6.10 removed the QHttpServer::listen() convenience overload. The
    // supported path is now QAbstractHttpServer::bind() with a QTcpServer that
    // is already listening. The Linux dev boxes still run Qt 6.4, so the old
    // listen() path is kept under the version guard below.
    auto* tcp = new TrackingTcpServer(this);
    tcp->onConnection = [this](QTcpSocket* socket) { activeSocket_ = socket; };
    if (!tcp->listen(QHostAddress::LocalHost, port) || !server_->bind(tcp)) {
        qCCritical(cRest).noquote()
            << "REST server failed to bind 127.0.0.1 port" << port;
        delete tcp;
        server_.reset();
        return 0;
    }
    const quint16 bound = tcp->serverPort();
#else
    const quint16 bound = server_->listen(QHostAddress::LocalHost, port);
    if (bound == 0) {
        qCCritical(cRest).noquote()
            << "REST server failed to bind 127.0.0.1 port" << port;
        server_.reset();
        return 0;
    }
#endif

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

    // ---- inspector panel readout ----------------------------------------

    // GET /inspector/tree -> the focused atom's Atom Info panel serialized as
    // [{field, value, tooltip?, children?}], so the curated primary set and its
    // provenance/status tooltips are REST-assertable (tooltips don't screenshot).
    server_->route(QStringLiteral("/inspector/tree"), [this]() {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader window not wired"), SC::ServiceUnavailable);
        return jsonResponse(readerWindow_->inspectorTreeJson());
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

    // GET /newman -> the focused residue's backbone phi/psi Newman projections
    // (torsion angle + substituent spokes) at the current frame. Drives and
    // verifies the dashboard Newman display; pure geometry, transform-invariant
    // so it reads the base conformation.
    server_->route(QStringLiteral("/newman"), [this]() {
        ASSERT_THREAD(this);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        const auto* conf    = loaded_ ? loaded_->conformation.get() : nullptr;
        if (!protein || !conf)
            return errorResponse(QStringLiteral("no protein loaded"), SC::ServiceUnavailable);
        if (!selection_ || !selection_->hasFocus())
            return errorResponse(QStringLiteral("no focused atom"), SC::ServiceUnavailable);
        const std::size_t atom = selection_->focus();
        if (atom >= protein->atomCount())
            return errorResponse(QStringLiteral("focus atom out of range"), SC::BadRequest);
        const int rRaw = protein->atom(atom).residueIndex;
        if (rRaw < 0 || static_cast<std::size_t>(rRaw) >= protein->residueCount())
            return errorResponse(QStringLiteral("focus atom has no residue"), SC::ServiceUnavailable);
        const std::size_t residue = static_cast<std::size_t>(rRaw);
        const int frame = playback_ ? playback_->currentFrame() : 0;
        const std::size_t f = (frame >= 0) ? static_cast<std::size_t>(frame) : 0u;

        const double kRadToDeg = 57.295779513082320876;
        auto encode = [&](const NewmanProjection& p) {
            QJsonObject o;
            o["valid"] = p.valid;
            if (!p.valid) { o["reason"] = p.invalidReason; return o; }
            o["angleDeg"] = p.torsionDeg;
            o["front"]    = p.frontLabel;
            o["back"]     = p.backLabel;
            QJsonArray spokes;
            for (const NewmanSpoke& s : p.spokes)
                spokes.append(QJsonObject{
                    {"label", s.label},
                    {"angleDeg", s.angleRad * kRadToDeg},
                    {"front", s.front},
                    {"reference", s.reference},
                });
            o["spokes"] = spokes;
            return o;
        };

        const NewmanProjection phiP = ComputeNewmanProjection(*protein, *conf, residue, f, NewmanKind::Phi);
        const NewmanProjection psiP = ComputeNewmanProjection(*protein, *conf, residue, f, NewmanKind::Psi);
        QJsonObject out;
        out["residue"]      = static_cast<qint64>(residue);
        out["residueLabel"] = phiP.residueLabel;
        out["frame"]        = frame;
        out["phi"]          = encode(phiP);
        out["psi"]          = encode(psiP);
        return jsonResponse(out);
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

    // ---- overlay visibility (automation / snapshot harness) -------------
    //
    // POST /overlay {"name": "ribbon"|"rings"|"butterfly"|"bfield"|"shadow",
    //                "visible": bool}
    // Drives the same toolbar toggle path a human click would, so the kernel
    // overlays (butterfly isosurfaces, B-field streamlines) get their per-frame
    // refresh and the toolbar checkbox stays in sync. Lets the headless
    // snapshot loop enable the field overlays, which default off and are not
    // persisted in QSettings.
    server_->route(QStringLiteral("/overlay"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("visible") || !body.value("visible").isBool())
            return errorResponse(
                QStringLiteral("body must be {\"name\": str, \"visible\": bool}"),
                SC::BadRequest);
        const QString name = body.value("name").toString();
        const bool visible = body.value("visible").toBool();
        if (!readerWindow_->setOverlayVisible(name, visible))
            return errorResponse(
                QStringLiteral("unknown overlay \"%1\" "
                               "(ribbon|rings|butterfly|bfield|shadow)").arg(name),
                SC::BadRequest);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- field-grid isosurface threshold (live tuning) ------------------
    //
    // POST /field/threshold {"ppm": number >= 0}
    // Sets the butterfly (field-grid) |T0| isosurface contour level in ppm and
    // re-renders. Lets the snapshot harness sweep the dominant-zone threshold
    // (default 0.30) without a rebuild while we tune the picture.
    server_->route(QStringLiteral("/field/threshold"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("ppm") || !body.value("ppm").isDouble())
            return errorResponse(QStringLiteral("body must be {\"ppm\": number}"),
                                 SC::BadRequest);
        const double ppm = body.value("ppm").toDouble();
        if (ppm < 0.0)
            return errorResponse(QStringLiteral("ppm must be >= 0"), SC::BadRequest);
        if (!readerWindow_->setFieldThreshold(ppm))
            return errorResponse(QStringLiteral("field-grid overlay not available"),
                                 SC::ServiceUnavailable);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- butterfly Gaussian taper: extent σ + peak amplitude ------------
    //
    // POST /field/gaussian {"extent": σ_Å >= 0.1 (optional),
    //                       "peak":   amplitude >= 0 (optional)}
    // The two orthogonal butterfly knobs: extent = radial Gaussian σ (lobe
    // reach), peak = amplitude gain (1.0 = true near-ring physics). Either or
    // both per call, so the tuning loop can sweep one at a time. (Contour level
    // is the separate /field/threshold knob.)
    server_->route(QStringLiteral("/field/gaussian"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);
        bool applied = false;
        if (body.value("extent").isDouble()) {
            const double sigma = body.value("extent").toDouble();
            if (sigma < 0.1)
                return errorResponse(QStringLiteral("extent (σ) must be >= 0.1"),
                                     SC::BadRequest);
            if (!readerWindow_->setFieldExtent(sigma))
                return errorResponse(QStringLiteral("field-grid overlay not available"),
                                     SC::ServiceUnavailable);
            applied = true;
        }
        if (body.value("peak").isDouble()) {
            const double peak = body.value("peak").toDouble();
            if (peak < 0.0)
                return errorResponse(QStringLiteral("peak must be >= 0"), SC::BadRequest);
            if (!readerWindow_->setFieldPeak(peak))
                return errorResponse(QStringLiteral("field-grid overlay not available"),
                                     SC::ServiceUnavailable);
            applied = true;
        }
        if (!applied)
            return errorResponse(
                QStringLiteral("body must set \"extent\" and/or \"peak\" (numbers)"),
                SC::BadRequest);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- UI state introspection -----------------------------------------
    //
    // GET /ui/state → full operating-state + per-control {enabled, checked}
    // snapshot. A real introspection API (a client can mirror the toolbar from
    // it), and how the selectability rules get asserted across states.
    server_->route(QStringLiteral("/ui/state"), Method::Get,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"),
                                 SC::ServiceUnavailable);
        return jsonResponse(readerWindow_->uiStateJson());
    });

    // ---- playback control (real transport API; mirrors the toolbar) ------
    //
    // POST /playback {"action": "play_forward"|"play_backward"|"pause"|
    //                 "step_forward"|"step_backward"|"toggle"}
    server_->route(QStringLiteral("/playback"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!playback_)
            return errorResponse(QStringLiteral("no playback (nothing loaded)"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        const QString a = (ok && body.contains("action"))
                              ? body.value("action").toString() : QString();
        if      (a == QStringLiteral("play_forward"))  playback_->playForward();
        else if (a == QStringLiteral("play_backward")) playback_->playBackward();
        else if (a == QStringLiteral("pause") ||
                 a == QStringLiteral("stop"))          playback_->pause();
        else if (a == QStringLiteral("step_forward"))  playback_->stepForward();
        else if (a == QStringLiteral("step_backward")) playback_->stepBackward();
        else if (a == QStringLiteral("toggle"))        playback_->togglePlayPause();
        else
            return errorResponse(
                QStringLiteral("action must be play_forward|play_backward|pause|"
                               "step_forward|step_backward|toggle"),
                SC::BadRequest);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- display isolation (filter mode) --------------------------------
    //
    // POST /filter {"residues": [int, ...]} — show only those residues' atoms in
    // the 3-D view; an empty/absent array restores the full structure. The
    // atom + nearby-residue view the camera modes used to approximate.
    server_->route(QStringLiteral("/filter"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);
        std::vector<std::size_t> residues;
        const QJsonArray arr = body.value(QStringLiteral("residues")).toArray();
        for (const QJsonValue v : arr)
            if (v.isDouble()) residues.push_back(static_cast<std::size_t>(v.toInt()));
        readerWindow_->setResidueFilter(residues);   // empty → restore full
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- transform (upstream data-layer transform — TransformedConformation) -
    //
    // POST /transform {"kind": "all_atom_fit"|"backbone_fit",
    //                   "reference_frame": int (default 0; mean seed + centroid anchor),
    //                   "subset_atoms": [int, ...] (fit_subset alias only),
    //                   "backbone_only": bool (fit_subset alias shorthand) }
    // Switches the wrapped Conformation's transform mode. Fire-and-forget;
    // the wrapper emits transformChanged() which is connected (in
    // ReaderMainWindow) to scene_->refreshCurrentFrame so the molecule
    // re-renders in the new frame without further client involvement.
    //
    // POST /transform/smoothing {"window": int} changes the symmetric
    // temporal rotation-smoothing half-width. window=0 disables smoothing;
    // translation is always re-derived from the current fit centroid.
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
            case model::TransformedConformation::Mode::FitReference: kind = QStringLiteral("all_atom_fit"); break;
            case model::TransformedConformation::Mode::FitSubset:    kind = QStringLiteral("backbone_fit"); break;
        }
        QJsonArray subsetArr;
        for (std::size_t a : transformed_->subsetAtoms())
            subsetArr.append(static_cast<qint64>(a));
        return jsonResponse(QJsonObject{
            {"kind", kind},
            {"reference_frame", static_cast<qint64>(transformed_->referenceFrame())},
            {"subset_atoms", subsetArr},
            {"subset_size", subsetArr.size()},
            {"window", transformed_->stabilisationWindow()},
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
        const QString kindStr = body.value("kind").toString().toLower();
        const std::size_t referenceFrame = body.contains("reference_frame")
            ? static_cast<std::size_t>(body.value("reference_frame").toInteger(0))
            : 0;

        using Mode = model::TransformedConformation::Mode;
        Mode mode = Mode::FitSubset;
        std::vector<std::size_t> subset;

        if (kindStr == QStringLiteral("all_atom_fit")
            || kindStr == QStringLiteral("fit_reference")) {
            mode = Mode::FitReference;
        } else if (kindStr == QStringLiteral("backbone_fit")) {
            mode = Mode::FitSubset;
            const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
            if (!protein)
                return errorResponse(QStringLiteral("no protein loaded for backbone_fit"),
                                     SC::ServiceUnavailable);
            subset = model::TransformedConformation::BackboneSubset(*protein);
            if (subset.size() < 3)
                return errorResponse(QStringLiteral("backbone subset has <3 atoms — "
                                                    "fit underdetermined"),
                                     SC::Conflict);
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
            return errorResponse(QStringLiteral("unknown transform kind: %1 "
                                                "(expected all_atom_fit or backbone_fit)")
                                     .arg(kindStr),
                                 SC::BadRequest);
        }

        transformed_->setMode(mode, referenceFrame, std::move(subset));
        // setMode emits transformChanged → ReaderMainWindow connects this
        // to scene_->refreshCurrentFrame. No explicit render here; the
        // connected slot handles it.
        return QHttpServerResponse(SC::NoContent);
    });

    server_->route(QStringLiteral("/transform/smoothing"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!transformed_)
            return errorResponse(QStringLiteral("transformed conformation not wired"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        const QJsonValue value = body.value(QStringLiteral("window"));
        if (!ok || !value.isDouble())
            return errorResponse(QStringLiteral("body must be {\"window\": int}"),
                                 SC::BadRequest);

        const double asDouble = value.toDouble(-1.0);
        if (!std::isfinite(asDouble)
            || std::floor(asDouble) != asDouble
            || asDouble < 0.0
            || asDouble > static_cast<double>(std::numeric_limits<int>::max())) {
            return errorResponse(QStringLiteral("window must be a non-negative integer"),
                                 SC::BadRequest);
        }

        transformed_->setStabilisationWindow(static_cast<int>(asDouble));
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
        scene_->syncCameraClippingRange();
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
        scene_->syncCameraClippingRange();
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
        scene_->syncCameraClippingRange();
        scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- camera inspect atom (free-camera framing to judge ONE atom) -----
    //
    // POST /camera/inspect_atom {"atom": int, "distance"?: double,
    //                            "azimuth"?: double, "elevation"?: double,
    //                            "roll"?: double}
    //
    // Frees the camera (so this absolute write STICKS -- a Free-mode per-frame
    // write is a no-op) and aims it at the atom's current display position from
    // `distance` Angstrom (default 12) along the current view direction,
    // optionally orbited by azimuth/elevation and rolled. Frames ONE atom close
    // enough, with its neighbourhood in view, so its orientation and its
    // appearance relative to nearby atoms can be judged. Returns the resulting
    // camera + the framed point. Pair with /selection/pick to highlight it and
    // (for tensor work) drive the CSA glyph.
    server_->route(QStringLiteral("/camera/inspect_atom"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->cameraComposer() || !scene_->Renderer())
            return errorResponse(QStringLiteral("scene not wired"), SC::ServiceUnavailable);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        if (!protein || !transformed_)
            return errorResponse(QStringLiteral("no protein loaded"), SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains("atom"))
            return errorResponse(QStringLiteral(
                "body must be {\"atom\": int, \"distance\"?, \"azimuth\"?, "
                "\"elevation\"?, \"roll\"?}"), SC::BadRequest);
        const qint64 rawAtom = body.value("atom").toInteger(-1);
        if (rawAtom < 0 || static_cast<std::size_t>(rawAtom) >= protein->atomCount())
            return errorResponse(QStringLiteral("atom out of range"), SC::BadRequest);
        const std::size_t atomIdx = static_cast<std::size_t>(rawAtom);

        double distance = body.value("distance").toDouble(12.0);
        if (!(distance > 0.5)) distance = 12.0;
        const double azimuth = body.value("azimuth").toDouble(0.0);
        const double elevation = body.value("elevation").toDouble(0.0);
        const double roll = body.value("roll").toDouble(0.0);

        const std::size_t frame = playback_
            ? static_cast<std::size_t>(playback_->currentFrame()) : 0;
        const model::Vec3 p = transformed_->atomPosition(frame, atomIdx);

        // Free the camera so the per-frame composer write does not clobber this.
        scene_->cameraComposer()->setMode(FreeMode(), DefaultPolicy(), frame);

        auto* cam = scene_->Renderer()->GetActiveCamera();
        if (!cam)
            return errorResponse(QStringLiteral("no active camera"), SC::ServiceUnavailable);
        double dop[3];
        cam->GetDirectionOfProjection(dop);  // unit vector, camera -> focal
        cam->SetFocalPoint(p.x(), p.y(), p.z());
        cam->SetPosition(p.x() - dop[0] * distance,
                         p.y() - dop[1] * distance,
                         p.z() - dop[2] * distance);
        cam->OrthogonalizeViewUp();
        if (azimuth != 0.0) cam->Azimuth(azimuth);
        if (elevation != 0.0) cam->Elevation(elevation);
        cam->OrthogonalizeViewUp();
        if (roll != 0.0) cam->Roll(roll);
        scene_->Renderer()->ResetCameraClippingRange();
        scene_->requestRender(MoleculeScene::RenderSource::Rest);

        double point[3] = {p.x(), p.y(), p.z()};
        double focal[3]{}, position[3]{}, viewUp[3]{};
        cam->GetFocalPoint(focal);
        cam->GetPosition(position);
        cam->GetViewUp(viewUp);
        return jsonResponse(QJsonObject{
            {"atom", static_cast<qint64>(atomIdx)},
            {"frame", static_cast<qint64>(frame)},
            {"distance", distance},
            {"point", vec3FromRaw(point)},
            {"focal", vec3FromRaw(focal)},
            {"position", vec3FromRaw(position)},
            {"view_up", vec3FromRaw(viewUp)},
        });
    });

    // ---- CSA probe (vet the shielding-tensor glyph numerically) ----------
    //
    // GET /csa          -> the focused atom's CSA result.
    // GET /csa?atom=N   -> atom N's CSA result (no selection side-effects).
    //
    // Returns {atom, dft_present, valid, framed, frame_kind, sigma_iso, eta,
    // span, skew, principal_values, pas_axes, molecular_axes}. Computed by the
    // SAME ComputeAtomCsa the glyph uses (via ReaderMainWindow::probeAtomCsa),
    // so the picture and these numbers can never disagree. pas_axes /
    // molecular_axes are 3 director columns (x,y,z) in display coordinates.
    server_->route(QStringLiteral("/csa"), [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader window not wired"), SC::ServiceUnavailable);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        if (!protein)
            return errorResponse(QStringLiteral("no protein loaded"), SC::ServiceUnavailable);

        std::size_t atom = 0;
        const QString atomQ = req.query().queryItemValue(QStringLiteral("atom"));
        if (!atomQ.isEmpty()) {
            bool okNum = false;
            const qint64 a = atomQ.toLongLong(&okNum);
            if (!okNum || a < 0 || static_cast<std::size_t>(a) >= protein->atomCount())
                return errorResponse(QStringLiteral("atom out of range"), SC::BadRequest);
            atom = static_cast<std::size_t>(a);
        } else if (selection_ && selection_->hasFocus()) {
            atom = selection_->focus();
        } else {
            return errorResponse(QStringLiteral("no atom: pass ?atom=N or focus one"),
                                 SC::BadRequest);
        }

        // Optional ?frame=N probes that frame's tensor straight from the DFT
        // store WITHOUT moving the live frame -- the per-frame basis the tensor
        // ghost trail samples. Default (or <0) is the live frame.
        int requestedFrame = -1;
        const QString frameQ = req.query().queryItemValue(QStringLiteral("frame"));
        if (!frameQ.isEmpty()) {
            bool okF = false;
            const qint64 f = frameQ.toLongLong(&okF);
            if (!okF || f < 0)
                return errorResponse(QStringLiteral("frame must be a non-negative integer"), SC::BadRequest);
            requestedFrame = static_cast<int>(f);
        }
        const model::AtomCsaResult r = readerWindow_->probeAtomCsa(atom, requestedFrame);
        const int resolvedFrame = requestedFrame >= 0
            ? requestedFrame : (playback_ ? playback_->currentFrame() : 0);
        QJsonObject out{
            {"atom", static_cast<qint64>(atom)},
            {"frame", resolvedFrame},
            {"dft_present", r.dftPresent},
            {"valid", r.valid},
            {"framed", r.framed},
            {"frame_kind", QString::fromLatin1(model::MolecularFrameKindName(r.frameKind))},
        };
        if (r.valid) {
            out.insert(QStringLiteral("sigma_iso"), r.shape.sigma_iso);
            out.insert(QStringLiteral("eta"), r.shape.eta);
            out.insert(QStringLiteral("span"), r.shape.span);
            out.insert(QStringLiteral("skew"), r.shape.skew);
            out.insert(QStringLiteral("principal_values"),
                       QJsonArray{r.shape.principal_values[0], r.shape.principal_values[1],
                                  r.shape.principal_values[2]});
            QJsonArray pas;
            for (int c = 0; c < 3; ++c) {
                double col[3] = {r.shape.pas_axes(0, c), r.shape.pas_axes(1, c),
                                 r.shape.pas_axes(2, c)};
                pas.append(vec3FromRaw(col));
            }
            out.insert(QStringLiteral("pas_axes"), pas);
            if (r.framed && r.molecularAxes) {
                QJsonArray mol;
                for (int c = 0; c < 3; ++c) {
                    double col[3] = {(*r.molecularAxes)(0, c), (*r.molecularAxes)(1, c),
                                     (*r.molecularAxes)(2, c)};
                    mol.append(vec3FromRaw(col));
                }
                out.insert(QStringLiteral("molecular_axes"), mol);
            }
        }
        return jsonResponse(out);
    });

    // ---- heroshot: transient molecule styling ---------------------------
    //
    // POST /heroshot/molecule_style
    //   {preset:"scaffold"|"sticks", atom_radius_scale?, bond_radius?,
    //    atom_color?, bond_color?, render_atoms?, render_bonds?, ...}
    //
    // Figure work often needs the molecule to recede so tensor/angle geometry
    // can read. This is deliberately transient and restored by /heroshot/clear.
    server_->route(QStringLiteral("/heroshot/molecule_style"), Method::Post,
                   [this](const QHttpServerRequest& request) {
        ASSERT_THREAD(this);
        if (!scene_)
            return errorResponse(QStringLiteral("scene not loaded"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(request, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("expected JSON object"), SC::BadRequest);

        if (!heroshotMoleculeStyleBefore_.has_value())
            heroshotMoleculeStyleBefore_ = scene_->moleculeStyle();
        MoleculeScene::MoleculeStyle style = scene_->moleculeStyle();
        applyMoleculeStylePreset(style, body.value(QStringLiteral("preset")).toString(
                                            QStringLiteral("scaffold")));

        auto styledBool = [&](const QString& key, bool fallback) {
            const QJsonValue v = body.value(key);
            return v.isBool() ? v.toBool() : fallback;
        };
        auto styledDouble = [&](const QString& key, double fallback, double lo, double hi) {
            const QJsonValue v = body.value(key);
            if (!v.isDouble()) return fallback;
            return std::clamp(v.toDouble(), lo, hi);
        };

        style.renderAtoms = styledBool(QStringLiteral("render_atoms"), style.renderAtoms);
        style.renderBonds = styledBool(QStringLiteral("render_bonds"), style.renderBonds);
        style.useMultiCylindersForBonds =
            styledBool(QStringLiteral("use_multi_bonds"), style.useMultiCylindersForBonds);
        style.atomicRadiusScaleFactor = static_cast<float>(
            styledDouble(QStringLiteral("atom_radius_scale"),
                         style.atomicRadiusScaleFactor, 0.0, 2.0));
        style.bondRadius = static_cast<float>(
            styledDouble(QStringLiteral("bond_radius"), style.bondRadius, 0.0, 0.5));
        style.atomicRadiusType = radiusTypeFromString(
            body.value(QStringLiteral("atomic_radius_type")).toString(), style.atomicRadiusType);
        style.atomColorMode = colorModeFromString(
            body.value(QStringLiteral("atom_color_mode")).toString(), style.atomColorMode);
        style.bondColorMode = colorModeFromString(
            body.value(QStringLiteral("bond_color_mode")).toString(), style.bondColorMode);
        style.atomColor = colorFromJson(body.value(QStringLiteral("atom_color")), style.atomColor);
        style.bondColor = colorFromJson(body.value(QStringLiteral("bond_color")), style.bondColor);

        scene_->applyMoleculeStyle(style);
        return jsonResponse(QJsonObject{
            {"style", moleculeStyleToJson(style)},
            {"will_restore_on_clear", true},
        });
    });

    // ---- heroshot: tensor ghost trail -----------------------------------
    //
    // POST /heroshot/ghost_trail  {atom?, n?, end_frame?, step?, frames?}
    //   Draw the focused (or {"atom":N}) atom's shielding tensor at its last
    //   `n` DFT-measured frames as a fading stack of glyphs -- newest opaque,
    //   oldest faint -- so the re-orientation reads "from the side". Walks
    //   backward from `end_frame` (default the live frame) by `step` (default
    //   2, the DFT cadence), keeping only frames whose tensor is valid; up to
    //   `n` (default 5, clamped 2..12) ghosts. Or pass `frames:[...]` to draw
    //   an explicit statistics-selected set such as rotamer wells. Each ghost
    //   is the REAL measured tensor at a REAL frame -- probeAtomCsa reads the
    //   DFT store directly, no interpolation and the live frame never moves.
    //   Heroshot layer only: the reader's own single live glyph is untouched.
    //   Echoes the frames + fades actually drawn so a script can verify what it
    //   sees.
    //   POST /heroshot/clear  removes the trail.
    server_->route(QStringLiteral("/heroshot/ghost_trail"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_ || !scene_)
            return errorResponse(QStringLiteral("reader window / scene not wired"),
                                 SC::ServiceUnavailable);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        if (!protein)
            return errorResponse(QStringLiteral("no protein loaded"), SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);

        std::size_t atom = 0;
        if (body.contains("atom")) {
            const int a = body.value("atom").toInt(-1);
            if (a < 0 || static_cast<std::size_t>(a) >= protein->atomCount())
                return errorResponse(QStringLiteral("atom out of range"), SC::BadRequest);
            atom = static_cast<std::size_t>(a);
        } else if (selection_ && selection_->hasFocus()) {
            atom = selection_->focus();
        } else {
            return errorResponse(QStringLiteral("no atom: pass {\"atom\":N} or focus one"),
                                 SC::BadRequest);
        }

        const int n = std::clamp(body.value("n").toInt(5), 2, 12);
        const int step = std::clamp(body.value("step").toInt(2), 1, 16);
        const bool includeCurrent =
            body.contains(QStringLiteral("include_current")) &&
                    body.value(QStringLiteral("include_current")).isBool()
                ? body.value(QStringLiteral("include_current")).toBool()
                : true;
        const bool hideSelectionMarker =
            body.contains(QStringLiteral("hide_selection_marker")) &&
                    body.value(QStringLiteral("hide_selection_marker")).isBool()
                ? body.value(QStringLiteral("hide_selection_marker")).toBool()
                : false;
        auto styledDouble = [&](const QString& key,
                                double fallback,
                                double lo,
                                double hi) {
            const QJsonValue v = body.value(key);
            if (!v.isDouble()) return fallback;
            return std::clamp(v.toDouble(), lo, hi);
        };
        auto styledBool = [&](const QString& key, bool fallback) {
            const QJsonValue v = body.value(key);
            return v.isBool() ? v.toBool() : fallback;
        };
        TensorGlyphActor::Style historyStyle;
        historyStyle.ovaloidScale =
            styledDouble(QStringLiteral("history_surface_scale"), 1.0, 0.0, 3.0);
        historyStyle.arrowLengthScale =
            styledDouble(QStringLiteral("history_arrow_scale"), 1.0, 0.0, 3.0);
        historyStyle.arrowWidthScale =
            styledDouble(QStringLiteral("history_arrow_width_scale"), 1.0, 0.0, 3.0);
        historyStyle.surfaceOpacity =
            styledDouble(QStringLiteral("history_surface_opacity"), 0.18, 0.0, 1.0);
        historyStyle.arrowOpacity =
            styledDouble(QStringLiteral("history_arrow_opacity"), 0.75, 0.0, 1.0);
        historyStyle.showSurface =
            styledBool(QStringLiteral("history_surface_visible"), false);
        historyStyle.showArrows =
            styledBool(QStringLiteral("history_arrows_visible"), true);
        TensorGlyphActor::Style currentStyle;
        currentStyle.ovaloidScale =
            styledDouble(QStringLiteral("current_surface_scale"), 1.0, 0.0, 3.0);
        currentStyle.arrowLengthScale =
            styledDouble(QStringLiteral("current_arrow_scale"), 1.0, 0.0, 3.0);
        currentStyle.arrowWidthScale =
            styledDouble(QStringLiteral("current_arrow_width_scale"), 1.0, 0.0, 3.0);
        currentStyle.surfaceOpacity =
            styledDouble(QStringLiteral("current_surface_opacity"), 0.50, 0.0, 1.0);
        currentStyle.arrowOpacity =
            styledDouble(QStringLiteral("current_arrow_opacity"), 1.0, 0.0, 1.0);
        currentStyle.showSurface =
            styledBool(QStringLiteral("current_surface_visible"), false);
        currentStyle.showArrows =
            styledBool(QStringLiteral("current_arrows_visible"), true);
        auto applyAxes = [](TensorGlyphActor::Style& style, const QString& mode) {
            if (mode == QStringLiteral("none")) {
                style.showAxes = {false, false, false};
            } else if (mode == QStringLiteral("sigma11") || mode == QStringLiteral("axis0")) {
                style.showAxes = {true, false, false};
            } else if (mode == QStringLiteral("sigma22") || mode == QStringLiteral("axis1")) {
                style.showAxes = {false, true, false};
            } else if (mode == QStringLiteral("sigma33") || mode == QStringLiteral("axis2")) {
                style.showAxes = {false, false, true};
            } else {
                style.showAxes = {true, true, true};
            }
        };
        const QString axesMode = body.value(QStringLiteral("axes")).toString(QStringLiteral("all"));
        applyAxes(historyStyle,
                  body.value(QStringLiteral("history_axes")).toString(axesMode));
        applyAxes(currentStyle,
                  body.value(QStringLiteral("current_axes")).toString(axesMode));

        const int liveFrame = playback_ ? playback_->currentFrame() : 0;
        int endFrame = body.contains("end_frame") ? body.value("end_frame").toInt(liveFrame)
                                                   : liveFrame;
        if (endFrame < 0) endFrame = 0;

        if (hideSelectionMarker && scene_->measurementOverlay()) {
            if (!heroshotMeasurementVisibleBefore_.has_value())
                heroshotMeasurementVisibleBefore_ = scene_->measurementOverlay()->isVisible();
            scene_->measurementOverlay()->setVisible(false);
        }

        std::vector<TensorGhostTrail::Sample> samples;
        std::vector<int> framesKept;
        const auto addFrameSample = [&](int f) {
            const model::AtomCsaResult r = readerWindow_->probeAtomCsa(atom, f);
            if (!r.valid) return false;
            TensorGhostTrail::Sample s;
            s.center = r.atomPos;
            s.principalValues = {r.shape.principal_values[0], r.shape.principal_values[1],
                                 r.shape.principal_values[2]};
            s.pasAxes = r.shape.pas_axes;
            s.iso = r.shape.sigma_iso;
            s.style = (f == endFrame) ? currentStyle : historyStyle;
            samples.push_back(s);
            framesKept.push_back(f);
            return true;
        };

        const QJsonArray explicitFrames =
            body.value(QStringLiteral("frames")).isArray()
                ? body.value(QStringLiteral("frames")).toArray()
                : QJsonArray{};
        if (!explicitFrames.isEmpty()) {
            for (const QJsonValue& v : explicitFrames) {
                if (static_cast<int>(samples.size()) >= 12) break;
                if (!v.isDouble()) continue;
                const int f = v.toInt(-1);
                if (f < 0) continue;
                addFrameSample(f);
            }
        } else {
            // Walk backward from endFrame keeping valid (DFT-present) frames
            // only. Bound the scan so a sparse DFT cadence cannot loop the whole
            // way back.
            const int scanLimit = n * step + 64;
            const int firstFrame = includeCurrent ? endFrame : endFrame - step;
            for (int f = firstFrame, scanned = 0;
                 f >= 0 && static_cast<int>(samples.size()) < n && scanned < scanLimit;
                 f -= step, ++scanned) {
                addFrameSample(f);
            }
        }
        if (samples.empty())
            return errorResponse(
                explicitFrames.isEmpty()
                    ? QStringLiteral("no DFT-valid frames at/under end_frame for this atom")
                    : QStringLiteral("none of the requested frames have valid DFT tensors for this atom"),
                SC::ServiceUnavailable);

        // Fade: newest (index 0 = endFrame) opaque -> oldest faint at kFloor.
        if (!explicitFrames.isEmpty()) {
            endFrame = framesKept.front();
            samples.front().style = currentStyle;
        }
        constexpr double kFloor = 0.20;  // faintest ghost opacity (oldest frame)
        const std::size_t m = samples.size();
        QJsonArray ghosts;
        for (std::size_t i = 0; i < m; ++i) {
            const double t = (m <= 1) ? 0.0
                                      : static_cast<double>(i) / static_cast<double>(m - 1);
            const double op = 1.0 - (1.0 - kFloor) * t;
            samples[i].opacity = op;
            ghosts.append(QJsonObject{
                {"frame", framesKept[i]},
                {"opacity", op},
                {"sigma_iso", samples[i].iso},
                {"role", framesKept[i] == endFrame ? QStringLiteral("current")
                                                   : QStringLiteral("history")},
            });
        }

        // Rebuild against the CURRENT scene renderer each call -- cheap, and it
        // sidesteps any stale-renderer hazard if the run was reloaded.
        heroshotTrail_ = std::make_unique<TensorGhostTrail>(
            vtkSmartPointer<vtkRenderer>(scene_->Renderer()));
        heroshotTrail_->show(samples);
        scene_->requestRender(MoleculeScene::RenderSource::Rest);

        return jsonResponse(QJsonObject{
            {"atom", static_cast<qint64>(atom)},
            {"end_frame", endFrame},
            {"reference_frame", endFrame},
            {"include_current", includeCurrent},
            {"mode", explicitFrames.isEmpty() ? QStringLiteral("cadence")
                                               : QStringLiteral("frames")},
            {"selection_marker_hidden", hideSelectionMarker},
            {"step", step},
            {"requested_n", n},
            {"kept", static_cast<qint64>(m)},
            {"ghosts", ghosts},
        });
    });

    // POST /heroshot/angle_collar
    //   Draw a transient, chart-derived dihedral collar in the molecule scene.
    //   Defaults are aimed at the Asp29 chi2 baton-pass: CA-CB-CG-OD1, with
    //   the signed angle read from the loaded H5 dihedral_time_series rather
    //   than guessed from labels. The collar sits around the CG->CB axis so the
    //   signed sweep matches residue.chi2's chart/stat convention.
    server_->route(QStringLiteral("/heroshot/angle_collar"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !loaded_ || !loaded_->protein || !loaded_->conformation || !transformed_)
            return errorResponse(QStringLiteral("scene / loaded run not wired"),
                                 SC::ServiceUnavailable);
        const model::QtProtein& protein = *loaded_->protein;
        const auto* traj = loaded_->conformation->asTrajectory();
        const auto* h5 = traj ? traj->h5() : nullptr;
        const auto* dihedrals = h5 ? h5->dihedrals() : nullptr;
        if (!dihedrals)
            return errorResponse(QStringLiteral("loaded run has no dihedral time series"),
                                 SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);

        const auto residueOpt = residueIndexFromRequest(protein, body, selection_);
        if (!residueOpt)
            return errorResponse(
                QStringLiteral("residue not found: pass {\"residue_number\":N}, {\"residue\":i}, or focus an atom"),
                SC::BadRequest);
        const std::size_t residue = *residueOpt;
        if (residue >= dihedrals->n_residues)
            return errorResponse(QStringLiteral("residue out of dihedral range"), SC::BadRequest);

        const int chi = std::clamp(body.value(QStringLiteral("chi")).toInt(2), 1, 4);
        const int chiIndex = chi - 1;
        const std::size_t chiExistsIndex = residue * 4u + static_cast<std::size_t>(chiIndex);
        if (chiExistsIndex >= dihedrals->chi_exists.size() || dihedrals->chi_exists[chiExistsIndex] == 0)
            return errorResponse(QStringLiteral("requested chi is not defined for this residue"),
                                 SC::BadRequest);

        const int liveFrame = playback_ ? playback_->currentFrame() : 0;
        const int maxFrame = static_cast<int>(std::min(transformed_->frameCount(), dihedrals->n_frames)) - 1;
        if (maxFrame < 0)
            return errorResponse(QStringLiteral("no frames available"), SC::ServiceUnavailable);
        const int frame = std::clamp(
            body.contains(QStringLiteral("frame")) ? body.value(QStringLiteral("frame")).toInt(liveFrame)
                                                   : liveFrame,
            0, maxFrame);
        const int previousFrame = std::clamp(
            body.contains(QStringLiteral("previous_frame"))
                ? body.value(QStringLiteral("previous_frame")).toInt(frame - 2)
                : frame - 2,
            0, maxFrame);

        const QString atomAName = body.value(QStringLiteral("a")).toString(QStringLiteral("CA"));
        const QString atomBName = body.value(QStringLiteral("b")).toString(QStringLiteral("CB"));
        const QString atomCName = body.value(QStringLiteral("c")).toString(QStringLiteral("CG"));
        const QString atomDName = body.value(QStringLiteral("d")).toString(QStringLiteral("OD1"));
        const auto atomA = findResidueAtomByName(protein, residue, atomAName);
        const auto atomB = findResidueAtomByName(protein, residue, atomBName);
        const auto atomC = findResidueAtomByName(protein, residue, atomCName);
        const auto atomD = findResidueAtomByName(protein, residue, atomDName);
        if (!atomA || !atomB || !atomC || !atomD)
            return errorResponse(QStringLiteral("could not resolve requested dihedral atoms in residue"),
                                 SC::BadRequest);

        auto styledDouble = [&](const QString& key, double fallback, double lo, double hi) {
            const QJsonValue v = body.value(key);
            if (!v.isDouble()) return fallback;
            return std::clamp(v.toDouble(), lo, hi);
        };
        const bool hideSelectionMarker =
            body.contains(QStringLiteral("hide_selection_marker")) &&
                    body.value(QStringLiteral("hide_selection_marker")).isBool()
                ? body.value(QStringLiteral("hide_selection_marker")).toBool()
                : false;

        if (hideSelectionMarker && scene_->measurementOverlay()) {
            if (!heroshotMeasurementVisibleBefore_.has_value())
                heroshotMeasurementVisibleBefore_ = scene_->measurementOverlay()->isVisible();
            scene_->measurementOverlay()->setVisible(false);
        }

        const double angle = dihedrals->chiAt(residue, static_cast<std::size_t>(frame), chiIndex);
        const double previousAngle =
            dihedrals->chiAt(residue, static_cast<std::size_t>(previousFrame), chiIndex);
        if (!std::isfinite(angle) || !std::isfinite(previousAngle))
            return errorResponse(QStringLiteral("requested frame has no finite chi angle"),
                                 SC::ServiceUnavailable);

        const model::Vec3 a = transformed_->atomPosition(static_cast<std::size_t>(frame), *atomA);
        const model::Vec3 b = transformed_->atomPosition(static_cast<std::size_t>(frame), *atomB);
        const model::Vec3 c = transformed_->atomPosition(static_cast<std::size_t>(frame), *atomC);

        AngleCollarActor::Style style;
        style.radius = styledDouble(QStringLiteral("radius"), 1.25, 0.1, 8.0);
        style.tubeRadius = styledDouble(QStringLiteral("tube_radius"), 0.035, 0.002, 0.25);
        style.axisPadding = styledDouble(QStringLiteral("axis_padding"), 0.35, 0.0, 3.0);
        style.coneLength = styledDouble(QStringLiteral("cone_length"), 0.0, 0.0, 8.0);
        style.neckRadius = styledDouble(QStringLiteral("neck_radius"), 0.0, 0.0, 4.0);
        style.rimRadius = styledDouble(QStringLiteral("rim_radius"), 0.0, 0.0, 8.0);
        style.coneOpacity = styledDouble(QStringLiteral("cone_opacity"), 0.26, 0.0, 1.0);
        style.coneDirection =
            body.contains(QStringLiteral("cone_flip")) &&
                    body.value(QStringLiteral("cone_flip")).isBool() &&
                    body.value(QStringLiteral("cone_flip")).toBool()
                ? 1.0
                : -1.0;
        style.arcSegments = std::clamp(body.value(QStringLiteral("segments")).toInt(96), 24, 240);

        std::vector<AngleCollarActor::Arc> arcs;
        arcs.push_back(AngleCollarActor::Arc{
            previousAngle,
            styledDouble(QStringLiteral("previous_opacity"), 0.42, 0.0, 1.0),
            styledDouble(QStringLiteral("previous_radius_scale"), 0.94, 0.2, 3.0),
            std::array<double, 3>{{0.95, 0.30, 0.38}},
        });
        arcs.push_back(AngleCollarActor::Arc{
            angle,
            styledDouble(QStringLiteral("current_opacity"), 0.92, 0.0, 1.0),
            styledDouble(QStringLiteral("current_radius_scale"), 1.04, 0.2, 3.0),
            std::array<double, 3>{{1.00, 0.72, 0.18}},
        });

        heroshotAngleCollar_ = std::make_unique<AngleCollarActor>(
            vtkSmartPointer<vtkRenderer>(scene_->Renderer()));
        // Axis start/end are CG->CB, not CB->CG: this matches the signed
        // residue.chi2 values used by nmrfigs/events_table.
        heroshotAngleCollar_->show(c, b, c, a - b, arcs, style);
        scene_->requestRender(MoleculeScene::RenderSource::Rest);

        const model::QtResidue& rr = protein.residue(residue);
        return jsonResponse(QJsonObject{
            {"residue", static_cast<qint64>(residue)},
            {"residue_number", rr.address.residueNumber},
            {"chi", chi},
            {"frame", frame},
            {"previous_frame", previousFrame},
            {"angle_degrees", radiansToDegrees(angle)},
            {"previous_angle_degrees", radiansToDegrees(previousAngle)},
            {"delta_degrees", radiansToDegrees(angle - previousAngle)},
            {"state", rotamerState(angle)},
            {"previous_state", rotamerState(previousAngle)},
            {"axis", QJsonArray{static_cast<qint64>(*atomC), static_cast<qint64>(*atomB)}},
            {"atoms", QJsonArray{static_cast<qint64>(*atomA), static_cast<qint64>(*atomB),
                                  static_cast<qint64>(*atomC), static_cast<qint64>(*atomD)}},
            {"atom_names", QJsonArray{atomAName, atomBName, atomCName, atomDName}},
            {"center", vec3FromEigen(c)},
            {"collar_shape", QStringLiteral("cone")},
            {"cone_opens_opposite_axis", style.coneDirection < 0.0},
            {"selection_marker_hidden", hideSelectionMarker},
        });
    });

    server_->route(QStringLiteral("/heroshot/clear"), Method::Post, [this]() {
        ASSERT_THREAD(this);
        if (heroshotTrail_) heroshotTrail_->clear();
        if (heroshotAngleCollar_) heroshotAngleCollar_->clear();
        if (scene_ && scene_->measurementOverlay()
            && heroshotMeasurementVisibleBefore_.has_value()) {
            scene_->measurementOverlay()->setVisible(*heroshotMeasurementVisibleBefore_);
            heroshotMeasurementVisibleBefore_.reset();
        }
        if (scene_ && heroshotMoleculeStyleBefore_.has_value()) {
            scene_->applyMoleculeStyle(*heroshotMoleculeStyleBefore_);
            heroshotMoleculeStyleBefore_.reset();
        }
        if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- catalog dump (the full "seeable list", for auditing) -----------
    //
    // GET /catalog -> every signal descriptor with its real state: provenance
    // (source), axis, value shape, storage path, temporal/static, the display
    // modes it offers, and -- crucially -- whether it is genuinely SHOWABLE for
    // the loaded run (field availability), not merely declared. Plus a summary
    // by source / availability / value-shape so the catalog's shape is visible
    // at a glance. The audit surface for the displayable-catalog cleanup.
    server_->route(QStringLiteral("/catalog"), [this]() {
        ASSERT_THREAD(this);
        if (!catalog_)
            return errorResponse(QStringLiteral("catalog not wired"), SC::ServiceUnavailable);
        const auto& descriptors = catalog_->allDescriptorList();
        QJsonArray arr;
        QJsonObject bySource, byAvail, byShape;
        int showable = 0;
        for (const model::SignalDescriptor& d : descriptors) {
            const QString source = model::ToString(d.sourceKind);
            const QString avail = availabilityString(catalog_, d.id);
            const QString shape = model::ToString(d.valueShape);
            bySource[source] = bySource.value(source).toInt() + 1;
            byAvail[avail] = byAvail.value(avail).toInt() + 1;
            byShape[shape] = byShape.value(shape).toInt() + 1;
            const bool isShowable = (avail == QStringLiteral("Available")
                                     || avail == QStringLiteral("AllZeroObserved"));
            if (isShowable) ++showable;
            QJsonArray modeArr;
            for (const QString& m : model::AllDisplayModes(d)) modeArr.append(m);
            QJsonArray tagArr;
            for (const QString& t : d.tags) tagArr.append(t);
            arr.append(QJsonObject{
                {"id", d.id},
                {"label", d.label},
                {"family", d.family},
                {"concept", d.conceptKey},
                {"import_set", d.importSet},
                {"source", source},
                {"axis", model::ToString(d.nativeAxis)},
                {"value_shape", shape},
                {"storage", d.storagePath},
                {"temporal", d.temporal},
                {"static", d.staticDisplay},
                {"availability", avail},
                {"showable", isShowable},
                {"sampling", model::ToString(d.samplingStatus)},
                {"modes", modeArr},
                {"tags", tagArr},
            });
        }
        return jsonResponse(QJsonObject{
            {"count", static_cast<qint64>(descriptors.size())},
            {"summary", QJsonObject{
                {"showable", showable},
                {"by_source", bySource},
                {"by_availability", byAvail},
                {"by_value_shape", byShape},
            }},
            {"descriptors", arr},
        });
    });

    // GET /catalog/tree -- the same catalog, organized mechanism -> concept ->
    // form (model/MetricTaxonomy.h). Sane hypothesis-first grouping: the four
    // shielding-contribution kernels (role=hypothesis), the DFT/ProCS15 reference
    // (role=reference), the conditioning inputs, dynamics, scaffold. The ~188 flat
    // rows fold to ~35 base concepts; each concept lists its FORMS (snapshot /
    // series / rollup / ...) as related items, not a hidden collapse.
    server_->route(QStringLiteral("/catalog/tree"), [this]() {
        ASSERT_THREAD(this);
        if (!catalog_)
            return errorResponse(QStringLiteral("catalog not wired"), SC::ServiceUnavailable);
        const QVector<model::MetricGroupNode> tree =
            model::GroupCatalog(catalog_->allDescriptorList());
        QJsonArray groups;
        int conceptTotal = 0, formTotal = 0;
        QJsonObject groupsByRole;
        for (const model::MetricGroupNode& g : tree) {
            QJsonArray concepts;
            for (const model::MetricConceptNode& c : g.concepts) {
                ++conceptTotal;
                QJsonArray forms;
                bool anyAvailable = false;
                for (const model::MetricFormEntry& f : c.forms) {
                    ++formTotal;
                    const QString avail = availabilityString(catalog_, f.descriptorId);
                    if (avail == QStringLiteral("Available")) anyAvailable = true;
                    forms.append(QJsonObject{
                        {"form", QString::fromLatin1(model::ToString(f.form))},
                        {"id", f.descriptorId},
                        {"availability", avail},
                    });
                }
                QJsonObject cobj{
                    {"concept", c.baseConcept},
                    {"label", c.label},
                    {"available", anyAvailable},
                    {"forms", forms},
                };
                if (!c.chargeModel.isEmpty()) cobj["charge_model"] = c.chargeModel;
                concepts.append(cobj);
            }
            const QString role = QString::fromLatin1(model::ToString(g.role));
            groupsByRole[role] = groupsByRole.value(role).toInt() + 1;
            groups.append(QJsonObject{
                {"group", QString::fromLatin1(model::ToString(g.group))},
                {"role", role},
                {"concept_count", static_cast<int>(concepts.size())},
                {"concepts", concepts},
            });
        }
        return jsonResponse(QJsonObject{
            {"descriptor_count", static_cast<qint64>(catalog_->allDescriptorList().size())},
            {"group_count", static_cast<int>(groups.size())},
            {"concept_count", conceptTotal},
            {"form_count", formTotal},
            {"groups_by_role", groupsByRole},
            {"groups", groups},
        });
    });

    // GET /catalog/display-audit -- coherence of the per-field DISPLAY OPTIONS.
    // For every descriptor: its value_shape, the modes it OFFERS (AllDisplayModes),
    // and the visualizations that actually SUPPORT it (registry.supporting). A
    // field whose supporting set is EMPTY can never be shown -- a dead option.
    // Aggregated by value_shape so incoherence is visible at a glance (e.g. a
    // shape that no viz supports, or a shape offering modes nothing renders).
    server_->route(QStringLiteral("/catalog/display-audit"), [this]() {
        ASSERT_THREAD(this);
        if (!catalog_)
            return errorResponse(QStringLiteral("catalog not wired"), SC::ServiceUnavailable);
        const model::VisualizationRegistry& reg = model::VisualizationRegistry::instance();
        QJsonArray fields;
        QJsonArray deadFields;            // intended as a signal but NO viz renders it (a real bug)
        QJsonArray nonDisplayableFields;  // intentionally not a dashboard signal (DisplayPolicy)
        QHash<QString, int> shapeCount, shapeDead;
        QHash<QString, QStringList> shapeViz;
        for (const model::SignalDescriptor& d : catalog_->allDescriptorList()) {
            const QString shape = model::ToString(d.valueShape);
            QJsonArray modeArr;
            for (const QString& m : model::AllDisplayModes(d)) modeArr.append(m);
            QJsonArray vizArr;
            for (const model::VisualizationDefinition* def : reg.supporting(d)) {
                const QString vt = model::ToString(def->type());
                vizArr.append(vt);
                if (!shapeViz[shape].contains(vt)) shapeViz[shape].append(vt);
            }
            const bool intended = model::IsDashboardDisplayable(d);
            const bool hasViz = !vizArr.isEmpty();
            const bool dead = intended && !hasViz;  // claims to be a signal but nothing renders it
            ++shapeCount[shape];
            if (dead) { ++shapeDead[shape]; deadFields.append(d.id); }
            if (!intended) nonDisplayableFields.append(d.id);
            fields.append(QJsonObject{
                {"id", d.id},
                {"group", QString::fromLatin1(model::ToString(model::ClassifyMetric(d).group))},
                {"value_shape", shape},
                {"displayable", intended},
                {"offered_modes", modeArr},
                {"supporting_viz", vizArr},
                {"dead", dead},
            });
        }
        QJsonObject byShape;
        for (auto it = shapeCount.constBegin(); it != shapeCount.constEnd(); ++it) {
            QJsonArray viz;
            for (const QString& v : shapeViz.value(it.key())) viz.append(v);
            byShape[it.key()] = QJsonObject{
                {"fields", it.value()},
                {"dead", shapeDead.value(it.key())},
                {"supporting_viz", viz},
            };
        }
        return jsonResponse(QJsonObject{
            {"total", static_cast<qint64>(catalog_->allDescriptorList().size())},
            {"dead_count", deadFields.size()},
            {"dead_fields", deadFields},
            {"non_displayable_count", nonDisplayableFields.size()},
            {"non_displayable_fields", nonDisplayableFields},
            {"by_value_shape", byShape},
            {"fields", fields},
        });
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

    server_->route(QStringLiteral("/dashboard/strip/series"), [this]() {
        ASSERT_THREAD(this);
        if (!dashboardController_)
            return errorResponse(QStringLiteral("dashboard display controller not wired"),
                                 SC::ServiceUnavailable);

        QJsonArray tracks;
        const QVector<DashboardDisplayController::StripTrack> stripTracks =
            dashboardController_->stripTracks();
        for (const DashboardDisplayController::StripTrack& track : stripTracks)
            tracks.append(stripTrackToJson(track));
        return jsonResponse(QJsonObject{{"tracks", tracks}});
    });

    // GET /dashboard/display -- the unified display MANIFEST: every active
    // display element + the data it actually plots. Temporal strip tracks (the
    // per-frame sampled series) PLUS the static panels (curve/spectrum/matrix/
    // fixed-freq/sequence-bar), which the strip-series path could not see. The
    // read-to-display test reads this to verify each metric reaches a sane,
    // renderable form -- closing the static-panel blind spot.
    server_->route(QStringLiteral("/dashboard/display"), [this]() {
        ASSERT_THREAD(this);
        if (!dashboardController_ || !readerWindow_)
            return errorResponse(QStringLiteral("dashboard not wired"), SC::ServiceUnavailable);
        QJsonArray tracks;
        for (const DashboardDisplayController::StripTrack& t : dashboardController_->stripTracks())
            tracks.append(stripTrackToJson(t));
        return jsonResponse(QJsonObject{
            {"strip_tracks", tracks},
            {"panels", readerWindow_->dashboardPanelManifest()},
        });
    });

    server_->route(QStringLiteral("/dashboard/state"), [this]() {
        ASSERT_THREAD(this);
        if (!signalModel_ || !panelModel_ || !selectionController_ || !catalog_)
            return errorResponse(QStringLiteral("dashboard models not wired"), SC::ServiceUnavailable);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"), SC::ServiceUnavailable);

        const QJsonArray selected = selectedStateArray(signalModel_.data(), catalog_);
        const DashboardSmokeSummary smoke =
            dashboardController_ ? dashboardController_->smokeSummary()
                                 : DashboardSmokeSummary{};
        return jsonResponse(QJsonObject{
            {"selected", selected},
            {"selected_count", signalModel_->selectedCount()},
            {"dock", QJsonObject{
                {"visible", readerWindow_->dashboardDockVisible()},
                {"width", readerWindow_->dashboardDockWidth()},
                {"raised", readerWindow_->dashboardDockRaised()},
            }},
            {"render", QJsonObject{
                {"owned_panel_count", readerWindow_->dashboardOwnedPanelCount()},
                {"strip_track_count", readerWindow_->dashboardStripTrackCount()},
                {"expected_empty", expectedEmptyArray(smoke)},
            }},
        });
    });

    server_->route(QStringLiteral("/dashboard/picker"), [this]() {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"), SC::ServiceUnavailable);
        return jsonResponse(readerWindow_->signalDisplayPickerState());
    });

    server_->route(QStringLiteral("/dashboard/picker/open"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"), SC::ServiceUnavailable);
        if (!selection_)
            return errorResponse(QStringLiteral("selection not wired"), SC::ServiceUnavailable);

        QJsonObject body;
        if (!req.body().isEmpty()) {
            bool ok = false;
            body = parseJsonBody(req, &ok);
            if (!ok)
                return errorResponse(QStringLiteral("body must be a JSON object"), SC::BadRequest);
        }

        if (body.contains(QStringLiteral("atom"))) {
            QString error;
            qint64 atom = -1;
            if (!readNonNegativeInteger(body, QStringLiteral("atom"), &atom, &error))
                return errorResponse(error, SC::BadRequest);
            const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
            if (!protein)
                return errorResponse(QStringLiteral("no protein loaded"), SC::ServiceUnavailable);
            if (!indexInRange(atom, protein->atomCount(), QStringLiteral("atom"), &error))
                return errorResponse(error, SC::BadRequest);
            selection_->applyPick(static_cast<std::size_t>(atom), Qt::NoModifier);
        }

        QString blockedReason;
        if (!readerWindow_->openSignalDisplayPicker(&blockedReason)) {
            return errorResponse(blockedReason.isEmpty()
                                     ? QStringLiteral("Metric Picker open was blocked")
                                     : blockedReason,
                                 SC::Conflict);
        }
        return jsonResponse(readerWindow_->signalDisplayPickerState());
    });

    server_->route(QStringLiteral("/dashboard/picker/add"), Method::Post,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"), SC::ServiceUnavailable);
        if (!signalModel_)
            return errorResponse(QStringLiteral("signal model not wired"), SC::ServiceUnavailable);

        const QJsonObject picker = readerWindow_->addSelectedSignalFromPicker();
        return jsonResponse(QJsonObject{
            {"picker", picker},
            {"selected_count", signalModel_->selectedCount()},
            {"dock", QJsonObject{
                {"visible", readerWindow_->dashboardDockVisible()},
                {"width", readerWindow_->dashboardDockWidth()},
                {"raised", readerWindow_->dashboardDockRaised()},
            }},
        });
    });

    server_->route(QStringLiteral("/dashboard/metric"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!signalModel_ || !panelModel_ || !catalog_)
            return errorResponse(QStringLiteral("dashboard models not wired"), SC::ServiceUnavailable);
        if (!loaded_)
            return errorResponse(QStringLiteral("loaded run not wired"), SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("body must be a JSON object"), SC::BadRequest);

        const QString descriptorId = body.value(QStringLiteral("descriptor_id")).toString().trimmed();
        if (descriptorId.isEmpty())
            return errorResponse(QStringLiteral("body must include non-empty string field \"descriptor_id\""),
                                 SC::BadRequest);
        const model::SignalDescriptor* descriptor = catalog_->findDescriptor(descriptorId);
        if (!descriptor) {
            if (findAnyDescriptor(catalog_, descriptorId)) {
                return errorResponse(QStringLiteral("descriptor not available: %1")
                                         .arg(availabilityString(catalog_, descriptorId)),
                                     SC::Conflict);
            }
            return errorResponse(QStringLiteral("descriptor not found: %1").arg(descriptorId),
                                 SC::NotFound);
        }

        QString error;
        ParsedDashboardAnchor parsedAnchor;
        if (!parseDashboardAnchor(body, loaded_, &parsedAnchor, &error))
            return errorResponse(error, SC::BadRequest);

        QStringList modes;
        if (!parseModeArray(body, &modes, &error))
            return errorResponse(error, SC::BadRequest);

        for (const QString& mode : modes) {
            model::DisplaySignalBinding binding;
            binding.sourceKind = descriptor->sourceKind;
            binding.descriptorId = descriptor->id;
            binding.conceptKey = descriptor->conceptKey;
            binding.displayModeId = mode;
            binding.anchor = parsedAnchor.anchor;
            binding.followsFocus = parsedAnchor.followsFocus;
            if (!model::SupportsDisplayMode(*descriptor, mode)) {
                return errorResponse(QStringLiteral("descriptor %1 does not support mode %2")
                                         .arg(descriptor->id, mode),
                                     SC::BadRequest);
            }
            if (!catalog_->canBind(binding)) {
                return errorResponse(QStringLiteral("descriptor %1 mode %2 cannot bind to anchor")
                                         .arg(descriptor->id, mode),
                                     SC::UnprocessableEntity);
            }
        }

        const QString label = QStringLiteral("%1 - %2").arg(descriptor->label, parsedAnchor.label);
        int addedRefs = 0;
        DashboardSelectionController::PanelTarget target;
        target.makeActive = true;
        const QUuid id = selectionController_->addMetric(*descriptor,
                                                         parsedAnchor.anchor,
                                                         modes,
                                                         target,
                                                         parsedAnchor.followsFocus,
                                                         label,
                                                         &addedRefs);

        return jsonResponse(QJsonObject{
            {"id", uuidToString(id)},
            {"added_refs", addedRefs},
        });
    });

    server_->route(QStringLiteral("/dashboard/metric/remove"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!signalModel_ || !selectionController_)
            return errorResponse(QStringLiteral("signal model not wired"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("body must be a JSON object"), SC::BadRequest);
        QString error;
        const std::optional<QUuid> id = uuidFromBody(body, &error);
        if (!id)
            return errorResponse(error, SC::BadRequest);
        if (!signalModel_->signalById(*id))
            return errorResponse(QStringLiteral("signal not found: %1").arg(uuidToString(*id)),
                                 SC::NotFound);
        if (!selectionController_->removeMetric(*id))
            return errorResponse(QStringLiteral("signal remove failed: %1").arg(uuidToString(*id)),
                                 SC::InternalServerError);
        return jsonResponse(QJsonObject{{"removed", true}});
    });

    server_->route(QStringLiteral("/dashboard/metric/mode"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!signalModel_ || !panelModel_ || !selectionController_ || !catalog_)
            return errorResponse(QStringLiteral("dashboard models not wired"), SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("body must be a JSON object"), SC::BadRequest);
        QString error;
        const std::optional<QUuid> id = uuidFromBody(body, &error);
        if (!id)
            return errorResponse(error, SC::BadRequest);
        const QString mode = body.value(QStringLiteral("mode")).toString().trimmed();
        if (mode.isEmpty())
            return errorResponse(QStringLiteral("body must include non-empty string field \"mode\""),
                                 SC::BadRequest);
        if (!body.contains(QStringLiteral("enabled")) || !body.value(QStringLiteral("enabled")).isBool())
            return errorResponse(QStringLiteral("body must include bool field \"enabled\""),
                                 SC::BadRequest);
        const bool enabled = body.value(QStringLiteral("enabled")).toBool();

        const model::DashboardSignal* before = signalModel_->signalById(*id);
        if (!before)
            return errorResponse(QStringLiteral("signal not found: %1").arg(uuidToString(*id)),
                                 SC::NotFound);
        const model::SignalDescriptor* descriptor = catalog_->findDescriptor(before->binding.descriptorId);
        if (!descriptor)
            return errorResponse(QStringLiteral("descriptor not available: %1")
                                     .arg(before->binding.descriptorId),
                                 SC::Conflict);
        if (!model::SupportsDisplayMode(*descriptor, mode)) {
            return errorResponse(QStringLiteral("descriptor %1 does not support mode %2")
                                     .arg(descriptor->id, mode),
                                 SC::BadRequest);
        }
        model::DisplaySignalBinding binding = before->binding;
        binding.displayModeId = mode;
        if (!catalog_->canBind(binding)) {
            return errorResponse(QStringLiteral("descriptor %1 mode %2 cannot bind to signal anchor")
                                     .arg(descriptor->id, mode),
                                 SC::UnprocessableEntity);
        }

        if (!selectionController_->setMetricMode(*id, mode, enabled)) {
            return errorResponse(QStringLiteral("display mode toggle failed"), SC::Conflict);
        }

        const model::DashboardSignal* after = signalModel_->signalById(*id);
        return jsonResponse(QJsonObject{
            {"modes", after ? stringListToJson(after->displayModeIds) : QJsonArray{}},
        });
    });

    server_->route(QStringLiteral("/dashboard/dock"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!readerWindow_)
            return errorResponse(QStringLiteral("reader main window not wired"), SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains(QStringLiteral("visible")) || !body.value(QStringLiteral("visible")).isBool())
            return errorResponse(QStringLiteral("body must be {\"visible\": bool}"), SC::BadRequest);
        readerWindow_->setDashboardDockVisible(body.value(QStringLiteral("visible")).toBool());
        return jsonResponse(QJsonObject{
            {"visible", readerWindow_->dashboardDockVisible()},
            {"width", readerWindow_->dashboardDockWidth()},
        });
    });

    // ---- shutdown -------------------------------------------------------

    // POST /shutdown — graceful exit so the operator (or a test harness)
    // doesn't need pkill. Returns 204, then exits on the socket's flush event
    // (bytesWritten -> bytesToWrite()==0), NOT a timer: the response write is
    // itself async, so any clock races it and truncates the reply. Quitting
    // synchronously in the handler would tear down before the bytes ship.
    server_->route(QStringLiteral("/shutdown"), Method::Post,
                   [this](const QHttpServerRequest&) {
        ASSERT_THREAD(this);
        QTcpSocket* sock = activeSocket_.data();
        if (sock) {
            qCInfo(cRest).noquote() << "REST /shutdown — exiting once the 204 has flushed";
            QObject::connect(sock, &QTcpSocket::bytesWritten, qApp, [sock]() {
                if (sock->bytesToWrite() == 0)
                    QCoreApplication::quit();
            });
            QObject::connect(sock, &QTcpSocket::disconnected,
                             qApp, &QCoreApplication::quit);
        } else {
            // Qt 6.4 listen() path: no socket handle captured — exit directly.
            qCInfo(cRest).noquote() << "REST /shutdown — quitting (no socket handle)";
            QCoreApplication::quit();
        }
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
        // Poster / print export: scale > 1 supersamples the scene capture for a
        // high-DPI, smoother image (target="scene" only). Clamp to 1..8.
        int scale = (ok && body.contains(QStringLiteral("scale")))
                        ? body.value(QStringLiteral("scale")).toInt(1)
                        : 1;
        if (scale < 1) scale = 1;
        else if (scale > 8) scale = 8;
        QByteArray png;
        if (target == QStringLiteral("window")) {
            png = captureWindowPng(mainWindow_.data());
        } else if (target == QStringLiteral("scene")) {
            png = captureScenePng(scene_.data(), forceRender, scale);
        } else if (target == QStringLiteral("widget")) {
            const QString objectName = body.value(QStringLiteral("object_name")).toString();
            if (objectName.isEmpty())
                return errorResponse(QStringLiteral("target \"widget\" requires \"object_name\""), SC::BadRequest);
            QWidget* widget = findWidgetByObjectName(mainWindow_.data(), objectName);
            if (!widget)
                return errorResponse(QStringLiteral("no live widget named \"%1\"").arg(objectName), SC::NotFound);
            png = captureWindowPng(widget);
        } else {
            return errorResponse(QStringLiteral("target must be \"scene\", \"window\", or \"widget\""), SC::BadRequest);
        }
        if (png.isEmpty())
            return errorResponse(QStringLiteral("screenshot capture failed"), SC::InternalServerError);
        return QHttpServerResponse(QByteArray(kMimePng), png);
    });
}

}  // namespace h5reader::app
