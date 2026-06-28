#include "RestServer.h"

#include "AngleCollarActor.h"
#include "AtomTrackOverlay.h"
#include "CameraAnchorHelper.h"
#include "CameraComposer.h"
#include "DashboardDisplayController.h"
#include "DashboardSelectionController.h"
#include "HeroshotButterflyOverlay.h"
#include "MeasurementOverlay.h"
#include "MoleculeScene.h"
#include "NewmanProjection.h"
#include "QtFieldGridOverlay.h"
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
#include "../model/ConformationGeometry.h"
#include "../model/Conformation.h"
#include "../model/DashboardPanelModel.h"
#include "../model/DashboardSignal.h"
#include "../model/DashboardSignalModel.h"
#include "../model/DisplayPolicy.h"
#include "../model/MetricTaxonomy.h"
#include "../model/QtProtein.h"
#include "../model/RingCurrentFaceCollar.h"
#include "../model/RingNullCollar.h"
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
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <utility>
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

QJsonValue finiteJson(double value) {
    return std::isfinite(value) ? QJsonValue(value) : QJsonValue(QJsonValue::Null);
}

QJsonObject restNamespace(const QString& prefix,
                          const QString& tier,
                          const QString& audience,
                          const QString& contract,
                          const QString& notes) {
    return QJsonObject{
        {"prefix", prefix},
        {"tier", tier},
        {"audience", audience},
        {"contract", contract},
        {"notes", notes},
    };
}

QJsonObject restRoute(const QString& method,
                      const QString& path,
                      const QString& tier,
                      const QString& summary) {
    return QJsonObject{
        {"method", method},
        {"path", path},
        {"tier", tier},
        {"summary", summary},
    };
}

QJsonObject restInterfaceDescription() {
    return QJsonObject{
        {"version", 1},
        {"namespaces", QJsonArray{
            restNamespace(QStringLiteral("/api"),
                          QStringLiteral("general"),
                          QStringLiteral("user/mcp"),
                          QStringLiteral("stable intent; suitable for automation clients"),
                          QStringLiteral("Promoted routes whose semantics are useful outside testing.")),
            restNamespace(QStringLiteral("/field"),
                          QStringLiteral("general"),
                          QStringLiteral("user/mcp"),
                          QStringLiteral("stable scene field display controls"),
                          QStringLiteral("Normal reader visualization controls, including ring-null context geometry.")),
            restNamespace(QStringLiteral("/resthero"),
                          QStringLiteral("figure_composition"),
                          QStringLiteral("operator/figure harness"),
                          QStringLiteral("transient scene styling; clear with /resthero/clear"),
                          QStringLiteral("Poster and review tools may move the scene away from normal UI state.")),
            restNamespace(QStringLiteral("/catalog"),
                          QStringLiteral("diagnostic"),
                          QStringLiteral("tests/development"),
                          QStringLiteral("display catalog audit"),
                          QStringLiteral("Inventory and coherence checks for the loaded signal catalog.")),
            restNamespace(QStringLiteral("/diagnostics"),
                          QStringLiteral("diagnostic"),
                          QStringLiteral("tests/development"),
                          QStringLiteral("inspection and snapshot utilities"),
                          QStringLiteral("Allowed to expose implementation detail; not the product API.")),
        }},
        {"routes", QJsonArray{
            restRoute(QStringLiteral("GET"), QStringLiteral("/api/interface"),
                      QStringLiteral("general"),
                      QStringLiteral("Machine-readable REST namespace map.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/api/screenshot"),
                      QStringLiteral("general"),
                      QStringLiteral("Scene/window screenshot capture for review, export, and automation.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/api/ring/null_crossings"),
                      QStringLiteral("general"),
                      QStringLiteral("Operational ring-null crossing collection.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/api/ring/current_face_collar"),
                      QStringLiteral("general"),
                      QStringLiteral("Ring-current weak-signal receiver and fit summary.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/field/null_cone"),
                      QStringLiteral("general"),
                      QStringLiteral("Regular-use ring null surface display.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/resthero/butterfly"),
                      QStringLiteral("figure_composition"),
                      QStringLiteral("High-resolution transient ring-current surfaces.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/resthero/atom_track"),
                      QStringLiteral("figure_composition"),
                      QStringLiteral("Transient point-cloud track for sampled atom positions.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/resthero/ghost_trail"),
                      QStringLiteral("figure_composition"),
                      QStringLiteral("Transient tensor ghost trail.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/resthero/clear"),
                      QStringLiteral("figure_composition"),
                      QStringLiteral("Clear transient resthero scene state.")),
            restRoute(QStringLiteral("POST"), QStringLiteral("/diagnostics/screenshot"),
                      QStringLiteral("diagnostic"),
                      QStringLiteral("Qt widget/framebuffer screenshot capture for tests.")),
            restRoute(QStringLiteral("GET"), QStringLiteral("/catalog"),
                      QStringLiteral("diagnostic"),
                      QStringLiteral("Full display catalog audit surface.")),
        }},
    };
}

QJsonArray mat3ToJson(const model::Mat3& m) {
    QJsonArray rows;
    for (int r = 0; r < 3; ++r)
        rows.append(QJsonArray{m(r, 0), m(r, 1), m(r, 2)});
    return rows;
}

QJsonArray array3ToJson(const std::array<double, 3>& values) {
    return QJsonArray{values[0], values[1], values[2]};
}

QJsonArray array5ToJson(const std::array<double, 5>& values) {
    return QJsonArray{values[0], values[1], values[2], values[3], values[4]};
}

QJsonArray intVectorToJson(const std::vector<int>& values) {
    QJsonArray out;
    for (int value : values)
        out.append(value);
    return out;
}

QJsonArray doubleVectorToJson(const std::vector<double>& values) {
    QJsonArray out;
    for (double value : values)
        out.append(finiteJson(value));
    return out;
}

struct RingLocalFrame {
    bool valid = false;
    model::RingGeometry geometry;
    model::Vec3 u = model::Vec3::Zero();
    model::Vec3 v = model::Vec3::Zero();
    model::Vec3 n = model::Vec3::Zero();
};

RingLocalFrame ringLocalFrameAt(const model::Conformation& conf,
                                std::size_t ringIdx,
                                std::size_t frame) {
    RingLocalFrame out;
    const std::vector<model::Vec3> verts = model::RingVertices(conf, ringIdx, frame);
    out.geometry = model::FitRingGeometry(verts);

    const double nNorm = out.geometry.normal.norm();
    if (verts.empty() || out.geometry.radius < 1e-9 || nNorm < 1e-12)
        return out;

    out.n = out.geometry.normal / nNorm;
    for (const model::Vec3& vert : verts) {
        model::Vec3 radial = vert - out.geometry.center;
        radial -= out.n * radial.dot(out.n);
        const double rNorm = radial.norm();
        if (rNorm > 1e-9) {
            out.u = radial / rNorm;
            break;
        }
    }
    if (out.u.norm() < 1e-9) {
        const model::RingOrthoBasis fallback = model::OrthoBasisFromNormal(out.n);
        out.u = fallback.u;
    }
    out.v = out.n.cross(out.u);
    const double vNorm = out.v.norm();
    if (vNorm < 1e-12)
        return out;
    out.v /= vNorm;
    out.u = out.v.cross(out.n).normalized();
    out.valid = true;
    return out;
}

model::Vec3 toRingLocal(const RingLocalFrame& frame, const model::Vec3& world) {
    const model::Vec3 delta = world - frame.geometry.center;
    return model::Vec3(delta.dot(frame.u), delta.dot(frame.v), delta.dot(frame.n));
}

model::Vec3 fromRingLocal(const RingLocalFrame& frame, const model::Vec3& local) {
    return frame.geometry.center
        + frame.u * local.x()
        + frame.v * local.y()
        + frame.n * local.z();
}

double ringCurrentExpectedValue(const model::RingNullMeasurement& m) {
    if (!m.valid || m.distanceA <= 1e-12)
        return std::numeric_limits<double>::quiet_NaN();
    return m.angularFactor / (m.distanceA * m.distanceA * m.distanceA);
}

QJsonObject sphericalTensorToJson(const model::SphericalTensor& tensor) {
    return QJsonObject{
        {"T0", finiteJson(tensor.T0)},
        {"T1", array3ToJson(tensor.T1)},
        {"T2", array5ToJson(tensor.T2)},
        {"T2_magnitude", finiteJson(tensor.T2Magnitude())},
    };
}

QJsonObject efgToJson(const model::QtEfg& efg) {
    return QJsonObject{
        {"T2", array5ToJson(efg.t2)},
        {"T2_magnitude", finiteJson(efg.t2Magnitude())},
    };
}

QJsonObject mopacScalarsToJson(const model::MopacScalars& scalars) {
    return QJsonObject{
        {"charge", finiteJson(scalars.charge)},
        {"sPop", finiteJson(scalars.sPop)},
        {"pPop", finiteJson(scalars.pPop)},
        {"valency", finiteJson(scalars.valency)},
    };
}

QJsonObject coulombScalarsToJson(const model::CoulombScalars& scalars) {
    return QJsonObject{
        {"E_magnitude", finiteJson(scalars.E_magnitude)},
        {"E_bond_proj", finiteJson(scalars.E_bond_proj)},
        {"E_backbone_signed_projection", finiteJson(scalars.E_backbone_frac)},
        {"aromatic_E_magnitude", finiteJson(scalars.aromatic_E_magnitude)},
    };
}

QJsonObject mcScalarsToJson(const model::McConnellScalars& scalars) {
    return QJsonObject{
        {"co_sum", finiteJson(scalars.co_sum)},
        {"cn_sum", finiteJson(scalars.cn_sum)},
        {"sidechain_sum", finiteJson(scalars.sidechain_sum)},
        {"aromatic_sum", finiteJson(scalars.aromatic_sum)},
        {"nearest_CO_dist", finiteJson(scalars.nearest_CO_dist)},
        {"nearest_CN_dist", finiteJson(scalars.nearest_CN_dist)},
        {"has_nearest_CO", scalars.hasNearestCO()},
        {"has_nearest_CN", scalars.hasNearestCN()},
    };
}

QJsonArray mcCategoryT2ToJson(const model::PerBondCategoryT2& categories) {
    static constexpr std::array<const char*, model::kMcConnellCategoryCount> kNames{
        "backbone_total",
        "sidechain_total",
        "aromatic_total",
        "co_nearest",
        "cn_nearest",
    };
    QJsonArray out;
    for (int i = 0; i < model::kMcConnellCategoryCount; ++i) {
        const auto& values = categories.byCategory[static_cast<std::size_t>(i)];
        double magnitude = 0.0;
        for (double value : values)
            magnitude += value * value;
        out.append(QJsonObject{
            {"category", QString::fromLatin1(kNames[static_cast<std::size_t>(i)])},
            {"T2", array5ToJson(values)},
            {"T2_magnitude", finiteJson(std::sqrt(magnitude))},
        });
    }
    return out;
}

QJsonObject mopacSignalsToJson(const model::RingNullMopacSignals& observed) {
    QJsonObject charge{{"present", observed.chargePresent}};
    if (observed.chargePresent)
        charge.insert(QStringLiteral("value"), finiteJson(observed.charge));

    QJsonObject core{{"present", observed.coreScalarsPresent}};
    if (observed.coreScalarsPresent)
        core.insert(QStringLiteral("scalars"), mopacScalarsToJson(observed.coreScalars));

    QJsonObject coulombE{{"present", observed.coulombEPresent}};
    if (observed.coulombEPresent) {
        coulombE.insert(QStringLiteral("vector"), vec3FromEigen(observed.coulombE));
        coulombE.insert(QStringLiteral("magnitude"), finiteJson(observed.coulombE.norm()));
    }

    QJsonObject coulombScalars{{"present", observed.coulombScalarsPresent}};
    if (observed.coulombScalarsPresent)
        coulombScalars.insert(QStringLiteral("scalars"),
                              coulombScalarsToJson(observed.coulombScalars));

    QJsonObject coulombShielding{{"present", observed.coulombShieldingPresent}};
    if (observed.coulombShieldingPresent)
        coulombShielding.insert(QStringLiteral("tensor"),
                                sphericalTensorToJson(observed.coulombShielding));

    QJsonObject efgBackbone{{"present", observed.coulombEfgBackbonePresent}};
    if (observed.coulombEfgBackbonePresent)
        efgBackbone.insert(QStringLiteral("efg"), efgToJson(observed.coulombEfgBackbone));

    QJsonObject efgAromatic{{"present", observed.coulombEfgAromaticPresent}};
    if (observed.coulombEfgAromaticPresent)
        efgAromatic.insert(QStringLiteral("efg"), efgToJson(observed.coulombEfgAromatic));

    QJsonObject mcShielding{{"present", observed.mcShieldingPresent}};
    if (observed.mcShieldingPresent)
        mcShielding.insert(QStringLiteral("tensor"),
                           sphericalTensorToJson(observed.mcShielding));

    QJsonObject mcCategory{{"present", observed.mcCategoryT2Present}};
    if (observed.mcCategoryT2Present)
        mcCategory.insert(QStringLiteral("categories"),
                          mcCategoryT2ToJson(observed.mcCategoryT2));

    QJsonObject mcScalars{{"present", observed.mcScalarsPresent}};
    if (observed.mcScalarsPresent)
        mcScalars.insert(QStringLiteral("scalars"), mcScalarsToJson(observed.mcScalars));

    return QJsonObject{
        {"present", observed.present},
        {"charge", charge},
        {"core", core},
        {"coulomb_E", coulombE},
        {"coulomb_scalars", coulombScalars},
        {"coulomb_shielding", coulombShielding},
        {"coulomb_efg_backbone", efgBackbone},
        {"coulomb_efg_aromatic", efgAromatic},
        {"mc_shielding", mcShielding},
        {"mc_category_T2", mcCategory},
        {"mc_scalars", mcScalars},
    };
}

QJsonObject csaShapeToJson(const model::CsaShape& shape) {
    QJsonObject out{{"valid", shape.valid}};
    if (!shape.valid)
        return out;
    out.insert(QStringLiteral("sigma_iso"), finiteJson(shape.sigma_iso));
    out.insert(QStringLiteral("span"), finiteJson(shape.span));
    out.insert(QStringLiteral("eta"), finiteJson(shape.eta));
    out.insert(QStringLiteral("skew"), finiteJson(shape.skew));
    out.insert(QStringLiteral("principal_values"), vec3FromEigen(shape.principal_values));
    out.insert(QStringLiteral("haeberlen_values"), vec3FromEigen(shape.haeberlen_values));
    return out;
}

QJsonObject dftShieldingToJson(const model::DftAtomShielding& shielding) {
    return QJsonObject{
        {"element", QString::fromLatin1(model::SymbolForElement(shielding.element))},
        {"total", sphericalTensorToJson(shielding.total)},
        {"diamagnetic", sphericalTensorToJson(shielding.dia)},
        {"paramagnetic", sphericalTensorToJson(shielding.para)},
        {"total_raw", mat3ToJson(shielding.total_raw)},
        {"diamagnetic_raw", mat3ToJson(shielding.dia_raw)},
        {"paramagnetic_raw", mat3ToJson(shielding.para_raw)},
        {"orca_coord_A", vec3FromEigen(shielding.orca_coord)},
    };
}

QJsonObject ringNullMeasurementToJson(const model::RingNullMeasurement& m) {
    QJsonObject out{
        {"valid", m.valid},
        {"side", QString::fromLatin1(model::RingNullSideName(m.side))},
        {"distance_A", finiteJson(m.distanceA)},
        {"axial_A", finiteJson(m.axialA)},
        {"abs_axial_A", finiteJson(m.absAxialA)},
        {"radial_A", finiteJson(m.radialA)},
        {"null_radius_A", finiteJson(m.nullRadiusA)},
        {"null_margin_A", finiteJson(m.nullMarginA)},
        {"angle_deg", finiteJson(m.angleDeg)},
        {"angular_factor", finiteJson(m.angularFactor)},
    };
    out.insert(QStringLiteral("atom_position"), vec3FromEigen(m.atomPosition));
    out.insert(QStringLiteral("ring_center"), vec3FromEigen(m.ring.center));
    out.insert(QStringLiteral("ring_normal"), vec3FromEigen(m.ring.normal));
    out.insert(QStringLiteral("ring_radius_A"), finiteJson(m.ring.radius));
    return out;
}

QJsonObject ringNullAtomIdentityToJson(const model::RingNullAtomIdentity& identity) {
    return QJsonObject{
        {"atom_label_amber", identity.atomLabelAmber},
        {"atom_label_iupac", identity.atomLabelIupac},
        {"atom_label_bmrb", identity.atomLabelBmrb},
        {"residue_index", identity.residueIndex},
        {"residue_number", identity.residueNumber},
        {"residue_label_amber", identity.residueLabelAmber},
        {"residue_label_iupac", identity.residueLabelIupac},
        {"residue_label_bmrb", identity.residueLabelBmrb},
    };
}

QJsonObject ringNullRingIdentityToJson(const model::RingNullRingIdentity& identity) {
    return QJsonObject{
        {"type_name", identity.typeName},
        {"type_index", identity.typeIndex},
        {"kind", identity.kind},
        {"parent_residue_index", identity.parentResidueIndex},
        {"parent_residue_number", identity.parentResidueNumber},
        {"parent_residue_label_amber", identity.parentResidueLabelAmber},
        {"parent_residue_label_iupac", identity.parentResidueLabelIupac},
        {"parent_residue_label_bmrb", identity.parentResidueLabelBmrb},
        {"fused_partner_ring_id", identity.fusedPartnerRingId},
        {"atom_indices", intVectorToJson(identity.atomIndices)},
    };
}

QJsonObject ringNullEventFrameToJson(const model::RingNullEventFrame& frame) {
    return QJsonObject{
        {"phase_coordinate", QStringLiteral("signed_null_margin_A")},
        {"zero_definition", QStringLiteral("linear interpolation to null_margin_A = 0")},
        {"from_signed_null_margin_A", finiteJson(frame.fromSignedNullMarginA)},
        {"to_signed_null_margin_A", finiteJson(frame.toSignedNullMarginA)},
        {"signed_null_margin_step_A", finiteJson(frame.signedNullMarginStepA)},
        {"zero_fraction", finiteJson(frame.zeroFraction)},
        {"zero_time_ps", finiteJson(frame.zeroTimePs)},
        {"zero_atom_position", vec3FromEigen(frame.zeroAtomPosition)},
    };
}

QJsonObject ringNullMotionToJson(const model::RingNullMotion& motion) {
    return QJsonObject{
        {"world_vector_A", vec3FromEigen(motion.worldVectorA)},
        {"distance_A", finiteJson(motion.distanceA)},
        {"time_step_ps", finiteJson(motion.timeStepPs)},
        {"radial_change_A", finiteJson(motion.radialChangeA)},
        {"abs_axial_change_A", finiteJson(motion.absAxialChangeA)},
        {"distance_change_A", finiteJson(motion.distanceChangeA)},
        {"angle_change_deg", finiteJson(motion.angleChangeDeg)},
        {"angular_factor_change", finiteJson(motion.angularFactorChange)},
    };
}

QJsonObject ringNullSnapshotToJson(const model::RingNullOrcaSnapshot& snapshot) {
    return QJsonObject{
        {"frame", snapshot.frameIndex},
        {"time_ps", finiteJson(snapshot.timePs)},
        {"null", ringNullMeasurementToJson(snapshot.null)},
        {"orca", dftShieldingToJson(snapshot.shielding)},
        {"shape", QJsonObject{
            {"total", csaShapeToJson(snapshot.totalShape)},
            {"diamagnetic", csaShapeToJson(snapshot.diaShape)},
            {"paramagnetic", csaShapeToJson(snapshot.paraShape)},
        }},
    };
}

QJsonObject ringNullSignalStampToJson(const model::RingNullSignalStamp& stamp) {
    QJsonObject orca{{"present", stamp.orcaPresent}};
    if (stamp.orcaPresent) {
        orca.insert(QStringLiteral("shielding"), dftShieldingToJson(stamp.orca));
        orca.insert(QStringLiteral("shape"), QJsonObject{
            {"total", csaShapeToJson(stamp.orcaTotalShape)},
            {"diamagnetic", csaShapeToJson(stamp.orcaDiaShape)},
            {"paramagnetic", csaShapeToJson(stamp.orcaParaShape)},
        });
    }
    return QJsonObject{
        {"frame", stamp.frameIndex},
        {"time_ps", finiteJson(stamp.timePs)},
        {"time_offset_from_zero_ps", finiteJson(stamp.timeOffsetFromZeroPs)},
        {"dft_ordinal_offset_from_zero", finiteJson(stamp.dftOrdinalOffsetFromZero)},
        {"null", ringNullMeasurementToJson(stamp.null)},
        {"orca", orca},
        {"snapshot_present", stamp.snapshotPresent},
        {"mopac", mopacSignalsToJson(stamp.mopac)},
    };
}

QJsonArray ringNullSignalStampsToJson(const std::vector<model::RingNullSignalStamp>& stamps) {
    QJsonArray out;
    for (const model::RingNullSignalStamp& stamp : stamps)
        out.append(ringNullSignalStampToJson(stamp));
    return out;
}

QJsonObject ringNullEntryToJson(const model::RingNullCollarEntry& entry) {
    return QJsonObject{
        {"atom", static_cast<qint64>(entry.atom)},
        {"ring", static_cast<qint64>(entry.ring)},
        {"atom_identity", ringNullAtomIdentityToJson(entry.atomIdentity)},
        {"ring_identity", ringNullRingIdentityToJson(entry.ringIdentity)},
        {"event_frame", ringNullEventFrameToJson(entry.eventFrame)},
        {"motion", ringNullMotionToJson(entry.motion)},
        {"from", ringNullSnapshotToJson(entry.from)},
        {"to", ringNullSnapshotToJson(entry.to)},
        {"signal_stamps", ringNullSignalStampsToJson(entry.signalStamps)},
    };
}

QJsonObject proteinAtomIdentityToJson(const model::QtProtein& protein, std::size_t atom) {
    QJsonObject out{
        {"atom_label_amber", protein.atomLabel(atom, model::NamingConvention::Amber)},
        {"atom_label_iupac", protein.atomLabel(atom, model::NamingConvention::Iupac)},
        {"atom_label_bmrb", protein.atomLabel(atom, model::NamingConvention::Bmrb)},
    };
    const model::QtAtom& qtAtom = protein.atom(atom);
    out.insert(QStringLiteral("residue_index"), qtAtom.residueIndex);
    if (qtAtom.residueIndex >= 0 &&
        static_cast<std::size_t>(qtAtom.residueIndex) < protein.residueCount()) {
        const std::size_t residue = static_cast<std::size_t>(qtAtom.residueIndex);
        out.insert(QStringLiteral("residue_number"),
                   protein.residue(residue).address.residueNumber);
        out.insert(QStringLiteral("residue_label_amber"),
                   protein.residueLabel(residue, model::NamingConvention::Amber,
                                        model::NamingSource::Verbatim));
        out.insert(QStringLiteral("residue_label_iupac"),
                   protein.residueLabel(residue, model::NamingConvention::Iupac,
                                        model::NamingSource::Verbatim));
        out.insert(QStringLiteral("residue_label_bmrb"),
                   protein.residueLabel(residue, model::NamingConvention::Bmrb,
                                        model::NamingSource::Verbatim));
    }
    return out;
}

QJsonObject proteinRingIdentityToJson(const model::QtProtein& protein, std::size_t ring) {
    const model::QtRing& qtRing = protein.ring(ring);
    QJsonObject out{
        {"type_name", QString::fromLatin1(qtRing.TypeName())},
        {"type_index", qtRing.TypeIndexAsInt()},
        {"kind", qtRing.IsAromatic() ? QStringLiteral("aromatic") : QStringLiteral("saturated")},
        {"parent_residue_index", qtRing.parentResidueIndex},
        {"parent_residue_number", qtRing.parentResidueNumber},
        {"fused_partner_ring_id", qtRing.fusedPartnerRingId},
    };
    QJsonArray atomIndices;
    for (int32_t atomIndex : qtRing.atomIndices)
        atomIndices.append(atomIndex);
    out.insert(QStringLiteral("atom_indices"), atomIndices);
    if (qtRing.parentResidueIndex >= 0 &&
        static_cast<std::size_t>(qtRing.parentResidueIndex) < protein.residueCount()) {
        const std::size_t residue = static_cast<std::size_t>(qtRing.parentResidueIndex);
        out.insert(QStringLiteral("parent_residue_label_amber"),
                   protein.residueLabel(residue, model::NamingConvention::Amber,
                                        model::NamingSource::Verbatim));
        out.insert(QStringLiteral("parent_residue_label_iupac"),
                   protein.residueLabel(residue, model::NamingConvention::Iupac,
                                        model::NamingSource::Verbatim));
        out.insert(QStringLiteral("parent_residue_label_bmrb"),
                   protein.residueLabel(residue, model::NamingConvention::Bmrb,
                                        model::NamingSource::Verbatim));
    }
    return out;
}

QJsonObject ringCurrentFitToJson(const model::RingCurrentLinearFit& fit) {
    return QJsonObject{
        {"valid", fit.valid},
        {"sample_count", fit.sampleCount},
        {"intercept", finiteJson(fit.intercept)},
        {"scale", finiteJson(fit.scale)},
        {"r2", finiteJson(fit.r2)},
        {"correlation", finiteJson(fit.correlation)},
        {"sse", finiteJson(fit.sse)},
        {"sst", finiteJson(fit.sst)},
        {"null_shift_count", fit.nullShiftCount},
        {"null_median_r2", finiteJson(fit.nullMedianR2)},
        {"null_max_r2", finiteJson(fit.nullMaxR2)},
        {"null_ge_real_fraction", finiteJson(fit.nullGeRealFraction)},
        {"null_r2", doubleVectorToJson(fit.nullR2)},
    };
}

QJsonObject ringCurrentSampleToJson(const model::RingCurrentFaceSample& sample) {
    return QJsonObject{
        {"frame", sample.frameIndex},
        {"time_ps", finiteJson(sample.timePs)},
        {"expected_relationship_value", finiteJson(sample.expectedRelationshipValue)},
        {"distance_only_value", finiteJson(sample.distanceOnlyValue)},
        {"angular_only_value", finiteJson(sample.angularOnlyValue)},
        {"biot_savart", sphericalTensorToJson(sample.biotSavart)},
        {"geometry", ringNullMeasurementToJson(sample.geometry)},
        {"orca", dftShieldingToJson(sample.orca)},
    };
}

QJsonArray ringCurrentSamplesToJson(const std::vector<model::RingCurrentFaceSample>& samples) {
    QJsonArray out;
    for (const model::RingCurrentFaceSample& sample : samples)
        out.append(ringCurrentSampleToJson(sample));
    return out;
}

QJsonObject ringCurrentEntryToJson(const model::QtProtein& protein,
                                   const model::RingCurrentFaceEntry& entry,
                                   bool includeSamples) {
    QJsonObject out{
        {"atom", static_cast<qint64>(entry.atom)},
        {"ring", static_cast<qint64>(entry.ring)},
        {"atom_identity", proteinAtomIdentityToJson(protein, entry.atom)},
        {"ring_identity", proteinRingIdentityToJson(protein, entry.ring)},
        {"sample_count", static_cast<int>(entry.samples.size())},
        {"hard_lobe_crossing", entry.hardLobeCrossing},
        {"positive_template_samples", entry.positiveTemplateSamples},
        {"negative_template_samples", entry.negativeTemplateSamples},
        {"template_sign_changes", entry.templateSignChanges},
        {"min_expected_relationship_value", finiteJson(entry.minTemplate)},
        {"max_expected_relationship_value", finiteJson(entry.maxTemplate)},
        {"expected_relationship_span", finiteJson(entry.templateSpan)},
        {"fits", QJsonObject{
            {"orca_total_T0", ringCurrentFitToJson(entry.orcaTotalT0)},
            {"orca_diamagnetic_T0", ringCurrentFitToJson(entry.orcaDiamagneticT0)},
            {"orca_paramagnetic_T0", ringCurrentFitToJson(entry.orcaParamagneticT0)},
            {"orca_total_T2_magnitude", ringCurrentFitToJson(entry.orcaTotalT2Magnitude)},
            {"biot_savart_T0_vs_orca_total_T0",
             ringCurrentFitToJson(entry.biotSavartOrcaTotalT0)},
        }},
        {"predictor_diagnostics", QJsonObject{
            {"expected_relationship_value_vs_biot_savart_T0",
             ringCurrentFitToJson(entry.expectedRelationshipBiotSavartT0)},
        }},
        {"confound_fits", QJsonObject{
            {"distance_only_orca_total_T0",
             ringCurrentFitToJson(entry.distanceOnlyOrcaTotalT0)},
            {"angular_only_orca_total_T0",
             ringCurrentFitToJson(entry.angularOnlyOrcaTotalT0)},
        }},
    };
    if (includeSamples)
        out.insert(QStringLiteral("samples"), ringCurrentSamplesToJson(entry.samples));
    return out;
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
    } else if (const auto* atomAnchor = std::get_if<AtomAnchor>(&anchor)) {
        out["kind"] = "atom"; out["atom"] = static_cast<qint64>(atomAnchor->atom);
    } else if (const auto* residueAnchor = std::get_if<ResidueAnchor>(&anchor)) {
        out["kind"] = "residue"; out["residue"] = static_cast<qint64>(residueAnchor->residue);
    } else if (const auto* tupleAnchor = std::get_if<AtomTupleAnchor>(&anchor)) {
        QJsonArray atoms;
        for (auto atom : tupleAnchor->atoms) atoms.append(static_cast<qint64>(atom));
        out["kind"] = "atom_tuple"; out["atoms"] = atoms;
    } else if (const auto* bondAnchor = std::get_if<BondAnchor>(&anchor)) {
        out["kind"] = "bond"; out["bond"] = static_cast<qint64>(bondAnchor->bond);
    } else if (const auto* vectorAnchor = std::get_if<BondVectorAnchor>(&anchor)) {
        out["kind"] = "bond_vector";
        out["residue"] = static_cast<qint64>(vectorAnchor->residue);
        out["kind_id"] = static_cast<qint64>(vectorAnchor->kind);
    } else if (const auto* ringAnchor = std::get_if<RingAnchor>(&anchor)) {
        out["kind"] = "ring"; out["ring"] = static_cast<qint64>(ringAnchor->ring);
    } else if (const auto* aromaticAnchor = std::get_if<AromaticRingAnchor>(&anchor)) {
        out["kind"] = "aromatic_ring"; out["ring"] = static_cast<qint64>(aromaticAnchor->ring);
    } else if (const auto* saturatedAnchor = std::get_if<SaturatedRingAnchor>(&anchor)) {
        out["kind"] = "saturated_ring"; out["ring"] = static_cast<qint64>(saturatedAnchor->ring);
    } else if (const auto* contributionAnchor = std::get_if<RingContributionPairAnchor>(&anchor)) {
        out["kind"] = "ring_contribution_pair";
        out["pair"] = static_cast<qint64>(contributionAnchor->pair);
    } else if (const auto* membershipAnchor = std::get_if<RingMembershipAnchor>(&anchor)) {
        out["kind"] = "ring_membership";
        out["membership"] = static_cast<qint64>(membershipAnchor->membership);
    } else if (const auto* mutationAnchor = std::get_if<MutationMatchPairAnchor>(&anchor)) {
        out["kind"] = "mutation_match_pair"; out["pair"] = static_cast<qint64>(mutationAnchor->pair);
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
    const QJsonArray modeArray = modesValue.toArray();
    for (const auto& value : std::as_const(modeArray)) {
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
    QHttpServerConfiguration config = server_->configuration();
    config.setKeepAliveTimeout(std::chrono::seconds(300));
    server_->setConfiguration(config);
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

    server_->route(QStringLiteral("/api/interface"), [this]() {
        ASSERT_THREAD(this);
        return jsonResponse(restInterfaceDescription());
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
        for (const auto& v : std::as_const(arr)) {
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
    // POST /overlay {"name": "ribbon"|"rings"|"butterfly"|"nullcone"|"bfield"
    //                         |"shadow", "visible": bool}
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
                               "(ribbon|rings|butterfly|nullcone|bfield|shadow)").arg(name),
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
    // ---- field-grid isosurface opacity (live tuning) --------------------
    //
    // POST /field/opacity {"opacity": number 0..1}
    // Sets the butterfly isosurface opacity independently from the null cone.
    // This is a resthero tuning knob for making the field a thin context layer.
    server_->route(QStringLiteral("/field/opacity"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->fieldGridOverlay())
            return errorResponse(QStringLiteral("field-grid overlay not available"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains(QStringLiteral("opacity")) ||
            !body.value(QStringLiteral("opacity")).isDouble()) {
            return errorResponse(QStringLiteral("body must be {\"opacity\": number}"),
                                 SC::BadRequest);
        }
        const double opacity = body.value(QStringLiteral("opacity")).toDouble();
        if (opacity < 0.0 || opacity > 1.0)
            return errorResponse(QStringLiteral("opacity must be between 0 and 1"),
                                 SC::BadRequest);
        scene_->fieldGridOverlay()->setOpacity(opacity);
        scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return QHttpServerResponse(SC::NoContent);
    });

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

    // ---- field null cone: magic-angle sign boundary ---------------------
    //
    // POST /field/null_cone {"visible"?: bool, "opacity"?: 0..1,
    //                         "length"?: Angstrom}
    // The null cone is lightweight geometry on the same ring-field overlay:
    // the dipolar magic-angle boundary where the ring-current sign flips.
    // It is intentionally independent of the expensive butterfly isosurface,
    // but shares selected-ring filtering.
    server_->route(QStringLiteral("/field/null_cone"), Method::Post,
                   [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->fieldGridOverlay())
            return errorResponse(QStringLiteral("field-grid overlay not available"),
                                 SC::ServiceUnavailable);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);

        bool applied = false;
        QtFieldGridOverlay* field = scene_->fieldGridOverlay();
        if (body.value(QStringLiteral("visible")).isBool()) {
            field->setNullConeVisible(body.value(QStringLiteral("visible")).toBool());
            applied = true;
        }
        if (body.value(QStringLiteral("opacity")).isDouble()) {
            const double opacity = body.value(QStringLiteral("opacity")).toDouble();
            if (opacity < 0.0 || opacity > 1.0)
                return errorResponse(QStringLiteral("opacity must be between 0 and 1"),
                                     SC::BadRequest);
            field->setNullConeOpacity(opacity);
            applied = true;
        }
        if (body.value(QStringLiteral("length")).isDouble()) {
            const double length = body.value(QStringLiteral("length")).toDouble();
            if (length <= 0.0)
                return errorResponse(QStringLiteral("length must be > 0"),
                                     SC::BadRequest);
            field->setNullConeLength(length);
            applied = true;
        }
        if (!applied)
            return errorResponse(
                QStringLiteral("body must set visible, opacity, and/or length"),
                SC::BadRequest);

        scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return QHttpServerResponse(SC::NoContent);
    });

    // ---- ring null collar collection ------------------------------------
    //
    // POST /api/ring/null_crossings {"atom"?, "ring"?, "start_frame"?, "end_frame"?}
    // Runs an explicit collar object over adjacent DFT frames and returns the
    // whole before/after per-atom ORCA shielding records for crossings.
    auto ringNullCrossingsHandler = [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        model::Conformation* conf = loaded_ ? loaded_->conformation.get() : nullptr;
        if (!protein || !conf)
            return errorResponse(QStringLiteral("protein/conformation not wired"),
                                 SC::ServiceUnavailable);
        if (!loaded_->manifest.dft)
            return errorResponse(QStringLiteral("calcset has no DFT frame manifest"),
                                 SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);

        model::RingNullCollarOptions options;
        if (body.contains(QStringLiteral("atom"))) {
            const qint64 raw = body.value(QStringLiteral("atom")).toInteger(-1);
            if (raw < 0 || static_cast<std::size_t>(raw) >= protein->atomCount())
                return errorResponse(QStringLiteral("atom out of range"), SC::BadRequest);
            options.atom = static_cast<std::size_t>(raw);
        }

        if (body.contains(QStringLiteral("ring"))) {
            const qint64 raw = body.value(QStringLiteral("ring")).toInteger(-1);
            if (raw < 0 || static_cast<std::size_t>(raw) >= protein->ringCount())
                return errorResponse(QStringLiteral("ring out of range"), SC::BadRequest);
            options.ring = static_cast<std::size_t>(raw);
        }

        const int frameCount = static_cast<int>(conf->frameCount());
        if (frameCount <= 0)
            return errorResponse(QStringLiteral("conformation has no frames"),
                                 SC::ServiceUnavailable);
        if (body.contains(QStringLiteral("start_frame")))
            options.startFrame = body.value(QStringLiteral("start_frame")).toInt(-1);
        if (body.contains(QStringLiteral("end_frame")))
            options.endFrame = body.value(QStringLiteral("end_frame")).toInt(-1);
        if (body.contains(QStringLiteral("surface_tolerance_A"))) {
            options.surfaceToleranceA =
                body.value(QStringLiteral("surface_tolerance_A")).toDouble(-1.0);
        } else if (body.contains(QStringLiteral("tolerance_A"))) {
            options.surfaceToleranceA =
                body.value(QStringLiteral("tolerance_A")).toDouble(-1.0);
        }

        if (!std::isfinite(options.surfaceToleranceA) || options.surfaceToleranceA < 0.0)
            return errorResponse(QStringLiteral("surface_tolerance_A must be finite and >= 0"),
                                 SC::BadRequest);
        options.includeSaturatedRings =
            body.value(QStringLiteral("include_saturated")).toBool(false);
        options.includeSignalStamps =
            body.value(QStringLiteral("include_signal_stamps")).toBool(true);
        options.stampRadiusDft =
            body.value(QStringLiteral("stamp_radius_dft")).toInt(2);
        if (options.stampRadiusDft < 0 || options.stampRadiusDft > 20)
            return errorResponse(QStringLiteral("stamp_radius_dft must be between 0 and 20"),
                                 SC::BadRequest);
        const bool explicitFullScan = body.value(QStringLiteral("full_scan")).toBool(false);
        if (!options.atom && !options.ring && !explicitFullScan) {
            return errorResponse(
                QStringLiteral("body must set atom and/or ring, or set full_scan=true"),
                SC::BadRequest);
        }

        model::RingNullCollar collar(options);
        QString collectError;
        if (!collar.collect(*protein, *conf, loaded_->manifest.dft->frames, &collectError))
            return errorResponse(collectError, SC::BadRequest);

        QJsonArray entries;
        for (const model::RingNullCollarEntry& entry : collar.entries())
            entries.append(ringNullEntryToJson(entry));

        const model::RingNullCollarSummary& summary = collar.summary();

        QJsonObject out{
            {"kind", QStringLiteral("ring_null_collar_collection")},
            {"collar", QJsonObject{
                {"mode", QStringLiteral("zero_width_surface_crossing")},
                {"physical_width_A", 0.0},
                {"surface_tolerance_A", options.surfaceToleranceA},
                {"phase_coordinate", QStringLiteral("signed_null_margin_A")},
                {"phase_zero", QStringLiteral("radial_A - sqrt(2) * abs(axial_A) = 0")},
                {"width_semantics",
                 QStringLiteral("surface_tolerance_A is numerical classification tolerance; it is not a finite physical aperture")},
            }},
            {"statistic", QJsonObject{
                {"name", QStringLiteral("ring_null_margin_A")},
                {"unit", QStringLiteral("A")},
                {"definition", QStringLiteral("radial_A - sqrt(2) * abs(axial_A)")},
                {"crossing_definition", QStringLiteral("sign change between adjacent DFT frames")},
                {"magic_angle_deg", model::RingNullMagicAngleDegrees()},
            }},
            {"options", QJsonObject{
                {"atom", options.atom ? QJsonValue(static_cast<qint64>(*options.atom))
                                      : QJsonValue(QJsonValue::Null)},
                {"ring", options.ring ? QJsonValue(static_cast<qint64>(*options.ring))
                                      : QJsonValue(QJsonValue::Null)},
                {"start_frame", options.startFrame ? QJsonValue(*options.startFrame)
                                                   : QJsonValue(QJsonValue::Null)},
                {"end_frame", options.endFrame ? QJsonValue(*options.endFrame)
                                               : QJsonValue(QJsonValue::Null)},
                {"surface_tolerance_A", options.surfaceToleranceA},
                {"include_saturated", options.includeSaturatedRings},
                {"full_scan", explicitFullScan},
                {"include_signal_stamps", options.includeSignalStamps},
                {"stamp_radius_dft", options.stampRadiusDft},
            }},
            {"summary", QJsonObject{
                {"complete", summary.complete},
                {"dft_frames_declared", summary.dftFramesDeclared},
                {"dft_frames_loaded", summary.dftFramesLoaded},
                {"dft_frames_skipped", summary.dftFramesSkipped},
                {"dft_pairs_scanned", summary.dftPairsScanned},
                {"atoms_scanned", summary.atomsScanned},
                {"rings_scanned", summary.ringsScanned},
                {"entry_count", summary.entryCount},
                {"signal_stamp_count", summary.signalStampCount},
            }},
            {"entries", entries},
        };
        return jsonResponse(out);
    };
    server_->route(QStringLiteral("/api/ring/null_crossings"), Method::Post,
                   ringNullCrossingsHandler);

    // ---- ring-current receiver: stash paths, fit ORCA -------------------
    //
    // POST /api/ring/current_face_collar {"atom"?, "ring"?, "start_frame"?,
    //                                     "end_frame"?, "min_samples"?}
    // Collects atom/ring paths whose expected ring-current relationship value
    // samples both lobes, then evaluates the stash as:
    // ORCA_component = intercept + scale * expected_relationship_value.
    auto ringCurrentFaceCollarHandler = [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        model::Conformation* conf = loaded_ ? loaded_->conformation.get() : nullptr;
        if (!protein || !conf)
            return errorResponse(QStringLiteral("protein/conformation not wired"),
                                 SC::ServiceUnavailable);
        if (!loaded_->manifest.dft)
            return errorResponse(QStringLiteral("calcset has no DFT frame manifest"),
                                 SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);

        model::RingCurrentFaceCollarOptions options;
        if (body.contains(QStringLiteral("atom"))) {
            const qint64 raw = body.value(QStringLiteral("atom")).toInteger(-1);
            if (raw < 0 || static_cast<std::size_t>(raw) >= protein->atomCount())
                return errorResponse(QStringLiteral("atom out of range"), SC::BadRequest);
            options.atom = static_cast<std::size_t>(raw);
        }
        if (body.contains(QStringLiteral("ring"))) {
            const qint64 raw = body.value(QStringLiteral("ring")).toInteger(-1);
            if (raw < 0 || static_cast<std::size_t>(raw) >= protein->ringCount())
                return errorResponse(QStringLiteral("ring out of range"), SC::BadRequest);
            options.ring = static_cast<std::size_t>(raw);
        }

        const int frameCount = static_cast<int>(conf->frameCount());
        if (frameCount <= 0)
            return errorResponse(QStringLiteral("conformation has no frames"),
                                 SC::ServiceUnavailable);
        if (body.contains(QStringLiteral("start_frame")))
            options.startFrame = body.value(QStringLiteral("start_frame")).toInt(-1);
        if (body.contains(QStringLiteral("end_frame")))
            options.endFrame = body.value(QStringLiteral("end_frame")).toInt(-1);
        if (body.contains(QStringLiteral("surface_tolerance_A"))) {
            options.surfaceToleranceA =
                body.value(QStringLiteral("surface_tolerance_A")).toDouble(-1.0);
        } else if (body.contains(QStringLiteral("tolerance_A"))) {
            options.surfaceToleranceA =
                body.value(QStringLiteral("tolerance_A")).toDouble(-1.0);
        }
        if (!std::isfinite(options.surfaceToleranceA) || options.surfaceToleranceA < 0.0)
            return errorResponse(QStringLiteral("surface_tolerance_A must be finite and >= 0"),
                                 SC::BadRequest);

        if (body.contains(QStringLiteral("template_zero_tolerance")))
            options.templateZeroTolerance =
                body.value(QStringLiteral("template_zero_tolerance")).toDouble(-1.0);
        if (!std::isfinite(options.templateZeroTolerance) ||
            options.templateZeroTolerance < 0.0) {
            return errorResponse(QStringLiteral("template_zero_tolerance must be finite and >= 0"),
                                 SC::BadRequest);
        }

        options.includeSaturatedRings =
            body.value(QStringLiteral("include_saturated")).toBool(false);
        options.minSamples = body.value(QStringLiteral("min_samples")).toInt(6);
        options.minSamplesPerLobe =
            body.value(QStringLiteral("min_samples_per_lobe")).toInt(3);
        options.minExpectedRelationshipSpan =
            body.value(QStringLiteral("min_expected_relationship_span")).toDouble(0.02);
        options.minAbsLobeExpectedValue =
            body.value(QStringLiteral("min_abs_lobe_expected_value")).toDouble(0.005);
        options.maxEntries = body.value(QStringLiteral("max_entries")).toInt(25);
        options.nullShiftCount = body.value(QStringLiteral("null_shift_count")).toInt(64);
        if (options.minSamples < 3 || options.minSamples > 1000)
            return errorResponse(QStringLiteral("min_samples must be between 3 and 1000"),
                                 SC::BadRequest);
        if (options.minSamplesPerLobe < 1 || options.minSamplesPerLobe > options.minSamples)
            return errorResponse(QStringLiteral("min_samples_per_lobe must be between 1 and min_samples"),
                                 SC::BadRequest);
        if (!std::isfinite(options.minExpectedRelationshipSpan) ||
            options.minExpectedRelationshipSpan < 0.0) {
            return errorResponse(QStringLiteral("min_expected_relationship_span must be finite and >= 0"),
                                 SC::BadRequest);
        }
        if (!std::isfinite(options.minAbsLobeExpectedValue) ||
            options.minAbsLobeExpectedValue < 0.0) {
            return errorResponse(QStringLiteral("min_abs_lobe_expected_value must be finite and >= 0"),
                                 SC::BadRequest);
        }
        if (options.maxEntries < 0 || options.maxEntries > 500)
            return errorResponse(QStringLiteral("max_entries must be between 0 and 500"),
                                 SC::BadRequest);
        if (options.nullShiftCount < 0 || options.nullShiftCount > 500)
            return errorResponse(QStringLiteral("null_shift_count must be between 0 and 500"),
                                 SC::BadRequest);

        const bool explicitFullScan = body.value(QStringLiteral("full_scan")).toBool(false);
        if (!options.atom && !options.ring && !explicitFullScan) {
            return errorResponse(
                QStringLiteral("body must set atom and/or ring, or set full_scan=true"),
                SC::BadRequest);
        }
        const bool includeSamples = body.value(QStringLiteral("include_samples")).toBool(true);

        model::RingCurrentFaceCollar collar(options);
        QString collectError;
        if (!collar.collect(*protein, *conf, loaded_->manifest.dft->frames, &collectError))
            return errorResponse(collectError, SC::BadRequest);

        QJsonArray entries;
        for (const model::RingCurrentFaceEntry& entry : collar.entries())
            entries.append(ringCurrentEntryToJson(*protein, entry, includeSamples));

        const model::RingCurrentFaceCollarSummary& summary = collar.summary();
        QJsonObject out{
            {"kind", QStringLiteral("ring_current_face_collar")},
            {"receiver", QJsonObject{
                {"collector", QStringLiteral("stash atom/ring DFT paths only when expected_relationship_value crosses both lobes")},
                {"expected_relationship_value",
                 QStringLiteral("(3*cos(theta)^2 - 1) / distance_A^3")},
                {"fit_model",
                 QStringLiteral("ORCA_component = intercept + scale * expected_relationship_value + residual")},
                {"biot_savart_fit_model",
                 QStringLiteral("ORCA_total_T0 = intercept + scale * recomputed_BS_T0 + residual")},
                {"biot_savart_relationship",
                 QStringLiteral("recomputed_BS_T0 is the finite Johnson-Bovey/Biot-Savart version of the same ring-current geometry; in the far field it tracks expected_relationship_value")},
                {"predictor_diagnostic_model",
                 QStringLiteral("recomputed_BS_T0 = intercept + scale * expected_relationship_value + residual")},
                {"biot_savart_source",
                 QStringLiteral("recomputed from current ring geometry via QtBiotSavartCalc and ring literature intensity; does not read ring_contributions.npy")},
                {"hard_crossing_requirement",
                 QStringLiteral("positive and negative expected_relationship_value samples with at least one sign change")},
                {"null_model",
                 QStringLiteral("circular shifts of the ORCA trace against the same expected relationship values")},
            }},
            {"options", QJsonObject{
                {"atom", options.atom ? QJsonValue(static_cast<qint64>(*options.atom))
                                      : QJsonValue(QJsonValue::Null)},
                {"ring", options.ring ? QJsonValue(static_cast<qint64>(*options.ring))
                                      : QJsonValue(QJsonValue::Null)},
                {"start_frame", options.startFrame ? QJsonValue(*options.startFrame)
                                                   : QJsonValue(QJsonValue::Null)},
                {"end_frame", options.endFrame ? QJsonValue(*options.endFrame)
                                               : QJsonValue(QJsonValue::Null)},
                {"surface_tolerance_A", options.surfaceToleranceA},
                {"template_zero_tolerance", options.templateZeroTolerance},
                {"include_saturated", options.includeSaturatedRings},
                {"full_scan", explicitFullScan},
                {"min_samples", options.minSamples},
                {"min_samples_per_lobe", options.minSamplesPerLobe},
                {"min_expected_relationship_span", options.minExpectedRelationshipSpan},
                {"min_abs_lobe_expected_value", options.minAbsLobeExpectedValue},
                {"max_entries", options.maxEntries},
                {"null_shift_count", options.nullShiftCount},
                {"include_samples", includeSamples},
            }},
            {"summary", QJsonObject{
                {"complete", summary.complete},
                {"dft_frames_declared", summary.dftFramesDeclared},
                {"dft_frames_loaded", summary.dftFramesLoaded},
                {"dft_frames_skipped", summary.dftFramesSkipped},
                {"atoms_scanned", summary.atomsScanned},
                {"rings_scanned", summary.ringsScanned},
                {"paths_considered", summary.pathsConsidered},
                {"paths_rejected_for_samples", summary.pathsRejectedForSamples},
                {"paths_rejected_for_hard_crossing", summary.pathsRejectedForHardCrossing},
                {"paths_rejected_for_weak_lobes", summary.pathsRejectedForWeakLobes},
                {"entry_count", summary.entryCount},
                {"truncated_by_max_entries", summary.truncatedByMaxEntries},
            }},
            {"entries", entries},
        };
        return jsonResponse(out);
    };
    server_->route(QStringLiteral("/api/ring/current_face_collar"), Method::Post,
                   ringCurrentFaceCollarHandler);

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
        for (const auto& v : std::as_const(arr))
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
                for (const auto& v : std::as_const(arr)) {
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
        for (const auto& v : std::as_const(arr))
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
            for (const auto& v : std::as_const(arr)) {
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
        scene_->cameraComposer()->setMode(result.mode,
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

    // ---- resthero: transient molecule styling ---------------------------
    //
    // POST /resthero/molecule_style
    //   {preset:"scaffold"|"sticks", atom_radius_scale?, bond_radius?,
    //    atom_color?, bond_color?, render_atoms?, render_bonds?, ...}
    //
    // Figure work often needs the molecule to recede so tensor/angle geometry
    // can read. This is deliberately transient and restored by /resthero/clear.
    auto restheroMoleculeStyleHandler = [this](const QHttpServerRequest& request) {
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
    };
    server_->route(QStringLiteral("/resthero/molecule_style"), Method::Post,
                   restheroMoleculeStyleHandler);

    // ---- resthero: isolate one aromatic ring field ----------------------
    //
    // POST /resthero/ring_field {"ring": N|null}
    // Narrows the butterfly isosurface to one ring for figure work. This is
    // deliberately resthero-only: the normal UI toggle still owns whether the
    // butterfly overlay is visible, and /resthero/clear restores the previous
    // all-rings/single-ring setting.
    auto restheroRingFieldHandler = [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !scene_->fieldGridOverlay())
            return errorResponse(QStringLiteral("field-grid overlay not available"),
                                 SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);
        if (!body.contains(QStringLiteral("ring")))
            return errorResponse(QStringLiteral("body must be {\"ring\": number|null}"),
                                 SC::BadRequest);

        QtFieldGridOverlay* field = scene_->fieldGridOverlay();
        std::optional<std::size_t> ring;
        const QJsonValue ringValue = body.value(QStringLiteral("ring"));
        if (ringValue.isNull()) {
            ring = std::nullopt;
        } else if (ringValue.isDouble()) {
            const int requested = ringValue.toInt(-1);
            if (requested < 0 ||
                static_cast<std::size_t>(requested) >= field->ringCount()) {
                return errorResponse(QStringLiteral("ring out of range"),
                                     SC::BadRequest);
            }
            ring = static_cast<std::size_t>(requested);
        } else {
            return errorResponse(QStringLiteral("ring must be a number or null"),
                                 SC::BadRequest);
        }

        if (!heroshotFieldRingWasSet_) {
            heroshotFieldRingBefore_ = field->visibleRing();
            heroshotFieldRingWasSet_ = true;
        }
        field->setVisibleRing(ring);
        scene_->requestRender(MoleculeScene::RenderSource::Rest);

        QJsonObject out{
            {"ring_count", static_cast<qint64>(field->ringCount())},
            {"will_restore_on_clear", true},
        };
        out.insert(QStringLiteral("ring"),
                   ring ? QJsonValue(static_cast<qint64>(*ring)) : QJsonValue());
        return jsonResponse(out);
    };
    server_->route(QStringLiteral("/resthero/ring_field"), Method::Post,
                   restheroRingFieldHandler);

    // ---- resthero: high-resolution ring-current butterfly --------------
    //
    // POST /resthero/butterfly {"ring": N, "frame"?, "dim"?,
    //                           "threshold_ppm"?, "opacity"?,
    //                           "extent"?, "peak"?, "mode"?}
    // A transient figure/export layer, separate from the normal
    // QtFieldGridOverlay used by the UI and playback. It samples the same
    // closed-form BS/HM scalar field more densely for one selected ring and
    // contours the same signed T0 isovalues. This is a numerical approximation
    // knob for resthero images, not a smoothing filter and not a mainline
    // butterfly change.
    auto restheroButterflyHandler = [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        const auto* protein = loaded_ ? loaded_->protein.get() : nullptr;
        model::Conformation* loadedConf = loaded_ ? loaded_->conformation.get() : nullptr;
        model::Conformation* conf = transformed_
            ? static_cast<model::Conformation*>(transformed_.data())
            : loadedConf;
        if (!scene_ || !protein || !conf)
            return errorResponse(QStringLiteral("scene / protein not wired"),
                                 SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok || !body.contains(QStringLiteral("ring")))
            return errorResponse(QStringLiteral("body must be {\"ring\": int, ...}"),
                                 SC::BadRequest);
        const int ringRaw = body.value(QStringLiteral("ring")).toInt(-1);
        if (ringRaw < 0 || static_cast<std::size_t>(ringRaw) >= protein->ringCount())
            return errorResponse(QStringLiteral("ring out of range"), SC::BadRequest);
        const std::size_t ring = static_cast<std::size_t>(ringRaw);

        const int frameCount = static_cast<int>(conf->frameCount());
        const int currentFrame = playback_ ? playback_->currentFrame() : 0;
        const int frame = std::clamp(
            body.contains(QStringLiteral("frame"))
                ? body.value(QStringLiteral("frame")).toInt(currentFrame)
                : currentFrame,
            0, std::max(0, frameCount - 1));

        auto styledDouble = [&](const QString& key, double fallback, double lo, double hi) {
            const QJsonValue v = body.value(key);
            if (!v.isDouble()) return fallback;
            return std::clamp(v.toDouble(), lo, hi);
        };

        HeroshotButterflyOverlay::Style style;
        style.gridDim = std::clamp(body.value(QStringLiteral("dim")).toInt(style.gridDim),
                                   16, 72);
        style.thresholdPpm = styledDouble(QStringLiteral("threshold_ppm"),
                                          styledDouble(QStringLiteral("ppm"),
                                                       style.thresholdPpm, 0.0, 1000.0),
                                          0.0, 1000.0);
        style.opacity = styledDouble(QStringLiteral("opacity"), style.opacity, 0.0, 1.0);
        style.gaussianExtentA = styledDouble(QStringLiteral("extent"),
                                             style.gaussianExtentA, 0.1, 30.0);
        style.gaussianPeak = styledDouble(QStringLiteral("peak"),
                                          style.gaussianPeak, 0.0, 1000.0);
        style.showShielded = body.value(QStringLiteral("show_shielded")).isBool()
            ? body.value(QStringLiteral("show_shielded")).toBool()
            : style.showShielded;
        style.showDeshielded = body.value(QStringLiteral("show_deshielded")).isBool()
            ? body.value(QStringLiteral("show_deshielded")).toBool()
            : style.showDeshielded;

        const QString mode = body.value(QStringLiteral("mode"))
                                 .toString(QStringLiteral("biot_savart")).toLower();
        if (mode == QStringLiteral("biot_savart") || mode == QStringLiteral("bs")) {
            style.mode = HeroshotButterflyOverlay::Mode::BiotSavart;
        } else if (mode == QStringLiteral("haigh_mallion") ||
                   mode == QStringLiteral("hm")) {
            style.mode = HeroshotButterflyOverlay::Mode::HaighMallion;
        } else if (mode == QStringLiteral("sum")) {
            style.mode = HeroshotButterflyOverlay::Mode::Sum;
        } else {
            return errorResponse(QStringLiteral("mode must be biot_savart, haigh_mallion, or sum"),
                                 SC::BadRequest);
        }

        heroshotButterfly_ = std::make_unique<HeroshotButterflyOverlay>(
            vtkSmartPointer<vtkRenderer>(scene_->Renderer()));
        if (!heroshotButterfly_->show(*protein, *conf, ring,
                                      static_cast<std::size_t>(frame), style)) {
            heroshotButterfly_.reset();
            return errorResponse(QStringLiteral("could not build resthero butterfly"),
                                 SC::Conflict);
        }
        const HeroshotButterflyOverlay::Stats stats = heroshotButterfly_->stats();
        scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return jsonResponse(QJsonObject{
            {"ring", static_cast<qint64>(ring)},
            {"frame", frame},
            {"mode", mode},
            {"grid_dim", style.gridDim},
            {"samples", static_cast<qint64>(style.gridDim * style.gridDim * style.gridDim)},
            {"threshold_ppm", finiteJson(style.thresholdPpm)},
            {"opacity", finiteJson(style.opacity)},
            {"extent", finiteJson(style.gaussianExtentA)},
            {"peak", finiteJson(style.gaussianPeak)},
            {"min_T0", finiteJson(stats.minT0)},
            {"max_T0", finiteJson(stats.maxT0)},
            {"shielded_points", static_cast<qint64>(stats.shieldedPoints)},
            {"shielded_cells", static_cast<qint64>(stats.shieldedCells)},
            {"deshielded_points", static_cast<qint64>(stats.deshieldedPoints)},
            {"deshielded_cells", static_cast<qint64>(stats.deshieldedCells)},
            {"will_clear_on_resthero_clear", true},
        });
    };
    server_->route(QStringLiteral("/resthero/butterfly"), Method::Post,
                   restheroButterflyHandler);

    // ---- resthero: atom track ------------------------------------------
    //
    // POST /resthero/atom_track
    //   {atom?, ring?, color_by?, color_mode?,
    //    frame_source:"dft"|"range", coordinate_space:"source_ring_local"|"world",
    //    frames? | start_frame?/end_frame?/step?/max_points?}
    //
    // Draws a sampled atom path as luminous signal marks plus optional hairline
    // connectors. This
    // is deliberately a distinct object from the selected-atom trajectory
    // envelope: the track says "these are the sampled positions", not "this is a
    // density volume". It also deliberately avoids molecular primitives: the
    // marks are screen-space points/halos, not atoms, and the connections are
    // lines, not bonds. With ring_current coloring, each mark is colored by the
    // signed analytic ring-current expected relationship for atom->ring at
    // that frame:
    // (3*cos(theta)^2 - 1) / r^3. The scene stays geometry-only; the response
    // echoes the sample table for audit and panels/scripts can carry text.
    //
    // coordinate_space="source_ring_local" answers the resthero question:
    // for each DFT graph sample, if the source ring were held stationary in
    // reference_frame, where would this atom's relative position land?
    auto restheroAtomTrackHandler = [this](const QHttpServerRequest& req) {
        ASSERT_THREAD(this);
        if (!scene_ || !loaded_ || !loaded_->protein)
            return errorResponse(QStringLiteral("scene / protein not wired"),
                                 SC::ServiceUnavailable);
        const model::Conformation* conf = transformed_
            ? static_cast<const model::Conformation*>(transformed_.data())
            : loaded_->conformation.get();
        if (!conf)
            return errorResponse(QStringLiteral("conformation not wired"),
                                 SC::ServiceUnavailable);

        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        if (!ok)
            return errorResponse(QStringLiteral("invalid JSON body"), SC::BadRequest);

        const model::QtProtein& protein = *loaded_->protein;
        std::size_t atom = 0;
        if (body.contains(QStringLiteral("atom"))) {
            const qint64 rawAtom = body.value(QStringLiteral("atom")).toInteger(-1);
            if (rawAtom < 0 || static_cast<std::size_t>(rawAtom) >= protein.atomCount())
                return errorResponse(QStringLiteral("atom out of range"), SC::BadRequest);
            atom = static_cast<std::size_t>(rawAtom);
        } else if (selection_ && selection_->hasFocus()) {
            atom = selection_->focus();
        } else {
            return errorResponse(QStringLiteral("no atom: pass {\"atom\":N} or focus one"),
                                 SC::BadRequest);
        }

        const int frameCount = static_cast<int>(conf->frameCount());
        if (frameCount <= 0)
            return errorResponse(QStringLiteral("no frames available"), SC::ServiceUnavailable);
        const int liveFrame = playback_ ? playback_->currentFrame() : 0;
        const int currentFrame = std::clamp(
            body.contains(QStringLiteral("current_frame"))
                ? body.value(QStringLiteral("current_frame")).toInt(liveFrame)
                : liveFrame,
            0, frameCount - 1);

        std::optional<std::size_t> ring;
        if (body.contains(QStringLiteral("ring"))) {
            const QJsonValue rv = body.value(QStringLiteral("ring"));
            if (!rv.isNull()) {
                const qint64 rawRing = rv.toInteger(-1);
                if (rawRing < 0 || static_cast<std::size_t>(rawRing) >= protein.ringCount())
                    return errorResponse(QStringLiteral("ring out of range"), SC::BadRequest);
                ring = static_cast<std::size_t>(rawRing);
            }
        }

        const QString colorBy = body.value(QStringLiteral("color_by")).toString(
            ring ? QStringLiteral("ring_current") : QStringLiteral("none")).toLower();
        if (colorBy != QStringLiteral("ring_current")
            && colorBy != QStringLiteral("null_margin")
            && colorBy != QStringLiteral("none")) {
            return errorResponse(
                QStringLiteral("color_by must be ring_current, null_margin, or none"),
                SC::BadRequest);
        }
        if ((colorBy == QStringLiteral("ring_current")
             || colorBy == QStringLiteral("null_margin"))
            && !ring) {
            return errorResponse(QStringLiteral("color_by requires {\"ring\": N}"),
                                 SC::BadRequest);
        }

        QString frameSource = body.value(QStringLiteral("frame_source"))
                                  .toString(QStringLiteral("range"))
                                  .toLower();
        if (frameSource != QStringLiteral("range") &&
            frameSource != QStringLiteral("dft")) {
            return errorResponse(QStringLiteral("frame_source must be range or dft"),
                                 SC::BadRequest);
        }

        const QJsonArray explicitFrames =
            body.value(QStringLiteral("frames")).isArray()
                ? body.value(QStringLiteral("frames")).toArray()
                : QJsonArray{};
        if (!explicitFrames.isEmpty())
            frameSource = QStringLiteral("explicit");

        const int maxPointDefault =
            frameSource == QStringLiteral("dft") ? 10000 : 120;
        const int maxPoints = std::clamp(
            body.value(QStringLiteral("max_points")).toInt(maxPointDefault),
            1, 10000);
        std::vector<int> frames;
        int candidateFrameCount = 0;
        bool decimatedByMaxPoints = false;
        if (!explicitFrames.isEmpty()) {
            frames.reserve(static_cast<std::size_t>(
                std::min(static_cast<int>(explicitFrames.size()), maxPoints)));
            for (const auto& v : std::as_const(explicitFrames)) {
                if (static_cast<int>(frames.size()) >= maxPoints)
                    break;
                const int f = v.toInt(-1);
                if (f < 0 || f >= frameCount)
                    return errorResponse(QStringLiteral("frame out of range"), SC::BadRequest);
                frames.push_back(f);
            }
            candidateFrameCount = static_cast<int>(explicitFrames.size());
        } else if (frameSource == QStringLiteral("dft")) {
            if (!loaded_->manifest.dft)
                return errorResponse(QStringLiteral("calcset has no DFT frame manifest"),
                                     SC::ServiceUnavailable);
            const int startFrame = std::clamp(
                body.contains(QStringLiteral("start_frame"))
                    ? body.value(QStringLiteral("start_frame")).toInt(0)
                    : 0,
                0, frameCount - 1);
            const int endFrame = std::clamp(
                body.contains(QStringLiteral("end_frame"))
                    ? body.value(QStringLiteral("end_frame")).toInt(frameCount - 1)
                    : frameCount - 1,
                0, frameCount - 1);
            const int lo = std::min(startFrame, endFrame);
            const int hi = std::max(startFrame, endFrame);
            std::vector<int> dftFrames;
            dftFrames.reserve(loaded_->manifest.dft->frames.size());
            for (const h5reader::io::DftFrame& declared : loaded_->manifest.dft->frames) {
                const int f = declared.frame_index;
                if (f < 0 || f >= frameCount)
                    continue;
                if (f < lo || f > hi)
                    continue;
                dftFrames.push_back(f);
            }
            std::sort(dftFrames.begin(), dftFrames.end());
            dftFrames.erase(std::unique(dftFrames.begin(), dftFrames.end()),
                            dftFrames.end());
            candidateFrameCount = static_cast<int>(dftFrames.size());
            if (candidateFrameCount > maxPoints) {
                decimatedByMaxPoints = true;
                const double stride = static_cast<double>(candidateFrameCount)
                                    / static_cast<double>(maxPoints);
                for (int i = 0; i < maxPoints; ++i) {
                    const int idx = std::clamp(
                        static_cast<int>(std::floor(static_cast<double>(i) * stride)),
                        0, candidateFrameCount - 1);
                    frames.push_back(dftFrames[static_cast<std::size_t>(idx)]);
                }
            } else {
                frames = std::move(dftFrames);
            }
        } else {
            const int startFrame = std::clamp(
                body.contains(QStringLiteral("start_frame"))
                    ? body.value(QStringLiteral("start_frame")).toInt(0)
                    : 0,
                0, frameCount - 1);
            const int endFrame = std::clamp(
                body.contains(QStringLiteral("end_frame"))
                    ? body.value(QStringLiteral("end_frame")).toInt(frameCount - 1)
                    : frameCount - 1,
                0, frameCount - 1);
            const int lo = std::min(startFrame, endFrame);
            const int hi = std::max(startFrame, endFrame);
            const int span = hi - lo + 1;
            const int requestedStep = body.value(QStringLiteral("step")).toInt(0);
            const int step = requestedStep > 0
                ? std::clamp(requestedStep, 1, frameCount)
                : std::max(1, static_cast<int>(
                      std::ceil(static_cast<double>(span) / static_cast<double>(maxPoints))));
            for (int f = lo; f <= hi && static_cast<int>(frames.size()) < maxPoints; f += step)
                frames.push_back(f);

            const bool includeCurrent =
                body.contains(QStringLiteral("include_current"))
                    ? body.value(QStringLiteral("include_current")).toBool(true)
                    : true;
            if (includeCurrent && currentFrame >= lo && currentFrame <= hi
                && std::find(frames.begin(), frames.end(), currentFrame) == frames.end()) {
                frames.push_back(currentFrame);
                std::sort(frames.begin(), frames.end());
                if (static_cast<int>(frames.size()) > maxPoints)
                    frames.erase(frames.begin());
            }
            candidateFrameCount = span;
        }
        frames.erase(std::unique(frames.begin(), frames.end()), frames.end());
        if (frames.empty())
            return errorResponse(QStringLiteral("no frames selected"), SC::BadRequest);

        const int referenceFrame = std::clamp(
            body.contains(QStringLiteral("reference_frame"))
                ? body.value(QStringLiteral("reference_frame")).toInt(currentFrame)
                : currentFrame,
            0, frameCount - 1);
        QString coordinateSpace = body.value(QStringLiteral("coordinate_space"))
                                      .toString(ring ? QStringLiteral("source_ring_local")
                                                     : QStringLiteral("world"))
                                      .toLower();
        if (coordinateSpace == QStringLiteral("frozen_ring"))
            coordinateSpace = QStringLiteral("source_ring_local");
        if (coordinateSpace != QStringLiteral("world") &&
            coordinateSpace != QStringLiteral("source_ring_local")) {
            return errorResponse(
                QStringLiteral("coordinate_space must be world or source_ring_local"),
                SC::BadRequest);
        }
        if (coordinateSpace == QStringLiteral("source_ring_local") && !ring) {
            return errorResponse(
                QStringLiteral("coordinate_space=source_ring_local requires {\"ring\": N}"),
                SC::BadRequest);
        }
        RingLocalFrame referenceRingFrame;
        if (coordinateSpace == QStringLiteral("source_ring_local")) {
            referenceRingFrame =
                ringLocalFrameAt(*conf, *ring, static_cast<std::size_t>(referenceFrame));
            if (!referenceRingFrame.valid)
                return errorResponse(QStringLiteral("reference ring geometry invalid"),
                                     SC::BadRequest);
        }

        auto styledDouble = [&](const QString& key, double fallback, double lo, double hi) {
            const QJsonValue v = body.value(key);
            if (!v.isDouble()) return fallback;
            return std::clamp(v.toDouble(), lo, hi);
        };
        auto styledBool = [&](const QString& key, bool fallback) {
            const QJsonValue v = body.value(key);
            return v.isBool() ? v.toBool() : fallback;
        };

        AtomTrackOverlay::Style style;
        style.pointSizePixels = styledDouble(QStringLiteral("point_size"),
                                             style.pointSizePixels, 1.0, 72.0);
        style.sphereRadiusA = styledDouble(QStringLiteral("dot_radius_A"),
                                           style.sphereRadiusA, 0.002, 0.20);
        style.currentPointScale = styledDouble(QStringLiteral("current_point_scale"),
                                               style.currentPointScale, 0.5, 6.0);
        style.pointOpacity = styledDouble(QStringLiteral("point_opacity"),
                                          style.pointOpacity, 0.0, 1.0);
        style.haloScale = styledDouble(QStringLiteral("halo_scale"),
                                       style.haloScale, 1.0, 10.0);
        style.haloOpacity = styledDouble(QStringLiteral("halo_opacity"),
                                         style.haloOpacity, 0.0, 1.0);
        style.lineWidthPixels = styledDouble(QStringLiteral("line_width"),
                                             style.lineWidthPixels, 1.0, 12.0);
        style.lineOpacity = styledDouble(QStringLiteral("line_opacity"), style.lineOpacity,
                                         0.0, 1.0);
        style.colorScale = styledDouble(QStringLiteral("color_scale"),
                                        style.colorScale, 0.0, 1e6);
        style.colorGamma = styledDouble(QStringLiteral("color_gamma"),
                                        style.colorGamma, 0.1, 4.0);
        style.minColorFraction = styledDouble(QStringLiteral("min_color_fraction"),
                                              style.minColorFraction, 0.0, 0.6);
        style.showPoints = styledBool(QStringLiteral("show_points"), style.showPoints);
        style.showHalos = styledBool(QStringLiteral("show_halos"), style.showHalos);
        style.showLines = styledBool(QStringLiteral("show_lines"), style.showLines);
        style.highlightCurrent = styledBool(QStringLiteral("highlight_current"),
                                            style.highlightCurrent);
        const QString pointShape = body.value(QStringLiteral("point_shape"))
                                       .toString(QStringLiteral("screen_point"))
                                       .toLower();
        if (pointShape == QStringLiteral("sphere") ||
            pointShape == QStringLiteral("spheres") ||
            pointShape == QStringLiteral("dot") ||
            pointShape == QStringLiteral("dots")) {
            style.pointShape = AtomTrackOverlay::PointShape::Sphere;
        } else if (pointShape == QStringLiteral("screen_point") ||
                   pointShape == QStringLiteral("screen_points") ||
                   pointShape == QStringLiteral("point") ||
                   pointShape == QStringLiteral("points")) {
            style.pointShape = AtomTrackOverlay::PointShape::ScreenPoint;
        } else {
            return errorResponse(
                QStringLiteral("point_shape must be screen_point or sphere"),
                SC::BadRequest);
        }
        const QString colorMode = body.value(QStringLiteral("color_mode"))
                                      .toString(QStringLiteral("signed")).toLower();
        if (colorMode == QStringLiteral("absolute")) {
            style.colorMode = AtomTrackOverlay::ColorMode::Absolute;
        } else if (colorMode == QStringLiteral("signed")) {
            style.colorMode = AtomTrackOverlay::ColorMode::Signed;
        } else {
            return errorResponse(QStringLiteral("color_mode must be signed or absolute"),
                                 SC::BadRequest);
        }
        std::vector<AtomTrackOverlay::Sample> samples;
        samples.reserve(frames.size());
        QJsonArray sampleJson;
        double minIntensity = std::numeric_limits<double>::infinity();
        double maxIntensity = -std::numeric_limits<double>::infinity();
        int ringValidCount = 0;
        int ringInvalidCount = 0;
        int vetComparedCount = 0;
        double vetMaxAbsKernelDelta = 0.0;
        double vetMaxAbsNullMarginDelta = 0.0;
        double vetMaxAbsDistanceDelta = 0.0;
        double vetMaxAbsLocalRoundtripDelta = 0.0;
        bool vetOk = true;
        for (int f : frames) {
            const model::Vec3 sourcePosition =
                conf->atomPosition(static_cast<std::size_t>(f), atom);
            model::Vec3 drawnPosition = sourcePosition;
            model::Vec3 sourceRingLocal = model::Vec3::Zero();
            model::RingNullMeasurement sourceMeasure;
            model::RingNullMeasurement projectedMeasure;
            double sourceKernel = std::numeric_limits<double>::quiet_NaN();
            double projectedKernel = std::numeric_limits<double>::quiet_NaN();
            bool ringMeasurementValid = false;

            AtomTrackOverlay::Sample sample;
            sample.current = (f == currentFrame);

            QJsonObject row{
                {"frame", f},
                {"original_frame", static_cast<qint64>(
                    conf->originalFrameIndex(static_cast<std::size_t>(f)))},
                {"current", sample.current},
            };
            if (ring) {
                if (coordinateSpace == QStringLiteral("source_ring_local")) {
                    const RingLocalFrame sourceRingFrame =
                        ringLocalFrameAt(*conf, *ring, static_cast<std::size_t>(f));
                    if (sourceRingFrame.valid) {
                        sourceRingLocal = toRingLocal(sourceRingFrame, sourcePosition);
                        drawnPosition = fromRingLocal(referenceRingFrame, sourceRingLocal);
                        sourceMeasure =
                            model::MeasureRingNull(sourceRingFrame.geometry, sourcePosition);
                        projectedMeasure =
                            model::MeasureRingNull(referenceRingFrame.geometry, drawnPosition);
                    }
                } else {
                    sourceMeasure =
                        model::MeasureRingNull(*conf, atom, *ring,
                                               static_cast<std::size_t>(f));
                    projectedMeasure = sourceMeasure;
                }

                ringMeasurementValid = sourceMeasure.valid;
                if (!ringMeasurementValid) {
                    row.insert(QStringLiteral("ring_measurement_valid"), false);
                    sample.intensity = 0.0;
                    ++ringInvalidCount;
                    if (frameSource == QStringLiteral("dft") &&
                        coordinateSpace == QStringLiteral("source_ring_local")) {
                        continue;
                    }
                } else {
                    ++ringValidCount;
                    sourceKernel = ringCurrentExpectedValue(sourceMeasure);
                    projectedKernel = ringCurrentExpectedValue(projectedMeasure);
                    if (colorBy == QStringLiteral("ring_current")) {
                        sample.intensity = sourceKernel;
                    } else if (colorBy == QStringLiteral("null_margin")) {
                        sample.intensity = sourceMeasure.nullMarginA;
                    } else {
                        sample.intensity = 0.0;
                    }
                    row.insert(QStringLiteral("ring_measurement_valid"), true);
                    row.insert(QStringLiteral("ring_current_kernel"), finiteJson(sourceKernel));
                    row.insert(QStringLiteral("expected_relationship_value"),
                               finiteJson(sourceKernel));
                    row.insert(QStringLiteral("null_margin_A"),
                               finiteJson(sourceMeasure.nullMarginA));
                    row.insert(QStringLiteral("distance_A"), finiteJson(sourceMeasure.distanceA));
                    row.insert(QStringLiteral("angle_degrees"), finiteJson(sourceMeasure.angleDeg));
                    row.insert(QStringLiteral("side"),
                               QString::fromLatin1(model::RingNullSideName(sourceMeasure.side)));
                    if (coordinateSpace == QStringLiteral("source_ring_local")) {
                        row.insert(QStringLiteral("projected_ring_current_kernel"),
                                   finiteJson(projectedKernel));
                        row.insert(QStringLiteral("projected_null_margin_A"),
                                   finiteJson(projectedMeasure.nullMarginA));
                        row.insert(QStringLiteral("projected_distance_A"),
                                   finiteJson(projectedMeasure.distanceA));
                        row.insert(QStringLiteral("projected_angle_degrees"),
                                   finiteJson(projectedMeasure.angleDeg));
                        row.insert(QStringLiteral("source_ring_local"),
                                   vec3FromEigen(sourceRingLocal));

                        const model::Vec3 roundtripLocal =
                            toRingLocal(referenceRingFrame, drawnPosition);
                        const double localRoundtripDelta =
                            (roundtripLocal - sourceRingLocal).norm();
                        const double kernelDelta =
                            std::abs(projectedKernel - sourceKernel);
                        const double nullMarginDelta =
                            std::abs(projectedMeasure.nullMarginA -
                                     sourceMeasure.nullMarginA);
                        const double distanceDelta =
                            std::abs(projectedMeasure.distanceA -
                                     sourceMeasure.distanceA);
                        vetMaxAbsKernelDelta =
                            std::max(vetMaxAbsKernelDelta, kernelDelta);
                        vetMaxAbsNullMarginDelta =
                            std::max(vetMaxAbsNullMarginDelta, nullMarginDelta);
                        vetMaxAbsDistanceDelta =
                            std::max(vetMaxAbsDistanceDelta, distanceDelta);
                        vetMaxAbsLocalRoundtripDelta =
                            std::max(vetMaxAbsLocalRoundtripDelta,
                                     localRoundtripDelta);
                        ++vetComparedCount;
                        if (!projectedMeasure.valid ||
                            kernelDelta > 1e-9 ||
                            nullMarginDelta > 1e-9 ||
                            distanceDelta > 1e-9 ||
                            localRoundtripDelta > 1e-9) {
                            vetOk = false;
                        }
                        row.insert(QStringLiteral("vet_kernel_delta"),
                                   finiteJson(projectedKernel - sourceKernel));
                        row.insert(QStringLiteral("vet_null_margin_delta_A"),
                                   finiteJson(projectedMeasure.nullMarginA -
                                              sourceMeasure.nullMarginA));
                        row.insert(QStringLiteral("vet_distance_delta_A"),
                                   finiteJson(projectedMeasure.distanceA -
                                              sourceMeasure.distanceA));
                        row.insert(QStringLiteral("vet_local_roundtrip_delta_A"),
                                   finiteJson(localRoundtripDelta));
                    }
                }
            } else {
                sample.intensity = 0.0;
            }
            sample.position = drawnPosition;
            row.insert(QStringLiteral("source_position"), vec3FromEigen(sourcePosition));
            row.insert(QStringLiteral("drawn_position"), vec3FromEigen(drawnPosition));
            row.insert(QStringLiteral("position"), vec3FromEigen(drawnPosition));
            row.insert(QStringLiteral("intensity"), finiteJson(sample.intensity));
            if (std::isfinite(sample.intensity)) {
                minIntensity = std::min(minIntensity, sample.intensity);
                maxIntensity = std::max(maxIntensity, sample.intensity);
            }
            samples.push_back(sample);
            sampleJson.append(row);
        }
        if (samples.empty())
            return errorResponse(QStringLiteral("no valid samples selected"), SC::BadRequest);

        heroshotAtomTrack_ = std::make_unique<AtomTrackOverlay>(
            vtkSmartPointer<vtkRenderer>(scene_->Renderer()));
        heroshotAtomTrack_->show(samples, style);
        scene_->requestRender(MoleculeScene::RenderSource::Rest);

        QJsonObject out{
            {"atom", static_cast<qint64>(atom)},
            {"current_frame", currentFrame},
            {"reference_frame", referenceFrame},
            {"frame_source", frameSource},
            {"coordinate_space", coordinateSpace},
            {"color_by", colorBy},
            {"color_mode", colorMode},
            {"kept", static_cast<qint64>(samples.size())},
            {"max_points", maxPoints},
            {"candidate_frame_count", candidateFrameCount},
            {"decimated_by_max_points", decimatedByMaxPoints},
            {"show_lines", style.showLines},
            {"show_points", style.showPoints},
            {"show_halos", style.showHalos},
            {"point_shape", pointShape},
            {"dot_radius_A", finiteJson(style.sphereRadiusA)},
            {"color_scale", finiteJson(style.colorScale)},
            {"color_gamma", finiteJson(style.colorGamma)},
            {"min_color_fraction", finiteJson(style.minColorFraction)},
            {"ring_measurements_valid", ringValidCount},
            {"ring_measurements_invalid", ringInvalidCount},
            {"min_intensity", finiteJson(minIntensity)},
            {"max_intensity", finiteJson(maxIntensity)},
            {"vet", QJsonObject{
                {"ok", coordinateSpace == QStringLiteral("source_ring_local") ? vetOk : true},
                {"compared", vetComparedCount},
                {"max_abs_kernel_delta", finiteJson(vetMaxAbsKernelDelta)},
                {"max_abs_null_margin_delta_A", finiteJson(vetMaxAbsNullMarginDelta)},
                {"max_abs_distance_delta_A", finiteJson(vetMaxAbsDistanceDelta)},
                {"max_abs_local_roundtrip_delta_A",
                 finiteJson(vetMaxAbsLocalRoundtripDelta)},
            }},
            {"samples", sampleJson},
        };
        out.insert(QStringLiteral("ring"),
                   ring ? QJsonValue(static_cast<qint64>(*ring)) : QJsonValue());
        return jsonResponse(out);
    };
    server_->route(QStringLiteral("/resthero/atom_track"), Method::Post,
                   restheroAtomTrackHandler);

    // ---- resthero: tensor ghost trail -----------------------------------
    //
    // POST /resthero/ghost_trail  {atom?, n?, end_frame?, step?, frames?}
    //   Draw the focused (or {"atom":N}) atom's shielding tensor at its last
    //   `n` DFT-measured frames as a fading stack of glyphs -- newest opaque,
    //   oldest faint -- so the re-orientation reads "from the side". Walks
    //   backward from `end_frame` (default the live frame) by `step` (default
    //   2, the DFT cadence), keeping only frames whose tensor is valid; up to
    //   `n` (default 5, clamped 2..12) ghosts. Or pass `frames:[...]` to draw
    //   an explicit statistics-selected set such as rotamer wells. Each ghost
    //   is the REAL measured tensor at a REAL frame -- probeAtomCsa reads the
    //   DFT store directly, no interpolation and the live frame never moves.
    //   Resthero layer only: the reader's own single live glyph is untouched.
    //   Echoes the frames + fades actually drawn so a script can verify what it
    //   sees.
    //   POST /resthero/clear  removes the trail.
    auto restheroGhostTrailHandler = [this](const QHttpServerRequest& req) {
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
            for (const auto& v : std::as_const(explicitFrames)) {
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
    };
    server_->route(QStringLiteral("/resthero/ghost_trail"), Method::Post,
                   restheroGhostTrailHandler);

    // POST /resthero/angle_collar
    //   Draw a transient, chart-derived dihedral collar in the molecule scene.
    //   Defaults are aimed at the Asp29 chi2 baton-pass: CA-CB-CG-OD1, with
    //   the signed angle read from the loaded H5 dihedral_time_series rather
    //   than guessed from labels. The collar sits around the CG->CB axis so the
    //   signed sweep matches residue.chi2's chart/stat convention.
    auto restheroAngleCollarHandler = [this](const QHttpServerRequest& req) {
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
    };
    server_->route(QStringLiteral("/resthero/angle_collar"), Method::Post,
                   restheroAngleCollarHandler);

    auto restheroClearHandler = [this]() {
        ASSERT_THREAD(this);
        if (heroshotAtomTrack_) heroshotAtomTrack_->clear();
        if (heroshotButterfly_) heroshotButterfly_->clear();
        if (heroshotTrail_) heroshotTrail_->clear();
        if (heroshotAngleCollar_) heroshotAngleCollar_->clear();
        if (scene_ && scene_->measurementOverlay()
            && heroshotMeasurementVisibleBefore_.has_value()) {
            scene_->measurementOverlay()->setVisible(*heroshotMeasurementVisibleBefore_);
            heroshotMeasurementVisibleBefore_.reset();
        }
        if (heroshotFieldRingWasSet_) {
            if (scene_ && scene_->fieldGridOverlay())
                scene_->fieldGridOverlay()->setVisibleRing(heroshotFieldRingBefore_);
            heroshotFieldRingBefore_.reset();
            heroshotFieldRingWasSet_ = false;
        }
        if (scene_ && heroshotMoleculeStyleBefore_.has_value()) {
            scene_->applyMoleculeStyle(*heroshotMoleculeStyleBefore_);
            heroshotMoleculeStyleBefore_.reset();
        }
        if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Rest);
        return QHttpServerResponse(SC::NoContent);
    };
    server_->route(QStringLiteral("/resthero/clear"), Method::Post, restheroClearHandler);

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
            const QStringList displayModes = model::AllDisplayModes(d);
            for (const QString& m : std::as_const(displayModes)) modeArr.append(m);
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
            const QStringList displayModes = model::AllDisplayModes(d);
            for (const QString& m : std::as_const(displayModes)) modeArr.append(m);
            QJsonArray vizArr;
            const auto supportingViz = reg.supporting(d);
            for (const model::VisualizationDefinition* def : std::as_const(supportingViz)) {
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
            const QStringList shapeVisualizations = shapeViz.value(it.key());
            for (const QString& v : std::as_const(shapeVisualizations)) viz.append(v);
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
        const QStringList activeNames = diagnostics::StructuredLogger::SymbolicNamesFromMask(mask);
        for (const QString& n : std::as_const(activeNames))
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
            const QJsonArray categories = body.value("categories").toArray();
            for (const auto& v : std::as_const(categories))
                names.append(v.toString());
            mask = diagnostics::StructuredLogger::MaskFromSymbolicNames(names);
        } else {
            return errorResponse(
                QStringLiteral("body must be {\"mask\": int} or {\"categories\": [...]}"),
                SC::BadRequest);
        }
        diagnostics::StructuredLogger::SetCategoryMask(mask);
        QJsonArray cats;
        const QStringList activeNames = diagnostics::StructuredLogger::SymbolicNamesFromMask(mask);
        for (const QString& n : std::as_const(activeNames))
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
        for (const auto& v : std::as_const(atomsArr)) {
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
        const auto stripTracks = dashboardController_->stripTracks();
        for (const DashboardDisplayController::StripTrack& t : std::as_const(stripTracks))
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

        for (const QString& mode : std::as_const(modes)) {
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

    auto screenshotHandler = [this](const QHttpServerRequest& req, bool diagnostic) {
        ASSERT_THREAD(this);
        bool ok = false;
        const QJsonObject body = parseJsonBody(req, &ok);
        const QString target = (ok && body.contains("target"))
                                   ? body.value("target").toString()
                                   : QStringLiteral("scene");
        if (!diagnostic && ok && body.contains(QStringLiteral("force_render"))) {
            return errorResponse(
                QStringLiteral("force_render is diagnostics-only; use /diagnostics/screenshot"),
                SC::BadRequest);
        }

        bool forceRender = true;
        if (diagnostic && ok && body.contains(QStringLiteral("force_render"))) {
            if (!body.value(QStringLiteral("force_render")).isBool())
                return errorResponse(QStringLiteral("force_render must be bool"),
                                     SC::BadRequest);
            // false skips vtkWindowToImageFilter::ShouldRerender so diagnostic
            // probes can read the current framebuffer instead of forcing a new
            // scene render. Only meaningful for target="scene".
            forceRender = body.value(QStringLiteral("force_render")).toBool();
        }
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
            if (!diagnostic) {
                return errorResponse(QStringLiteral("target must be \"scene\" or \"window\""),
                                     SC::BadRequest);
            }
            const QString objectName = body.value(QStringLiteral("object_name")).toString();
            if (objectName.isEmpty())
                return errorResponse(QStringLiteral("target \"widget\" requires \"object_name\""), SC::BadRequest);
            QWidget* widget = findWidgetByObjectName(mainWindow_.data(), objectName);
            if (!widget)
                return errorResponse(QStringLiteral("no live widget named \"%1\"").arg(objectName), SC::NotFound);
            png = captureWindowPng(widget);
        } else {
            return errorResponse(
                diagnostic
                    ? QStringLiteral("target must be \"scene\", \"window\", or \"widget\"")
                    : QStringLiteral("target must be \"scene\" or \"window\""),
                SC::BadRequest);
        }
        if (png.isEmpty())
            return errorResponse(QStringLiteral("screenshot capture failed"), SC::InternalServerError);
        return QHttpServerResponse(QByteArray(kMimePng), png);
    };
    server_->route(QStringLiteral("/api/screenshot"), Method::Post,
                   [screenshotHandler](const QHttpServerRequest& req) {
        return screenshotHandler(req, false);
    });
    server_->route(QStringLiteral("/diagnostics/screenshot"), Method::Post,
                   [screenshotHandler](const QHttpServerRequest& req) {
        return screenshotHandler(req, true);
    });
}

}  // namespace h5reader::app
