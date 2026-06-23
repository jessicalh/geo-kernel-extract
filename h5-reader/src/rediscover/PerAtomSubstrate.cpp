#include "PerAtomSubstrate.h"

#include "CaseHunter.h"
#include "Catalog.h"
#include "LiteratureConstants.h"
#include "LocalFrameBasis.h"
#include "McConnellLiteratureKernel.h"
#include "RelationshipEngine.h"
#include "RingCurrentKernel.h"
#include "SphericalBasis.h"
#include "Verbs.h"

#include "../io/QtTrajectoryH5.h"
#include "../io/QtFieldCatalog.gen.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtConformationSnapshot.h"
#include "../model/QtPerResidueBuffers.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtRing.h"
#include "../model/QtTimeSeriesBuffers.h"
#include "../model/QtTopology.h"

#include <QByteArray>
#include <QDir>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLoggingCategory>
#include <QSaveFile>
#include <QTextStream>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cPerAtom, "h5reader.rediscover.per_atom_substrate")

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr float kNaNF32 = std::numeric_limits<float>::quiet_NaN();

QString num(double v) { return QString::number(v, 'g', 12); }

template <typename E>
int ord(E e) {
    return static_cast<int>(e);
}

bool finiteVec3(const Vec3& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

bool finiteT2(const std::array<double, 5>& t2) {
    for (double v : t2)
        if (!std::isfinite(v)) return false;
    return true;
}

bool finiteRaw(const double* values, std::size_t n) {
    if (!values) return false;
    for (std::size_t i = 0; i < n; ++i)
        if (!std::isfinite(values[i])) return false;
    return true;
}

double t2Magnitude(const std::array<double, 5>& t2) {
    if (!finiteT2(t2)) return kNaN;
    double s = 0.0;
    for (double v : t2) s += v * v;
    return std::sqrt(s);
}

double vecMagnitude(const Vec3& v, bool present) {
    return present && finiteVec3(v) ? v.norm() : kNaN;
}

double finiteOrZero(double v) {
    return std::isfinite(v) ? v : 0.0;
}

std::array<double, 5> nanT2() {
    std::array<double, 5> out;
    out.fill(kNaN);
    return out;
}

std::array<double, 3> nanT1() {
    std::array<double, 3> out;
    out.fill(kNaN);
    return out;
}

void addT2(std::array<double, 5>& dst, const std::array<double, 5>& src) {
    for (std::size_t i = 0; i < dst.size(); ++i) dst[i] += src[i];
}

void fillNaN(std::array<double, kPerAtomClassicalCols>& v) { v.fill(kNaN); }
void fillNaN(std::array<double, kPerAtomConditioningCols>& v) { v.fill(kNaN); }
void fillNaN(std::array<double, kPerAtomDominanceCols>& v) { v.fill(kNaN); }
void fillNaN(std::array<double, kPerAtomBackboneAuditCols>& v) { v.fill(kNaN); }
void fillNaN(std::array<double, kPerAtomRingPathCols>& v) { v.fill(kNaN); }
void fillNaN(std::array<double, kPerAtomMethodPathCols>& v) { v.fill(kNaN); }
void fillNaN(std::array<double, kPerAtomHbondConditioningCols>& v) { v.fill(kNaN); }

double dominanceFraction(const std::vector<double>& vals) {
    double sum = 0.0;
    double maxv = 0.0;
    for (double v : vals) {
        if (!std::isfinite(v)) continue;
        const double a = std::abs(v);
        sum += a;
        maxv = std::max(maxv, a);
    }
    return sum > 0.0 ? maxv / sum : kNaN;
}

const model::NpyColumn* indexedColumn(const model::QtConformationSnapshot* snapshot,
                                      io::FieldKind kind,
                                      std::size_t rowIndex,
                                      int minCols) {
    if (!snapshot || !snapshot->has(kind)) return nullptr;
    const model::NpyColumn& col = snapshot->column(kind);
    if (rowIndex >= static_cast<std::size_t>(std::max(0, col.rows))) return nullptr;
    if (col.cols < minCols) return nullptr;
    return &col;
}

const model::NpyColumn* atomColumn(const model::QtConformationSnapshot* snapshot,
                                   io::FieldKind kind,
                                   std::size_t atom,
                                   int minCols) {
    return indexedColumn(snapshot, kind, atom, minCols);
}

template <std::size_t N>
bool copyAtomField(std::array<double, N>& out,
                   std::size_t& offset,
                   const model::QtConformationSnapshot* snapshot,
                   io::FieldKind kind,
                   std::size_t atom,
                   int count) {
    const model::NpyColumn* col = atomColumn(snapshot, kind, atom, count);
    bool present = false;
    if (col) {
        const double* row = col->row(atom);
        present = finiteRaw(row, static_cast<std::size_t>(count));
        if (present) {
            for (int i = 0; i < count; ++i) out[offset + static_cast<std::size_t>(i)] = row[i];
        }
    }
    offset += static_cast<std::size_t>(count);
    return present;
}

template <std::size_t N>
bool copyIndexedField(std::array<double, N>& out,
                      std::size_t& offset,
                      const model::QtConformationSnapshot* snapshot,
                      io::FieldKind kind,
                      std::size_t rowIndex,
                      int count) {
    const model::NpyColumn* col = indexedColumn(snapshot, kind, rowIndex, count);
    bool present = false;
    if (col) {
        const double* row = col->row(rowIndex);
        present = finiteRaw(row, static_cast<std::size_t>(count));
        if (present) {
            for (int i = 0; i < count; ++i) out[offset + static_cast<std::size_t>(i)] = row[i];
        }
    }
    offset += static_cast<std::size_t>(count);
    return present;
}

template <std::size_t N>
bool copyAtomTensorT2(std::array<double, N>& out,
                      std::size_t& offset,
                      const model::QtConformationSnapshot* snapshot,
                      io::FieldKind kind,
                      std::size_t atom) {
    const model::NpyColumn* col = atomColumn(snapshot, kind, atom, 9);
    bool present = false;
    if (col) {
        const double* row = col->row(atom);
        present = finiteRaw(row + 4, 5);
        if (present) {
            for (int i = 0; i < 5; ++i) out[offset + static_cast<std::size_t>(i)] = row[4 + i];
        }
    }
    offset += 5;
    return present;
}

struct ChargeScalar {
    double value = kNaN;
    bool present = false;
};

struct RowChargeScalars {
    ChargeScalar ff14sb;
    ChargeScalar mopac_welford_mean;
    ChargeScalar eeq_charge;
    ChargeScalar eeq_coordination_number;
    bool charge_complete = false;
};

ChargeScalar catalogChargeScalar(const Body& body, ArrayId id, std::size_t atom, std::size_t row) {
    ChargeScalar out;
    if (!body.catalog.present(body, id, atom, row)) return out;
    const double value = body.catalog.value(body, id, atom, row);
    if (!std::isfinite(value)) return out;
    out.value = value;
    out.present = true;
    return out;
}

QString csvScalar(const ChargeScalar& scalar) {
    return scalar.present ? num(scalar.value) : QStringLiteral("NaN");
}

bool validResidue(const model::QtProtein& p, int32_t r) {
    return r >= 0 && static_cast<std::size_t>(r) < p.residueCount();
}

std::vector<std::size_t> allAtomsStratum(const Body& body) {
    std::vector<std::size_t> atoms;
    if (!body.run.protein) return atoms;
    atoms.resize(body.run.protein->atomCount());
    std::iota(atoms.begin(), atoms.end(), std::size_t{0});
    return atoms;
}

FrameResult labFrameFn(const Body&, std::size_t, std::size_t) {
    return {};  // default axes, invalid local-frame flag; primary payload is lab frame
}

enum class RoleOrd : int {
    BackboneN = 0,
    BackboneCA = 1,
    BackboneC = 2,
    BackboneO = 3,
    BackboneHN = 4,
    BackboneHA = 5,
    AromaticH = 6,
    PolarSidechainH = 7,
    AliphaticH = 8,
    AromaticHeavy = 9,
    SidechainHeavy = 10,
    TerminalOrCap = 11,
    Other = 12,
};

enum class StratumOrd : int {
    N = 0,
    CA = 1,
    C = 2,
    O = 3,
    HN = 4,
    HA = 5,
    AromaticH = 6,
    PolarH = 7,
    AliphaticH = 8,
    AromaticHeavy = 9,
    ChargedHeavy = 10,
    CarbonylSidechain = 11,
    AmideSidechain = 12,
    Sulfur = 13,
    OtherHeavy = 14,
    Other = 15,
};

QString roleName(RoleOrd r) {
    switch (r) {
    case RoleOrd::BackboneN: return QStringLiteral("backbone_N");
    case RoleOrd::BackboneCA: return QStringLiteral("backbone_CA");
    case RoleOrd::BackboneC: return QStringLiteral("backbone_C");
    case RoleOrd::BackboneO: return QStringLiteral("backbone_O");
    case RoleOrd::BackboneHN: return QStringLiteral("backbone_HN");
    case RoleOrd::BackboneHA: return QStringLiteral("backbone_HA");
    case RoleOrd::AromaticH: return QStringLiteral("aromatic_H");
    case RoleOrd::PolarSidechainH: return QStringLiteral("polar_sidechain_H");
    case RoleOrd::AliphaticH: return QStringLiteral("aliphatic_H");
    case RoleOrd::AromaticHeavy: return QStringLiteral("aromatic_heavy");
    case RoleOrd::SidechainHeavy: return QStringLiteral("sidechain_heavy");
    case RoleOrd::TerminalOrCap: return QStringLiteral("terminal_or_cap");
    case RoleOrd::Other: return QStringLiteral("other");
    }
    return QStringLiteral("other");
}

QString stratumName(StratumOrd s) {
    switch (s) {
    case StratumOrd::N: return QStringLiteral("N");
    case StratumOrd::CA: return QStringLiteral("CA");
    case StratumOrd::C: return QStringLiteral("C");
    case StratumOrd::O: return QStringLiteral("O");
    case StratumOrd::HN: return QStringLiteral("HN");
    case StratumOrd::HA: return QStringLiteral("HA");
    case StratumOrd::AromaticH: return QStringLiteral("aromatic_H");
    case StratumOrd::PolarH: return QStringLiteral("polar_H");
    case StratumOrd::AliphaticH: return QStringLiteral("aliphatic_H");
    case StratumOrd::AromaticHeavy: return QStringLiteral("aromatic_heavy");
    case StratumOrd::ChargedHeavy: return QStringLiteral("charged_heavy");
    case StratumOrd::CarbonylSidechain: return QStringLiteral("carbonyl_sidechain");
    case StratumOrd::AmideSidechain: return QStringLiteral("amide_sidechain");
    case StratumOrd::Sulfur: return QStringLiteral("sulfur");
    case StratumOrd::OtherHeavy: return QStringLiteral("other_heavy");
    case StratumOrd::Other: return QStringLiteral("other");
    }
    return QStringLiteral("other");
}

RoleOrd roleForAtom(const model::QtAtom& a) {
    if (a.IsBackboneNitrogen()) return RoleOrd::BackboneN;
    if (a.IsBackboneAlphaCarbon()) return RoleOrd::BackboneCA;
    if (a.IsBackboneCarbonylCarbon()) return RoleOrd::BackboneC;
    if (a.IsBackboneCarbonylOxygen()) return RoleOrd::BackboneO;
    if (a.IsBackboneAmideHydrogen()) return RoleOrd::BackboneHN;
    if (a.IsAnyAlphaHydrogen()) return RoleOrd::BackboneHA;
    if (a.residueIndex < 0) return RoleOrd::TerminalOrCap;
    if (a.element == model::Element::H) {
        if (a.IsAromaticRingHydrogen()) return RoleOrd::AromaticH;
        if (a.IsPolarH() || a.isExchangeable) return RoleOrd::PolarSidechainH;
        return RoleOrd::AliphaticH;
    }
    if (a.aromatic) return RoleOrd::AromaticHeavy;
    if (a.element != model::Element::Unknown) return RoleOrd::SidechainHeavy;
    return RoleOrd::Other;
}

StratumOrd stratumForAtom(const model::QtAtom& a) {
    if (a.IsBackboneNitrogen()) return StratumOrd::N;
    if (a.IsBackboneAlphaCarbon()) return StratumOrd::CA;
    if (a.IsBackboneCarbonylCarbon()) return StratumOrd::C;
    if (a.IsBackboneCarbonylOxygen()) return StratumOrd::O;
    if (a.IsBackboneAmideHydrogen()) return StratumOrd::HN;
    if (a.IsAnyAlphaHydrogen()) return StratumOrd::HA;
    if (a.element == model::Element::H) {
        if (a.IsAromaticRingHydrogen()) return StratumOrd::AromaticH;
        if (a.IsPolarH() || a.isExchangeable) return StratumOrd::PolarH;
        return StratumOrd::AliphaticH;
    }
    if (a.aromatic) return StratumOrd::AromaticHeavy;
    if (a.formalCharge != 0) return StratumOrd::ChargedHeavy;
    if (a.planarGroup == model::PlanarGroupKind::Carboxylate) return StratumOrd::CarbonylSidechain;
    if (a.planarGroup == model::PlanarGroupKind::SidechainAmide) return StratumOrd::AmideSidechain;
    if (a.element == model::Element::S) return StratumOrd::Sulfur;
    if (a.element != model::Element::H && a.element != model::Element::Unknown)
        return StratumOrd::OtherHeavy;
    return StratumOrd::Other;
}

bool isBackboneBroadAtom(const model::QtAtom& a) {
    return a.IsBackbone() || a.IsAnyAlphaHydrogen();
}

FrameResult backboneAuditFrame(const Body& body, std::size_t atom, std::size_t frame) {
    FrameResult fr;
    const model::QtProtein& p = *body.run.protein;
    const model::QtAtom& a = p.atom(atom);
    if (!validResidue(p, a.residueIndex)) return fr;
    const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
    auto posOf = [&](int32_t ai) { return verbs::pos(body, static_cast<std::size_t>(ai), frame); };

    if (a.IsBackboneNitrogen()) {
        if (r.N == model::QtResidue::NONE || r.CA == model::QtResidue::NONE) return fr;
        bool cPrevValid = false;
        Vec3 cRef = Vec3::Zero();
        if (validResidue(p, r.prevResidueIndex)) {
            const model::QtResidue& prev = p.residue(static_cast<std::size_t>(r.prevResidueIndex));
            if (prev.C != model::QtResidue::NONE) {
                cRef = posOf(prev.C);
                cPrevValid = true;
            }
        }
        if (!cPrevValid && r.C != model::QtResidue::NONE) cRef = posOf(r.C);
        fr.frame = BuildBackboneNFrame(posOf(r.N), posOf(r.CA), cRef, cPrevValid);
        fr.anchor_atom_index = cPrevValid && validResidue(p, r.prevResidueIndex)
                                   ? p.residue(static_cast<std::size_t>(r.prevResidueIndex)).C
                                   : r.C;
        return fr;
    }
    if (a.IsBackboneAlphaCarbon()) {
        if (r.CA == model::QtResidue::NONE || r.N == model::QtResidue::NONE || r.C == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCaFrame(posOf(r.CA), posOf(r.N), posOf(r.C));
        fr.anchor_atom_index = r.N;
        return fr;
    }
    if (a.IsBackboneCarbonylCarbon()) {
        if (r.C == model::QtResidue::NONE || r.O == model::QtResidue::NONE || r.CA == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCarbonylCFrame(posOf(r.C), posOf(r.O), posOf(r.CA));
        fr.anchor_atom_index = r.CA;
        return fr;
    }
    if (a.IsBackboneCarbonylOxygen()) {
        if (r.O == model::QtResidue::NONE || r.C == model::QtResidue::NONE || r.CA == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCarbonylOFrame(posOf(r.O), posOf(r.C), posOf(r.CA));
        fr.anchor_atom_index = r.CA;
        return fr;
    }
    if (a.IsBackboneAmideHydrogen()) {
        if (r.N == model::QtResidue::NONE || r.CA == model::QtResidue::NONE) return fr;
        bool cPrevValid = false;
        Vec3 cPrev = Vec3::Zero();
        if (validResidue(p, r.prevResidueIndex)) {
            const model::QtResidue& prev = p.residue(static_cast<std::size_t>(r.prevResidueIndex));
            if (prev.C != model::QtResidue::NONE) {
                cPrev = posOf(prev.C);
                cPrevValid = true;
            }
        }
        fr.frame = BuildHNFrame(posOf(r.N), verbs::pos(body, atom, frame), posOf(r.CA), cPrev, cPrevValid);
        return fr;
    }
    if (a.IsAnyAlphaHydrogen()) {
        if (r.CA == model::QtResidue::NONE || r.N == model::QtResidue::NONE) return fr;
        fr.frame = BuildBackboneHaFrame(verbs::pos(body, atom, frame), posOf(r.CA), posOf(r.N));
        fr.anchor_atom_index = r.N;
        return fr;
    }
    return fr;
}

class StreamingNpy {
public:
    StreamingNpy() = default;

    bool open(const QString& path, const std::vector<std::size_t>& shape, const QByteArray& descr) {
        file_ = std::make_unique<QSaveFile>(path);
        if (!file_->open(QIODevice::WriteOnly)) return false;
        expected_ = 1;
        for (std::size_t dim : shape) expected_ *= dim;
        written_ = 0;

        QByteArray header;
        header += "{'descr': '";
        header += descr;
        header += "', 'fortran_order': False, 'shape': (";
        for (std::size_t i = 0; i < shape.size(); ++i) {
            if (i) header += ", ";
            header += QByteArray::number(static_cast<qulonglong>(shape[i]));
        }
        if (shape.size() == 1) header += ",";
        header += "), }";

        constexpr int kPreambleBytes = 10;
        const int newlineBytes = 1;
        int pad = 16 - ((kPreambleBytes + header.size() + newlineBytes) % 16);
        if (pad == 16) pad = 0;
        header += QByteArray(pad, ' ');
        header += '\n';
        if (header.size() > 65535) return false;

        QByteArray prefix;
        prefix.append("\x93NUMPY", 6);
        prefix.append(char(1));
        prefix.append(char(0));
        const quint16 headerLen = static_cast<quint16>(header.size());
        prefix.append(char(headerLen & 0xff));
        prefix.append(char((headerLen >> 8) & 0xff));
        if (file_->write(prefix) != prefix.size()) return false;
        return file_->write(header) == header.size();
    }

    template <typename T>
    bool writeScalar(T value) {
        if (!file_) return false;
        const qsizetype n = static_cast<qsizetype>(sizeof(T));
        if (file_->write(reinterpret_cast<const char*>(&value), n) != n) return false;
        ++written_;
        return true;
    }

    template <typename T, std::size_t N>
    bool writeArray(const std::array<T, N>& values) {
        for (T v : values)
            if (!writeScalar<T>(v)) return false;
        return true;
    }

    bool commit() {
        return file_ && written_ == expected_ && file_->commit();
    }

private:
    std::unique_ptr<QSaveFile> file_;
    std::size_t expected_ = 0;
    std::size_t written_ = 0;
};

template <typename T>
bool writeNpy(const QString& path, const std::vector<std::size_t>& shape,
              const std::vector<T>& data, const QByteArray& descr) {
    StreamingNpy writer;
    if (!writer.open(path, shape, descr)) return false;
    for (T v : data)
        if (!writer.writeScalar<T>(v)) return false;
    return writer.commit();
}

struct WelfordCell {
    std::size_t n = 0;
    double mean = 0.0;
    double m2 = 0.0;

    void push(double x) {
        if (!std::isfinite(x)) return;
        ++n;
        const double delta = x - mean;
        mean += delta / static_cast<double>(n);
        const double delta2 = x - mean;
        m2 += delta * delta2;
    }

    double sd() const {
        return n > 1 ? std::sqrt(m2 / static_cast<double>(n - 1)) : 0.0;
    }
};

struct MechanismAggregate {
    bool ring_present = false;
    double ring_T0 = 0.0;
    std::array<double, 5> ring_T2 = {};
    int ring_n = 0;
    int ring_valid_n = 0;
    int ring_self_or_bonded_n = 0;

    bool charge_present = false;
    std::array<double, 5> charge_T2 = {};
    int charge_n = 0;
    int charge_excluded_same_residue_n = 0;
    bool ff14sb_field_present = false;
    Vec3 ff14sb_field = Vec3::Zero();

    bool mc_lit_valid_present = false;
    model::SphericalTensor mc_lit_valid;
    int bond_n = 0;
    int bond_n_valid = 0;
    int bond_self_or_bonded_n = 0;

    PerAtomIsolationScalars isolation;
    std::array<double, kPerAtomBackboneAuditCols> backbone_audit = {};
};

bool ringSourceIsSelfOrBonded(const Body& body, std::size_t targetAtom, const model::QtRing& sourceRing) {
    const std::vector<int32_t> ownRings = verbs::ringsOf(body, targetAtom);
    for (int32_t r : ownRings)
        if (r == sourceRing.ringId) return true;
    const std::vector<int32_t> ownAtoms = verbs::ownRingAtoms(body, targetAtom);
    for (int32_t ra : sourceRing.atomIndices)
        for (int32_t oa : ownAtoms)
            if (oa == ra) return true;
    const std::size_t heavy = verbs::heavyParent(body, targetAtom);
    return std::find(sourceRing.atomIndices.begin(), sourceRing.atomIndices.end(),
                     static_cast<int32_t>(heavy)) != sourceRing.atomIndices.end();
}

PairContribution makeRingContribution(const Body& body, std::size_t targetAtom,
                                      std::size_t row, const SourceRef& ref,
                                      const LocalFrame& frame) {
    PairContribution out;
    out.kernel_T2.fill(kNaN);
    if (ref.entity_index < 0 || !body.run.protein) return out;
    const model::QtProtein& p = *body.run.protein;
    const std::size_t ringIdx = static_cast<std::size_t>(ref.entity_index);
    if (ringIdx >= p.topology().ringCount()) return out;
    const model::QtRing& ring = p.topology().ringAt(ringIdx);
    const model::RingGeometry& g = verbs::ringGeom(body, ringIdx, row);
    const Vec3 targetPos = verbs::pos(body, targetAtom, row);
    const verbs::Displacement d = verbs::displacement(targetPos, g.center, g.normal);
    const model::SphericalTensor unit =
        JohnsonBoveySourceUnitKernelLocal(body, frame, targetPos, ringIdx, ring, row);
    const model::SphericalTensor fixed = ScaleSphericalTensor(unit, ring.LiteratureIntensity());
    const bool self = ringSourceIsSelfOrBonded(body, targetAtom, ring);

    out.mechanism = QStringLiteral("ring_jb");
    out.source_kind = QStringLiteral("ring_center");
    out.source_id = static_cast<int32_t>(ringIdx);
    out.source_cloud_index = ref.cloud_index;
    out.source_category_ord = ring.TypeIndexAsInt();
    out.pointer_flags = PresentFlag | (self ? SelfOrBondedFlag : ProducerValidFlag);
    out.disp = d.disp;
    out.r = d.r;
    out.inv_r3 = d.inv_r3;
    out.cos_theta = d.cos_theta;
    out.dipolar = d.dipolar;
    out.kernel_T0 = fixed.T0;
    out.kernel_T2 = fixed.T2;
    out.contribution = RingForwardContributionPpm(fixed);
    return out;
}

bool chargeRejectedSameResidue(const model::QtProtein& p, std::size_t targetAtom, std::size_t sourceAtom) {
    const model::QtAtom& t = p.atom(targetAtom);
    const model::QtAtom& s = p.atom(sourceAtom);
    return t.residueIndex >= 0 && s.residueIndex == t.residueIndex;
}

bool chargeSourceAcceptedForContribution(const Body& body, const model::QtProtein& p,
                                         std::size_t targetAtom, std::size_t sourceAtom,
                                         std::size_t row, double distance) {
    if (sourceAtom >= p.atomCount() || sourceAtom == targetAtom) return false;
    if (!(distance > 1e-9)) return false;
    if (chargeRejectedSameResidue(p, targetAtom, sourceAtom)) return false;
    if (!body.catalog.has(ArrayId::Ff14sbCharge)) return false;
    if (!body.catalog.present(body, ArrayId::Ff14sbCharge, sourceAtom, row)) return false;
    return std::isfinite(body.catalog.value(body, ArrayId::Ff14sbCharge, sourceAtom, row));
}

PairContribution makeChargeContribution(const Body& body, std::size_t targetAtom,
                                        std::size_t row, const SourceRef& ref,
                                        bool excludeSameResidue) {
    PairContribution out;
    out.kernel_T2.fill(kNaN);
    if (!body.run.protein || ref.entity_index < 0 || !body.catalog.has(ArrayId::Ff14sbCharge))
        return out;
    const model::QtProtein& p = *body.run.protein;
    const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
    if (sourceAtom >= p.atomCount() || sourceAtom == targetAtom) return out;
    const bool sameResidue = chargeRejectedSameResidue(p, targetAtom, sourceAtom);
    if (excludeSameResidue && sameResidue) return out;
    if (!body.catalog.present(body, ArrayId::Ff14sbCharge, sourceAtom, row)) return out;
    const double q = body.catalog.value(body, ArrayId::Ff14sbCharge, sourceAtom, row);
    if (!std::isfinite(q)) return out;
    const Vec3 disp = verbs::pos(body, sourceAtom, row) - verbs::pos(body, targetAtom, row);
    const double r = disp.norm();
    if (!(r > 1e-9)) return out;
    const double r3 = r * r * r;
    const double r5 = r3 * r * r;
    Mat3 efg = q * (3.0 * disp * disp.transpose() / r5 - Mat3::Identity() / r3);
    efg *= CoulombKeVA();
    efg -= (efg.trace() / 3.0) * Mat3::Identity();
    const model::SphericalTensor st = DecomposeLibrary(efg);

    out.mechanism = QStringLiteral("charge_q_over_r3");
    out.source_kind = QStringLiteral("ff14sb_charge_site");
    out.source_id = static_cast<int32_t>(sourceAtom);
    out.source_atom_index = static_cast<int32_t>(sourceAtom);
    out.source_cloud_index = ref.cloud_index;
    out.source_category_ord = 0;
    out.pointer_flags = PresentFlag | ProducerValidFlag;
    out.disp = disp;
    out.r = r;
    out.inv_r3 = 1.0 / r3;
    out.kernel_T0 = st.T0;
    out.kernel_T2 = st.T2;
    out.contribution = t2Magnitude(st.T2);
    return out;
}

bool mopacFieldSourceAcceptedForContribution(const Body& body, const model::QtProtein& p,
                                             std::size_t targetAtom, std::size_t sourceAtom,
                                             std::size_t row, double distance) {
    if (sourceAtom >= p.atomCount() || sourceAtom == targetAtom) return false;
    if (!(distance > 1e-9)) return false;
    if (chargeRejectedSameResidue(p, targetAtom, sourceAtom)) return false;
    if (!body.catalog.has(ArrayId::MopacChargeWelfordMean)) return false;
    if (!body.catalog.present(body, ArrayId::MopacChargeWelfordMean, sourceAtom, row)) return false;
    return std::isfinite(body.catalog.value(body, ArrayId::MopacChargeWelfordMean, sourceAtom, row));
}

PairContribution makeMopacFieldContribution(const Body& body, std::size_t targetAtom,
                                            std::size_t row, const SourceRef& ref) {
    PairContribution out;
    out.kernel_T2.fill(kNaN);
    if (!body.run.protein || ref.entity_index < 0) return out;
    const model::QtProtein& p = *body.run.protein;
    const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
    if (sourceAtom >= p.atomCount() || sourceAtom == targetAtom) return out;
    const Vec3 disp = verbs::pos(body, sourceAtom, row) - verbs::pos(body, targetAtom, row);
    const double r = disp.norm();
    if (!mopacFieldSourceAcceptedForContribution(body, p, targetAtom, sourceAtom, row, r)) return out;
    const double q = body.catalog.value(body, ArrayId::MopacChargeWelfordMean, sourceAtom, row);
    const double r3 = r * r * r;
    const Vec3 field = (-q / r3) * disp;

    out.mechanism = QStringLiteral("field_mopac_coulomb");
    out.source_kind = QStringLiteral("mopac_welford_charge_site");
    out.source_id = static_cast<int32_t>(sourceAtom);
    out.source_atom_index = static_cast<int32_t>(sourceAtom);
    out.source_cloud_index = ref.cloud_index;
    out.source_category_ord = 0;
    out.pointer_flags = PresentFlag | ProducerValidFlag;
    out.disp = disp;
    out.r = r;
    out.inv_r3 = 1.0 / r3;
    out.contribution = field.norm();
    return out;
}

PairContribution makeBondContribution(const Body& body, std::size_t targetAtom,
                                      std::size_t row, const SourceRef& ref,
                                      double nearFieldRatio, const LocalFrame& frame) {
    PairContribution out;
    out.kernel_T2.fill(kNaN);
    if (!body.run.protein || ref.entity_index < 0) return out;
    const model::QtProtein& p = *body.run.protein;
    const std::size_t bondIdx = static_cast<std::size_t>(ref.entity_index);
    if (bondIdx >= p.topology().bondCount()) return out;
    const model::QtBond& b = p.topology().bondAt(bondIdx);
    if (b.atomIndexA < 0 || b.atomIndexB < 0) return out;
    const Vec3 posA = verbs::pos(body, static_cast<std::size_t>(b.atomIndexA), row);
    const Vec3 posB = verbs::pos(body, static_cast<std::size_t>(b.atomIndexB), row);
    const Vec3 axis = posB - posA;
    const double axisNorm = axis.norm();
    if (!(axisNorm > 1e-9)) return out;
    const Vec3 axisU = axis / axisNorm;
    const Vec3 targetPos = verbs::pos(body, targetAtom, row);
    const Vec3 midpoint = 0.5 * (posA + posB);
    const verbs::Displacement d = verbs::displacement(targetPos, midpoint, axis);
    if (!(d.r > 1e-6) || !std::isfinite(d.dipolar)) return out;
    const bool endpointSelf = static_cast<int32_t>(targetAtom) == b.atomIndexA
                              || static_cast<int32_t>(targetAtom) == b.atomIndexB;
    const bool nearField = d.r <= axisNorm * nearFieldRatio;

    SourceSlot slot;
    slot.kind = SourceKind::Bond;
    slot.disp_local = frame.is_valid ? frame.ToLocal(d.disp) : d.disp;
    slot.r = d.r;
    slot.cos_theta = d.cos_theta;
    slot.dipolar = d.dipolar;
    slot.bond_category = static_cast<int>(b.category);
    slot.bond_order = static_cast<int>(b.order);
    slot.bond_index = b.bondIndex;
    slot.bond_atom_a = b.atomIndexA;
    slot.bond_atom_b = b.atomIndexB;
    slot.bond_axis_local = frame.is_valid ? frame.ToLocal(axisU) : axisU;
    slot.mc_source_is_self_or_bonded = endpointSelf || nearField;
    bool litPresent = false;
    const model::SphericalTensor lit = McConnellSourceLiteratureKernelLocal(slot, &litPresent);

    out.mechanism = QStringLiteral("mc_lit_valid");
    out.source_kind = QStringLiteral("bond_midpoint");
    out.source_id = b.bondIndex;
    out.source_cloud_index = ref.cloud_index;
    out.source_category_ord = static_cast<int>(b.category);
    out.pointer_flags = PresentFlag
                        | (slot.mc_source_is_self_or_bonded ? SelfOrBondedFlag : ProducerValidFlag)
                        | (nearField ? NearFieldFlag : 0);
    out.disp = d.disp;
    out.r = d.r;
    out.inv_r3 = d.inv_r3;
    out.cos_theta = d.cos_theta;
    out.dipolar = d.dipolar;
    if (litPresent && !slot.mc_source_is_self_or_bonded) {
        out.kernel_T0 = lit.T0;
        out.kernel_T2 = lit.T2;
        out.contribution = t2Magnitude(lit.T2);
    }
    return out;
}

}  // namespace

std::vector<PairContribution> PerAtomRowPairContributions(const Body& body, std::size_t atom,
                                                          std::size_t row,
                                                          const PerAtomSubstrateConfig& cfg,
                                                          const LocalFrame& frame) {
    std::vector<PairContribution> out;
    for (const SourceRef& ref : verbs::near(body, CloudKind::RingCenters, atom, row, cfg.ring_cutoff_A))
        out.push_back(makeRingContribution(body, atom, row, ref, frame));
    for (const SourceRef& ref : verbs::near(body, CloudKind::ChargeSites, atom, row, cfg.charge_cutoff_A)) {
        PairContribution c = makeChargeContribution(body, atom, row, ref, true);
        if (c.pointer_flags & PresentFlag) out.push_back(std::move(c));
        PairContribution f = makeMopacFieldContribution(body, atom, row, ref);
        if (f.pointer_flags & PresentFlag) out.push_back(std::move(f));
    }
    for (const SourceRef& ref : verbs::near(body, CloudKind::BondMidpoints, atom, row, cfg.bond_cutoff_A))
        out.push_back(makeBondContribution(body, atom, row, ref, cfg.mc_near_field_ratio, frame));
    return out;
}

namespace {

std::array<double, 5> chargeKernelT2Local(const Body& body, std::size_t atom, std::size_t row,
                                          const PerAtomSubstrateConfig& cfg,
                                          const LocalFrame& frame) {
    if (!body.run.protein || !frame.is_valid || !body.catalog.has(ArrayId::Ff14sbCharge))
        return nanT2();
    const model::QtProtein& p = *body.run.protein;
    const Vec3 targetPos = verbs::pos(body, atom, row);
    Mat3 efg = Mat3::Zero();
    bool any = false;
    for (const SourceRef& ref : verbs::near(body, CloudKind::ChargeSites, atom, row, cfg.charge_cutoff_A)) {
        if (ref.entity_index < 0) continue;
        const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
        if (sourceAtom >= p.atomCount() || sourceAtom == atom) continue;
        if (chargeRejectedSameResidue(p, atom, sourceAtom)) continue;
        if (!body.catalog.present(body, ArrayId::Ff14sbCharge, sourceAtom, row)) continue;
        const double q = body.catalog.value(body, ArrayId::Ff14sbCharge, sourceAtom, row);
        if (!std::isfinite(q)) continue;
        const Vec3 disp = frame.ToLocal(verbs::pos(body, sourceAtom, row) - targetPos);
        const double r = disp.norm();
        if (!(r > 1e-9)) continue;
        const double r3 = r * r * r;
        const double r5 = r3 * r * r;
        efg += q * (3.0 * disp * disp.transpose() / r5 - Mat3::Identity() / r3);
        any = true;
    }
    if (!any) return std::array<double, 5>{};
    efg *= CoulombKeVA();
    efg -= (efg.trace() / 3.0) * Mat3::Identity();
    const model::SphericalTensor st = DecomposeLibrary(efg);
    return finiteT2(st.T2) ? st.T2 : nanT2();
}

double gapToSecondDistance(const Body& body, CloudKind kind, std::size_t atom, std::size_t row,
                           double cutoff, bool heavyOnly = false);
double gapToSecondMopacFieldDistance(const Body& body, std::size_t atom, std::size_t row,
                                     double cutoff);

MechanismAggregate reduceMechanisms(const Body& body, std::size_t atom, std::size_t row,
                                    const PerAtomSubstrateConfig& cfg) {
    MechanismAggregate agg;
    fillNaN(agg.backbone_audit);
    const LocalFrame lab;
    const std::vector<PairContribution> labPairs = PerAtomRowPairContributions(body, atom, row, cfg, lab);
    std::vector<double> ringVals;
    std::vector<double> chargeVals;
    std::vector<double> mcVals;
    std::vector<double> fieldVals;
    for (const PairContribution& p : labPairs) {
        if (p.mechanism == QStringLiteral("ring_jb") && finiteT2(p.kernel_T2)) {
            agg.ring_present = true;
            agg.ring_T0 += p.kernel_T0;
            addT2(agg.ring_T2, p.kernel_T2);
            ++agg.ring_n;
            if (p.pointer_flags & SelfOrBondedFlag)
                ++agg.ring_self_or_bonded_n;
            else
                ++agg.ring_valid_n;
            ringVals.push_back(p.contribution);
        } else if (p.mechanism == QStringLiteral("charge_q_over_r3") && finiteT2(p.kernel_T2)) {
            agg.charge_present = true;
            addT2(agg.charge_T2, p.kernel_T2);
            ++agg.charge_n;
            chargeVals.push_back(p.contribution);
            if (p.r > 1e-9 && std::isfinite(p.inv_r3)) {
                // FIELD E = sum q_i (r_atom - r_i) / r^3. The per-source q is
                // reconstructed from q/r^3 only for the vector field by reading
                // the catalog again below, preserving the same source filter.
            }
        } else if (p.mechanism == QStringLiteral("mc_lit_valid") && finiteT2(p.kernel_T2)) {
            agg.mc_lit_valid_present = true;
            agg.mc_lit_valid.T0 += p.kernel_T0;
            for (std::size_t i = 0; i < 5; ++i) agg.mc_lit_valid.T2[i] += p.kernel_T2[i];
            mcVals.push_back(p.contribution);
        } else if (p.mechanism == QStringLiteral("field_mopac_coulomb")
                   && std::isfinite(p.contribution)) {
            fieldVals.push_back(p.contribution);
        }
        if (p.mechanism == QStringLiteral("mc_lit_valid")) {
            ++agg.bond_n;
            if (p.pointer_flags & SelfOrBondedFlag)
                ++agg.bond_self_or_bonded_n;
            else
                ++agg.bond_n_valid;
        }
    }
    agg.isolation.dominant_fraction_ring = dominanceFraction(ringVals);
    agg.isolation.dominant_fraction_charge = dominanceFraction(chargeVals);
    agg.isolation.dominant_fraction_mc = dominanceFraction(mcVals);
    agg.isolation.dominant_fraction_field = dominanceFraction(fieldVals);
    agg.isolation.gap_to_2nd_ring_r =
        gapToSecondDistance(body, CloudKind::RingCenters, atom, row, 1000.0);
    agg.isolation.gap_to_2nd_charge_r =
        gapToSecondDistance(body, CloudKind::ChargeSites, atom, row, 1000.0);
    agg.isolation.gap_to_2nd_bond_r =
        gapToSecondDistance(body, CloudKind::BondMidpoints, atom, row, 1000.0);
    agg.isolation.gap_to_2nd_field_r =
        gapToSecondMopacFieldDistance(body, atom, row, cfg.charge_cutoff_A);

    // The FF14SB field uses the same included charge source set as the q/r^3
    // T2 reducer. Count excluded same-residue charge sites separately.
    if (body.run.protein && body.catalog.has(ArrayId::Ff14sbCharge)) {
        const model::QtProtein& p = *body.run.protein;
        const Vec3 targetPos = verbs::pos(body, atom, row);
        for (const SourceRef& ref : verbs::near(body, CloudKind::ChargeSites, atom, row, cfg.charge_cutoff_A)) {
            if (ref.entity_index < 0) continue;
            const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
            if (sourceAtom >= p.atomCount() || sourceAtom == atom) continue;
            if (chargeRejectedSameResidue(p, atom, sourceAtom)) {
                ++agg.charge_excluded_same_residue_n;
                continue;
            }
            if (!body.catalog.present(body, ArrayId::Ff14sbCharge, sourceAtom, row)) continue;
            const double q = body.catalog.value(body, ArrayId::Ff14sbCharge, sourceAtom, row);
            const Vec3 disp = verbs::pos(body, sourceAtom, row) - targetPos;
            const double r = disp.norm();
            if (!(r > 1e-9) || !std::isfinite(q)) continue;
            agg.ff14sb_field += (-q / (r * r * r)) * disp;
            agg.ff14sb_field_present = true;
        }
    }

    // Backbone-local audit sidecar: same reducer choices and filters, but
    // expressed in the broad-backbone local frame so validation can compare
    // current all-atom backbone rows against the existing broad substrate.
    if (body.run.protein && isBackboneBroadAtom(body.run.protein->atom(atom))) {
        const FrameResult fr = backboneAuditFrame(body, atom, row);
        if (fr.frame.is_valid) {
            const std::vector<PairContribution> localPairs =
                PerAtomRowPairContributions(body, atom, row, cfg, fr.frame);
            const std::array<double, 5> chargeLocal = chargeKernelT2Local(body, atom, row, cfg, fr.frame);
            Vec3 fieldLocal = Vec3::Zero();
            model::SphericalTensor mcLocal;
            bool mcLocalPresent = false;
            if (body.run.protein && body.catalog.has(ArrayId::Ff14sbCharge)) {
                const model::QtProtein& p = *body.run.protein;
                const Vec3 targetPos = verbs::pos(body, atom, row);
                for (const SourceRef& ref : verbs::near(body, CloudKind::ChargeSites, atom, row, cfg.charge_cutoff_A)) {
                    if (ref.entity_index < 0) continue;
                    const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
                    if (sourceAtom >= p.atomCount() || sourceAtom == atom) continue;
                    if (chargeRejectedSameResidue(p, atom, sourceAtom)) continue;
                    if (!body.catalog.present(body, ArrayId::Ff14sbCharge, sourceAtom, row)) continue;
                    const double q = body.catalog.value(body, ArrayId::Ff14sbCharge, sourceAtom, row);
                    const Vec3 disp = verbs::pos(body, sourceAtom, row) - targetPos;
                    const double r = disp.norm();
                    if (!(r > 1e-9) || !std::isfinite(q)) continue;
                    fieldLocal += (-q / (r * r * r)) * fr.frame.ToLocal(disp);
                }
            }
            for (const PairContribution& pc : localPairs) {
                if (pc.mechanism != QStringLiteral("mc_lit_valid") || !finiteT2(pc.kernel_T2))
                    continue;
                mcLocalPresent = true;
                mcLocal.T0 += pc.kernel_T0;
                for (std::size_t i = 0; i < 5; ++i) mcLocal.T2[i] += pc.kernel_T2[i];
            }
            std::size_t c = 0;
            for (double v : chargeLocal) agg.backbone_audit[c++] = v;
            agg.backbone_audit[c++] = fieldLocal.x();
            agg.backbone_audit[c++] = fieldLocal.y();
            agg.backbone_audit[c++] = fieldLocal.z();
            agg.backbone_audit[c++] = mcLocalPresent ? mcLocal.T0 : kNaN;
            for (double v : mcLocal.T2) agg.backbone_audit[c++] = mcLocalPresent ? v : kNaN;
        }
    }

    if (!agg.ring_present) agg.ring_T2 = nanT2();
    if (!agg.charge_present) agg.charge_T2 = nanT2();
    return agg;
}

double nearestDistance(const Body& body, CloudKind kind, std::size_t atom, std::size_t row,
                       double cutoff, bool heavyOnly = false) {
    if (!body.run.protein) return kNaN;
    const model::QtProtein& p = *body.run.protein;
    const Vec3 query = verbs::pos(body, atom, row);
    double best = std::numeric_limits<double>::infinity();
    for (const SourceRef& ref : verbs::near(body, kind, atom, row, cutoff)) {
        if (ref.cloud_index < 0) continue;
        const Vec3 point = body.idx.spatial.tree(kind, row).pointAt(static_cast<std::size_t>(ref.cloud_index));
        const double d = (point - query).norm();
        if (!(d > 1e-9) || !std::isfinite(d)) continue;
        if (kind == CloudKind::Atoms) {
            if (ref.entity_index < 0) continue;
            const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
            if (sourceAtom == atom || sourceAtom >= p.atomCount()) continue;
            if (heavyOnly && p.atom(sourceAtom).element == model::Element::H) continue;
        } else if (kind == CloudKind::ChargeSites) {
            if (ref.entity_index < 0) continue;
            const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
            if (!chargeSourceAcceptedForContribution(body, p, atom, sourceAtom, row, d)) continue;
        }
        best = std::min(best, d);
    }
    return std::isfinite(best) ? best : kNaN;
}

double gapToSecondDistance(const Body& body, CloudKind kind, std::size_t atom, std::size_t row,
                           double cutoff, bool heavyOnly) {
    if (!body.run.protein) return kNaN;
    const model::QtProtein& p = *body.run.protein;
    const Vec3 query = verbs::pos(body, atom, row);
    std::array<double, 2> best = {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
    };
    for (const SourceRef& ref : verbs::near(body, kind, atom, row, cutoff)) {
        if (ref.cloud_index < 0) continue;
        const Vec3 point = body.idx.spatial.tree(kind, row).pointAt(static_cast<std::size_t>(ref.cloud_index));
        const double d = (point - query).norm();
        if (!(d > 1e-9) || !std::isfinite(d)) continue;
        if (kind == CloudKind::Atoms) {
            if (ref.entity_index < 0) continue;
            const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
            if (sourceAtom == atom || sourceAtom >= p.atomCount()) continue;
            if (heavyOnly && p.atom(sourceAtom).element == model::Element::H) continue;
        } else if (kind == CloudKind::ChargeSites) {
            if (ref.entity_index < 0) continue;
            const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
            if (!chargeSourceAcceptedForContribution(body, p, atom, sourceAtom, row, d)) continue;
        }
        if (d < best[0]) {
            best[1] = best[0];
            best[0] = d;
        } else if (d < best[1]) {
            best[1] = d;
        }
    }
    return std::isfinite(best[0]) && std::isfinite(best[1]) ? best[1] - best[0] : kNaN;
}

double gapToSecondMopacFieldDistance(const Body& body, std::size_t atom, std::size_t row,
                                     double cutoff) {
    if (!body.run.protein) return kNaN;
    const model::QtProtein& p = *body.run.protein;
    const Vec3 query = verbs::pos(body, atom, row);
    std::array<double, 2> best = {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
    };
    for (const SourceRef& ref : verbs::near(body, CloudKind::ChargeSites, atom, row, cutoff)) {
        if (ref.entity_index < 0 || ref.cloud_index < 0) continue;
        const std::size_t sourceAtom = static_cast<std::size_t>(ref.entity_index);
        const Vec3 point = body.idx.spatial.tree(CloudKind::ChargeSites, row)
                               .pointAt(static_cast<std::size_t>(ref.cloud_index));
        const double d = (point - query).norm();
        if (!mopacFieldSourceAcceptedForContribution(body, p, atom, sourceAtom, row, d)) continue;
        if (d < best[0]) {
            best[1] = best[0];
            best[0] = d;
        } else if (d < best[1]) {
            best[1] = d;
        }
    }
    return std::isfinite(best[0]) && std::isfinite(best[1]) ? best[1] - best[0] : kNaN;
}

}  // namespace

PerAtomIsolationScalars PerAtomIsolationScalarsForRow(const Body& body,
                                                      std::size_t atom,
                                                      std::size_t row,
                                                      const PerAtomSubstrateConfig& cfg) {
    PerAtomIsolationScalars out;
    out.gap_to_2nd_ring_r = gapToSecondDistance(body, CloudKind::RingCenters, atom, row, 1000.0);
    out.gap_to_2nd_charge_r = gapToSecondDistance(body, CloudKind::ChargeSites, atom, row, 1000.0);
    out.gap_to_2nd_bond_r = gapToSecondDistance(body, CloudKind::BondMidpoints, atom, row, 1000.0);
    out.gap_to_2nd_field_r = gapToSecondMopacFieldDistance(body, atom, row, cfg.charge_cutoff_A);
    std::vector<double> ringVals;
    std::vector<double> chargeVals;
    std::vector<double> mcVals;
    std::vector<double> fieldVals;
    for (const PairContribution& p : PerAtomRowPairContributions(body, atom, row, cfg, LocalFrame{})) {
        if (p.mechanism == QStringLiteral("ring_jb")) ringVals.push_back(p.contribution);
        if (p.mechanism == QStringLiteral("charge_q_over_r3")) chargeVals.push_back(p.contribution);
        if (p.mechanism == QStringLiteral("mc_lit_valid")) mcVals.push_back(p.contribution);
        if (p.mechanism == QStringLiteral("field_mopac_coulomb")) fieldVals.push_back(p.contribution);
    }
    out.dominant_fraction_ring = dominanceFraction(ringVals);
    out.dominant_fraction_charge = dominanceFraction(chargeVals);
    out.dominant_fraction_mc = dominanceFraction(mcVals);
    out.dominant_fraction_field = dominanceFraction(fieldVals);
    return out;
}

namespace {

int countNear(const Body& body, CloudKind kind, std::size_t atom, std::size_t row, double cutoff) {
    return static_cast<int>(verbs::near(body, kind, atom, row, cutoff).size());
}

bool isPeptidePlaneAtom(const model::QtAtom& a) {
    return a.IsBackboneNitrogen() || a.IsBackboneAmideHydrogen()
           || a.IsBackboneCarbonylOxygen() || a.IsBackboneCarbonylCarbon();
}

template <std::size_t N>
bool copyDsspRawBackup(std::array<double, N>& out,
                       std::size_t& offset,
                       const Body& body,
                       std::size_t atom,
                       std::size_t h5Row) {
    bool present = false;
    std::array<double, 10> values;
    values.fill(kNaN);
    if (body.run.protein) {
        const model::QtProtein& p = *body.run.protein;
        const model::QtAtom& a = p.atom(atom);
        const io::QtTrajectoryH5* h5 = body.run.h5();
        const model::QtDssp8TimeSeries* dssp = h5 ? h5->dssp8() : nullptr;
        if (dssp && validResidue(p, a.residueIndex) && isPeptidePlaneAtom(a)
            && h5Row < dssp->n_frames && dssp->sourceAttachedAt(h5Row)) {
            const std::size_t res = static_cast<std::size_t>(a.residueIndex);
            const std::size_t base = (res * dssp->n_frames + h5Row) * 2;
            if (base + 1 < dssp->hbond_acceptor_partner.size()
                && base + 1 < dssp->hbond_acceptor_energy.size()
                && base + 1 < dssp->hbond_donor_partner.size()
                && base + 1 < dssp->hbond_donor_energy.size()) {
                int accCount = 0;
                int donCount = 0;
                std::size_t c = 0;
                for (int slot = 0; slot < 2; ++slot) {
                    const std::size_t i = base + static_cast<std::size_t>(slot);
                    const int32_t partner = dssp->hbond_acceptor_partner[i];
                    const double energy = dssp->hbond_acceptor_energy[i];
                    if (partner >= 0 && std::isfinite(energy)) ++accCount;
                    values[c++] = static_cast<double>(partner);
                    values[c++] = std::isfinite(energy) ? energy : 0.0;
                }
                values[c++] = static_cast<double>(accCount);
                for (int slot = 0; slot < 2; ++slot) {
                    const std::size_t i = base + static_cast<std::size_t>(slot);
                    const int32_t partner = dssp->hbond_donor_partner[i];
                    const double energy = dssp->hbond_donor_energy[i];
                    if (partner >= 0 && std::isfinite(energy)) ++donCount;
                    values[c++] = static_cast<double>(partner);
                    values[c++] = std::isfinite(energy) ? energy : 0.0;
                }
                values[c++] = static_cast<double>(donCount);
                present = finiteRaw(values.data(), values.size());
            }
        }
    }
    if (present) {
        for (std::size_t i = 0; i < values.size(); ++i) out[offset + i] = values[i];
    }
    offset += values.size();
    return present;
}

template <std::size_t N>
bool copyPrimaryRingGeometry(std::array<double, N>& out,
                             std::size_t& offset,
                             const Body& body,
                             const model::QtConformationSnapshot* snapshot,
                             std::size_t atom) {
    bool present = false;
    if (body.run.protein && snapshot && snapshot->has(io::FieldKind::RingGeometry)) {
        const model::QtTopology& topo = body.run.protein->topology();
        const model::NpyColumn& col = snapshot->column(io::FieldKind::RingGeometry);
        for (int32_t mi : topo.ringMembershipsForAtom(atom)) {
            if (mi < 0) continue;
            const model::QtRingMembership& membership =
                topo.ringMembershipAt(static_cast<std::size_t>(mi));
            if (membership.ringId < 0
                || static_cast<std::size_t>(membership.ringId) >= topo.ringCount())
                continue;
            const model::QtRing& ring = topo.ringAt(static_cast<std::size_t>(membership.ringId));
            if (!ring.IsAromatic()) continue;
            const int32_t native = ring.nativeAxisIndex >= 0 ? ring.nativeAxisIndex : ring.ringId;
            if (native < 0 || static_cast<std::size_t>(native) >= static_cast<std::size_t>(std::max(0, col.rows))
                || col.cols < 10)
                continue;
            const double* row = col.row(static_cast<std::size_t>(native));
            present = finiteRaw(row, 10);
            if (present) {
                for (int i = 0; i < 10; ++i) out[offset + static_cast<std::size_t>(i)] = row[i];
            }
            break;
        }
    }
    offset += 10;
    return present;
}

struct DirectFeatures {
    DftTarget target;
    bool dft_present = false;
    bool apbs_E_present = false;
    Vec3 apbs_E = Vec3::Zero();
    bool apbs_efg_present = false;
    std::array<double, 5> apbs_efg_T2 = {};
    bool aimnet2_charge_present = false;
    double aimnet2_charge = kNaN;
    bool aimnet2_crg_present = false;
    double aimnet2_crg_scalar = kNaN;
    Vec3 aimnet2_crg = Vec3::Zero();
    bool aimnet2_embedding_present = false;
    const float* aimnet2_embedding = nullptr;
    std::size_t aimnet2_embedding_dims = 0;
    bool mopac_coulomb_present = false;
    std::array<double, 5> mopac_coulomb_T2 = {};
    bool mopac_mc_present = false;
    std::array<double, 5> mopac_mc_T2 = {};
    bool hbond_nearest_direction_present = false;
    Vec3 hbond_nearest_direction = Vec3::Zero();
    bool hbond_flags_present = false;
    std::array<double, 3> hbond_flags = {};
    bool hbond_count_present = false;
    double hbond_count = kNaN;
    bool hm_shielding_present = false;
    std::array<double, 5> hm_T2 = {};
    bool water_efield_present = false;
    Vec3 water_efield = Vec3::Zero();
    bool water_n_first_present = false;
    double water_n_first = kNaN;
    bool water_n_second_present = false;
    double water_n_second = kNaN;
    bool hydration_half_shell_present = false;
    double hydration_half_shell_asymmetry = kNaN;
    bool hydration_dipole_cos_present = false;
    double hydration_dipole_cos = kNaN;
    bool sasa_present = false;
    double sasa = kNaN;
    bool sasa_normal_present = false;
    Vec3 sasa_normal = Vec3::Zero();
    std::array<double, kPerAtomRingPathCols> ring_paths = {};
    bool bs_shielding_path_present = false;
    bool bs_per_type_T0_present = false;
    bool bs_per_type_T2_present = false;
    bool bs_total_B_present = false;
    bool bs_ring_counts_present = false;
    bool hm_shielding_path_present = false;
    bool hm_per_type_T0_present = false;
    bool hm_per_type_T2_present = false;
    bool pq_per_type_T0_present = false;
    bool disp_per_type_T0_present = false;
    std::array<double, kPerAtomMethodPathCols> method_paths = {};
    bool mc_scalars_present = false;
    bool mopac_mc_scalars_present = false;
    bool mopac_bond_orders_aggregate_present = false;
    bool mopac_coulomb_E_present = false;
    bool mopac_coulomb_efg_backbone_present = false;
    bool mopac_coulomb_efg_aromatic_present = false;
    bool mopac_coulomb_scalars_present = false;
    bool aimnet2_efg_present = false;
    bool aimnet2_efg_aromatic_present = false;
    bool aimnet2_efg_backbone_present = false;
    bool water_efg_present = false;
    bool water_efield_first_present = false;
    bool eeq_cn_present = false;
    bool mopac_scalars_present = false;
    std::array<double, kPerAtomHbondConditioningCols> hbond_conditioning = {};
    bool larsen_1pHB_T2_present = false;
    bool larsen_2pHB_T2_present = false;
    bool larsen_1pHaB_T2_present = false;
    bool larsen_2pHaB_T2_present = false;
    bool larsen_water_term_present = false;
    bool hbond_scalars_present = false;
    bool dssp_chemical_flags_present = false;
    bool dssp_hbond_energy_present = false;
    bool dssp_raw_backup_present = false;
    bool dssp_ss8_present = false;
    bool dssp_chi_present = false;
    bool omega_actual_present = false;
    bool pyramidalization_present = false;
    bool ring_geometry_present = false;
};

DirectFeatures directFeatures(const Body& body, std::size_t atom, std::size_t row,
                              std::size_t orig, const LocalFrame& frame,
                              const model::QtConformationSnapshot* snapshot = nullptr,
                              const std::array<double, 4>* mopacBondAggregate = nullptr) {
    DirectFeatures out;
    fillNaN(out.ring_paths);
    fillNaN(out.method_paths);
    fillNaN(out.hbond_conditioning);
    out.target = BuildTarget(body.run, atom, orig, frame);
    out.dft_present = out.target.present && finiteT2(out.target.total_decomp.T2);
    out.apbs_E_present = body.catalog.present(body, ArrayId::ApbsEfield, atom, row);
    if (out.apbs_E_present) out.apbs_E = body.catalog.valueVec3(body, ArrayId::ApbsEfield, atom, row);
    out.apbs_E_present = out.apbs_E_present && finiteVec3(out.apbs_E);
    out.apbs_efg_present = body.catalog.present(body, ArrayId::ApbsEfg, atom, row);
    if (out.apbs_efg_present) out.apbs_efg_T2 = body.catalog.valueT2(body, ArrayId::ApbsEfg, atom, row);
    out.apbs_efg_present = out.apbs_efg_present && finiteT2(out.apbs_efg_T2);
    out.aimnet2_charge_present = body.catalog.present(body, ArrayId::Aimnet2Charge, atom, row);
    if (out.aimnet2_charge_present) out.aimnet2_charge = body.catalog.value(body, ArrayId::Aimnet2Charge, atom, row);
    out.aimnet2_charge_present = out.aimnet2_charge_present && std::isfinite(out.aimnet2_charge);
    out.aimnet2_crg_present =
        body.catalog.present(body, ArrayId::Aimnet2ChargeRespScalar, atom, row)
        && body.catalog.present(body, ArrayId::Aimnet2ChargeRespVector, atom, row);
    if (out.aimnet2_crg_present) {
        out.aimnet2_crg_scalar = body.catalog.value(body, ArrayId::Aimnet2ChargeRespScalar, atom, row);
        out.aimnet2_crg = body.catalog.valueVec3(body, ArrayId::Aimnet2ChargeRespVector, atom, row);
    }
    out.aimnet2_crg_present = out.aimnet2_crg_present && std::isfinite(out.aimnet2_crg_scalar)
                              && finiteVec3(out.aimnet2_crg);
    out.aimnet2_embedding =
        body.catalog.valueEmbedding(body, ArrayId::Aimnet2Embedding, atom, row, out.aimnet2_embedding_dims);
    out.aimnet2_embedding_present =
        out.aimnet2_embedding && body.catalog.present(body, ArrayId::Aimnet2Embedding, atom, row);
    out.mopac_coulomb_present = body.catalog.present(body, ArrayId::MopacCoulombShielding, atom, row);
    if (out.mopac_coulomb_present)
        out.mopac_coulomb_T2 = body.catalog.valueT2(body, ArrayId::MopacCoulombShielding, atom, row);
    out.mopac_coulomb_present = out.mopac_coulomb_present && finiteT2(out.mopac_coulomb_T2);
    out.mopac_mc_present = body.catalog.present(body, ArrayId::MopacMcShielding, atom, row);
    if (out.mopac_mc_present)
        out.mopac_mc_T2 = body.catalog.valueT2(body, ArrayId::MopacMcShielding, atom, row);
    out.mopac_mc_present = out.mopac_mc_present && finiteT2(out.mopac_mc_T2);
    if (const model::NpyColumn* hdir =
            atomColumn(snapshot, io::FieldKind::HBondNearestDir, atom, 3)) {
        const double* v = hdir->row(atom);
        out.hbond_nearest_direction = Vec3(v[0], v[1], v[2]);
        out.hbond_nearest_direction_present = finiteRaw(v, 3);
    } else {
        out.hbond_nearest_direction_present =
            body.catalog.present(body, ArrayId::HbondNearestDirection, atom, row);
        if (out.hbond_nearest_direction_present)
            out.hbond_nearest_direction =
                body.catalog.valueVec3(body, ArrayId::HbondNearestDirection, atom, row);
        out.hbond_nearest_direction_present =
            out.hbond_nearest_direction_present && finiteVec3(out.hbond_nearest_direction);
    }
    if (const model::NpyColumn* hflags =
            atomColumn(snapshot, io::FieldKind::HBondFlags, atom, 3)) {
        const double* v = hflags->row(atom);
        for (std::size_t i = 0; i < out.hbond_flags.size(); ++i) out.hbond_flags[i] = v[i];
        out.hbond_flags_present = finiteRaw(v, 3);
    } else {
        out.hbond_flags_present = body.catalog.present(body, ArrayId::HbondFlags, atom, row);
        if (out.hbond_flags_present) {
            for (std::size_t i = 0; i < out.hbond_flags.size(); ++i)
                out.hbond_flags[i] =
                    body.catalog.value(body, ArrayId::HbondFlags, atom, row, -1, static_cast<int>(i));
        }
        out.hbond_flags_present =
            out.hbond_flags_present && finiteRaw(out.hbond_flags.data(), out.hbond_flags.size());
    }
    out.hbond_count_present = body.catalog.present(body, ArrayId::HbondCount, atom, row);
    if (out.hbond_count_present)
        out.hbond_count = body.catalog.value(body, ArrayId::HbondCount, atom, row);
    out.hbond_count_present = out.hbond_count_present && std::isfinite(out.hbond_count);
    out.hm_shielding_present = body.catalog.present(body, ArrayId::HmShielding, atom, row);
    if (out.hm_shielding_present)
        out.hm_T2 = body.catalog.valueT2(body, ArrayId::HmShielding, atom, row);
    out.hm_shielding_present = out.hm_shielding_present && finiteT2(out.hm_T2);
    out.water_efield_present = body.catalog.present(body, ArrayId::WaterEfield, atom, row);
    if (out.water_efield_present)
        out.water_efield = body.catalog.valueVec3(body, ArrayId::WaterEfield, atom, row);
    out.water_efield_present = out.water_efield_present && finiteVec3(out.water_efield);
    out.water_n_first_present = body.catalog.present(body, ArrayId::WaterNFirst, atom, row);
    if (out.water_n_first_present)
        out.water_n_first = body.catalog.value(body, ArrayId::WaterNFirst, atom, row);
    out.water_n_first_present = out.water_n_first_present && std::isfinite(out.water_n_first);
    out.water_n_second_present = body.catalog.present(body, ArrayId::WaterNSecond, atom, row);
    if (out.water_n_second_present)
        out.water_n_second = body.catalog.value(body, ArrayId::WaterNSecond, atom, row);
    out.water_n_second_present = out.water_n_second_present && std::isfinite(out.water_n_second);
    out.hydration_half_shell_present =
        body.catalog.present(body, ArrayId::HydrationHalfShellAsymmetry, atom, row);
    if (out.hydration_half_shell_present)
        out.hydration_half_shell_asymmetry =
            body.catalog.value(body, ArrayId::HydrationHalfShellAsymmetry, atom, row);
    out.hydration_half_shell_present =
        out.hydration_half_shell_present && std::isfinite(out.hydration_half_shell_asymmetry);
    out.hydration_dipole_cos_present =
        body.catalog.present(body, ArrayId::HydrationDipoleCos, atom, row);
    if (out.hydration_dipole_cos_present)
        out.hydration_dipole_cos = body.catalog.value(body, ArrayId::HydrationDipoleCos, atom, row);
    out.hydration_dipole_cos_present =
        out.hydration_dipole_cos_present && std::isfinite(out.hydration_dipole_cos);
    out.sasa_present = body.catalog.present(body, ArrayId::Sasa, atom, row);
    if (out.sasa_present) out.sasa = body.catalog.value(body, ArrayId::Sasa, atom, row);
    out.sasa_present = out.sasa_present && std::isfinite(out.sasa);
    out.sasa_normal_present = body.catalog.present(body, ArrayId::SasaNormal, atom, row);
    if (out.sasa_normal_present)
        out.sasa_normal = body.catalog.valueVec3(body, ArrayId::SasaNormal, atom, row);
    out.sasa_normal_present = out.sasa_normal_present && finiteVec3(out.sasa_normal);

    std::size_t rp = 0;
    out.bs_shielding_path_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::BSShielding, atom, 9);
    out.bs_per_type_T0_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::BSPerTypeT0, atom, 8);
    out.bs_per_type_T2_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::BSPerTypeT2, atom, 40);
    out.bs_total_B_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::BSTotalB, atom, 3);
    out.bs_ring_counts_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::BSRingCounts, atom, 4);
    out.hm_shielding_path_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::HMShielding, atom, 9);
    out.hm_per_type_T0_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::HMPerTypeT0, atom, 8);
    out.hm_per_type_T2_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::HMPerTypeT2, atom, 40);
    out.pq_per_type_T0_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::PQPerTypeT0, atom, 8);
    out.disp_per_type_T0_present =
        copyAtomField(out.ring_paths, rp, snapshot, io::FieldKind::DispPerTypeT0, atom, 8);

    std::size_t mp = 0;
    out.mc_scalars_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::McScalars, atom, 6);
    out.mopac_mc_scalars_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::MOPACMcScalars, atom, 6);
    if (mopacBondAggregate && finiteRaw(mopacBondAggregate->data(), mopacBondAggregate->size())) {
        for (std::size_t i = 0; i < mopacBondAggregate->size(); ++i)
            out.method_paths[mp + i] = (*mopacBondAggregate)[i];
        out.mopac_bond_orders_aggregate_present = true;
    }
    mp += 4;
    out.mopac_coulomb_E_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::MOPACCoulombE, atom, 3);
    out.mopac_coulomb_efg_backbone_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::MOPACCoulombEFGBackbone, atom, 5);
    out.mopac_coulomb_efg_aromatic_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::MOPACCoulombEFGAromatic, atom, 5);
    out.mopac_coulomb_scalars_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::MOPACCoulombScalars, atom, 4);
    out.aimnet2_efg_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::AIMNet2EFG, atom, 5);
    out.aimnet2_efg_aromatic_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::AIMNet2EFGAromatic, atom, 5);
    out.aimnet2_efg_backbone_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::AIMNet2EFGBackbone, atom, 5);
    out.water_efg_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::WaterEFG, atom, 5);
    out.water_efield_first_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::WaterEfieldFirst, atom, 3);
    out.eeq_cn_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::EEQCN, atom, 1);
    out.mopac_scalars_present =
        copyAtomField(out.method_paths, mp, snapshot, io::FieldKind::MOPACScalars, atom, 4);

    std::size_t hc = 0;
    out.larsen_1pHB_T2_present =
        copyAtomTensorT2(out.hbond_conditioning, hc, snapshot,
                         io::FieldKind::LarsenHBond1pHBShielding, atom);
    out.larsen_2pHB_T2_present =
        copyAtomTensorT2(out.hbond_conditioning, hc, snapshot,
                         io::FieldKind::LarsenHBond2pHBShielding, atom);
    out.larsen_1pHaB_T2_present =
        copyAtomTensorT2(out.hbond_conditioning, hc, snapshot,
                         io::FieldKind::LarsenHBond1pHaBShielding, atom);
    out.larsen_2pHaB_T2_present =
        copyAtomTensorT2(out.hbond_conditioning, hc, snapshot,
                         io::FieldKind::LarsenHBond2pHaBShielding, atom);
    out.larsen_water_term_present =
        copyAtomField(out.hbond_conditioning, hc, snapshot, io::FieldKind::LarsenHBondWaterTerm, atom, 1);
    out.hbond_scalars_present =
        copyAtomField(out.hbond_conditioning, hc, snapshot, io::FieldKind::HBondScalars, atom, 4);
    if (body.run.protein) {
        const model::QtAtom& a = body.run.protein->atom(atom);
        out.hbond_conditioning[hc++] =
            (a.IsBackboneNitrogen() || a.IsBackboneAmideHydrogen()) ? 1.0 : 0.0;
        out.hbond_conditioning[hc++] =
            (a.IsBackboneCarbonylOxygen() || a.IsBackboneCarbonylCarbon()) ? 1.0 : 0.0;
        out.dssp_chemical_flags_present = true;
    } else {
        hc += 2;
    }
    out.dssp_hbond_energy_present =
        copyAtomField(out.hbond_conditioning, hc, snapshot, io::FieldKind::DSSPHBondEnergy, atom, 4);
    out.dssp_raw_backup_present = copyDsspRawBackup(out.hbond_conditioning, hc, body, atom, row);
    out.dssp_ss8_present =
        copyAtomField(out.hbond_conditioning, hc, snapshot, io::FieldKind::DSSPSs8, atom, 8);
    out.dssp_chi_present =
        copyAtomField(out.hbond_conditioning, hc, snapshot, io::FieldKind::DSSPChi, atom, 12);
    if (body.run.protein && validResidue(*body.run.protein, body.run.protein->atom(atom).residueIndex)) {
        const std::size_t res = static_cast<std::size_t>(body.run.protein->atom(atom).residueIndex);
        out.omega_actual_present =
            copyIndexedField(out.hbond_conditioning, hc, snapshot, io::FieldKind::OmegaActual, res, 1);
    } else {
        hc += 1;
    }
    out.pyramidalization_present =
        copyAtomField(out.hbond_conditioning, hc, snapshot, io::FieldKind::Pyramidalization, atom, 1);
    out.ring_geometry_present =
        copyPrimaryRingGeometry(out.hbond_conditioning, hc, body, snapshot, atom);
    return out;
}

RowChargeScalars rowChargeScalars(const Body& body, std::size_t atom, std::size_t row,
                                  const DirectFeatures& direct) {
    RowChargeScalars out;
    out.ff14sb = catalogChargeScalar(body, ArrayId::Ff14sbCharge, atom, row);
    out.mopac_welford_mean =
        catalogChargeScalar(body, ArrayId::MopacChargeWelfordMean, atom, row);
    out.eeq_charge = catalogChargeScalar(body, ArrayId::EeqChargeMean, atom, row);
    out.eeq_coordination_number =
        catalogChargeScalar(body, ArrayId::EeqCoordinationNumber, atom, row);
    out.charge_complete = out.ff14sb.present && out.mopac_welford_mean.present
                          && direct.aimnet2_charge_present;
    return out;
}

void auditRange(PerAtomSubstrateStats& stats, const QString& name, bool present,
                const std::vector<double>& values) {
    PerAtomChannelAudit& audit = stats.new_channel_audit[name];
    if (!present) return;
    ++audit.present;
    for (double v : values) {
        if (!std::isfinite(v)) continue;
        if (!audit.has_range) {
            audit.min = v;
            audit.max = v;
            audit.has_range = true;
        } else {
            audit.min = std::min(audit.min, v);
            audit.max = std::max(audit.max, v);
        }
    }
}

void auditScalar(PerAtomSubstrateStats& stats, const QString& name, bool present, double value) {
    auditRange(stats, name, present, {value});
}

void auditVec3(PerAtomSubstrateStats& stats, const QString& name, bool present, const Vec3& value) {
    auditRange(stats, name, present, {value.x(), value.y(), value.z()});
}

void auditT2(PerAtomSubstrateStats& stats, const QString& name, bool present,
             const std::array<double, 5>& value) {
    auditRange(stats, name, present,
               {value[0], value[1], value[2], value[3], value[4]});
}

void auditT1(PerAtomSubstrateStats& stats, const QString& name, bool present,
             const std::array<double, 3>& value) {
    auditRange(stats, name, present, {value[0], value[1], value[2]});
}

template <std::size_t N>
void auditArraySegment(PerAtomSubstrateStats& stats,
                       const QString& name,
                       bool present,
                       const std::array<double, N>& values,
                       std::size_t offset,
                       std::size_t count) {
    std::vector<double> segment;
    segment.reserve(count);
    for (std::size_t i = 0; i < count; ++i) segment.push_back(values[offset + i]);
    auditRange(stats, name, present, segment);
}

void auditRingPathFeatures(PerAtomSubstrateStats& stats, const DirectFeatures& direct) {
    std::size_t c = 0;
    auditArraySegment(stats, QStringLiteral("bs_shielding"), direct.bs_shielding_path_present,
                      direct.ring_paths, c, 9);
    c += 9;
    auditArraySegment(stats, QStringLiteral("bs_per_type_T0"), direct.bs_per_type_T0_present,
                      direct.ring_paths, c, 8);
    c += 8;
    auditArraySegment(stats, QStringLiteral("bs_per_type_T2"), direct.bs_per_type_T2_present,
                      direct.ring_paths, c, 40);
    c += 40;
    auditArraySegment(stats, QStringLiteral("bs_total_B"), direct.bs_total_B_present,
                      direct.ring_paths, c, 3);
    c += 3;
    auditArraySegment(stats, QStringLiteral("bs_ring_counts"), direct.bs_ring_counts_present,
                      direct.ring_paths, c, 4);
    c += 4;
    auditArraySegment(stats, QStringLiteral("hm_shielding"), direct.hm_shielding_path_present,
                      direct.ring_paths, c, 9);
    c += 9;
    auditArraySegment(stats, QStringLiteral("hm_per_type_T0"), direct.hm_per_type_T0_present,
                      direct.ring_paths, c, 8);
    c += 8;
    auditArraySegment(stats, QStringLiteral("hm_per_type_T2"), direct.hm_per_type_T2_present,
                      direct.ring_paths, c, 40);
    c += 40;
    auditArraySegment(stats, QStringLiteral("pq_per_type_T0"), direct.pq_per_type_T0_present,
                      direct.ring_paths, c, 8);
    c += 8;
    auditArraySegment(stats, QStringLiteral("disp_per_type_T0"), direct.disp_per_type_T0_present,
                      direct.ring_paths, c, 8);
}

void auditMethodPathFeatures(PerAtomSubstrateStats& stats, const DirectFeatures& direct) {
    std::size_t c = 0;
    auditArraySegment(stats, QStringLiteral("mc_scalars"), direct.mc_scalars_present,
                      direct.method_paths, c, 6);
    c += 6;
    auditArraySegment(stats, QStringLiteral("mopac_mc_scalars"),
                      direct.mopac_mc_scalars_present, direct.method_paths, c, 6);
    c += 6;
    auditArraySegment(stats, QStringLiteral("mopac_bond_orders_aggregate"),
                      direct.mopac_bond_orders_aggregate_present, direct.method_paths, c, 4);
    c += 4;
    auditArraySegment(stats, QStringLiteral("mopac_coulomb_E"),
                      direct.mopac_coulomb_E_present, direct.method_paths, c, 3);
    c += 3;
    auditArraySegment(stats, QStringLiteral("mopac_coulomb_efg_backbone"),
                      direct.mopac_coulomb_efg_backbone_present, direct.method_paths, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("mopac_coulomb_efg_aromatic"),
                      direct.mopac_coulomb_efg_aromatic_present, direct.method_paths, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("mopac_coulomb_scalars"),
                      direct.mopac_coulomb_scalars_present, direct.method_paths, c, 4);
    c += 4;
    auditArraySegment(stats, QStringLiteral("aimnet2_efg"), direct.aimnet2_efg_present,
                      direct.method_paths, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("aimnet2_efg_aromatic"),
                      direct.aimnet2_efg_aromatic_present, direct.method_paths, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("aimnet2_efg_backbone"),
                      direct.aimnet2_efg_backbone_present, direct.method_paths, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("water_efg"), direct.water_efg_present,
                      direct.method_paths, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("water_efield_first"),
                      direct.water_efield_first_present, direct.method_paths, c, 3);
    c += 3;
    auditArraySegment(stats, QStringLiteral("eeq_cn"), direct.eeq_cn_present,
                      direct.method_paths, c, 1);
    c += 1;
    auditArraySegment(stats, QStringLiteral("mopac_scalars"), direct.mopac_scalars_present,
                      direct.method_paths, c, 4);
}

void auditHbondConditioningFeatures(PerAtomSubstrateStats& stats, const DirectFeatures& direct) {
    std::size_t c = 0;
    auditArraySegment(stats, QStringLiteral("larsen_hbond_1pHB_T2"),
                      direct.larsen_1pHB_T2_present, direct.hbond_conditioning, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("larsen_hbond_2pHB_T2"),
                      direct.larsen_2pHB_T2_present, direct.hbond_conditioning, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("larsen_hbond_1pHaB_T2"),
                      direct.larsen_1pHaB_T2_present, direct.hbond_conditioning, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("larsen_hbond_2pHaB_T2"),
                      direct.larsen_2pHaB_T2_present, direct.hbond_conditioning, c, 5);
    c += 5;
    auditArraySegment(stats, QStringLiteral("larsen_hbond_water_term"),
                      direct.larsen_water_term_present, direct.hbond_conditioning, c, 1);
    c += 1;
    auditArraySegment(stats, QStringLiteral("hbond_scalars"), direct.hbond_scalars_present,
                      direct.hbond_conditioning, c, 4);
    c += 4;
    auditArraySegment(stats, QStringLiteral("dssp_chemical_flags"),
                      direct.dssp_chemical_flags_present, direct.hbond_conditioning, c, 2);
    c += 2;
    auditArraySegment(stats, QStringLiteral("dssp_hbond_energy"),
                      direct.dssp_hbond_energy_present, direct.hbond_conditioning, c, 4);
    c += 4;
    auditArraySegment(stats, QStringLiteral("dssp_raw_hbond_backup"),
                      direct.dssp_raw_backup_present, direct.hbond_conditioning, c, 10);
    c += 10;
    auditArraySegment(stats, QStringLiteral("dssp_ss8"), direct.dssp_ss8_present,
                      direct.hbond_conditioning, c, 8);
    c += 8;
    auditArraySegment(stats, QStringLiteral("dssp_chi"), direct.dssp_chi_present,
                      direct.hbond_conditioning, c, 12);
    c += 12;
    auditArraySegment(stats, QStringLiteral("omega_actual"), direct.omega_actual_present,
                      direct.hbond_conditioning, c, 1);
    c += 1;
    auditArraySegment(stats, QStringLiteral("pyramidalization"),
                      direct.pyramidalization_present, direct.hbond_conditioning, c, 1);
    c += 1;
    auditArraySegment(stats, QStringLiteral("ring_geometry"), direct.ring_geometry_present,
                      direct.hbond_conditioning, c, 10);
}

void auditIsolationFeatures(PerAtomSubstrateStats& stats, const PerAtomIsolationScalars& iso) {
    auditScalar(stats, QStringLiteral("gap_to_2nd_ring_r"),
                std::isfinite(iso.gap_to_2nd_ring_r), iso.gap_to_2nd_ring_r);
    auditScalar(stats, QStringLiteral("gap_to_2nd_charge_r"),
                std::isfinite(iso.gap_to_2nd_charge_r), iso.gap_to_2nd_charge_r);
    auditScalar(stats, QStringLiteral("gap_to_2nd_bond_r"),
                std::isfinite(iso.gap_to_2nd_bond_r), iso.gap_to_2nd_bond_r);
    auditScalar(stats, QStringLiteral("gap_to_2nd_field_r"),
                std::isfinite(iso.gap_to_2nd_field_r), iso.gap_to_2nd_field_r);
    auditScalar(stats, QStringLiteral("dominant_fraction_ring"),
                std::isfinite(iso.dominant_fraction_ring), iso.dominant_fraction_ring);
    auditScalar(stats, QStringLiteral("dominant_fraction_charge"),
                std::isfinite(iso.dominant_fraction_charge), iso.dominant_fraction_charge);
    auditScalar(stats, QStringLiteral("dominant_fraction_mc"),
                std::isfinite(iso.dominant_fraction_mc), iso.dominant_fraction_mc);
    auditScalar(stats, QStringLiteral("dominant_fraction_field"),
                std::isfinite(iso.dominant_fraction_field), iso.dominant_fraction_field);
}

void auditDirectFeatures(PerAtomSubstrateStats& stats, const DirectFeatures& direct,
                         const RowChargeScalars& charges) {
    const bool hbondGeometryPresent = direct.hbond_nearest_direction_present
                                      && direct.hbond_flags_present;

    auditScalar(stats, QStringLiteral("hbond_count"), direct.hbond_count_present, direct.hbond_count);
    auditVec3(stats, QStringLiteral("hbond_nearest_dir"), hbondGeometryPresent,
              direct.hbond_nearest_direction);
    auditScalar(stats, QStringLiteral("hbond_flags_is_backbone"), hbondGeometryPresent,
                direct.hbond_flags[0]);
    auditScalar(stats, QStringLiteral("hbond_flags_is_donor"), hbondGeometryPresent,
                direct.hbond_flags[1]);
    auditScalar(stats, QStringLiteral("hbond_flags_is_acceptor"), hbondGeometryPresent,
                direct.hbond_flags[2]);
    auditT2(stats, QStringLiteral("hm_shielding_T2"), direct.hm_shielding_present,
            direct.hm_T2);
    auditVec3(stats, QStringLiteral("water_efield"), direct.water_efield_present,
              direct.water_efield);
    auditScalar(stats, QStringLiteral("water_n_first"), direct.water_n_first_present,
                direct.water_n_first);
    auditScalar(stats, QStringLiteral("water_n_second"), direct.water_n_second_present,
                direct.water_n_second);
    auditScalar(stats, QStringLiteral("water_half_shell_asymmetry"),
                direct.hydration_half_shell_present,
                direct.hydration_half_shell_asymmetry);
    auditScalar(stats, QStringLiteral("water_dipole_cos"), direct.hydration_dipole_cos_present,
                direct.hydration_dipole_cos);
    auditScalar(stats, QStringLiteral("sasa"), direct.sasa_present, direct.sasa);
    auditVec3(stats, QStringLiteral("sasa_normal"), direct.sasa_normal_present,
              direct.sasa_normal);
    auditScalar(stats, QStringLiteral("eeq_charge"), charges.eeq_charge.present,
                charges.eeq_charge.value);
    auditScalar(stats, QStringLiteral("eeq_coordination_number"),
                charges.eeq_coordination_number.present,
                charges.eeq_coordination_number.value);
    auditT1(stats, QStringLiteral("target_T1"), direct.dft_present,
            direct.target.total_decomp.T1);
    auditScalar(stats, QStringLiteral("target_dia_T0"), direct.dft_present,
                direct.target.dia_decomp.T0);
    auditT1(stats, QStringLiteral("target_dia_T1"), direct.dft_present,
            direct.target.dia_decomp.T1);
    auditT2(stats, QStringLiteral("target_dia_T2"), direct.dft_present,
            direct.target.dia_decomp.T2);
    auditScalar(stats, QStringLiteral("target_para_T0"), direct.dft_present,
                direct.target.para_decomp.T0);
    auditT1(stats, QStringLiteral("target_para_T1"), direct.dft_present,
            direct.target.para_decomp.T1);
    auditT2(stats, QStringLiteral("target_para_T2"), direct.dft_present,
            direct.target.para_decomp.T2);
    auditRingPathFeatures(stats, direct);
    auditMethodPathFeatures(stats, direct);
    auditHbondConditioningFeatures(stats, direct);
}

void recordAbsentNewChannelSlabs(const Body& body, PerAtomSubstrateStats& stats) {
    const std::vector<std::pair<QString, ArrayId>> slabs = {
        {QStringLiteral("hbond_scalars"), ArrayId::HbondScalars},
        {QStringLiteral("hbond_nearest_dir"), ArrayId::HbondNearestDirection},
        {QStringLiteral("hbond_flags"), ArrayId::HbondFlags},
        {QStringLiteral("larsen_hbond_count_time_series"), ArrayId::HbondCount},
        {QStringLiteral("hm_shielding_time_series"), ArrayId::HmShielding},
        {QStringLiteral("water_field_time_series/efield"), ArrayId::WaterEfield},
        {QStringLiteral("water_field_time_series/n_first"), ArrayId::WaterNFirst},
        {QStringLiteral("water_field_time_series/n_second"), ArrayId::WaterNSecond},
        {QStringLiteral("hydration_shell_time_series/half_shell_asymmetry"),
         ArrayId::HydrationHalfShellAsymmetry},
        {QStringLiteral("hydration_shell_time_series/mean_water_dipole_cos"),
         ArrayId::HydrationDipoleCos},
        {QStringLiteral("sasa_time_series/sasa"), ArrayId::Sasa},
        {QStringLiteral("hydration_geometry_time_series/surface_normal"), ArrayId::SasaNormal},
        {QStringLiteral("eeq_welford/charge_mean"), ArrayId::EeqChargeMean},
    };
    stats.absent_new_channel_slabs.clear();
    for (const auto& slab : slabs) {
        if (!body.catalog.has(slab.second))
            stats.absent_new_channel_slabs << slab.first;
    }
}

std::array<double, kPerAtomClassicalCols> classicalFeatures(const MechanismAggregate& agg,
                                                            const DirectFeatures& direct) {
    std::array<double, kPerAtomClassicalCols> f;
    fillNaN(f);
    std::size_t c = 0;
    f[c++] = agg.ring_present ? agg.ring_T0 : kNaN;
    for (double v : agg.ring_T2) f[c++] = agg.ring_present ? v : kNaN;
    for (double v : agg.charge_T2) f[c++] = agg.charge_present ? v : kNaN;
    f[c++] = agg.mc_lit_valid_present ? agg.mc_lit_valid.T0 : kNaN;
    for (double v : agg.mc_lit_valid.T2) f[c++] = agg.mc_lit_valid_present ? v : kNaN;
    for (double v : direct.mopac_coulomb_T2) f[c++] = direct.mopac_coulomb_present ? v : kNaN;
    for (double v : direct.mopac_mc_T2) f[c++] = direct.mopac_mc_present ? v : kNaN;
    f[c++] = agg.ff14sb_field_present ? agg.ff14sb_field.x() : kNaN;
    f[c++] = agg.ff14sb_field_present ? agg.ff14sb_field.y() : kNaN;
    f[c++] = agg.ff14sb_field_present ? agg.ff14sb_field.z() : kNaN;
    f[c++] = agg.ff14sb_field_present ? agg.ff14sb_field.norm() : kNaN;
    f[c++] = direct.apbs_E_present ? direct.apbs_E.x() : kNaN;
    f[c++] = direct.apbs_E_present ? direct.apbs_E.y() : kNaN;
    f[c++] = direct.apbs_E_present ? direct.apbs_E.z() : kNaN;
    f[c++] = direct.apbs_E_present ? direct.apbs_E.norm() : kNaN;
    for (double v : direct.apbs_efg_T2) f[c++] = direct.apbs_efg_present ? v : kNaN;
    f[c++] = direct.aimnet2_charge_present ? direct.aimnet2_charge : kNaN;
    f[c++] = direct.aimnet2_crg_present ? direct.aimnet2_crg_scalar : kNaN;
    f[c++] = direct.aimnet2_crg_present ? direct.aimnet2_crg.x() : kNaN;
    f[c++] = direct.aimnet2_crg_present ? direct.aimnet2_crg.y() : kNaN;
    f[c++] = direct.aimnet2_crg_present ? direct.aimnet2_crg.z() : kNaN;
    f[c++] = direct.hbond_flags_present ? direct.hbond_flags[0] : kNaN;
    f[c++] = direct.hbond_nearest_direction_present ? direct.hbond_nearest_direction.x() : kNaN;
    f[c++] = direct.hbond_nearest_direction_present ? direct.hbond_nearest_direction.y() : kNaN;
    f[c++] = direct.hbond_nearest_direction_present ? direct.hbond_nearest_direction.z() : kNaN;
    f[c++] = direct.hbond_count_present ? direct.hbond_count : kNaN;
    f[c++] = direct.hbond_flags_present ? direct.hbond_flags[1] : kNaN;
    f[c++] = direct.hbond_flags_present ? direct.hbond_flags[2] : kNaN;
    for (double v : direct.hm_T2) f[c++] = direct.hm_shielding_present ? v : kNaN;
    f[c++] = direct.water_efield_present ? direct.water_efield.x() : kNaN;
    f[c++] = direct.water_efield_present ? direct.water_efield.y() : kNaN;
    f[c++] = direct.water_efield_present ? direct.water_efield.z() : kNaN;
    f[c++] = direct.water_efield_present ? direct.water_efield.norm() : kNaN;
    f[c++] = direct.water_n_first_present ? direct.water_n_first : kNaN;
    f[c++] = direct.water_n_second_present ? direct.water_n_second : kNaN;
    f[c++] = direct.hydration_half_shell_present ? direct.hydration_half_shell_asymmetry : kNaN;
    f[c++] = direct.hydration_dipole_cos_present ? direct.hydration_dipole_cos : kNaN;
    f[c++] = direct.sasa_present ? direct.sasa : kNaN;
    f[c++] = direct.sasa_normal_present ? direct.sasa_normal.x() : kNaN;
    f[c++] = direct.sasa_normal_present ? direct.sasa_normal.y() : kNaN;
    f[c++] = direct.sasa_normal_present ? direct.sasa_normal.z() : kNaN;
    return f;
}

std::array<double, kPerAtomDriverMagnitudeCols> driverMagnitudes(const MechanismAggregate& agg,
                                                                 const DirectFeatures& direct) {
    return {
        agg.ring_present ? t2Magnitude(agg.ring_T2) : kNaN,
        agg.charge_present ? t2Magnitude(agg.charge_T2) : kNaN,
        agg.mc_lit_valid_present ? t2Magnitude(agg.mc_lit_valid.T2) : kNaN,
        direct.mopac_coulomb_present ? t2Magnitude(direct.mopac_coulomb_T2) : kNaN,
        direct.mopac_mc_present ? t2Magnitude(direct.mopac_mc_T2) : kNaN,
        direct.apbs_efg_present ? t2Magnitude(direct.apbs_efg_T2) : kNaN,
        vecMagnitude(direct.apbs_E, direct.apbs_E_present),
        vecMagnitude(agg.ff14sb_field, agg.ff14sb_field_present),
        vecMagnitude(direct.aimnet2_crg, direct.aimnet2_crg_present),
    };
}

std::array<double, kPerAtomConditioningCols> conditioningFeatures(
    const Body& body, std::size_t atom, std::size_t row, const MechanismAggregate& agg,
    const std::array<double, kPerAtomDriverMagnitudeCols>& mag) {
    std::array<double, kPerAtomConditioningCols> c;
    fillNaN(c);
    std::size_t i = 0;
    c[i++] = nearestDistance(body, CloudKind::RingCenters, atom, row, 1000.0);
    c[i++] = nearestDistance(body, CloudKind::ChargeSites, atom, row, 1000.0);
    c[i++] = nearestDistance(body, CloudKind::BondMidpoints, atom, row, 1000.0);
    c[i++] = nearestDistance(body, CloudKind::Atoms, atom, row, 1000.0, true);
    c[i++] = countNear(body, CloudKind::RingCenters, atom, row, 4.0);
    c[i++] = countNear(body, CloudKind::RingCenters, atom, row, 6.0);
    c[i++] = countNear(body, CloudKind::RingCenters, atom, row, 8.0);
    c[i++] = countNear(body, CloudKind::ChargeSites, atom, row, 4.0);
    c[i++] = countNear(body, CloudKind::ChargeSites, atom, row, 6.0);
    c[i++] = countNear(body, CloudKind::ChargeSites, atom, row, 10.0);
    c[i++] = countNear(body, CloudKind::BondMidpoints, atom, row, 4.0);
    c[i++] = countNear(body, CloudKind::BondMidpoints, atom, row, 8.0);
    c[i++] = countNear(body, CloudKind::BondMidpoints, atom, row, 10.0);
    c[i++] = agg.ring_self_or_bonded_n;
    c[i++] = agg.bond_self_or_bonded_n;
    c[i++] = agg.charge_excluded_same_residue_n;
    c[i++] = (agg.ring_self_or_bonded_n + agg.bond_self_or_bonded_n
              + agg.charge_excluded_same_residue_n) > 0 ? 1.0 : 0.0;
    for (double v : mag) c[i++] = finiteOrZero(v);
    c[i++] = agg.isolation.gap_to_2nd_ring_r;
    c[i++] = agg.isolation.gap_to_2nd_charge_r;
    c[i++] = agg.isolation.gap_to_2nd_bond_r;
    c[i++] = agg.isolation.dominant_fraction_ring;
    c[i++] = agg.isolation.dominant_fraction_charge;
    c[i++] = agg.isolation.dominant_fraction_mc;
    return c;
}

std::array<double, kPerAtomDominanceCols> dominanceFeatures(const PerAtomIsolationScalars& iso) {
    std::array<double, kPerAtomDominanceCols> out;
    fillNaN(out);
    std::size_t i = 0;
    out[i++] = iso.gap_to_2nd_ring_r;
    out[i++] = iso.gap_to_2nd_charge_r;
    out[i++] = iso.gap_to_2nd_bond_r;
    out[i++] = iso.gap_to_2nd_field_r;
    out[i++] = iso.dominant_fraction_ring;
    out[i++] = iso.dominant_fraction_charge;
    out[i++] = iso.dominant_fraction_mc;
    out[i++] = iso.dominant_fraction_field;
    return out;
}

void appendTensor9Names(QStringList& cols, const QString& prefix) {
    cols << QStringLiteral("%1_T0").arg(prefix);
    for (int i = 0; i < 3; ++i) cols << QStringLiteral("%1_T1_%2").arg(prefix).arg(i);
    for (int i = 0; i < 5; ++i) cols << QStringLiteral("%1_T2_%2").arg(prefix).arg(i);
}

void appendIndexedNames(QStringList& cols, const QString& prefix, int count) {
    for (int i = 0; i < count; ++i) cols << QStringLiteral("%1_%2").arg(prefix).arg(i);
}

QStringList makeRingPathColumns() {
    QStringList cols;
    appendTensor9Names(cols, QStringLiteral("bs_shielding"));
    appendIndexedNames(cols, QStringLiteral("bs_per_type_T0"), 8);
    appendIndexedNames(cols, QStringLiteral("bs_per_type_T2"), 40);
    cols << QStringLiteral("bs_total_B_x") << QStringLiteral("bs_total_B_y")
         << QStringLiteral("bs_total_B_z");
    appendIndexedNames(cols, QStringLiteral("bs_ring_counts"), 4);
    appendTensor9Names(cols, QStringLiteral("hm_shielding"));
    appendIndexedNames(cols, QStringLiteral("hm_per_type_T0"), 8);
    appendIndexedNames(cols, QStringLiteral("hm_per_type_T2"), 40);
    appendIndexedNames(cols, QStringLiteral("pq_per_type_T0"), 8);
    appendIndexedNames(cols, QStringLiteral("disp_per_type_T0"), 8);
    return cols;
}

QStringList makeMethodPathColumns() {
    QStringList cols;
    appendIndexedNames(cols, QStringLiteral("mc_scalars"), 6);
    appendIndexedNames(cols, QStringLiteral("mopac_mc_scalars"), 6);
    cols << QStringLiteral("mopac_bond_orders_incident_count")
         << QStringLiteral("mopac_bond_orders_incident_sum")
         << QStringLiteral("mopac_bond_orders_incident_mean")
         << QStringLiteral("mopac_bond_orders_incident_max");
    cols << QStringLiteral("mopac_coulomb_E_x")
         << QStringLiteral("mopac_coulomb_E_y")
         << QStringLiteral("mopac_coulomb_E_z");
    appendIndexedNames(cols, QStringLiteral("mopac_coulomb_efg_backbone"), 5);
    appendIndexedNames(cols, QStringLiteral("mopac_coulomb_efg_aromatic"), 5);
    appendIndexedNames(cols, QStringLiteral("mopac_coulomb_scalars"), 4);
    appendIndexedNames(cols, QStringLiteral("aimnet2_efg"), 5);
    appendIndexedNames(cols, QStringLiteral("aimnet2_efg_aromatic"), 5);
    appendIndexedNames(cols, QStringLiteral("aimnet2_efg_backbone"), 5);
    appendIndexedNames(cols, QStringLiteral("water_efg"), 5);
    cols << QStringLiteral("water_efield_first_x")
         << QStringLiteral("water_efield_first_y")
         << QStringLiteral("water_efield_first_z");
    cols << QStringLiteral("eeq_cn");
    appendIndexedNames(cols, QStringLiteral("mopac_scalars"), 4);
    return cols;
}

QStringList makeHbondConditioningColumns() {
    QStringList cols;
    appendIndexedNames(cols, QStringLiteral("larsen_hbond_1pHB_T2"), 5);
    appendIndexedNames(cols, QStringLiteral("larsen_hbond_2pHB_T2"), 5);
    appendIndexedNames(cols, QStringLiteral("larsen_hbond_1pHaB_T2"), 5);
    appendIndexedNames(cols, QStringLiteral("larsen_hbond_2pHaB_T2"), 5);
    cols << QStringLiteral("larsen_hbond_water_term");
    appendIndexedNames(cols, QStringLiteral("hbond_scalars"), 4);
    cols << QStringLiteral("dssp_chem_donor_flag")
         << QStringLiteral("dssp_chem_acceptor_flag");
    appendIndexedNames(cols, QStringLiteral("dssp_hbond_energy"), 4);
    cols << QStringLiteral("dssp_acceptor_partner_0")
         << QStringLiteral("dssp_acceptor_energy_0")
         << QStringLiteral("dssp_acceptor_partner_1")
         << QStringLiteral("dssp_acceptor_energy_1")
         << QStringLiteral("dssp_acceptor_count")
         << QStringLiteral("dssp_donor_partner_0")
         << QStringLiteral("dssp_donor_energy_0")
         << QStringLiteral("dssp_donor_partner_1")
         << QStringLiteral("dssp_donor_energy_1")
         << QStringLiteral("dssp_donor_count");
    appendIndexedNames(cols, QStringLiteral("dssp_ss8"), 8);
    appendIndexedNames(cols, QStringLiteral("dssp_chi"), 12);
    cols << QStringLiteral("omega_actual");
    cols << QStringLiteral("pyramidalization");
    appendIndexedNames(cols, QStringLiteral("ring_geometry"), 10);
    return cols;
}

const QStringList kClassicalColumns = {
    QStringLiteral("ring_jb_T0"),
    QStringLiteral("ring_jb_T2_0"), QStringLiteral("ring_jb_T2_1"),
    QStringLiteral("ring_jb_T2_2"), QStringLiteral("ring_jb_T2_3"),
    QStringLiteral("ring_jb_T2_4"),
    QStringLiteral("charge_q_over_r3_T2_0"), QStringLiteral("charge_q_over_r3_T2_1"),
    QStringLiteral("charge_q_over_r3_T2_2"), QStringLiteral("charge_q_over_r3_T2_3"),
    QStringLiteral("charge_q_over_r3_T2_4"),
    QStringLiteral("mc_lit_T0_valid"),
    QStringLiteral("mc_lit_T2_valid_0"), QStringLiteral("mc_lit_T2_valid_1"),
    QStringLiteral("mc_lit_T2_valid_2"), QStringLiteral("mc_lit_T2_valid_3"),
    QStringLiteral("mc_lit_T2_valid_4"),
    QStringLiteral("mopac_coulomb_shielding_T2_0"), QStringLiteral("mopac_coulomb_shielding_T2_1"),
    QStringLiteral("mopac_coulomb_shielding_T2_2"), QStringLiteral("mopac_coulomb_shielding_T2_3"),
    QStringLiteral("mopac_coulomb_shielding_T2_4"),
    QStringLiteral("mopac_mc_shielding_T2_0"), QStringLiteral("mopac_mc_shielding_T2_1"),
    QStringLiteral("mopac_mc_shielding_T2_2"), QStringLiteral("mopac_mc_shielding_T2_3"),
    QStringLiteral("mopac_mc_shielding_T2_4"),
    QStringLiteral("ff14sb_field_x"), QStringLiteral("ff14sb_field_y"),
    QStringLiteral("ff14sb_field_z"), QStringLiteral("ff14sb_field_mag"),
    QStringLiteral("apbs_E_x"), QStringLiteral("apbs_E_y"),
    QStringLiteral("apbs_E_z"), QStringLiteral("apbs_E_mag"),
    QStringLiteral("apbs_efg_T2_0"), QStringLiteral("apbs_efg_T2_1"),
    QStringLiteral("apbs_efg_T2_2"), QStringLiteral("apbs_efg_T2_3"),
    QStringLiteral("apbs_efg_T2_4"),
    QStringLiteral("aimnet2_charge"), QStringLiteral("aimnet2_crg_scalar"),
    QStringLiteral("aimnet2_crg_x"), QStringLiteral("aimnet2_crg_y"),
    QStringLiteral("aimnet2_crg_z"),
    QStringLiteral("hbond_flags_is_backbone"), QStringLiteral("hbond_nearest_dir_x"),
    QStringLiteral("hbond_nearest_dir_y"), QStringLiteral("hbond_nearest_dir_z"),
    QStringLiteral("hbond_count"), QStringLiteral("hbond_flags_is_donor"),
    QStringLiteral("hbond_flags_is_acceptor"),
    QStringLiteral("hm_shielding_T2_0"), QStringLiteral("hm_shielding_T2_1"),
    QStringLiteral("hm_shielding_T2_2"), QStringLiteral("hm_shielding_T2_3"),
    QStringLiteral("hm_shielding_T2_4"),
    QStringLiteral("water_efield_x"), QStringLiteral("water_efield_y"),
    QStringLiteral("water_efield_z"), QStringLiteral("water_efield_mag"),
    QStringLiteral("water_n_first"), QStringLiteral("water_n_second"),
    QStringLiteral("water_half_shell_asymmetry"), QStringLiteral("water_dipole_cos"),
    QStringLiteral("sasa"), QStringLiteral("sasa_normal_x"),
    QStringLiteral("sasa_normal_y"), QStringLiteral("sasa_normal_z"),
};

const QStringList kConditioningColumns = {
    QStringLiteral("nearest_ring_r"), QStringLiteral("nearest_charge_r"),
    QStringLiteral("nearest_bond_midpoint_r"), QStringLiteral("nearest_heavy_atom_r"),
    QStringLiteral("ring_count_4A"), QStringLiteral("ring_count_6A"), QStringLiteral("ring_count_8A"),
    QStringLiteral("charge_count_4A"), QStringLiteral("charge_count_6A"), QStringLiteral("charge_count_10A"),
    QStringLiteral("bond_count_4A"), QStringLiteral("bond_count_8A"), QStringLiteral("bond_count_10A"),
    QStringLiteral("ring_self_or_bonded_n"), QStringLiteral("bond_self_or_bonded_n"),
    QStringLiteral("charge_excluded_same_residue_n"), QStringLiteral("has_self_or_bonded_driver"),
    QStringLiteral("abs_ring_jb_T2"), QStringLiteral("abs_charge_T2"), QStringLiteral("abs_mc_lit_T2"),
    QStringLiteral("abs_mopac_coulomb_T2"), QStringLiteral("abs_mopac_mc_T2"),
    QStringLiteral("abs_apbs_efg_T2"), QStringLiteral("abs_apbs_E"),
    QStringLiteral("abs_ff14sb_E"), QStringLiteral("abs_aimnet2_crg"),
    QStringLiteral("gap_to_2nd_ring_r"), QStringLiteral("gap_to_2nd_charge_r"),
    QStringLiteral("gap_to_2nd_bond_r"), QStringLiteral("dominant_fraction_ring"),
    QStringLiteral("dominant_fraction_charge"), QStringLiteral("dominant_fraction_mc"),
};

const QStringList kDominanceColumns = {
    QStringLiteral("gap_to_2nd_ring_r"),
    QStringLiteral("gap_to_2nd_charge_r"),
    QStringLiteral("gap_to_2nd_bond_r"),
    QStringLiteral("gap_to_2nd_field_r"),
    QStringLiteral("dominant_fraction_ring"),
    QStringLiteral("dominant_fraction_charge"),
    QStringLiteral("dominant_fraction_mc"),
    QStringLiteral("dominant_fraction_field"),
};

const QStringList kMagnitudeColumns = {
    QStringLiteral("sd_ring_jb_T2_by_atom"), QStringLiteral("sd_charge_T2_by_atom"),
    QStringLiteral("sd_mc_lit_T2_by_atom"), QStringLiteral("sd_mopac_coulomb_T2_by_atom"),
    QStringLiteral("sd_mopac_mc_T2_by_atom"), QStringLiteral("sd_apbs_efg_T2_by_atom"),
    QStringLiteral("sd_apbs_E_by_atom"), QStringLiteral("sd_ff14sb_E_by_atom"),
    QStringLiteral("sd_aimnet2_crg_by_atom"),
};

const QStringList kPartitionBinColumns = {
    QStringLiteral("bin_nearest_ring_r_4_6_8_10"),
    QStringLiteral("bin_nearest_charge_r_4_6_8_10"),
    QStringLiteral("bin_nearest_bond_midpoint_r_4_6_8_10"),
    QStringLiteral("bin_nearest_heavy_atom_r_4_6_8_10"),
    QStringLiteral("bin_gap_to_2nd_ring_r_4_6_8_10"),
    QStringLiteral("bin_gap_to_2nd_charge_r_4_6_8_10"),
    QStringLiteral("bin_gap_to_2nd_bond_r_4_6_8_10"),
    QStringLiteral("bin_abs_ring_jb_T2_quintile"),
    QStringLiteral("bin_abs_charge_T2_quintile"),
    QStringLiteral("bin_abs_mc_lit_T2_quintile"),
    QStringLiteral("bin_abs_mopac_coulomb_T2_quintile"),
    QStringLiteral("bin_abs_mopac_mc_T2_quintile"),
    QStringLiteral("bin_abs_apbs_efg_T2_quintile"),
    QStringLiteral("bin_abs_apbs_E_quintile"),
    QStringLiteral("bin_abs_ff14sb_E_quintile"),
    QStringLiteral("bin_abs_aimnet2_crg_quintile"),
    QStringLiteral("bin_sd_ring_jb_T2_by_atom_quintile"),
    QStringLiteral("bin_sd_charge_T2_by_atom_quintile"),
    QStringLiteral("bin_sd_mc_lit_T2_by_atom_quintile"),
    QStringLiteral("bin_sd_mopac_coulomb_T2_by_atom_quintile"),
    QStringLiteral("bin_sd_mopac_mc_T2_by_atom_quintile"),
    QStringLiteral("bin_sd_apbs_efg_T2_by_atom_quintile"),
    QStringLiteral("bin_sd_apbs_E_by_atom_quintile"),
    QStringLiteral("bin_sd_ff14sb_E_by_atom_quintile"),
    QStringLiteral("bin_sd_aimnet2_crg_by_atom_quintile"),
};

const QStringList kDominanceBinColumns = {
    QStringLiteral("bin_dominant_fraction_ring_quintile"),
    QStringLiteral("bin_dominant_fraction_charge_quintile"),
    QStringLiteral("bin_dominant_fraction_mc_quintile"),
    QStringLiteral("bin_dominant_fraction_field_quintile"),
};

const QStringList kBackboneAuditColumns = {
    QStringLiteral("broad_charge_literature_kernel_T2_0"),
    QStringLiteral("broad_charge_literature_kernel_T2_1"),
    QStringLiteral("broad_charge_literature_kernel_T2_2"),
    QStringLiteral("broad_charge_literature_kernel_T2_3"),
    QStringLiteral("broad_charge_literature_kernel_T2_4"),
    QStringLiteral("broad_field_local_x"),
    QStringLiteral("broad_field_local_y"),
    QStringLiteral("broad_field_local_z"),
    QStringLiteral("broad_mc_lit_T0_valid"),
    QStringLiteral("broad_mc_lit_T2_valid_0"),
    QStringLiteral("broad_mc_lit_T2_valid_1"),
    QStringLiteral("broad_mc_lit_T2_valid_2"),
    QStringLiteral("broad_mc_lit_T2_valid_3"),
    QStringLiteral("broad_mc_lit_T2_valid_4"),
};

const QStringList kRingPathColumns = makeRingPathColumns();
const QStringList kMethodPathColumns = makeMethodPathColumns();
const QStringList kHbondConditioningColumns = makeHbondConditioningColumns();

void validateColumnCounts() {
    auto check = [](const QString& name, qsizetype actual, std::size_t expected) {
        if (actual != static_cast<qsizetype>(expected)) {
            throw std::runtime_error(QStringLiteral("%1 column count %2 != expected %3")
                                         .arg(name)
                                         .arg(actual)
                                         .arg(static_cast<qulonglong>(expected))
                                         .toStdString());
        }
    };
    check(QStringLiteral("per_atom_substrate_features_classical"),
          kClassicalColumns.size(), kPerAtomClassicalCols);
    check(QStringLiteral("per_atom_substrate_features_conditioning"),
          kConditioningColumns.size(), kPerAtomConditioningCols);
    check(QStringLiteral("per_atom_substrate_features_dominance"),
          kDominanceColumns.size(), kPerAtomDominanceCols);
    check(QStringLiteral("per_atom_substrate_driver_modulation_by_atom"),
          kMagnitudeColumns.size(), kPerAtomDriverMagnitudeCols);
    check(QStringLiteral("per_atom_substrate_partition_bins"),
          kPartitionBinColumns.size(), kPerAtomPartitionBinCols);
    check(QStringLiteral("per_atom_substrate_dominance_bins"),
          kDominanceBinColumns.size(), kPerAtomDominanceBinCols);
    check(QStringLiteral("per_atom_substrate_backbone_audit"),
          kBackboneAuditColumns.size(), kPerAtomBackboneAuditCols);
    check(QStringLiteral("per_atom_substrate_features_ring_paths"),
          kRingPathColumns.size(), kPerAtomRingPathCols);
    check(QStringLiteral("per_atom_substrate_features_method_paths"),
          kMethodPathColumns.size(), kPerAtomMethodPathCols);
    check(QStringLiteral("per_atom_substrate_features_hbond_conditioning"),
	          kHbondConditioningColumns.size(), kPerAtomHbondConditioningCols);
}

struct RowPartitionScratch {
    std::size_t atom = 0;
    int stratum = 0;
    std::array<double, 7> geometry = {};
    std::array<double, kPerAtomDriverMagnitudeCols> magnitude = {};
};

struct RowDominanceScratch {
    int stratum = 0;
    std::array<double, kPerAtomDominanceBinCols> fraction = {};
};

struct QuintileEdges {
    std::array<double, 4> edge = {kNaN, kNaN, kNaN, kNaN};
    bool valid = false;
};

int distanceBandId(double v) {
    if (!std::isfinite(v)) return -1;
    if (v < 4.0) return 0;
    if (v < 6.0) return 1;
    if (v < 8.0) return 2;
    if (v < 10.0) return 3;
    return 4;
}

int chargeGapBandId(double v) {
    if (!std::isfinite(v)) return -1;
    if (v < 0.25) return 0;
    if (v < 0.5) return 1;
    if (v < 1.0) return 2;
    if (v < 1.5) return 3;
    return 4;
}

QuintileEdges computeQuintileEdges(std::vector<double> values) {
    QuintileEdges out;
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double v) { return !std::isfinite(v); }),
                 values.end());
    if (values.empty()) return out;
    std::sort(values.begin(), values.end());
    const std::array<double, 4> probs = {0.2, 0.4, 0.6, 0.8};
    for (std::size_t i = 0; i < probs.size(); ++i) {
        const double pos = probs[i] * static_cast<double>(values.size() - 1);
        out.edge[i] = values[static_cast<std::size_t>(std::llround(pos))];
    }
    out.valid = true;
    return out;
}

int quintileBinId(double v, const QuintileEdges& edges) {
    if (!edges.valid || !std::isfinite(v)) return -1;
    for (std::size_t i = 0; i < edges.edge.size(); ++i) {
        if (v <= edges.edge[i]) return static_cast<int>(i);
    }
    return 4;
}

QJsonArray jsonArrayDoubles(const std::array<double, 4>& values, bool valid) {
    QJsonArray arr;
    for (double v : values)
        arr.append(valid && std::isfinite(v) ? QJsonValue(v) : QJsonValue(QJsonValue::Null));
    return arr;
}

bool writePartitionBins(const QString& outDir,
                        const std::vector<RowPartitionScratch>& rows,
                        const std::vector<WelfordCell>& modulation,
                        std::size_t atoms,
                        PerAtomSubstrateStats& stats) {
    constexpr int kStratumCount = static_cast<int>(StratumOrd::Other) + 1;
    using FeatureBuckets = std::array<std::vector<double>, kPerAtomDriverMagnitudeCols>;
    std::array<FeatureBuckets, kStratumCount> magnitudeBuckets;
    std::array<FeatureBuckets, kStratumCount> modulationBuckets;
    std::vector<int> atomStrata(atoms, -1);

    for (const RowPartitionScratch& row : rows) {
        const int s = std::clamp(row.stratum, 0, kStratumCount - 1);
        if (row.atom < atomStrata.size() && atomStrata[row.atom] < 0) atomStrata[row.atom] = s;
        for (std::size_t c = 0; c < kPerAtomDriverMagnitudeCols; ++c)
            if (std::isfinite(row.magnitude[c])) magnitudeBuckets[s][c].push_back(row.magnitude[c]);
    }
    for (std::size_t atom = 0; atom < atoms; ++atom) {
        const int s = atomStrata[atom] >= 0 ? atomStrata[atom] : kStratumCount - 1;
        for (std::size_t c = 0; c < kPerAtomDriverMagnitudeCols; ++c) {
            const double v = modulation[atom * kPerAtomDriverMagnitudeCols + c].sd();
            if (std::isfinite(v)) modulationBuckets[s][c].push_back(v);
        }
    }

    std::array<std::array<QuintileEdges, kPerAtomDriverMagnitudeCols>, kStratumCount> magnitudeEdges;
    std::array<std::array<QuintileEdges, kPerAtomDriverMagnitudeCols>, kStratumCount> modulationEdges;
    for (int s = 0; s < kStratumCount; ++s) {
        for (std::size_t c = 0; c < kPerAtomDriverMagnitudeCols; ++c) {
            magnitudeEdges[s][c] = computeQuintileEdges(std::move(magnitudeBuckets[s][c]));
            modulationEdges[s][c] = computeQuintileEdges(std::move(modulationBuckets[s][c]));
        }
    }

    std::vector<int32_t> out;
    out.reserve(rows.size() * kPerAtomPartitionBinCols);
    std::array<int32_t, kPerAtomPartitionBinCols> minBin;
    std::array<int32_t, kPerAtomPartitionBinCols> maxBin;
    std::array<std::size_t, kPerAtomPartitionBinCols> present;
    minBin.fill(std::numeric_limits<int32_t>::max());
    maxBin.fill(std::numeric_limits<int32_t>::min());
    present.fill(0);
    auto pushBin = [&](std::size_t col, int32_t bin) {
        out.push_back(bin);
        if (bin < 0) return;
        ++present[col];
        minBin[col] = std::min(minBin[col], bin);
        maxBin[col] = std::max(maxBin[col], bin);
    };

    for (const RowPartitionScratch& row : rows) {
        std::size_t col = 0;
        const int s = std::clamp(row.stratum, 0, kStratumCount - 1);
        for (std::size_t g = 0; g < row.geometry.size(); ++g) {
            const int bin = g == 5 ? chargeGapBandId(row.geometry[g])
                                   : distanceBandId(row.geometry[g]);
            pushBin(col++, bin);
        }
        for (std::size_t c = 0; c < kPerAtomDriverMagnitudeCols; ++c)
            pushBin(col++, quintileBinId(row.magnitude[c], magnitudeEdges[s][c]));
        for (std::size_t c = 0; c < kPerAtomDriverMagnitudeCols; ++c) {
            const double v = row.atom < atoms ? modulation[row.atom * kPerAtomDriverMagnitudeCols + c].sd()
                                              : kNaN;
            pushBin(col++, quintileBinId(v, modulationEdges[s][c]));
        }
    }
    if (!writeNpy<int32_t>(QStringLiteral("%1/per_atom_substrate_partition_bins.npy").arg(outDir),
                           {rows.size(), kPerAtomPartitionBinCols}, out, QByteArray("<i4"))) {
        return false;
    }

    for (int i = 0; i < kPartitionBinColumns.size(); ++i) {
        PerAtomChannelAudit& audit = stats.new_channel_audit[kPartitionBinColumns[i]];
        audit.present = present[static_cast<std::size_t>(i)];
        if (present[static_cast<std::size_t>(i)] > 0) {
            audit.has_range = true;
            audit.min = minBin[static_cast<std::size_t>(i)];
            audit.max = maxBin[static_cast<std::size_t>(i)];
        }
    }

    QJsonObject manifest;
    manifest.insert(QStringLiteral("array"), QStringLiteral("per_atom_substrate_partition_bins"));
    manifest.insert(QStringLiteral("dtype"), QStringLiteral("int32"));
    manifest.insert(QStringLiteral("missing_bin_id"), -1);
    manifest.insert(QStringLiteral("distance_band_edges_A"), QJsonArray{4.0, 6.0, 8.0, 10.0});
    manifest.insert(QStringLiteral("charge_gap_band_edges_A"), QJsonArray{0.25, 0.5, 1.0, 1.5});
    QJsonArray shape;
    shape.append(static_cast<qint64>(rows.size()));
    shape.append(static_cast<qint64>(kPerAtomPartitionBinCols));
    manifest.insert(QStringLiteral("shape"), shape);
    QJsonArray columns;
    for (const QString& name : kPartitionBinColumns) columns.append(name);
    manifest.insert(QStringLiteral("columns"), columns);

    QJsonObject strata;
    for (int s = 0; s < kStratumCount; ++s) {
        QJsonObject stratum;
        QJsonObject mag;
        QJsonObject mod;
        for (std::size_t c = 0; c < kPerAtomDriverMagnitudeCols; ++c) {
            mag.insert(kConditioningColumns[17 + static_cast<int>(c)],
                       jsonArrayDoubles(magnitudeEdges[s][c].edge, magnitudeEdges[s][c].valid));
            mod.insert(kMagnitudeColumns[static_cast<int>(c)],
                       jsonArrayDoubles(modulationEdges[s][c].edge, modulationEdges[s][c].valid));
        }
        stratum.insert(QStringLiteral("driver_magnitude_quintiles"), mag);
        stratum.insert(QStringLiteral("driver_modulation_quintiles"), mod);
        strata.insert(stratumName(static_cast<StratumOrd>(s)), stratum);
    }
    manifest.insert(QStringLiteral("quintile_edges_by_stratum"), strata);
    stats.partition_bin_manifest = manifest;
    return true;
}

bool writeDominanceBins(const QString& outDir,
                        const std::vector<RowDominanceScratch>& rows,
                        PerAtomSubstrateStats& stats) {
    constexpr int kStratumCount = static_cast<int>(StratumOrd::Other) + 1;
    using DominanceBuckets = std::array<std::vector<double>, kPerAtomDominanceBinCols>;
    std::array<DominanceBuckets, kStratumCount> buckets;
    for (const RowDominanceScratch& row : rows) {
        const int s = std::clamp(row.stratum, 0, kStratumCount - 1);
        for (std::size_t c = 0; c < kPerAtomDominanceBinCols; ++c) {
            if (std::isfinite(row.fraction[c])) buckets[s][c].push_back(row.fraction[c]);
        }
    }

    std::array<std::array<QuintileEdges, kPerAtomDominanceBinCols>, kStratumCount> edges;
    for (int s = 0; s < kStratumCount; ++s)
        for (std::size_t c = 0; c < kPerAtomDominanceBinCols; ++c)
            edges[s][c] = computeQuintileEdges(std::move(buckets[s][c]));

    std::vector<int32_t> out;
    out.reserve(rows.size() * kPerAtomDominanceBinCols);
    std::array<int32_t, kPerAtomDominanceBinCols> minBin;
    std::array<int32_t, kPerAtomDominanceBinCols> maxBin;
    std::array<std::size_t, kPerAtomDominanceBinCols> present;
    minBin.fill(std::numeric_limits<int32_t>::max());
    maxBin.fill(std::numeric_limits<int32_t>::min());
    present.fill(0);
    auto pushBin = [&](std::size_t col, int32_t bin) {
        out.push_back(bin);
        if (bin < 0) return;
        ++present[col];
        minBin[col] = std::min(minBin[col], bin);
        maxBin[col] = std::max(maxBin[col], bin);
    };

    for (const RowDominanceScratch& row : rows) {
        const int s = std::clamp(row.stratum, 0, kStratumCount - 1);
        for (std::size_t c = 0; c < kPerAtomDominanceBinCols; ++c)
            pushBin(c, quintileBinId(row.fraction[c], edges[s][c]));
    }
    if (!writeNpy<int32_t>(QStringLiteral("%1/per_atom_substrate_dominance_bins.npy").arg(outDir),
                           {rows.size(), kPerAtomDominanceBinCols}, out, QByteArray("<i4"))) {
        return false;
    }

    for (int i = 0; i < kDominanceBinColumns.size(); ++i) {
        PerAtomChannelAudit& audit = stats.new_channel_audit[kDominanceBinColumns[i]];
        audit.present = present[static_cast<std::size_t>(i)];
        if (present[static_cast<std::size_t>(i)] > 0) {
            audit.has_range = true;
            audit.min = minBin[static_cast<std::size_t>(i)];
            audit.max = maxBin[static_cast<std::size_t>(i)];
        }
    }

    QJsonObject manifest;
    manifest.insert(QStringLiteral("array"), QStringLiteral("per_atom_substrate_dominance_bins"));
    manifest.insert(QStringLiteral("dtype"), QStringLiteral("int32"));
    manifest.insert(QStringLiteral("missing_bin_id"), -1);
    manifest.insert(QStringLiteral("source_scalar_array"),
                    QStringLiteral("per_atom_substrate_features_dominance"));
    QJsonArray shape;
    shape.append(static_cast<qint64>(rows.size()));
    shape.append(static_cast<qint64>(kPerAtomDominanceBinCols));
    manifest.insert(QStringLiteral("shape"), shape);
    QJsonArray columns;
    for (const QString& name : kDominanceBinColumns) columns.append(name);
    manifest.insert(QStringLiteral("columns"), columns);
    QJsonArray sourceCols;
    for (int i = static_cast<int>(kPerAtomDominanceGapCols); i < kDominanceColumns.size(); ++i)
        sourceCols.append(kDominanceColumns[i]);
    manifest.insert(QStringLiteral("source_scalar_columns"), sourceCols);

    QJsonObject strata;
    for (int s = 0; s < kStratumCount; ++s) {
        QJsonObject stratum;
        for (std::size_t c = 0; c < kPerAtomDominanceBinCols; ++c) {
            stratum.insert(kDominanceColumns[static_cast<int>(kPerAtomDominanceGapCols + c)],
                           jsonArrayDoubles(edges[s][c].edge, edges[s][c].valid));
        }
        strata.insert(stratumName(static_cast<StratumOrd>(s)), stratum);
    }
    manifest.insert(QStringLiteral("quintile_edges_by_stratum"), strata);
    stats.dominance_bin_manifest = manifest;
    return true;
}

class PerAtomWriter {
public:
    PerAtomWriter(const QString& outDir, std::size_t rows, std::size_t atoms,
                  const PerAtomSubstrateConfig& cfg)
        : outDir_(outDir), rows_(rows), atoms_(atoms), cfg_(cfg) {
        QDir().mkpath(outDir_);
        rowsFile_ = std::make_unique<QSaveFile>(QStringLiteral("%1/per_atom_substrate_rows.csv").arg(outDir_));
        ringIdentityFile_ =
            std::make_unique<QSaveFile>(QStringLiteral("%1/per_atom_substrate_ring_identity.csv").arg(outDir_));
        if (!rowsFile_->open(QIODevice::WriteOnly | QIODevice::Text)) return;
        if (!ringIdentityFile_->open(QIODevice::WriteOnly | QIODevice::Text)) return;
        rowsOut_ = std::make_unique<QTextStream>(rowsFile_.get());
        ringIdentityOut_ = std::make_unique<QTextStream>(ringIdentityFile_.get());
        *rowsOut_ << rowHeader() << '\n';
        *ringIdentityOut_ << "atom_index,membership_slot,ring_id,ring_type_index_ord,ring_kind_ord,"
                              "fused_partner_ring_id,ring_atom_count,ring_atom_order,is_vertex,is_substituent\n";

        ok_ = targetT2_.open(path(QStringLiteral("per_atom_substrate_target_T2.npy")),
                             {rows_, 5}, QByteArray("<f8"))
              && targetT0_.open(path(QStringLiteral("per_atom_substrate_target_T0.npy")),
                                {rows_}, QByteArray("<f8"))
              && targetT1_.open(path(QStringLiteral("per_atom_substrate_target_T1.npy")),
                                {rows_, 3}, QByteArray("<f8"))
              && targetDiaT0_.open(path(QStringLiteral("per_atom_substrate_target_dia_T0.npy")),
                                   {rows_}, QByteArray("<f8"))
              && targetDiaT1_.open(path(QStringLiteral("per_atom_substrate_target_dia_T1.npy")),
                                   {rows_, 3}, QByteArray("<f8"))
              && targetDiaT2_.open(path(QStringLiteral("per_atom_substrate_target_dia_T2.npy")),
                                   {rows_, 5}, QByteArray("<f8"))
              && targetParaT0_.open(path(QStringLiteral("per_atom_substrate_target_para_T0.npy")),
                                    {rows_}, QByteArray("<f8"))
              && targetParaT1_.open(path(QStringLiteral("per_atom_substrate_target_para_T1.npy")),
                                    {rows_, 3}, QByteArray("<f8"))
              && targetParaT2_.open(path(QStringLiteral("per_atom_substrate_target_para_T2.npy")),
                                    {rows_, 5}, QByteArray("<f8"))
              && classical_.open(path(QStringLiteral("per_atom_substrate_features_classical.npy")),
                                 {rows_, kPerAtomClassicalCols}, QByteArray("<f8"))
              && ringPaths_.open(path(QStringLiteral("per_atom_substrate_features_ring_paths.npy")),
                                 {rows_, kPerAtomRingPathCols}, QByteArray("<f8"))
              && methodPaths_.open(path(QStringLiteral("per_atom_substrate_features_method_paths.npy")),
                                   {rows_, kPerAtomMethodPathCols}, QByteArray("<f8"))
              && hbondConditioning_.open(path(QStringLiteral("per_atom_substrate_features_hbond_conditioning.npy")),
                                         {rows_, kPerAtomHbondConditioningCols}, QByteArray("<f8"))
              && conditioning_.open(path(QStringLiteral("per_atom_substrate_features_conditioning.npy")),
                                    {rows_, kPerAtomConditioningCols}, QByteArray("<f8"))
              && dominance_.open(path(QStringLiteral("per_atom_substrate_features_dominance.npy")),
                                 {rows_, kPerAtomDominanceCols}, QByteArray("<f8"))
              && backboneAudit_.open(path(QStringLiteral("per_atom_substrate_backbone_audit.npy")),
                                     {rows_, kPerAtomBackboneAuditCols}, QByteArray("<f8"));
        if (ok_ && cfg_.emit_embedding) {
            ok_ = embedding_.open(path(QStringLiteral("per_atom_substrate_aimnet2_embedding.npy")),
                                  {rows_, cfg_.embedding_dims}, QByteArray("<f4"));
        }
    }

    bool ok() const { return ok_; }

    void writeRingIdentity(const Body& body) {
        if (!ok_ || !body.run.protein) return;
        const model::QtTopology& topo = body.run.protein->topology();
        for (std::size_t atom = 0; atom < body.run.protein->atomCount(); ++atom) {
            int slot = 0;
            for (int32_t mi : topo.ringMembershipsForAtom(atom)) {
                const model::QtRingMembership& m = topo.ringMembershipAt(static_cast<std::size_t>(mi));
                if (m.ringId < 0 || static_cast<std::size_t>(m.ringId) >= topo.ringCount()) continue;
                const model::QtRing& r = topo.ringAt(static_cast<std::size_t>(m.ringId));
                *ringIdentityOut_ << static_cast<qulonglong>(atom) << ',' << slot++ << ','
                                  << m.ringId << ',' << r.TypeIndexAsInt() << ','
                                  << static_cast<int>(r.ringKind) << ',' << r.fusedPartnerRingId
                                  << ',' << static_cast<int>(r.atomIndices.size()) << ','
                                  << static_cast<int>(m.ringAtomOrder) << ','
                                  << (m.isVertex ? 1 : 0) << ',' << (m.isSubstituent ? 1 : 0)
                                  << '\n';
            }
        }
    }

    bool writeRow(const Body& body, std::size_t atom, std::size_t row, std::size_t frameSlot,
                  std::size_t orig, const DirectFeatures& direct,
                  const MechanismAggregate& agg,
                  const RowChargeScalars& charges,
                  const std::array<double, kPerAtomClassicalCols>& classical,
                  const std::array<double, kPerAtomConditioningCols>& conditioning,
                  const std::array<double, kPerAtomDominanceCols>& dominance) {
        if (!ok_ || !body.run.protein) return false;
        const model::QtProtein& p = *body.run.protein;
        const model::QtAtom& a = p.atom(atom);
        const RoleOrd role = roleForAtom(a);
        const StratumOrd stratum = stratumForAtom(a);
        int residueNumber = 0;
        int amino = -1;
        if (validResidue(p, a.residueIndex)) {
            const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
            residueNumber = r.address.residueNumber;
            amino = static_cast<int>(r.aminoAcid);
        }
        const int64_t rowId = static_cast<int64_t>(rowsWritten_);
        *rowsOut_ << rowId << ',' << a.atomIndex << ',' << static_cast<qint64>(row) << ','
                  << static_cast<qint64>(frameSlot) << ',' << static_cast<qint64>(orig) << ','
                  << num(body.run.timePs(row)) << ','
                  << static_cast<int>(a.element) << ',' << static_cast<int>(a.ffAtomType) << ','
                  << p.atomLabel(atom, model::NamingConvention::Bmrb) << ','
                  << p.atomLabel(atom, model::NamingConvention::Iupac) << ','
                  << a.residueIndex << ',' << residueNumber << ',' << amino << ','
                  << static_cast<int>(a.backboneRole) << ',' << static_cast<int>(a.locant) << ','
                  << static_cast<int>(a.branch.outer) << ',' << static_cast<int>(a.branch.inner) << ','
                  << static_cast<int>(a.diIndex) << ','
                  << static_cast<int>(a.ringPositionPrimary) << ','
                  << static_cast<int>(a.ringPositionSecondary) << ','
                  << static_cast<int>(a.planarGroup) << ','
                  << static_cast<int>(a.polarH) << ','
                  << static_cast<int>(a.prochiral) << ','
                  << static_cast<int>(a.equivalenceClass) << ','
                  << (a.aromatic ? 1 : 0) << ',' << static_cast<int>(a.formalCharge) << ','
                  << (a.isExchangeable ? 1 : 0) << ','
                  << static_cast<int>(role) << ',' << roleName(role) << ','
                  << static_cast<int>(stratum) << ',' << stratumName(stratum) << ','
                  << (direct.dft_present ? 1 : 0) << ','
                  << (agg.ring_present ? 1 : 0) << ',' << agg.ring_n << ',' << agg.ring_valid_n << ','
                  << (agg.charge_present ? 1 : 0) << ',' << agg.charge_n << ','
                  << agg.charge_excluded_same_residue_n << ','
                  << (agg.mc_lit_valid_present ? 1 : 0) << ',' << agg.bond_n << ','
                  << agg.bond_n_valid << ',' << (agg.ff14sb_field_present ? 1 : 0) << ','
                  << (direct.apbs_E_present ? 1 : 0) << ',' << (direct.apbs_efg_present ? 1 : 0) << ','
                  << (direct.mopac_coulomb_present ? 1 : 0) << ','
                  << (direct.mopac_mc_present ? 1 : 0) << ','
                  << (direct.aimnet2_charge_present ? 1 : 0) << ','
                  << (direct.aimnet2_crg_present ? 1 : 0) << ','
                  << (direct.aimnet2_embedding_present ? 1 : 0) << ','
                  << agg.ring_self_or_bonded_n << ',' << agg.bond_self_or_bonded_n << ','
                  << ((agg.ring_self_or_bonded_n + agg.bond_self_or_bonded_n
                       + agg.charge_excluded_same_residue_n) > 0 ? 1 : 0)
                  << ',' << csvScalar(charges.ff14sb)
                  << ',' << (charges.ff14sb.present ? 1 : 0)
                  << ',' << csvScalar(charges.mopac_welford_mean)
                  << ',' << (charges.mopac_welford_mean.present ? 1 : 0)
                  << ',' << (direct.hbond_count_present ? 1 : 0)
                  << ',' << ((direct.hbond_nearest_direction_present
                               && direct.hbond_flags_present) ? 1 : 0)
                  << ',' << (direct.hm_shielding_present ? 1 : 0)
                  << ',' << ((direct.water_efield_present
                               && direct.water_n_first_present
                               && direct.water_n_second_present) ? 1 : 0)
                  << ',' << ((direct.hydration_half_shell_present
                               && direct.hydration_dipole_cos_present) ? 1 : 0)
                  << ',' << (direct.sasa_present ? 1 : 0)
                  << ',' << (direct.sasa_normal_present ? 1 : 0)
                  << ',' << csvScalar(charges.eeq_charge)
                  << ',' << (charges.eeq_charge.present ? 1 : 0)
                  << ',' << csvScalar(charges.eeq_coordination_number)
                  << ',' << (charges.eeq_coordination_number.present ? 1 : 0)
                  << '\n';

        const std::array<double, 5> targetT2 =
            direct.dft_present ? direct.target.total_decomp.T2 : nanT2();
        const std::array<double, 3> targetT1 =
            direct.dft_present ? direct.target.total_decomp.T1 : nanT1();
        const std::array<double, 3> targetDiaT1 =
            direct.dft_present ? direct.target.dia_decomp.T1 : nanT1();
        const std::array<double, 5> targetDiaT2 =
            direct.dft_present ? direct.target.dia_decomp.T2 : nanT2();
        const std::array<double, 3> targetParaT1 =
            direct.dft_present ? direct.target.para_decomp.T1 : nanT1();
        const std::array<double, 5> targetParaT2 =
            direct.dft_present ? direct.target.para_decomp.T2 : nanT2();
        if (!targetT2_.writeArray(targetT2)) return false;
        if (!targetT0_.writeScalar<double>(direct.dft_present ? direct.target.total_decomp.T0 : kNaN))
            return false;
        if (!targetT1_.writeArray(targetT1)) return false;
        if (!targetDiaT0_.writeScalar<double>(direct.dft_present ? direct.target.dia_decomp.T0 : kNaN))
            return false;
        if (!targetDiaT1_.writeArray(targetDiaT1)) return false;
        if (!targetDiaT2_.writeArray(targetDiaT2)) return false;
        if (!targetParaT0_.writeScalar<double>(direct.dft_present ? direct.target.para_decomp.T0 : kNaN))
            return false;
        if (!targetParaT1_.writeArray(targetParaT1)) return false;
        if (!targetParaT2_.writeArray(targetParaT2)) return false;
        if (!classical_.writeArray(classical)) return false;
        if (!ringPaths_.writeArray(direct.ring_paths)) return false;
        if (!methodPaths_.writeArray(direct.method_paths)) return false;
        if (!hbondConditioning_.writeArray(direct.hbond_conditioning)) return false;
        if (!conditioning_.writeArray(conditioning)) return false;
        if (!dominance_.writeArray(dominance)) return false;
        if (!backboneAudit_.writeArray(agg.backbone_audit)) return false;
        if (cfg_.emit_embedding) {
            const bool embOk = direct.aimnet2_embedding_present && direct.aimnet2_embedding
                               && direct.aimnet2_embedding_dims == cfg_.embedding_dims;
            for (std::size_t i = 0; i < cfg_.embedding_dims; ++i) {
                const float v = embOk ? direct.aimnet2_embedding[i] : kNaNF32;
                if (!embedding_.writeScalar<float>(v)) return false;
            }
        }
        ++rowsWritten_;
        return true;
    }

    bool commit(const std::vector<WelfordCell>& modulation) {
        if (!ok_) return false;
        if (rowsOut_) rowsOut_->flush();
        if (ringIdentityOut_) ringIdentityOut_->flush();
        std::vector<double> mod;
        mod.reserve(atoms_ * kPerAtomDriverMagnitudeCols);
        for (std::size_t atom = 0; atom < atoms_; ++atom) {
            for (std::size_t c = 0; c < kPerAtomDriverMagnitudeCols; ++c)
                mod.push_back(modulation[atom * kPerAtomDriverMagnitudeCols + c].sd());
        }
        const bool modOk =
            writeNpy<double>(path(QStringLiteral("per_atom_substrate_driver_modulation_by_atom.npy")),
                             {atoms_, kPerAtomDriverMagnitudeCols}, mod, QByteArray("<f8"));
        return rowsFile_->commit() && ringIdentityFile_->commit()
               && targetT2_.commit() && targetT0_.commit() && targetT1_.commit()
               && targetDiaT0_.commit() && targetDiaT1_.commit() && targetDiaT2_.commit()
               && targetParaT0_.commit() && targetParaT1_.commit() && targetParaT2_.commit()
               && classical_.commit() && ringPaths_.commit() && methodPaths_.commit()
               && hbondConditioning_.commit()
               && conditioning_.commit() && dominance_.commit() && backboneAudit_.commit()
               && (!cfg_.emit_embedding || embedding_.commit()) && modOk;
    }

    std::size_t rowsWritten() const { return rowsWritten_; }

private:
    QString path(const QString& name) const { return QStringLiteral("%1/%2").arg(outDir_, name); }

    QString rowHeader() const {
        return QStringLiteral(
            "row_id,atom_index,h5_row,frame_slot,original_frame_index,time_ps,"
            "element_ord,ff_atom_type_ord,atom_name,iupac_atom_name,residue_index,residue_number,"
            "amino_acid_ord,backbone_role_ord,locant_ord,branch_outer_ord,branch_inner_ord,"
            "di_index_ord,ring_position_primary_ord,ring_position_secondary_ord,planar_group_ord,"
            "polar_h_kind_ord,prochiral_ord,equivalence_class,aromatic,formal_charge,is_exchangeable,"
            "role_ord,role,stratum_ord,stratum,dft_present,ring_present,ring_n,ring_valid_n,"
            "charge_present,charge_n,charge_excluded_same_residue_n,mc_lit_valid_present,"
            "bond_n,bond_n_valid,ff14sb_field_present,apbs_E_present,apbs_efg_present,"
            "mopac_coulomb_shielding_present,mopac_mc_shielding_present,aimnet2_charge_present,"
            "aimnet2_crg_present,aimnet2_embedding_present,ring_self_or_bonded_n,"
            "bond_self_or_bonded_n,has_self_or_bonded_driver,"
            "ff14sb_charge,ff14sb_charge_present,"
            "mopac_welford_mean_charge,mopac_welford_mean_charge_present,"
            "hbond_count_present,hbond_geometry_present,hm_shielding_present,"
            "water_field_present,hydration_shell_present,"
            "sasa_present,sasa_normal_present,eeq_charge,eeq_charge_present,"
            "eeq_coordination_number,eeq_coordination_number_present");
    }

    QString outDir_;
    std::size_t rows_ = 0;
    std::size_t atoms_ = 0;
    PerAtomSubstrateConfig cfg_;
    bool ok_ = false;
    std::size_t rowsWritten_ = 0;
    std::unique_ptr<QSaveFile> rowsFile_;
    std::unique_ptr<QSaveFile> ringIdentityFile_;
    std::unique_ptr<QTextStream> rowsOut_;
    std::unique_ptr<QTextStream> ringIdentityOut_;
    StreamingNpy targetT2_;
    StreamingNpy targetT0_;
    StreamingNpy targetT1_;
    StreamingNpy targetDiaT0_;
    StreamingNpy targetDiaT1_;
    StreamingNpy targetDiaT2_;
    StreamingNpy targetParaT0_;
    StreamingNpy targetParaT1_;
    StreamingNpy targetParaT2_;
    StreamingNpy classical_;
    StreamingNpy ringPaths_;
    StreamingNpy methodPaths_;
    StreamingNpy hbondConditioning_;
    StreamingNpy conditioning_;
    StreamingNpy dominance_;
    StreamingNpy backboneAudit_;
    StreamingNpy embedding_;
};

void addColumnSpec(QJsonArray& arr, const QString& array, const QString& name,
                   int index, const QString& units, const QString& irreps,
                   const QString& mechanism,
                   const QString& sign = QString(), const QString& nativeAxis = QStringLiteral("rediscover_target_row")) {
    QJsonObject o;
    o.insert(QStringLiteral("array"), array);
    o.insert(QStringLiteral("column"), name);
    o.insert(QStringLiteral("index"), index);
    o.insert(QStringLiteral("units"), units);
    o.insert(QStringLiteral("irreps"), irreps);
    o.insert(QStringLiteral("mechanism"), mechanism);
    o.insert(QStringLiteral("sign_convention"), sign);
    o.insert(QStringLiteral("native_axis"), nativeAxis);
    arr.append(o);
}

bool writeColumnSpecs(const QString& outDir, const PerAtomSubstrateConfig& cfg) {
    QJsonObject root;
    root.insert(QStringLiteral("relationship"), QStringLiteral("per_atom_substrate"));
    root.insert(QStringLiteral("feature_dtype"), QStringLiteral("float64"));
    root.insert(QStringLiteral("t2_order"), QStringLiteral("[xy,yz,zz,xz,xx-yy]"));
    QJsonArray cols;
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_target_T0"),
                  QStringLiteral("target_T0"), 0, QStringLiteral("ppm"), QStringLiteral("1x0e"),
                  QStringLiteral("quantum_reference"));
    for (int i = 0; i < 3; ++i) {
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_target_T1"),
                      QStringLiteral("target_T1_%1").arg(i), i, QStringLiteral("ppm"),
                      QStringLiteral("1x1e"), QStringLiteral("quantum_reference_diagnostic"));
    }
    for (int i = 0; i < 5; ++i) {
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_target_T2"),
                      QStringLiteral("target_T2_%1").arg(i), i, QStringLiteral("ppm"),
                      QStringLiteral("1x2e"), QStringLiteral("quantum_reference"));
    }
    const QStringList targetSplits = {
        QStringLiteral("dia"),
        QStringLiteral("para"),
    };
    for (const QString& split : targetSplits) {
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_target_%1_T0").arg(split),
                      QStringLiteral("target_%1_T0").arg(split), 0, QStringLiteral("ppm"),
                      QStringLiteral("1x0e"), QStringLiteral("quantum_reference_split"));
        for (int i = 0; i < 3; ++i) {
            addColumnSpec(cols, QStringLiteral("per_atom_substrate_target_%1_T1").arg(split),
                          QStringLiteral("target_%1_T1_%2").arg(split).arg(i), i,
                          QStringLiteral("ppm"), QStringLiteral("1x1e"),
                          QStringLiteral("quantum_reference_split_diagnostic"));
        }
        for (int i = 0; i < 5; ++i) {
            addColumnSpec(cols, QStringLiteral("per_atom_substrate_target_%1_T2").arg(split),
                          QStringLiteral("target_%1_T2_%2").arg(split).arg(i), i,
                          QStringLiteral("ppm"), QStringLiteral("1x2e"),
                          QStringLiteral("quantum_reference_split"));
        }
    }
    for (int i = 0; i < kClassicalColumns.size(); ++i) {
        const QString name = kClassicalColumns[i];
        QString units;
        QString irreps = QStringLiteral("0e");
        QString mechanism = QStringLiteral("geometry");
        QString sign;
        if (name.contains(QStringLiteral("ring_jb"))) {
            units = QStringLiteral("ppm");
            irreps = name.contains(QStringLiteral("T2")) ? QStringLiteral("1x2e") : QStringLiteral("1x0e");
            mechanism = QStringLiteral("ring_current");
            sign = QStringLiteral("sigma_ab=-dB_sec_a/dB0_b");
        } else if (name.contains(QStringLiteral("hbond_nearest_dir"))) {
            irreps = QStringLiteral("1o");
            mechanism = QStringLiteral("hbond");
        } else if (name.contains(QStringLiteral("hbond_count"))
                   || name.contains(QStringLiteral("hbond_flags_"))) {
            units = name.contains(QStringLiteral("count")) ? QStringLiteral("count") : QString();
            irreps = QStringLiteral("0e");
            mechanism = QStringLiteral("hbond");
        } else if (name.contains(QStringLiteral("hm_shielding"))) {
            units = QStringLiteral("Angstrom^-1");
            irreps = QStringLiteral("1x2e");
            mechanism = QStringLiteral("ring_current_alt");
            sign = QStringLiteral("sigma_ab=-dB_sec_a/dB0_b");
        } else if (name.contains(QStringLiteral("water_efield"))) {
            units = QStringLiteral("V/A");
            irreps = name.endsWith(QStringLiteral("_mag")) ? QStringLiteral("0e") : QStringLiteral("1o");
            mechanism = QStringLiteral("water_field");
        } else if (name.contains(QStringLiteral("water_n_"))) {
            units = QStringLiteral("count");
            irreps = QStringLiteral("0e");
            mechanism = QStringLiteral("water_field");
        } else if (name.contains(QStringLiteral("water_half_shell"))
                   || name.contains(QStringLiteral("water_dipole_cos"))) {
            irreps = QStringLiteral("0e");
            mechanism = QStringLiteral("water_field");
        } else if (name == QStringLiteral("sasa")) {
            units = QStringLiteral("A^2");
            irreps = QStringLiteral("0e");
            mechanism = QStringLiteral("sasa");
        } else if (name.contains(QStringLiteral("sasa_normal"))) {
            irreps = QStringLiteral("1o");
            mechanism = QStringLiteral("sasa");
        } else if (name.contains(QStringLiteral("charge_q_over_r3")) || name.contains(QStringLiteral("apbs_efg"))) {
            units = name.contains(QStringLiteral("apbs")) ? QStringLiteral("V/A^2") : QStringLiteral("CoulombKe*e/A^3");
            irreps = QStringLiteral("1x2e");
            mechanism = QStringLiteral("electrostatic_efg");
        } else if (name.contains(QStringLiteral("mc_lit")) || name.contains(QStringLiteral("mopac_mc"))) {
            units = QStringLiteral("ppm");
            irreps = name.contains(QStringLiteral("T0")) ? QStringLiteral("1x0e") : QStringLiteral("1x2e");
            mechanism = QStringLiteral("bond_anisotropy");
            sign = QStringLiteral("sigma_ab=-dB_sec_a/dB0_b");
        } else if (name.contains(QStringLiteral("mopac_coulomb"))) {
            units = QStringLiteral("ppm");
            irreps = QStringLiteral("1x2e");
            mechanism = QStringLiteral("electrostatic_efg");
            sign = QStringLiteral("sigma_ab=-dB_sec_a/dB0_b");
        } else if (name.contains(QStringLiteral("field")) || name.contains(QStringLiteral("apbs_E"))) {
            units = name.contains(QStringLiteral("apbs")) ? QStringLiteral("V/A") : QStringLiteral("e/A^2");
            irreps = name.endsWith(QStringLiteral("_mag")) ? QStringLiteral("0e") : QStringLiteral("1o");
            mechanism = name.contains(QStringLiteral("ff14sb")) ? QStringLiteral("charges")
                                                                 : QStringLiteral("electrostatic_efg");
        } else if (name.contains(QStringLiteral("aimnet2"))) {
            units = name.contains(QStringLiteral("charge")) ? QStringLiteral("e") : QStringLiteral("e^2/A");
            irreps = name.endsWith(QStringLiteral("_x")) || name.endsWith(QStringLiteral("_y"))
                         || name.endsWith(QStringLiteral("_z")) ? QStringLiteral("1o") : QStringLiteral("0e");
            mechanism = QStringLiteral("aimnet2");
        }
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_features_classical"), name, i,
                      units, irreps, mechanism, sign);
    }
    for (int i = 0; i < kRingPathColumns.size(); ++i) {
        const QString name = kRingPathColumns[i];
        QString units;
        QString irreps = QStringLiteral("0e");
        QString mechanism = QStringLiteral("ring_current");
        QString sign;
        if (name.contains(QStringLiteral("_T1_")) || name.contains(QStringLiteral("total_B"))) {
            irreps = QStringLiteral("1x1e");
        } else if (name.contains(QStringLiteral("_T2_"))) {
            irreps = QStringLiteral("1x2e");
        }
        if (name.contains(QStringLiteral("bs_"))) {
            units = name.contains(QStringLiteral("total_B")) ? QStringLiteral("T")
                                                             : QStringLiteral("ppm_T_per_nA");
        } else if (name.contains(QStringLiteral("hm_"))) {
            units = QStringLiteral("Angstrom^-1");
        } else if (name.contains(QStringLiteral("ringchi_"))) {
            units = QStringLiteral("Angstrom^-3");
        } else if (name.contains(QStringLiteral("pq_"))) {
            units = name.contains(QStringLiteral("_T0_")) ? QStringLiteral("Angstrom^-4")
                                                          : QStringLiteral("Angstrom^-5");
            mechanism = QStringLiteral("ring_efg");
        } else if (name.contains(QStringLiteral("disp_"))) {
            units = QStringLiteral("Angstrom^-6");
            mechanism = QStringLiteral("ring_dispersion");
        }
        if (name.contains(QStringLiteral("shielding")) || name.contains(QStringLiteral("per_type_T")))
            sign = QStringLiteral("sigma_ab=-dB_sec_a/dB0_b");
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_features_ring_paths"),
                      name, i, units, irreps, mechanism, sign);
    }
    for (int i = 0; i < kMethodPathColumns.size(); ++i) {
        const QString name = kMethodPathColumns[i];
        QString units;
        QString irreps = QStringLiteral("0e");
        QString mechanism = QStringLiteral("comparison_method");
        QString sign;
        if (name.contains(QStringLiteral("_T2_")) || name.contains(QStringLiteral("efg"))) {
            irreps = QStringLiteral("1x2e");
        } else if (name.endsWith(QStringLiteral("_x")) || name.endsWith(QStringLiteral("_y"))
                   || name.endsWith(QStringLiteral("_z"))) {
            irreps = QStringLiteral("1o");
        }
        if (name.contains(QStringLiteral("mc_"))) {
            units = QStringLiteral("Angstrom^-3");
            mechanism = QStringLiteral("bond_anisotropy");
            sign = QStringLiteral("sigma_ab=-dB_sec_a/dB0_b");
        } else if (name.contains(QStringLiteral("bond_orders"))) {
            mechanism = QStringLiteral("mopac_bond_order");
        } else if (name.contains(QStringLiteral("mopac_coulomb_E"))) {
            units = QStringLiteral("V/A");
            mechanism = QStringLiteral("electrostatic_efg");
        } else if (name.contains(QStringLiteral("mopac_coulomb_efg"))
                   || name.contains(QStringLiteral("aimnet2_efg"))
                   || name.contains(QStringLiteral("water_efg"))) {
            units = QStringLiteral("V/A^2");
            mechanism = name.contains(QStringLiteral("aimnet2")) ? QStringLiteral("aimnet2")
                       : name.contains(QStringLiteral("water")) ? QStringLiteral("water_field")
                                                                 : QStringLiteral("electrostatic_efg");
        } else if (name.contains(QStringLiteral("water_efield_first"))) {
            units = QStringLiteral("V/A");
            mechanism = QStringLiteral("water_field");
        } else if (name.contains(QStringLiteral("eeq_cn"))) {
            mechanism = QStringLiteral("eeq");
        } else if (name.contains(QStringLiteral("mopac_scalars"))) {
            mechanism = QStringLiteral("charges");
        }
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_features_method_paths"),
                      name, i, units, irreps, mechanism, sign);
    }
    for (int i = 0; i < kHbondConditioningColumns.size(); ++i) {
        const QString name = kHbondConditioningColumns[i];
        QString units;
        QString irreps = QStringLiteral("0e");
        QString mechanism = QStringLiteral("conditioning");
        QString sign;
        if (name.contains(QStringLiteral("larsen_hbond")) && name.contains(QStringLiteral("_T2_"))) {
            units = QStringLiteral("ppm");
            irreps = QStringLiteral("1x2e");
            mechanism = QStringLiteral("larsen_hbond");
            sign = QStringLiteral("sigma_ab=-dB_sec_a/dB0_b");
        } else if (name == QStringLiteral("larsen_hbond_water_term")) {
            units = QStringLiteral("ppm");
            mechanism = QStringLiteral("larsen_hbond");
        } else if (name.contains(QStringLiteral("hbond_scalars"))) {
            mechanism = QStringLiteral("hbond_geometry");
        } else if (name.contains(QStringLiteral("dssp_hbond"))
                   || name.contains(QStringLiteral("dssp_chem"))
                   || name.contains(QStringLiteral("dssp_acceptor"))
                   || name.contains(QStringLiteral("dssp_donor"))) {
            mechanism = QStringLiteral("dssp_hbond_backup");
            units = name.contains(QStringLiteral("energy")) ? QStringLiteral("kcal/mol") : QString();
        } else if (name.contains(QStringLiteral("dssp_ss8")) || name.contains(QStringLiteral("dssp_chi"))) {
            mechanism = QStringLiteral("secondary_structure");
        } else if (name == QStringLiteral("omega_actual")) {
            units = QStringLiteral("radians");
            mechanism = QStringLiteral("planar_geometry");
        } else if (name == QStringLiteral("pyramidalization")) {
            units = QStringLiteral("radians");
            mechanism = QStringLiteral("planar_geometry");
        } else if (name.contains(QStringLiteral("ring_geometry"))) {
            units = QStringLiteral("A");
            mechanism = QStringLiteral("ring_geometry");
        }
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_features_hbond_conditioning"),
                      name, i, units, irreps, mechanism, sign);
    }
    for (int i = 0; i < kConditioningColumns.size(); ++i) {
        const QString name = kConditioningColumns[i];
        const QString units = name.endsWith(QStringLiteral("_r")) ? QStringLiteral("A") : QString();
        const QString mechanism = name.contains(QStringLiteral("dominant_fraction"))
                                      || name.contains(QStringLiteral("gap_to_2nd"))
                                  ? QStringLiteral("isolation")
                                  : QStringLiteral("conditioning");
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_features_conditioning"),
                      name, i, units, QStringLiteral("0e"),
                      mechanism);
    }
    for (int i = 0; i < kDominanceColumns.size(); ++i) {
        const QString name = kDominanceColumns[i];
        const QString units = name.startsWith(QStringLiteral("gap_to_2nd")) ? QStringLiteral("A") : QString();
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_features_dominance"),
                      name, i, units, QStringLiteral("0e"),
                      QStringLiteral("isolation"));
    }
    for (int i = 0; i < kMagnitudeColumns.size(); ++i) {
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_driver_modulation_by_atom"),
                      kMagnitudeColumns[i], i, QString(), QStringLiteral("0e"),
                      QStringLiteral("conditioning"), QString(),
                      QStringLiteral("atom"));
    }
    for (int i = 0; i < kPartitionBinColumns.size(); ++i) {
        QString family = QStringLiteral("partition_bin");
        if (kPartitionBinColumns[i].contains(QStringLiteral("_4_6_8_10")))
            family = QStringLiteral("geometry_distance_band");
        else if (kPartitionBinColumns[i].contains(QStringLiteral("_sd_")))
            family = QStringLiteral("driver_modulation_quintile");
        else
            family = QStringLiteral("driver_magnitude_quintile");
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_partition_bins"),
                      kPartitionBinColumns[i], i, QStringLiteral("bin_id"),
                      QStringLiteral("0e"), family);
    }
    for (int i = 0; i < kDominanceBinColumns.size(); ++i) {
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_dominance_bins"),
                      kDominanceBinColumns[i], i, QStringLiteral("bin_id"),
                      QStringLiteral("0e"), QStringLiteral("dominance_quintile"));
    }
    for (int i = 0; i < kBackboneAuditColumns.size(); ++i) {
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_backbone_audit"),
                      kBackboneAuditColumns[i], i, QString(), QStringLiteral("audit"),
                      QStringLiteral("provenance_qc"));
    }
    if (cfg.emit_embedding) {
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_aimnet2_embedding"),
                      QStringLiteral("embedding_000..255"), 0, QString(), QStringLiteral("256x0e"),
                      QStringLiteral("aimnet2"));
    }
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                  QStringLiteral("ff14sb_charge"), 52, QStringLiteral("e"),
                  QStringLiteral("0e"), QStringLiteral("charges"));
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                  QStringLiteral("ff14sb_charge_present"), 53, QString(),
                  QStringLiteral("0e"), QStringLiteral("provenance_qc"));
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                  QStringLiteral("mopac_welford_mean_charge"), 54, QStringLiteral("e"),
                  QStringLiteral("0e"), QStringLiteral("charges"));
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                  QStringLiteral("mopac_welford_mean_charge_present"), 55, QString(),
                  QStringLiteral("0e"), QStringLiteral("provenance_qc"));
    const QStringList rowPresentCols = {
        QStringLiteral("hbond_count_present"),
        QStringLiteral("hbond_geometry_present"),
        QStringLiteral("hm_shielding_present"),
        QStringLiteral("water_field_present"),
        QStringLiteral("hydration_shell_present"),
        QStringLiteral("sasa_present"),
        QStringLiteral("sasa_normal_present"),
    };
    for (int i = 0; i < rowPresentCols.size(); ++i) {
        addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                      rowPresentCols[i], 56 + i, QString(), QStringLiteral("0e"),
                      QStringLiteral("provenance_qc"));
    }
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                  QStringLiteral("eeq_charge"), 67, QStringLiteral("e"),
                  QStringLiteral("0e"), QStringLiteral("eeq"));
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                  QStringLiteral("eeq_charge_present"), 68, QString(),
                  QStringLiteral("0e"), QStringLiteral("provenance_qc"));
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                  QStringLiteral("eeq_coordination_number"), 69, QString(),
                  QStringLiteral("0e"), QStringLiteral("eeq"));
    addColumnSpec(cols, QStringLiteral("per_atom_substrate_rows"),
                  QStringLiteral("eeq_coordination_number_present"), 70, QString(),
                  QStringLiteral("0e"), QStringLiteral("provenance_qc"));
    root.insert(QStringLiteral("columns"), cols);
    QSaveFile f(QStringLiteral("%1/per_atom_substrate_column_specs.json").arg(outDir));
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    f.write(QJsonDocument(root).toJson(QJsonDocument::Indented));
    return f.commit();
}

QJsonObject labelVocab() {
    QJsonObject root;
    QJsonObject roles;
    for (int i = 0; i <= static_cast<int>(RoleOrd::Other); ++i)
        roles.insert(QString::number(i), roleName(static_cast<RoleOrd>(i)));
    QJsonObject strata;
    for (int i = 0; i <= static_cast<int>(StratumOrd::Other); ++i)
        strata.insert(QString::number(i), stratumName(static_cast<StratumOrd>(i)));
    root.insert(QStringLiteral("role"), roles);
    root.insert(QStringLiteral("stratum"), strata);
    root.insert(QStringLiteral("source"), QStringLiteral("C++ typed projection over QtAtom fields"));
    return root;
}

bool writePerAtomManifest(const QString& outDir, const PerAtomSubstrateStats& stats,
                          const PerAtomSubstrateConfig& cfg, const DftFrameAlignment& alignment) {
    QJsonObject root;
    root.insert(QStringLiteral("relationship"), QStringLiteral("per_atom_substrate"));
    root.insert(QStringLiteral("relationship_kind"), QStringLiteral("per_atom_aggregate"));
    root.insert(QStringLiteral("normalization"), QStringLiteral("raw_lab_frame"));
    root.insert(QStringLiteral("row_alignment"),
                QStringLiteral("row_id == frame_slot * n_atoms + atom_index; sidecar row i == rows.csv row_id i"));
    root.insert(QStringLiteral("n_atoms"), static_cast<qint64>(stats.atom_count));
    root.insert(QStringLiteral("n_dft_frames"), static_cast<qint64>(stats.dft_rows));
    root.insert(QStringLiteral("rows"), static_cast<qint64>(stats.rows));
    root.insert(QStringLiteral("embedding_separable"), cfg.emit_embedding);
    root.insert(QStringLiteral("label_vocabulary"), labelVocab());
    QJsonObject cutoffs;
    cutoffs.insert(QStringLiteral("ring_cutoff_A"), cfg.ring_cutoff_A);
    cutoffs.insert(QStringLiteral("bond_cutoff_A"), cfg.bond_cutoff_A);
    cutoffs.insert(QStringLiteral("charge_cutoff_A"), cfg.charge_cutoff_A);
    cutoffs.insert(QStringLiteral("mc_near_field_ratio"), cfg.mc_near_field_ratio);
    root.insert(QStringLiteral("cutoffs"), cutoffs);
    QJsonObject align;
    align.insert(QStringLiteral("n_frames"), alignment.n_frames);
    align.insert(QStringLiteral("max_angle_deg"), alignment.max_angle_deg);
    align.insert(QStringLiteral("mean_angle_deg"), alignment.mean_angle_deg);
    align.insert(QStringLiteral("max_rmsd_A"), alignment.max_rmsd_A);
    align.insert(QStringLiteral("mean_rmsd_A"), alignment.mean_rmsd_A);
    align.insert(QStringLiteral("t2_components"),
                 alignment.max_angle_deg < 1.0 ? QStringLiteral("FRAME-ALIGNED")
                                                : QStringLiteral("ROTATED"));
    align.insert(QStringLiteral("target_T1_frame"),
                 alignment.max_angle_deg < 1.0 ? QStringLiteral("frame_verified")
                                                : QStringLiteral("t1_frame_unverified"));
    root.insert(QStringLiteral("dft_frame_alignment"), align);
    QJsonObject sizeGate;
    constexpr std::size_t kPiece3BAppendFloat64Cols =
        kPerAtomTargetDecompositionCols + kPerAtomRingPathCols + kPerAtomMethodPathCols
        + kPerAtomHbondConditioningCols;
    const double appendGiB =
        static_cast<double>(stats.rows) * static_cast<double>(kPiece3BAppendFloat64Cols)
        * 8.0 / 1073741824.0;
    sizeGate.insert(QStringLiteral("status"), appendGiB <= 10.0 ? QStringLiteral("pass")
                                                                : QStringLiteral("fail"));
    sizeGate.insert(QStringLiteral("appended_float64_columns"),
                    static_cast<qint64>(kPiece3BAppendFloat64Cols));
    sizeGate.insert(QStringLiteral("appended_payload_GiB"), appendGiB);
    sizeGate.insert(QStringLiteral("axis_contract"), QStringLiteral("(atom,frame) only"));
    root.insert(QStringLiteral("piece3b_size_gate"), sizeGate);
    QJsonArray sidecars;
    for (const QString& s : PerAtomSubstrateSidecars(cfg)) sidecars.append(s);
    root.insert(QStringLiteral("sidecars"), sidecars);
    QJsonObject queries;
    queries.insert(QStringLiteral("pairs_rows"), static_cast<qint64>(stats.pair_query_rows));
    queries.insert(QStringLiteral("top_sources_rows"), static_cast<qint64>(stats.top_source_query_rows));
    queries.insert(QStringLiteral("dominance_fractions_rows"), static_cast<qint64>(stats.dominance_query_rows));
    queries.insert(QStringLiteral("dominance_fractions_build4_rows"),
                   static_cast<qint64>(stats.dominance_build4_query_rows));
    queries.insert(QStringLiteral("reader_pairs_rows"), static_cast<qint64>(stats.reader_pair_query_rows));
    root.insert(QStringLiteral("pair_index_named_queries"), queries);
    QJsonObject dominanceDefinitions;
    QJsonObject fieldDef;
    fieldDef.insert(QStringLiteral("mechanism"), QStringLiteral("field"));
    fieldDef.insert(QStringLiteral("source"), QStringLiteral("MOPAC-Coulomb"));
    fieldDef.insert(QStringLiteral("charge_array"), QStringLiteral("MopacChargeWelfordMean"));
    fieldDef.insert(QStringLiteral("charge_scope"), QStringLiteral("resident atom-axis mean MOPAC charge"));
    fieldDef.insert(QStringLiteral("contribution"), QStringLiteral("|q_mopac * rhat / r^2|"));
    fieldDef.insert(QStringLiteral("source_set"), QStringLiteral("charge_sites within charge_cutoff_A"));
    fieldDef.insert(QStringLiteral("exclusions"), QStringLiteral("self and same residue"));
    dominanceDefinitions.insert(QStringLiteral("field"), fieldDef);
    root.insert(QStringLiteral("dominance_definitions"), dominanceDefinitions);
    QJsonObject hunter;
    hunter.insert(QStringLiteral("anti_circular_assertion"), stats.hunter_anti_circular_assertion);
    QJsonObject hunterCounts;
    for (auto it = stats.hunter_candidate_counts.constBegin();
         it != stats.hunter_candidate_counts.constEnd(); ++it) {
        hunterCounts.insert(it.key(), static_cast<qint64>(it.value()));
    }
    hunter.insert(QStringLiteral("candidate_counts"), hunterCounts);
    hunter.insert(QStringLiteral("selection_predicate"),
                  QStringLiteral("typed_habitat && isolation && motion && quiet; DFT measured after selection only"));
    root.insert(QStringLiteral("case_hunter"), hunter);
    QJsonObject support;
    support.insert(QStringLiteral("ff14sb_charge_present_rows"),
                   static_cast<qint64>(stats.ff14sb_charge_present));
    support.insert(QStringLiteral("mopac_welford_mean_charge_present_rows"),
                   static_cast<qint64>(stats.mopac_welford_mean_charge_present));
    support.insert(QStringLiteral("charge_complete_rows"),
                   static_cast<qint64>(stats.charge_complete));
    support.insert(QStringLiteral("hbond_count_present_rows"),
                   static_cast<qint64>(stats.hbond_count_present));
    support.insert(QStringLiteral("hbond_geometry_present_rows"),
                   static_cast<qint64>(stats.hbond_geometry_present));
    support.insert(QStringLiteral("hm_shielding_present_rows"),
                   static_cast<qint64>(stats.hm_shielding_present));
    support.insert(QStringLiteral("water_field_present_rows"),
                   static_cast<qint64>(stats.water_field_present));
    support.insert(QStringLiteral("hydration_shell_present_rows"),
                   static_cast<qint64>(stats.hydration_shell_present));
    support.insert(QStringLiteral("sasa_present_rows"),
                   static_cast<qint64>(stats.sasa_present));
    support.insert(QStringLiteral("sasa_normal_present_rows"),
                   static_cast<qint64>(stats.sasa_normal_present));
    support.insert(QStringLiteral("eeq_charge_present_rows"),
                   static_cast<qint64>(stats.eeq_charge_present));
    support.insert(QStringLiteral("eeq_coordination_number_present_rows"),
                   static_cast<qint64>(stats.eeq_coordination_number_present));
    root.insert(QStringLiteral("feature_support_rows"), support);
    QJsonObject audit;
    for (auto it = stats.new_channel_audit.constBegin(); it != stats.new_channel_audit.constEnd(); ++it) {
        QJsonObject o;
        o.insert(QStringLiteral("present_rows"), static_cast<qint64>(it.value().present));
        if (it.value().has_range) {
            o.insert(QStringLiteral("min"), it.value().min);
            o.insert(QStringLiteral("max"), it.value().max);
        } else {
            o.insert(QStringLiteral("min"), QJsonValue::Null);
            o.insert(QStringLiteral("max"), QJsonValue::Null);
        }
        audit.insert(it.key(), o);
    }
    root.insert(QStringLiteral("new_channel_audit"), audit);
    root.insert(QStringLiteral("partition_bins"), stats.partition_bin_manifest);
    root.insert(QStringLiteral("dominance_bins"), stats.dominance_bin_manifest);
    QJsonArray absent;
    for (const QString& s : stats.absent_new_channel_slabs) absent.append(s);
    root.insert(QStringLiteral("absent_new_channel_slabs"), absent);
    QSaveFile f(QStringLiteral("%1/per_atom_substrate_manifest.json").arg(outDir));
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    f.write(QJsonDocument(root).toJson(QJsonDocument::Indented));
    return f.commit();
}

bool writeQueryManifest(const QString& outDir, const QString& name, const QString& mechanism,
                        std::size_t rows, const QJsonArray& columns) {
    QDir().mkpath(QStringLiteral("%1/query_results").arg(outDir));
    QJsonObject root;
    root.insert(QStringLiteral("query_name"), name);
    root.insert(QStringLiteral("relationship_or_mechanism"), mechanism);
    root.insert(QStringLiteral("source"), QStringLiteral("pair_index_query"));
    root.insert(QStringLiteral("delivery"), QStringLiteral("transient_named_query_output"));
    root.insert(QStringLiteral("row_count"), static_cast<qint64>(rows));
    root.insert(QStringLiteral("lazy_columns_computed"), columns);
    root.insert(QStringLiteral("cutoff_policy"), QStringLiteral("aggregate_reducer_cutoffs"));
    QSaveFile f(QStringLiteral("%1/query_results/%2_manifest.json").arg(outDir, name));
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    f.write(QJsonDocument(root).toJson(QJsonDocument::Indented));
    return f.commit();
}

void writePairRows(QTextStream& out, std::size_t targetRowId, std::size_t atom, std::size_t frameSlot,
                   std::size_t h5Row, std::size_t orig, const PairContribution& p) {
    out << static_cast<qint64>(targetRowId) << ',' << static_cast<qint64>(atom) << ','
        << static_cast<qint64>(frameSlot) << ',' << static_cast<qint64>(h5Row) << ','
        << static_cast<qint64>(orig) << ',' << p.mechanism << ',' << p.source_kind << ','
        << p.source_id << ',' << p.source_atom_index << ',' << p.source_cloud_index << ','
        << p.source_category_ord << ',' << p.pointer_flags << '\n';
}

bool isLegacyPairMechanism(const PairContribution& p) {
    return p.mechanism == QStringLiteral("ring_jb")
           || p.mechanism == QStringLiteral("charge_q_over_r3")
           || p.mechanism == QStringLiteral("mc_lit_valid");
}

bool emitPairQueries(const Body& body, const QString& outDir, const PerAtomSubstrateConfig& cfg,
                     PerAtomSubstrateStats& stats) {
    if (!body.run.protein) return true;
    const std::size_t nAtoms = body.run.protein->atomCount();
    const std::vector<std::size_t>& rows = body.run.frameMap.dftRows();
    QDir().mkpath(QStringLiteral("%1/query_results").arg(outDir));

    auto openText = [&](const QString& name, const QString& header, std::unique_ptr<QSaveFile>& file,
                        std::unique_ptr<QTextStream>& out) {
        file = std::make_unique<QSaveFile>(QStringLiteral("%1/query_results/%2_rows.csv").arg(outDir, name));
        if (!file->open(QIODevice::WriteOnly | QIODevice::Text)) return false;
        out = std::make_unique<QTextStream>(file.get());
        *out << header << '\n';
        return true;
    };

    const QString pairHeader =
        QStringLiteral("target_row_id,target_atom_index,frame_slot,h5_row,original_frame_index,"
                       "mechanism,source_kind,source_id,source_atom_index,source_cloud_index,"
                       "source_category_ord,pointer_flags");
    std::unique_ptr<QSaveFile> pairsFile;
    std::unique_ptr<QTextStream> pairsOut;
    if (!openText(QStringLiteral("pairs"), pairHeader, pairsFile, pairsOut)) return false;
    std::vector<double> geom;
    std::vector<double> kt2;
    std::vector<double> kt0;
    std::vector<double> contrib;

    std::unique_ptr<QSaveFile> topFile;
    std::unique_ptr<QTextStream> topOut;
    if (!openText(QStringLiteral("top_sources"),
                  QStringLiteral("target_row_id,target_atom_index,frame_slot,mechanism,rank,source_kind,source_id,contribution"),
                  topFile, topOut)) return false;

    std::unique_ptr<QSaveFile> domFile;
    std::unique_ptr<QTextStream> domOut;
    if (!openText(QStringLiteral("dominance_fractions"),
                  QStringLiteral("target_row_id,target_atom_index,frame_slot,ring_abs_T2_dom_frac,charge_abs_T2_dom_frac,mc_abs_T2_dom_frac,mopac_abs_T2_dom_frac,aimnet2_embedding_norm_rank_frac"),
                  domFile, domOut)) return false;

    std::unique_ptr<QSaveFile> domBuild4File;
    std::unique_ptr<QTextStream> domBuild4Out;
    if (!openText(QStringLiteral("dominance_fractions_build4"),
                  QStringLiteral("target_row_id,target_atom_index,frame_slot,field_mopac_coulomb_dom_frac"),
                  domBuild4File, domBuild4Out)) return false;

    std::unique_ptr<QSaveFile> readerFile;
    std::unique_ptr<QTextStream> readerOut;
    if (!openText(QStringLiteral("reader_pairs"), pairHeader, readerFile, readerOut)) return false;

    auto appendPairPayload = [&](const PairContribution& p) {
        geom.insert(geom.end(), {p.disp.x(), p.disp.y(), p.disp.z(), p.r, p.inv_r3, p.cos_theta, p.dipolar});
        for (double v : p.kernel_T2) kt2.push_back(v);
        kt0.push_back(p.kernel_T0);
        contrib.push_back(p.contribution);
    };

    const std::size_t frameSlots = std::min(cfg.query_frame_slots, rows.size());
    for (std::size_t fs = 0; fs < frameSlots; ++fs) {
        const std::size_t h5Row = rows[fs];
        const std::size_t orig = body.run.frameMap.originalIndex(h5Row);
        for (std::size_t atom = 0; atom < nAtoms; ++atom) {
            const std::size_t targetRowId = fs * nAtoms + atom;
            const std::vector<PairContribution> pairs = PerAtomRowPairContributions(body, atom, h5Row, cfg, LocalFrame{});
            std::vector<double> ringVals, chargeVals, mcVals, fieldVals;
            std::vector<PairContribution> legacyPairs;
            legacyPairs.reserve(pairs.size());
            for (const PairContribution& p : pairs) {
                if (p.mechanism == QStringLiteral("field_mopac_coulomb")) fieldVals.push_back(p.contribution);
                if (isLegacyPairMechanism(p)) legacyPairs.push_back(p);
            }
            std::vector<PairContribution> top = legacyPairs;
            std::sort(top.begin(), top.end(), [](const PairContribution& a, const PairContribution& b) {
                const double av = std::isfinite(a.contribution) ? a.contribution : -1.0;
                const double bv = std::isfinite(b.contribution) ? b.contribution : -1.0;
                return av > bv;
            });
            int rank = 0;
            for (const PairContribution& p : legacyPairs) {
                writePairRows(*pairsOut, targetRowId, atom, fs, h5Row, orig, p);
                appendPairPayload(p);
                ++stats.pair_query_rows;
                if (p.mechanism == QStringLiteral("ring_jb")) ringVals.push_back(p.contribution);
                if (p.mechanism == QStringLiteral("charge_q_over_r3")) chargeVals.push_back(p.contribution);
                if (p.mechanism == QStringLiteral("mc_lit_valid")) mcVals.push_back(p.contribution);
            }
            for (const PairContribution& p : top) {
                if (rank >= cfg.top_k) break;
                if (!std::isfinite(p.contribution)) continue;
                *topOut << static_cast<qint64>(targetRowId) << ',' << static_cast<qint64>(atom)
                        << ',' << static_cast<qint64>(fs) << ',' << p.mechanism << ','
                        << rank << ',' << p.source_kind << ',' << p.source_id << ','
                        << num(p.contribution) << '\n';
                ++rank;
                ++stats.top_source_query_rows;
            }
            DirectFeatures direct = directFeatures(body, atom, h5Row, orig, LocalFrame{});
            *domOut << static_cast<qint64>(targetRowId) << ',' << static_cast<qint64>(atom)
                    << ',' << static_cast<qint64>(fs) << ',' << num(dominanceFraction(ringVals))
                    << ',' << num(dominanceFraction(chargeVals)) << ',' << num(dominanceFraction(mcVals))
                    << ',' << (direct.mopac_coulomb_present ? "1" : "nan")
                    << ',' << (direct.aimnet2_embedding_present ? "1" : "nan") << '\n';
            ++stats.dominance_query_rows;
            *domBuild4Out << static_cast<qint64>(targetRowId) << ',' << static_cast<qint64>(atom)
                          << ',' << static_cast<qint64>(fs) << ','
                          << num(dominanceFraction(fieldVals)) << '\n';
            ++stats.dominance_build4_query_rows;
        }
    }

    const std::size_t readerFrames = std::min(cfg.reader_pair_frame_slots, rows.size());
    const std::size_t readerAtom = std::min(cfg.reader_pair_atom, nAtoms == 0 ? std::size_t{0} : nAtoms - 1);
    for (std::size_t fs = 0; fs < readerFrames; ++fs) {
        const std::size_t h5Row = rows[fs];
        const std::size_t orig = body.run.frameMap.originalIndex(h5Row);
        const std::size_t targetRowId = fs * nAtoms + readerAtom;
        for (const PairContribution& p : PerAtomRowPairContributions(body, readerAtom, h5Row, cfg, LocalFrame{})) {
            if (!isLegacyPairMechanism(p)) continue;
            writePairRows(*readerOut, targetRowId, readerAtom, fs, h5Row, orig, p);
            ++stats.reader_pair_query_rows;
        }
    }

    pairsOut->flush();
    topOut->flush();
    domOut->flush();
    domBuild4Out->flush();
    readerOut->flush();
    const bool filesOk = pairsFile->commit() && topFile->commit() && domFile->commit()
                         && domBuild4File->commit() && readerFile->commit();
    const bool arraysOk =
        writeNpy<double>(QStringLiteral("%1/query_results/pairs_geometry.npy").arg(outDir),
                         {stats.pair_query_rows, 7}, geom, QByteArray("<f8"))
        && writeNpy<double>(QStringLiteral("%1/query_results/pairs_kernel_T2.npy").arg(outDir),
                            {stats.pair_query_rows, 5}, kt2, QByteArray("<f8"))
        && writeNpy<double>(QStringLiteral("%1/query_results/pairs_kernel_T0.npy").arg(outDir),
                            {stats.pair_query_rows}, kt0, QByteArray("<f8"))
        && writeNpy<double>(QStringLiteral("%1/query_results/pairs_contribution.npy").arg(outDir),
                            {stats.pair_query_rows}, contrib, QByteArray("<f8"));
    QJsonArray pairCols;
    for (const QString& c : {QStringLiteral("disp_x"), QStringLiteral("disp_y"), QStringLiteral("disp_z"),
                             QStringLiteral("r"), QStringLiteral("inv_r3"), QStringLiteral("cos_theta"),
                             QStringLiteral("dipolar"), QStringLiteral("kernel_T0"),
                             QStringLiteral("kernel_T2_0..4"), QStringLiteral("contribution")})
        pairCols.append(c);
    const bool manifestsOk =
        writeQueryManifest(outDir, QStringLiteral("pairs"), QStringLiteral("mixed"), stats.pair_query_rows, pairCols)
        && writeQueryManifest(outDir, QStringLiteral("top_sources"), QStringLiteral("mixed"),
                              stats.top_source_query_rows, QJsonArray{QStringLiteral("contribution")})
        && writeQueryManifest(outDir, QStringLiteral("dominance_fractions"), QStringLiteral("mixed"),
                              stats.dominance_query_rows, QJsonArray{QStringLiteral("dominance_fraction")})
        && writeQueryManifest(outDir, QStringLiteral("dominance_fractions_build4"),
                              QStringLiteral("field"), stats.dominance_build4_query_rows,
                              QJsonArray{QStringLiteral("field_mopac_coulomb_dom_frac")})
        && writeQueryManifest(outDir, QStringLiteral("reader_pairs"), QStringLiteral("mixed"),
                              stats.reader_pair_query_rows, pairCols);
    return filesOk && arraysOk && manifestsOk;
}

class SnapshotCache {
public:
    explicit SnapshotCache(const Body& body) : body_(&body) {
        if (body.run.protein) atomCount_ = body.run.protein->atomCount();
    }

    std::shared_ptr<const model::QtConformationSnapshot> get(std::size_t h5Row) {
        if (!body_ || !body_->run.conformation) return nullptr;
        if (!have_ || h5Row != h5Row_) {
            body_->run.conformation->requestSnapshot(h5Row);
            snapshot_ = body_->run.conformation->snapshot(h5Row);
            rebuildMopacBondAggregates();
            h5Row_ = h5Row;
            have_ = true;
        }
        return snapshot_;
    }

    const std::array<double, 4>* mopacBondAggregate(std::size_t atom) const {
        if (atom >= mopacBondAggregates_.size()) return nullptr;
        return &mopacBondAggregates_[atom];
    }

private:
    void rebuildMopacBondAggregates() {
        mopacBondAggregates_.assign(atomCount_, {kNaN, kNaN, kNaN, kNaN});
        if (!snapshot_ || !snapshot_->has(io::FieldKind::MOPACBondOrders)) return;
        const model::NpyColumn& col = snapshot_->column(io::FieldKind::MOPACBondOrders);
        if (col.cols < 3) return;
        std::vector<double> sums(atomCount_, 0.0);
        std::vector<double> maxes(atomCount_, -std::numeric_limits<double>::infinity());
        std::vector<int> counts(atomCount_, 0);
        for (int r = 0; r < col.rows; ++r) {
            const double* row = col.row(static_cast<std::size_t>(r));
            if (!finiteRaw(row, 3)) continue;
            const std::size_t a = static_cast<std::size_t>(row[0]);
            const std::size_t b = static_cast<std::size_t>(row[1]);
            const double order = row[2];
            auto push = [&](std::size_t atom) {
                if (atom >= atomCount_) return;
                ++counts[atom];
                sums[atom] += order;
                maxes[atom] = std::max(maxes[atom], order);
            };
            push(a);
            push(b);
        }
        for (std::size_t atom = 0; atom < atomCount_; ++atom) {
            if (counts[atom] <= 0) continue;
            mopacBondAggregates_[atom] = {
                static_cast<double>(counts[atom]),
                sums[atom],
                sums[atom] / static_cast<double>(counts[atom]),
                maxes[atom],
            };
        }
    }

    const Body* body_ = nullptr;
    bool have_ = false;
    std::size_t h5Row_ = 0;
    std::size_t atomCount_ = 0;
    std::shared_ptr<const model::QtConformationSnapshot> snapshot_;
    std::vector<std::array<double, 4>> mopacBondAggregates_;
};

}  // namespace

QStringList PerAtomSubstrateSidecars(const PerAtomSubstrateConfig& config) {
    QStringList out;
    out.reserve(static_cast<int>(kPerAtomSubstrateAlwaysSidecars.size() + 1));
    for (const std::string_view stem : kPerAtomSubstrateAlwaysSidecars)
        out << QString::fromLatin1(stem.data(), static_cast<qsizetype>(stem.size()));
    if (config.emit_embedding)
        out << QString::fromLatin1(kPerAtomSubstrateEmbeddingSidecar.data(),
                                   static_cast<qsizetype>(kPerAtomSubstrateEmbeddingSidecar.size()));
    return out;
}

QMap<QString, std::size_t> PerAtomSubstrateFeatureSupport(const PerAtomSubstrateStats& stats) {
    return {
        {QStringLiteral("dft_present_rows"), stats.dft_present},
        {QStringLiteral("ring_jb_present_rows"), stats.ring_present},
        {QStringLiteral("charge_q_over_r3_present_rows"), stats.charge_present},
        {QStringLiteral("mc_lit_valid_present_rows"), stats.mc_lit_valid_present},
        {QStringLiteral("ff14sb_field_present_rows"), stats.ff14sb_field_present},
        {QStringLiteral("apbs_efield_present_rows"), stats.apbs_efield_present},
        {QStringLiteral("apbs_efg_present_rows"), stats.apbs_efg_present},
        {QStringLiteral("aimnet2_charge_present_rows"), stats.aimnet2_charge_present},
        {QStringLiteral("aimnet2_crg_present_rows"), stats.aimnet2_crg_present},
        {QStringLiteral("aimnet2_embedding_present_rows"), stats.aimnet2_embedding_present},
        {QStringLiteral("ff14sb_charge_present_rows"), stats.ff14sb_charge_present},
        {QStringLiteral("mopac_welford_mean_charge_present_rows"), stats.mopac_welford_mean_charge_present},
        {QStringLiteral("charge_complete_rows"), stats.charge_complete},
        {QStringLiteral("mopac_coulomb_shielding_present_rows"), stats.mopac_coulomb_shielding_present},
        {QStringLiteral("mopac_mc_shielding_present_rows"), stats.mopac_mc_shielding_present},
        {QStringLiteral("hbond_count_present_rows"), stats.hbond_count_present},
        {QStringLiteral("hbond_geometry_present_rows"), stats.hbond_geometry_present},
        {QStringLiteral("hm_shielding_present_rows"), stats.hm_shielding_present},
        {QStringLiteral("water_field_present_rows"), stats.water_field_present},
        {QStringLiteral("hydration_shell_present_rows"), stats.hydration_shell_present},
        {QStringLiteral("sasa_present_rows"), stats.sasa_present},
        {QStringLiteral("sasa_normal_present_rows"), stats.sasa_normal_present},
        {QStringLiteral("eeq_charge_present_rows"), stats.eeq_charge_present},
        {QStringLiteral("eeq_coordination_number_present_rows"), stats.eeq_coordination_number_present},
    };
}

PerAtomSubstrateStats RunPerAtomSubstrateEmit(const Body& body,
                                              const QString& outDir,
                                              const PerAtomSubstrateConfig& config,
                                              const DftFrameAlignment& alignment) {
    validateColumnCounts();
    PerAtomSubstrateStats stats;
    if (!body.run.protein || !body.run.trajectory()) return stats;
    const model::QtProtein& p = *body.run.protein;
    stats.atom_count = p.atomCount();
    stats.dft_rows = body.run.frameMap.dftRows().size();
    const std::size_t expectedRows = stats.atom_count * stats.dft_rows;
    constexpr std::size_t kPiece3BAppendFloat64Cols =
        kPerAtomTargetDecompositionCols + kPerAtomRingPathCols + kPerAtomMethodPathCols
        + kPerAtomHbondConditioningCols;
    const double appendGiB =
        static_cast<double>(expectedRows) * static_cast<double>(kPiece3BAppendFloat64Cols)
        * 8.0 / 1073741824.0;
    if (appendGiB > 10.0) {
        throw std::runtime_error(QStringLiteral("per_atom_substrate SIZE-GATE failed: appended flat payload %1 GiB")
                                     .arg(appendGiB, 0, 'f', 3)
                                     .toStdString());
    }
    recordAbsentNewChannelSlabs(body, stats);

    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (!body.catalog.has(ArrayId::MopacChargeWelfordMean)) {
        throw std::runtime_error("Build4 field dominance requires resident MopacChargeWelfordMean; looked in Catalog/trajectory mopacChargeWelford");
    }
    const model::QtEmbeddingTimeSeries* emb = h5 ? h5->aimnet2Embedding() : nullptr;
    if (config.emit_embedding && emb && emb->n_dims != config.embedding_dims) {
        throw std::runtime_error(QStringLiteral("aimnet2_embedding dimension %1 != expected %2")
                                     .arg(static_cast<qulonglong>(emb->n_dims))
                                     .arg(static_cast<qulonglong>(config.embedding_dims))
                                     .toStdString());
    }

    qCInfo(cPerAtom).noquote()
        << "per_atom_substrate | atoms=" << stats.atom_count
        << "| dft rows=" << stats.dft_rows
        << "| ring_cut=" << config.ring_cutoff_A
        << "| bond_cut=" << config.bond_cutoff_A
        << "| charge_cut=" << config.charge_cutoff_A
        << "| frame=raw_lab_frame"
        << "| emit_embedding=" << config.emit_embedding
        << "| piece3b_append_float64_cols=" << kPiece3BAppendFloat64Cols
        << "| piece3b_append_GiB=" << appendGiB;
    if (!stats.absent_new_channel_slabs.isEmpty()) {
        qCWarning(cPerAtom).noquote()
            << "per_atom_substrate absent new-channel slabs |"
            << stats.absent_new_channel_slabs.join(QStringLiteral(", "));
    }

    PerAtomWriter writer(outDir, expectedRows, stats.atom_count, config);
    if (!writer.ok())
        throw std::runtime_error("per_atom_substrate writer open failed");
    writer.writeRingIdentity(body);
    std::vector<WelfordCell> modulation(stats.atom_count * kPerAtomDriverMagnitudeCols);
    std::vector<RowPartitionScratch> partitionRows;
    partitionRows.reserve(expectedRows);
    std::vector<RowDominanceScratch> dominanceRows;
    dominanceRows.reserve(expectedRows);
    SnapshotCache snapshotCache(body);

    std::size_t frameSlot = 0;
    RunTraversal(
        body,
        allAtomsStratum,
        labFrameFn,
        [&](const Body& b, std::size_t atom, std::size_t row, std::size_t orig,
            const FrameResult& fr) {
            const std::shared_ptr<const model::QtConformationSnapshot> snapshot =
                snapshotCache.get(row);
            return directFeatures(b, atom, row, orig, fr.frame, snapshot.get(),
                                  snapshotCache.mopacBondAggregate(atom));
        },
        [&](const Body& b, const AtomState& st, const FrameResult&,
            const DirectFeatures&) {
            return reduceMechanisms(b, st.atom, st.frame, config);
        },
        [&](std::size_t atom, std::size_t row, std::size_t orig, const FrameResult&,
            const DirectFeatures& direct, const MechanismAggregate& agg) {
            const std::size_t thisFrameSlot = writer.rowsWritten() / stats.atom_count;
            const std::array<double, kPerAtomClassicalCols> classical =
                classicalFeatures(agg, direct);
            const std::array<double, kPerAtomDriverMagnitudeCols> mag =
                driverMagnitudes(agg, direct);
            const RowChargeScalars charges = rowChargeScalars(body, atom, row, direct);
            const bool hbondGeometryPresent = direct.hbond_nearest_direction_present
                                              && direct.hbond_flags_present;
            const bool waterFieldPresent = direct.water_efield_present
                                           && direct.water_n_first_present
                                           && direct.water_n_second_present;
            const bool hydrationShellPresent = direct.hydration_half_shell_present
                                               && direct.hydration_dipole_cos_present;
            auditDirectFeatures(stats, direct, charges);
            auditIsolationFeatures(stats, agg.isolation);
            for (std::size_t c = 0; c < kPerAtomDriverMagnitudeCols; ++c)
                modulation[atom * kPerAtomDriverMagnitudeCols + c].push(mag[c]);
            const std::array<double, kPerAtomConditioningCols> conditioning =
                conditioningFeatures(body, atom, row, agg, mag);
            const std::array<double, kPerAtomDominanceCols> dominance =
                dominanceFeatures(agg.isolation);
            RowPartitionScratch scratch;
            scratch.atom = atom;
            scratch.stratum = static_cast<int>(stratumForAtom(body.run.protein->atom(atom)));
            scratch.geometry = {
                conditioning[0],
                conditioning[1],
                conditioning[2],
                conditioning[3],
                agg.isolation.gap_to_2nd_ring_r,
                agg.isolation.gap_to_2nd_charge_r,
                agg.isolation.gap_to_2nd_bond_r,
            };
            scratch.magnitude = mag;
            partitionRows.push_back(scratch);
            RowDominanceScratch domScratch;
            domScratch.stratum = scratch.stratum;
            domScratch.fraction = {
                agg.isolation.dominant_fraction_ring,
                agg.isolation.dominant_fraction_charge,
                agg.isolation.dominant_fraction_mc,
                agg.isolation.dominant_fraction_field,
            };
            dominanceRows.push_back(domScratch);

            if (!writer.writeRow(body, atom, row, thisFrameSlot, orig, direct, agg, charges,
                                 classical, conditioning, dominance)) {
                throw std::runtime_error("per_atom_substrate row write failed");
            }
            ++stats.rows;
            if (direct.dft_present) ++stats.dft_present;
            if (agg.ring_present) ++stats.ring_present;
            if (agg.charge_present) ++stats.charge_present;
            if (agg.mc_lit_valid_present) ++stats.mc_lit_valid_present;
            if (agg.ff14sb_field_present) ++stats.ff14sb_field_present;
            if (direct.apbs_E_present) ++stats.apbs_efield_present;
            if (direct.apbs_efg_present) ++stats.apbs_efg_present;
            if (direct.aimnet2_charge_present) ++stats.aimnet2_charge_present;
            if (direct.aimnet2_crg_present) ++stats.aimnet2_crg_present;
            if (direct.aimnet2_embedding_present) ++stats.aimnet2_embedding_present;
            if (charges.ff14sb.present) ++stats.ff14sb_charge_present;
            if (charges.mopac_welford_mean.present) ++stats.mopac_welford_mean_charge_present;
            if (charges.charge_complete) ++stats.charge_complete;
            if (direct.mopac_coulomb_present) ++stats.mopac_coulomb_shielding_present;
            if (direct.mopac_mc_present) ++stats.mopac_mc_shielding_present;
            if (direct.hbond_count_present) ++stats.hbond_count_present;
            if (hbondGeometryPresent) ++stats.hbond_geometry_present;
            if (direct.hm_shielding_present) ++stats.hm_shielding_present;
            if (waterFieldPresent) ++stats.water_field_present;
            if (hydrationShellPresent) ++stats.hydration_shell_present;
            if (direct.sasa_present) ++stats.sasa_present;
            if (direct.sasa_normal_present) ++stats.sasa_normal_present;
            if (charges.eeq_charge.present) ++stats.eeq_charge_present;
            if (charges.eeq_coordination_number.present) ++stats.eeq_coordination_number_present;
        });
    (void)frameSlot;

    if (stats.rows != expectedRows) {
        throw std::runtime_error(QStringLiteral("per_atom_substrate rows %1 != expected %2")
                                     .arg(static_cast<qulonglong>(stats.rows))
                                     .arg(static_cast<qulonglong>(expectedRows))
                                     .toStdString());
    }
    if (!writer.commit(modulation))
        throw std::runtime_error("per_atom_substrate commit failed");
    if (!writePartitionBins(outDir, partitionRows, modulation, stats.atom_count, stats))
        throw std::runtime_error("per_atom_substrate partition-bin write failed");
    if (!writeDominanceBins(outDir, dominanceRows, stats))
        throw std::runtime_error("per_atom_substrate dominance-bin write failed");
    if (!writeColumnSpecs(outDir, config))
        throw std::runtime_error("per_atom_substrate column spec write failed");
    if (!emitPairQueries(body, outDir, config, stats))
        throw std::runtime_error("per_atom_substrate pair query write failed");
    const CaseHunterStats hunterStats = RunCaseHunter(body, config, outDir);
    stats.hunter_candidate_counts = hunterStats.candidate_counts;
    stats.hunter_anti_circular_assertion = hunterStats.anti_circular_assertion;
    if (!writePerAtomManifest(outDir, stats, config, alignment))
        throw std::runtime_error("per_atom_substrate manifest write failed");

    qCInfo(cPerAtom).noquote()
        << "per_atom_substrate rows | rows=" << stats.rows
        << "| dft_present=" << stats.dft_present
        << "| ring_present=" << stats.ring_present
        << "| charge_present=" << stats.charge_present
        << "| mc_lit_valid_present=" << stats.mc_lit_valid_present
        << "| water_field_present=" << stats.water_field_present
        << "| sasa_present=" << stats.sasa_present
        << "| eeq_charge_present=" << stats.eeq_charge_present
        << "| pair_rows=" << stats.pair_query_rows
        << "| top_source_rows=" << stats.top_source_query_rows
        << "| dominance_rows=" << stats.dominance_query_rows
        << "| reader_pair_rows=" << stats.reader_pair_query_rows;
    return stats;
}

}  // namespace h5reader::rediscover
