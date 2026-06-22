#include "AnalysisAtom.h"

#include "ClassicalSourceMath.h"
#include "DistributionSummary.h"
#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "LiteratureConstants.h"
#include "McConnellLiteratureKernel.h"
#include "RamaRegion.h"
#include "RowDesign.h"
#include "SphericalBasis.h"
#include "SubspaceCompare.h"
#include "TensorConventionGuard.h"
#include "Verbs.h"

#include "../io/QtFieldCatalog.gen.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtConformationSnapshot.h"
#include "../model/QtMopacCoreGroup.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtRing.h"
#include "../model/QtTopology.h"

#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/statistics/linear_regression.hpp>
#include <boost/math/statistics/ljung_box.hpp>
#include <boost/math/statistics/runs_test.hpp>
#include <boost/math/statistics/signal_statistics.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>

#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QLoggingCategory>
#include <QProcess>
#include <QStringList>
#include <QTextStream>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace h5reader::rediscover {
namespace {

Q_LOGGING_CATEGORY(cAnalysisObject, "h5reader.rediscover.analysis_object")

namespace bmstats = boost::math::statistics;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr std::uint32_t kNullSeed = 0xC0DE5EEDu;
constexpr int kRollingWindow = 32;
constexpr int kNullBlock = 16;
// kNullShifts raised 200 -> 2000 for parity with the inventory's circular-shift
// null floor (ACCUMULATOR_RESPEC §4.6: the reads used 2000 shifts -> the 1/2001
// floor). Cost is bounded by the fixed characterization catalogue.
constexpr int kNullShifts = 2000;
constexpr int kMaxChannels = 16;          // structural catalogue cap, not a significance cap.
constexpr int kMaxCollinearityChains = 4;  // bounded named driver-driver chains.
constexpr int kMaxMediationChains = 1;     // the named salt-bridge mediation chain.
constexpr double kDeltaSurvivalRMin = 0.2;  // minimum signed-delta support carried in the gauntlet.
constexpr int kMinContextFrames = 30;    // sparse stratified contexts are marked (§2.4B).
constexpr int kMinThirdFrames = 10;      // segment-thirds edge policy (§4.6.1).
constexpr int kMinDwell = 3;             // debounce a well transition (§4.6.1).
constexpr int kLagWindow = 4;          // lead_lag window, within the sigma-Nyquist (§4.5).
constexpr double kCollinearThreshold = 0.7;  // driver-driver collinearity report floor (§4.6).

bool finite(double v) { return std::isfinite(v); }

double tensorFrobenius(const Mat3& m) {
    return std::sqrt((m.array() * m.array()).sum());
}

double norm5(const std::array<double, 5>& values) {
    double s = 0.0;
    for (double v : values) s += v * v;
    return std::sqrt(s);
}

QJsonValue jd(double v) {
    return finite(v) ? QJsonValue(v) : QJsonValue(QJsonValue::Null);
}

QString csvEscape(QString s) {
    s.replace(QLatin1Char('"'), QStringLiteral("\"\""));
    return QStringLiteral("\"%1\"").arg(s);
}

QString csvNum(double v) {
    return finite(v) ? QString::number(v, 'g', 17) : QString();
}

QString csvBool(bool v) {
    return v ? QStringLiteral("true") : QStringLiteral("false");
}

QString csvSemiDoubles(const std::vector<double>& values) {
    QStringList out;
    for (double v : values) out << csvNum(v);
    return csvEscape(out.join(QLatin1Char(';')));
}

QString csvSemiStrings(const QStringList& values) {
    return csvEscape(values.join(QLatin1Char(';')));
}

QString frameVariantName(FrameVariant v) {
    switch (v) {
    case FrameVariant::Invalid: return QStringLiteral("invalid");
    case FrameVariant::HN_Standard: return QStringLiteral("HN_Standard");
    case FrameVariant::HN_NTerminus: return QStringLiteral("HN_NTerminus");
    case FrameVariant::AromaticHRing: return QStringLiteral("AromaticHRing");
    case FrameVariant::BackboneN: return QStringLiteral("BackboneN");
    case FrameVariant::BackboneN_NTerminus: return QStringLiteral("BackboneN_NTerminus");
    case FrameVariant::BackboneCA: return QStringLiteral("BackboneCA");
    case FrameVariant::BackboneCarbonylC: return QStringLiteral("BackboneCarbonylC");
    case FrameVariant::BackboneCarbonylO: return QStringLiteral("BackboneCarbonylO");
    case FrameVariant::BackboneHA: return QStringLiteral("BackboneHA");
    }
    return QStringLiteral("unknown");
}

QJsonArray vec3Json(const Vec3& v, bool present = true) {
    QJsonArray a;
    for (int i = 0; i < 3; ++i) a.append(present ? jd(v[i]) : QJsonValue(QJsonValue::Null));
    return a;
}

QJsonArray vec5Json(const std::array<double, 5>& v, bool present = true) {
    QJsonArray a;
    for (double x : v) a.append(present ? jd(x) : QJsonValue(QJsonValue::Null));
    return a;
}

template <std::size_t N>
QJsonArray arrayJson(const std::array<double, N>& v, bool present = true) {
    QJsonArray a;
    for (double x : v) a.append(present ? jd(x) : QJsonValue(QJsonValue::Null));
    return a;
}

QJsonArray mat3Json(const Mat3& m, bool present = true) {
    QJsonArray rows;
    for (int r = 0; r < 3; ++r) {
        QJsonArray row;
        for (int c = 0; c < 3; ++c)
            row.append(present ? jd(m(r, c)) : QJsonValue(QJsonValue::Null));
        rows.append(row);
    }
    return rows;
}

QJsonArray boolArrayJson(const std::vector<bool>& values) {
    QJsonArray a;
    for (bool v : values) a.append(v);
    return a;
}

QJsonArray sizeArrayJson(const std::vector<std::size_t>& values) {
    QJsonArray a;
    for (std::size_t v : values) a.append(static_cast<qint64>(v));
    return a;
}

model::SphericalTensor nanTensor() {
    model::SphericalTensor t;
    t.T0 = kNaN;
    t.T1 = {kNaN, kNaN, kNaN};
    t.T2 = {kNaN, kNaN, kNaN, kNaN, kNaN};
    return t;
}

bool tensorFinite(const model::SphericalTensor& t) {
    bool ok = finite(t.T0);
    for (double v : t.T1) ok = ok && finite(v);
    for (double v : t.T2) ok = ok && finite(v);
    return ok;
}

QJsonObject tensorJson(const model::SphericalTensor& t, bool present) {
    QJsonObject o;
    o.insert(QStringLiteral("T0"), present ? jd(t.T0) : QJsonValue(QJsonValue::Null));
    QJsonArray t1;
    for (double v : t.T1) t1.append(present ? jd(v) : QJsonValue(QJsonValue::Null));
    QJsonArray t2;
    for (double v : t.T2) t2.append(present ? jd(v) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("T1"), t1);
    o.insert(QStringLiteral("T2"), t2);
    return o;
}

QString elementName(model::Element e) {
    switch (e) {
    case model::Element::H: return QStringLiteral("H");
    case model::Element::C: return QStringLiteral("C");
    case model::Element::N: return QStringLiteral("N");
    case model::Element::O: return QStringLiteral("O");
    case model::Element::S: return QStringLiteral("S");
    case model::Element::Unknown: return QStringLiteral("Unknown");
    }
    return QStringLiteral("Unknown");
}

QString hybridName(model::Hybridisation h) {
    switch (h) {
    case model::Hybridisation::sp: return QStringLiteral("sp");
    case model::Hybridisation::sp2: return QStringLiteral("sp2");
    case model::Hybridisation::sp3: return QStringLiteral("sp3");
    case model::Hybridisation::Unassigned: return QStringLiteral("Unassigned");
    }
    return QStringLiteral("Unassigned");
}

QString backboneRoleName(model::BackboneRole r) {
    switch (r) {
    case model::BackboneRole::None: return QStringLiteral("None");
    case model::BackboneRole::Nitrogen: return QStringLiteral("Nitrogen");
    case model::BackboneRole::AlphaCarbon: return QStringLiteral("AlphaCarbon");
    case model::BackboneRole::CarbonylCarbon: return QStringLiteral("CarbonylCarbon");
    case model::BackboneRole::CarbonylOxygen: return QStringLiteral("CarbonylOxygen");
    case model::BackboneRole::AmideHydrogen: return QStringLiteral("AmideHydrogen");
    case model::BackboneRole::AlphaHydrogen: return QStringLiteral("AlphaHydrogen");
    }
    return QStringLiteral("None");
}

QString terminalStateName(model::TerminalState s) {
    switch (s) {
    case model::TerminalState::Internal: return QStringLiteral("Internal");
    case model::TerminalState::NTerminus: return QStringLiteral("NTerminus");
    case model::TerminalState::CTerminus: return QStringLiteral("CTerminus");
    case model::TerminalState::NAndCTerminus: return QStringLiteral("NAndCTerminus");
    case model::TerminalState::Unknown: return QStringLiteral("Unknown");
    }
    return QStringLiteral("Unknown");
}

QString aaName(model::AminoAcid aa) {
    static constexpr std::array<const char*, 21> names = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "Unknown",
    };
    const int i = static_cast<int>(aa);
    if (i >= 0 && i < static_cast<int>(names.size())) return QString::fromLatin1(names[static_cast<std::size_t>(i)]);
    return QStringLiteral("Unknown");
}

QString ffAtomTypeName(model::QtFfAtomType t) {
    switch (t) {
    case model::QtFfAtomType::C: return QStringLiteral("C");
    case model::QtFfAtomType::CA: return QStringLiteral("CA");
    case model::QtFfAtomType::CB: return QStringLiteral("CB");
    case model::QtFfAtomType::CC: return QStringLiteral("CC");
    case model::QtFfAtomType::CD: return QStringLiteral("CD");
    case model::QtFfAtomType::CK: return QStringLiteral("CK");
    case model::QtFfAtomType::CM: return QStringLiteral("CM");
    case model::QtFfAtomType::CN: return QStringLiteral("CN");
    case model::QtFfAtomType::CO: return QStringLiteral("CO");
    case model::QtFfAtomType::CP: return QStringLiteral("CP");
    case model::QtFfAtomType::CQ: return QStringLiteral("CQ");
    case model::QtFfAtomType::CR: return QStringLiteral("CR");
    case model::QtFfAtomType::CT: return QStringLiteral("CT");
    case model::QtFfAtomType::CV: return QStringLiteral("CV");
    case model::QtFfAtomType::CW: return QStringLiteral("CW");
    case model::QtFfAtomType::CX: return QStringLiteral("CX");
    case model::QtFfAtomType::CY: return QStringLiteral("CY");
    case model::QtFfAtomType::CZ: return QStringLiteral("CZ");
    case model::QtFfAtomType::Cstar: return QStringLiteral("C*");
    case model::QtFfAtomType::N: return QStringLiteral("N");
    case model::QtFfAtomType::N2: return QStringLiteral("N2");
    case model::QtFfAtomType::N3: return QStringLiteral("N3");
    case model::QtFfAtomType::NA: return QStringLiteral("NA");
    case model::QtFfAtomType::NB: return QStringLiteral("NB");
    case model::QtFfAtomType::NC: return QStringLiteral("NC");
    case model::QtFfAtomType::NP: return QStringLiteral("NP");
    case model::QtFfAtomType::NT: return QStringLiteral("NT");
    case model::QtFfAtomType::NY: return QStringLiteral("NY");
    case model::QtFfAtomType::Nstar: return QStringLiteral("N*");
    case model::QtFfAtomType::H: return QStringLiteral("H");
    case model::QtFfAtomType::H1: return QStringLiteral("H1");
    case model::QtFfAtomType::H2: return QStringLiteral("H2");
    case model::QtFfAtomType::H3: return QStringLiteral("H3");
    case model::QtFfAtomType::H4: return QStringLiteral("H4");
    case model::QtFfAtomType::H5: return QStringLiteral("H5");
    case model::QtFfAtomType::HA: return QStringLiteral("HA");
    case model::QtFfAtomType::HC: return QStringLiteral("HC");
    case model::QtFfAtomType::HO: return QStringLiteral("HO");
    case model::QtFfAtomType::HP: return QStringLiteral("HP");
    case model::QtFfAtomType::HS: return QStringLiteral("HS");
    case model::QtFfAtomType::HW: return QStringLiteral("HW");
    case model::QtFfAtomType::HZ: return QStringLiteral("HZ");
    case model::QtFfAtomType::O: return QStringLiteral("O");
    case model::QtFfAtomType::O2: return QStringLiteral("O2");
    case model::QtFfAtomType::OH: return QStringLiteral("OH");
    case model::QtFfAtomType::OS: return QStringLiteral("OS");
    case model::QtFfAtomType::OW: return QStringLiteral("OW");
    case model::QtFfAtomType::OP: return QStringLiteral("OP");
    case model::QtFfAtomType::OD: return QStringLiteral("OD");
    case model::QtFfAtomType::S: return QStringLiteral("S");
    case model::QtFfAtomType::SH: return QStringLiteral("SH");
    case model::QtFfAtomType::Unknown: return QStringLiteral("Unknown");
    }
    return QStringLiteral("Unknown");
}

QString bondOrderName(model::BondOrder o) {
    switch (o) {
    case model::BondOrder::Single: return QStringLiteral("Single");
    case model::BondOrder::Double: return QStringLiteral("Double");
    case model::BondOrder::Triple: return QStringLiteral("Triple");
    case model::BondOrder::Aromatic: return QStringLiteral("Aromatic");
    case model::BondOrder::Peptide: return QStringLiteral("Peptide");
    case model::BondOrder::Unknown: return QStringLiteral("Unknown");
    }
    return QStringLiteral("Unknown");
}

QString bondCategoryName(model::BondCategory c) {
    switch (c) {
    case model::BondCategory::PeptideCO: return QStringLiteral("PeptideCO");
    case model::BondCategory::PeptideCN: return QStringLiteral("PeptideCN");
    case model::BondCategory::BackboneOther: return QStringLiteral("BackboneOther");
    case model::BondCategory::SidechainCO: return QStringLiteral("SidechainCO");
    case model::BondCategory::Aromatic: return QStringLiteral("Aromatic");
    case model::BondCategory::Disulfide: return QStringLiteral("Disulfide");
    case model::BondCategory::SidechainOther: return QStringLiteral("SidechainOther");
    case model::BondCategory::Unknown: return QStringLiteral("Unknown");
    }
    return QStringLiteral("Unknown");
}

QString ringTypeName(model::RingTypeIndex t) {
    switch (t) {
    case model::RingTypeIndex::PheBenzene: return QStringLiteral("PheBenzene");
    case model::RingTypeIndex::TyrPhenol: return QStringLiteral("TyrPhenol");
    case model::RingTypeIndex::TrpBenzene: return QStringLiteral("TrpBenzene");
    case model::RingTypeIndex::TrpPyrrole: return QStringLiteral("TrpPyrrole");
    case model::RingTypeIndex::TrpPerimeter: return QStringLiteral("TrpPerimeter");
    case model::RingTypeIndex::HisImidazole: return QStringLiteral("HisImidazole");
    case model::RingTypeIndex::HidImidazole: return QStringLiteral("HidImidazole");
    case model::RingTypeIndex::HieImidazole: return QStringLiteral("HieImidazole");
    case model::RingTypeIndex::ProPyrrolidine: return QStringLiteral("ProPyrrolidine");
    }
    return QStringLiteral("Unknown");
}

QString aromaticityName(model::RingAromaticity a) {
    switch (a) {
    case model::RingAromaticity::Full: return QStringLiteral("Full");
    case model::RingAromaticity::Reduced: return QStringLiteral("Reduced");
    case model::RingAromaticity::Weak: return QStringLiteral("Weak");
    case model::RingAromaticity::None: return QStringLiteral("None");
    }
    return QStringLiteral("Unknown");
}

QString ss3ContextName(SecondaryStructure3 ss3) {
    switch (ss3) {
    case SecondaryStructure3::Helix: return QStringLiteral("helix");
    case SecondaryStructure3::Sheet: return QStringLiteral("sheet");
    case SecondaryStructure3::Coil: return QStringLiteral("coil");
    case SecondaryStructure3::Unknown: return QStringLiteral("unknown");
    }
    return QStringLiteral("unknown");
}

QString dihedralBinContextName(int8_t bin) {
    switch (bin) {
    case 0: return QStringLiteral("neg");
    case 1: return QStringLiteral("mid");
    case 2: return QStringLiteral("pos");
    default: return QStringLiteral("unknown");
    }
}

QString dihedralContextName(const ResidentIndexes& idx, std::size_t atom, std::size_t step) {
    const DihedralState phi = idx.dihedrals.state(DihedralKind::Phi, atom, step);
    const DihedralState psi = idx.dihedrals.state(DihedralKind::Psi, atom, step);
    if (phi.present && psi.present)
        return QString::fromLatin1(NameForRowRama(ClassifyRowRama(phi.radians, psi.radians)));
    if (psi.present) return QStringLiteral("psi_%1").arg(dihedralBinContextName(psi.fixed_bin));
    if (phi.present) return QStringLiteral("phi_%1").arg(dihedralBinContextName(phi.fixed_bin));
    return QStringLiteral("unknown");
}

QString nearestContactClassContext(const Body& body, std::size_t atom, std::size_t step) {
    if (!body.run.protein || atom >= body.run.protein->atomCount())
        return QStringLiteral("unknown");
    const model::QtProtein& protein = *body.run.protein;
    const model::QtAtom& target = protein.atom(atom);
    if (target.residueIndex < 0) return QStringLiteral("unknown");
    const Vec3 p = verbs::pos(body, atom, step);
    double best = std::numeric_limits<double>::infinity();
    int bestResidue = -1;
    for (std::size_t other = 0; other < protein.atomCount(); ++other) {
        const model::QtAtom& oa = protein.atom(other);
        if (oa.residueIndex < 0 || oa.residueIndex == target.residueIndex) continue;
        const Vec3 q = verbs::pos(body, other, step);
        const double d2 = (p - q).squaredNorm();
        if (d2 < best) {
            best = d2;
            bestResidue = oa.residueIndex;
        }
    }
    if (bestResidue < 0 || best > 36.0) return QStringLiteral("no_contact");
    return QString::fromLatin1(
        NameForResidueClass(ClassifyResidue(protein.residue(static_cast<std::size_t>(bestResidue)).aminoAcid)));
}

QString uid(const QString& type, std::size_t index) {
    const int width = type == QStringLiteral("atom") ? 5 : 3;
    return QStringLiteral("%1:%2").arg(type).arg(index, width, 10, QLatin1Char('0'));
}

std::pair<std::size_t, std::size_t> canonicalPair(std::size_t a, std::size_t b) {
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
}

std::uint64_t pairKey(std::size_t a, std::size_t b) {
    const auto p = canonicalPair(a, b);
    return (static_cast<std::uint64_t>(p.first) << 32u) | static_cast<std::uint64_t>(p.second);
}

int mechanismOrd(const QString& mechanism) {
    if (mechanism == QStringLiteral("field_mopac_coulomb")) return 0;
    if (mechanism == QStringLiteral("charge_q_over_r3")) return 1;
    if (mechanism == QStringLiteral("mc_lit_valid")) return 2;
    if (mechanism == QStringLiteral("ring_jb")) return 3;
    if (mechanism.contains(QStringLiteral("larsen"))) return 4;
    if (mechanism.contains(QStringLiteral("efg")) || mechanism.contains(QStringLiteral("gradient"))) return 5;
    return -1;
}

QString mechanismName(int ord, const QString& fallback = {}) {
    switch (ord) {
    case 0: return QStringLiteral("field_mopac_coulomb");
    case 1: return QStringLiteral("charge_q_over_r3");
    case 2: return QStringLiteral("mc_lit_valid");
    case 3: return QStringLiteral("ring_jb");
    case 4: return QStringLiteral("hbond_larsen");
    case 5: return QStringLiteral("efg_coulomb_gradient");
    default: return fallback.isEmpty() ? QStringLiteral("unknown") : fallback;
    }
}

int sourceKindOrd(const QString& source) {
    if (source == QStringLiteral("ff14sb_charge_site")) return 0;
    if (source == QStringLiteral("mopac_welford_charge_site")) return 1;
    if (source == QStringLiteral("bond_midpoint")) return 2;
    if (source == QStringLiteral("ring_center")) return 3;
    if (source == QStringLiteral("self")) return 4;
    if (source.contains(QStringLiteral("charge_site"))) return 0;
    return -1;
}

QString sourceKindName(int ord, const QString& fallback = {}) {
    switch (ord) {
    case 0: return QStringLiteral("ff14sb_charge_site");
    case 1: return QStringLiteral("mopac_welford_charge_site");
    case 2: return QStringLiteral("bond_midpoint");
    case 3: return QStringLiteral("ring_center");
    case 4: return QStringLiteral("self");
    default: return fallback.isEmpty() ? QStringLiteral("unknown") : fallback;
    }
}

int scopeOrd(const QString& scope) {
    if (scope == QStringLiteral("self")) return 0;
    if (scope == QStringLiteral("bonded")) return 1;
    if (scope == QStringLiteral("bonded_near_field")) return 2;
    if (scope == QStringLiteral("through_space")) return 3;
    return -1;
}

QString scopeName(int ord) {
    switch (ord) {
    case 0: return QStringLiteral("self");
    case 1: return QStringLiteral("bonded");
    case 2: return QStringLiteral("bonded_near_field");
    case 3: return QStringLiteral("through_space");
    default: return QStringLiteral("unknown");
    }
}

struct ProducerTensorField {
    ArrayId id;
    const char* key;
};

static constexpr std::array<ProducerTensorField, 17> kMcTensorFields = {{
    {ArrayId::McPeptideCoFixed, "mc_peptide_co_fixed"},
    {ArrayId::McPeptideCoBo, "mc_peptide_co_bo"},
    {ArrayId::McPeptideCoRhombic, "mc_peptide_co_rhombic"},
    {ArrayId::McPeptideCnFixed, "mc_peptide_cn_fixed"},
    {ArrayId::McPeptideCnBo, "mc_peptide_cn_bo"},
    {ArrayId::McBackboneOtherFixed, "mc_backbone_other_fixed"},
    {ArrayId::McBackboneOtherBo, "mc_backbone_other_bo"},
    {ArrayId::McSidechainCoFixed, "mc_sidechain_co_fixed"},
    {ArrayId::McSidechainCoBo, "mc_sidechain_co_bo"},
    {ArrayId::McSidechainOtherFixed, "mc_sidechain_other_fixed"},
    {ArrayId::McSidechainOtherBo, "mc_sidechain_other_bo"},
    {ArrayId::McDisulfideFixed, "mc_disulfide_fixed"},
    {ArrayId::McDisulfideBo, "mc_disulfide_bo"},
    {ArrayId::McAromaticZeroedFixed, "mc_aromatic_zeroed_fixed"},
    {ArrayId::McAromaticZeroedBo, "mc_aromatic_zeroed_bo"},
    {ArrayId::McNearestCoT2, "mc_nearest_co_T2"},
    {ArrayId::McNearestCnT2, "mc_nearest_cn_T2"},
}};

struct ScalarSeries {
    explicit ScalarSeries(std::size_t n = 0) : values(n, kNaN) {}
    void reset(std::size_t n) { values.assign(n, kNaN); }
    void set(std::size_t step, double v) {
        if (step < values.size()) values[step] = v;
    }
    QJsonArray json() const {
        QJsonArray a;
        for (double v : values) a.append(jd(v));
        return a;
    }
    std::vector<double> values;
};

struct IntSeries {
    static constexpr int kMissing = std::numeric_limits<int>::min();
    explicit IntSeries(std::size_t n = 0) : values(n, kMissing) {}
    void reset(std::size_t n) { values.assign(n, kMissing); }
    void set(std::size_t step, int v) {
        if (step < values.size()) values[step] = v;
    }
    QJsonArray json() const {
        QJsonArray a;
        for (int v : values)
            a.append(v == kMissing ? QJsonValue(QJsonValue::Null) : QJsonValue(v));
        return a;
    }
    std::vector<int> values;
};

struct BoolTriSeries {
    explicit BoolTriSeries(std::size_t n = 0) : values(n, -1) {}
    void reset(std::size_t n) { values.assign(n, -1); }
    void set(std::size_t step, bool v) {
        if (step < values.size()) values[step] = v ? 1 : 0;
    }
    QJsonArray json() const {
        QJsonArray a;
        for (int v : values)
            a.append(v < 0 ? QJsonValue(QJsonValue::Null) : QJsonValue(v != 0));
        return a;
    }
    std::vector<int> values;
};

struct Vec3Series {
    explicit Vec3Series(std::size_t n = 0) : values(n, Vec3::Zero()), present(n, false) {}
    void reset(std::size_t n) {
        values.assign(n, Vec3::Zero());
        present.assign(n, false);
    }
    void set(std::size_t step, const Vec3& v) {
        if (step >= values.size()) return;
        values[step] = v;
        present[step] = finite(v.x()) && finite(v.y()) && finite(v.z());
    }
    QJsonArray json() const {
        QJsonArray a;
        for (std::size_t i = 0; i < values.size(); ++i)
            a.append(present[i] ? QJsonValue(vec3Json(values[i])) : QJsonValue(QJsonValue::Null));
        return a;
    }
    std::vector<Vec3> values;
    std::vector<bool> present;
};

struct Mat3Series {
    explicit Mat3Series(std::size_t n = 0) : values(n, Mat3::Constant(kNaN)), present(n, false) {}
    void reset(std::size_t n) {
        values.assign(n, Mat3::Constant(kNaN));
        present.assign(n, false);
    }
    void set(std::size_t step, const Mat3& m) {
        if (step >= values.size()) return;
        values[step] = m;
        present[step] = m.allFinite();
    }
    QJsonArray json() const {
        QJsonArray a;
        for (std::size_t i = 0; i < values.size(); ++i)
            a.append(present[i] ? QJsonValue(mat3Json(values[i])) : QJsonValue(QJsonValue::Null));
        return a;
    }
    std::vector<Mat3> values;
    std::vector<bool> present;
};

template <std::size_t N>
struct FixedArraySeries {
    explicit FixedArraySeries(std::size_t n = 0) : values(n), present(n, false) {
        for (auto& v : values) v.fill(kNaN);
    }
    void reset(std::size_t n) {
        values.assign(n, {});
        for (auto& v : values) v.fill(kNaN);
        present.assign(n, false);
    }
    void set(std::size_t step, const std::array<double, N>& v) {
        if (step >= values.size()) return;
        values[step] = v;
        bool ok = true;
        for (double x : v) ok = ok && finite(x);
        present[step] = ok;
    }
    QJsonArray json() const {
        QJsonArray a;
        for (std::size_t i = 0; i < values.size(); ++i)
            a.append(present[i] ? QJsonValue(arrayJson(values[i])) : QJsonValue(QJsonValue::Null));
        return a;
    }
    std::vector<std::array<double, N>> values;
    std::vector<bool> present;
};

struct TensorSeries {
    explicit TensorSeries(std::size_t n = 0) : values(n, nanTensor()), present(n, false) {}
    void reset(std::size_t n) {
        values.assign(n, nanTensor());
        present.assign(n, false);
    }
    void set(std::size_t step, const model::SphericalTensor& t) {
        if (step >= values.size()) return;
        values[step] = t;
        present[step] = finite(t.T0);
    }
    QJsonArray json() const {
        QJsonArray a;
        for (std::size_t i = 0; i < values.size(); ++i)
            a.append(present[i] ? QJsonValue(tensorJson(values[i], true)) : QJsonValue(QJsonValue::Null));
        return a;
    }
    std::vector<model::SphericalTensor> values;
    std::vector<bool> present;
};

struct EfgValue {
    double t0 = kNaN;
    std::array<double, 5> t2 = {kNaN, kNaN, kNaN, kNaN, kNaN};
};

struct EfgSeries {
    explicit EfgSeries(std::size_t n = 0) : values(n), present(n, false) {}
    void reset(std::size_t n) {
        values.assign(n, {});
        present.assign(n, false);
    }
    void set(std::size_t step, EfgValue v) {
        if (step >= values.size()) return;
        values[step] = v;
        bool any = finite(v.t0);
        for (double x : v.t2) any = any || finite(x);
        present[step] = any;
    }
    QJsonArray json() const {
        QJsonArray a;
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (!present[i]) {
                a.append(QJsonValue(QJsonValue::Null));
                continue;
            }
            QJsonObject o;
            o.insert(QStringLiteral("T0"), jd(values[i].t0));
            o.insert(QStringLiteral("T2"), vec5Json(values[i].t2, true));
            a.append(o);
        }
        return a;
    }
    std::vector<EfgValue> values;
    std::vector<bool> present;
};

struct CsaDescriptor {
    bool valid = false;
    Vec3 principal_values = Vec3::Constant(kNaN);   // value order: sigma11 <= sigma22 <= sigma33
    Mat3 pas_axes = Mat3::Constant(kNaN);           // columns follow principal_values
    Vec3 haeberlen_values = Vec3::Constant(kNaN);   // columns/values ordered sigma_xx, sigma_yy, sigma_zz
    Mat3 haeberlen_axes = Mat3::Constant(kNaN);
    double sigma_iso = kNaN;
    double haeberlen_eta = kNaN;
    double haeberlen_span = kNaN;
    double haeberlen_skew = kNaN;
};

struct CsaDescriptorSeries {
    explicit CsaDescriptorSeries(std::size_t n = 0) : values(n) {}
    void reset(std::size_t n) { values.assign(n, {}); }
    void set(std::size_t step, const CsaDescriptor& v) {
        if (step < values.size()) values[step] = v;
    }
    QJsonObject json() const {
        QJsonObject o;
        QJsonArray valid;
        QJsonArray iso;
        QJsonArray eta;
        QJsonArray span;
        QJsonArray skew;
        QJsonArray principal;
        QJsonArray pasAxes;
        QJsonArray haeberlenValues;
        QJsonArray haeberlenAxes;
        for (const CsaDescriptor& v : values) {
            valid.append(v.valid);
            iso.append(v.valid ? jd(v.sigma_iso) : QJsonValue(QJsonValue::Null));
            eta.append(v.valid ? jd(v.haeberlen_eta) : QJsonValue(QJsonValue::Null));
            span.append(v.valid ? jd(v.haeberlen_span) : QJsonValue(QJsonValue::Null));
            skew.append(v.valid ? jd(v.haeberlen_skew) : QJsonValue(QJsonValue::Null));
            principal.append(v.valid ? QJsonValue(vec3Json(v.principal_values))
                                     : QJsonValue(QJsonValue::Null));
            pasAxes.append(v.valid ? QJsonValue(mat3Json(v.pas_axes))
                                   : QJsonValue(QJsonValue::Null));
            haeberlenValues.append(v.valid ? QJsonValue(vec3Json(v.haeberlen_values))
                                           : QJsonValue(QJsonValue::Null));
            haeberlenAxes.append(v.valid ? QJsonValue(mat3Json(v.haeberlen_axes))
                                         : QJsonValue(QJsonValue::Null));
        }
        o.insert(QStringLiteral("valid"), valid);
        o.insert(QStringLiteral("sigma_iso"), iso);
        o.insert(QStringLiteral("haeberlen_eta"), eta);
        o.insert(QStringLiteral("haeberlen_span"), span);
        o.insert(QStringLiteral("haeberlen_skew"), skew);
        o.insert(QStringLiteral("principal_values"), principal);
        o.insert(QStringLiteral("pas_axes"), pasAxes);
        o.insert(QStringLiteral("haeberlen_values"), haeberlenValues);
        o.insert(QStringLiteral("haeberlen_axes"), haeberlenAxes);
        o.insert(QStringLiteral("layout"),
                 QStringLiteral("per-step arrays; axes are 3x3 row-major matrices with lab-frame director columns"));
        o.insert(QStringLiteral("source"), QStringLiteral("sigma.total total_raw symmetric part"));
        return o;
    }
    std::vector<CsaDescriptor> values;
};

int molecularFrameKindOrd(const QString& kind) {
    if (kind == QStringLiteral("none")) return 0;
    if (kind == QStringLiteral("backbone_carbonyl")) return 1;
    if (kind == QStringLiteral("backbone_amide_n")) return 2;
    if (kind == QStringLiteral("aromatic_ring_local")) return 3;
    if (kind == QStringLiteral("met_sd_cg_s_ce")) return 4;
    if (kind == QStringLiteral("sidechain_carboxylate")) return 5;
    if (kind == QStringLiteral("sidechain_guanidinium")) return 6;
    if (kind == QStringLiteral("sidechain_carboxamide")) return 7;
    if (kind == QStringLiteral("backbone_amide_h")) return 8;
    if (kind == QStringLiteral("backbone_ca")) return 9;
    if (kind == QStringLiteral("aliphatic_carbon")) return 10;
    if (kind == QStringLiteral("hydroxyl_oxygen")) return 11;
    return -1;
}

struct MolecularFrameValue {
    bool valid = false;
    Mat3 axes = Mat3::Constant(kNaN);  // columns are x, y, z in lab coordinates
};

struct MolecularFrameSeries {
    explicit MolecularFrameSeries(std::size_t n = 0) : values(n) {}
    void reset(std::size_t n) { values.assign(n, {}); }
    void setMetadata(QString kindValue,
                     std::vector<int32_t> anchorValues,
                     FrameVariant variantValue = FrameVariant::Invalid,
                     QString ringUidValue = {},
                     QString ringTypeValue = {}) {
        kind = std::move(kindValue);
        kind_ord = molecularFrameKindOrd(kind);
        anchors = std::move(anchorValues);
        variant = frameVariantName(variantValue);
        variant_ord = static_cast<int>(variantValue);
        ring_uid = std::move(ringUidValue);
        ring_type = std::move(ringTypeValue);
    }
    void set(std::size_t step, const Mat3& axesValue) {
        if (step >= values.size()) return;
        MolecularFrameValue v;
        v.axes = axesValue;
        v.valid = axesValue.allFinite();
        values[step] = v;
    }
    QJsonObject json() const {
        QJsonObject o;
        o.insert(QStringLiteral("frame_kind"), kind);
        o.insert(QStringLiteral("frame_kind_ord"), kind_ord);
        o.insert(QStringLiteral("frame_variant"), variant);
        o.insert(QStringLiteral("frame_variant_ord"), variant_ord);
        QJsonArray anchorArray;
        for (int32_t a : anchors) anchorArray.append(a);
        o.insert(QStringLiteral("anchors"), anchorArray);
        if (!ring_uid.isEmpty()) o.insert(QStringLiteral("ring_uid"), ring_uid);
        if (!ring_type.isEmpty()) o.insert(QStringLiteral("ring_type"), ring_type);

        QJsonArray valid;
        QJsonArray axesArray;
        for (const MolecularFrameValue& v : values) {
            valid.append(v.valid);
            axesArray.append(v.valid ? QJsonValue(mat3Json(v.axes))
                                     : QJsonValue(QJsonValue::Null));
        }
        o.insert(QStringLiteral("valid"), valid);
        o.insert(QStringLiteral("axes"), axesArray);
        o.insert(QStringLiteral("layout"),
                 QStringLiteral("per-step arrays; axes[t][row][col] with columns x,y,z in lab coordinates"));
        o.insert(QStringLiteral("source"), QStringLiteral("ported pas_projection.molecular_frame() geometry"));
        return o;
    }

    QString kind = QStringLiteral("none");
    int kind_ord = 0;
    QString variant = QStringLiteral("invalid");
    int variant_ord = static_cast<int>(FrameVariant::Invalid);
    std::vector<int32_t> anchors;
    QString ring_uid;
    QString ring_type;
    std::vector<MolecularFrameValue> values;
};

std::optional<Vec3> normalizeFrameVec(const Vec3& v, double eps = 1e-12) {
    const double norm = v.norm();
    if (!finite(norm) || !(norm > eps)) return std::nullopt;
    return v / norm;
}

std::optional<Mat3> frameFromColumns(const Vec3& x, const Vec3& y, const Vec3& z) {
    Mat3 frame;
    frame.col(0) = x;
    frame.col(1) = y;
    frame.col(2) = z;
    if (!frame.allFinite()) return std::nullopt;
    return frame;
}

std::optional<Mat3> frameFromXAndPlane(const Vec3& xVec, const Vec3& planeVec) {
    const auto xOpt = normalizeFrameVec(xVec);
    if (!xOpt) return std::nullopt;
    const Vec3& x = *xOpt;
    const Vec3 plane = planeVec - planeVec.dot(x) * x;
    const auto y0Opt = normalizeFrameVec(plane);
    if (!y0Opt) return std::nullopt;
    const auto zOpt = normalizeFrameVec(x.cross(*y0Opt));
    if (!zOpt) return std::nullopt;
    const auto yOpt = normalizeFrameVec(zOpt->cross(x));
    if (!yOpt) return std::nullopt;
    return frameFromColumns(x, *yOpt, *zOpt);
}

std::optional<Mat3> frameFromXAndZ(const Vec3& xVec, const Vec3& zVec) {
    const auto zOpt = normalizeFrameVec(zVec);
    if (!zOpt) return std::nullopt;
    const Vec3& z = *zOpt;
    const Vec3 xProj = xVec - xVec.dot(z) * z;
    const auto xOpt = normalizeFrameVec(xProj);
    if (!xOpt) return std::nullopt;
    const auto yOpt = normalizeFrameVec(z.cross(*xOpt));
    if (!yOpt) return std::nullopt;
    return frameFromColumns(*xOpt, *yOpt, z);
}

std::optional<Mat3> frameFromXAndPlaneLocked(const Vec3& xVec,
                                             const Vec3& planeVec,
                                             std::optional<Vec3>& prevX,
                                             std::optional<Vec3>& prevZ) {
    auto xOpt = normalizeFrameVec(xVec);
    if (!xOpt) return std::nullopt;
    Vec3 x = *xOpt;
    const Vec3 plane = planeVec - planeVec.dot(x) * x;
    const auto y0Opt = normalizeFrameVec(plane);
    if (!y0Opt) return std::nullopt;
    auto zOpt = normalizeFrameVec(x.cross(*y0Opt));
    if (!zOpt) return std::nullopt;
    Vec3 z = *zOpt;

    if (prevX && prevX->dot(x) < 0.0) x *= -1.0;
    if (prevZ && prevZ->dot(z) < 0.0) z *= -1.0;

    const auto yOpt = normalizeFrameVec(z.cross(x));
    if (!yOpt) return std::nullopt;
    const auto frame = frameFromColumns(x, *yOpt, z);
    if (!frame) return std::nullopt;
    prevX = x;
    prevZ = z;
    return frame;
}

double sigmaRoundTripMaxAbs(const Mat3& axes, const Mat3& raw) {
    long double x[3] = {
        static_cast<long double>(axes(0, 0)),
        static_cast<long double>(axes(1, 0)),
        static_cast<long double>(axes(2, 0)),
    };
    long double z[3] = {
        static_cast<long double>(axes(0, 2)),
        static_cast<long double>(axes(1, 2)),
        static_cast<long double>(axes(2, 2)),
    };
    auto normalize = [](long double (&v)[3]) -> bool {
        const long double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (!(n > 0.0L) || !std::isfinite(static_cast<double>(n))) return false;
        v[0] /= n;
        v[1] /= n;
        v[2] /= n;
        return true;
    };
    if (!normalize(x)) return kNaN;
    const long double zx = z[0] * x[0] + z[1] * x[1] + z[2] * x[2];
    z[0] -= zx * x[0];
    z[1] -= zx * x[1];
    z[2] -= zx * x[2];
    if (!normalize(z)) return kNaN;
    long double y[3] = {
        z[1] * x[2] - z[2] * x[1],
        z[2] * x[0] - z[0] * x[2],
        z[0] * x[1] - z[1] * x[0],
    };
    if (!normalize(y)) return kNaN;
    long double a[3][3] = {
        {x[0], y[0], z[0]},
        {x[1], y[1], z[1]},
        {x[2], y[2], z[2]},
    };
    long double local[3][3]{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            long double sum = 0.0L;
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    sum += a[k][i]
                           * static_cast<long double>(raw(k, l))
                           * a[l][j];
                }
            }
            local[i][j] = sum;
        }
    }
    long double maxAbs = 0.0L;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            long double sum = 0.0L;
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    sum += a[i][k]
                           * local[k][l]
                           * a[j][l];
                }
            }
            const long double err = std::abs(sum - static_cast<long double>(raw(i, j)));
            if (err > maxAbs) maxAbs = err;
        }
    }
    return static_cast<double>(maxAbs);
}

std::vector<double> componentSeries(const ScalarSeries& s) { return s.values; }

std::vector<double> componentSeries(const Vec3Series& s, int c) {
    std::vector<double> out(s.values.size(), kNaN);
    for (std::size_t i = 0; i < s.values.size(); ++i)
        if (s.present[i]) out[i] = s.values[i][c];
    return out;
}

std::vector<double> componentSeriesT0(const TensorSeries& s) {
    std::vector<double> out(s.values.size(), kNaN);
    for (std::size_t i = 0; i < s.values.size(); ++i)
        if (s.present[i]) out[i] = s.values[i].T0;
    return out;
}

struct RunningSeriesRef {
    QString series_ref;
    QJsonValue component_ref = QJsonValue(QJsonValue::Null);
    QJsonValue relationship_index = QJsonValue(QJsonValue::Null);
    QJsonValue relationship_facet = QJsonValue(QJsonValue::Null);
    std::vector<double> values;
};

struct OlsResult {
    double slope = kNaN;
    double intercept = kNaN;
};

OlsResult ols(const std::vector<double>& x, const std::vector<double>& y) {
    OlsResult r;
    if (x.size() != y.size() || x.size() < 2) return r;
    const double mx = bmstats::mean(x);
    const double my = bmstats::mean(y);
    double sxx = 0.0;
    double sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - mx;
        sxx += dx * dx;
        sxy += dx * (y[i] - my);
    }
    if (!(sxx > 0.0)) return r;
    r.slope = sxy / sxx;
    r.intercept = my - r.slope * mx;
    return r;
}

QJsonObject coverageJson(const std::vector<double>& x) {
    QJsonObject o;
    if (x.empty()) {
        o.insert(QStringLiteral("facet_min"), QJsonValue(QJsonValue::Null));
        o.insert(QStringLiteral("facet_max"), QJsonValue(QJsonValue::Null));
        o.insert(QStringLiteral("sampled_range_frac"), QJsonValue(QJsonValue::Null));
        return o;
    }
    const auto [mn, mx] = std::minmax_element(x.begin(), x.end());
    o.insert(QStringLiteral("facet_min"), jd(*mn));
    o.insert(QStringLiteral("facet_max"), jd(*mx));
    if (*mx == *mn) {
        o.insert(QStringLiteral("sampled_range_frac"), QJsonValue(QJsonValue::Null));
        return o;
    }
    const int bins = std::min<int>(16, std::max<int>(1, static_cast<int>(x.size())));
    std::vector<bool> occupied(static_cast<std::size_t>(bins), false);
    for (double v : x) {
        int b = static_cast<int>(std::floor((v - *mn) / (*mx - *mn) * bins));
        if (b == bins) --b;
        if (b >= 0 && b < bins) occupied[static_cast<std::size_t>(b)] = true;
    }
    const int occ = static_cast<int>(std::count(occupied.begin(), occupied.end(), true));
    o.insert(QStringLiteral("sampled_range_frac"), static_cast<double>(occ) / bins);
    return o;
}

double leverageTop1(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return kNaN;
    const double mx = bmstats::mean(x);
    const double my = bmstats::mean(y);
    std::vector<double> contrib;
    contrib.reserve(x.size());
    double sumAbs = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double c = (x[i] - mx) * (y[i] - my);
        contrib.push_back(std::abs(c));
        sumAbs += std::abs(c);
    }
    if (!(sumAbs > 0.0)) return kNaN;
    return *std::max_element(contrib.begin(), contrib.end()) / sumAbs;
}

double pchipMidpoint(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return kNaN;
    std::vector<std::pair<double, double>> pairs;
    pairs.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) pairs.emplace_back(x[i], y[i]);
    std::sort(pairs.begin(), pairs.end());
    std::vector<double> ux;
    std::vector<double> uy;
    for (const auto& p : pairs) {
        if (!ux.empty() && std::abs(p.first - ux.back()) < 1e-12) {
            uy.back() = 0.5 * (uy.back() + p.second);
            continue;
        }
        ux.push_back(p.first);
        uy.push_back(p.second);
    }
    if (ux.size() < 2) return kNaN;
    try {
        boost::math::interpolators::pchip<std::vector<double>> interp(std::move(ux), std::move(uy));
        const double minX = pairs.front().first;
        const double maxX = pairs.back().first;
        return interp(0.5 * (minX + maxX));
    } catch (...) {
        return kNaN;
    }
}

double chatterjeeXi(std::vector<double> x, std::vector<double> y) {
    if (x.size() != y.size() || x.size() < 2) return kNaN;
    std::vector<std::pair<double, double>> pairs;
    pairs.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) pairs.emplace_back(x[i], y[i]);
    std::sort(pairs.begin(), pairs.end());
    for (std::size_t i = 0; i < pairs.size(); ++i) y[i] = pairs[i].second;

    std::vector<std::pair<double, std::size_t>> ranked;
    ranked.reserve(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) ranked.emplace_back(y[i], i);
    std::sort(ranked.begin(), ranked.end());

    std::vector<std::size_t> ranks(y.size(), 0);
    std::size_t i = 0;
    while (i < ranked.size()) {
        std::size_t j = i + 1;
        while (j < ranked.size() && ranked[j].first == ranked[i].first) ++j;
        for (std::size_t k = i; k < j; ++k) ranks[ranked[k].second] = i;
        i = j;
    }

    std::size_t sum = 0;
    for (std::size_t k = 1; k < ranks.size(); ++k)
        sum += ranks[k] > ranks[k - 1] ? ranks[k] - ranks[k - 1] : ranks[k - 1] - ranks[k];

    const double denom = static_cast<double>(ranks.size() * ranks.size() - 1);
    if (!(denom > 0.0)) return kNaN;
    const double xi = 1.0 - (3.0 * static_cast<double>(sum) / denom);
    return std::abs(xi - 1.0) < std::numeric_limits<double>::epsilon() ? kNaN : xi;
}

std::tuple<double, double, double, double> nullShiftStats(const std::vector<double>& x,
                                                          const std::vector<double>& y,
                                                          double observedSlope) {
    if (x.size() != y.size() || x.size() < 2) return {kNaN, kNaN, kNaN, kNaN};
    std::mt19937 rng(kNullSeed);
    const int n = static_cast<int>(x.size());
    const int blocks = std::max(1, (n + kNullBlock - 1) / kNullBlock);
    std::uniform_int_distribution<int> dist(1, std::max(1, blocks - 1));
    std::vector<double> slopes;
    slopes.reserve(kNullShifts);
    std::vector<double> shifted(x.size(), kNaN);
    int ge = 0;
    for (int i = 0; i < kNullShifts; ++i) {
        const int shift = (dist(rng) * kNullBlock) % n;
        for (int j = 0; j < n; ++j) shifted[static_cast<std::size_t>(j)] = x[static_cast<std::size_t>((j + shift) % n)];
        const OlsResult r = ols(shifted, y);
        if (!finite(r.slope)) continue;
        slopes.push_back(r.slope);
        if (finite(observedSlope) && std::abs(r.slope) >= std::abs(observedSlope)) ++ge;
    }
    if (slopes.empty()) return {kNaN, kNaN, kNaN, kNaN};
    const double mean = bmstats::mean(slopes);
    double variance = 0.0;
    if (slopes.size() > 1) {
        for (double slope : slopes) {
            const double d = slope - mean;
            variance += d * d;
        }
        variance /= static_cast<double>(slopes.size() - 1);
    }
    const double sd = slopes.size() > 1 ? std::sqrt(variance) : 0.0;
    const double z = (sd > 0.0 && finite(observedSlope)) ? (observedSlope - mean) / sd : kNaN;
    const double p = finite(observedSlope)
        ? static_cast<double>(ge + 1) / static_cast<double>(slopes.size() + 1)
        : kNaN;
    return {mean, sd, z, p};
}

// ---------------------------------------------------------------------------
// NEW statistics (ACCUMULATOR_RESPEC §4.6). Each carries its inventory
// cross-check in a comment. They reuse the certified bmstats primitives
// (mean/variance/covariance) but are otherwise new code (the engine's bmstats
// has only mean/var/cov/ljung_box/runs/snr).
// ---------------------------------------------------------------------------

// Pull the finite paired (driver, sigma) samples on the sigma rows, in time
// order. This is the single gather every response-law/gauntlet statistic runs
// on, so the right-frame pairing happens once.
void pairedOnSigmaRows(const std::vector<double>& sigma,
                       const std::vector<double>& driver,
                       const std::vector<std::size_t>& sigmaRows,
                       std::vector<double>& xOut,
                       std::vector<double>& yOut) {
    xOut.clear();
    yOut.clear();
    xOut.reserve(sigmaRows.size());
    yOut.reserve(sigmaRows.size());
    for (std::size_t row : sigmaRows) {
        if (row >= sigma.size() || row >= driver.size()) continue;
        if (!finite(sigma[row]) || !finite(driver[row])) continue;
        xOut.push_back(driver[row]);  // x = driver
        yOut.push_back(sigma[row]);   // y = sigma
    }
}

// Pearson r on a paired series. r = cov / sqrt(var_x * var_y).
// Cross-check (§4.6): ASP7 CG efg_xz->sigma_xz r=0.952; LYS30 C sigma_yy via
// efg_zz r~=0.83; carbonyl C lab 0.14 -> mol 0.75.
double pearsonR(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return kNaN;
    const double mx = bmstats::mean(x);
    const double my = bmstats::mean(y);
    double sxx = 0.0;
    double syy = 0.0;
    double sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    const double denom = std::sqrt(sxx * syy);
    if (!(denom > 0.0)) return kNaN;
    return sxy / denom;
}

// Segment-thirds weakest-third |r|. Split the paired series into 3 equal
// contiguous thirds; r in each; report the min-|r|.
// Cross-check (§4.6): carbonyl C' weakest-third 0.63; ASP7 CG seg-min 0.940;
// the ASP29 OD2 demotion fires here.
double segMinR(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 6) return kNaN;
    const std::size_t n = x.size();
    const std::size_t first = (n + 2) / 3;  // ceil(n/3)
    const std::size_t second = first;
    if (n <= first + second) return kNaN;
    const std::array<std::size_t, 3> sizes = {first, second, n - first - second};
    double minAbs = std::numeric_limits<double>::infinity();
    double minSigned = kNaN;
    std::size_t lo = 0;
    for (int seg = 0; seg < 3; ++seg) {
        const std::size_t hi = lo + sizes[static_cast<std::size_t>(seg)];
        if (hi > n || hi <= lo || (hi - lo) < static_cast<std::size_t>(kMinThirdFrames)) return kNaN;
        std::vector<double> sx(x.begin() + static_cast<std::ptrdiff_t>(lo),
                               x.begin() + static_cast<std::ptrdiff_t>(hi));
        std::vector<double> sy(y.begin() + static_cast<std::ptrdiff_t>(lo),
                               y.begin() + static_cast<std::ptrdiff_t>(hi));
        const double r = pearsonR(sx, sy);
        if (!finite(r)) return kNaN;
        if (std::abs(r) < minAbs) {
            minAbs = std::abs(r);
            minSigned = r;
        }
        lo = hi;
    }
    return finite(minSigned) ? minSigned : kNaN;
}

// First-difference r on the paired series (the differenced-pair Pearson r).
double deltaR(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 3) return kNaN;
    std::vector<double> dx;
    std::vector<double> dy;
    dx.reserve(x.size() - 1);
    dy.reserve(y.size() - 1);
    for (std::size_t i = 1; i < x.size(); ++i) {
        dx.push_back(x[i] - x[i - 1]);
        dy.push_back(y[i] - y[i - 1]);
    }
    return pearsonR(dx, dy);
}

// delta_survives: the first-difference relationship retains the level sign AND
// |r_delta| >= floor (the co-drift demotion). Cross-check (§4.6): the 206/274
// co-drift nulls collapse (ASP29 OD1 level 0.73 -> delta 0.04 = FALSE); the
// carboxylate delta >= level = TRUE.
bool deltaSurvives(const std::vector<double>& x, const std::vector<double>& y, double levelR) {
    const double rd = deltaR(x, y);
    if (!finite(rd) || !finite(levelR)) return false;
    const bool signKept = (rd > 0.0) == (levelR > 0.0) || levelR == 0.0;
    return signKept && std::abs(rd) >= kDeltaSurvivalRMin;
}

// Plain autocorrelation function at a lag.
double acfAtLag(const std::vector<double>& v, int lag) {
    const int n = static_cast<int>(v.size());
    if (lag <= 0 || lag >= n) return kNaN;
    const double m = bmstats::mean(v);
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; ++i) {
        const double d = v[static_cast<std::size_t>(i)] - m;
        den += d * d;
        if (i + lag < n) num += d * (v[static_cast<std::size_t>(i + lag)] - m);
    }
    if (!(den > 0.0)) return kNaN;
    return num / den;
}

struct DwellResult {
    double mean_dwell_frames = kNaN;
    int n_transitions = -1;
    double autocorr_time = kNaN;  // first lag where ACF < 1/e.
};

// dwell / n_transitions / autocorr_time on a scalar series. A "well" is the
// sign of (value - median); a transition is a sign flip; dwell = mean run
// length; autocorr_time = first lag where ACF < 1/e (MEASURED, never a deflator).
// Cross-check (§4.6): MET23 chi3 dwell; TYR46 2 transitions; HIS18 chi2 102.
DwellResult dwellStats(const std::vector<double>& series) {
    DwellResult out;
    std::vector<double> v;
    v.reserve(series.size());
    for (double x : series)
        if (finite(x)) v.push_back(x);
    if (v.size() < 3) return out;
    std::vector<double> sorted = v;
    std::sort(sorted.begin(), sorted.end());
    const double median = sorted[sorted.size() / 2];
    std::vector<int> wells;
    wells.reserve(v.size());
    for (double x : v) wells.push_back(x >= median ? 1 : 0);
    std::vector<int> rawRunLengths;
    int run = 1;
    for (std::size_t i = 1; i < wells.size(); ++i) {
        if (wells[i] != wells[i - 1]) {
            rawRunLengths.push_back(run);
            run = 1;
        } else {
            ++run;
        }
    }
    rawRunLengths.push_back(run);

    std::vector<int> runLengths;
    for (int len : rawRunLengths) {
        if (!runLengths.empty() && len < kMinDwell)
            runLengths.back() += len;
        else
            runLengths.push_back(len);
    }
    const int transitions = runLengths.empty() ? 0 : static_cast<int>(runLengths.size()) - 1;
    out.n_transitions = transitions;
    double sum = 0.0;
    for (int r : runLengths) sum += r;
    out.mean_dwell_frames = runLengths.empty() ? kNaN : sum / static_cast<double>(runLengths.size());
    const double oneOverE = 1.0 / std::exp(1.0);
    for (int lag = 1; lag < static_cast<int>(v.size()); ++lag) {
        const double acf = acfAtLag(v, lag);
        if (finite(acf) && acf < oneOverE) {
            out.autocorr_time = static_cast<double>(lag);
            break;
        }
    }
    return out;
}

// Circular wells for a dihedral series: the standard 3-rotamer bins by the
// nearest of {+60, 180, -60} degrees (radians). Returns -1 where absent.
std::vector<int> circularWells(const std::vector<double>& radians) {
    std::vector<int> wells(radians.size(), -1);
    constexpr double pi = 3.14159265358979323846264338327950288;
    const std::array<double, 3> centres = {pi / 3.0, pi, -pi / 3.0};  // +60, 180, -60
    for (std::size_t i = 0; i < radians.size(); ++i) {
        if (!finite(radians[i])) continue;
        int best = 0;
        double bestDist = std::numeric_limits<double>::infinity();
        for (int w = 0; w < 3; ++w) {
            double d = std::abs(radians[i] - centres[static_cast<std::size_t>(w)]);
            d = std::min(d, 2.0 * pi - d);  // circular distance
            if (d < bestDist) {
                bestDist = d;
                best = w;
            }
        }
        wells[i] = best;
    }
    return wells;
}

// eta^2 (one-way between-well SS / total SS) of a response series stratified by
// circular wells of a driver, both gathered on the sigma rows.
// Cross-check (§4.6): MET23 chi3->SD T2norm eta^2=0.58; PRO9 pucker 0.63/0.70/0.81.
double etaSquaredByWell(const std::vector<double>& response,
                        const std::vector<int>& wells,
                        const std::vector<std::size_t>& sigmaRows) {
    std::vector<double> vals;
    std::vector<int> grp;
    for (std::size_t row : sigmaRows) {
        if (row >= response.size() || row >= wells.size()) continue;
        if (!finite(response[row]) || wells[row] < 0) continue;
        vals.push_back(response[row]);
        grp.push_back(wells[row]);
    }
    if (vals.size() < 4) return kNaN;
    const double grand = bmstats::mean(vals);
    std::array<double, 3> sum = {0, 0, 0};
    std::array<int, 3> cnt = {0, 0, 0};
    for (std::size_t i = 0; i < vals.size(); ++i) {
        sum[static_cast<std::size_t>(grp[i])] += vals[i];
        ++cnt[static_cast<std::size_t>(grp[i])];
    }
    int nonEmpty = 0;
    for (int c : cnt) nonEmpty += (c > 0) ? 1 : 0;
    if (nonEmpty < 2) return kNaN;  // need >=2 occupied wells for between-well SS.
    double ssBetween = 0.0;
    double ssTotal = 0.0;
    for (int w = 0; w < 3; ++w) {
        if (cnt[static_cast<std::size_t>(w)] == 0) continue;
        const double mw = sum[static_cast<std::size_t>(w)] / cnt[static_cast<std::size_t>(w)];
        ssBetween += cnt[static_cast<std::size_t>(w)] * (mw - grand) * (mw - grand);
    }
    for (double v : vals) ssTotal += (v - grand) * (v - grand);
    if (!(ssTotal > 0.0)) return kNaN;
    return ssBetween / ssTotal;
}

// Partial r: corr(sigma, driver | mediator) via residualised series.
// Cross-check (§4.6): Arg43 CZ 91% mediation-collapse (raw r ~0.50 -> partial ~0.045).
double partialR(const std::vector<double>& y,
                const std::vector<double>& x,
                const std::vector<double>& m) {
    if (y.size() != x.size() || y.size() != m.size() || y.size() < 3) return kNaN;
    const OlsResult yOnM = ols(m, y);
    const OlsResult xOnM = ols(m, x);
    if (!finite(yOnM.slope) || !finite(xOnM.slope)) return kNaN;
    std::vector<double> yres;
    std::vector<double> xres;
    yres.reserve(y.size());
    xres.reserve(x.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        yres.push_back(y[i] - (yOnM.intercept + yOnM.slope * m[i]));
        xres.push_back(x[i] - (xOnM.intercept + xOnM.slope * m[i]));
    }
    return pearsonR(xres, yres);
}

// Block-bootstrap CI for partial r: resample contiguous blocks (block =
// kNullBlock), recompute partial r, take the 2.5/97.5 percentiles. NOT
// n_eff-deflated (progression-not-mean). The honest small-n CI (§4.6).
std::pair<double, double> blockBootstrapPartialR(const std::vector<double>& y,
                                                 const std::vector<double>& x,
                                                 const std::vector<double>& m) {
    if (y.size() != x.size() || y.size() != m.size() || y.size() < kNullBlock * 2)
        return {kNaN, kNaN};
    std::mt19937 rng(kNullSeed);
    const int n = static_cast<int>(y.size());
    const int nBlocks = (n + kNullBlock - 1) / kNullBlock;
    std::uniform_int_distribution<int> startDist(0, std::max(0, n - kNullBlock));
    constexpr int kBootstrap = 2000;
    std::vector<double> estimates;
    estimates.reserve(kBootstrap);
    for (int b = 0; b < kBootstrap; ++b) {
        std::vector<double> by;
        std::vector<double> bx;
        std::vector<double> bm;
        by.reserve(static_cast<std::size_t>(n));
        bx.reserve(static_cast<std::size_t>(n));
        bm.reserve(static_cast<std::size_t>(n));
        for (int blk = 0; blk < nBlocks; ++blk) {
            const int start = startDist(rng);
            for (int j = 0; j < kNullBlock && static_cast<int>(by.size()) < n; ++j) {
                const int idx = start + j;
                if (idx >= n) break;
                by.push_back(y[static_cast<std::size_t>(idx)]);
                bx.push_back(x[static_cast<std::size_t>(idx)]);
                bm.push_back(m[static_cast<std::size_t>(idx)]);
            }
        }
        const double pr = partialR(by, bx, bm);
        if (finite(pr)) estimates.push_back(pr);
    }
    if (estimates.size() < 20) return {kNaN, kNaN};
    std::sort(estimates.begin(), estimates.end());
    auto pct = [&](double q) {
        const double pos = q * static_cast<double>(estimates.size() - 1);
        const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
        const std::size_t hi = std::min(lo + 1, estimates.size() - 1);
        const double frac = pos - static_cast<double>(lo);
        return estimates[lo] * (1.0 - frac) + estimates[hi] * frac;
    };
    return {pct(0.025), pct(0.975)};
}

std::pair<double, double> olsSlopeCi95(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const OlsResult& fit) {
    if (x.size() != y.size() || x.size() < 3 || !finite(fit.slope) || !finite(fit.intercept))
        return {kNaN, kNaN};
    const double mx = bmstats::mean(x);
    double sxx = 0.0;
    double rss = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - mx;
        sxx += dx * dx;
        const double r = y[i] - (fit.intercept + fit.slope * x[i]);
        rss += r * r;
    }
    if (!(sxx > 0.0)) return {kNaN, kNaN};
    const double sigma2 = rss / static_cast<double>(x.size() - 2);
    const double se = std::sqrt(sigma2 / sxx);
    if (!finite(se)) return {kNaN, kNaN};
    constexpr double z95 = 1.959963984540054;
    return {fit.slope - z95 * se, fit.slope + z95 * se};
}

QString supportClass(std::size_t finiteN) {
    if (finiteN == 0) return QStringLiteral("empty");
    if (finiteN == 1) return QStringLiteral("singleton");
    if (finiteN < 10) return QStringLiteral("very_low");
    if (finiteN < 30) return QStringLiteral("low");
    return QStringLiteral("standard");
}

bool summaryPresent(const DistributionSummary& s) {
    return s.finite_n > 0;
}

void appendDistributionHeader(QStringList& cols, const QString& prefix) {
    cols << prefix + QStringLiteral("_n")
         << prefix + QStringLiteral("_finite_n")
         << prefix + QStringLiteral("_finite_frac")
         << prefix + QStringLiteral("_mean")
         << prefix + QStringLiteral("_sd")
         << prefix + QStringLiteral("_min")
         << prefix + QStringLiteral("_p05")
         << prefix + QStringLiteral("_p25")
         << prefix + QStringLiteral("_median")
         << prefix + QStringLiteral("_p75")
         << prefix + QStringLiteral("_p95")
         << prefix + QStringLiteral("_max");
}

void appendDistributionCells(QStringList& cells, const DistributionSummary& s) {
    cells << QString::number(static_cast<qulonglong>(s.n))
          << QString::number(static_cast<qulonglong>(s.finite_n))
          << csvNum(s.finite_frac)
          << csvNum(summaryPresent(s) ? s.mean : kNaN)
          << csvNum(summaryPresent(s) ? s.sd : kNaN)
          << csvNum(summaryPresent(s) ? s.min : kNaN)
          << csvNum(summaryPresent(s) ? s.p05 : kNaN)
          << csvNum(summaryPresent(s) ? s.p25 : kNaN)
          << csvNum(summaryPresent(s) ? s.median : kNaN)
          << csvNum(summaryPresent(s) ? s.p75 : kNaN)
          << csvNum(summaryPresent(s) ? s.p95 : kNaN)
          << csvNum(summaryPresent(s) ? s.max : kNaN);
}

QStringList classicalSourceHeaderColumns() {
    QStringList cols = {
        QStringLiteral("schema_version"),
        QStringLiteral("granularity"),
        QStringLiteral("dataset_id"),
        QStringLiteral("protein_id"),
        QStringLiteral("atom_uid"),
        QStringLiteral("atom_index"),
        QStringLiteral("residue_index"),
        QStringLiteral("residue_number"),
        QStringLiteral("residue_type"),
        QStringLiteral("atom_name"),
        QStringLiteral("element"),
        QStringLiteral("backbone_role"),
        QStringLiteral("frame_kind"),
        QStringLiteral("frame_kind_ord"),
        QStringLiteral("n_instances"),
        QStringLiteral("finite_n"),
        QStringLiteral("finite_frac"),
        QStringLiteral("support_class"),
        QStringLiteral("singleton_flag"),
        QStringLiteral("missing_n"),
        QStringLiteral("missing_reason"),
        QStringLiteral("sigma_cl_formula"),
        QStringLiteral("formula_terms_present"),
        QStringLiteral("sigma0_term_present"),
        QStringLiteral("sigma0_constant_status"),
        QStringLiteral("buckingham_linear_term_present"),
        QStringLiteral("buckingham_linear_constant_status"),
        QStringLiteral("buckingham_linear_constant_key"),
        QStringLiteral("buckingham_linear_constant_value"),
        QStringLiteral("buckingham_linear_constant_units"),
        QStringLiteral("buckingham_quadratic_term_present"),
        QStringLiteral("buckingham_quadratic_constant_status"),
        QStringLiteral("buckingham_quadratic_constant_key"),
        QStringLiteral("buckingham_quadratic_constant_value"),
        QStringLiteral("buckingham_quadratic_constant_units"),
        QStringLiteral("ring_term_present"),
        QStringLiteral("ring_constant_status"),
        QStringLiteral("mcconnell_term_present"),
        QStringLiteral("mcconnell_constant_status"),
        QStringLiteral("larsen_term_present"),
        QStringLiteral("larsen_constant_status"),
    };
    for (const QString& prefix : {
             QStringLiteral("sigma_qm"),
             QStringLiteral("sigma0"),
             QStringLiteral("contrib_buckingham_linear"),
             QStringLiteral("contrib_buckingham_quadratic"),
             QStringLiteral("contrib_ring"),
             QStringLiteral("contrib_mcconnell"),
             QStringLiteral("contrib_larsen"),
             QStringLiteral("sigma_cl"),
             QStringLiteral("residual"),
         }) {
        appendDistributionHeader(cols, prefix);
    }
    cols << QStringLiteral("slope_cl_qm")
         << QStringLiteral("slope_ci_low")
         << QStringLiteral("slope_ci_high")
         << QStringLiteral("r_cl_qm")
         << QStringLiteral("rmsd_ppm")
         << QStringLiteral("scale_factor")
         << QStringLiteral("tracking_label")
         << QStringLiteral("constant_placeholder_n");
    return cols;
}

QString termsPresentString(bool sigma0Present,
                           bool buckLinearPresent,
                           bool buckQuadraticPresent,
                           bool ringPresent,
                           bool mcPresent,
                           bool larsenPresent) {
    QStringList terms;
    if (sigma0Present) terms << QStringLiteral("sigma0");
    if (buckLinearPresent) terms << QStringLiteral("buckingham_linear");
    if (buckQuadraticPresent) terms << QStringLiteral("buckingham_quadratic");
    if (ringPresent) terms << QStringLiteral("ring");
    if (mcPresent) terms << QStringLiteral("mcconnell");
    if (larsenPresent) terms << QStringLiteral("larsen");
    return terms.join(QLatin1Char('+'));
}

double rmsd(const std::vector<double>& residuals) {
    double ss = 0.0;
    std::size_t n = 0;
    for (double v : residuals) {
        if (!finite(v)) continue;
        ss += v * v;
        ++n;
    }
    return n > 0 ? std::sqrt(ss / static_cast<double>(n)) : kNaN;
}

ClassicalSourceTermRecord MergeClassicalSourceRecords(
        const ClassicalSourceTermRecord& seed,
        const std::vector<ClassicalSourceTermRecord>& records) {
    ClassicalSourceTermRecord out = seed;
    out.schema_version = QStringLiteral("classical_source_terms_v2");
    out.granularity = QStringLiteral("residue_iupac");
    out.atom_uid.clear();
    out.atom_index = 0;
    out.residue_index = -1;
    out.residue_number = 0;
    out.backbone_role = seed.backbone_role;
    out.frame_kind = seed.frame_kind;
    out.frame_kind_ord = seed.frame_kind_ord;
    out.sigma_qm.clear();
    out.sigma0.clear();
    out.buckingham_linear.clear();
    out.buckingham_quadratic.clear();
    out.ring.clear();
    out.mcconnell.clear();
    out.larsen.clear();
    out.sigma_cl.clear();
    out.residual.clear();
    out.constant_placeholder_n = 0;

    auto append = [](std::vector<double>& dst, const std::vector<double>& src) {
        dst.insert(dst.end(), src.begin(), src.end());
    };
    QStringList sigma0Statuses;
    QStringList buckLinStatuses;
    QStringList buckQuadStatuses;
    QStringList ringStatuses;
    QStringList mcStatuses;
    QStringList larsenStatuses;
    auto addStatus = [](QStringList& statuses, const QString& status) {
        if (!status.isEmpty() && !statuses.contains(status)) statuses.append(status);
    };

    for (const ClassicalSourceTermRecord& r : records) {
        append(out.sigma_qm, r.sigma_qm);
        append(out.sigma0, r.sigma0);
        append(out.buckingham_linear, r.buckingham_linear);
        append(out.buckingham_quadratic, r.buckingham_quadratic);
        append(out.ring, r.ring);
        append(out.mcconnell, r.mcconnell);
        append(out.larsen, r.larsen);
        append(out.sigma_cl, r.sigma_cl);
        append(out.residual, r.residual);
        out.constant_placeholder_n += r.constant_placeholder_n;
        addStatus(sigma0Statuses, r.sigma0_status);
        addStatus(buckLinStatuses, r.buckingham_linear_status);
        addStatus(buckQuadStatuses, r.buckingham_quadratic_status);
        addStatus(ringStatuses, r.ring_status);
        addStatus(mcStatuses, r.mcconnell_status);
        addStatus(larsenStatuses, r.larsen_status);
    }
    out.sigma0_status = sigma0Statuses.join(QLatin1Char(';'));
    out.buckingham_linear_status = buckLinStatuses.join(QLatin1Char(';'));
    out.buckingham_quadratic_status = buckQuadStatuses.join(QLatin1Char(';'));
    out.ring_status = ringStatuses.join(QLatin1Char(';'));
    out.mcconnell_status = mcStatuses.join(QLatin1Char(';'));
    out.larsen_status = larsenStatuses.join(QLatin1Char(';'));
    return out;
}

void writeClassicalSourceRecord(QTextStream& out, const ClassicalSourceTermRecord& r) {
    const DistributionSummary sigmaQm = SummarizeDistribution(r.sigma_qm);
    const DistributionSummary sigma0 = SummarizeDistribution(r.sigma0);
    const DistributionSummary buckLinear = SummarizeDistribution(r.buckingham_linear);
    const DistributionSummary buckQuadratic = SummarizeDistribution(r.buckingham_quadratic);
    const DistributionSummary ring = SummarizeDistribution(r.ring);
    const DistributionSummary mc = SummarizeDistribution(r.mcconnell);
    const DistributionSummary larsen = SummarizeDistribution(r.larsen);
    const DistributionSummary sigmaCl = SummarizeDistribution(r.sigma_cl);
    const DistributionSummary residual = SummarizeDistribution(r.residual);

    std::vector<double> x;
    std::vector<double> y;
    const std::size_t n = std::min(r.sigma_qm.size(), r.sigma_cl.size());
    x.reserve(n);
    y.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!finite(r.sigma_qm[i]) || !finite(r.sigma_cl[i])) continue;
        x.push_back(r.sigma_cl[i]);
        y.push_back(r.sigma_qm[i]);
    }
    const OlsResult fit = ols(x, y);
    const auto [slopeLo, slopeHi] = olsSlopeCi95(x, y, fit);
    const double corr = pearsonR(x, y);
    const double scale = SdRatioScaleFactor(y, x);
    const double rmsdPpm = rmsd(r.residual);
    const QString trackingLabel =
        finite(corr) && std::abs(corr) >= 0.7 ? QStringLiteral("tracking")
        : finite(corr) && std::abs(corr) >= 0.3 ? QStringLiteral("partial_tracking")
        : finite(corr) ? QStringLiteral("weak_or_flat")
                       : QStringLiteral("insufficient_finite_pairs");

    const bool sigma0Present = summaryPresent(sigma0);
    const bool buckLinearPresent = summaryPresent(buckLinear);
    const bool buckQuadraticPresent = summaryPresent(buckQuadratic);
    const bool ringPresent = summaryPresent(ring);
    const bool mcPresent = summaryPresent(mc);
    const bool larsenPresent = summaryPresent(larsen);
    const QString termsPresent = termsPresentString(sigma0Present, buckLinearPresent,
                                                    buckQuadraticPresent, ringPresent,
                                                    mcPresent, larsenPresent);
    const std::size_t finiteN = x.size();
    const std::size_t missingN = n > finiteN ? n - finiteN : 0;
    const QString missingReason = missingN == 0 ? QString()
                                : QStringLiteral("nonfinite_sigma_qm_or_sigma_cl");
    QStringList cells;
    cells << csvEscape(r.schema_version)
          << csvEscape(r.granularity)
          << csvEscape(r.dataset_id)
          << csvEscape(r.protein_id)
          << csvEscape(r.atom_uid)
          << (r.granularity == QStringLiteral("atom_leaf")
                  ? QString::number(static_cast<qulonglong>(r.atom_index))
                  : QString())
          << QString::number(r.residue_index)
          << QString::number(r.residue_number)
          << csvEscape(r.residue_type)
          << csvEscape(r.atom_name)
          << csvEscape(r.element)
          << csvEscape(r.backbone_role)
          << csvEscape(r.frame_kind)
          << csvEscape(r.frame_kind_ord)
          << QString::number(static_cast<qulonglong>(n))
          << QString::number(static_cast<qulonglong>(finiteN))
          << csvNum(n > 0 ? static_cast<double>(finiteN) / static_cast<double>(n) : kNaN)
          << csvEscape(supportClass(finiteN))
          << csvBool(finiteN == 1)
          << QString::number(static_cast<qulonglong>(missingN))
          << csvEscape(missingReason)
          << csvEscape(QStringLiteral("sigma_cl = sigma0 + buckingham_linear + buckingham_quadratic + ring + mcconnell + larsen (present terms only)"))
          << csvEscape(termsPresent)
          << csvBool(sigma0Present)
          << csvEscape(r.sigma0_status)
          << csvBool(buckLinearPresent)
          << csvEscape(r.buckingham_linear_status)
          << csvEscape(r.buckingham_linear_key)
          << csvNum(r.buckingham_linear_constant_value)
          << csvEscape(r.buckingham_linear_units)
          << csvBool(buckQuadraticPresent)
          << csvEscape(r.buckingham_quadratic_status)
          << csvEscape(r.buckingham_quadratic_key)
          << csvNum(r.buckingham_quadratic_constant_value)
          << csvEscape(r.buckingham_quadratic_units)
          << csvBool(ringPresent)
          << csvEscape(r.ring_status)
          << csvBool(mcPresent)
          << csvEscape(r.mcconnell_status)
          << csvBool(larsenPresent)
          << csvEscape(r.larsen_status);
    appendDistributionCells(cells, sigmaQm);
    appendDistributionCells(cells, sigma0);
    appendDistributionCells(cells, buckLinear);
    appendDistributionCells(cells, buckQuadratic);
    appendDistributionCells(cells, ring);
    appendDistributionCells(cells, mc);
    appendDistributionCells(cells, larsen);
    appendDistributionCells(cells, sigmaCl);
    appendDistributionCells(cells, residual);
    cells << csvNum(fit.slope)
          << csvNum(slopeLo)
          << csvNum(slopeHi)
          << csvNum(corr)
          << csvNum(rmsdPpm)
          << csvNum(scale)
          << csvEscape(trackingLabel)
          << QString::number(r.constant_placeholder_n);
    out << cells.join(QLatin1Char(',')) << '\n';
}

struct LeadLagResult {
    int best_lag = 0;
    double lead_r = kNaN;
};

// Lagged cross-correlation: driver(t) vs sigma(t+k) over |k| <= kLagWindow.
// best_lag is the k maximising |r|; positive k = driver leads sigma (§4.5).
// Both gathered on the sigma rows (the driver is read at the sigma cadence).
LeadLagResult leadLag(const std::vector<double>& x, const std::vector<double>& y) {
    LeadLagResult out;
    if (x.size() != y.size() || x.size() < static_cast<std::size_t>(2 * kLagWindow + 2))
        return out;
    const int n = static_cast<int>(x.size());
    double bestAbs = -1.0;
    for (int k = -kLagWindow; k <= kLagWindow; ++k) {
        std::vector<double> dx;
        std::vector<double> sy;
        for (int i = 0; i < n; ++i) {
            const int j = i + k;  // sigma at t+k
            if (j < 0 || j >= n) continue;
            dx.push_back(x[static_cast<std::size_t>(i)]);
            sy.push_back(y[static_cast<std::size_t>(j)]);
        }
        const double r = pearsonR(dx, sy);
        if (finite(r) && std::abs(r) > bestAbs) {
            bestAbs = std::abs(r);
            out.best_lag = k;
            out.lead_r = r;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Molecular-component extraction (ACCUMULATOR_RESPEC §4.7) + ContactedClass
// (§2.4). Concrete, named, version-pinned.
// ---------------------------------------------------------------------------

// The 6 symmetrized components of a (generally asymmetric) molecular-frame
// tensor, in the fixed order/naming the molecular-frame reads use (§4.7.1):
// {xx, yy, xy, xz, yz, zz} of S = 1/2 (M + M^T).
enum class MolComp { xx = 0, yy, xy, xz, yz, zz };
constexpr std::array<const char*, 6> kMolCompNames = {"xx", "yy", "xy", "xz", "yz", "zz"};

std::array<double, 6> symmetricComponents(const Mat3& m) {
    const Mat3 s = 0.5 * (m + m.transpose());
    return {s(0, 0), s(1, 1), s(0, 1), s(0, 2), s(1, 2), s(2, 2)};
}

double antisymmetricNorm(const Mat3& m) {
    const Mat3 a = 0.5 * (m - m.transpose());
    return tensorFrobenius(a);
}

struct TensorInvariants {
    double iso = kNaN;
    double span = kNaN;     // Omega = s33 - s11
    double aniso = kNaN;    // ||T2|| of the symmetric-traceless part (the size axis)
    double eta_H = kNaN;    // proper distance-from-isotropic asymmetry in [0,1]
    double skew = kNaN;     // Maryland skew 3(s22-iso)/Omega
    double frobenius = kNaN;
};

// Invariants from the symmetric part of a tensor, by the SAME proper
// distance-from-isotropic eta_H discipline as foldCsa (§4.7.2, the D1 fix).
TensorInvariants tensorInvariants(const Mat3& raw) {
    TensorInvariants inv;
    const Mat3 sym = 0.5 * (raw + raw.transpose());
    inv.frobenius = tensorFrobenius(raw);
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym);
    if (solver.info() != Eigen::Success) return inv;
    const Vec3 ev = solver.eigenvalues();  // ascending: s11 <= s22 <= s33
    if (!ev.allFinite()) return inv;
    const double s11 = ev[0];
    const double s22 = ev[1];
    const double s33 = ev[2];
    inv.iso = (s11 + s22 + s33) / 3.0;
    inv.span = s33 - s11;
    // aniso = ||symmetric-traceless T2||
    const Mat3 trless = sym - inv.iso * Mat3::Identity();
    inv.aniso = tensorFrobenius(trless);
    if (std::abs(inv.span) > 1e-12) inv.skew = 3.0 * (s22 - inv.iso) / inv.span;
    // proper eta_H: order by distance from isotropic (zz farthest), xx the next
    // farthest of the remaining, yy the closest.
    int zz = 0;
    double farthest = -1.0;
    for (int i = 0; i < 3; ++i) {
        const double d = std::abs(ev[i] - inv.iso);
        if (d > farthest) { farthest = d; zz = i; }
    }
    std::array<int, 2> rem{};
    int k = 0;
    for (int i = 0; i < 3; ++i) if (i != zz) rem[static_cast<std::size_t>(k++)] = i;
    int xx = rem[0];
    int yy = rem[1];
    if (std::abs(ev[xx] - inv.iso) < std::abs(ev[yy] - inv.iso)) std::swap(xx, yy);
    const double denom = ev[zz] - inv.iso;
    if (std::abs(denom) > 1e-12) {
        double eta = (ev[yy] - ev[xx]) / denom;
        if (eta < 0.0 && eta > -1e-9) eta = 0.0;
        if (eta > 1.0 && eta < 1.0 + 1e-9) eta = 1.0;
        if (eta >= 0.0 && eta <= 1.0) inv.eta_H = eta;
    }
    return inv;
}

struct EfgEigen {
    double v_zz_abs = kNaN;  // |V_zz|, NMR convention |V_zz| >= |V_yy| >= |V_xx|
    double eta = kNaN;       // (V_xx - V_yy)/V_zz, the internal asymmetry
};

// EFG eigenframe from a library-basis T2 5-vector: reconstruct the symmetric
// 3x3, eigendecompose, order by |eigenvalue| (NMR convention), eta from the
// two smaller-magnitude principal values (§4.7.2).
EfgEigen efgEigen(const std::array<double, 5>& t2) {
    EfgEigen out;
    const Mat3 m = ReconstructLibraryT2Matrix(0.0, t2);
    Eigen::SelfAdjointEigenSolver<Mat3> solver(m);
    if (solver.info() != Eigen::Success) return out;
    Vec3 ev = solver.eigenvalues();
    if (!ev.allFinite()) return out;
    std::array<double, 3> vals = {ev[0], ev[1], ev[2]};
    std::sort(vals.begin(), vals.end(),
              [](double a, double b) { return std::abs(a) < std::abs(b); });
    const double vxx = vals[0];
    const double vyy = vals[1];
    const double vzz = vals[2];
    out.v_zz_abs = std::abs(vzz);
    if (std::abs(vzz) > 1e-12) out.eta = (vxx - vyy) / vzz;
    return out;
}

// ContactedClass — the six-value CHARGE/AROMATICITY enum for the context key
// (§2.4 / D-S3). NONE is the through-space-isolated baseline (no contact).
enum class ContactedClass { NONE = 0, ACIDIC, BASIC, AROMATIC, POLAR, APOLAR };

constexpr std::array<ContactedClass, 5> kPresentContactClasses = {
    ContactedClass::ACIDIC,
    ContactedClass::BASIC,
    ContactedClass::AROMATIC,
    ContactedClass::POLAR,
    ContactedClass::APOLAR,
};

ContactedClass classOf(model::AminoAcid aa) {
    switch (aa) {
    case model::AminoAcid::ASP:
    case model::AminoAcid::GLU:
        return ContactedClass::ACIDIC;
    case model::AminoAcid::ARG:
    case model::AminoAcid::LYS:
    case model::AminoAcid::HIS:  // titratable salt-bridge/pi partner -> BASIC (§2.4)
        return ContactedClass::BASIC;
    case model::AminoAcid::PHE:
    case model::AminoAcid::TYR:  // AROMATIC for its ring (the TYR default, §2.4)
    case model::AminoAcid::TRP:
        return ContactedClass::AROMATIC;
    case model::AminoAcid::ASN:
    case model::AminoAcid::GLN:
    case model::AminoAcid::SER:
    case model::AminoAcid::THR:
    case model::AminoAcid::CYS:
    case model::AminoAcid::GLY:
        return ContactedClass::POLAR;
    case model::AminoAcid::ALA:
    case model::AminoAcid::VAL:
    case model::AminoAcid::LEU:
    case model::AminoAcid::ILE:
    case model::AminoAcid::MET:
    case model::AminoAcid::PRO:
        return ContactedClass::APOLAR;
    case model::AminoAcid::Unknown:
        return ContactedClass::NONE;
    }
    return ContactedClass::NONE;
}

QString contactedClassName(ContactedClass c) {
    switch (c) {
    case ContactedClass::NONE: return QStringLiteral("NONE");
    case ContactedClass::ACIDIC: return QStringLiteral("ACIDIC");
    case ContactedClass::BASIC: return QStringLiteral("BASIC");
    case ContactedClass::AROMATIC: return QStringLiteral("AROMATIC");
    case ContactedClass::POLAR: return QStringLiteral("POLAR");
    case ContactedClass::APOLAR: return QStringLiteral("APOLAR");
    }
    return QStringLiteral("NONE");
}

std::uint8_t contactedClassBit(ContactedClass c) {
    const int ord = static_cast<int>(c);
    if (ord <= static_cast<int>(ContactedClass::NONE)) return 0;
    return static_cast<std::uint8_t>(1u << (ord - 1));
}

QJsonArray contactedClassMaskOrdsJson(std::uint8_t mask) {
    QJsonArray a;
    for (ContactedClass c : kPresentContactClasses) {
        if ((mask & contactedClassBit(c)) != 0)
            a.append(static_cast<int>(c));
    }
    return a;
}

QJsonArray contactedClassMaskNamesJson(std::uint8_t mask) {
    QJsonArray a;
    for (ContactedClass c : kPresentContactClasses) {
        if ((mask & contactedClassBit(c)) != 0)
            a.append(contactedClassName(c));
    }
    return a;
}

QString contactedClassMaskName(std::uint8_t mask) {
    QStringList names;
    for (ContactedClass c : kPresentContactClasses) {
        if ((mask & contactedClassBit(c)) != 0)
            names.append(contactedClassName(c));
    }
    return names.isEmpty() ? QStringLiteral("NONE") : names.join(QStringLiteral("+"));
}

int contactedClassSingleOrd(std::uint8_t mask) {
    int ord = static_cast<int>(ContactedClass::NONE);
    for (ContactedClass c : kPresentContactClasses) {
        if ((mask & contactedClassBit(c)) == 0) continue;
        if (ord != static_cast<int>(ContactedClass::NONE)) return -1;
        ord = static_cast<int>(c);
    }
    return ord;
}

// ---------------------------------------------------------------------------
// THE CHARACTERIZATION CATALOGUE (ACCUMULATOR_RESPEC §2.2A).
// The fixed, version-pinned set
// of named (driver-family -> sigma-target -> frame) relationships any context
// may instantiate. The catalogue is the ONLY door: every emitted responses[]
// entry carries a CatalogueId from this enum -> the row-mill is structurally
// impossible (there is no slot for an arbitrary (driver, sigma-target) pair).
// 16 rows; each gated by a membership predicate evaluated against the atom.
// ---------------------------------------------------------------------------
enum class CatalogueId {
    FIELD_ISO = 0,
    FIELD_SPAN,
    FIELD_MOLCOMP,
    FIELD_SIGNED_BOND,
    EFG_ISO,
    EFG_MOLCOMP,
    EFG_EIGEN,
    RING_HEAVY,
    RING_H,
    MCCONNELL,
    SCHAR_ISO,
    WIBERG,
    BONDLEN_ISO,
    DIHEDRAL,
    HBOND,
    STRUCT,
};

// The right-frame tag carried on each response (§2.1 frame).
enum class FrameTag { Invariant = 0, Molecular, SignedBond, Scalar, Circular };

QString frameTagName(FrameTag f) {
    switch (f) {
    case FrameTag::Invariant: return QStringLiteral("invariant");
    case FrameTag::Molecular: return QStringLiteral("molecular");
    case FrameTag::SignedBond: return QStringLiteral("signed-bond");
    case FrameTag::Scalar: return QStringLiteral("scalar");
    case FrameTag::Circular: return QStringLiteral("circular");
    }
    return QStringLiteral("scalar");
}

struct CatalogueRow {
    CatalogueId id;
    const char* id_name;
    const char* driver_family;  // the mechanism family tag (category, not a filter)
    FrameTag frame;
};

static constexpr std::array<CatalogueRow, 16> kCharacterizationCatalogue = {{
    {CatalogueId::FIELD_ISO, "FIELD_ISO", "Electric-field", FrameTag::Invariant},
    {CatalogueId::FIELD_SPAN, "FIELD_SPAN", "Electric-field", FrameTag::Invariant},
    {CatalogueId::FIELD_MOLCOMP, "FIELD_MOLCOMP", "Electric-field", FrameTag::Molecular},
    {CatalogueId::FIELD_SIGNED_BOND, "FIELD_SIGNED_BOND", "Electric-field", FrameTag::SignedBond},
    {CatalogueId::EFG_ISO, "EFG_ISO", "EFG", FrameTag::Invariant},
    {CatalogueId::EFG_MOLCOMP, "EFG_MOLCOMP", "EFG", FrameTag::Molecular},
    {CatalogueId::EFG_EIGEN, "EFG_EIGEN", "EFG", FrameTag::Molecular},
    {CatalogueId::RING_HEAVY, "RING_HEAVY", "Ring-current", FrameTag::Invariant},
    {CatalogueId::RING_H, "RING_H", "Ring-current", FrameTag::Molecular},
    {CatalogueId::MCCONNELL, "MCCONNELL", "McConnell", FrameTag::Molecular},
    {CatalogueId::SCHAR_ISO, "SCHAR_ISO", "Local-electronic", FrameTag::Scalar},
    {CatalogueId::WIBERG, "WIBERG", "Local-electronic", FrameTag::Scalar},
    {CatalogueId::BONDLEN_ISO, "BONDLEN_ISO", "Local-electronic", FrameTag::Scalar},
    {CatalogueId::DIHEDRAL, "DIHEDRAL", "Geometry-dihedral", FrameTag::Circular},
    {CatalogueId::HBOND, "HBOND", "Structural-environment", FrameTag::Scalar},
    {CatalogueId::STRUCT, "STRUCT", "Structural-environment", FrameTag::Scalar},
}};

const CatalogueRow& catalogueRow(CatalogueId id) {
    return kCharacterizationCatalogue[static_cast<std::size_t>(id)];
}

QJsonObject couplingRecord(const QString& sigmaComponent,
                           const std::vector<double>& sigma,
                           const RunningSeriesRef& series,
                           const std::vector<std::size_t>& sigmaRows) {
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(sigmaRows.size());
    y.reserve(sigmaRows.size());
    for (std::size_t row : sigmaRows) {
        if (row >= sigma.size() || row >= series.values.size()) continue;
        if (!finite(sigma[row]) || !finite(series.values[row])) continue;
        x.push_back(series.values[row]);
        y.push_back(sigma[row]);
    }

    QJsonObject o;
    o.insert(QStringLiteral("sigma_component"), sigmaComponent);
    o.insert(QStringLiteral("series_ref"), series.series_ref);
    o.insert(QStringLiteral("component_ref"), series.component_ref);
    o.insert(QStringLiteral("relationship_index"), series.relationship_index);
    o.insert(QStringLiteral("relationship_facet"), series.relationship_facet);
    o.insert(QStringLiteral("n"), static_cast<int>(x.size()));
    o.insert(QStringLiteral("mean_facet"), x.empty() ? QJsonValue(QJsonValue::Null) : jd(bmstats::mean(x)));
    o.insert(QStringLiteral("var_facet"), x.empty() ? QJsonValue(QJsonValue::Null)
                                                     : jd(x.size() > 1 ? bmstats::variance(x) : 0.0));
    double covariance = kNaN;
    if (x.size() > 1) {
        try { covariance = bmstats::covariance(x, y); } catch (...) {}
    }
    const OlsResult fit = ols(x, y);
    o.insert(QStringLiteral("covariance"), jd(covariance));
    o.insert(QStringLiteral("ols_slope"), jd(fit.slope));
    o.insert(QStringLiteral("ols_intercept"), jd(fit.intercept));
    if (x.size() >= static_cast<std::size_t>(kRollingWindow)) {
        const auto xb = x.end() - kRollingWindow;
        const auto yb = y.end() - kRollingWindow;
        std::vector<double> rx(xb, x.end());
        std::vector<double> ry(yb, y.end());
        const OlsResult rolling = ols(rx, ry);
        o.insert(QStringLiteral("ols_slope_rolling"), jd(rolling.slope));
    } else {
        o.insert(QStringLiteral("ols_slope_rolling"), QJsonValue(QJsonValue::Null));
    }
    o.insert(QStringLiteral("chatterjee_xi"), jd(chatterjeeXi(x, y)));
    double deltaSlope = kNaN;
    if (x.size() >= 2) {
        std::vector<double> dx;
        std::vector<double> dy;
        dx.reserve(x.size() - 1);
        dy.reserve(y.size() - 1);
        for (std::size_t i = 1; i < x.size(); ++i) {
            dx.push_back(x[i] - x[i - 1]);
            dy.push_back(y[i] - y[i - 1]);
        }
        deltaSlope = ols(dx, dy).slope;
    }
    o.insert(QStringLiteral("delta_slope"), jd(deltaSlope));
    auto [nullMean, nullSd, obsZ, circP] = nullShiftStats(x, y, fit.slope);
    (void)circP;
    o.insert(QStringLiteral("null_slope_mean"), jd(nullMean));
    o.insert(QStringLiteral("null_slope_sd"), jd(nullSd));
    o.insert(QStringLiteral("obs_slope_z"), jd(obsZ));
    o.insert(QStringLiteral("pchip_midpoint"), jd(pchipMidpoint(x, y)));
    o.insert(QStringLiteral("coverage"), coverageJson(x));
    QJsonObject lev;
    lev.insert(QStringLiteral("top1_frame_slope_frac"), jd(leverageTop1(x, y)));
    o.insert(QStringLiteral("leverage"), lev);
    return o;
}

QJsonObject serialRecord(const RunningSeriesRef& series) {
    std::vector<double> x;
    x.reserve(series.values.size());
    for (double v : series.values)
        if (finite(v)) x.push_back(v);
    QJsonObject o;
    o.insert(QStringLiteral("series_ref"), series.series_ref);
    o.insert(QStringLiteral("component_ref"), series.component_ref);
    o.insert(QStringLiteral("relationship_index"), series.relationship_index);
    o.insert(QStringLiteral("relationship_facet"), series.relationship_facet);
    o.insert(QStringLiteral("n"), static_cast<int>(x.size()));
    const int lags = x.empty() ? 0 : std::min<int>(20, static_cast<int>(x.size() / 5));
    o.insert(QStringLiteral("ljung_box_lags"), lags);
    double ljung = kNaN;
    if (lags > 0 && x.size() > static_cast<std::size_t>(lags)) {
        try { ljung = bmstats::ljung_box(x, lags).first; } catch (...) {}
    }
    o.insert(QStringLiteral("ljung_box"), jd(ljung));
    double runs = kNaN;
    if (x.size() >= 2) {
        try { runs = bmstats::runs_above_and_below_median(x).first; } catch (...) {}
    }
    o.insert(QStringLiteral("runs_test"), jd(runs));
    double snr = kNaN;
    if (x.size() >= 4) {
        try { snr = bmstats::m2m4_snr_estimator(x); } catch (...) {}
    }
    o.insert(QStringLiteral("snr"), jd(snr));
    double period = kNaN;
    double amplitude = kNaN;
    if (x.size() >= 4) {
        try {
            const double mean = bmstats::mean(x);
            const std::size_t n = x.size();
            double bestAmp = 0.0;
            std::size_t bestK = 0;
            constexpr double pi = 3.14159265358979323846264338327950288;
            for (std::size_t k = 1; k <= n / 2; ++k) {
                double re = 0.0;
                double im = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    const double angle = 2.0 * pi * static_cast<double>(k * i) / static_cast<double>(n);
                    re += (x[i] - mean) * std::cos(angle);
                    im -= (x[i] - mean) * std::sin(angle);
                }
                const double amp = 2.0 * std::sqrt(re * re + im * im) / static_cast<double>(n);
                if (amp > bestAmp) {
                    bestAmp = amp;
                    bestK = k;
                }
            }
            if (bestK > 0 && bestAmp > 0.0) {
                period = static_cast<double>(n) / static_cast<double>(bestK);
                amplitude = bestAmp;
            }
        } catch (...) {}
    }
    o.insert(QStringLiteral("period_dominant"), jd(period));
    o.insert(QStringLiteral("period_amplitude"), jd(amplitude));
    return o;
}

void addScalarSeriesRef(std::vector<RunningSeriesRef>& refs, const QString& key, const ScalarSeries& s) {
    refs.push_back({key, QJsonValue(QJsonValue::Null), QJsonValue(QJsonValue::Null),
                    QJsonValue(QJsonValue::Null), componentSeries(s)});
}

void addVec3SeriesRefs(std::vector<RunningSeriesRef>& refs, const QString& key, const Vec3Series& s) {
    static constexpr std::array<const char*, 3> comps = {"x", "y", "z"};
    for (int c = 0; c < 3; ++c)
        refs.push_back({key, QString::fromLatin1(comps[static_cast<std::size_t>(c)]),
                        QJsonValue(QJsonValue::Null), QJsonValue(QJsonValue::Null),
                        componentSeries(s, c)});
}

QJsonObject boostJson(const std::vector<std::pair<QString, std::vector<double>>>& sigmas,
                      const std::vector<RunningSeriesRef>& serialRefs,
                      const std::vector<std::size_t>& sigmaRows,
                      bool coupling) {
    QJsonObject boost;
    QJsonArray couplingArray;
    if (coupling) {
        for (const auto& s : sigmas) {
            for (const RunningSeriesRef& ref : serialRefs) {
                couplingArray.append(couplingRecord(s.first, s.second, ref, sigmaRows));
            }
        }
    }
    QJsonArray serialArray;
    for (const RunningSeriesRef& ref : serialRefs) serialArray.append(serialRecord(ref));
    boost.insert(QStringLiteral("coupling"), couplingArray);
    boost.insert(QStringLiteral("serial"), serialArray);
    return boost;
}

bool validResidue(const model::QtProtein& protein, int32_t residue) {
    return residue >= 0 && static_cast<std::size_t>(residue) < protein.residueCount();
}

int residueNumber(const model::QtProtein& protein, int32_t residueIndex) {
    if (!validResidue(protein, residueIndex)) return -1;
    return protein.residue(static_cast<std::size_t>(residueIndex)).address.residueNumber;
}

model::AminoAcid residueAa(const model::QtProtein& protein, int32_t residueIndex) {
    if (!validResidue(protein, residueIndex)) return model::AminoAcid::Unknown;
    return protein.residue(static_cast<std::size_t>(residueIndex)).aminoAcid;
}

struct RelationshipKey {
    int hybridisation_ord = static_cast<int>(model::Hybridisation::Unassigned);
    int element_ord = static_cast<int>(model::Element::Unknown);
    int contacted_residue_number = -1;
    int contacted_amino_acid_ord = -1;
    int scope_ord = -1;
    int mechanism_ord = -1;
    int source_kind_ord = -1;
    int32_t source_id = -1;
    int32_t source_atom_index = -1;
    int source_category_ord = -1;

    bool operator<(const RelationshipKey& o) const {
        return std::tie(hybridisation_ord, element_ord, contacted_residue_number,
                        contacted_amino_acid_ord, scope_ord, mechanism_ord, source_kind_ord,
                        source_id, source_atom_index, source_category_ord)
               < std::tie(o.hybridisation_ord, o.element_ord, o.contacted_residue_number,
                          o.contacted_amino_acid_ord, o.scope_ord, o.mechanism_ord,
                          o.source_kind_ord, o.source_id, o.source_atom_index, o.source_category_ord);
    }
};

struct RelationshipSeries {
    explicit RelationshipSeries(std::size_t n = 0)
        : step_count(n) {}

    void reset(std::size_t n) {
        step_count = n;
        present_steps.clear();
        distance.samples.clear();
        cos_theta.samples.clear();
        inv_r3.samples.clear();
        dipolar.samples.clear();
        kernel_T0.samples.clear();
        contribution.samples.clear();
        kernel_T2.samples.clear();
        kernel_mol.samples.clear();
    }

    struct SparseScalar {
        void set(std::size_t step, double v) {
            if (!samples.empty() && samples.back().first == step) {
                samples.back().second = v;
                return;
            }
            samples.push_back({step, v});
        }
        std::vector<double> dense(std::size_t n) const {
            std::vector<double> out(n, kNaN);
            for (const auto& sample : samples)
                if (sample.first < n) out[sample.first] = sample.second;
            return out;
        }
        QJsonArray json(std::size_t n) const {
            QJsonArray a;
            std::size_t sample = 0;
            for (std::size_t step = 0; step < n; ++step) {
                while (sample < samples.size() && samples[sample].first < step) ++sample;
                if (sample < samples.size() && samples[sample].first == step)
                    a.append(jd(samples[sample].second));
                else
                    a.append(QJsonValue(QJsonValue::Null));
            }
            return a;
        }
        std::vector<std::pair<std::size_t, double>> samples;
    };

    struct SparseT2 {
        void set(std::size_t step, const std::array<double, 5>& v) {
            if (!std::any_of(v.begin(), v.end(), [](double x) { return finite(x); })) return;
            if (!samples.empty() && samples.back().first == step) {
                samples.back().second = v;
                return;
            }
            samples.push_back({step, v});
        }
        std::vector<double> componentDense(std::size_t n, int component) const {
            std::vector<double> out(n, kNaN);
            for (const auto& sample : samples)
                if (sample.first < n) out[sample.first] = sample.second[static_cast<std::size_t>(component)];
            return out;
        }
        QJsonArray json(std::size_t n) const {
            QJsonArray a;
            std::size_t sample = 0;
            for (std::size_t step = 0; step < n; ++step) {
                while (sample < samples.size() && samples[sample].first < step) ++sample;
                if (sample < samples.size() && samples[sample].first == step)
                    a.append(vec5Json(samples[sample].second));
                else
                    a.append(QJsonValue(QJsonValue::Null));
            }
            return a;
        }
        std::vector<std::pair<std::size_t, std::array<double, 5>>> samples;
    };

    struct SparseMat3 {
        void set(std::size_t step, const Mat3& v) {
            if (!v.allFinite()) return;
            if (!samples.empty() && samples.back().first == step) {
                samples.back().second = v;
                return;
            }
            samples.push_back({step, v});
        }
        std::vector<double> componentDense(std::size_t n, MolComp comp) const {
            std::vector<double> out(n, kNaN);
            const std::size_t idx = static_cast<std::size_t>(comp);
            for (const auto& sample : samples)
                if (sample.first < n)
                    out[sample.first] = symmetricComponents(sample.second)[idx];
            return out;
        }
        Mat3Series dense(std::size_t n) const {
            Mat3Series out(n);
            for (const auto& sample : samples)
                if (sample.first < n) out.set(sample.first, sample.second);
            return out;
        }
        std::vector<std::pair<std::size_t, Mat3>> samples;
    };

    void markPresent(std::size_t step, bool present) {
        if (!present) return;
        if (present_steps.empty() || present_steps.back() != step) present_steps.push_back(step);
    }

    QJsonArray presentJson() const {
        QJsonArray a;
        std::size_t sample = 0;
        for (std::size_t step = 0; step < step_count; ++step) {
            while (sample < present_steps.size() && present_steps[sample] < step) ++sample;
            a.append(sample < present_steps.size() && present_steps[sample] == step);
        }
        return a;
    }

    std::size_t step_count = 0;
    std::vector<std::size_t> present_steps;
    SparseScalar distance;
    SparseScalar cos_theta;
    SparseScalar inv_r3;
    SparseScalar dipolar;
    SparseScalar kernel_T0;
    SparseScalar contribution;
    SparseT2 kernel_T2;
    SparseMat3 kernel_mol;
};

QString relationshipScope(const Body& body, std::size_t target_atom, const PairContribution& pair) {
    const bool selfOrBonded = (pair.pointer_flags & SelfOrBondedFlag) != 0;
    if (!selfOrBonded) return QStringLiteral("through_space");
    if (!body.run.protein) return QStringLiteral("bonded");
    const model::QtProtein& protein = *body.run.protein;
    const bool self =
        pair.source_atom_index == static_cast<int32_t>(target_atom)
        || (pair.source_kind == QStringLiteral("ring_center") && pair.source_id >= 0
            && static_cast<std::size_t>(pair.source_id) < protein.topology().ringCount()
            && std::find(protein.topology().ringAt(static_cast<std::size_t>(pair.source_id)).atomIndices.begin(),
                         protein.topology().ringAt(static_cast<std::size_t>(pair.source_id)).atomIndices.end(),
                         static_cast<int32_t>(verbs::heavyParent(body, target_atom)))
                   != protein.topology().ringAt(static_cast<std::size_t>(pair.source_id)).atomIndices.end())
        || (pair.source_kind == QStringLiteral("bond_midpoint") && pair.source_id >= 0
            && static_cast<std::size_t>(pair.source_id) < protein.topology().bondCount()
            && (protein.topology().bondAt(static_cast<std::size_t>(pair.source_id)).atomIndexA
                    == static_cast<int32_t>(target_atom)
                || protein.topology().bondAt(static_cast<std::size_t>(pair.source_id)).atomIndexB
                       == static_cast<int32_t>(target_atom)));
    if (self) return QStringLiteral("self");
    if (pair.pointer_flags & NearFieldFlag) return QStringLiteral("bonded_near_field");
    return QStringLiteral("bonded");
}

std::pair<int, model::AminoAcid> contactedResidue(const Body& body, const PairContribution& pair) {
    if (!body.run.protein) return {-1, model::AminoAcid::Unknown};
    const model::QtProtein& p = *body.run.protein;
    if (pair.source_atom_index >= 0
        && static_cast<std::size_t>(pair.source_atom_index) < p.atomCount()) {
        const int32_t ri = p.atom(static_cast<std::size_t>(pair.source_atom_index)).residueIndex;
        return {residueNumber(p, ri), residueAa(p, ri)};
    }
    if (pair.source_kind == QStringLiteral("ring_center") && pair.source_id >= 0
        && static_cast<std::size_t>(pair.source_id) < p.topology().ringCount()) {
        const int32_t ri = p.topology().ringAt(static_cast<std::size_t>(pair.source_id)).parentResidueIndex;
        return {residueNumber(p, ri), residueAa(p, ri)};
    }
    if (pair.source_kind == QStringLiteral("bond_midpoint") && pair.source_id >= 0
        && static_cast<std::size_t>(pair.source_id) < p.topology().bondCount()) {
        const model::QtBond& b = p.topology().bondAt(static_cast<std::size_t>(pair.source_id));
        if (b.atomIndexA >= 0 && static_cast<std::size_t>(b.atomIndexA) < p.atomCount()) {
            const int32_t ri = p.atom(static_cast<std::size_t>(b.atomIndexA)).residueIndex;
            return {residueNumber(p, ri), residueAa(p, ri)};
        }
    }
    return {-1, model::AminoAcid::Unknown};
}

EfgValue readT2Field(const Body& body, io::FieldKind kind, ArrayId id, std::size_t atom, std::size_t row) {
    EfgValue out;
    if (body.catalog.present(body, id, atom, row)) {
        out.t2 = body.catalog.valueT2(body, id, atom, row);
    } else {
        for (int i = 0; i < 5; ++i) {
            if (const std::optional<double> v = body.catalog.value(body, kind, atom, row, i))
                out.t2[static_cast<std::size_t>(i)] = *v;
        }
    }
    return out;
}

}  // namespace

AnalysisCadence::AnalysisCadence(const FrameMap& frameMap) {
    original_by_step_.reserve(frameMap.frameCount());
    for (std::size_t row = 0; row < frameMap.frameCount(); ++row)
        original_by_step_.push_back(frameMap.originalIndex(row));
    sigma_present_.assign(frameMap.frameCount(), false);
    sigma_rows_ = frameMap.dftRows();
    for (std::size_t row : sigma_rows_)
        if (row < sigma_present_.size()) sigma_present_[row] = true;
    if (sigma_rows_.size() >= 2) {
        const std::size_t stride = sigma_rows_[1] - sigma_rows_[0];
        bool regular = stride > 0;
        for (std::size_t i = 2; i < sigma_rows_.size(); ++i)
            regular = regular && (sigma_rows_[i] - sigma_rows_[i - 1] == stride);
        if (regular) nominal_stride_ = stride;
    }
}

AnalysisObjectContext::AnalysisObjectContext(const Body& body, const AnalysisCadence& cadence)
    : body_(body), cadence_(cadence) {}

std::shared_ptr<const model::QtConformationSnapshot> AnalysisObjectContext::snapshot(std::size_t step) const {
    if (!body_.run.conformation) return nullptr;
    if (cached_step_ != step) rebuildMopacCache(step);
    return cached_snapshot_;
}

void AnalysisObjectContext::rebuildMopacCache(std::size_t step) const {
    cached_step_ = step;
    cached_snapshot_.reset();
    cached_mopac_bonds_.clear();
    if (!body_.run.conformation) return;
    body_.run.conformation->requestSnapshot(step);
    cached_snapshot_ = body_.run.conformation->snapshot(step);
    if (!cached_snapshot_ || !cached_snapshot_->has(io::FieldKind::MOPACBondOrders)) return;
    model::QtMopacCoreGroup mopac(*cached_snapshot_);
    cached_mopac_bonds_.reserve(mopac.bondOrderCount());
    for (std::size_t i = 0; i < mopac.bondOrderCount(); ++i) {
        const std::optional<model::MopacBondOrder> b = mopac.bondOrder(i);
        if (!b || !finite(b->wibergOrder)) continue;
        const auto p = canonicalPair(b->atomA, b->atomB);
        cached_mopac_bonds_.push_back({p.first, p.second, b->wibergOrder});
        std::vector<double>& series = mopac_wiberg_by_pair_[pairKey(p.first, p.second)];
        if (series.size() != cadence_.stepCount())
            series.assign(cadence_.stepCount(), kNaN);
        if (step < series.size()) series[step] = b->wibergOrder;
    }
}

const std::vector<MopacFrameBond>& AnalysisObjectContext::mopacBonds(std::size_t step) const {
    if (cached_step_ != step) rebuildMopacCache(step);
    return cached_mopac_bonds_;
}

std::vector<MopacFrameBond> AnalysisObjectContext::mopacBondsForAtom(std::size_t step, std::size_t atom) const {
    std::vector<MopacFrameBond> out;
    for (const MopacFrameBond& b : mopacBonds(step))
        if (b.atom_a == atom || b.atom_b == atom) out.push_back(b);
    return out;
}

std::optional<double> AnalysisObjectContext::mopacWiberg(std::size_t step,
                                                         std::size_t atomA,
                                                         std::size_t atomB) const {
    const auto p = canonicalPair(atomA, atomB);
    const auto sit = mopac_wiberg_by_pair_.find(pairKey(p.first, p.second));
    if (sit != mopac_wiberg_by_pair_.end()
        && step < sit->second.size()
        && finite(sit->second[step])) {
        return sit->second[step];
    }
    for (const MopacFrameBond& b : mopacBonds(step)) {
        if (b.atom_a == p.first && b.atom_b == p.second) return b.wiberg_order;
    }
    return std::nullopt;
}

std::vector<double> AnalysisObjectContext::mopacWibergSeries(std::size_t atomA,
                                                             std::size_t atomB) const {
    const auto p = canonicalPair(atomA, atomB);
    const auto it = mopac_wiberg_by_pair_.find(pairKey(p.first, p.second));
    if (it != mopac_wiberg_by_pair_.end()) return it->second;
    return std::vector<double>(cadence_.stepCount(), kNaN);
}

AnalysisElement::AnalysisElement(const AnalysisObjectContext& context,
                                 QString objectType,
                                 std::size_t modelIndex)
    : context_(context),
      body_(context.body()),
      cadence_(context.cadence()),
      object_type_(std::move(objectType)),
      model_index_(modelIndex) {}

AnalysisStructure::AnalysisStructure(const AnalysisObjectContext& context,
                                     QString objectType,
                                     std::size_t modelIndex)
    : AnalysisElement(context, std::move(objectType), modelIndex) {}

class AnalysisAtom::Impl {
public:
    Impl(const AnalysisObjectContext& context, std::size_t atom, PerAtomSubstrateConfig config)
        : context_(context),
          body_(context.body()),
          cadence_(context.cadence()),
          atom_(atom),
          config_(std::move(config)),
          coord_(cadence_.stepCount()),
          molecular_frame_(cadence_.stepCount()),
          sigma_total_(cadence_.stepCount()),
          sigma_dia_(cadence_.stepCount()),
          sigma_para_(cadence_.stepCount()),
          sigma_total_raw_(cadence_.stepCount()),
          sigma_dia_raw_(cadence_.stepCount()),
          sigma_para_raw_(cadence_.stepCount()),
          sigma_mol_total_(cadence_.stepCount()),
          sigma_mol_dia_(cadence_.stepCount()),
          sigma_mol_para_(cadence_.stepCount()),
          sigma_frob_(cadence_.stepCount()),
          sigma_csa_(cadence_.stepCount()),
          field_mopac_(cadence_.stepCount()),
          field_ff14sb_(cadence_.stepCount()),
          field_apbs_(cadence_.stepCount()),
          field_charge_ff14sb_(cadence_.stepCount()),
          field_mol_mopac_(cadence_.stepCount()),
          field_mol_apbs_(cadence_.stepCount()),
          field_mol_charge_ff14sb_(cadence_.stepCount()),
          field_E_z_mopac_(cadence_.stepCount()),
          field_E_z_apbs_(cadence_.stepCount()),
          field_E_z_charge_ff14sb_(cadence_.stepCount()),
          field_E2_mopac_(cadence_.stepCount()),
          field_E2_apbs_(cadence_.stepCount()),
          field_E2_charge_ff14sb_(cadence_.stepCount()),
          field_abs_E_mopac_(cadence_.stepCount()),
          field_abs_E_apbs_(cadence_.stepCount()),
          field_abs_E_charge_ff14sb_(cadence_.stepCount()),
          field_abs_E2_mopac_(cadence_.stepCount()),
          field_abs_E2_apbs_(cadence_.stepCount()),
          field_abs_E2_charge_ff14sb_(cadence_.stepCount()),
          charge_field_source_count_(cadence_.stepCount()),
          charge_field_rejected_self_(cadence_.stepCount()),
          charge_field_rejected_degenerate_(cadence_.stepCount()),
          charge_field_rejected_cutoff_(cadence_.stepCount()),
          charge_field_missing_charge_(cadence_.stepCount()),
          efg_apbs_(cadence_.stepCount()),
          efg_aimnet2_(cadence_.stepCount()),
          efg_mol_apbs_(cadence_.stepCount()),
          efg_mol_aimnet2_(cadence_.stepCount()),
          efg_abs_apbs_(cadence_.stepCount()),
          efg_abs_aimnet2_(cadence_.stepCount()),
          shielding_mopac_coulomb_(cadence_.stepCount()),
          shielding_abs_mopac_coulomb_(cadence_.stepCount()),
          shielding_mol_mopac_(cadence_.stepCount()),
          bs_per_type_T2_(cadence_.stepCount()),
          hm_per_type_T2_(cadence_.stepCount()),
          mc_tensor_series_(kMcTensorFields.size(), TensorSeries(cadence_.stepCount())),
          bs_per_type_mol_(8, Mat3Series(cadence_.stepCount())),
          hm_per_type_mol_(8, Mat3Series(cadence_.stepCount())),
          mc_tensor_mol_series_(kMcTensorFields.size(), Mat3Series(cadence_.stepCount())),
          tripeptide_bb_shielding_(cadence_.stepCount()),
          larsen_hbond_shielding_(cadence_.stepCount()),
          tripeptide_bb_shielding_mol_(cadence_.stepCount()),
          larsen_hbond_shielding_mol_(cadence_.stepCount()),
          sasa_(cadence_.stepCount()),
          hbond_nearest_dir_(cadence_.stepCount()),
          hbond_nearest_dir_mol_(cadence_.stepCount()),
          hbond_count_(cadence_.stepCount()),
          ff_pb_radius_(cadence_.stepCount()),
          phi_(cadence_.stepCount()),
          psi_(cadence_.stepCount()),
          omega_(cadence_.stepCount()),
          chi1_(cadence_.stepCount()),
          chi2_(cadence_.stepCount()),
          chi3_(cadence_.stepCount()),
          chi4_(cadence_.stepCount()),
          mopac_charge_(cadence_.stepCount()),
          mopac_s_pop_(cadence_.stepCount()),
          mopac_p_pop_(cadence_.stepCount()),
          mopac_valency_(cadence_.stepCount()),
          mopac_s_character_(cadence_.stepCount()),
          hybridisation_(cadence_.stepCount()),
          contacted_class_membership_(cadence_.stepCount(), 0),
          pi_character_(cadence_.stepCount()) {
        buildStaticBonds();
        readStaticHybridisation();
        configureMolecularFrame();
    }

    void calculate(std::size_t step) {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return;
        coord_.set(step, verbs::pos(body_, atom_, step));
        foldMolecularFrame(step);
        foldMopacSelfState(step);
        foldTopologyJoin(step);
        foldGeometry(step);
        foldDihedrals(step);
        foldFieldsAndRelationships(step);
        foldEfg(step);
        foldProducerSubstrate(step);
        if (cadence_.sigmaPresent(step)) foldSigma(step);
    }

    QJsonObject truth() const {
        QJsonObject root;
        root.insert(QStringLiteral("object_type"), QStringLiteral("atom"));
        root.insert(QStringLiteral("identity"), identityJson());
        root.insert(QStringLiteral("model"), modelJson());
        // ACCUMULATOR_RESPEC §5.1.1: the all-pairs boost grid is GONE; the
        // certified methods feed the compact per-context accumulator instead.
        // There is no `boost` key in the emitted atom object (§2.2A anti-boost).
        root.insert(QStringLiteral("accumulator"), buildAccumulatorJson());
        return root;
    }

    std::size_t sigmaFolds() const { return sigma_folds_; }
    std::size_t relationshipCount() const { return relationships_.size(); }
    std::size_t mappedBondCount() const {
        std::size_t n = 0;
        for (const auto& item : bond_series_)
            if (static_bonds_.count(item.first) != 0) ++n;
        return n;
    }
    std::size_t mismatchEventCount() const { return mismatchCount(); }
    bool oxygenGatePassed() const { return !oxygenGateVerdict().suspect_sp3; }
    std::size_t writeBoundedSigmaRows(QTextStream& out,
                                      const QString& datasetId,
                                      const QString& proteinId) const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return 0;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        QString residueType;
        int residueNumber = 0;
        if (validResidue(p, a.residueIndex)) {
            const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
            residueType = aaName(r.aminoAcid);
            residueNumber = r.address.residueNumber;
        }
        const QString atomUid = uid(QStringLiteral("atom"), atom_);
        const QString atomName = p.atomLabel(atom_, model::NamingConvention::Iupac);
        const QString frameVariant =
            molecular_frame_.variant_ord == static_cast<int>(FrameVariant::Invalid)
                ? molecular_frame_.kind
                : molecular_frame_.variant;
        std::size_t rows = 0;
        std::size_t sigmaOrdinal = 0;
        for (std::size_t step : cadence_.sigmaRows()) {
            const bool sigmaPresent =
                step < sigma_total_raw_.present.size() && sigma_total_raw_.present[step];
            const bool frameValid =
                step < molecular_frame_.values.size() && molecular_frame_.values[step].valid;
            QString frameMissingReason;
            if (molecular_frame_.kind_ord == 0)
                frameMissingReason = QStringLiteral("frame_not_configured");
            else if (!frameValid)
                frameMissingReason = QStringLiteral("frame_degenerate_or_missing_anchor");

            TensorInvariants inv;
            double antisym = kNaN;
            double t1Fraction = kNaN;
            if (sigmaPresent) {
                const Mat3& raw = sigma_total_raw_.values[step];
                inv = tensorInvariants(raw);
                antisym = antisymmetricNorm(raw);
                const double frob = tensorFrobenius(raw);
                if (finite(frob) && frob > 1e-12) t1Fraction = antisym / frob;
            }

            std::array<double, 6> mol{};
            mol.fill(kNaN);
            double molFrob = kNaN;
            double roundtrip = kNaN;
            if (step < sigma_mol_total_.present.size() && sigma_mol_total_.present[step]) {
                mol = symmetricComponents(sigma_mol_total_.values[step]);
                molFrob = tensorFrobenius(sigma_mol_total_.values[step]);
                if (sigmaPresent && step < molecular_frame_.values.size())
                    roundtrip = sigmaRoundTripMaxAbs(molecular_frame_.values[step].axes,
                                                     sigma_total_raw_.values[step]);
            }
            const std::size_t finiteSupportN = sigmaPresent ? 1u : 0u;
            const std::size_t missingSupportN = sigmaPresent ? 0u : 1u;
            const QString missingSupportReason =
                sigmaPresent ? QString() : QStringLiteral("missing_sigma");
            QString hybContext = QStringLiteral("unknown");
            if (static_hybridisation_ord_ >= static_cast<int>(model::Hybridisation::sp)
                && static_hybridisation_ord_ <= static_cast<int>(model::Hybridisation::Unassigned)) {
                hybContext = hybridName(static_cast<model::Hybridisation>(static_hybridisation_ord_)).toLower();
                if (hybContext == QStringLiteral("unassigned")) hybContext = QStringLiteral("unknown");
            }
            const SecondaryStructureState ssContext =
                body_.idx.secondaryStructure.state(atom_, step);
            const QString contactContext = nearestContactClassContext(body_, atom_, step);
            const QString dihedralContext = dihedralContextName(body_.idx, atom_, step);
            const QString ssContextName = ssContext.observed
                                              ? ss3ContextName(ssContext.ss3)
                                              : QStringLiteral("unknown");

            QStringList cells;
            cells << QStringLiteral("bounded_sigma_v1")
                  << csvEscape(datasetId)
                  << csvEscape(proteinId)
                  << csvEscape(atomUid)
                  << QString::number(static_cast<qulonglong>(atom_))
                  << QString::number(a.residueIndex)
                  << QString::number(residueNumber)
                  << csvEscape(residueType)
                  << csvEscape(atomName)
                  << csvEscape(elementName(a.element))
                  << csvEscape(backboneRoleName(a.backboneRole))
                  << csvEscape(molecular_frame_.kind)
                  << QString::number(molecular_frame_.kind_ord)
                  << csvEscape(frameVariant)
                  << QString::number(molecular_frame_.variant_ord)
                  << csvBool(frameValid)
                  << csvEscape(frameMissingReason)
                  << QString::number(static_cast<qulonglong>(cadence_.originalIndex(step)))
                  << QString::number(static_cast<qulonglong>(sigmaOrdinal))
                  << QString::number(static_cast<qulonglong>(step))
                  << csvBool(sigmaPresent)
                  << csvNum(sigmaPresent && step < sigma_total_.present.size() && sigma_total_.present[step]
                                ? sigma_total_.values[step].T0
                                : inv.iso);
            for (double v : mol) cells << csvNum(v);
            cells << csvNum(inv.span)
                  << csvNum(inv.aniso)
                  << csvNum(inv.eta_H)
                  << csvNum(inv.skew)
                  << csvNum(inv.frobenius)
                  << csvNum(antisym)
                  << csvNum(t1Fraction)
                  << csvNum(molFrob)
                  << csvNum(roundtrip)
                  << csvEscape(supportClass(finiteSupportN))
                  << csvNum(static_cast<double>(finiteSupportN))
                  << csvBool(finiteSupportN == 1)
                  << QString::number(static_cast<qulonglong>(missingSupportN))
                  << csvEscape(missingSupportReason)
                  << csvEscape(hybContext)
                  << csvEscape(contactContext)
                  << csvEscape(dihedralContext)
                  << csvEscape(ssContextName);
            out << cells.join(QLatin1Char(',')) << '\n';
            ++rows;
            ++sigmaOrdinal;
        }
        return rows;
    }

    ClassicalSourceTermRecord classicalSourceTermRecord(const QString& datasetId,
                                                        const QString& proteinId) const {
        ClassicalSourceTermRecord record;
        record.dataset_id = datasetId;
        record.protein_id = proteinId;
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return record;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        QString residueType;
        int residueNumber = 0;
        if (validResidue(p, a.residueIndex)) {
            const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
            residueType = aaName(r.aminoAcid);
            residueNumber = r.address.residueNumber;
        }

        auto sumMechanism = [&](const QString& mechanism) {
            std::vector<double> outv(cadence_.stepCount(), kNaN);
            const int ord = mechanismOrd(mechanism);
            for (const auto& item : relationships_) {
                if (item.first.mechanism_ord != ord) continue;
                const std::vector<double> v = item.second.contribution.dense(cadence_.stepCount());
                for (std::size_t i = 0; i < outv.size() && i < v.size(); ++i) {
                    if (!finite(v[i])) continue;
                    outv[i] = finite(outv[i]) ? outv[i] + v[i] : v[i];
                }
            }
            return outv;
        };

        const QString atomName = p.atomLabel(atom_, model::NamingConvention::Iupac);
        record.atom_uid = uid(QStringLiteral("atom"), atom_);
        record.atom_index = atom_;
        record.residue_index = a.residueIndex;
        record.residue_number = residueNumber;
        record.residue_type = residueType;
        record.atom_name = atomName;
        record.element = elementName(a.element);
        record.backbone_role = backboneRoleName(a.backboneRole);
        record.frame_kind = molecular_frame_.kind;
        record.frame_kind_ord = QString::number(molecular_frame_.kind_ord);

        const std::string residueStd = residueType.toStdString();
        const std::string atomStd = atomName.toStdString();
        const std::string frameKindStd = molecular_frame_.kind.toStdString();
        const LiteratureConstant sigma0 = Sigma0(a.element, residueStd, atomStd);
        const LiteratureConstant buckinghamA = BuckinghamA(a.element, residueStd, atomStd, frameKindStd);
        const LiteratureConstant buckinghamB = BuckinghamB(a.element, residueStd, atomStd, frameKindStd);

        record.sigma_qm = sigmaIsoSeries();
        record.sigma0.assign(cadence_.stepCount(), finite(sigma0.value) ? sigma0.value : kNaN);
        auto statusName = [](const LiteratureConstant& c) {
            return QString::fromLatin1(LiteratureStatusName(c.status));
        };
        record.sigma0_status = statusName(sigma0);
        record.buckingham_linear_status = statusName(buckinghamA);
        record.buckingham_linear_key = QString::fromLatin1(buckinghamA.key);
        record.buckingham_linear_constant_value = buckinghamA.value;
        record.buckingham_linear_units = QString::fromLatin1(buckinghamA.units);
        record.buckingham_quadratic_status = statusName(buckinghamB);
        record.buckingham_quadratic_key = QString::fromLatin1(buckinghamB.key);
        record.buckingham_quadratic_constant_value = buckinghamB.value;
        record.buckingham_quadratic_units = QString::fromLatin1(buckinghamB.units);
        for (const LiteratureConstant* c : {&sigma0, &buckinghamA, &buckinghamB})
            if (c->status == LiteratureStatus::Placeholder) ++record.constant_placeholder_n;

        std::vector<double> eParallel = signedBondFieldSeries();
        const bool hasSignedBondField =
            std::any_of(eParallel.begin(), eParallel.end(), [](double v) { return finite(v); });
        if (!hasSignedBondField) eParallel = field_E_z_mopac_.values;
        record.buckingham_linear.assign(cadence_.stepCount(), kNaN);
        record.buckingham_quadratic.assign(cadence_.stepCount(), kNaN);
        for (std::size_t i = 0; i < cadence_.stepCount(); ++i) {
            if (i < eParallel.size() && finite(eParallel[i])) {
                record.buckingham_linear[i] = -buckinghamA.value * eParallel[i];
                record.buckingham_quadratic[i] = -buckinghamB.value * eParallel[i] * eParallel[i];
            }
        }
        record.ring = sumMechanism(QStringLiteral("ring_jb"));
        record.mcconnell = sumMechanism(QStringLiteral("mc_lit_valid"));
        record.larsen = componentSeriesT0(larsen_hbond_shielding_);
        record.sigma_cl.assign(cadence_.stepCount(), kNaN);
        record.residual.assign(cadence_.stepCount(), kNaN);
        for (std::size_t i = 0; i < cadence_.stepCount(); ++i) {
            bool anyPresent = false;
            const double sigmaCl = FoldClassicalSigma(
                sigma0.value,
                {
                    {i < record.buckingham_linear.size() ? record.buckingham_linear[i] : kNaN,
                     i < record.buckingham_linear.size() && finite(record.buckingham_linear[i])},
                    {i < record.buckingham_quadratic.size() ? record.buckingham_quadratic[i] : kNaN,
                     i < record.buckingham_quadratic.size() && finite(record.buckingham_quadratic[i])},
                    {i < record.ring.size() ? record.ring[i] : kNaN,
                     i < record.ring.size() && finite(record.ring[i])},
                    {i < record.mcconnell.size() ? record.mcconnell[i] : kNaN,
                     i < record.mcconnell.size() && finite(record.mcconnell[i])},
                    {i < record.larsen.size() ? record.larsen[i] : kNaN,
                     i < record.larsen.size() && finite(record.larsen[i])},
                },
                &anyPresent);
            if (!anyPresent || !finite(sigmaCl)) continue;
            record.sigma_cl[i] = sigmaCl;
            if (i < record.sigma_qm.size() && finite(record.sigma_qm[i]))
                record.residual[i] = record.sigma_qm[i] - sigmaCl;
        }
        return record;
    }

    std::size_t writeClassicalSourceTermRows(QTextStream& out,
                                             const QString& datasetId,
                                             const QString& proteinId) const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return 0;
        writeClassicalSourceRecord(out, classicalSourceTermRecord(datasetId, proteinId));
        return 1;
    }

    std::size_t writeSourceFamilyMatrixRows(QTextStream& out,
                                            const QString& datasetId,
                                            const QString& proteinId) const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return 0;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        QString residueType;
        int residueNumber = 0;
        if (validResidue(p, a.residueIndex)) {
            const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
            residueType = aaName(r.aminoAcid);
            residueNumber = r.address.residueNumber;
        }
        const QString atomName = p.atomLabel(atom_, model::NamingConvention::Iupac);
        const QString atomUid = uid(QStringLiteral("atom"), atom_);
        const QString identityKey = QStringLiteral("%1:%2:%3").arg(datasetId, proteinId, atomUid);
        const QString contextKey =
            QStringLiteral("%1:%2:%3:%4").arg(proteinId, residueType, atomName, elementName(a.element));

        const std::vector<SubspaceFamily> families = {
            fieldMopacFamily(),
            fieldExternalFamily(),
            efgNodeFamily(),
            ringTensorFamily(QStringLiteral("ring_bs_tensor"), bs_per_type_mol_),
            ringTensorFamily(QStringLiteral("ring_hm_tensor"), hm_per_type_mol_),
            mcAxialFamily(),
            mcRhombicFamily(),
            localElectronicPopulationFamily(),
            localElectronicBondingFamily(),
            larsenHbondFamily(),
        };

        const std::vector<std::size_t> rows = cadence_.sigmaRows();
        std::size_t emitted = 0;
        for (const SubspaceFamily& f : families) {
            QStringList channelNames;
            std::vector<double> means;
            std::vector<double> sds;
            int completeRows = 0;
            for (const SubspaceChannel& ch : f.channels) {
                channelNames << ch.name;
                std::vector<double> finiteVals;
                finiteVals.reserve(rows.size());
                for (std::size_t row : rows)
                    if (row < ch.values.size() && finite(ch.values[row]))
                        finiteVals.push_back(ch.values[row]);
                const DistributionSummary s = SummarizeDistribution(finiteVals);
                means.push_back(s.hasFinite() ? s.mean : kNaN);
                sds.push_back(s.hasFinite() ? s.sd : kNaN);
            }
            for (std::size_t row : rows) {
                bool ok = !f.channels.empty();
                for (const SubspaceChannel& ch : f.channels) {
                    if (row >= ch.values.size() || !finite(ch.values[row])) {
                        ok = false;
                        break;
                    }
                }
                if (ok) ++completeRows;
            }
            const std::size_t finiteRows = completeRows > 0
                ? static_cast<std::size_t>(completeRows)
                : 0u;
            const std::size_t missingRows = rows.size() > finiteRows ? rows.size() - finiteRows : 0u;
            const QString missingReason =
                missingRows == 0 ? QString()
                                 : f.channels.empty()
                                       ? QStringLiteral("no_channels")
                                       : QStringLiteral("nonfinite_channel_row");

            QStringList cells;
            cells << QStringLiteral("source_family_matrix_v1")
                  << csvEscape(datasetId)
                  << csvEscape(proteinId)
                  << csvEscape(atomUid)
                  << QString::number(static_cast<qulonglong>(atom_))
                  << QString::number(a.residueIndex)
                  << QString::number(residueNumber)
                  << csvEscape(residueType)
                  << csvEscape(atomName)
                  << csvEscape(elementName(a.element))
                  << csvEscape(backboneRoleName(a.backboneRole))
                  << csvEscape(identityKey)
                  << csvEscape(contextKey)
                  << csvEscape(QStringLiteral("axis1"))
                  << csvEscape(QStringLiteral("axis1_trajectory"))
                  << csvEscape(f.name)
                  << csvEscape(f.name)
                  << QString::number(static_cast<qulonglong>(f.channels.size()))
                  << QString::number(completeRows)
                  << csvEscape(supportClass(finiteRows))
                  << csvNum(rows.empty() ? kNaN
                                         : static_cast<double>(finiteRows)
                                               / static_cast<double>(rows.size()))
                  << csvBool(finiteRows == 1)
                  << QString::number(static_cast<qulonglong>(missingRows))
                  << csvEscape(missingReason)
                  << csvSemiStrings(channelNames)
                  << csvSemiDoubles(means)
                  << csvSemiDoubles(sds)
                  << csvEscape(QStringLiteral("axis1_ready;cross_axis_join_blocked_pending_axis2"));
            out << cells.join(QLatin1Char(',')) << '\n';
            ++emitted;
        }
        return emitted;
    }

    std::size_t writeSubspaceOverlapRows(QTextStream& out,
                                         const QString& datasetId,
                                         const QString& proteinId) const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return 0;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        QString residueType;
        int residueNumber = 0;
        if (validResidue(p, a.residueIndex)) {
            const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
            residueType = aaName(r.aminoAcid);
            residueNumber = r.address.residueNumber;
        }
        const QString atomName = p.atomLabel(atom_, model::NamingConvention::Iupac);
        const QString atomUid = uid(QStringLiteral("atom"), atom_);
        const QString identityKey = QStringLiteral("%1:%2:%3").arg(datasetId, proteinId, atomUid);
        const QString contextKey =
            QStringLiteral("%1:%2:%3:%4").arg(proteinId, residueType, atomName, elementName(a.element));

        struct ReadSpec {
            QString id;
            SubspaceFamily a;
            SubspaceFamily b;
        };
        const std::vector<ReadSpec> reads = {
            {QStringLiteral("G2_field_sources"),
             fieldMopacFamily(),
             fieldExternalFamily()},
            {QStringLiteral("G6_field_vs_efg"),
             fieldSourcesFamily(),
             efgNodeFamily()},
            {QStringLiteral("G4_bs_vs_hm_tensor_components"),
             ringTensorFamily(QStringLiteral("ring_bs_tensor"), bs_per_type_mol_),
             ringTensorFamily(QStringLiteral("ring_hm_tensor"), hm_per_type_mol_)},
            {QStringLiteral("G6_efg_node"),
             efgAimnetApbsFamily(),
             efgMopacFamily()},
            {QStringLiteral("G6_local_electronic_bundle"),
             localElectronicPopulationFamily(),
             localElectronicBondingFamily()},
            {QStringLiteral("G6_mcconnell_axial_vs_rhombic"),
             mcAxialFamily(),
             mcRhombicFamily()},
        };

        std::size_t emitted = 0;
        const std::vector<std::size_t> rows = cadence_.sigmaRows();
        for (const ReadSpec& read : reads) {
            const SubspaceCompareResult r = CompareSubspaces(read.a, read.b, rows);
            const std::size_t finiteRows = r.finite_n >= 0 ? static_cast<std::size_t>(r.finite_n) : 0u;
            const std::size_t missingRows = rows.size() > finiteRows ? rows.size() - finiteRows : 0u;
            const QString supportMissingReason =
                missingRows == 0 ? QString()
                                 : (r.missing_reason.isEmpty()
                                        ? QStringLiteral("nonfinite_subspace_row")
                                        : r.missing_reason);
            const bool bsHmRead = read.id == QStringLiteral("G4_bs_vs_hm_tensor_components");
            const QString independenceVerdict =
                bsHmRead
                    ? QStringLiteral("independent_forms_checked_distinct_kernel_forms_shared_geometry_projection_decomposition")
                    : QString();
            const QString independenceBasis =
                bsHmRead
                    ? QStringLiteral("Producer finding CODEX_BS_HM_RESEARCH_20260622: BS is a Biot-Savart current-loop kernel and HM is a Haigh-Mallion surface-integral kernel under shared geometry/projection/decomposition; near-null agreement is model robustness, not a relabel")
                    : QString();
            QStringList cells;
            cells << QStringLiteral("subspace_overlap_v1")
                  << csvEscape(datasetId)
                  << csvEscape(proteinId)
                  << csvEscape(atomUid)
                  << QString::number(static_cast<qulonglong>(atom_))
                  << QString::number(a.residueIndex)
                  << QString::number(residueNumber)
                  << csvEscape(residueType)
                  << csvEscape(atomName)
                  << csvEscape(elementName(a.element))
                  << csvEscape(backboneRoleName(a.backboneRole))
                  << csvEscape(identityKey)
                  << csvEscape(contextKey)
                  << csvEscape(QStringLiteral("axis1"))
                  << csvEscape(QStringLiteral("axis1_trajectory"))
                  << csvEscape(read.id)
                  << csvEscape(read.a.name)
                  << csvEscape(read.b.name)
                  << QString::number(r.finite_n)
                  << QString::number(r.input_dim_a)
                  << QString::number(r.input_dim_b)
                  << QString::number(r.active_dim_a)
                  << QString::number(r.active_dim_b)
                  << (r.computed ? QString::number(r.basis_dim_a) : QString())
                  << (r.computed ? QString::number(r.basis_dim_b) : QString())
                  << csvNum(r.explained_fraction_a)
                  << csvNum(r.explained_fraction_b)
                  << csvNum(r.condition_number_a)
                  << csvNum(r.condition_number_b)
                  << csvNum(r.max_canonical_corr)
                  << csvNum(r.mean_canonical_corr)
                  << QString::number(r.n_cc_ge_0_80)
                  << QString::number(r.n_cc_ge_0_95)
                  << csvNum(r.min_angle_deg)
                  << csvSemiDoubles(r.canonical_corrs)
                  << csvSemiDoubles(r.principal_angles_deg)
                  << csvEscape(r.overlap_label)
                  << csvEscape(r.provenance)
                  << csvEscape(r.computed ? QStringLiteral("computed") : QStringLiteral("missing"))
                  << csvEscape(supportMissingReason)
                  << csvEscape(supportClass(finiteRows))
                  << csvNum(rows.empty() ? kNaN
                                         : static_cast<double>(finiteRows)
                                               / static_cast<double>(rows.size()))
                  << csvBool(finiteRows == 1)
                  << QString::number(static_cast<qulonglong>(missingRows))
                  << QString::number(static_cast<qulonglong>(missingRows))
                  << csvSemiStrings(r.dropped_channels_a)
                  << csvSemiStrings(r.dropped_channels_b)
                  << csvEscape(independenceVerdict)
                  << csvEscape(independenceBasis)
                  << csvEscape(QStringLiteral("axis1_ready;cross_axis_join_blocked_pending_axis2"));
            out << cells.join(QLatin1Char(',')) << '\n';
            ++emitted;
        }
        return emitted;
    }

    std::size_t writeEtaByWellRows(QTextStream& out,
                                   const QString& datasetId,
                                   const QString& proteinId) const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return 0;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        QString residueType;
        int residueNumber = 0;
        if (validResidue(p, a.residueIndex)) {
            const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
            residueType = aaName(r.aminoAcid);
            residueNumber = r.address.residueNumber;
        }
        const QString atomName = p.atomLabel(atom_, model::NamingConvention::Iupac);
        const QString atomUid = uid(QStringLiteral("atom"), atom_);

        auto ternaryWells = [&](const std::vector<double>& series) {
            std::vector<int> wells(series.size(), -1);
            std::vector<double> finiteVals;
            finiteVals.reserve(series.size());
            for (double v : series) {
                if (finite(v)) finiteVals.push_back(v);
            }
            if (finiteVals.size() < 4) return wells;
            std::sort(finiteVals.begin(), finiteVals.end());
            const double q1 = finiteVals[finiteVals.size() / 3];
            const double q2 = finiteVals[(2 * finiteVals.size()) / 3];
            for (std::size_t i = 0; i < series.size(); ++i) {
                if (!finite(series[i])) continue;
                wells[i] = series[i] <= q1 ? 0 : (series[i] <= q2 ? 1 : 2);
            }
            return wells;
        };
        auto contactedClassWells = [&]() {
            std::vector<int> wells(contacted_class_membership_.size(), -1);
            const std::uint8_t charged = static_cast<std::uint8_t>(
                contactedClassBit(ContactedClass::ACIDIC) | contactedClassBit(ContactedClass::BASIC));
            const std::uint8_t piOrApolar = static_cast<std::uint8_t>(
                contactedClassBit(ContactedClass::AROMATIC) | contactedClassBit(ContactedClass::APOLAR));
            for (std::size_t i = 0; i < contacted_class_membership_.size(); ++i) {
                const std::uint8_t mask = contacted_class_membership_[i];
                wells[i] = (mask & charged) != 0 ? 2
                         : (mask & piOrApolar) != 0 ? 1
                         : 0;
            }
            return wells;
        };
        struct WellSpec {
            QString family;
            QString source;
            std::vector<int> wells;
        };
        const std::vector<WellSpec> wells = {
            {QStringLiteral("dihedral_rotamer"),
             QStringLiteral("rotamerRelevantDihedralWells"),
             rotamerRelevantDihedralWells()},
            {QStringLiteral("contacted_class"),
             QStringLiteral("contacted_class_membership"),
             contactedClassWells()},
            {QStringLiteral("hbond_count_tertiles"),
             QStringLiteral("hbond.count"),
             ternaryWells(hbond_count_.values)},
            {QStringLiteral("sasa_tertiles"),
             QStringLiteral("sasa"),
             ternaryWells(sasa_.values)},
            {QStringLiteral("mopac_charge_tertiles"),
             QStringLiteral("mopac.charge"),
             ternaryWells(mopac_charge_.values)},
            {QStringLiteral("field_abs_E_mopac_tertiles"),
             QStringLiteral("field.mopac_coulomb.abs_E"),
             ternaryWells(field_abs_E_mopac_.values)},
        };

        struct TargetSpec {
            QString key;
            QString frame;
            std::vector<double> values;
        };
        std::vector<TargetSpec> targets = {
            {QStringLiteral("sigma.iso"), QStringLiteral("lab"), sigmaIsoSeries()},
            {QStringLiteral("sigma.invariant.span"), QStringLiteral("pas_invariant"),
             sigmaInvariantSeries(QStringLiteral("span"))},
            {QStringLiteral("sigma.invariant.aniso"), QStringLiteral("pas_invariant"),
             sigmaInvariantSeries(QStringLiteral("aniso"))},
            {QStringLiteral("sigma.invariant.eta_H"), QStringLiteral("pas_invariant"),
             sigmaInvariantSeries(QStringLiteral("eta_H"))},
            {QStringLiteral("sigma.shape.t1_fraction"), QStringLiteral("lab_shape"),
             sigmaT1FractionSeries()},
        };
        for (int sc = 0; sc < 6; ++sc) {
            const auto comp = static_cast<MolComp>(sc);
            targets.push_back({
                QStringLiteral("sigma.total.molcomp:%1")
                    .arg(QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(sc)])),
                QStringLiteral("molecular"),
                sigmaMolCompSeries(sigma_mol_total_, comp),
            });
        }

        const std::vector<std::size_t> rows = cadence_.sigmaRows();
        std::size_t emitted = 0;
        for (const WellSpec& well : wells) {
            for (const TargetSpec& target : targets) {
                std::array<int, 3> counts = {0, 0, 0};
                std::size_t finiteN = 0;
                for (std::size_t row : rows) {
                    if (row >= target.values.size() || row >= well.wells.size()) continue;
                    if (!finite(target.values[row]) || well.wells[row] < 0 || well.wells[row] > 2)
                        continue;
                    ++counts[static_cast<std::size_t>(well.wells[row])];
                    ++finiteN;
                }
                const std::size_t missingN = rows.size() > finiteN ? rows.size() - finiteN : 0u;
                QStringList countStrings;
                for (int c : counts) countStrings << QString::number(c);
                QStringList cells;
                cells << QStringLiteral("eta2_by_well_v1")
                      << csvEscape(datasetId)
                      << csvEscape(proteinId)
                      << csvEscape(atomUid)
                      << QString::number(static_cast<qulonglong>(atom_))
                      << QString::number(a.residueIndex)
                      << QString::number(residueNumber)
                      << csvEscape(residueType)
                      << csvEscape(atomName)
                      << csvEscape(elementName(a.element))
                      << csvEscape(backboneRoleName(a.backboneRole))
                      << csvEscape(well.family)
                      << csvEscape(well.source)
                      << csvEscape(target.key)
                      << csvEscape(target.frame)
                      << QString::number(static_cast<qulonglong>(finiteN))
                      << csvEscape(supportClass(finiteN))
                      << csvNum(rows.empty() ? kNaN
                                             : static_cast<double>(finiteN)
                                                   / static_cast<double>(rows.size()))
                      << csvBool(finiteN == 1)
                      << QString::number(static_cast<qulonglong>(missingN))
                      << csvEscape(missingN == 0 ? QString()
                                                 : QStringLiteral("nonfinite_target_or_well"))
                      << csvNum(etaSquaredByWell(target.values, well.wells, rows))
                      << csvEscape(countStrings.join(QLatin1Char(';')));
                out << cells.join(QLatin1Char(',')) << '\n';
                ++emitted;
            }
        }
        return emitted;
    }
    mutable std::size_t last_accumulator_response_count = 0;
    mutable std::size_t last_accumulator_context_count = 0;

private:
    struct StaticBondInfo {
        int index = -1;
        model::BondOrder order = model::BondOrder::Unknown;
        model::BondCategory category = model::BondCategory::Unknown;
        std::size_t other = 0;
    };

    struct OxygenGateVerdict {
        bool checked = false;
        bool suspect_sp3 = false;
    };

    static bool oxygenExpectedSp2(const model::QtAtom& a) {
        return a.IsBackboneCarbonylOxygen()
               || a.IsSidechainCarboxylateOxygen()
               || a.IsSidechainAmideOxygen()
               || a.planarGroup == model::PlanarGroupKind::Carboxylate
               || a.planarGroup == model::PlanarGroupKind::PeptideAmide
               || a.planarGroup == model::PlanarGroupKind::SidechainAmide;
    }

    static bool oxygenExpectedSp3(const model::QtAtom& a) {
        return a.ffAtomType == model::QtFfAtomType::OH
               || a.ffAtomType == model::QtFfAtomType::OS
               || a.planarGroup == model::PlanarGroupKind::None
               || a.planarGroup == model::PlanarGroupKind::AromaticHydroxyl;
    }

    OxygenGateVerdict oxygenGateVerdict() const {
        OxygenGateVerdict verdict;
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return verdict;
        const model::QtAtom& a = body_.run.protein->atom(atom_);
        if (a.element != model::Element::O) return verdict;

        verdict.checked = oxygen_gate_checked_;
        if (!verdict.checked) return verdict;

        const bool expectedSp2 = oxygenExpectedSp2(a);
        const bool expectedSp3 = oxygenExpectedSp3(a) && !expectedSp2;
        if (!expectedSp2 && !expectedSp3) return verdict;

        std::size_t n = 0;
        double sumPi = 0.0;
        double maxPi = -std::numeric_limits<double>::infinity();
        std::size_t sp2Like = 0;
        std::size_t sp3Like = 0;
        for (std::size_t i = 0; i < pi_character_.values.size(); ++i) {
            const double pi = pi_character_.values[i];
            if (!finite(pi)) continue;
            ++n;
            sumPi += pi;
            maxPi = std::max(maxPi, pi);
            const int h = i < hybridisation_.values.size() ? hybridisation_.values[i] : IntSeries::kMissing;
            if (h == static_cast<int>(model::Hybridisation::sp2)
                || h == static_cast<int>(model::Hybridisation::sp)) {
                ++sp2Like;
            } else if (h == static_cast<int>(model::Hybridisation::sp3)) {
                ++sp3Like;
            }
        }
        if (n == 0) {
            verdict.suspect_sp3 = true;
            return verdict;
        }

        const double meanPi = sumPi / static_cast<double>(n);
        constexpr double kSp2Boundary = 0.4;
        constexpr double kDeadband = 0.05;
        if (expectedSp2) {
            verdict.suspect_sp3 = (meanPi < (kSp2Boundary - kDeadband))
                                  || (sp2Like == 0 && maxPi < kSp2Boundary);
        } else if (expectedSp3) {
            verdict.suspect_sp3 = (meanPi > (kSp2Boundary + kDeadband))
                                  && (sp2Like > sp3Like);
        }
        return verdict;
    }

    QJsonObject identityJson() const {
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        QJsonObject o;
        o.insert(QStringLiteral("uid"), uid(QStringLiteral("atom"), atom_));
        o.insert(QStringLiteral("object_type"), QStringLiteral("atom"));
        o.insert(QStringLiteral("model_index"), static_cast<int>(atom_));
        o.insert(QStringLiteral("atom_index"), static_cast<int>(atom_));
        o.insert(QStringLiteral("residue_index"), a.residueIndex);
        o.insert(QStringLiteral("parent_atom_index"), a.parentAtomIndex);
        o.insert(QStringLiteral("element_ord"), static_cast<int>(a.element));
        o.insert(QStringLiteral("element"), elementName(a.element));
        o.insert(QStringLiteral("atom_name"), p.atomLabel(atom_, model::NamingConvention::Iupac));
        o.insert(QStringLiteral("backbone_role_ord"), static_cast<int>(a.backboneRole));
        o.insert(QStringLiteral("backbone_role"), backboneRoleName(a.backboneRole));
        o.insert(QStringLiteral("ff_atom_type_ord"), static_cast<int>(a.ffAtomType));
        o.insert(QStringLiteral("ff_atom_type"), ffAtomTypeName(a.ffAtomType));
        o.insert(QStringLiteral("formal_charge"), a.formalCharge);
        o.insert(QStringLiteral("aromatic"), a.aromatic);
        o.insert(QStringLiteral("is_exchangeable"), a.isExchangeable);
        return o;
    }

    QJsonObject typingJson() const {
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        QJsonObject t;
        t.insert(QStringLiteral("element"), elementName(a.element));
        t.insert(QStringLiteral("element_ord"), static_cast<int>(a.element));
        t.insert(QStringLiteral("atom_name"), p.atomLabel(atom_, model::NamingConvention::Iupac));
        t.insert(QStringLiteral("backbone_role"), backboneRoleName(a.backboneRole));
        t.insert(QStringLiteral("backbone_role_ord"), static_cast<int>(a.backboneRole));
        return t;
    }

    QJsonObject modelJson() const {
        QJsonObject m;
        m.insert(QStringLiteral("typing"), typingJson());
        m.insert(QStringLiteral("category"), categoryJson());
        m.insert(QStringLiteral("topology_join"), topologyJoinJson());
        m.insert(QStringLiteral("substrate_parameters"), substrateParametersJson());
        m.insert(QStringLiteral("graph_topology"), graphTopologyJson());
        return m;
    }

    QJsonObject substrateParametersJson() const {
        QJsonObject m;
        QJsonObject charge;
        charge.insert(QStringLiteral("field_key"), QStringLiteral("field.charge_ff14sb"));
        charge.insert(QStringLiteral("cutoff_A"), config_.charge_cutoff_A);
        charge.insert(QStringLiteral("same_residue_sources_included"), true);
        charge.insert(QStringLiteral("self_sources_rejected_with_flag"), true);
        charge.insert(QStringLiteral("degenerate_geometry_rejected_with_flag"), true);
        charge.insert(QStringLiteral("unit_conversion_ke_V_A_per_e_over_A2"), CoulombKeVA());
        m.insert(QStringLiteral("charge_field"), charge);

        QJsonObject ring;
        ring.insert(QStringLiteral("cutoff_A"), config_.ring_cutoff_A);
        QJsonArray intensities;
        auto addRing = [&](model::RingTypeIndex type) {
            const LiteratureConstant c = RingIntensity(type);
            QJsonObject o;
            o.insert(QStringLiteral("ring_type"), ringTypeName(type));
            o.insert(QStringLiteral("literature_intensity_nA_per_T"), c.value);
            o.insert(QStringLiteral("constant_status"),
                     QString::fromLatin1(LiteratureStatusName(c.status)));
            intensities.append(o);
        };
        for (int type = 0; type < model::kAromaticRingTypeCount; ++type)
            addRing(static_cast<model::RingTypeIndex>(type));
        ring.insert(QStringLiteral("literature_intensities"), intensities);
        m.insert(QStringLiteral("ring_current"), ring);

        QJsonObject mc;
        mc.insert(QStringLiteral("near_field_ratio"), config_.mc_near_field_ratio);
        mc.insert(QStringLiteral("molar_prefactor"), McConnellMolarPrefactor());
        QJsonArray dchi;
        for (model::BondCategory category : {model::BondCategory::PeptideCO,
                                             model::BondCategory::PeptideCN,
                                             model::BondCategory::SidechainCO,
                                             model::BondCategory::Aromatic}) {
            QJsonObject o;
            o.insert(QStringLiteral("bond_category"), bondCategoryName(category));
            o.insert(QStringLiteral("delta_chi_q_1e_minus_6_cm3_per_mol"),
                     McConnellDeltaChiQ(category));
            dchi.append(o);
        }
        mc.insert(QStringLiteral("delta_chi_literature"), dchi);
        m.insert(QStringLiteral("mcconnell"), mc);
        return m;
    }

    QJsonObject categoryJson() const {
        const model::QtAtom& a = body_.run.protein->atom(atom_);
        QJsonObject c;
        c.insert(QStringLiteral("element"), elementName(a.element));
        c.insert(QStringLiteral("element_ord"), static_cast<int>(a.element));
        c.insert(QStringLiteral("hybridisation_series"), hybridisation_.json());
        c.insert(QStringLiteral("hybridisation_pi_character_series"), pi_character_.json());
        c.insert(QStringLiteral("s_character_series"), mopac_s_character_.json());
        c.insert(QStringLiteral("contacted_class_membership_series"),
                 contactedClassMembershipSeriesJson());
        c.insert(QStringLiteral("static_hybridisation_ord"), static_hybridisation_ord_);
        c.insert(QStringLiteral("static_hybridisation"), hybridName(static_cast<model::Hybridisation>(static_hybridisation_ord_)));
        QJsonArray events;
        int prev = IntSeries::kMissing;
        for (std::size_t i = 0; i < hybridisation_.values.size(); ++i) {
            const int cur = hybridisation_.values[i];
            if (cur == IntSeries::kMissing) continue;
            if (prev != IntSeries::kMissing && cur != prev) {
                QJsonObject e;
                e.insert(QStringLiteral("step_row"), static_cast<int>(i));
                e.insert(QStringLiteral("from_ord"), prev);
                e.insert(QStringLiteral("from"), hybridName(static_cast<model::Hybridisation>(prev)));
                e.insert(QStringLiteral("to_ord"), cur);
                e.insert(QStringLiteral("to"), hybridName(static_cast<model::Hybridisation>(cur)));
                events.append(e);
            }
            prev = cur;
        }
        c.insert(QStringLiteral("hybridisation_change_events"), events);
        c.insert(QStringLiteral("contacted_class_change_events"), contactedClassChangeEventsJson());
        QJsonObject contactSelection;
        contactSelection.insert(QStringLiteral("per_frame_policy"),
                                QStringLiteral("MEMBERSHIP_OF_ALL_PRESENT_CONTACT_CLASSES"));
        contactSelection.insert(QStringLiteral("multiple_classes_policy"),
                                QStringLiteral("frame contributes to every present contacted class; contacted-class regimes overlap"));
        contactSelection.insert(QStringLiteral("none_policy"),
                                QStringLiteral("NONE when no classified present contact is folded for the frame"));
        contactSelection.insert(QStringLiteral("regime_occupancy_policy"),
                                QStringLiteral("overlapping contacted-class memberships; occupancy sums may exceed the frame count"));
        c.insert(QStringLiteral("contacted_class_selection"), contactSelection);
        QJsonArray contacted;
        for (const auto& r : contacted_residues_) {
            QJsonObject o;
            o.insert(QStringLiteral("residue_number"), r.first);
            o.insert(QStringLiteral("amino_acid_ord"), r.second);
            o.insert(QStringLiteral("amino_acid"), aaName(static_cast<model::AminoAcid>(r.second)));
            contacted.append(o);
        }
        c.insert(QStringLiteral("contacted_residues"), contacted);
        if (a.element == model::Element::O) {
            const OxygenGateVerdict oxygenGate = oxygenGateVerdict();
            QJsonObject gate;
            gate.insert(QStringLiteral("checked"), oxygenGate.checked);
            gate.insert(QStringLiteral("suspect_sp3"), oxygenGate.suspect_sp3);
            c.insert(QStringLiteral("oxygen_label_gate"), gate);
        }
        return c;
    }

    QJsonObject topologyJoinJson() const {
        QJsonObject t;
        QJsonArray mapped;
        for (const auto& item : bond_series_) {
            const auto sit = static_bonds_.find(item.first);
            if (sit == static_bonds_.end()) continue;
            const auto p = pairFromKey(item.first);
            QJsonObject o;
            o.insert(QStringLiteral("atom_a"), static_cast<int>(p.first));
            o.insert(QStringLiteral("atom_b"), static_cast<int>(p.second));
            o.insert(QStringLiteral("static_bond_index"), sit->second.index);
            o.insert(QStringLiteral("static_order_ord"), static_cast<int>(sit->second.order));
            o.insert(QStringLiteral("static_order"), bondOrderName(sit->second.order));
            o.insert(QStringLiteral("static_category_ord"), static_cast<int>(sit->second.category));
            o.insert(QStringLiteral("static_category"), bondCategoryName(sit->second.category));
            o.insert(QStringLiteral("mopac_wiberg_order_series"), item.second.json());
            double sum = 0.0;
            int n = 0;
            for (double v : item.second.values)
                if (finite(v)) {
                    sum += v;
                    ++n;
                }
            o.insert(QStringLiteral("mopac_wiberg_order_mean"), n > 0 ? QJsonValue(sum / n) : QJsonValue(QJsonValue::Null));
            o.insert(QStringLiteral("is_partner_bond"), true);
            mapped.append(o);
        }
        t.insert(QStringLiteral("mapped_bonds"), mapped);

        QJsonArray mismatch;
        for (const auto& item : bond_series_) {
            if (static_bonds_.count(item.first) != 0) continue;
            const auto p = pairFromKey(item.first);
            QJsonObject o;
            o.insert(QStringLiteral("kind"), QStringLiteral("mopac_only"));
            o.insert(QStringLiteral("atom_a"), static_cast<int>(p.first));
            o.insert(QStringLiteral("atom_b"), static_cast<int>(p.second));
            o.insert(QStringLiteral("mopac_wiberg_order_series"), item.second.json());
            o.insert(QStringLiteral("static_bond_index"), -1);
            o.insert(QStringLiteral("static_order_ord"), -1);
            int n = 0;
            for (double v : item.second.values)
                if (finite(v)) ++n;
            o.insert(QStringLiteral("n_steps_present"), n);
            mismatch.append(o);
        }
        for (const auto& item : static_bonds_) {
            const auto bit = bond_series_.find(item.first);
            int missing = static_cast<int>(cadence_.stepCount());
            if (bit != bond_series_.end()) {
                missing = 0;
                for (double v : bit->second.values)
                    if (!finite(v)) ++missing;
            }
            if (bit != bond_series_.end() && missing == 0) continue;
            const auto p = pairFromKey(item.first);
            QJsonObject o;
            o.insert(QStringLiteral("kind"), QStringLiteral("static_only"));
            o.insert(QStringLiteral("atom_a"), static_cast<int>(p.first));
            o.insert(QStringLiteral("atom_b"), static_cast<int>(p.second));
            o.insert(QStringLiteral("mopac_wiberg_order_series"), QJsonValue(QJsonValue::Null));
            o.insert(QStringLiteral("static_bond_index"), item.second.index);
            o.insert(QStringLiteral("static_order_ord"), static_cast<int>(item.second.order));
            o.insert(QStringLiteral("n_steps_present"), missing);
            mismatch.append(o);
        }
        t.insert(QStringLiteral("mismatch_events"), mismatch);
        return t;
    }

    template <typename Predicate>
    int minGraphHopsTo(Predicate pred) const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return -1;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtTopology& topo = p.topology();
        std::vector<int> dist(p.atomCount(), -1);
        std::vector<std::size_t> queue;
        queue.reserve(p.atomCount());
        dist[atom_] = 0;
        queue.push_back(atom_);
        for (std::size_t qi = 0; qi < queue.size(); ++qi) {
            const std::size_t cur = queue[qi];
            if (pred(p.atom(cur))) return dist[cur];
            for (int32_t bi : topo.bondIndicesForAtom(cur)) {
                if (bi < 0 || static_cast<std::size_t>(bi) >= topo.bondCount()) continue;
                const model::QtBond& b = topo.bondAt(static_cast<std::size_t>(bi));
                const int32_t next = b.atomIndexA == static_cast<int32_t>(cur) ? b.atomIndexB : b.atomIndexA;
                if (next < 0 || static_cast<std::size_t>(next) >= p.atomCount()) continue;
                const std::size_t ni = static_cast<std::size_t>(next);
                if (dist[ni] >= 0) continue;
                dist[ni] = dist[cur] + 1;
                queue.push_back(ni);
            }
        }
        return -1;
    }

    QJsonObject graphTopologyJson() const {
        QJsonObject g;
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return g;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        const model::QtTopology& topo = p.topology();
        const auto& bondIndices = topo.bondIndicesForAtom(atom_);
        int degreeHeavy = 0;
        int degreeHydrogen = 0;
        for (int32_t bi : bondIndices) {
            if (bi < 0 || static_cast<std::size_t>(bi) >= topo.bondCount()) continue;
            const model::QtBond& b = topo.bondAt(static_cast<std::size_t>(bi));
            const int32_t other = b.atomIndexA == static_cast<int32_t>(atom_) ? b.atomIndexB : b.atomIndexA;
            if (other < 0 || static_cast<std::size_t>(other) >= p.atomCount()) continue;
            if (p.atom(static_cast<std::size_t>(other)).element == model::Element::H)
                ++degreeHydrogen;
            else
                ++degreeHeavy;
        }
        g.insert(QStringLiteral("definition"),
                 QStringLiteral("covalent topology graph; hop counts are bond-edge counts; ring centroid hops use the nearest ring-member atom as the graph anchor"));
        g.insert(QStringLiteral("degree_total"), static_cast<int>(bondIndices.size()));
        g.insert(QStringLiteral("degree_heavy"), degreeHeavy);
        g.insert(QStringLiteral("degree_hydrogen"), degreeHydrogen);
        g.insert(QStringLiteral("is_hydrogen"), a.element == model::Element::H);
        g.insert(QStringLiteral("is_heavy_atom"), a.element != model::Element::H);
        g.insert(QStringLiteral("is_ring_atom"), a.IsInAnyRing());
        g.insert(QStringLiteral("is_aromatic"), a.aromatic);
        g.insert(QStringLiteral("parent_heavy_atom_index"), a.parentAtomIndex);
        g.insert(QStringLiteral("nearest_backbone_hops"),
                 minGraphHopsTo([](const model::QtAtom& x) { return x.IsBackbone(); }));
        const int ringHops = minGraphHopsTo([](const model::QtAtom& x) { return x.IsInAnyRing(); });
        g.insert(QStringLiteral("nearest_ring_atom_hops"), ringHops);
        g.insert(QStringLiteral("nearest_ring_centroid_hops"), ringHops);
        g.insert(QStringLiteral("electronegativity_scale"), QStringLiteral("Pauling"));
        g.insert(QStringLiteral("electronegativity_pauling"), a.Electronegativity());
        g.insert(QStringLiteral("covalent_radius_A"), a.CovalentRadius());
        g.insert(QStringLiteral("atomic_number"), a.AtomicNumber());
        if (validResidue(p, a.residueIndex)) {
            g.insert(QStringLiteral("residue_atom_count"),
                     static_cast<int>(p.residue(static_cast<std::size_t>(a.residueIndex)).atomIndices.size()));
        }
        return g;
    }

    QJsonArray bondLengthsJson() const {
        QJsonArray out;
        if (!body_.run.protein) return out;
        const model::QtProtein& p = *body_.run.protein;
        for (const auto& item : static_bonds_) {
            const std::uint64_t key = item.first;
            const StaticBondInfo& info = item.second;
            const auto pair = pairFromKey(key);
            QJsonObject o;
            o.insert(QStringLiteral("atom_a"), static_cast<int>(pair.first));
            o.insert(QStringLiteral("atom_b"), static_cast<int>(pair.second));
            o.insert(QStringLiteral("other_atom_index"), static_cast<int>(info.other));
            if (info.other < p.atomCount())
                o.insert(QStringLiteral("other_atom_name"),
                         p.atomLabel(info.other, model::NamingConvention::Iupac));
            o.insert(QStringLiteral("static_bond_index"), info.index);
            o.insert(QStringLiteral("static_order_ord"), static_cast<int>(info.order));
            o.insert(QStringLiteral("static_order"), bondOrderName(info.order));
            o.insert(QStringLiteral("static_category_ord"), static_cast<int>(info.category));
            o.insert(QStringLiteral("static_category"), bondCategoryName(info.category));
            bool disulfideSS = false;
            if (pair.first < p.atomCount() && pair.second < p.atomCount()) {
                disulfideSS = info.category == model::BondCategory::Disulfide
                              && p.atom(pair.first).element == model::Element::S
                              && p.atom(pair.second).element == model::Element::S;
            }
            o.insert(QStringLiteral("is_disulfide_ss"), disulfideSS);
            const auto lenIt = bond_length_series_.find(key);
            const auto rejIt = bond_length_rejected_degenerate_.find(key);
            o.insert(QStringLiteral("length_A"),
                     lenIt != bond_length_series_.end() ? QJsonValue(lenIt->second.json())
                                                        : QJsonValue(QJsonValue::Null));
            o.insert(QStringLiteral("rejected_degenerate"),
                     rejIt != bond_length_rejected_degenerate_.end()
                         ? QJsonValue(rejIt->second.json())
                         : QJsonValue(QJsonValue::Null));
            out.append(o);
        }
        return out;
    }

    std::size_t mismatchCount() const {
        std::size_t n = 0;
        for (const auto& item : bond_series_)
            if (static_bonds_.count(item.first) == 0) ++n;
        for (const auto& item : static_bonds_) {
            const auto bit = bond_series_.find(item.first);
            if (bit == bond_series_.end()) {
                ++n;
                continue;
            }
            if (std::any_of(bit->second.values.begin(), bit->second.values.end(),
                            [](double v) { return !finite(v); }))
                ++n;
        }
        return n;
    }

    QJsonObject relationshipJson(const RelationshipKey& key, const RelationshipSeries& rel, int) const {
        QJsonObject o;
        QJsonObject k;
        const auto hyb = static_cast<model::Hybridisation>(key.hybridisation_ord);
        const auto elem = static_cast<model::Element>(key.element_ord);
        k.insert(QStringLiteral("chemistry"),
                 QStringLiteral("hyb=%1|element=%2").arg(hybridName(hyb), elementName(elem)));
        k.insert(QStringLiteral("hybridisation_ord"), key.hybridisation_ord);
        k.insert(QStringLiteral("element_ord"), key.element_ord);
        k.insert(QStringLiteral("contacted_residue"),
                 QStringLiteral("residue=%1|aa=%2")
                     .arg(key.contacted_residue_number)
                     .arg(aaName(static_cast<model::AminoAcid>(key.contacted_amino_acid_ord))));
        k.insert(QStringLiteral("contacted_residue_number"), key.contacted_residue_number);
        k.insert(QStringLiteral("contacted_amino_acid_ord"), key.contacted_amino_acid_ord);
        k.insert(QStringLiteral("scope_ord"), key.scope_ord);
        k.insert(QStringLiteral("scope"), scopeName(key.scope_ord));
        k.insert(QStringLiteral("mechanism_ord"), key.mechanism_ord);
        k.insert(QStringLiteral("mechanism"), mechanismName(key.mechanism_ord));
        k.insert(QStringLiteral("source_kind_ord"), key.source_kind_ord);
        k.insert(QStringLiteral("source_kind"), sourceKindName(key.source_kind_ord));
        k.insert(QStringLiteral("source_id"), key.source_id);
        k.insert(QStringLiteral("source_atom_index"), key.source_atom_index);
        k.insert(QStringLiteral("source_category_ord"), key.source_category_ord);
        o.insert(QStringLiteral("key"), k);
        // The partner census is FOLDED, not retained: every partner stays enumerated
        // (the key above), but each facet is a characterization (summary + shape), not a
        // per-frame series — "who I danced with, via which mechanism, how often, and what
        // each did," as character rather than rows.
        const std::size_t n = cadence_.stepCount();
        QJsonObject present;
        present.insert(QStringLiteral("n"), static_cast<int>(rel.present_steps.size()));
        if (!rel.present_steps.empty()) {
            present.insert(QStringLiteral("first"), static_cast<int>(rel.present_steps.front()));
            present.insert(QStringLiteral("last"), static_cast<int>(rel.present_steps.back()));
        }
        o.insert(QStringLiteral("present"), present);
        QJsonObject facets;
        facets.insert(QStringLiteral("distance"), scalarSummaryJson(rel.distance.dense(n)));
        facets.insert(QStringLiteral("cos_theta"), scalarSummaryJson(rel.cos_theta.dense(n)));
        facets.insert(QStringLiteral("inv_r3"), scalarSummaryJson(rel.inv_r3.dense(n)));
        facets.insert(QStringLiteral("dipolar"), scalarSummaryJson(rel.dipolar.dense(n)));
        facets.insert(QStringLiteral("kernel_T0"), scalarSummaryJson(rel.kernel_T0.dense(n)));
        QJsonObject kernelMol;
        for (int c = 0; c < 6; ++c) {
            kernelMol.insert(QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(c)]),
                             scalarSummaryJson(rel.kernel_mol.componentDense(n, static_cast<MolComp>(c))));
        }
        facets.insert(QStringLiteral("kernel_mol_components"), kernelMol);
        facets.insert(QStringLiteral("contribution"), scalarSummaryJson(rel.contribution.dense(n)));
        o.insert(QStringLiteral("facets"), facets);
        return o;
    }

    // =======================================================================
    // THE COMPACT PER-CONTEXT ACCUMULATOR (ACCUMULATOR_RESPEC PART 2).
    // The atom folds its own response by context and tells its truth. The
    // certified methods (ols/pearsonR/chatterjeeXi/circular-shift null/foldCsa
    // eta_H/serial dynamics) are re-pointed at the MOLECULAR-frame quantities
    // and the invariants; the fixed characterization catalogue (§2.2A) is the only door; the cap
    // bounded by the fixed characterization catalogue; there is NO boost key.
    // =======================================================================

    // -- molecular-frame sigma component series (symmetrized), on every step.
    std::vector<double> sigmaMolCompSeries(const Mat3Series& mol, MolComp comp) const {
        std::vector<double> out(mol.values.size(), kNaN);
        for (std::size_t i = 0; i < mol.values.size(); ++i) {
            if (!mol.present[i]) continue;
            out[i] = symmetricComponents(mol.values[i])[static_cast<std::size_t>(comp)];
        }
        return out;
    }

    // -- sigma invariant series from sigma_mol.total (frame-independent, so
    // computed from the lab raw tensor where present; molecular not required).
    std::vector<double> sigmaInvariantSeries(const QString& which) const {
        std::vector<double> out(sigma_total_raw_.values.size(), kNaN);
        for (std::size_t i = 0; i < sigma_total_raw_.values.size(); ++i) {
            if (!sigma_total_raw_.present[i]) continue;
            const TensorInvariants inv = tensorInvariants(sigma_total_raw_.values[i]);
            if (which == QStringLiteral("iso")) out[i] = inv.iso;
            else if (which == QStringLiteral("span")) out[i] = inv.span;
            else if (which == QStringLiteral("aniso")) out[i] = inv.aniso;
            else if (which == QStringLiteral("eta_H")) out[i] = inv.eta_H;
            else if (which == QStringLiteral("frobenius")) out[i] = inv.frobenius;
        }
        return out;
    }

    std::vector<double> sigmaT1FractionSeries() const {
        std::vector<double> out(sigma_total_raw_.values.size(), kNaN);
        for (std::size_t i = 0; i < sigma_total_raw_.values.size(); ++i) {
            if (!sigma_total_raw_.present[i]) continue;
            const double antisym = antisymmetricNorm(sigma_total_raw_.values[i]);
            const double frob = tensorFrobenius(sigma_total_raw_.values[i]);
            if (finite(antisym) && finite(frob) && frob > 1e-12) out[i] = antisym / frob;
        }
        return out;
    }

    // -- sigma iso T0 series (the scalar target).
    std::vector<double> sigmaIsoSeries() const { return componentSeriesT0(sigma_total_); }

    // -- driver invariant magnitude series helpers.
    std::vector<double> fieldAbsSeries() const { return field_abs_E_mopac_.values; }
    std::vector<double> efgAbsSeries() const {
        // best of the available |EFG| channels (aimnet2 primary, shielding fallback).
        std::vector<double> out = efg_abs_aimnet2_.values;
        bool any = std::any_of(out.begin(), out.end(), [](double v) { return finite(v); });
        if (!any) out = shielding_abs_mopac_coulomb_.values;
        return out;
    }

    // -- molecular field component series.
    std::vector<double> fieldMolCompSeries(int c) const {
        return componentSeries(field_mol_mopac_, c);
    }

    // -- EFG molecular component series from a specific source tensor series.
    std::vector<double> molCompFrom(const Mat3Series& src, MolComp comp) const {
        std::vector<double> out(src.values.size(), kNaN);
        for (std::size_t i = 0; i < src.values.size(); ++i)
            if (src.present[i])
                out[i] = symmetricComponents(src.values[i])[static_cast<std::size_t>(comp)];
        return out;
    }

    bool t2Complete(const std::array<double, 5>& t2) const {
        for (double v : t2)
            if (!finite(v)) return false;
        return true;
    }

    std::optional<Mat3> projectLabTensorToMolecular(std::size_t step, const Mat3& lab) const {
        if (!lab.allFinite()) return std::nullopt;
        const auto axes = molecularAxesAt(step);
        if (!axes) return std::nullopt;
        const Mat3 local = axes->transpose() * lab * (*axes);
        if (!local.allFinite()) return std::nullopt;
        return local;
    }

    std::optional<Mat3> projectLibraryT2ToMolecular(std::size_t step,
                                                    double t0,
                                                    const std::array<double, 5>& t2) const {
        if (!t2Complete(t2)) return std::nullopt;
        const Mat3 lab = ReconstructLibraryT2Matrix(finite(t0) ? t0 : 0.0, t2);
        return projectLabTensorToMolecular(step, lab);
    }

    std::optional<Mat3> projectSphericalTensorToMolecular(std::size_t step,
                                                          const model::SphericalTensor& t) const {
        if (!finite(t.T0) || !t2Complete(t.T2)) return std::nullopt;
        return projectLibraryT2ToMolecular(step, t.T0, t.T2);
    }

    void projectFixedT2Array(std::size_t step,
                             const std::array<double, 40>& values,
                             std::vector<Mat3Series>& out) {
        for (int type = 0; type < 8; ++type) {
            std::array<double, 5> t2{};
            for (int c = 0; c < 5; ++c)
                t2[static_cast<std::size_t>(c)] =
                    values[static_cast<std::size_t>(type * 5 + c)];
            if (const auto local = projectLibraryT2ToMolecular(step, 0.0, t2))
                out[static_cast<std::size_t>(type)].set(step, *local);
        }
    }

    static double tensorDoubleDot(const Mat3& a, const Mat3& b) {
        return (a.array() * b.array()).sum();
    }

    Mat3 pairedMeanSigmaTensor(const Mat3Series& sigma,
                               const Mat3Series& mechanism,
                               const std::vector<std::size_t>& rows) const {
        Mat3 sum = Mat3::Zero();
        int n = 0;
        for (std::size_t row : rows) {
            if (row >= sigma.values.size() || row >= mechanism.values.size()) continue;
            if (!sigma.present[row] || !mechanism.present[row]) continue;
            sum += sigma.values[row];
            ++n;
        }
        if (n == 0) return Mat3::Constant(kNaN);
        return sum / static_cast<double>(n);
    }

    std::vector<double> tensorContractionSeries(const Mat3Series& sigma,
                                                const Mat3Series& mechanism,
                                                const std::vector<std::size_t>& rows,
                                                bool centeredSigma,
                                                bool cosineAlignment) const {
        std::vector<double> out(sigma.values.size(), kNaN);
        const Mat3 mean = centeredSigma ? pairedMeanSigmaTensor(sigma, mechanism, rows)
                                        : Mat3::Zero();
        if (centeredSigma && !mean.allFinite()) return out;
        for (std::size_t row : rows) {
            if (row >= sigma.values.size() || row >= mechanism.values.size()) continue;
            if (!sigma.present[row] || !mechanism.present[row]) continue;
            const Mat3 s = centeredSigma ? sigma.values[row] - mean : sigma.values[row];
            const Mat3& m = mechanism.values[row];
            double v = tensorDoubleDot(s, m);
            if (cosineAlignment) {
                const double sn = tensorFrobenius(s);
                const double mn = tensorFrobenius(m);
                if (!(sn > 0.0) || !(mn > 0.0)) continue;
                v /= (sn * mn);
            }
            if (finite(v)) out[row] = v;
        }
        return out;
    }

    std::vector<double> centeredScalarSeries(const std::vector<double>& values,
                                             const std::vector<std::size_t>& rows) const {
        std::vector<double> out(values.size(), kNaN);
        double sum = 0.0;
        int n = 0;
        for (std::size_t row : rows) {
            if (row >= values.size() || !finite(values[row])) continue;
            sum += values[row];
            ++n;
        }
        if (n == 0) return out;
        const double mean = sum / static_cast<double>(n);
        for (std::size_t row : rows) {
            if (row < values.size() && finite(values[row])) out[row] = values[row] - mean;
        }
        return out;
    }

    // -- EFG |V_zz| and eta series from the (preferred) EFG T2.
    std::vector<double> efgEigenSeries(bool wantVzz) const {
        std::vector<double> out(efg_aimnet2_.values.size(), kNaN);
        for (std::size_t i = 0; i < efg_aimnet2_.values.size(); ++i) {
            const EfgValue* src = nullptr;
            if (efg_aimnet2_.present[i]) src = &efg_aimnet2_.values[i];
            else if (i < shielding_mopac_coulomb_.values.size() && shielding_mopac_coulomb_.present[i])
                src = &shielding_mopac_coulomb_.values[i];
            if (!src) continue;
            const EfgEigen e = efgEigen(src->t2);
            out[i] = wantVzz ? e.v_zz_abs : e.eta;
        }
        return out;
    }

    // -- signed X-H bond field E_|| series (§4.8): the field projected onto the
    // unit X-H bond vector, per frame. Signed (negative for a donor). Only for
    // an exchangeable/polar H with a single heavy parent.
    std::vector<double> signedBondFieldSeries() const {
        std::vector<double> out(cadence_.stepCount(), kNaN);
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return out;
        const model::QtAtom& a = body_.run.protein->atom(atom_);
        if (a.element != model::Element::H || a.parentAtomIndex < 0) return out;
        for (std::size_t step = 0; step < cadence_.stepCount(); ++step) {
            if (!field_mopac_.present[step]) continue;
            const auto h = coordAt(static_cast<int32_t>(atom_), step);
            const auto x = coordAt(a.parentAtomIndex, step);
            if (!h || !x) continue;
            const auto u = normalizeFrameVec(*h - *x);  // bond axis X->H
            if (!u) continue;
            out[step] = field_mopac_.values[step].dot(*u);  // signed E_||
        }
        return out;
    }

    static void appendChannel(SubspaceFamily& f,
                              const QString& name,
                              const std::vector<double>& values) {
        f.channels.push_back({name, values});
    }

    void appendMat3Channels(SubspaceFamily& f,
                            const QString& prefix,
                            const Mat3Series& src) const {
        for (int c = 0; c < 6; ++c) {
            const QString comp = QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(c)]);
            appendChannel(f, QStringLiteral("%1.%2").arg(prefix, comp),
                          molCompFrom(src, static_cast<MolComp>(c)));
        }
    }

    std::vector<double> intSeriesAsDouble(const IntSeries& src) const {
        std::vector<double> out(src.values.size(), kNaN);
        for (std::size_t i = 0; i < src.values.size(); ++i) {
            if (src.values[i] != IntSeries::kMissing)
                out[i] = static_cast<double>(src.values[i]);
        }
        return out;
    }

    SubspaceFamily fieldMopacFamily() const {
        SubspaceFamily f{QStringLiteral("field_mopac"), {}};
        appendChannel(f, QStringLiteral("field.mopac_coulomb.abs_E"), field_abs_E_mopac_.values);
        appendChannel(f, QStringLiteral("field.mopac_coulomb.E2"), field_abs_E2_mopac_.values);
        appendChannel(f, QStringLiteral("field.mopac_coulomb.E_parallel_mol_z"), field_E_z_mopac_.values);
        appendChannel(f, QStringLiteral("field.mopac_coulomb.mol_x"), componentSeries(field_mol_mopac_, 0));
        appendChannel(f, QStringLiteral("field.mopac_coulomb.mol_y"), componentSeries(field_mol_mopac_, 1));
        appendChannel(f, QStringLiteral("field.mopac_coulomb.mol_z"), componentSeries(field_mol_mopac_, 2));
        return f;
    }

    SubspaceFamily fieldExternalFamily() const {
        SubspaceFamily f{QStringLiteral("field_external"), {}};
        appendChannel(f, QStringLiteral("field.apbs.abs_E"), field_abs_E_apbs_.values);
        appendChannel(f, QStringLiteral("field.apbs.E2"), field_abs_E2_apbs_.values);
        appendChannel(f, QStringLiteral("field.apbs.E_parallel_mol_z"), field_E_z_apbs_.values);
        appendChannel(f, QStringLiteral("field.apbs.mol_x"), componentSeries(field_mol_apbs_, 0));
        appendChannel(f, QStringLiteral("field.apbs.mol_y"), componentSeries(field_mol_apbs_, 1));
        appendChannel(f, QStringLiteral("field.apbs.mol_z"), componentSeries(field_mol_apbs_, 2));
        appendChannel(f, QStringLiteral("field.charge_ff14sb.abs_E"), field_abs_E_charge_ff14sb_.values);
        appendChannel(f, QStringLiteral("field.charge_ff14sb.E2"), field_abs_E2_charge_ff14sb_.values);
        appendChannel(f, QStringLiteral("field.charge_ff14sb.E_parallel_mol_z"), field_E_z_charge_ff14sb_.values);
        appendChannel(f, QStringLiteral("field.charge_ff14sb.mol_x"), componentSeries(field_mol_charge_ff14sb_, 0));
        appendChannel(f, QStringLiteral("field.charge_ff14sb.mol_y"), componentSeries(field_mol_charge_ff14sb_, 1));
        appendChannel(f, QStringLiteral("field.charge_ff14sb.mol_z"), componentSeries(field_mol_charge_ff14sb_, 2));
        return f;
    }

    SubspaceFamily fieldSourcesFamily() const {
        SubspaceFamily f{QStringLiteral("field_sources"), {}};
        for (const SubspaceChannel& ch : fieldMopacFamily().channels) f.channels.push_back(ch);
        for (const SubspaceChannel& ch : fieldExternalFamily().channels) f.channels.push_back(ch);
        return f;
    }

    SubspaceFamily efgNodeFamily() const {
        SubspaceFamily f{QStringLiteral("efg_node"), {}};
        appendChannel(f, QStringLiteral("efg.aimnet2.abs_T2"), efg_abs_aimnet2_.values);
        appendChannel(f, QStringLiteral("efg.apbs.abs_T2"), efg_abs_apbs_.values);
        appendChannel(f, QStringLiteral("efg.mopac_coulomb.abs_T2"), shielding_abs_mopac_coulomb_.values);
        appendChannel(f, QStringLiteral("efg.best.Vzz_abs"), efgEigenSeries(true));
        appendChannel(f, QStringLiteral("efg.best.eta"), efgEigenSeries(false));
        appendMat3Channels(f, QStringLiteral("efg.aimnet2.mol"), efg_mol_aimnet2_);
        appendMat3Channels(f, QStringLiteral("efg.apbs.mol"), efg_mol_apbs_);
        appendMat3Channels(f, QStringLiteral("efg.mopac_coulomb.mol"), shielding_mol_mopac_);
        return f;
    }

    SubspaceFamily efgAimnetApbsFamily() const {
        SubspaceFamily f{QStringLiteral("efg_aimnet2_apbs"), {}};
        appendChannel(f, QStringLiteral("efg.aimnet2.abs_T2"), efg_abs_aimnet2_.values);
        appendChannel(f, QStringLiteral("efg.apbs.abs_T2"), efg_abs_apbs_.values);
        appendMat3Channels(f, QStringLiteral("efg.aimnet2.mol"), efg_mol_aimnet2_);
        appendMat3Channels(f, QStringLiteral("efg.apbs.mol"), efg_mol_apbs_);
        return f;
    }

    SubspaceFamily efgMopacFamily() const {
        SubspaceFamily f{QStringLiteral("efg_mopac_coulomb"), {}};
        appendChannel(f, QStringLiteral("efg.mopac_coulomb.abs_T2"), shielding_abs_mopac_coulomb_.values);
        appendMat3Channels(f, QStringLiteral("efg.mopac_coulomb.mol"), shielding_mol_mopac_);
        return f;
    }

    SubspaceFamily ringTensorFamily(const QString& name,
                                    const std::vector<Mat3Series>& tensors) const {
        SubspaceFamily f{name, {}};
        for (std::size_t type = 0; type < tensors.size(); ++type)
            appendMat3Channels(f, QStringLiteral("%1.type%2").arg(name).arg(type), tensors[type]);
        return f;
    }

    SubspaceFamily bsTensorFamilyForType(int type) const {
        SubspaceFamily f{QStringLiteral("bs_tensor_type%1").arg(type), {}};
        if (type >= 0 && static_cast<std::size_t>(type) < bs_per_type_mol_.size())
            appendMat3Channels(f, QStringLiteral("bs.type%1").arg(type),
                               bs_per_type_mol_[static_cast<std::size_t>(type)]);
        return f;
    }

    SubspaceFamily hmTensorFamilyForType(int type) const {
        SubspaceFamily f{QStringLiteral("hm_tensor_type%1").arg(type), {}};
        if (type >= 0 && static_cast<std::size_t>(type) < hm_per_type_mol_.size())
            appendMat3Channels(f, QStringLiteral("hm.type%1").arg(type),
                               hm_per_type_mol_[static_cast<std::size_t>(type)]);
        return f;
    }

    SubspaceFamily localElectronicPopulationFamily() const {
        SubspaceFamily f{QStringLiteral("local_electronic_population"), {}};
        appendChannel(f, QStringLiteral("mopac.charge"), mopac_charge_.values);
        appendChannel(f, QStringLiteral("mopac.s_population"), mopac_s_pop_.values);
        appendChannel(f, QStringLiteral("mopac.p_population"), mopac_p_pop_.values);
        appendChannel(f, QStringLiteral("mopac.valency"), mopac_valency_.values);
        return f;
    }

    SubspaceFamily localElectronicBondingFamily() const {
        SubspaceFamily f{QStringLiteral("local_electronic_bonding"), {}};
        appendChannel(f, QStringLiteral("mopac.s_character"), mopac_s_character_.values);
        appendChannel(f, QStringLiteral("hybridisation.ordinal"), intSeriesAsDouble(hybridisation_));
        appendChannel(f, QStringLiteral("hybridisation.pi_character"), pi_character_.values);
        appendChannel(f, QStringLiteral("dominant_partner_wiberg"), partnerWibergSeries());
        return f;
    }

    SubspaceFamily localElectronicBundleFamily() const {
        SubspaceFamily f{QStringLiteral("local_electronic_bundle"), {}};
        for (const SubspaceChannel& ch : localElectronicPopulationFamily().channels) f.channels.push_back(ch);
        for (const SubspaceChannel& ch : localElectronicBondingFamily().channels) f.channels.push_back(ch);
        return f;
    }

    SubspaceFamily larsenHbondFamily() const {
        SubspaceFamily f{QStringLiteral("hbond_larsen"), {}};
        appendChannel(f, QStringLiteral("larsen_hbond.T0"), componentSeriesT0(larsen_hbond_shielding_));
        appendChannel(f, QStringLiteral("hbond.count"), hbond_count_.values);
        appendChannel(f, QStringLiteral("hbond.nearest_dir.mol_x"), componentSeries(hbond_nearest_dir_mol_, 0));
        appendChannel(f, QStringLiteral("hbond.nearest_dir.mol_y"), componentSeries(hbond_nearest_dir_mol_, 1));
        appendChannel(f, QStringLiteral("hbond.nearest_dir.mol_z"), componentSeries(hbond_nearest_dir_mol_, 2));
        appendMat3Channels(f, QStringLiteral("larsen_hbond.mol"), larsen_hbond_shielding_mol_);
        return f;
    }

    std::vector<double> mcAxialSeries(const Mat3Series& src) const {
        std::vector<double> xx = molCompFrom(src, MolComp::xx);
        std::vector<double> yy = molCompFrom(src, MolComp::yy);
        std::vector<double> zz = molCompFrom(src, MolComp::zz);
        std::vector<double> out(std::max({xx.size(), yy.size(), zz.size()}), kNaN);
        for (std::size_t i = 0; i < out.size(); ++i) {
            const double x = i < xx.size() ? xx[i] : kNaN;
            const double y = i < yy.size() ? yy[i] : kNaN;
            const double z = i < zz.size() ? zz[i] : kNaN;
            if (finite(x) && finite(y) && finite(z)) out[i] = z - 0.5 * (x + y);
        }
        return out;
    }

    std::vector<double> mcRhombicSeries(const Mat3Series& src) const {
        std::vector<double> xx = molCompFrom(src, MolComp::xx);
        std::vector<double> yy = molCompFrom(src, MolComp::yy);
        std::vector<double> xy = molCompFrom(src, MolComp::xy);
        std::vector<double> xz = molCompFrom(src, MolComp::xz);
        std::vector<double> yz = molCompFrom(src, MolComp::yz);
        std::vector<double> out(std::max({xx.size(), yy.size(), xy.size(), xz.size(), yz.size()}), kNaN);
        for (std::size_t i = 0; i < out.size(); ++i) {
            const double x = i < xx.size() ? xx[i] : kNaN;
            const double y = i < yy.size() ? yy[i] : kNaN;
            const double a = i < xy.size() ? xy[i] : kNaN;
            const double b = i < xz.size() ? xz[i] : kNaN;
            const double c = i < yz.size() ? yz[i] : kNaN;
            if (finite(x) && finite(y) && finite(a) && finite(b) && finite(c))
                out[i] = std::sqrt((x - y) * (x - y) + 4.0 * (a * a + b * b + c * c));
        }
        return out;
    }

    SubspaceFamily mcAxialFamily() const {
        SubspaceFamily f{QStringLiteral("mcconnell_axial"), {}};
        for (std::size_t i = 0; i < mc_tensor_mol_series_.size(); ++i)
            appendChannel(f, QStringLiteral("%1.axial")
                                .arg(QString::fromLatin1(kMcTensorFields[i].key)),
                          mcAxialSeries(mc_tensor_mol_series_[i]));
        return f;
    }

    SubspaceFamily mcRhombicFamily() const {
        SubspaceFamily f{QStringLiteral("mcconnell_rhombic"), {}};
        for (std::size_t i = 0; i < mc_tensor_mol_series_.size(); ++i)
            appendChannel(f, QStringLiteral("%1.rhombic")
                                .arg(QString::fromLatin1(kMcTensorFields[i].key)),
                          mcRhombicSeries(mc_tensor_mol_series_[i]));
        return f;
    }

    QJsonObject subspaceCompareObject(const QString& readId,
                                      const SubspaceFamily& a,
                                      const SubspaceFamily& b,
                                      const std::vector<std::size_t>& rows) const {
        const SubspaceCompareResult r = CompareSubspaces(a, b, rows);
        QJsonObject o = SubspaceCompareJson(r);
        o.insert(QStringLiteral("read_id"), readId);
        o.insert(QStringLiteral("family_a"), a.name);
        o.insert(QStringLiteral("family_b"), b.name);
        return o;
    }

    bool isDonorH() const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return false;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);
        if (a.element != model::Element::H || a.parentAtomIndex < 0
            || static_cast<std::size_t>(a.parentAtomIndex) >= p.atomCount())
            return false;
        const model::Element parentEl = p.atom(static_cast<std::size_t>(a.parentAtomIndex)).element;
        return a.isExchangeable || parentEl == model::Element::N || parentEl == model::Element::O;
    }

    bool withinRingCutoff() const {
        // a heavy/H atom that has any ring-current kernel signal at all.
        for (std::size_t i = 0; i < bs_per_type_T2_.values.size(); ++i)
            if (bs_per_type_T2_.present[i]) return true;
        return false;
    }

    bool hasMcPartner() const {
        for (const TensorSeries& s : mc_tensor_series_) {
            for (std::size_t i = 0; i < s.values.size(); ++i) {
                if (!s.present[i]) continue;
                for (double v : s.values[i].T2)
                    if (finite(v) && v != 0.0) return true;
            }
        }
        return false;
    }

    bool hasMappedBond() const {
        for (const auto& item : bond_series_)
            if (static_bonds_.count(item.first) != 0) return true;
        return false;
    }

    bool hasFramedAtom() const {
        for (const MolecularFrameValue& v : molecular_frame_.values)
            if (v.valid) return true;
        return false;
    }

    int graphStratum() const {
        // coarse hops-bucket on nearest_ring_atom_hops (the spatial-locality
        // moderator); a single static value per atom (does not multiply contexts).
        const int hops = minGraphHopsTo([](const model::QtAtom& x) { return x.IsInAnyRing(); });
        if (hops < 0) return -1;
        if (hops == 0) return 0;     // is a ring atom
        if (hops <= 2) return 1;     // proximal
        return 2;                    // distal
    }

    // Does the catalogue row's membership predicate hold for this atom (§2.2A)?
    bool membershipHolds(CatalogueId id) const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return false;
        const model::QtAtom& a = body_.run.protein->atom(atom_);
        const bool heavy = a.element != model::Element::H;
        const bool framed = hasFramedAtom();
        switch (id) {
        case CatalogueId::FIELD_ISO:
        case CatalogueId::FIELD_SPAN:
        case CatalogueId::STRUCT:
            return true;  // all atoms
        case CatalogueId::FIELD_MOLCOMP:
            return framed;
        case CatalogueId::FIELD_SIGNED_BOND:
            return isDonorH();
        case CatalogueId::EFG_ISO:
            return heavy;
        case CatalogueId::EFG_MOLCOMP:
            return heavy && framed;
        case CatalogueId::EFG_EIGEN:
            return heavy && framed;  // framed heavy atom (sp2 groups carry the cleanest eigen)
        case CatalogueId::RING_HEAVY:
            return heavy && withinRingCutoff();
        case CatalogueId::RING_H:
            return a.element == model::Element::H && withinRingCutoff();
        case CatalogueId::MCCONNELL:
            return hasMcPartner();
        case CatalogueId::SCHAR_ISO:
            return a.element == model::Element::N || a.element == model::Element::C;
        case CatalogueId::WIBERG:
            return hasMappedBond();
        case CatalogueId::BONDLEN_ISO:
            return a.element == model::Element::S || !static_bonds_.empty();
        case CatalogueId::DIHEDRAL:
            return hasAnyDihedral();
        case CatalogueId::HBOND:
            return hasHbondSignal();
        }
        return false;
    }

    bool hasAnyDihedral() const {
        for (const ScalarSeries* s : {&phi_, &psi_, &chi1_, &chi2_, &chi3_, &chi4_})
            for (double v : s->values)
                if (finite(v)) return true;
        return false;
    }

    bool hasHbondSignal() const {
        for (double v : hbond_count_.values)
            if (finite(v) && v > 0.0) return true;
        return false;
    }

    // The candidate (driver-series, sigma-target-series, sigma_target label,
    // driver component label) pairs a catalogue row evaluates. The builder
    // groups these by sigma_target: one primary source per target, every
    // physically-relevant target carried in targets[].
    struct Candidate {
        std::vector<double> driver;
        std::vector<double> sigma;
        QString sigma_target;     // iso | invariant:span | molcomp:yy | sigma_total.molcomp:yy | ...
        QString driver_channel;   // the channel + component
    };

    void appendTensorComponentCandidates(std::vector<Candidate>& out,
                                         const Mat3Series& mechanism,
                                         const QString& channel,
                                         bool crossProduct,
                                         bool includePara) const {
        for (int dc = 0; dc < 6; ++dc) {
            const auto dcomp = static_cast<MolComp>(dc);
            const std::vector<double> driver = molCompFrom(mechanism, dcomp);
            const QString dcName = QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(dc)]);
            for (int sc = 0; sc < 6; ++sc) {
                if (!crossProduct && sc != dc) continue;
                const auto scomp = static_cast<MolComp>(sc);
                const QString scName = QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(sc)]);
                out.push_back({driver, sigmaMolCompSeries(sigma_mol_total_, scomp),
                               QStringLiteral("sigma_total.molcomp:%1").arg(scName),
                               QStringLiteral("%1|molcomp:%2").arg(channel, dcName)});
                if (includePara) {
                    out.push_back({driver, sigmaMolCompSeries(sigma_mol_para_, scomp),
                                   QStringLiteral("sigma_para.molcomp:%1").arg(scName),
                                   QStringLiteral("%1|molcomp:%2").arg(channel, dcName)});
                }
            }
        }
    }

    void appendRingTensorCandidates(std::vector<Candidate>& out, bool crossProduct) const {
        for (int type = 0; type < 8; ++type) {
            appendTensorComponentCandidates(
                out, bs_per_type_mol_[static_cast<std::size_t>(type)],
                QStringLiteral("bs.type%1").arg(type), crossProduct, true);
            appendTensorComponentCandidates(
                out, hm_per_type_mol_[static_cast<std::size_t>(type)],
                QStringLiteral("hm.type%1").arg(type), crossProduct, true);
        }
    }

    void appendMcTensorCandidates(std::vector<Candidate>& out, bool crossProduct) const {
        for (std::size_t i = 0; i < mc_tensor_mol_series_.size(); ++i) {
            appendTensorComponentCandidates(out, mc_tensor_mol_series_[i],
                                            QString::fromLatin1(kMcTensorFields[i].key),
                                            crossProduct, true);
        }
    }

    std::vector<double> summedMolCompFrom(const std::vector<Mat3Series>& sources, MolComp comp) const {
        std::vector<double> out(cadence_.stepCount(), kNaN);
        for (const Mat3Series& src : sources) {
            const std::vector<double> c = molCompFrom(src, comp);
            for (std::size_t i = 0; i < out.size() && i < c.size(); ++i) {
                if (!finite(c[i])) continue;
                out[i] = finite(out[i]) ? out[i] + c[i] : c[i];
            }
        }
        return out;
    }

    std::vector<double> sumSeries(const std::vector<double>& a, const std::vector<double>& b) const {
        std::vector<double> out(std::max(a.size(), b.size()), kNaN);
        for (std::size_t i = 0; i < out.size(); ++i) {
            const double av = i < a.size() ? a[i] : kNaN;
            const double bv = i < b.size() ? b[i] : kNaN;
            if (finite(av) && finite(bv)) out[i] = av + bv;
            else if (finite(av)) out[i] = av;
            else if (finite(bv)) out[i] = bv;
        }
        return out;
    }

    void appendCatalogueSignedComponentCandidate(std::vector<Candidate>& out,
                                                 const std::vector<double>& driver,
                                                 const QString& channel) const {
        out.push_back({driver, sigmaMolCompSeries(sigma_mol_total_, MolComp::zz),
                       QStringLiteral("sigma_total.molcomp:zz"), channel});
        out.push_back({driver, sigmaMolCompSeries(sigma_mol_para_, MolComp::zz),
                       QStringLiteral("sigma_para.molcomp:zz"), channel});
    }

    std::vector<Candidate> candidatesFor(CatalogueId id) const {
        std::vector<Candidate> out;
        const std::vector<double> iso = sigmaIsoSeries();
        switch (id) {
        case CatalogueId::FIELD_ISO:
            out.push_back({fieldAbsSeries(), iso, QStringLiteral("iso"),
                           QStringLiteral("field.mopac_coulomb|abs_E")});
            break;
        case CatalogueId::FIELD_SPAN:
            out.push_back({fieldAbsSeries(), sigmaInvariantSeries(QStringLiteral("span")),
                           QStringLiteral("invariant:span"), QStringLiteral("field.mopac_coulomb|abs_E")});
            break;
        case CatalogueId::FIELD_MOLCOMP:
            // cross-product (field axis x/y/z) x (sigma molcomp) so a
            // cross-component coupling (field_y -> sigma_xy etc.) is found, not
            // only same-component pairs. Sources are summarized per sigma target.
            for (int axis = 0; axis < 3; ++axis) {
                const std::vector<double> fieldComp = fieldMolCompSeries(axis);
                const QString axisName = QStringLiteral("xyz").mid(axis, 1);
                for (int sc = 0; sc < 6; ++sc) {
                    const auto smc = static_cast<MolComp>(sc);
                    out.push_back({fieldComp, sigmaMolCompSeries(sigma_mol_total_, smc),
                                   QStringLiteral("molcomp:%1").arg(QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(sc)])),
                                   QStringLiteral("field_mol|%1").arg(axisName)});
                }
            }
            break;
        case CatalogueId::FIELD_SIGNED_BOND:
            out.push_back({signedBondFieldSeries(), iso, QStringLiteral("iso"),
                           QStringLiteral("field.mopac_coulomb|E_parallel_XH")});
            break;
        case CatalogueId::EFG_ISO:
            out.push_back({efgAbsSeries(), iso, QStringLiteral("iso"),
                           QStringLiteral("efg|abs_T2")});
            break;
        case CatalogueId::EFG_MOLCOMP: {
            // cross-product (efg molcomp) x (sigma molcomp), over BOTH EFG
            // sources: the AIMNet2 EFG and the MOPAC-Coulomb EFG (the inventory's
            // strongest heavy-atom |EFG| correlate, now shielding_mol). The
            // keystone cross-component coupling efg_zz -> sigma_yy (LYS30 C) and
            // the ASP7 CG efg_xz -> sigma_xz both live here; primary source is
            // summarized per sigma target.
            struct EfgSrc { const Mat3Series* s; const char* tag; };
            const std::array<EfgSrc, 2> srcs = {{
                {&efg_mol_aimnet2_, "efg_mol.aimnet2"},
                {&shielding_mol_mopac_, "efg_mol.mopac_coulomb"},
            }};
            for (const EfgSrc& src : srcs) {
                for (int dc = 0; dc < 6; ++dc) {
                    const std::vector<double> efgComp = molCompFrom(*src.s, static_cast<MolComp>(dc));
                    const QString dcName = QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(dc)]);
                    for (int sc = 0; sc < 6; ++sc) {
                        const auto smc = static_cast<MolComp>(sc);
                        out.push_back({efgComp, sigmaMolCompSeries(sigma_mol_total_, smc),
                                       QStringLiteral("molcomp:%1").arg(QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(sc)])),
                                       QStringLiteral("%1|%2").arg(QString::fromLatin1(src.tag), dcName)});
                    }
                }
            }
            break;
        }
        case CatalogueId::EFG_EIGEN:
            out.push_back({efgEigenSeries(true), sigmaInvariantSeries(QStringLiteral("span")),
                           QStringLiteral("invariant:span"), QStringLiteral("efg|V_zz_abs")});
            out.push_back({efgEigenSeries(false), sigmaInvariantSeries(QStringLiteral("eta_H")),
                           QStringLiteral("invariant:eta_H"), QStringLiteral("efg|eta")});
            break;
        case CatalogueId::RING_HEAVY:
            appendCatalogueSignedComponentCandidate(
                out,
                sumSeries(summedMolCompFrom(bs_per_type_mol_, MolComp::zz),
                          summedMolCompFrom(hm_per_type_mol_, MolComp::zz)),
                QStringLiteral("ring_tensor_signed_sum|molcomp:zz"));
            break;
        case CatalogueId::RING_H:
            appendCatalogueSignedComponentCandidate(
                out,
                sumSeries(summedMolCompFrom(bs_per_type_mol_, MolComp::zz),
                          summedMolCompFrom(hm_per_type_mol_, MolComp::zz)),
                QStringLiteral("ring_tensor_signed_sum|molcomp:zz"));
            break;
        case CatalogueId::MCCONNELL:
            appendCatalogueSignedComponentCandidate(
                out,
                summedMolCompFrom(mc_tensor_mol_series_, MolComp::zz),
                QStringLiteral("mc_tensor_signed_sum|molcomp:zz"));
            break;
        case CatalogueId::SCHAR_ISO:
            out.push_back({mopac_s_character_.values, iso, QStringLiteral("iso"),
                           QStringLiteral("mopac.s_character")});
            break;
        case CatalogueId::WIBERG:
            // WIBERG is a structured accumulator object (§2.1), not a
            // responses[] row. buildWibergJson() evaluates its member series.
            break;
        case CatalogueId::BONDLEN_ISO:
            out.push_back({disulfideOrFirstBondLengthSeries(), iso, QStringLiteral("iso"),
                           QStringLiteral("bond.length")});
            break;
        case CatalogueId::DIHEDRAL:
            // the dihedral with the most signal (chi-first), correlated as
            // |sin| projection to bounded scalar sigma targets. Use eta^2 over
            // wells in the response_law; the linear pair here uses cos(dihedral).
            out.push_back({dihedralCosSeries(), iso, QStringLiteral("iso"),
                           QStringLiteral("dihedral|cos")});
            out.push_back({dihedralCosSeries(), sigmaInvariantSeries(QStringLiteral("aniso")),
                           QStringLiteral("invariant:aniso"), QStringLiteral("dihedral|cos")});
            out.push_back({dihedralCosSeries(), sigmaInvariantSeries(QStringLiteral("span")),
                           QStringLiteral("invariant:span"), QStringLiteral("dihedral|cos")});
            out.push_back({dihedralCosSeries(), sigmaInvariantSeries(QStringLiteral("eta_H")),
                           QStringLiteral("invariant:eta_H"), QStringLiteral("dihedral|cos")});
            out.push_back({dihedralCosSeries(), sigmaInvariantSeries(QStringLiteral("frobenius")),
                           QStringLiteral("invariant:frobenius"), QStringLiteral("dihedral|cos")});
            out.push_back({dihedralCosSeries(), sigmaT1FractionSeries(),
                           QStringLiteral("shape:t1_fraction"), QStringLiteral("dihedral|cos")});
            if (hasFramedAtom()) {
                for (int sc = 0; sc < 6; ++sc) {
                    const auto smc = static_cast<MolComp>(sc);
                    const QString scName = QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(sc)]);
                    out.push_back({dihedralCosSeries(), sigmaMolCompSeries(sigma_mol_total_, smc),
                                   QStringLiteral("molcomp:%1").arg(scName),
                                   QStringLiteral("dihedral|cos")});
                }
            }
            break;
        case CatalogueId::HBOND:
            out.push_back({hbond_count_.values, iso, QStringLiteral("iso"),
                           QStringLiteral("hbond.count")});
            break;
        case CatalogueId::STRUCT:
            out.push_back({sasa_.values, iso, QStringLiteral("iso"),
                           QStringLiteral("sasa")});
            break;
        }
        return out;
    }


    bool hasFiniteSample(const ScalarSeries& s) const {
        for (double v : s.values)
            if (finite(v)) return true;
        return false;
    }

    const ScalarSeries* rotamerRelevantDihedralSeries() const {
        if (body_.run.protein && atom_ < body_.run.protein->atomCount()) {
            const model::QtProtein& p = *body_.run.protein;
            const model::QtAtom& a = p.atom(atom_);
            if (validResidue(p, a.residueIndex)) {
                const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
                const QString atomName = p.atomLabel(atom_, model::NamingConvention::Iupac);
                if (r.aminoAcid == model::AminoAcid::MET && atomName == QStringLiteral("SD")
                    && hasFiniteSample(chi3_))
                    return &chi3_;
            }
        }

        // Physical side-chain priority, not "most populated": use the deepest
        // available chi first, then backbone torsions as the fallback context.
        for (const ScalarSeries* s : {&chi4_, &chi3_, &chi2_, &chi1_, &phi_, &psi_})
            if (hasFiniteSample(*s)) return s;
        return nullptr;
    }

    std::vector<double> dihedralCosSeries() const {
        const ScalarSeries* best = rotamerRelevantDihedralSeries();
        std::vector<double> out(cadence_.stepCount(), kNaN);
        if (!best) return out;
        for (std::size_t i = 0; i < best->values.size(); ++i)
            if (finite(best->values[i])) out[i] = std::cos(best->values[i]);
        return out;
    }

    // The rotamer-relevant dihedral's circular-well series (for eta^2 in DIHEDRAL).
    std::vector<int> rotamerRelevantDihedralWells() const {
        const ScalarSeries* best = rotamerRelevantDihedralSeries();
        if (!best) return std::vector<int>(cadence_.stepCount(), -1);
        return circularWells(best->values);
    }

    std::vector<double> partnerWibergSeries() const {
        // the partner-bond Wiberg order with the most present frames.
        const ScalarSeries* best = nullptr;
        int bestN = 0;
        for (const auto& item : bond_series_) {
            if (static_bonds_.count(item.first) == 0) continue;
            int n = 0;
            for (double v : item.second.values) if (finite(v)) ++n;
            if (n > bestN) { bestN = n; best = &item.second; }
        }
        return best ? best->values : std::vector<double>(cadence_.stepCount(), kNaN);
    }

    std::vector<double> disulfideOrFirstBondLengthSeries() const {
        // prefer a disulfide S-S length; else the first static bond length.
        const ScalarSeries* chosen = nullptr;
        for (const auto& item : static_bonds_) {
            if (item.second.category == model::BondCategory::Disulfide) {
                const auto it = bond_length_series_.find(item.first);
                if (it != bond_length_series_.end()) { chosen = &it->second; break; }
            }
        }
        if (!chosen && !bond_length_series_.empty()) chosen = &bond_length_series_.begin()->second;
        return chosen ? chosen->values : std::vector<double>(cadence_.stepCount(), kNaN);
    }

    std::uint8_t contactedClassMaskAtStep(std::size_t step) const {
        if (step >= contacted_class_membership_.size()) return 0;
        return contacted_class_membership_[step];
    }

    QJsonObject contactedClassMembershipSeriesJson() const {
        QJsonObject root;
        QJsonArray classes;
        QJsonArray ords;
        QJsonObject presence;
        for (ContactedClass c : kPresentContactClasses) {
            const QString name = contactedClassName(c);
            classes.append(name);
            ords.append(static_cast<int>(c));
            QJsonArray values;
            const std::uint8_t bit = contactedClassBit(c);
            for (std::uint8_t mask : contacted_class_membership_)
                values.append((mask & bit) != 0);
            presence.insert(name, values);
        }
        root.insert(QStringLiteral("classes"), classes);
        root.insert(QStringLiteral("class_ords"), ords);
        root.insert(QStringLiteral("presence_by_class"), presence);
        return root;
    }

    QJsonArray contactedClassChangeEventsJson() const {
        QJsonArray events;
        bool havePrev = false;
        std::uint8_t prev = 0;
        for (std::size_t i = 0; i < contacted_class_membership_.size(); ++i) {
            const std::uint8_t cur = contacted_class_membership_[i];
            if (havePrev && cur != prev) {
                const std::uint8_t entered = static_cast<std::uint8_t>(cur & ~prev);
                const std::uint8_t exited = static_cast<std::uint8_t>(prev & ~cur);
                QJsonObject e;
                e.insert(QStringLiteral("step_row"), static_cast<int>(i));
                e.insert(QStringLiteral("from_ord"), contactedClassSingleOrd(prev));
                e.insert(QStringLiteral("from"), contactedClassMaskName(prev));
                e.insert(QStringLiteral("to_ord"), contactedClassSingleOrd(cur));
                e.insert(QStringLiteral("to"), contactedClassMaskName(cur));
                e.insert(QStringLiteral("from_ords"), contactedClassMaskOrdsJson(prev));
                e.insert(QStringLiteral("from_classes"), contactedClassMaskNamesJson(prev));
                e.insert(QStringLiteral("to_ords"), contactedClassMaskOrdsJson(cur));
                e.insert(QStringLiteral("to_classes"), contactedClassMaskNamesJson(cur));
                e.insert(QStringLiteral("entered_ords"), contactedClassMaskOrdsJson(entered));
                e.insert(QStringLiteral("entered_classes"), contactedClassMaskNamesJson(entered));
                e.insert(QStringLiteral("exited_ords"), contactedClassMaskOrdsJson(exited));
                e.insert(QStringLiteral("exited_classes"), contactedClassMaskNamesJson(exited));
                events.append(e);
            }
            havePrev = true;
            prev = cur;
        }
        return events;
    }

    void foldContactClass(std::size_t step,
                          const std::pair<int, model::AminoAcid>& contacted,
                          bool present,
                          double /*distance*/) {
        if (!present || step >= contacted_class_membership_.size() || contacted.first < 0) return;
        const ContactedClass c = classOf(contacted.second);
        if (c == ContactedClass::NONE) return;
        contacted_class_membership_[step] =
            static_cast<std::uint8_t>(contacted_class_membership_[step] | contactedClassBit(c));
    }

    // One characterized sigma target inside a catalogue channel. A channel may
    // summarize sources to a primary, but it carries every physically-relevant
    // sigma target in response.targets[].
    struct EvaluatedTarget {
        QString sigma_target;
        QString driver_channel;
        FrameTag frame;
        // response_law
        double slope = kNaN;
        double intercept = kNaN;
        double r = kNaN;
        double null_z = kNaN;
        double xi = kNaN;
        double pchip_mid = kNaN;
        double leverage = kNaN;
        QJsonObject coverage;
        // gauntlet
        bool delta_survives = false;
        double seg_min_r = kNaN;
        double circshift_p = kNaN;
        // lead_lag
        int best_lag = 0;
        double lead_r = kNaN;
        // dihedral eta^2 (only for DIHEDRAL)
        double eta_squared = kNaN;
        std::vector<double> driver_series;
        std::vector<double> sigma_series;
    };

    struct EvaluatedResponse {
        CatalogueId id;
        QString driver_channel;
        FrameTag frame;
        std::vector<EvaluatedTarget> targets;
        QStringList considerations;
    };

    FrameTag frameForTarget(CatalogueId id, const QString& sigmaTarget) const {
        if (sigmaTarget.startsWith(QStringLiteral("molcomp:"))
            || sigmaTarget.contains(QStringLiteral(".molcomp:")))
            return FrameTag::Molecular;
        if (sigmaTarget.startsWith(QStringLiteral("invariant:"))) return FrameTag::Invariant;
        return catalogueRow(id).frame;
    }

    double targetCharacterizationScore(CatalogueId id, const EvaluatedTarget& t) const {
        if (id == CatalogueId::DIHEDRAL
            && t.sigma_target == QStringLiteral("invariant:aniso")
            && finite(t.eta_squared))
            return t.eta_squared;
        return finite(t.r) ? std::abs(t.r) : -1.0;
    }

    const EvaluatedTarget* primaryTarget(const EvaluatedResponse& e) const {
        const EvaluatedTarget* best = nullptr;
        double bestScore = -1.0;
        for (const EvaluatedTarget& t : e.targets) {
            const double s = targetCharacterizationScore(e.id, t);
            if (s > bestScore) {
                bestScore = s;
                best = &t;
            }
        }
        return best;
    }

    EvaluatedTarget nullTargetFor(const Candidate& c, CatalogueId id) const {
        EvaluatedTarget t;
        t.sigma_target = c.sigma_target;
        t.driver_channel = c.driver_channel;
        t.frame = frameForTarget(id, c.sigma_target);
        t.driver_series = c.driver;
        t.sigma_series = c.sigma;
        return t;
    }

    std::optional<EvaluatedTarget> evaluateCandidateTarget(CatalogueId id,
                                                           const Candidate& c,
                                                           const std::vector<std::size_t>& rows) const {
        std::vector<double> x;
        std::vector<double> y;
        pairedOnSigmaRows(c.sigma, c.driver, rows, x, y);
        if (x.size() < 3) return std::nullopt;

        EvaluatedTarget t;
        t.sigma_target = c.sigma_target;
        t.driver_channel = c.driver_channel;
        t.frame = frameForTarget(id, c.sigma_target);
        t.driver_series = c.driver;
        t.sigma_series = c.sigma;
        const OlsResult fit = ols(x, y);
        t.slope = fit.slope;
        t.intercept = fit.intercept;
        t.r = pearsonR(x, y);
        const auto [nm, nsd, nz, cp] = nullShiftStats(x, y, fit.slope);
        (void)nm; (void)nsd;
        t.null_z = nz;
        t.xi = chatterjeeXi(x, y);
        t.pchip_mid = pchipMidpoint(x, y);
        t.leverage = leverageTop1(x, y);
        t.coverage = coverageJson(x);
        t.delta_survives = deltaSurvives(x, y, t.r);
        t.seg_min_r = segMinR(x, y);
        t.circshift_p = cp;
        const LeadLagResult ll = leadLag(x, y);
        t.best_lag = ll.best_lag;
        t.lead_r = ll.lead_r;
        if (id == CatalogueId::DIHEDRAL)
            t.eta_squared = etaSquaredByWell(c.sigma, rotamerRelevantDihedralWells(), rows);
        return t;
    }

    std::optional<EvaluatedResponse> evaluateCatalogueRowFromCandidates(
            CatalogueId id,
            const std::vector<Candidate>& cands,
            const std::vector<std::size_t>& rows) const {
        if (cands.empty()) return std::nullopt;

        EvaluatedResponse e;
        e.id = id;
        e.frame = catalogueRow(id).frame;

        std::vector<QString> order;
        std::map<QString, Candidate> firstCandidate;
        std::map<QString, EvaluatedTarget> bestByTarget;
        std::map<QString, QStringList> sourcesByTarget;

        for (const Candidate& c : cands) {
            if (!firstCandidate.count(c.sigma_target)) {
                firstCandidate.emplace(c.sigma_target, c);
                order.push_back(c.sigma_target);
            }
            if (!sourcesByTarget[c.sigma_target].contains(c.driver_channel))
                sourcesByTarget[c.sigma_target].append(c.driver_channel);

            const auto t = evaluateCandidateTarget(id, c, rows);
            if (!t) continue;
            const auto it = bestByTarget.find(c.sigma_target);
            if (it == bestByTarget.end()
                || targetCharacterizationScore(id, *t) > targetCharacterizationScore(id, it->second))
                bestByTarget[c.sigma_target] = *t;
        }

        for (const QString& targetName : order) {
            const auto hit = bestByTarget.find(targetName);
            if (hit != bestByTarget.end()) {
                e.targets.push_back(hit->second);
            } else {
                e.targets.push_back(nullTargetFor(firstCandidate.at(targetName), id));
            }
        }

        const EvaluatedTarget* primary = primaryTarget(e);
        if (primary) {
            e.driver_channel = primary->driver_channel;
        } else {
            e.driver_channel = QString::fromLatin1(catalogueRow(id).driver_family);
        }

        for (const EvaluatedTarget& t : e.targets) {
            const QStringList srcs = sourcesByTarget[t.sigma_target];
            for (const QString& src : srcs) {
                if (src == t.driver_channel) continue;
                e.considerations.append(QStringLiteral("candidate_source[%1]: %2")
                                            .arg(t.sigma_target, src));
            }
        }

        return e;
    }

    std::optional<EvaluatedResponse> evaluateCatalogueRow(CatalogueId id,
                                                          const std::vector<std::size_t>& rows) const {
        return evaluateCatalogueRowFromCandidates(id, candidatesFor(id), rows);
    }

    std::optional<int32_t> residueIndexForNumber(int residueNumberValue,
                                                 model::AminoAcid aa = model::AminoAcid::Unknown) const {
        if (!body_.run.protein) return std::nullopt;
        const model::QtProtein& p = *body_.run.protein;
        for (std::size_t ri = 0; ri < p.residueCount(); ++ri) {
            const model::QtResidue& r = p.residue(ri);
            if (r.address.residueNumber != residueNumberValue) continue;
            if (aa != model::AminoAcid::Unknown && r.aminoAcid != aa) continue;
            return static_cast<int32_t>(ri);
        }
        return std::nullopt;
    }

    // The salt-bridge partner = the NEAREST acidic contact, chosen by geometry
    // (smallest mean donor↔carboxylate-O distance) from the atom's own contact
    // set. No residue is pinned — the selection carries no 1P9J knowledge, so the
    // pass is protein-agnostic (feed any protein, it finds that protein's partner).
    std::optional<int32_t> acidicContactResidueIndex() const {
        std::optional<int32_t> best;
        double bestMean = 0.0;
        for (const auto& r : contacted_residues_) {
            const auto aa = static_cast<model::AminoAcid>(r.second);
            if (classOf(aa) != ContactedClass::ACIDIC) continue;
            const auto ri = residueIndexForNumber(r.first, aa);
            if (!ri) continue;
            const std::vector<double> d = saltBridgeDistanceSeriesFor(*ri);
            double sum = 0.0;
            std::size_t n = 0;
            for (double v : d) if (finite(v)) { sum += v; ++n; }
            if (n == 0) continue;
            const double mean = sum / static_cast<double>(n);
            if (!best || mean < bestMean) { bestMean = mean; best = ri; }
        }
        return best;
    }

    std::vector<int32_t> guanidiniumDonorAtoms(int32_t residueIndex) const {
        std::vector<int32_t> atoms;
        const std::array<QString, 3> names = {
            QStringLiteral("NH1"), QStringLiteral("NH2"), QStringLiteral("NE")
        };
        for (const QString& name : names) {
            if (const auto ai = atomByResidueName(residueIndex, name)) atoms.push_back(*ai);
        }
        return atoms;
    }

    std::vector<int32_t> acidicSidechainAtoms(int32_t residueIndex) const {
        std::vector<int32_t> atoms;
        if (!body_.run.protein || !validResidue(*body_.run.protein, residueIndex)) return atoms;
        const model::QtResidue& r =
            body_.run.protein->residue(static_cast<std::size_t>(residueIndex));
        const std::array<QString, 2> asp = {QStringLiteral("OD1"), QStringLiteral("OD2")};
        const std::array<QString, 2> glu = {QStringLiteral("OE1"), QStringLiteral("OE2")};
        const auto& names = r.aminoAcid == model::AminoAcid::GLU ? glu : asp;
        for (const QString& name : names) {
            if (const auto ai = atomByResidueName(residueIndex, name)) atoms.push_back(*ai);
        }
        if (!atoms.empty()) return atoms;
        for (int32_t ai : r.atomIndices) {
            if (!validAtomIndex(ai)) continue;
            if (body_.run.protein->atom(static_cast<std::size_t>(ai)).element != model::Element::H)
                atoms.push_back(ai);
        }
        return atoms;
    }

    // Per-step min donor↔carboxylate-O distance to a GIVEN acidic residue. The
    // donor is the atom's guanidinium N's (NH1/NH2/NE) for an ARG, else the atom
    // itself. Parameterised by the acidic residue so the partner selection can
    // score every acidic candidate — no residue pinned.
    std::vector<double> saltBridgeDistanceSeriesFor(int32_t acidRi) const {
        std::vector<double> out(cadence_.stepCount(), kNaN);
        const std::vector<int32_t> acidAtoms = acidicSidechainAtoms(acidRi);
        if (acidAtoms.empty()) return out;
        std::vector<int32_t> donorAtoms;
        if (body_.run.protein && atom_ < body_.run.protein->atomCount()) {
            const model::QtProtein& p = *body_.run.protein;
            const model::QtAtom& a = p.atom(atom_);
            if (validResidue(p, a.residueIndex)) {
                const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
                if (r.aminoAcid == model::AminoAcid::ARG)
                    donorAtoms = guanidiniumDonorAtoms(a.residueIndex);
            }
        }
        if (donorAtoms.empty()) donorAtoms.push_back(static_cast<int32_t>(atom_));
        for (std::size_t step = 0; step < cadence_.stepCount(); ++step) {
            double minR = kNaN;
            for (int32_t donor : donorAtoms) {
                const auto dpos = coordAt(donor, step);
                if (!dpos) continue;
                for (int32_t ai : acidAtoms) {
                    const auto apos = coordAt(ai, step);
                    if (!apos) continue;
                    const double r = (*dpos - *apos).norm();
                    if (!finite(minR) || r < minR) minR = r;
                }
            }
            if (finite(minR)) out[step] = minR;
        }
        return out;
    }

    std::vector<double> saltBridgeDistanceSeries() const {
        const auto acidRi = acidicContactResidueIndex();
        if (!acidRi) return std::vector<double>(cadence_.stepCount(), kNaN);
        return saltBridgeDistanceSeriesFor(*acidRi);
    }

    bool saltBridgeMediationApplies() const {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return false;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = body_.run.protein->atom(atom_);
        return a.element == model::Element::C
               && validResidue(p, a.residueIndex)
               && p.residue(static_cast<std::size_t>(a.residueIndex)).aminoAcid == model::AminoAcid::ARG
               && p.atomLabel(atom_, model::NamingConvention::Iupac) == QStringLiteral("CZ")
               && acidicContactResidueIndex().has_value();
    }

    QJsonArray mediationChainsJson(const std::vector<EvaluatedResponse>& characterized,
                                   const std::vector<std::size_t>& rows) const {
        QJsonArray out;
        if (!saltBridgeMediationApplies()) return out;

        const EvaluatedResponse* chosen = nullptr;
        const EvaluatedTarget* target = nullptr;
        const std::array<CatalogueId, 4> preference = {{
            CatalogueId::EFG_MOLCOMP,
            CatalogueId::EFG_ISO,
            CatalogueId::FIELD_MOLCOMP,
            CatalogueId::FIELD_ISO,
        }};
        for (CatalogueId id : preference) {
            for (const EvaluatedResponse& e : characterized) {
                if (e.id != id) continue;
                const EvaluatedTarget* primary = primaryTarget(e);
                if (!primary) continue;
                chosen = &e;
                target = primary;
                break;
            }
            if (chosen && target) break;
        }
        if (!chosen || !target) return out;

        std::vector<double> sigma;
        std::vector<double> distance;
        std::vector<double> mediator;
        const std::vector<double> distanceSeries = saltBridgeDistanceSeries();
        for (std::size_t row : rows) {
            if (row >= target->sigma_series.size()
                || row >= target->driver_series.size()
                || row >= distanceSeries.size())
                continue;
            const double sig = target->sigma_series[row];
            const double med = target->driver_series[row];
            const double dist = distanceSeries[row];
            if (!finite(sig) || !finite(med) || !finite(dist)) continue;
            sigma.push_back(sig);
            mediator.push_back(med);
            distance.push_back(dist);
        }
        if (sigma.size() < kNullBlock * 2) return out;

        QJsonObject chain;
        chain.insert(QStringLiteral("driver"), QStringLiteral("salt_bridge.distance_A"));
        chain.insert(QStringLiteral("driver_geometry"),
                     QStringLiteral("ARG guanidinium N(NH1/NH2/NE) <-> acidic carboxylate O(OD/OE)"));
        chain.insert(QStringLiteral("mediator"), target->driver_channel);
        chain.insert(QStringLiteral("catalogue_id"),
                     QString::fromLatin1(catalogueRow(chosen->id).id_name));
        chain.insert(QStringLiteral("sigma_target"), target->sigma_target);
        chain.insert(QStringLiteral("n"), static_cast<int>(sigma.size()));

        QJsonObject links;
        links.insert(QStringLiteral("dist_to_field_r"), jd(pearsonR(distance, mediator)));
        links.insert(QStringLiteral("field_to_sigma_r"), jd(pearsonR(mediator, sigma)));
        links.insert(QStringLiteral("dist_to_sigma_r"), jd(pearsonR(distance, sigma)));
        chain.insert(QStringLiteral("links"), links);
        chain.insert(QStringLiteral("partial_r"), jd(partialR(sigma, distance, mediator)));
        const auto [lo, hi] = blockBootstrapPartialR(sigma, distance, mediator);
        chain.insert(QStringLiteral("bootstrap_lo"), jd(lo));
        chain.insert(QStringLiteral("bootstrap_hi"), jd(hi));
        out.append(chain);
        return out;
    }

    std::vector<Candidate> wibergTensorCandidates(const QString& channel,
                                                  const std::vector<double>& driver) const {
        std::vector<Candidate> cands;
        cands.push_back({driver, sigmaIsoSeries(), QStringLiteral("iso"), channel});
        cands.push_back({driver, sigmaInvariantSeries(QStringLiteral("aniso")),
                         QStringLiteral("invariant:aniso"), channel});
        cands.push_back({driver, sigmaInvariantSeries(QStringLiteral("span")),
                         QStringLiteral("invariant:span"), channel});
        for (int sc = 0; sc < 6; ++sc) {
            const auto smc = static_cast<MolComp>(sc);
            cands.push_back({driver, sigmaMolCompSeries(sigma_mol_total_, smc),
                             QStringLiteral("molcomp:%1").arg(QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(sc)])),
                             channel});
        }
        return cands;
    }

    std::optional<EvaluatedResponse> evaluateWibergMember(const QString& channel,
                                                          const std::vector<double>& driver,
                                                          const std::vector<std::size_t>& rows) const {
        return evaluateCatalogueRowFromCandidates(CatalogueId::WIBERG,
                                                  wibergTensorCandidates(channel, driver),
                                                  rows);
    }

    QJsonObject scalarSummaryJson(const std::vector<double>& values) const {
        const DistributionSummary distribution = SummarizeDistribution(values);
        QJsonObject o = DistributionSummaryJson(distribution);
        std::vector<double> finiteValues;
        finiteValues.reserve(values.size());
        for (double v : values)
            if (finite(v)) finiteValues.push_back(v);
        if (finiteValues.empty()) {
            o.insert(QStringLiteral("skewness"), QJsonValue(QJsonValue::Null));
            o.insert(QStringLiteral("excess_kurtosis"), QJsonValue(QJsonValue::Null));
            o.insert(QStringLiteral("bimodality"), QJsonValue(QJsonValue::Null));
            o.insert(QStringLiteral("constant"), false);
            return o;
        }
        const auto [mn, mx] = std::minmax_element(finiteValues.begin(), finiteValues.end());
        // A (near-)constant series has no shape, and the flatness is itself the signal
        // (a static partner / locked geometry). boost special-cases EXACTLY zero variance
        // ("a constant dataset has no skewness" -> 0); near-constant values (sd at the
        // floating-point-noise floor relative to magnitude) would otherwise yield noise
        // skewness/kurtosis -- so detect them here, mark `constant`, and null the shape.
        constexpr double kRelSpreadFloor = 1e-9;  // numerical-validity floor, NOT a data cutoff
        const double scale = std::max({std::abs(*mn), std::abs(*mx), 1e-300});
        const double sd = distribution.sd;
        const bool hasSpread = finiteValues.size() >= 2 && sd > kRelSpreadFloor * scale;
        o.insert(QStringLiteral("constant"), !hasSpread);
        // Shape-signal for descriptive shape detection (priors to confirm, never imposed):
        // skewness, excess kurtosis, and Sarle's bimodality coefficient
        // BC = (skew^2 + 1) / kurtosis (raw 4th moment = excess + 3).
        QJsonValue skewVal(QJsonValue::Null), exKurtVal(QJsonValue::Null), bcVal(QJsonValue::Null);
        if (hasSpread && finiteValues.size() >= 3) {
            try {
                const double skew = bmstats::skewness(finiteValues);
                const double kurt = bmstats::kurtosis(finiteValues);
                skewVal = jd(skew);
                exKurtVal = jd(kurt - 3.0);
                if (kurt > 0.0) bcVal = jd((skew * skew + 1.0) / kurt);
            } catch (...) {}
        }
        o.insert(QStringLiteral("skewness"), skewVal);
        o.insert(QStringLiteral("excess_kurtosis"), exKurtVal);
        o.insert(QStringLiteral("bimodality"), bcVal);
        return o;
    }

    QJsonObject wibergResponseOrNull(const std::optional<EvaluatedResponse>& e) const {
        if (e) return responseJson(*e);
        QJsonObject o;
        o.insert(QStringLiteral("catalogue_id"), QStringLiteral("WIBERG"));
        QJsonArray considerations;
        considerations.append(QStringLiteral("ZERO_VARIANCE_OR_INSUFFICIENT_SUPPORT"));
        o.insert(QStringLiteral("considerations"), considerations);
        o.insert(QStringLiteral("targets"), QJsonArray{});
        return o;
    }

    std::vector<double> wibergSeriesForBond(std::uint64_t key) const {
        const auto it = bond_series_.find(key);
        if (it != bond_series_.end()) return it->second.values;
        const auto p = pairFromKey(key);
        return context_.mopacWibergSeries(p.first, p.second);
    }

    double finiteMeanOrNan(const std::vector<double>& values) const {
        double sum = 0.0;
        int n = 0;
        for (double v : values) {
            if (!finite(v)) continue;
            sum += v;
            ++n;
        }
        return n > 0 ? sum / static_cast<double>(n) : kNaN;
    }

    QJsonObject buildWibergJson() const {
        QJsonObject w;
        w.insert(QStringLiteral("catalogue_id"), QStringLiteral("WIBERG"));
        w.insert(QStringLiteral("slot_kind"), QStringLiteral("STRUCTURED"));
        const std::vector<std::size_t> rows = cadence_.sigmaRows();
        if (!hasMappedBond()) {
            w.insert(QStringLiteral("present"), false);
            return w;
        }
        w.insert(QStringLiteral("present"), true);

        const auto valencyScore =
            evaluateWibergMember(QStringLiteral("mopac.valency"), mopac_valency_.values, rows);
        QJsonObject valency;
        valency.insert(QStringLiteral("order_series_summary"), scalarSummaryJson(mopac_valency_.values));
        valency.insert(QStringLiteral("response"), wibergResponseOrNull(valencyScore));
        w.insert(QStringLiteral("valency"), valency);

        std::uint64_t dominantKey = 0;
        double dominantMean = -std::numeric_limits<double>::infinity();
        for (const auto& item : bond_series_) {
            if (static_bonds_.count(item.first) == 0) continue;
            const double mean = finiteMeanOrNan(item.second.values);
            if (finite(mean) && mean > dominantMean) {
                dominantMean = mean;
                dominantKey = item.first;
            }
        }

        QJsonArray own;
        std::vector<std::pair<QString, std::vector<double>>> ownSeries;
        for (const auto& item : bond_series_) {
            const auto sit = static_bonds_.find(item.first);
            if (sit == static_bonds_.end()) continue;
            const auto p = pairFromKey(item.first);
            const std::size_t partner = p.first == atom_ ? p.second : p.first;
            const QString channel = QStringLiteral("mopac_wiberg_order|atom%1").arg(partner);
            const auto score = evaluateWibergMember(channel, item.second.values, rows);
            QJsonObject o;
            o.insert(QStringLiteral("partner_atom_index"), static_cast<int>(partner));
            if (partner < body_.run.protein->atomCount())
                o.insert(QStringLiteral("partner_element"),
                         elementName(body_.run.protein->atom(partner).element));
            o.insert(QStringLiteral("bond_type"), bondOrderName(sit->second.order));
            o.insert(QStringLiteral("bond_category"), bondCategoryName(sit->second.category));
            o.insert(QStringLiteral("order_summary"), scalarSummaryJson(item.second.values));
            o.insert(QStringLiteral("is_dominant"), item.first == dominantKey);
            o.insert(QStringLiteral("response"), wibergResponseOrNull(score));
            own.append(o);
            ownSeries.push_back({channel, item.second.values});
        }
        w.insert(QStringLiteral("own_bonds"), own);

        QJsonArray neighborhood;
        if (body_.run.protein) {
            const model::QtTopology& topo = body_.run.protein->topology();
            std::set<std::uint64_t> seen;
            for (const auto& item : static_bonds_) {
                const std::size_t other = item.second.other;
                if (other >= body_.run.protein->atomCount()) continue;
                for (int32_t bi : topo.bondIndicesForAtom(other)) {
                    if (bi < 0 || static_cast<std::size_t>(bi) >= topo.bondCount()) continue;
                    const model::QtBond& b = topo.bondAt(static_cast<std::size_t>(bi));
                    if (b.atomIndexA < 0 || b.atomIndexB < 0) continue;
                    const std::size_t a = static_cast<std::size_t>(b.atomIndexA);
                    const std::size_t bb = static_cast<std::size_t>(b.atomIndexB);
                    if (a == atom_ || bb == atom_) continue;
                    const std::uint64_t key = pairKey(a, bb);
                    if (!seen.insert(key).second) continue;
                    const std::vector<double> series = wibergSeriesForBond(key);
                    const QString channel =
                        QStringLiteral("mopac_wiberg_neighbor|atom%1-atom%2").arg(a).arg(bb);
                    const auto score = evaluateWibergMember(channel, series, rows);
                    QJsonObject o;
                    o.insert(QStringLiteral("resonance_partner_bond"),
                             QStringLiteral("atom%1-atom%2").arg(a).arg(bb));
                    o.insert(QStringLiteral("reach"), QStringLiteral("hops<=1"));
                    o.insert(QStringLiteral("order_summary"), scalarSummaryJson(series));
                    o.insert(QStringLiteral("response"), wibergResponseOrNull(score));
                    neighborhood.append(o);
                }
            }
        }
        w.insert(QStringLiteral("neighborhood"), neighborhood);

        QJsonObject collinearity;
        double bestAbs = -1.0;
        QString bestPair;
        double bestR = kNaN;
        for (std::size_t i = 0; i < ownSeries.size(); ++i) {
            for (std::size_t j = i + 1; j < ownSeries.size(); ++j) {
                std::vector<double> x;
                std::vector<double> y;
                pairedOnSigmaRows(ownSeries[i].second, ownSeries[j].second, rows, x, y);
                const double r = pearsonR(x, y);
                if (finite(r) && std::abs(r) > bestAbs) {
                    bestAbs = std::abs(r);
                    bestR = r;
                    bestPair = QStringLiteral("%1 <-> %2").arg(ownSeries[i].first, ownSeries[j].first);
                }
            }
        }
        collinearity.insert(QStringLiteral("best_pair"), bestPair);
        collinearity.insert(QStringLiteral("r"), jd(bestR));
        collinearity.insert(QStringLiteral("reported_not_precut"), true);
        w.insert(QStringLiteral("collinearity"), collinearity);
        return w;
    }

    QJsonObject fieldReversalTailDiagnosticsJson() const {
        QJsonObject o;
        const std::vector<std::size_t> rows = cadence_.sigmaRows();
        const std::vector<double> sigma = sigmaIsoSeries();
        const std::vector<double> signedField = field_E_z_mopac_.values;
        const std::vector<double> e2 = field_abs_E2_mopac_.values;
        std::vector<double> xSigned;
        std::vector<double> ySigma;
        pairedOnSigmaRows(sigma, signedField, rows, xSigned, ySigma);
        std::vector<double> xE2;
        std::vector<double> yE2;
        pairedOnSigmaRows(sigma, e2, rows, xE2, yE2);

        int pos = 0;
        int neg = 0;
        int zero = 0;
        for (double v : xSigned) {
            if (v > 0.0) ++pos;
            else if (v < 0.0) ++neg;
            else ++zero;
        }
        o.insert(QStringLiteral("signed_field_positive_n"), pos);
        o.insert(QStringLiteral("signed_field_negative_n"), neg);
        o.insert(QStringLiteral("signed_field_zero_n"), zero);
        o.insert(QStringLiteral("field_sign_balance"),
                 (pos + neg) > 0
                     ? jd(static_cast<double>(std::min(pos, neg)) / static_cast<double>(pos + neg))
                     : jd(kNaN));

        const DistributionSummary e2Summary = SummarizeDistribution(xE2);
        o.insert(QStringLiteral("E2_summary"), DistributionSummaryJson(e2Summary));
        int highTail = 0;
        for (double v : xE2)
            if (e2Summary.hasFinite() && finite(v) && v >= e2Summary.p95) ++highTail;
        o.insert(QStringLiteral("high_field_tail_threshold_E2"), e2Summary.hasFinite() ? jd(e2Summary.p95) : jd(kNaN));
        o.insert(QStringLiteral("high_field_tail_n"), highTail);

        const OlsResult linear = ols(xSigned, ySigma);
        std::vector<double> residual;
        residual.reserve(std::min(xSigned.size(), ySigma.size()));
        for (std::size_t i = 0; i < xSigned.size() && i < ySigma.size(); ++i) {
            if (!finite(linear.slope) || !finite(linear.intercept)) continue;
            residual.push_back(ySigma[i] - (linear.intercept + linear.slope * xSigned[i]));
        }
        std::vector<double> e2ForResidual;
        e2ForResidual.reserve(rows.size());
        std::vector<double> residualAligned;
        residualAligned.reserve(rows.size());
        for (std::size_t row : rows) {
            if (row >= sigma.size() || row >= signedField.size() || row >= e2.size()) continue;
            if (!finite(sigma[row]) || !finite(signedField[row]) || !finite(e2[row])) continue;
            if (!finite(linear.slope) || !finite(linear.intercept)) continue;
            e2ForResidual.push_back(e2[row]);
            residualAligned.push_back(sigma[row] - (linear.intercept + linear.slope * signedField[row]));
        }
        const OlsResult curvature = ols(e2ForResidual, residualAligned);
        o.insert(QStringLiteral("linear_signed_field_slope"), jd(linear.slope));
        o.insert(QStringLiteral("linear_signed_field_r"), jd(pearsonR(xSigned, ySigma)));
        o.insert(QStringLiteral("residual_curvature_slope_E2"), jd(curvature.slope));
        o.insert(QStringLiteral("residual_curvature_r_E2"), jd(pearsonR(e2ForResidual, residualAligned)));
        o.insert(QStringLiteral("quadratic_raw_slope_E2"), jd(ols(xE2, yE2).slope));
        o.insert(QStringLiteral("quadratic_raw_r_E2"), jd(pearsonR(xE2, yE2)));
        o.insert(QStringLiteral("reversal_testable"), pos > 0 && neg > 0);
        o.insert(QStringLiteral("structural_reason"),
                 (pos > 0 && neg > 0) ? QString() : QStringLiteral("not testable for reversal"));
        o.insert(QStringLiteral("emit_policy"), QStringLiteral("bounded_summary_no_raw_field_sidecar"));
        return o;
    }

    QJsonObject targetJson(CatalogueId id, const EvaluatedTarget& t) const {
        QJsonObject o;
        o.insert(QStringLiteral("sigma_target"), t.sigma_target);
        o.insert(QStringLiteral("frame"), frameTagName(t.frame));
        o.insert(QStringLiteral("driver"), t.driver_channel);
        QJsonObject law;
        law.insert(QStringLiteral("slope"), jd(t.slope));
        law.insert(QStringLiteral("intercept"), jd(t.intercept));
        law.insert(QStringLiteral("r"), jd(t.r));
        law.insert(QStringLiteral("null_z"), jd(t.null_z));
        law.insert(QStringLiteral("xi"), jd(t.xi));
        law.insert(QStringLiteral("pchip_mid"), jd(t.pchip_mid));
        law.insert(QStringLiteral("leverage"), jd(t.leverage));
        law.insert(QStringLiteral("coverage"), t.coverage);
        o.insert(QStringLiteral("response_law"), law);
        QJsonObject g;
        g.insert(QStringLiteral("delta_survives"), t.delta_survives);
        g.insert(QStringLiteral("seg_min_r"), jd(t.seg_min_r));
        g.insert(QStringLiteral("circshift_p"), jd(t.circshift_p));
        o.insert(QStringLiteral("gauntlet"), g);
        if (id == CatalogueId::DIHEDRAL && finite(t.eta_squared)) {
            QJsonObject change;
            change.insert(QStringLiteral("eta_sq"), jd(t.eta_squared));
            change.insert(QStringLiteral("change_point_ratio"), QJsonValue(QJsonValue::Null));
            change.insert(QStringLiteral("change_point_p"), QJsonValue(QJsonValue::Null));
            o.insert(QStringLiteral("change_metrics"), change);
        } else {
            o.insert(QStringLiteral("change_metrics"), QJsonValue(QJsonValue::Null));
        }
        return o;
    }

    QJsonObject evaluatedTargetJson(CatalogueId id,
                                    const Candidate& c,
                                    const std::vector<std::size_t>& rows) const {
        std::vector<double> x;
        std::vector<double> y;
        pairedOnSigmaRows(c.sigma, c.driver, rows, x, y);
        if (x.size() < 3) {
            QJsonObject miss = targetJson(id, nullTargetFor(c, id));
            miss.insert(QStringLiteral("bounded_reveal_response"), true);
            return miss;
        }

        EvaluatedTarget t;
        t.sigma_target = c.sigma_target;
        t.driver_channel = c.driver_channel;
        t.frame = frameForTarget(id, c.sigma_target);
        const OlsResult fit = ols(x, y);
        t.slope = fit.slope;
        t.intercept = fit.intercept;
        t.r = pearsonR(x, y);
        t.leverage = leverageTop1(x, y);
        t.coverage = coverageJson(x);
        t.delta_survives = deltaSurvives(x, y, t.r);
        t.seg_min_r = segMinR(x, y);
        const LeadLagResult ll = leadLag(x, y);
        t.best_lag = ll.best_lag;
        t.lead_r = ll.lead_r;

        QJsonObject o = targetJson(id, t);
        o.insert(QStringLiteral("bounded_reveal_response"), true);
        return o;
    }

    int presentCount(const Mat3Series& s) const {
        return static_cast<int>(std::count(s.present.begin(), s.present.end(), true));
    }

    QJsonObject tensorComponentSummariesJson(const Mat3Series& s) const {
        QJsonObject o;
        for (int c = 0; c < 6; ++c) {
            o.insert(QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(c)]),
                     scalarSummaryJson(molCompFrom(s, static_cast<MolComp>(c))));
        }
        return o;
    }

    QJsonArray componentResponsesJson(const Mat3Series& mechanism,
                                      const QString& channel,
                                      const std::vector<std::size_t>& rows,
                                      CatalogueId id,
                                      bool crossProduct) const {
        QJsonArray out;
        for (int dc = 0; dc < 6; ++dc) {
            const auto dcomp = static_cast<MolComp>(dc);
            const std::vector<double> driver = molCompFrom(mechanism, dcomp);
            const QString dcName = QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(dc)]);
            for (int sc = 0; sc < 6; ++sc) {
                if (!crossProduct && sc != dc) continue;
                const auto scomp = static_cast<MolComp>(sc);
                const QString scName = QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(sc)]);
                const std::array<std::pair<const Mat3Series*, const char*>, 2> sigmaSources = {{
                    {&sigma_mol_total_, "sigma_total"},
                    {&sigma_mol_para_, "sigma_para"},
                }};
                for (const auto& sigmaSource : sigmaSources) {
                    const QString basis = QString::fromLatin1(sigmaSource.second);
                    Candidate c{driver,
                                sigmaMolCompSeries(*sigmaSource.first, scomp),
                                QStringLiteral("%1.molcomp:%2").arg(basis, scName),
                                QStringLiteral("%1|molcomp:%2").arg(channel, dcName)};
                    QJsonObject o = evaluatedTargetJson(id, c, rows);
                    o.insert(QStringLiteral("sigma_basis"), basis);
                    o.insert(QStringLiteral("gauge_flag"), basis == QStringLiteral("sigma_para"));
                    o.insert(QStringLiteral("driver_component"), dcName);
                    o.insert(QStringLiteral("sigma_component"), scName);
                    o.insert(QStringLiteral("component_pair"),
                             QStringLiteral("%1->%2").arg(dcName, scName));
                    out.append(o);
                }
            }
        }
        return out;
    }

    std::vector<double> sigmaIsoForBasis(const QString& basis) const {
        if (basis == QStringLiteral("sigma_para")) return componentSeriesT0(sigma_para_);
        return sigmaIsoSeries();
    }

    QJsonArray contractionResponsesJson(const Mat3Series& mechanism,
                                        const QString& channel,
                                        const std::vector<std::size_t>& rows,
                                        CatalogueId id) const {
        QJsonArray out;
        const std::array<std::pair<const Mat3Series*, const char*>, 2> sigmaSources = {{
            {&sigma_mol_total_, "sigma_total"},
            {&sigma_mol_para_, "sigma_para"},
        }};
        for (const auto& sigmaSource : sigmaSources) {
            const QString basis = QString::fromLatin1(sigmaSource.second);
            for (bool centered : {false, true}) {
                for (bool cosine : {false, true}) {
                    const QString form = cosine ? QStringLiteral("cosine_alignment")
                                                : QStringLiteral("raw_inner");
                    const std::vector<double> driver =
                        tensorContractionSeries(*sigmaSource.first, mechanism, rows, centered, cosine);
                    std::vector<double> sigma = sigmaIsoForBasis(basis);
                    if (centered) sigma = centeredScalarSeries(sigma, rows);
                    Candidate c{driver,
                                sigma,
                                centered ? QStringLiteral("%1.delta_iso").arg(basis)
                                         : QStringLiteral("%1.iso").arg(basis),
                                QStringLiteral("%1|contraction:%2:%3")
                                    .arg(channel, form,
                                         centered ? QStringLiteral("delta_sigma")
                                                  : QStringLiteral("sigma"))};
                    QJsonObject o;
                    o.insert(QStringLiteral("sigma_basis"), basis);
                    o.insert(QStringLiteral("gauge_flag"), basis == QStringLiteral("sigma_para"));
                    o.insert(QStringLiteral("form"), form);
                    o.insert(QStringLiteral("sigma_delta"), centered);
                    o.insert(QStringLiteral("series_summary"), scalarSummaryJson(driver));
                    o.insert(QStringLiteral("response"), evaluatedTargetJson(id, c, rows));
                    out.append(o);
                }
            }
        }
        return out;
    }

    QJsonObject bsHmDivergenceJson(int type, const std::vector<std::size_t>& rows) const {
        const Mat3Series& bs = bs_per_type_mol_[static_cast<std::size_t>(type)];
        const Mat3Series& hm = hm_per_type_mol_[static_cast<std::size_t>(type)];
        const std::vector<double> ringContribution = ringContributionMagnitudeSeries(type);
        QJsonObject o;
        o.insert(QStringLiteral("level"), QStringLiteral("per_type"));
        o.insert(QStringLiteral("family"), QStringLiteral("bs_hm_divergence"));
        o.insert(QStringLiteral("source_key"), QStringLiteral("ring.type%1").arg(type));
        o.insert(QStringLiteral("frame"), QStringLiteral("molecular"));
        o.insert(QStringLiteral("projection"), QStringLiteral("T_mol = R^T * T_lab * R"));

        std::vector<double> divergence(bs.values.size(), kNaN);
        std::vector<double> divergenceRatio(bs.values.size(), kNaN);
        for (std::size_t row : rows) {
            if (row >= bs.values.size() || row >= hm.values.size()) continue;
            if (!bs.present[row] || !hm.present[row]) continue;
            divergence[row] = tensorFrobenius(bs.values[row] - hm.values[row]);
            const double scale =
                row < ringContribution.size() ? ringContribution[row] : kNaN;
            if (finite(scale) && scale > 1.0e-15)
                divergenceRatio[row] = divergence[row] / scale;
        }
        o.insert(QStringLiteral("bs_hm_divergence"), scalarSummaryJson(divergence));
        o.insert(QStringLiteral("bs_hm_divergence_to_ring_contribution_ratio"),
                 scalarSummaryJson(divergenceRatio));

        QJsonObject componentDeltas;
        for (int c = 0; c < 6; ++c) {
            const auto comp = static_cast<MolComp>(c);
            const std::vector<double> bsc = molCompFrom(bs, comp);
            const std::vector<double> hmc = molCompFrom(hm, comp);
            std::vector<double> delta(std::max(bsc.size(), hmc.size()), kNaN);
            for (std::size_t i = 0; i < delta.size(); ++i) {
                const double bv = i < bsc.size() ? bsc[i] : kNaN;
                const double hv = i < hmc.size() ? hmc[i] : kNaN;
                if (finite(bv) && finite(hv)) delta[i] = bv - hv;
            }
            componentDeltas.insert(QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(c)]),
                                   scalarSummaryJson(delta));
        }
        o.insert(QStringLiteral("component_deltas"), componentDeltas);
        o.insert(QStringLiteral("subspace_compare"),
                 subspaceCompareObject(QStringLiteral("bs_vs_hm_tensor_components_type%1").arg(type),
                                       bsTensorFamilyForType(type),
                                       hmTensorFamilyForType(type),
                                       rows));

        QJsonArray contractions;
        const std::array<std::pair<const Mat3Series*, const char*>, 2> sigmaSources = {{
            {&sigma_mol_total_, "sigma_total"},
            {&sigma_mol_para_, "sigma_para"},
        }};
        for (const auto& sigmaSource : sigmaSources) {
            const QString basis = QString::fromLatin1(sigmaSource.second);
            for (bool centered : {false, true}) {
                for (bool cosine : {false, true}) {
                    const std::vector<double> bsDriver =
                        tensorContractionSeries(*sigmaSource.first, bs, rows, centered, cosine);
                    const std::vector<double> hmDriver =
                        tensorContractionSeries(*sigmaSource.first, hm, rows, centered, cosine);
                    std::vector<double> x;
                    std::vector<double> y;
                    pairedOnSigmaRows(bsDriver, hmDriver, rows, x, y);
                    const double r = pearsonR(x, y);
                    QJsonObject c;
                    c.insert(QStringLiteral("sigma_basis"), basis);
                    c.insert(QStringLiteral("gauge_flag"), basis == QStringLiteral("sigma_para"));
                    c.insert(QStringLiteral("form"),
                             cosine ? QStringLiteral("cosine_alignment")
                                    : QStringLiteral("raw_inner"));
                    c.insert(QStringLiteral("sigma_delta"), centered);
                    c.insert(QStringLiteral("finite_n"), static_cast<int>(x.size()));
                    c.insert(QStringLiteral("r_bs_hm_contraction"), jd(r));
                    contractions.append(c);
                }
            }
        }
        o.insert(QStringLiteral("contraction_pearson_correlations"), contractions);
        return o;
    }

    std::vector<double> ringContributionMagnitudeSeries(int type) const {
        std::vector<double> out(cadence_.stepCount(), kNaN);
        const int ringOrd = mechanismOrd(QStringLiteral("ring_jb"));
        for (const auto& item : relationships_) {
            const RelationshipKey& key = item.first;
            if (key.mechanism_ord != ringOrd || key.source_category_ord != type) continue;
            const std::vector<double> contribution =
                item.second.contribution.dense(cadence_.stepCount());
            for (std::size_t i = 0; i < out.size() && i < contribution.size(); ++i) {
                if (!finite(contribution[i])) continue;
                const double v = std::abs(contribution[i]);
                out[i] = finite(out[i]) ? out[i] + v : v;
            }
        }
        return out;
    }

    static int signClass(double v) {
        if (!finite(v) || v == 0.0) return 0;
        return v > 0.0 ? 1 : -1;
    }

    int crossingCount(const std::vector<double>& values) const {
        int crossings = 0;
        int prev = 0;
        for (double v : values) {
            const int s = signClass(v);
            if (s == 0) continue;
            if (prev != 0 && s != prev) ++crossings;
            prev = s;
        }
        return crossings;
    }

    QJsonObject magicAngleSummaryJson(const std::vector<double>& cosTheta,
                                      const std::vector<double>& dipolar,
                                      int crossingOverride = -1) const {
        int pos = 0;
        int neg = 0;
        int zero = 0;
        for (double v : dipolar) {
            const int s = signClass(v);
            if (s > 0) ++pos;
            else if (s < 0) ++neg;
            else if (finite(v)) ++zero;
        }
        const int n = pos + neg + zero;
        QJsonObject o;
        o.insert(QStringLiteral("cos_theta"), scalarSummaryJson(cosTheta));
        o.insert(QStringLiteral("dipolar"), scalarSummaryJson(dipolar));
        o.insert(QStringLiteral("magic_angle_abs_cos"), jd(1.0 / std::sqrt(3.0)));
        o.insert(QStringLiteral("crossing_count"),
                 crossingOverride >= 0 ? crossingOverride : crossingCount(dipolar));
        o.insert(QStringLiteral("above_magic_fraction"), n > 0 ? jd(static_cast<double>(pos) / n) : jd(kNaN));
        o.insert(QStringLiteral("below_magic_fraction"), n > 0 ? jd(static_cast<double>(neg) / n) : jd(kNaN));
        o.insert(QStringLiteral("on_crossing_fraction"), n > 0 ? jd(static_cast<double>(zero) / n) : jd(kNaN));
        return o;
    }

    QJsonObject angularContextForMechanism(int mechanismOrdValue,
                                           std::optional<int> sourceCategory = std::nullopt) const {
        std::vector<double> cosTheta;
        std::vector<double> dipolar;
        int relationships = 0;
        int crossings = 0;
        const std::size_t n = cadence_.stepCount();
        for (const auto& item : relationships_) {
            const RelationshipKey& key = item.first;
            if (key.mechanism_ord != mechanismOrdValue) continue;
            if (sourceCategory && key.source_category_ord != *sourceCategory) continue;
            const RelationshipSeries& rel = item.second;
            ++relationships;
            const std::vector<double> relCos = rel.cos_theta.dense(n);
            const std::vector<double> relDip = rel.dipolar.dense(n);
            crossings += crossingCount(relDip);
            for (double v : relCos)
                if (finite(v)) cosTheta.push_back(v);
            for (double v : relDip)
                if (finite(v)) dipolar.push_back(v);
        }
        QJsonObject o = magicAngleSummaryJson(cosTheta, dipolar, crossings);
        o.insert(QStringLiteral("relationship_count"), relationships);
        o.insert(QStringLiteral("mechanism"), mechanismName(mechanismOrdValue));
        if (sourceCategory) o.insert(QStringLiteral("source_category_ord"), *sourceCategory);
        return o;
    }

    QJsonArray sigmaFlipFollowJson(const RelationshipSeries& rel) const {
        QJsonArray out;
        const std::size_t n = cadence_.stepCount();
        const std::vector<double> dip = rel.dipolar.dense(n);
        const std::array<std::pair<const Mat3Series*, const char*>, 2> sigmaSources = {{
            {&sigma_mol_total_, "sigma_total"},
            {&sigma_mol_para_, "sigma_para"},
        }};
        for (const auto& sigmaSource : sigmaSources) {
            const QString basis = QString::fromLatin1(sigmaSource.second);
            for (int sc = 0; sc < 6; ++sc) {
                const auto comp = static_cast<MolComp>(sc);
                const QString compName = QString::fromLatin1(kMolCompNames[static_cast<std::size_t>(sc)]);
                const std::vector<double> sigma = sigmaMolCompSeries(*sigmaSource.first, comp);
                std::vector<int> prevFinite(n, -1);
                std::vector<int> nextFinite(n, -1);
                int last = -1;
                for (std::size_t i = 0; i < n; ++i) {
                    if (i < sigma.size() && finite(sigma[i])) last = static_cast<int>(i);
                    prevFinite[i] = last;
                }
                last = -1;
                for (std::size_t i = n; i-- > 0;) {
                    if (i < sigma.size() && finite(sigma[i])) last = static_cast<int>(i);
                    nextFinite[i] = last;
                }
                int crossings = 0;
                int compared = 0;
                int followed = 0;
                for (std::size_t i = 1; i < n; ++i) {
                    const int s0 = signClass(dip[i - 1]);
                    const int s1 = signClass(dip[i]);
                    if (s0 == 0 || s1 == 0 || s0 == s1) continue;
                    ++crossings;
                    const int before = prevFinite[i - 1];
                    const int after = nextFinite[i];
                    if (before < 0 || after < 0 || before == after) continue;
                    const int ds = signClass(sigma[static_cast<std::size_t>(after)]
                                             - sigma[static_cast<std::size_t>(before)]);
                    const int dd = signClass(dip[i] - dip[i - 1]);
                    if (ds == 0 || dd == 0) continue;
                    ++compared;
                    if (ds == dd) ++followed;
                }
                QJsonObject o;
                o.insert(QStringLiteral("sigma_basis"), basis);
                o.insert(QStringLiteral("gauge_flag"), basis == QStringLiteral("sigma_para"));
                o.insert(QStringLiteral("sigma_component"), compName);
                o.insert(QStringLiteral("crossing_count"), crossings);
                o.insert(QStringLiteral("compared_crossings"), compared);
                o.insert(QStringLiteral("followed_crossings"), followed);
                o.insert(QStringLiteral("opposed_crossings"), compared - followed);
                o.insert(QStringLiteral("follow_fraction"),
                         compared > 0 ? jd(static_cast<double>(followed) / compared) : jd(kNaN));
                o.insert(QStringLiteral("sigma_flips_on_driver_crossing"),
                         compared > 0 ? QJsonValue(followed > (compared - followed))
                                      : QJsonValue(QJsonValue::Null));
                out.append(o);
            }
        }
        return out;
    }

    QJsonObject typeTensorRevealJson(const QString& family,
                                     const QString& sourceKey,
                                     const Mat3Series& mechanism,
                                     const std::vector<std::size_t>& rows,
                                     CatalogueId id,
                                     const QJsonValue& angularContext = QJsonValue(QJsonValue::Null)) const {
        QJsonObject o;
        o.insert(QStringLiteral("level"), QStringLiteral("per_type"));
        o.insert(QStringLiteral("family"), family);
        o.insert(QStringLiteral("source_key"), sourceKey);
        o.insert(QStringLiteral("frame"), QStringLiteral("molecular"));
        o.insert(QStringLiteral("projection"), QStringLiteral("T_mol = R^T * T_lab * R"));
        o.insert(QStringLiteral("present_frames"), presentCount(mechanism));
        o.insert(QStringLiteral("component_summaries"), tensorComponentSummariesJson(mechanism));
        o.insert(QStringLiteral("component_responses"),
                 componentResponsesJson(mechanism, sourceKey, rows, id, true));
        o.insert(QStringLiteral("contractions"),
                 contractionResponsesJson(mechanism, sourceKey, rows, id));
        if (!angularContext.isNull()) o.insert(QStringLiteral("angular_context"), angularContext);
        return o;
    }

    QJsonArray perTypeTensorRevealsJson(const std::vector<std::size_t>& rows) const {
        QJsonArray out;
        const QJsonObject mcAngular =
            angularContextForMechanism(mechanismOrd(QStringLiteral("mc_lit_valid")));
        for (std::size_t i = 0; i < mc_tensor_mol_series_.size(); ++i) {
            out.append(typeTensorRevealJson(QStringLiteral("mc"),
                                            QString::fromLatin1(kMcTensorFields[i].key),
                                            mc_tensor_mol_series_[i], rows,
                                            CatalogueId::MCCONNELL, mcAngular));
        }
        const int ringOrd = mechanismOrd(QStringLiteral("ring_jb"));
        for (int type = 0; type < 8; ++type) {
            const QJsonObject ringAngular = angularContextForMechanism(ringOrd, type);
            out.append(typeTensorRevealJson(QStringLiteral("bs"),
                                            QStringLiteral("bs.type%1").arg(type),
                                            bs_per_type_mol_[static_cast<std::size_t>(type)],
                                            rows, CatalogueId::RING_H, ringAngular));
            out.append(typeTensorRevealJson(QStringLiteral("hm"),
                                            QStringLiteral("hm.type%1").arg(type),
                                            hm_per_type_mol_[static_cast<std::size_t>(type)],
                                            rows, CatalogueId::RING_H, ringAngular));
            out.append(bsHmDivergenceJson(type, rows));
        }
        out.append(typeTensorRevealJson(QStringLiteral("tripeptide"),
                                        QStringLiteral("tripeptide_bb_shielding"),
                                        tripeptide_bb_shielding_mol_, rows,
                                        CatalogueId::MCCONNELL));
        out.append(typeTensorRevealJson(QStringLiteral("larsen"),
                                        QStringLiteral("larsen_hbond_shielding"),
                                        larsen_hbond_shielding_mol_, rows,
                                        CatalogueId::MCCONNELL));
        return out;
    }

    bool isAngularPartner(const RelationshipKey& key) const {
        return key.mechanism_ord == mechanismOrd(QStringLiteral("mc_lit_valid"))
               || key.mechanism_ord == mechanismOrd(QStringLiteral("ring_jb"));
    }

    QJsonArray partnerIsoResponsesJson(const RelationshipSeries& rel,
                                       const QString& channel,
                                       const std::vector<std::size_t>& rows,
                                       CatalogueId id) const {
        QJsonArray out;
        const std::vector<double> contribution = rel.contribution.dense(cadence_.stepCount());
        const std::array<std::pair<std::vector<double>, const char*>, 2> sigmaSources = {{
            {sigmaIsoSeries(), "sigma_total"},
            {componentSeriesT0(sigma_para_), "sigma_para"},
        }};
        for (const auto& sigmaSource : sigmaSources) {
            const QString basis = QString::fromLatin1(sigmaSource.second);
            Candidate c{contribution,
                        sigmaSource.first,
                        QStringLiteral("%1.iso").arg(basis),
                        QStringLiteral("%1|signed_contribution").arg(channel)};
            QJsonObject o = evaluatedTargetJson(id, c, rows);
            o.insert(QStringLiteral("sigma_basis"), basis);
            o.insert(QStringLiteral("gauge_flag"), basis == QStringLiteral("sigma_para"));
            o.insert(QStringLiteral("sigma_component"), QStringLiteral("iso"));
            o.insert(QStringLiteral("driver_component"), QStringLiteral("signed_contribution"));
            out.append(o);
        }
        return out;
    }

    QJsonObject partnerTensorRevealJson(const RelationshipKey& key,
                                        const RelationshipSeries& rel,
                                        int index,
                                        const std::vector<std::size_t>& rows) const {
        const Mat3Series kernel = rel.kernel_mol.dense(cadence_.stepCount());
        const QString channel =
            QStringLiteral("partner%1.%2.source%3")
                .arg(index)
                .arg(mechanismName(key.mechanism_ord))
                .arg(key.source_id);
        const CatalogueId id = key.mechanism_ord == mechanismOrd(QStringLiteral("ring_jb"))
                                   ? CatalogueId::RING_H
                                   : CatalogueId::MCCONNELL;
        QJsonObject o;
        o.insert(QStringLiteral("level"), QStringLiteral("per_partner"));
        o.insert(QStringLiteral("frame"), QStringLiteral("molecular"));
        o.insert(QStringLiteral("projection"), QStringLiteral("T_mol = R^T * T_lab * R"));
        o.insert(QStringLiteral("present_frames"), presentCount(kernel));
        o.insert(QStringLiteral("kernel_component_summaries"), tensorComponentSummariesJson(kernel));
        o.insert(QStringLiteral("component_responses"),
                 componentResponsesJson(kernel, channel, rows, id, false));
        o.insert(QStringLiteral("iso_responses"),
                 partnerIsoResponsesJson(rel, channel, rows, id));
        o.insert(QStringLiteral("magic_angle_crossing"),
                 magicAngleSummaryJson(rel.cos_theta.dense(cadence_.stepCount()),
                                       rel.dipolar.dense(cadence_.stepCount())));
        o.insert(QStringLiteral("sigma_follow_flip"), sigmaFlipFollowJson(rel));
        return o;
    }

    QJsonArray perPartnerTensorRevealsJson(const std::vector<std::size_t>& rows) const {
        QJsonArray out;
        int idx = 0;
        for (const auto& item : relationships_) {
            if (!isAngularPartner(item.first)) continue;
            QJsonObject o = relationshipJson(item.first, item.second, idx);
            o.insert(QStringLiteral("tensorial_reveal"),
                     partnerTensorRevealJson(item.first, item.second, idx, rows));
            out.append(o);
            ++idx;
        }
        return out;
    }

    QJsonObject buildTensorialRevealsJson() const {
        const std::vector<std::size_t> rows = cadence_.sigmaRows();
        QJsonObject o;
        o.insert(QStringLiteral("projection"), QStringLiteral("T_mol = R^T * T_lab * R"));
        o.insert(QStringLiteral("sigma_targets"),
                 QJsonArray{QStringLiteral("sigma_total.molcomp:*"),
                            QStringLiteral("sigma_para.molcomp:*")});
        o.insert(QStringLiteral("para_gauge_flagged"), true);
        o.insert(QStringLiteral("per_type"), perTypeTensorRevealsJson(rows));
        o.insert(QStringLiteral("per_partner"), perPartnerTensorRevealsJson(rows));
        return o;
    }

    QJsonObject responseJson(const EvaluatedResponse& e) const {
        QJsonObject o;
        o.insert(QStringLiteral("catalogue_id"), QString::fromLatin1(catalogueRow(e.id).id_name));
        o.insert(QStringLiteral("driver"), QStringLiteral("%1 [%2]")
                                               .arg(e.driver_channel, QString::fromLatin1(catalogueRow(e.id).driver_family)));
        o.insert(QStringLiteral("frame"), frameTagName(e.frame));
        QJsonArray targets;
        for (const EvaluatedTarget& t : e.targets) targets.append(targetJson(e.id, t));
        o.insert(QStringLiteral("targets"), targets);
        const EvaluatedTarget* primary = primaryTarget(e);
        QJsonObject ll;
        ll.insert(QStringLiteral("best_lag"), primary ? primary->best_lag : 0);
        ll.insert(QStringLiteral("lead_r"), primary ? jd(primary->lead_r) : QJsonValue(QJsonValue::Null));
        o.insert(QStringLiteral("lead_lag"), ll);
        QJsonArray considerations;
        for (const QString& note : e.considerations) considerations.append(note);
        o.insert(QStringLiteral("considerations"), considerations);
        return o;
    }

    // Build a single context entry's bounded responses[] plus its within-atom chains.
    QJsonObject buildContextJson(int hybOrd,
                                 const QString& contactedCategory,
                                 const QString& graphStratum,
                                 const std::vector<std::size_t>& rows,
                                 std::size_t& characterizedOut,
                                 std::set<QString>& serialDriversOut,
                                 std::set<QString>& serialSigmaOut,
                                 bool forceComputed) const {
        QJsonObject ctx;
        QJsonObject key;
        key.insert(QStringLiteral("hybridisation"),
                   hybridName(static_cast<model::Hybridisation>(hybOrd)));
        key.insert(QStringLiteral("hybridisation_ord"), hybOrd);
        const model::QtAtom& a = body_.run.protein->atom(atom_);
        key.insert(QStringLiteral("element"), elementName(a.element));
        key.insert(QStringLiteral("contacted_category"), contactedCategory);
        key.insert(QStringLiteral("graph_stratum"), graphStratum);
        key.insert(QStringLiteral("n_frames_in_context"), static_cast<int>(rows.size()));
        ctx.insert(QStringLiteral("context_key"), key);

        const bool lowSupport =
            !forceComputed && rows.size() < static_cast<std::size_t>(kMinContextFrames);
        ctx.insert(QStringLiteral("context_status"),
                   lowSupport ? QStringLiteral("COMPUTED_LOW_SUPPORT") : QStringLiteral("COMPUTED"));

        std::vector<EvaluatedResponse> characterized;
        int applicableChannels = 0;
        for (const CatalogueRow& row : kCharacterizationCatalogue) {
            if (row.id == CatalogueId::WIBERG) continue;  // structured object slot, not responses[].
            if (!membershipHolds(row.id)) continue;
            ++applicableChannels;
            auto e = evaluateCatalogueRow(row.id, rows);
            if (!e) {
                // Membership held but no computable candidate: recorded as a flat
                // characterization, not dropped.
                EvaluatedResponse miss;
                miss.id = row.id;
                miss.frame = row.frame;
                miss.driver_channel = QString::fromLatin1(row.driver_family);
                miss.considerations.append(QStringLiteral("ZERO_VARIANCE_OR_INSUFFICIENT_SUPPORT"));
                characterized.push_back(miss);
                continue;
            }
            characterized.push_back(*e);
        }
        const int nOverCap =
            std::max(0, static_cast<int>(characterized.size()) - applicableChannels);
        QJsonArray responses;
        for (const EvaluatedResponse& e : characterized) {
            responses.append(responseJson(e));
            for (const EvaluatedTarget& t : e.targets) serialSigmaOut.insert(t.sigma_target);
            serialDriversOut.insert(e.driver_channel);
        }
        ctx.insert(QStringLiteral("responses"), responses);
        ctx.insert(QStringLiteral("n_over_cap"), nOverCap);
        characterizedOut = characterized.size();

        // considerations[] -- notes that ride with the characterization.
        QJsonArray considerations;
        if (lowSupport) {
            considerations.append(
                QStringLiteral("LOW_SUPPORT: n_frames_in_context < kMinContextFrames; responses and chains computed with existing support guards"));
        }
        considerations.append(QStringLiteral("gauge(dia/para): GIAO total is gauge-invariant; dia/para is context colour only"));
        ctx.insert(QStringLiteral("considerations"), considerations);

        QJsonObject chains;
        chains.insert(QStringLiteral("collinearity"), driverCollinearityChains(characterized, rows));
        chains.insert(QStringLiteral("mediations"), mediationChainsJson(characterized, rows));
        ctx.insert(QStringLiteral("chains"), chains);
        return ctx;
    }

    QJsonArray driverCollinearityChains(const std::vector<EvaluatedResponse>& characterized,
                                        const std::vector<std::size_t>& rows) const {
        auto hasAny = [&](std::initializer_list<CatalogueId> ids) {
            for (const EvaluatedResponse& e : characterized)
                for (CatalogueId id : ids)
                    if (e.id == id) return true;
            return false;
        };
        const bool hasField = hasAny({CatalogueId::FIELD_ISO,
                                      CatalogueId::FIELD_SPAN,
                                      CatalogueId::FIELD_MOLCOMP});
        const bool hasEfg = hasAny({CatalogueId::EFG_ISO,
                                    CatalogueId::EFG_MOLCOMP,
                                    CatalogueId::EFG_EIGEN});
        QJsonArray out;
        if (!hasField || !hasEfg) return out;

        std::vector<double> x;
        std::vector<double> y;
        pairedOnSigmaRows(fieldAbsSeries(), efgAbsSeries(), rows, x, y);
        if (x.size() < 4) return out;
        const double r = pearsonR(x, y);
        if (!finite(r)) return out;
        QJsonObject o;
        o.insert(QStringLiteral("driver_a"), QStringLiteral("field.mopac_coulomb|abs_E"));
        o.insert(QStringLiteral("driver_b"), QStringLiteral("efg|abs_T2"));
        o.insert(QStringLiteral("field_efg_pearson_r_unthresholded"), jd(r));
        o.insert(QStringLiteral("finite_n"), static_cast<int>(x.size()));
        o.insert(QStringLiteral("family_a"), QStringLiteral("field_sources"));
        o.insert(QStringLiteral("family_b"), QStringLiteral("efg_sources"));
        o.insert(QStringLiteral("subspace_compare"),
                 subspaceCompareObject(QStringLiteral("field_sources_vs_efg_node"),
                                       fieldSourcesFamily(),
                                       efgNodeFamily(),
                                       rows));
        o.insert(QStringLiteral("threshold"), kCollinearThreshold);
        o.insert(QStringLiteral("above_threshold"), std::abs(r) > kCollinearThreshold);
        o.insert(QStringLiteral("emit_policy"), QStringLiteral("unthresholded"));
        out.append(o);
        return out;
    }

    QJsonObject buildAccumulatorJson() const {
        QJsonObject acc;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& a = p.atom(atom_);

        // atom_key
        QJsonObject atomKey;
        atomKey.insert(QStringLiteral("uid"), uid(QStringLiteral("atom"), atom_));
        atomKey.insert(QStringLiteral("atom_index"), static_cast<int>(atom_));
        atomKey.insert(QStringLiteral("atom_name"), p.atomLabel(atom_, model::NamingConvention::Iupac));
        atomKey.insert(QStringLiteral("element"), elementName(a.element));
        atomKey.insert(QStringLiteral("residue_index"), a.residueIndex);
        if (validResidue(p, a.residueIndex)) {
            const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
            atomKey.insert(QStringLiteral("residue_number"), r.address.residueNumber);
            atomKey.insert(QStringLiteral("amino_acid"), aaName(r.aminoAcid));
        }
        atomKey.insert(QStringLiteral("backbone_role"), backboneRoleName(a.backboneRole));
        acc.insert(QStringLiteral("atom_key"), atomKey);

        // sigma0_ref -- the per-IUPAC zero-field reference handle (mean iso).
        double isoMean = kNaN;
        {
            std::vector<double> iso;
            for (double v : sigmaIsoSeries()) if (finite(v)) iso.push_back(v);
            if (!iso.empty()) isoMean = bmstats::mean(iso);
        }
        acc.insert(QStringLiteral("sigma0_ref"), jd(isoMean));

        // contexts[] -- emit the whole-trajectory BASE context first (§2.4B),
        // then additional hybridisation/contact-class stratified contexts where
        // the trajectory actually splits into regimes.
        const int stratum = graphStratum();
        std::map<std::pair<int, int>, std::vector<std::size_t>> rowsByHybContact;
        for (std::size_t row : cadence_.sigmaRows()) {
            int hyb = static_cast<int>(model::Hybridisation::Unassigned);
            if (row < hybridisation_.values.size() && hybridisation_.values[row] != IntSeries::kMissing)
                hyb = hybridisation_.values[row];
            else
                hyb = static_hybridisation_ord_;
            const std::uint8_t membership = contactedClassMaskAtStep(row);
            bool emittedContactClass = false;
            for (ContactedClass c : kPresentContactClasses) {
                if ((membership & contactedClassBit(c)) == 0) continue;
                rowsByHybContact[{hyb, static_cast<int>(c)}].push_back(row);
                emittedContactClass = true;
            }
            if (!emittedContactClass)
                rowsByHybContact[{hyb, static_cast<int>(ContactedClass::NONE)}].push_back(row);
        }
        if (rowsByHybContact.empty())
            rowsByHybContact[{static_hybridisation_ord_, static_cast<int>(ContactedClass::NONE)}] =
                cadence_.sigmaRows();

        QJsonArray contexts;
        std::set<QString> serialDrivers;
        std::set<QString> serialSigma;
        std::size_t totalResponses = 0;
        {
            std::size_t characterizedCount = 0;
            QJsonObject ctx = buildContextJson(static_hybridisation_ord_,
                                               QStringLiteral("WHOLE_TRAJECTORY"),
                                               QStringLiteral("ALL"),
                                               cadence_.sigmaRows(), characterizedCount,
                                               serialDrivers, serialSigma, true);
            contexts.append(ctx);
            totalResponses += characterizedCount;
        }
        for (const auto& kv : rowsByHybContact) {
            std::size_t characterizedCount = 0;
            QJsonObject ctx = buildContextJson(kv.first.first,
                                               contactedClassName(
                                                   static_cast<ContactedClass>(kv.first.second)),
                                               QString::number(stratum), kv.second,
                                               characterizedCount,
                                               serialDrivers, serialSigma, false);
            contexts.append(ctx);
            totalResponses += characterizedCount;
        }
        acc.insert(QStringLiteral("contexts"), contexts);
        acc.insert(QStringLiteral("field_reversal_tail_diagnostics"),
                   fieldReversalTailDiagnosticsJson());

        // serial[] -- one entry per (sigma_target ∪ driver) a characterized
        // response touched (§4.2; bounded by the catalogue, never all 85 channels).
        QJsonArray serial;
        auto addSerial = [&](const QString& ref, const QString& kind, const std::vector<double>& series) {
            RunningSeriesRef sref;
            sref.series_ref = ref;
            sref.values = series;
            QJsonObject o = serialRecord(sref);
            o.insert(QStringLiteral("target_ref_kind"), kind);
            const DwellResult d = dwellStats(series);
            QJsonObject dwell;
            dwell.insert(QStringLiteral("mean_dwell_frames"), jd(d.mean_dwell_frames));
            dwell.insert(QStringLiteral("n_transitions"), d.n_transitions);
            dwell.insert(QStringLiteral("autocorr_time"), jd(d.autocorr_time));
            o.insert(QStringLiteral("dwell"), dwell);
            serial.append(o);
        };
        // the sigma iso/aniso targets (a sigma-response always present).
        addSerial(QStringLiteral("sigma.iso"), QStringLiteral("SIGMA"), sigmaIsoSeries());
        if (serialSigma.count(QStringLiteral("invariant:aniso")) || serialSigma.empty())
            addSerial(QStringLiteral("sigma.aniso"), QStringLiteral("SIGMA"),
                      sigmaInvariantSeries(QStringLiteral("aniso")));
        // the characterized drivers' own dynamics (bounded set).
        if (serialDrivers.count(QStringLiteral("field.mopac_coulomb|abs_E")))
            addSerial(QStringLiteral("field.abs_E"), QStringLiteral("DRIVER"), fieldAbsSeries());
        bool anyEfgDriver = false;
        for (const QString& d : serialDrivers) if (d.startsWith(QStringLiteral("efg"))) anyEfgDriver = true;
        if (anyEfgDriver) addSerial(QStringLiteral("efg.abs_T2"), QStringLiteral("DRIVER"), efgAbsSeries());
        acc.insert(QStringLiteral("serial"), serial);

        // dia_para_split[] -- per characterized FIELD/EFG driver; gauge colour only.
        QJsonArray diaPara;
        {
            // for the headline field/EFG drivers, does the response live in dia/para?
            const std::vector<double> diaIso = componentSeriesT0(sigma_dia_);
            const std::vector<double> paraIso = componentSeriesT0(sigma_para_);
            const std::vector<double> field = fieldAbsSeries();
            std::vector<double> x;
            std::vector<double> y;
            pairedOnSigmaRows(diaIso, field, cadence_.sigmaRows(), x, y);
            const double rDia = pearsonR(x, y);
            pairedOnSigmaRows(paraIso, field, cadence_.sigmaRows(), x, y);
            const double rPara = pearsonR(x, y);
            if (finite(rDia) || finite(rPara)) {
                QJsonObject o;
                o.insert(QStringLiteral("driver"), QStringLiteral("field.mopac_coulomb|abs_E"));
                o.insert(QStringLiteral("moves_dia"), jd(rDia));
                o.insert(QStringLiteral("moves_para"), jd(rPara));
                o.insert(QStringLiteral("gauge_caveat"), true);
                diaPara.append(o);
            }
        }
        acc.insert(QStringLiteral("dia_para_split"), diaPara);

        // recurrence[] -- per hybridisation_change_event: the sigma response and
        // whether it recurred (§6.1.5). Most atoms carry [].
        QJsonArray recurrence;
        {
            std::vector<std::size_t> eventRows;
            int prev = IntSeries::kMissing;
            for (std::size_t i = 0; i < hybridisation_.values.size(); ++i) {
                const int cur = hybridisation_.values[i];
                if (cur == IntSeries::kMissing) continue;
                if (prev != IntSeries::kMissing && cur != prev) eventRows.push_back(i);
                prev = cur;
            }
            const std::vector<double> iso = sigmaIsoSeries();
            for (std::size_t k = 0; k < eventRows.size(); ++k) {
                QJsonObject o;
                o.insert(QStringLiteral("event_frame"), static_cast<int>(eventRows[k]));
                // sigma response = the iso step magnitude across the nearest sigma rows.
                double before = kNaN;
                double after = kNaN;
                for (std::size_t r = eventRows[k]; r-- > 0;) if (finite(iso[r])) { before = iso[r]; break; }
                for (std::size_t r = eventRows[k]; r < iso.size(); ++r) if (finite(iso[r])) { after = iso[r]; break; }
                const double step = (finite(before) && finite(after)) ? (after - before) : kNaN;
                o.insert(QStringLiteral("sigma_response"), jd(step));
                o.insert(QStringLiteral("recurred"), eventRows.size() > 1);
                recurrence.append(o);
            }
        }
        acc.insert(QStringLiteral("recurrence"), recurrence);

        acc.insert(QStringLiteral("molecular_frame_audit"), molecularFrameAuditJson());
        acc.insert(QStringLiteral("tensorial_reveals"), buildTensorialRevealsJson());

        // Wiberg is a structured slot (§2.1), not a flat responses[] row.
        acc.insert(QStringLiteral("wiberg"), buildWibergJson());

        // companion: the individual contacted residues (§2.4), un-folded.
        QJsonArray contactedCompanion;
        for (const auto& r : contacted_residues_) {
            QJsonObject o;
            o.insert(QStringLiteral("residue_number"), r.first);
            o.insert(QStringLiteral("amino_acid"), aaName(static_cast<model::AminoAcid>(r.second)));
            o.insert(QStringLiteral("contacted_class"),
                     contactedClassName(classOf(static_cast<model::AminoAcid>(r.second))));
            contactedCompanion.append(o);
        }
        acc.insert(QStringLiteral("contacted_residues"), contactedCompanion);

        last_accumulator_response_count = totalResponses;
        last_accumulator_context_count = static_cast<std::size_t>(contexts.size());

        // -- the four structural asserts (§5.1.3), fail-loud.
        assertAccumulatorStructure(acc);
        return acc;
    }

    // ACCUMULATOR_RESPEC §5.1.3 anti-regression asserts (fail-loud on the
    // emitted object). These guard the doctrine structurally.
    void assertAccumulatorStructure(const QJsonObject& acc) const {
        const auto containsKey = [](const QJsonValue& root, const QString& key) {
            std::function<bool(const QJsonValue&)> walk = [&](const QJsonValue& v) -> bool {
                if (v.isObject()) {
                    const QJsonObject o = v.toObject();
                    if (o.contains(key)) return true;
                    for (auto it = o.begin(); it != o.end(); ++it)
                        if (walk(it.value())) return true;
                } else if (v.isArray()) {
                    const QJsonArray a = v.toArray();
                    for (const QJsonValue& item : a)
                        if (walk(item)) return true;
                }
                return false;
            };
            return walk(root);
        };

        // 1. no `boost` key anywhere.
        if (containsKey(acc, QStringLiteral("boost")))
            qFatal("accumulator regression: atom %lld emitted a 'boost' key",
                   static_cast<long long>(atom_));
        if (containsKey(acc, QStringLiteral("grade")))
            qFatal("accumulator regression: atom %lld emitted a 'grade' key",
                   static_cast<long long>(atom_));

        std::set<QString> expectedIds;
        for (const CatalogueRow& row : kCharacterizationCatalogue) {
            if (row.id == CatalogueId::WIBERG) continue;
            if (membershipHolds(row.id))
                expectedIds.insert(QString::fromLatin1(row.id_name));
        }

        const QJsonArray contexts = acc.value(QStringLiteral("contexts")).toArray();
        for (const QJsonValue& cv : contexts) {
            const QJsonObject ctx = cv.toObject();
            const QJsonArray responses = ctx.value(QStringLiteral("responses")).toArray();
            const QString status = ctx.value(QStringLiteral("context_status")).toString();
            if (responses.size() > kMaxChannels)
                qFatal("accumulator regression: atom %lld context has %d responses > structural cap %d",
                       static_cast<long long>(atom_), static_cast<int>(responses.size()), kMaxChannels);
            const int nOverCap = ctx.value(QStringLiteral("n_over_cap")).toInt(0);
            if (nOverCap != 0)
                qFatal("accumulator regression: atom %lld context has n_over_cap=%d",
                       static_cast<long long>(atom_), nOverCap);
            if (!ctx.contains(QStringLiteral("chains")) || !ctx.value(QStringLiteral("chains")).isObject())
                qFatal("accumulator regression: atom %lld context missing bounded chains slot",
                       static_cast<long long>(atom_));
            const QJsonObject chains = ctx.value(QStringLiteral("chains")).toObject();
            if (!chains.value(QStringLiteral("collinearity")).isArray()
                || !chains.value(QStringLiteral("mediations")).isArray())
                qFatal("accumulator regression: atom %lld context chains malformed",
                       static_cast<long long>(atom_));
            const int nCollinear = chains.value(QStringLiteral("collinearity")).toArray().size();
            const int nMediations = chains.value(QStringLiteral("mediations")).toArray().size();
            if (nCollinear > kMaxCollinearityChains || nMediations > kMaxMediationChains)
                qFatal("accumulator regression: atom %lld context chains are unbounded (%d collinearity, %d mediations)",
                       static_cast<long long>(atom_), nCollinear, nMediations);

            std::set<QString> seenIds;
            for (const QJsonValue& rv : responses) {
                const QJsonObject r = rv.toObject();
                if (!r.contains(QStringLiteral("catalogue_id"))
                    || r.value(QStringLiteral("catalogue_id")).toString().isEmpty())
                    qFatal("accumulator regression: atom %lld response missing catalogue_id",
                           static_cast<long long>(atom_));
                if (r.contains(QStringLiteral("kept_reason")))
                    qFatal("accumulator regression: atom %lld response still emits kept_reason",
                           static_cast<long long>(atom_));
                if (r.contains(QStringLiteral("mediation")))
                    qFatal("accumulator regression: atom %lld response still emits per-entry mediation",
                           static_cast<long long>(atom_));
                const QString cid = r.value(QStringLiteral("catalogue_id")).toString();
                if (!expectedIds.count(cid))
                    qFatal("accumulator regression: atom %lld response has non-applicable catalogue_id '%s'",
                           static_cast<long long>(atom_), qPrintable(cid));
                if (!seenIds.insert(cid).second)
                    qFatal("accumulator regression: atom %lld duplicate catalogue_id '%s'",
                           static_cast<long long>(atom_), qPrintable(cid));
                const QJsonArray responseConsiderations =
                    r.value(QStringLiteral("considerations")).toArray();
                bool insufficientSupport = false;
                for (const QJsonValue& note : responseConsiderations) {
                    if (note.toString() == QStringLiteral("ZERO_VARIANCE_OR_INSUFFICIENT_SUPPORT")) {
                        insufficientSupport = true;
                        break;
                    }
                }
                if ((!r.contains(QStringLiteral("targets")) || !r.value(QStringLiteral("targets")).isArray()
                     || r.value(QStringLiteral("targets")).toArray().isEmpty())
                    && !insufficientSupport)
                    qFatal("accumulator regression: atom %lld response '%s' missing targets[]",
                           static_cast<long long>(atom_), qPrintable(cid));
            }
            if ((status == QStringLiteral("COMPUTED")
                 || status == QStringLiteral("COMPUTED_LOW_SUPPORT"))
                && seenIds != expectedIds)
                qFatal("accumulator regression: atom %lld context has %d catalogue entries, expected %d",
                       static_cast<long long>(atom_), static_cast<int>(seenIds.size()),
                       static_cast<int>(expectedIds.size()));
        }
    }

    void buildStaticBonds() {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return;
        const model::QtTopology& topo = body_.run.protein->topology();
        for (int32_t bi : topo.bondIndicesForAtom(atom_)) {
            if (bi < 0 || static_cast<std::size_t>(bi) >= topo.bondCount()) continue;
            const model::QtBond& b = topo.bondAt(static_cast<std::size_t>(bi));
            if (b.atomIndexA < 0 || b.atomIndexB < 0) continue;
            const std::size_t a = static_cast<std::size_t>(b.atomIndexA);
            const std::size_t bb = static_cast<std::size_t>(b.atomIndexB);
            const std::size_t other = a == atom_ ? bb : a;
            const std::uint64_t key = pairKey(a, bb);
            static_bonds_[key] = {b.bondIndex, b.order, b.category, other};
            bond_length_series_.try_emplace(key, cadence_.stepCount());
            bond_length_rejected_degenerate_.try_emplace(key, cadence_.stepCount());
        }
    }

    void readStaticHybridisation() {
        if (body_.catalog.present(body_, io::FieldKind::EnrichmentHybridisation, atom_, 0, -1)) {
            const std::optional<double> h =
                body_.catalog.value(body_, io::FieldKind::EnrichmentHybridisation, atom_, 0, -1);
            if (h && finite(*h)) static_hybridisation_ord_ = static_cast<int>(*h);
        }
    }

    enum class MolFrameKind {
        None,
        BackboneCarbonyl,
        BackboneAmideN,
        BackboneAmideH,
        AromaticRingLocal,
        MetSd,
        SidechainCarboxylate,
        SidechainGuanidinium,
        SidechainCarboxamide,
        BackboneCA,
        AliphaticCarbon,
        HydroxylOxygen,
    };

    struct MolFrameSpec {
        MolFrameKind kind = MolFrameKind::None;
        QString label = QStringLiteral("none");
        std::vector<int32_t> anchors;
        int32_t origin = -1;
        int32_t x_anchor = -1;
        int32_t plane_anchor = -1;
        int32_t second_anchor = -1;
        int32_t ring = -1;
        int32_t heavy = -1;
        QString ring_type;
        FrameVariant variant = FrameVariant::Invalid;
    };

    bool validAtomIndex(int32_t atomIndex) const {
        return atomIndex >= 0 && body_.run.protein
               && static_cast<std::size_t>(atomIndex) < body_.run.protein->atomCount();
    }

    std::optional<int32_t> atomByResidueName(int32_t residueIndex, const QString& name) const {
        if (!body_.run.protein || !validResidue(*body_.run.protein, residueIndex)) return std::nullopt;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtResidue& r = p.residue(static_cast<std::size_t>(residueIndex));
        for (int32_t atomIndex : r.atomIndices) {
            if (!validAtomIndex(atomIndex)) continue;
            if (p.atomLabel(static_cast<std::size_t>(atomIndex), model::NamingConvention::Iupac) == name)
                return atomIndex;
        }
        return std::nullopt;
    }

    std::optional<int32_t> ringContainingAtom(int32_t atomIndex) const {
        if (!validAtomIndex(atomIndex) || !body_.run.protein) return std::nullopt;
        const model::QtTopology& topo = body_.run.protein->topology();
        for (std::size_t ring = 0; ring < topo.ringCount(); ++ring) {
            const model::QtRing& r = topo.ringAt(ring);
            if (ringTypeName(r.TypeIndex()).endsWith(QStringLiteral("Perimeter"))) continue;
            if (std::find(r.atomIndices.begin(), r.atomIndices.end(), atomIndex) != r.atomIndices.end())
                return static_cast<int32_t>(ring);
        }
        return std::nullopt;
    }

    std::optional<int32_t> molecularFrameRingForAtom(const model::QtAtom& atom) const {
        if (const auto direct = ringContainingAtom(static_cast<int32_t>(atom_))) return direct;
        if (atom.parentAtomIndex >= 0) return ringContainingAtom(atom.parentAtomIndex);
        return std::nullopt;
    }

    bool hasResidueAtom(int32_t residueIndex, const QString& name) const {
        return atomByResidueName(residueIndex, name).has_value();
    }

    std::vector<int32_t> bondedAtomIndices(bool heavyOnly = false) const {
        std::vector<int32_t> out;
        if (!body_.run.protein) return out;
        for (const auto& item : static_bonds_) {
            const int32_t other = static_cast<int32_t>(item.second.other);
            if (!validAtomIndex(other)) continue;
            if (heavyOnly
                && body_.run.protein->atom(static_cast<std::size_t>(other)).element
                       == model::Element::H)
                continue;
            if (std::find(out.begin(), out.end(), other) == out.end()) out.push_back(other);
        }
        std::sort(out.begin(), out.end());
        return out;
    }

    std::optional<int32_t> firstBondedAtomByElement(model::Element element) const {
        if (!body_.run.protein) return std::nullopt;
        for (int32_t ai : bondedAtomIndices(false)) {
            if (body_.run.protein->atom(static_cast<std::size_t>(ai)).element == element)
                return ai;
        }
        return std::nullopt;
    }

    bool isHydroxylOxygen(const model::QtAtom& atom, const QString& name) const {
        if (atom.element != model::Element::O) return false;
        if (oxygenExpectedSp2(atom)) return false;
        if (atom.ffAtomType == model::QtFfAtomType::OH
            || atom.planarGroup == model::PlanarGroupKind::AromaticHydroxyl)
            return true;
        return name == QStringLiteral("OG") || name == QStringLiteral("OG1")
               || name == QStringLiteral("OH");
    }

    bool isAliphaticCarbon(const model::QtAtom& atom, const QString& name) const {
        if (atom.element != model::Element::C) return false;
        if (name == QStringLiteral("C") || name == QStringLiteral("CA")) return false;
        if (atom.backboneRole == model::BackboneRole::CarbonylCarbon) return false;
        if (atom.IsInAnyRing()) return false;
        return atom.planarGroup == model::PlanarGroupKind::None
               || static_hybridisation_ord_ == static_cast<int>(model::Hybridisation::sp3)
               || name.startsWith(QStringLiteral("C"));
    }

    std::optional<std::array<QString, 3>> carboxylateGroup(const QString& name, int32_t residueIndex) const {
        auto present = [&](const QString& atomName) { return hasResidueAtom(residueIndex, atomName); };
        if ((name == QStringLiteral("CG") || name == QStringLiteral("OD1") || name == QStringLiteral("OD2"))
            && present(QStringLiteral("CG")) && present(QStringLiteral("OD1")) && present(QStringLiteral("OD2")))
            return std::array<QString, 3>{QStringLiteral("CG"), QStringLiteral("OD1"), QStringLiteral("OD2")};
        if ((name == QStringLiteral("CD") || name == QStringLiteral("OE1") || name == QStringLiteral("OE2"))
            && present(QStringLiteral("CD")) && present(QStringLiteral("OE1")) && present(QStringLiteral("OE2")))
            return std::array<QString, 3>{QStringLiteral("CD"), QStringLiteral("OE1"), QStringLiteral("OE2")};
        return std::nullopt;
    }

    std::optional<std::array<QString, 4>> guanidiniumGroup(const QString& name, int32_t residueIndex) const {
        auto present = [&](const QString& atomName) { return hasResidueAtom(residueIndex, atomName); };
        if ((name == QStringLiteral("CZ") || name == QStringLiteral("NE") || name == QStringLiteral("NH1")
             || name == QStringLiteral("NH2"))
            && present(QStringLiteral("CZ")) && present(QStringLiteral("NE")) && present(QStringLiteral("NH1"))
            && present(QStringLiteral("NH2")))
            return std::array<QString, 4>{QStringLiteral("CZ"), QStringLiteral("NE"), QStringLiteral("NH1"),
                                          QStringLiteral("NH2")};
        return std::nullopt;
    }

    std::optional<std::array<QString, 3>> carboxamideGroup(const QString& name, int32_t residueIndex) const {
        auto present = [&](const QString& atomName) { return hasResidueAtom(residueIndex, atomName); };
        if (name == QStringLiteral("ND2") && present(QStringLiteral("CG")) && present(QStringLiteral("OD1"))
            && present(QStringLiteral("ND2")) && !present(QStringLiteral("OD2")))
            return std::array<QString, 3>{QStringLiteral("CG"), QStringLiteral("OD1"), QStringLiteral("ND2")};
        if (name == QStringLiteral("NE2") && present(QStringLiteral("CD")) && present(QStringLiteral("OE1"))
            && present(QStringLiteral("NE2")) && !present(QStringLiteral("OE2")))
            return std::array<QString, 3>{QStringLiteral("CD"), QStringLiteral("OE1"), QStringLiteral("NE2")};
        if (name == QStringLiteral("CG") && present(QStringLiteral("OD1")) && present(QStringLiteral("ND2"))
            && !present(QStringLiteral("OD2")))
            return std::array<QString, 3>{QStringLiteral("CG"), QStringLiteral("OD1"), QStringLiteral("ND2")};
        if (name == QStringLiteral("OD1") && present(QStringLiteral("CG")) && present(QStringLiteral("ND2"))
            && !present(QStringLiteral("OD2")))
            return std::array<QString, 3>{QStringLiteral("CG"), QStringLiteral("OD1"), QStringLiteral("ND2")};
        if (name == QStringLiteral("CD") && present(QStringLiteral("OE1")) && present(QStringLiteral("NE2"))
            && !present(QStringLiteral("OE2")))
            return std::array<QString, 3>{QStringLiteral("CD"), QStringLiteral("OE1"), QStringLiteral("NE2")};
        if (name == QStringLiteral("OE1") && present(QStringLiteral("CD")) && present(QStringLiteral("NE2"))
            && !present(QStringLiteral("OE2")))
            return std::array<QString, 3>{QStringLiteral("CD"), QStringLiteral("OE1"), QStringLiteral("NE2")};
        return std::nullopt;
    }

    void applyMolecularFrameSpec(const MolFrameSpec& spec) {
        molecular_frame_spec_ = spec;
        molecular_frame_.setMetadata(spec.label, spec.anchors,
                                     spec.variant,
                                     spec.ring >= 0 ? uid(QStringLiteral("ring"), static_cast<std::size_t>(spec.ring))
                                                    : QString(),
                                     spec.ring_type);
    }

    void configureMolecularFrame() {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return;
        const model::QtProtein& p = *body_.run.protein;
        const model::QtAtom& atom = p.atom(atom_);
        const int32_t ri = atom.residueIndex;
        if (!validResidue(p, ri)) return;
        const QString name = p.atomLabel(atom_, model::NamingConvention::Iupac);
        auto recAt = [&](const QString& atomName, int32_t residueIndex) {
            return atomByResidueName(residueIndex, atomName);
        };
        auto rec = [&](const QString& atomName) { return recAt(atomName, ri); };

        if (name == QStringLiteral("CA")) {
            const auto ca = rec(QStringLiteral("CA"));
            const auto n = rec(QStringLiteral("N"));
            const auto c = rec(QStringLiteral("C"));
            if (ca && n && c) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::BackboneCA;
                spec.label = QStringLiteral("backbone_ca");
                spec.origin = *ca;
                spec.x_anchor = *n;
                spec.plane_anchor = *c;
                spec.anchors = std::vector<int32_t>{*ca, *n, *c};
                spec.variant = FrameVariant::BackboneCA;
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (atom.backboneRole == model::BackboneRole::CarbonylCarbon
            || atom.backboneRole == model::BackboneRole::CarbonylOxygen) {
            const auto c = rec(QStringLiteral("C"));
            const auto o = rec(QStringLiteral("O"));
            const auto ca = rec(QStringLiteral("CA"));
            const auto nNext = recAt(QStringLiteral("N"), ri + 1);
            if (c && o && ca) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::BackboneCarbonyl;
                spec.label = QStringLiteral("backbone_carbonyl");
                spec.origin = *c;
                spec.x_anchor = *o;
                spec.plane_anchor = nNext ? *nNext : *ca;
                spec.anchors = std::vector<int32_t>{*c, *o, *ca};
                if (nNext) spec.anchors.push_back(*nNext);
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (atom.backboneRole == model::BackboneRole::Nitrogen || name == QStringLiteral("N")) {
            const auto n = rec(QStringLiteral("N"));
            const auto h = rec(QStringLiteral("H"));
            const auto ca = rec(QStringLiteral("CA"));
            const auto cPrev = recAt(QStringLiteral("C"), ri - 1);
            if (n && ca && (h || cPrev)) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::BackboneAmideN;
                spec.label = QStringLiteral("backbone_amide_n");
                spec.origin = *n;
                spec.x_anchor = h ? *h : *ca;
                spec.plane_anchor = cPrev ? *cPrev : *ca;
                spec.anchors = std::vector<int32_t>{*n, *ca};
                if (h) spec.anchors.push_back(*h);
                if (cPrev) spec.anchors.push_back(*cPrev);
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (atom.element == model::Element::H
            && atom.backboneRole == model::BackboneRole::AmideHydrogen
            && atom.parentAtomIndex >= 0
            && validAtomIndex(atom.parentAtomIndex)
            && p.atom(static_cast<std::size_t>(atom.parentAtomIndex)).backboneRole
                   == model::BackboneRole::Nitrogen) {
            const auto n = rec(QStringLiteral("N"));
            const auto ca = rec(QStringLiteral("CA"));
            const auto cPrev = recAt(QStringLiteral("C"), ri - 1);
            if (n && ca && *n == atom.parentAtomIndex) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::BackboneAmideH;
                spec.label = QStringLiteral("backbone_amide_h");
                spec.origin = *n;
                spec.x_anchor = static_cast<int32_t>(atom_);
                spec.second_anchor = *ca;
                spec.plane_anchor = cPrev ? *cPrev : -1;
                spec.anchors = std::vector<int32_t>{*n, static_cast<int32_t>(atom_), *ca};
                if (cPrev) spec.anchors.push_back(*cPrev);
                spec.variant = cPrev ? FrameVariant::HN_Standard : FrameVariant::HN_NTerminus;
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        const bool ringElement = atom.element == model::Element::C || atom.element == model::Element::H
                                 || atom.element == model::Element::N || atom.element == model::Element::O;
        if (ringElement) {
            const auto ring = molecularFrameRingForAtom(atom);
            if (ring && *ring >= 0 && static_cast<std::size_t>(*ring) < p.topology().ringCount()) {
                const model::QtRing& r = p.topology().ringAt(static_cast<std::size_t>(*ring));
                MolFrameSpec spec;
                spec.kind = MolFrameKind::AromaticRingLocal;
                spec.label = QStringLiteral("aromatic_ring_local");
                spec.ring = *ring;
                spec.heavy = atom.element == model::Element::H && atom.parentAtomIndex >= 0
                                 ? atom.parentAtomIndex
                                 : static_cast<int32_t>(atom_);
                spec.anchors = r.atomIndices;
                if (std::find(spec.anchors.begin(), spec.anchors.end(), spec.heavy) == spec.anchors.end())
                    spec.anchors.push_back(spec.heavy);
                spec.ring_type = ringTypeName(r.TypeIndex());
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (atom.element == model::Element::S && name == QStringLiteral("SD")) {
            const auto sd = rec(QStringLiteral("SD"));
            const auto cg = rec(QStringLiteral("CG"));
            const auto ce = rec(QStringLiteral("CE"));
            if (sd && cg && ce) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::MetSd;
                spec.label = QStringLiteral("met_sd_cg_s_ce");
                spec.origin = *sd;
                spec.x_anchor = *ce;
                spec.plane_anchor = *cg;
                spec.anchors = std::vector<int32_t>{*cg, *sd, *ce};
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (const auto group = carboxylateGroup(name, ri)) {
            const auto c = rec((*group)[0]);
            const auto o1 = rec((*group)[1]);
            const auto o2 = rec((*group)[2]);
            if (c && o1 && o2) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::SidechainCarboxylate;
                spec.label = QStringLiteral("sidechain_carboxylate");
                spec.origin = *c;
                spec.x_anchor = *o1;
                spec.second_anchor = *o2;
                spec.plane_anchor = *o1;
                spec.anchors = std::vector<int32_t>{*c, *o1, *o2};
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (const auto group = guanidiniumGroup(name, ri)) {
            const auto cz = rec((*group)[0]);
            const auto ne = rec((*group)[1]);
            const auto nh1 = rec((*group)[2]);
            const auto nh2 = rec((*group)[3]);
            if (cz && ne && nh1) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::SidechainGuanidinium;
                spec.label = QStringLiteral("sidechain_guanidinium");
                spec.origin = *cz;
                spec.x_anchor = *ne;
                spec.plane_anchor = *nh1;
                spec.anchors = std::vector<int32_t>{*cz, *ne, *nh1};
                if (nh2) spec.anchors.push_back(*nh2);
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (const auto group = carboxamideGroup(name, ri)) {
            const auto c = rec((*group)[0]);
            const auto o = rec((*group)[1]);
            const auto n = rec((*group)[2]);
            if (c && o && n) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::SidechainCarboxamide;
                spec.label = QStringLiteral("sidechain_carboxamide");
                spec.origin = *c;
                spec.x_anchor = *o;
                spec.plane_anchor = *n;
                spec.anchors = std::vector<int32_t>{*c, *o, *n};
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (isHydroxylOxygen(atom, name)) {
            const std::vector<int32_t> bonded = bondedAtomIndices(false);
            const auto h = firstBondedAtomByElement(model::Element::H);
            int32_t heavy = atom.parentAtomIndex;
            if (!validAtomIndex(heavy)) {
                const std::vector<int32_t> heavyBonded = bondedAtomIndices(true);
                if (!heavyBonded.empty()) heavy = heavyBonded.front();
            }
            int32_t plane = h.value_or(-1);
            if (!validAtomIndex(plane)) {
                for (int32_t ai : bonded) {
                    if (ai != heavy) { plane = ai; break; }
                }
            }
            if (validAtomIndex(heavy) && validAtomIndex(plane) && heavy != plane) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::HydroxylOxygen;
                spec.label = QStringLiteral("hydroxyl_oxygen");
                spec.origin = static_cast<int32_t>(atom_);
                spec.x_anchor = heavy;
                spec.plane_anchor = plane;
                spec.anchors = std::vector<int32_t>{static_cast<int32_t>(atom_), heavy, plane};
                applyMolecularFrameSpec(spec);
                return;
            }
        }

        if (isAliphaticCarbon(atom, name)) {
            const std::vector<int32_t> bonded = bondedAtomIndices(false);
            const std::vector<int32_t> heavyBonded = bondedAtomIndices(true);
            int32_t xAnchor = atom.parentAtomIndex;
            if (!validAtomIndex(xAnchor) && !heavyBonded.empty()) xAnchor = heavyBonded.front();
            if (!validAtomIndex(xAnchor) && !bonded.empty()) xAnchor = bonded.front();
            int32_t planeAnchor = -1;
            for (int32_t ai : heavyBonded) {
                if (ai != xAnchor) { planeAnchor = ai; break; }
            }
            if (!validAtomIndex(planeAnchor)) {
                for (int32_t ai : bonded) {
                    if (ai != xAnchor) { planeAnchor = ai; break; }
                }
            }
            if (validAtomIndex(xAnchor) && validAtomIndex(planeAnchor) && xAnchor != planeAnchor) {
                MolFrameSpec spec;
                spec.kind = MolFrameKind::AliphaticCarbon;
                spec.label = QStringLiteral("aliphatic_carbon");
                spec.origin = static_cast<int32_t>(atom_);
                spec.x_anchor = xAnchor;
                spec.plane_anchor = planeAnchor;
                spec.anchors = std::vector<int32_t>{static_cast<int32_t>(atom_), xAnchor, planeAnchor};
                applyMolecularFrameSpec(spec);
                return;
            }
        }
    }

    std::optional<Vec3> coordAt(int32_t atomIndex, std::size_t step) const {
        if (!validAtomIndex(atomIndex)) return std::nullopt;
        const Vec3 v = verbs::pos(body_, static_cast<std::size_t>(atomIndex), step);
        if (!v.allFinite()) return std::nullopt;
        return v;
    }

    void foldMolecularFrame(std::size_t step) {
        const MolFrameSpec& spec = molecular_frame_spec_;
        std::optional<Mat3> frame;
        switch (spec.kind) {
        case MolFrameKind::None:
            return;
        case MolFrameKind::BackboneCA: {
            const auto ca = coordAt(spec.origin, step);
            const auto n = coordAt(spec.x_anchor, step);
            const auto c = coordAt(spec.plane_anchor, step);
            if (ca && n && c) {
                const LocalFrame caFrame = BuildBackboneCaFrame(*ca, *n, *c);
                if (caFrame.is_valid) {
                    Mat3 axes;
                    axes.col(0) = caFrame.x;
                    axes.col(1) = caFrame.y;
                    axes.col(2) = caFrame.z;
                    frame = axes;
                }
            }
            break;
        }
        case MolFrameKind::BackboneCarbonyl:
        case MolFrameKind::BackboneAmideN: {
            const auto origin = coordAt(spec.origin, step);
            const auto xAnchor = coordAt(spec.x_anchor, step);
            const auto planeAnchor = coordAt(spec.plane_anchor, step);
            if (origin && xAnchor && planeAnchor)
                frame = frameFromXAndPlane(*xAnchor - *origin, *planeAnchor - *origin);
            break;
        }
        case MolFrameKind::BackboneAmideH: {
            const auto n = coordAt(spec.origin, step);
            const auto h = coordAt(spec.x_anchor, step);
            const auto ca = coordAt(spec.second_anchor, step);
            const bool hasPrevC = spec.plane_anchor >= 0;
            const auto cPrev = hasPrevC ? coordAt(spec.plane_anchor, step) : std::optional<Vec3>{};
            if (n && h && ca && (!hasPrevC || cPrev)) {
                const LocalFrame hn = BuildHNFrame(*n, *h, *ca,
                                                   hasPrevC ? *cPrev : Vec3::Zero(),
                                                   hasPrevC);
                if (hn.is_valid) {
                    Mat3 axes;
                    axes.col(0) = hn.x;
                    axes.col(1) = hn.y;
                    axes.col(2) = hn.z;
                    frame = axes;
                }
            }
            break;
        }
        case MolFrameKind::MetSd: {
            const auto origin = coordAt(spec.origin, step);
            const auto xAnchor = coordAt(spec.x_anchor, step);
            const auto planeAnchor = coordAt(spec.plane_anchor, step);
            if (origin && xAnchor && planeAnchor)
                frame = frameFromXAndPlaneLocked(*xAnchor - *origin, *planeAnchor - *origin,
                                                 prev_mol_locked_x_, prev_mol_locked_z_);
            break;
        }
        case MolFrameKind::AromaticRingLocal: {
            if (spec.ring < 0 || spec.heavy < 0) break;
            const auto heavy = coordAt(spec.heavy, step);
            if (!heavy) break;
            const model::RingGeometry& g = body_.idx.ringGeometry.at(static_cast<std::size_t>(spec.ring), step);
            frame = frameFromXAndZ(*heavy - g.center, g.normal);
            break;
        }
        case MolFrameKind::SidechainCarboxylate: {
            const auto c = coordAt(spec.origin, step);
            const auto o1 = coordAt(spec.x_anchor, step);
            const auto o2 = coordAt(spec.second_anchor, step);
            if (c && o1 && o2) {
                const auto d1 = normalizeFrameVec(*o1 - *c);
                const auto d2 = normalizeFrameVec(*o2 - *c);
                if (d1 && d2)
                    frame = frameFromXAndPlaneLocked(*d1 + *d2, *o1 - *c,
                                                     prev_mol_locked_x_, prev_mol_locked_z_);
            }
            break;
        }
        case MolFrameKind::SidechainGuanidinium:
        case MolFrameKind::SidechainCarboxamide:
        case MolFrameKind::AliphaticCarbon:
        case MolFrameKind::HydroxylOxygen: {
            const auto origin = coordAt(spec.origin, step);
            const auto xAnchor = coordAt(spec.x_anchor, step);
            const auto planeAnchor = coordAt(spec.plane_anchor, step);
            if (origin && xAnchor && planeAnchor)
                frame = frameFromXAndPlaneLocked(*xAnchor - *origin, *planeAnchor - *origin,
                                                 prev_mol_locked_x_, prev_mol_locked_z_);
            break;
        }
        }
        if (frame) molecular_frame_.set(step, *frame);
    }

    std::optional<Mat3> molecularAxesAt(std::size_t step) const {
        if (step >= molecular_frame_.values.size()) return std::nullopt;
        const MolecularFrameValue& f = molecular_frame_.values[step];
        if (!f.valid || !f.axes.allFinite()) return std::nullopt;
        return f.axes;
    }

    QJsonObject molecularFrameAuditJson() const {
        QJsonObject audit;
        audit.insert(QStringLiteral("frame_kind"), molecular_frame_.kind);
        audit.insert(QStringLiteral("frame_kind_ord"), molecular_frame_.kind_ord);
        audit.insert(QStringLiteral("frame_variant"), molecular_frame_.variant);
        audit.insert(QStringLiteral("frame_variant_ord"), molecular_frame_.variant_ord);
        audit.insert(QStringLiteral("checked"), molecular_frame_.kind_ord != 0);
        int validFrames = 0;
        int normalFlips = 0;
        double maxOrthonormalAbs = 0.0;
        double maxDetAbs = 0.0;
        double maxSigmaRoundtripAbs = 0.0;
        std::optional<Vec3> prevZ;
        for (std::size_t step = 0; step < molecular_frame_.values.size(); ++step) {
            const MolecularFrameValue& f = molecular_frame_.values[step];
            if (!f.valid) continue;
            ++validFrames;
            const Mat3 gram = f.axes.transpose() * f.axes;
            const Mat3 err = gram - Mat3::Identity();
            maxOrthonormalAbs = std::max(maxOrthonormalAbs, err.cwiseAbs().maxCoeff());
            maxDetAbs = std::max(maxDetAbs, std::abs(f.axes.determinant() - 1.0));
            const Vec3 z = f.axes.col(2);
            if (prevZ && prevZ->dot(z) < 0.0) ++normalFlips;
            prevZ = z;
            if (step < sigma_total_raw_.values.size() && sigma_total_raw_.present[step]) {
                maxSigmaRoundtripAbs =
                    std::max(maxSigmaRoundtripAbs,
                             sigmaRoundTripMaxAbs(f.axes, sigma_total_raw_.values[step]));
            }
        }
        audit.insert(QStringLiteral("valid_frames"), validFrames);
        audit.insert(QStringLiteral("max_orthonormal_abs"), jd(validFrames > 0 ? maxOrthonormalAbs : kNaN));
        audit.insert(QStringLiteral("max_det_abs"), jd(validFrames > 0 ? maxDetAbs : kNaN));
        audit.insert(QStringLiteral("normal_flip_count"), normalFlips);
        audit.insert(QStringLiteral("max_sigma_roundtrip_abs"),
                     jd(validFrames > 0 ? maxSigmaRoundtripAbs : kNaN));
        audit.insert(QStringLiteral("roundtrip_tolerance"), 1.0e-13);
        audit.insert(QStringLiteral("roundtrip_method"),
                     QStringLiteral("long_double_reorthonormalized_axes"));
        return audit;
    }

    void setScalarIfPresent(ScalarSeries& out, ArrayId id, std::size_t step) {
        if (!body_.catalog.present(body_, id, atom_, step)) return;
        const double v = body_.catalog.value(body_, id, atom_, step);
        if (finite(v)) out.set(step, v);
    }

    void setVec3IfPresent(Vec3Series& out, ArrayId id, std::size_t step) {
        if (!body_.catalog.present(body_, id, atom_, step)) return;
        const Vec3 v = body_.catalog.valueVec3(body_, id, atom_, step);
        if (v.allFinite()) out.set(step, v);
    }

    void foldGeometry(std::size_t step) {
        for (const auto& item : static_bonds_) {
            const std::uint64_t key = item.first;
            const StaticBondInfo& info = item.second;
            const auto a = coordAt(static_cast<int32_t>(atom_), step);
            const auto b = coordAt(static_cast<int32_t>(info.other), step);
            if (!a || !b) continue;
            const double length = (*b - *a).norm();
            if (finite(length) && length > 1e-12) {
                bond_length_series_[key].set(step, length);
                bond_length_rejected_degenerate_[key].set(step, false);
            } else {
                bond_length_rejected_degenerate_[key].set(step, true);
            }
        }

        setScalarIfPresent(sasa_, ArrayId::Sasa, step);
        setVec3IfPresent(hbond_nearest_dir_, ArrayId::HbondNearestDirection, step);
        if (step < hbond_nearest_dir_.present.size() && hbond_nearest_dir_.present[step]) {
            if (const auto axes = molecularAxesAt(step)) {
                const Vec3 local = axes->transpose() * hbond_nearest_dir_.values[step];
                if (local.allFinite()) hbond_nearest_dir_mol_.set(step, local);
            }
        }
        setScalarIfPresent(hbond_count_, ArrayId::HbondCount, step);
        setScalarIfPresent(ff_pb_radius_, ArrayId::FfPbRadius, step);
    }

    void foldMopacSelfState(std::size_t step) {
        double pi = 0.0;
        bool anyBond = false;
        for (const MopacFrameBond& b : context_.mopacBondsForAtom(step, atom_)) {
            pi += std::max(0.0, b.wiberg_order - 1.0);
            anyBond = true;
        }
        if (anyBond) {
            pi_character_.set(step, pi);
            model::Hybridisation h = model::Hybridisation::Unassigned;
            if (pi < 0.4)
                h = model::Hybridisation::sp3;
            else if (pi < 1.5)
                h = model::Hybridisation::sp2;
            else
                h = model::Hybridisation::sp;
            hybridisation_.set(step, static_cast<int>(h));
            updateOxygenGate(h);
        }

        auto read = [&](int comp) -> std::optional<double> {
            return body_.catalog.value(body_, io::FieldKind::MOPACScalars, atom_, step, comp);
        };
        const std::optional<double> charge = read(0);
        const std::optional<double> s = read(1);
        const std::optional<double> p = read(2);
        const std::optional<double> val = read(3);
        if (charge) mopac_charge_.set(step, *charge);
        if (s) mopac_s_pop_.set(step, *s);
        if (p) mopac_p_pop_.set(step, *p);
        if (val) mopac_valency_.set(step, *val);
        if (s && p && finite(*s) && finite(*p) && (*s + *p) > 1e-12)
            mopac_s_character_.set(step, *s / (*s + *p));
    }

    void updateOxygenGate(model::Hybridisation h) {
        (void)h;
        if (!body_.run.protein) return;
        const model::QtAtom& a = body_.run.protein->atom(atom_);
        if (a.element != model::Element::O) return;
        oxygen_gate_checked_ = true;
    }

    void foldTopologyJoin(std::size_t step) {
        for (const MopacFrameBond& b : context_.mopacBondsForAtom(step, atom_)) {
            const std::uint64_t key = pairKey(b.atom_a, b.atom_b);
            auto it = bond_series_.find(key);
            if (it == bond_series_.end()) it = bond_series_.emplace(key, ScalarSeries(cadence_.stepCount())).first;
            it->second.set(step, b.wiberg_order);
        }
    }

    void foldDihedrals(std::size_t step) {
        static constexpr std::array<DihedralKind, 7> kinds = {
            DihedralKind::Phi, DihedralKind::Psi, DihedralKind::Omega, DihedralKind::Chi1,
            DihedralKind::Chi2, DihedralKind::Chi3, DihedralKind::Chi4,
        };
        std::array<ScalarSeries*, 7> series = {&phi_, &psi_, &omega_, &chi1_, &chi2_, &chi3_, &chi4_};
        for (std::size_t i = 0; i < kinds.size(); ++i) {
            const DihedralState st = body_.idx.dihedrals.state(kinds[i], atom_, step);
            if (st.present && finite(st.radians)) series[i]->set(step, st.radians);
        }
    }

    void setFieldDerived(std::size_t step,
                         const Vec3& field,
                         Vec3Series& mol,
                         ScalarSeries& ez,
                         ScalarSeries& e2,
                         ScalarSeries& absE,
                         ScalarSeries& absE2) {
        if (!field.allFinite()) return;
        const double e2Value = field.squaredNorm();
        if (!finite(e2Value)) return;
        absE.set(step, std::sqrt(e2Value));
        absE2.set(step, e2Value);

        const auto axes = molecularAxesAt(step);
        if (!axes) return;
        const Vec3 local = axes->transpose() * field;
        if (!local.allFinite()) return;
        mol.set(step, local);
        ez.set(step, local.z());
        e2.set(step, e2Value);
    }

    void foldDirectChargeField(std::size_t step) {
        if (!body_.run.protein || !body_.catalog.has(ArrayId::Ff14sbCharge)) return;
        const auto target = coordAt(static_cast<int32_t>(atom_), step);
        if (!target) return;

        Vec3 field = Vec3::Zero();
        int sourceCount = 0;
        int rejectedSelf = 0;
        int rejectedDegenerate = 0;
        int rejectedCutoff = 0;
        int missingCharge = 0;
        const std::size_t nAtoms = body_.run.protein->atomCount();
        for (std::size_t source = 0; source < nAtoms; ++source) {
            if (source == atom_) {
                ++rejectedSelf;
                continue;
            }
            if (!body_.catalog.present(body_, ArrayId::Ff14sbCharge, source, step)) {
                ++missingCharge;
                continue;
            }
            const auto sourcePos = coordAt(static_cast<int32_t>(source), step);
            if (!sourcePos) {
                ++rejectedDegenerate;
                continue;
            }
            const Vec3 disp = *sourcePos - *target;
            const double r = disp.norm();
            if (!finite(r) || !(r > 1e-12)) {
                ++rejectedDegenerate;
                continue;
            }
            if (r > config_.charge_cutoff_A) {
                ++rejectedCutoff;
                continue;
            }
            const double q = body_.catalog.value(body_, ArrayId::Ff14sbCharge, source, step);
            if (!finite(q)) {
                ++missingCharge;
                continue;
            }
            field += (-CoulombKeVA() * q / (r * r * r)) * disp;
            ++sourceCount;
        }

        charge_field_source_count_.set(step, sourceCount);
        charge_field_rejected_self_.set(step, rejectedSelf);
        charge_field_rejected_degenerate_.set(step, rejectedDegenerate);
        charge_field_rejected_cutoff_.set(step, rejectedCutoff);
        charge_field_missing_charge_.set(step, missingCharge);
        if (field.allFinite()) field_charge_ff14sb_.set(step, field);
    }

    void foldFieldsAndRelationships(std::size_t step) {
        Vec3 ff = Vec3::Zero();
        bool ffPresent = false;
        const auto pairs = PerAtomRowPairContributions(body_, atom_, step, config_, LocalFrame{});
        for (const PairContribution& pair : pairs) {
            foldRelationship(step, pair);
            if (pair.mechanism == QStringLiteral("charge_q_over_r3")
                && pair.source_atom_index >= 0 && pair.r > 1e-12
                && static_cast<std::size_t>(pair.source_atom_index) < body_.run.protein->atomCount()
                && body_.catalog.present(body_, ArrayId::Ff14sbCharge,
                                         static_cast<std::size_t>(pair.source_atom_index), step)) {
                const double q = body_.catalog.value(body_, ArrayId::Ff14sbCharge,
                                                     static_cast<std::size_t>(pair.source_atom_index), step);
                if (finite(q)) {
                    ff += (-q / (pair.r * pair.r * pair.r)) * pair.disp;
                    ffPresent = true;
                }
            }
        }
        if (ffPresent) field_ff14sb_.set(step, ff);
        if (body_.catalog.present(body_, ArrayId::MopacCoulombEfield, atom_, step)) {
            const Vec3 v = body_.catalog.valueVec3(body_, ArrayId::MopacCoulombEfield, atom_, step);
            if (v.allFinite()) field_mopac_.set(step, v);
        }
        if (body_.catalog.present(body_, ArrayId::ApbsEfield, atom_, step)) {
            const Vec3 v = body_.catalog.valueVec3(body_, ArrayId::ApbsEfield, atom_, step);
            if (v.allFinite()) field_apbs_.set(step, v);
        }
        foldDirectChargeField(step);

        if (field_mopac_.present[step])
            setFieldDerived(step, field_mopac_.values[step], field_mol_mopac_, field_E_z_mopac_,
                            field_E2_mopac_, field_abs_E_mopac_, field_abs_E2_mopac_);
        if (field_apbs_.present[step])
            setFieldDerived(step, field_apbs_.values[step], field_mol_apbs_, field_E_z_apbs_,
                            field_E2_apbs_, field_abs_E_apbs_, field_abs_E2_apbs_);
        if (field_charge_ff14sb_.present[step])
            setFieldDerived(step, field_charge_ff14sb_.values[step], field_mol_charge_ff14sb_,
                            field_E_z_charge_ff14sb_, field_E2_charge_ff14sb_,
                            field_abs_E_charge_ff14sb_, field_abs_E2_charge_ff14sb_);
    }

    template <std::size_t N>
    bool setFixedArrayIfPresent(FixedArraySeries<N>& out, ArrayId id, std::size_t step) {
        if (!body_.catalog.has(id)) return false;
        const ArraySpec& spec = body_.catalog.spec(id);
        if (spec.axes.comp_count < static_cast<int>(N)) return false;
        if (!body_.catalog.present(body_, id, atom_, step)) return false;
        std::array<double, N> values{};
        bool ok = true;
        for (std::size_t c = 0; c < N; ++c) {
            const double v = body_.catalog.value(body_, id, atom_, step, -1, static_cast<int>(c));
            values[c] = v;
            ok = ok && finite(v);
        }
        if (!ok) return false;
        out.set(step, values);
        return true;
    }

    bool setTensorIfPresent(TensorSeries& out, ArrayId id, std::size_t step) {
        if (!body_.catalog.has(id)) return false;
        if (body_.catalog.spec(id).axes.comp_count < 9) return false;
        if (!body_.catalog.present(body_, id, atom_, step)) return false;
        const model::SphericalTensor t = body_.catalog.valueTensor(body_, id, atom_, step);
        if (!tensorFinite(t)) return false;
        out.set(step, t);
        return true;
    }

    void foldProducerSubstrate(std::size_t step) {
        if (setFixedArrayIfPresent(bs_per_type_T2_, ArrayId::BSPerTypeT2, step))
            projectFixedT2Array(step, bs_per_type_T2_.values[step], bs_per_type_mol_);
        if (setFixedArrayIfPresent(hm_per_type_T2_, ArrayId::HMPerTypeT2, step))
            projectFixedT2Array(step, hm_per_type_T2_.values[step], hm_per_type_mol_);
        for (std::size_t i = 0; i < kMcTensorFields.size(); ++i) {
            if (setTensorIfPresent(mc_tensor_series_[i], kMcTensorFields[i].id, step)) {
                if (const auto local =
                        projectSphericalTensorToMolecular(step, mc_tensor_series_[i].values[step]))
                    mc_tensor_mol_series_[i].set(step, *local);
            }
        }
        if (setTensorIfPresent(tripeptide_bb_shielding_, ArrayId::TripeptideBBShielding, step)) {
            if (const auto local =
                    projectSphericalTensorToMolecular(step, tripeptide_bb_shielding_.values[step]))
                tripeptide_bb_shielding_mol_.set(step, *local);
        }
        if (setTensorIfPresent(larsen_hbond_shielding_, ArrayId::LarsenHBondShielding, step)) {
            if (const auto local =
                    projectSphericalTensorToMolecular(step, larsen_hbond_shielding_.values[step]))
                larsen_hbond_shielding_mol_.set(step, *local);
        }
    }

    void foldRelationship(std::size_t step, const PairContribution& pair) {
        const int hyb = hybridisation_.values[step] == IntSeries::kMissing
                            ? static_cast<int>(model::Hybridisation::Unassigned)
                            : hybridisation_.values[step];
        const model::QtAtom& atom = body_.run.protein->atom(atom_);
        const auto cr = contactedResidue(body_, pair);
        if (cr.first >= 0) contacted_residues_.insert({cr.first, static_cast<int>(cr.second)});
        const bool present = (pair.pointer_flags & PresentFlag) != 0;
        const QString scope = relationshipScope(body_, atom_, pair);
        RelationshipKey key;
        key.hybridisation_ord = hyb;
        key.element_ord = static_cast<int>(atom.element);
        key.contacted_residue_number = cr.first;
        key.contacted_amino_acid_ord = static_cast<int>(cr.second);
        key.scope_ord = scopeOrd(scope);
        key.mechanism_ord = mechanismOrd(pair.mechanism);
        key.source_kind_ord = sourceKindOrd(pair.source_kind);
        key.source_id = pair.source_id;
        key.source_atom_index = pair.source_atom_index;
        key.source_category_ord = pair.source_category_ord;
        auto it = relationships_.find(key);
        if (it == relationships_.end())
            it = relationships_.emplace(key, RelationshipSeries(cadence_.stepCount())).first;
        RelationshipSeries& r = it->second;
        if (step >= cadence_.stepCount()) return;
        foldContactClass(step, cr, present, pair.r);
        r.markPresent(step, present);
        r.distance.set(step, pair.r);
        r.cos_theta.set(step, pair.cos_theta);
        r.inv_r3.set(step, pair.inv_r3);
        r.dipolar.set(step, pair.dipolar);
        r.kernel_T0.set(step, pair.kernel_T0);
        r.contribution.set(step, pair.contribution);
        r.kernel_T2.set(step, pair.kernel_T2);
        if (const auto local = projectLibraryT2ToMolecular(step, 0.0, pair.kernel_T2))
            r.kernel_mol.set(step, *local);
    }

    static bool efgT2Complete(const EfgValue& v) {
        for (double x : v.t2)
            if (!finite(x)) return false;
        return true;
    }

    void setEfgDerived(std::size_t step, const EfgValue& v, Mat3Series& mol, ScalarSeries& absSeries) {
        if (!efgT2Complete(v)) return;
        absSeries.set(step, norm5(v.t2));
        const auto axes = molecularAxesAt(step);
        if (!axes) return;
        const Mat3 lab = ReconstructLibraryT2Matrix(finite(v.t0) ? v.t0 : 0.0, v.t2);
        const Mat3 local = axes->transpose() * lab * (*axes);
        if (local.allFinite()) mol.set(step, local);
    }

    void foldEfg(std::size_t step) {
        const EfgValue apbs = readT2Field(body_, io::FieldKind::APBSEFG, ArrayId::ApbsEfg, atom_, step);
        const EfgValue aimnet = readT2Field(body_, io::FieldKind::AIMNet2EFG, ArrayId::Aimnet2Efg, atom_, step);
        efg_apbs_.set(step, apbs);
        efg_aimnet2_.set(step, aimnet);
        setEfgDerived(step, apbs, efg_mol_apbs_, efg_abs_apbs_);
        setEfgDerived(step, aimnet, efg_mol_aimnet2_, efg_abs_aimnet2_);

        const EfgValue shielding = readT2Field(body_, io::FieldKind::MOPACCoulombEFG,
                                              ArrayId::MopacCoulombShielding, atom_, step);
        shielding_mopac_coulomb_.set(step, shielding);
        // §4.7.3: |T2| + molecular-frame projection (it is an EFG, projected
        // exactly as foldEfg does for efg_mol). setEfgDerived sets the abs norm
        // and, on framed atoms, the molecular matrix.
        setEfgDerived(step, shielding, shielding_mol_mopac_, shielding_abs_mopac_coulomb_);
    }

    void foldSigma(std::size_t step) {
        const std::size_t original = cadence_.originalIndex(step);
        const DftTarget target = BuildTarget(body_.run, atom_, original, LocalFrame{});
        if (!target.present || !finite(target.total_decomp.T0)) return;
        sigma_total_.set(step, target.total_decomp);
        sigma_dia_.set(step, target.dia_decomp);
        sigma_para_.set(step, target.para_decomp);
        sigma_total_raw_.set(step, target.total_raw);
        sigma_dia_raw_.set(step, target.dia_raw);
        sigma_para_raw_.set(step, target.para_raw);
        if (const auto axes = molecularAxesAt(step)) {
            sigma_mol_total_.set(step, axes->transpose() * target.total_raw * (*axes));
            sigma_mol_dia_.set(step, axes->transpose() * target.dia_raw * (*axes));
            sigma_mol_para_.set(step, axes->transpose() * target.para_raw * (*axes));
        }
        sigma_frob_.set(step, tensorFrobenius(target.total_raw));
        foldCsa(step, target.total_raw);
        ++sigma_folds_;
    }

    void foldCsa(std::size_t step, const Mat3& totalRaw) {
        const Mat3 sym = 0.5 * (totalRaw + totalRaw.transpose());
        Eigen::SelfAdjointEigenSolver<Mat3> solver(sym);
        if (solver.info() != Eigen::Success) return;

        CsaDescriptor d;
        d.principal_values = solver.eigenvalues();
        d.pas_axes = solver.eigenvectors();
        if (!d.principal_values.allFinite() || !d.pas_axes.allFinite()) return;

        if (prev_pas_axes_) {
            for (int c = 0; c < 3; ++c) {
                if (prev_pas_axes_->col(c).dot(d.pas_axes.col(c)) < 0.0)
                    d.pas_axes.col(c) *= -1.0;
            }
        }
        if (d.pas_axes.determinant() < 0.0) {
            if (prev_pas_axes_) {
                int flipCol = 0;
                double smallest = std::numeric_limits<double>::infinity();
                for (int c = 0; c < 3; ++c) {
                    const double continuity = std::abs(prev_pas_axes_->col(c).dot(d.pas_axes.col(c)));
                    if (continuity < smallest) {
                        smallest = continuity;
                        flipCol = c;
                    }
                }
                d.pas_axes.col(flipCol) *= -1.0;
            } else {
                d.pas_axes.col(2) *= -1.0;
            }
        }

        const double s11 = d.principal_values[0];
        const double s22 = d.principal_values[1];
        const double s33 = d.principal_values[2];
        d.sigma_iso = (s11 + s22 + s33) / 3.0;
        d.haeberlen_span = s33 - s11;
        if (!(std::abs(d.haeberlen_span) > 1e-12)) return;

        int zz = 0;
        double farthest = -1.0;
        for (int i = 0; i < 3; ++i) {
            const double dist = std::abs(d.principal_values[i] - d.sigma_iso);
            if (dist > farthest) {
                farthest = dist;
                zz = i;
            }
        }
        std::array<int, 2> rem{};
        int k = 0;
        for (int i = 0; i < 3; ++i)
            if (i != zz) rem[static_cast<std::size_t>(k++)] = i;

        int xx = rem[0];
        int yy = rem[1];
        const double dx = std::abs(d.principal_values[xx] - d.sigma_iso);
        const double dy = std::abs(d.principal_values[yy] - d.sigma_iso);
        if (dx < dy) std::swap(xx, yy);  // xx is the remaining value farther from isotropic.

        const double denomEta = d.principal_values[zz] - d.sigma_iso;
        if (!(std::abs(denomEta) > 1e-12)) return;
        d.haeberlen_eta = (d.principal_values[yy] - d.principal_values[xx]) / denomEta;
        if (d.haeberlen_eta < 0.0 && d.haeberlen_eta > -1e-12) d.haeberlen_eta = 0.0;
        if (d.haeberlen_eta > 1.0 && d.haeberlen_eta < 1.0 + 1e-12) d.haeberlen_eta = 1.0;
        if (!(d.haeberlen_eta >= 0.0 && d.haeberlen_eta <= 1.0)) return;

        d.haeberlen_skew = 3.0 * (s22 - d.sigma_iso) / d.haeberlen_span;
        d.haeberlen_values = Vec3(d.principal_values[xx], d.principal_values[yy], d.principal_values[zz]);
        d.haeberlen_axes.col(0) = d.pas_axes.col(xx);
        d.haeberlen_axes.col(1) = d.pas_axes.col(yy);
        d.haeberlen_axes.col(2) = d.pas_axes.col(zz);
        d.valid = finite(d.sigma_iso) && finite(d.haeberlen_eta) && finite(d.haeberlen_span)
                  && finite(d.haeberlen_skew) && d.haeberlen_axes.allFinite();
        if (!d.valid) return;
        sigma_csa_.set(step, d);
        prev_pas_axes_ = d.pas_axes;
    }

    std::pair<std::size_t, std::size_t> pairFromKey(std::uint64_t k) const {
        return {static_cast<std::size_t>(k >> 32u), static_cast<std::size_t>(k & 0xffffffffu)};
    }

    const AnalysisObjectContext& context_;
    const Body& body_;
    const AnalysisCadence& cadence_;
    std::size_t atom_ = 0;
    PerAtomSubstrateConfig config_;

    Vec3Series coord_;
    MolecularFrameSeries molecular_frame_;
    TensorSeries sigma_total_;
    TensorSeries sigma_dia_;
    TensorSeries sigma_para_;
    Mat3Series sigma_total_raw_;
    Mat3Series sigma_dia_raw_;
    Mat3Series sigma_para_raw_;
    Mat3Series sigma_mol_total_;
    Mat3Series sigma_mol_dia_;
    Mat3Series sigma_mol_para_;
    ScalarSeries sigma_frob_;
    CsaDescriptorSeries sigma_csa_;
    Vec3Series field_mopac_;
    Vec3Series field_ff14sb_;
    Vec3Series field_apbs_;
    Vec3Series field_charge_ff14sb_;
    Vec3Series field_mol_mopac_;
    Vec3Series field_mol_apbs_;
    Vec3Series field_mol_charge_ff14sb_;
    ScalarSeries field_E_z_mopac_;
    ScalarSeries field_E_z_apbs_;
    ScalarSeries field_E_z_charge_ff14sb_;
    ScalarSeries field_E2_mopac_;
    ScalarSeries field_E2_apbs_;
    ScalarSeries field_E2_charge_ff14sb_;
    ScalarSeries field_abs_E_mopac_;
    ScalarSeries field_abs_E_apbs_;
    ScalarSeries field_abs_E_charge_ff14sb_;
    ScalarSeries field_abs_E2_mopac_;
    ScalarSeries field_abs_E2_apbs_;
    ScalarSeries field_abs_E2_charge_ff14sb_;
    IntSeries charge_field_source_count_;
    IntSeries charge_field_rejected_self_;
    IntSeries charge_field_rejected_degenerate_;
    IntSeries charge_field_rejected_cutoff_;
    IntSeries charge_field_missing_charge_;
    EfgSeries efg_apbs_;
    EfgSeries efg_aimnet2_;
    Mat3Series efg_mol_apbs_;
    Mat3Series efg_mol_aimnet2_;
    ScalarSeries efg_abs_apbs_;
    ScalarSeries efg_abs_aimnet2_;
    EfgSeries shielding_mopac_coulomb_;
    ScalarSeries shielding_abs_mopac_coulomb_;
    Mat3Series shielding_mol_mopac_;  // §4.7.3 molecular projection of the MOPAC-Coulomb EFG
    FixedArraySeries<40> bs_per_type_T2_;
    FixedArraySeries<40> hm_per_type_T2_;
    std::vector<TensorSeries> mc_tensor_series_;
    std::vector<Mat3Series> bs_per_type_mol_;
    std::vector<Mat3Series> hm_per_type_mol_;
    std::vector<Mat3Series> mc_tensor_mol_series_;
    TensorSeries tripeptide_bb_shielding_;
    TensorSeries larsen_hbond_shielding_;
    Mat3Series tripeptide_bb_shielding_mol_;
    Mat3Series larsen_hbond_shielding_mol_;
    ScalarSeries sasa_;
    Vec3Series hbond_nearest_dir_;
    Vec3Series hbond_nearest_dir_mol_;
    ScalarSeries hbond_count_;
    ScalarSeries ff_pb_radius_;
    ScalarSeries phi_;
    ScalarSeries psi_;
    ScalarSeries omega_;
    ScalarSeries chi1_;
    ScalarSeries chi2_;
    ScalarSeries chi3_;
    ScalarSeries chi4_;
    ScalarSeries mopac_charge_;
    ScalarSeries mopac_s_pop_;
    ScalarSeries mopac_p_pop_;
    ScalarSeries mopac_valency_;
    ScalarSeries mopac_s_character_;
    IntSeries hybridisation_;
    std::vector<std::uint8_t> contacted_class_membership_;
    ScalarSeries pi_character_;
    std::map<RelationshipKey, RelationshipSeries> relationships_;
    std::map<std::uint64_t, StaticBondInfo> static_bonds_;
    std::map<std::uint64_t, ScalarSeries> bond_series_;
    std::map<std::uint64_t, ScalarSeries> bond_length_series_;
    std::map<std::uint64_t, BoolTriSeries> bond_length_rejected_degenerate_;
    std::set<std::pair<int, int>> contacted_residues_;
    int static_hybridisation_ord_ = static_cast<int>(model::Hybridisation::Unassigned);
    bool oxygen_gate_checked_ = false;
    std::size_t sigma_folds_ = 0;
    std::optional<Mat3> prev_pas_axes_;
    MolFrameSpec molecular_frame_spec_;
    std::optional<Vec3> prev_mol_locked_x_;
    std::optional<Vec3> prev_mol_locked_z_;
};

AnalysisAtom::AnalysisAtom(const AnalysisObjectContext& context,
                           std::size_t atomIndex,
                           PerAtomSubstrateConfig config)
    : AnalysisElement(context, QStringLiteral("atom"), atomIndex),
      impl_(std::make_unique<Impl>(context, atomIndex, std::move(config))) {}

AnalysisAtom::~AnalysisAtom() = default;

void AnalysisAtom::Calculate(std::size_t step) {
    impl_->calculate(step);
}

QJsonObject AnalysisAtom::Truth() const {
    return impl_->truth();
}

std::size_t AnalysisAtom::sigmaFolds() const { return impl_->sigmaFolds(); }
std::size_t AnalysisAtom::relationshipCount() const { return impl_->relationshipCount(); }
std::size_t AnalysisAtom::mappedBondCount() const { return impl_->mappedBondCount(); }
std::size_t AnalysisAtom::mismatchEventCount() const { return impl_->mismatchEventCount(); }
std::size_t AnalysisAtom::accumulatorResponseCount() const { return impl_->last_accumulator_response_count; }
std::size_t AnalysisAtom::accumulatorContextCount() const { return impl_->last_accumulator_context_count; }
bool AnalysisAtom::oxygenGatePassed() const { return impl_->oxygenGatePassed(); }
std::size_t AnalysisAtom::writeBoundedSigmaRows(QTextStream& out,
                                                const QString& datasetId,
                                                const QString& proteinId) const {
    return impl_->writeBoundedSigmaRows(out, datasetId, proteinId);
}
std::size_t AnalysisAtom::writeClassicalSourceTermRows(QTextStream& out,
                                                       const QString& datasetId,
                                                       const QString& proteinId) const {
    return impl_->writeClassicalSourceTermRows(out, datasetId, proteinId);
}
ClassicalSourceTermRecord AnalysisAtom::classicalSourceTermRecord(const QString& datasetId,
                                                                  const QString& proteinId) const {
    return impl_->classicalSourceTermRecord(datasetId, proteinId);
}
std::size_t AnalysisAtom::writeSourceFamilyMatrixRows(QTextStream& out,
                                                      const QString& datasetId,
                                                      const QString& proteinId) const {
    return impl_->writeSourceFamilyMatrixRows(out, datasetId, proteinId);
}
std::size_t AnalysisAtom::writeSubspaceOverlapRows(QTextStream& out,
                                                   const QString& datasetId,
                                                   const QString& proteinId) const {
    return impl_->writeSubspaceOverlapRows(out, datasetId, proteinId);
}
std::size_t AnalysisAtom::writeEtaByWellRows(QTextStream& out,
                                             const QString& datasetId,
                                             const QString& proteinId) const {
    return impl_->writeEtaByWellRows(out, datasetId, proteinId);
}

void AnalysisAtom::WriteBoundedSigmaHeader(QTextStream& out) {
    out << "schema_version,dataset_id,protein_id,atom_uid,atom_index,residue_index,"
           "residue_number,residue_type,atom_name,element,backbone_role,frame_kind,"
           "frame_kind_ord,frame_variant,frame_variant_ord,frame_valid,"
           "frame_missing_reason,original_frame_index,sigma_ordinal,sigma_step,"
           "sigma_valid,sigma_iso,molcomp_xx,molcomp_yy,molcomp_xy,molcomp_xz,"
           "molcomp_yz,molcomp_zz,invariant_span,invariant_aniso,invariant_eta_H,"
           "invariant_skew,invariant_frobenius,antisymmetric_norm,t1_fraction,"
           "mol_frobenius,sigma_roundtrip_abs,support_class,finite_frac,"
           "singleton_flag,missing_n,missing_reason,hyb,contact_class,"
           "dihedral_region,SS\n";
}

void AnalysisAtom::WriteClassicalSourceTermHeader(QTextStream& out) {
    out << classicalSourceHeaderColumns().join(QLatin1Char(',')) << '\n';
}

void AnalysisAtom::WriteSourceFamilyMatrixHeader(QTextStream& out) {
    out << "schema_version,dataset_id,protein_id,atom_uid,atom_index,residue_index,"
           "residue_number,residue_type,atom_name,element,backbone_role,identity_key,"
           "context_key,axis,target_view,family_key,family_label,channel_count,"
           "finite_row_count,support_class,finite_frac,singleton_flag,missing_n,"
           "missing_reason,channel_names,channel_means,channel_sds,cross_axis_join_status\n";
}

void AnalysisAtom::WriteSubspaceOverlapHeader(QTextStream& out) {
    out << "schema_version,dataset_id,protein_id,atom_uid,atom_index,residue_index,"
           "residue_number,residue_type,atom_name,element,backbone_role,identity_key,"
           "context_key,axis,target_view,read_id,family_a,family_b,finite_n,input_dim_a,"
           "input_dim_b,active_dim_a,active_dim_b,basis_dim_a,basis_dim_b,"
           "explained_fraction_a,explained_fraction_b,condition_number_a,"
           "condition_number_b,max_canonical_corr,mean_canonical_corr,n_cc_ge_0_80,"
           "n_cc_ge_0_95,min_angle_deg,canonical_corrs,principal_angles_deg,"
           "overlap_label,provenance,status,missing_reason,support_class,finite_frac,"
           "singleton_flag,missing_n,missing_reason_count,dropped_channels_a,"
           "dropped_channels_b,independence_verdict,independence_basis,"
           "cross_axis_join_status\n";
}

void AnalysisAtom::WriteEtaByWellHeader(QTextStream& out) {
    out << "schema_version,dataset_id,protein_id,atom_uid,atom_index,residue_index,"
           "residue_number,residue_type,atom_name,element,backbone_role,well_family,"
           "well_source,target_key,target_frame,finite_n,support_class,finite_frac,"
           "singleton_flag,missing_n,missing_reason,eta2,well_counts\n";
}

bool AnalysisAtom::AssertMolCompOrder(QString* errOut) {
    constexpr std::array<const char*, 6> expected = {"xx", "yy", "xy", "xz", "yz", "zz"};
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (QString::fromLatin1(kMolCompNames[i]) != QString::fromLatin1(expected[i])
            || static_cast<int>(static_cast<MolComp>(i)) != static_cast<int>(i)) {
            if (errOut) {
                *errOut = QStringLiteral("MolComp order mismatch; expected xx,yy,xy,xz,yz,zz");
            }
            return false;
        }
    }
    Mat3 probe;
    probe << 1.0, 3.0, 4.0,
             3.0, 2.0, 5.0,
             4.0, 5.0, 6.0;
    const std::array<double, 6> c = symmetricComponents(probe);
    const std::array<double, 6> values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (std::abs(c[i] - values[i]) > 1e-12) {
            if (errOut) {
                *errOut = QStringLiteral("MolComp symmetricComponents order mismatch; expected xx,yy,xy,xz,yz,zz");
            }
            return false;
        }
    }
    return true;
}

bool AnalysisAtom::AssertPasShapeConvention(QString* errOut) {
    return AssertPasShapeConventionEnv(errOut);
}

class AnalysisRing::Impl {
public:
    Impl(const AnalysisObjectContext& context, std::size_t ring, double cutoff)
        : context_(context),
          body_(context.body()),
          cadence_(context.cadence()),
          ring_(ring),
          cutoff_(cutoff),
          center_(cadence_.stepCount()),
          normal_(cadence_.stepCount()),
          radius_(cadence_.stepCount()),
          flip_(cadence_.stepCount()),
          pucker_q_(cadence_.stepCount()),
          pucker_theta_(cadence_.stepCount()),
          near_center_(cadence_.stepCount()) {}

    void calculate(std::size_t step) {
        if (!body_.run.protein || ring_ >= body_.run.protein->topology().ringCount()) return;
        const model::QtRing& ring = body_.run.protein->topology().ringAt(ring_);
        const model::RingGeometry& g = body_.idx.ringGeometry.at(ring_, step);
        center_.set(step, g.center);
        normal_.set(step, g.normal);
        radius_.set(step, g.radius);
        if (ring.IsAromatic()) {
            const std::size_t native = ring.nativeAxisIndex >= 0 ? static_cast<std::size_t>(ring.nativeAxisIndex) : ring_;
            const std::optional<double> v =
                body_.catalog.value(body_, io::FieldKind::AromaticChi2, native, step, 0);
            if (v) flip_.set(step, *v);
        } else {
            const std::size_t native = ring.nativeAxisIndex >= 0 ? static_cast<std::size_t>(ring.nativeAxisIndex) : 0;
            if (const std::optional<double> q =
                    body_.catalog.value(body_, io::FieldKind::PuckerQ, native, step, 0))
                pucker_q_.set(step, *q);
            if (const std::optional<double> th =
                    body_.catalog.value(body_, io::FieldKind::PuckerTheta, native, step, 0))
                pucker_theta_.set(step, *th);
        }
        for (const SourceRef& ref : verbs::nearPoint(body_, CloudKind::Atoms, g.center, step, cutoff_)) {
            if (ref.entity_index < 0 || ref.cloud_index < 0) continue;
            const Vec3 pos = body_.idx.spatial.tree(CloudKind::Atoms, step)
                                 .pointAt(static_cast<std::size_t>(ref.cloud_index));
            const Vec3 d = pos - g.center;
            const double dist = d.norm();
            const double nNorm = g.normal.norm();
            const double axial = nNorm > 1e-12 ? d.dot(g.normal / nNorm) : kNaN;
            const double radial = finite(axial) ? std::sqrt(std::max(0.0, dist * dist - axial * axial)) : kNaN;
            near_center_[step].push_back({ref.entity_index, dist, axial, radial});
        }
    }

    QJsonObject truth() const {
        QJsonObject root;
        root.insert(QStringLiteral("object_type"), QStringLiteral("ring"));
        root.insert(QStringLiteral("identity"), identityJson());
        root.insert(QStringLiteral("model"), modelJson());
        root.insert(QStringLiteral("series"), seriesJson());
        root.insert(QStringLiteral("boost"), buildBoostJson());
        return root;
    }

    std::size_t nearHits() const {
        std::size_t n = 0;
        for (const auto& step : near_center_) n += step.size();
        return n;
    }

    mutable std::size_t last_boost_serial_count = 0;

private:
    struct Hit {
        int32_t atom = -1;
        double distance = kNaN;
        double axial = kNaN;
        double radial = kNaN;
    };

    QJsonObject identityJson() const {
        const model::QtRing& r = body_.run.protein->topology().ringAt(ring_);
        QJsonObject o;
        o.insert(QStringLiteral("uid"), uid(QStringLiteral("ring"), ring_));
        o.insert(QStringLiteral("object_type"), QStringLiteral("ring"));
        o.insert(QStringLiteral("model_index"), static_cast<int>(ring_));
        o.insert(QStringLiteral("ring_id"), static_cast<int>(ring_));
        o.insert(QStringLiteral("parent_residue_index"), r.parentResidueIndex);
        o.insert(QStringLiteral("parent_residue_number"), r.parentResidueNumber);
        o.insert(QStringLiteral("ring_kind_ord"), static_cast<int>(r.ringKind));
        o.insert(QStringLiteral("ring_kind"), r.ringKind == model::RingKind::Aromatic ? QStringLiteral("Aromatic")
                                                                                       : QStringLiteral("Saturated"));
        o.insert(QStringLiteral("ring_type_ord"), static_cast<int>(r.TypeIndex()));
        o.insert(QStringLiteral("ring_type"), ringTypeName(r.TypeIndex()));
        o.insert(QStringLiteral("native_axis_index"), r.nativeAxisIndex);
        o.insert(QStringLiteral("fused_partner_ring_id"), r.fusedPartnerRingId);
        QJsonArray atoms;
        for (int32_t a : r.atomIndices) atoms.append(a);
        o.insert(QStringLiteral("atom_indices"), atoms);
        QJsonObject phys;
        phys.insert(QStringLiteral("literature_intensity"), r.LiteratureIntensity());
        phys.insert(QStringLiteral("jb_lobe_offset"), r.JohnsonBoveyLobeOffset());
        phys.insert(QStringLiteral("nitrogen_count"), r.NitrogenCount());
        phys.insert(QStringLiteral("aromaticity_ord"), static_cast<int>(r.Aromaticity()));
        phys.insert(QStringLiteral("aromaticity"), aromaticityName(r.Aromaticity()));
        phys.insert(QStringLiteral("ring_size"), r.RingSizeValue());
        o.insert(QStringLiteral("physics"), phys);
        return o;
    }

    QJsonObject modelJson() const {
        const model::QtRing& r = body_.run.protein->topology().ringAt(ring_);
        QJsonObject m;
        QJsonObject typing;
        typing.insert(QStringLiteral("ring_type"), ringTypeName(r.TypeIndex()));
        typing.insert(QStringLiteral("ring_type_ord"), static_cast<int>(r.TypeIndex()));
        typing.insert(QStringLiteral("ring_kind"), r.ringKind == model::RingKind::Aromatic ? QStringLiteral("Aromatic")
                                                                                           : QStringLiteral("Saturated"));
        m.insert(QStringLiteral("typing"), typing);
        QJsonObject category;
        category.insert(QStringLiteral("ring_type"), ringTypeName(r.TypeIndex()));
        category.insert(QStringLiteral("ring_kind"), r.ringKind == model::RingKind::Aromatic ? QStringLiteral("Aromatic")
                                                                                             : QStringLiteral("Saturated"));
        m.insert(QStringLiteral("category"), category);
        return m;
    }

    QJsonObject seriesJson() const {
        QJsonObject s;
        s.insert(QStringLiteral("sigma_present"), boolArrayJson(cadence_.sigmaMask()));
        s.insert(QStringLiteral("ring.center"), center_.json());
        s.insert(QStringLiteral("ring.normal"), normal_.json());
        s.insert(QStringLiteral("ring.radius"), radius_.json());
        s.insert(QStringLiteral("ring.flip_state"), flip_.json());
        s.insert(QStringLiteral("ring.pucker_Q"), pucker_q_.json());
        s.insert(QStringLiteral("ring.pucker_theta"), pucker_theta_.json());
        QJsonArray steps;
        for (const auto& hits : near_center_) {
            QJsonArray h;
            for (const Hit& hit : hits) {
                QJsonObject o;
                o.insert(QStringLiteral("atom_index"), hit.atom);
                o.insert(QStringLiteral("distance"), jd(hit.distance));
                o.insert(QStringLiteral("axial_offset"), jd(hit.axial));
                o.insert(QStringLiteral("radial_offset"), jd(hit.radial));
                h.append(o);
            }
            steps.append(h);
        }
        s.insert(QStringLiteral("ring.near_center"), steps);
        return s;
    }

    QJsonObject buildBoostJson() const {
        std::vector<RunningSeriesRef> refs;
        addVec3SeriesRefs(refs, QStringLiteral("ring.center"), center_);
        addVec3SeriesRefs(refs, QStringLiteral("ring.normal"), normal_);
        addScalarSeriesRef(refs, QStringLiteral("ring.radius"), radius_);
        addScalarSeriesRef(refs, QStringLiteral("ring.flip_state"), flip_);
        addScalarSeriesRef(refs, QStringLiteral("ring.pucker_Q"), pucker_q_);
        addScalarSeriesRef(refs, QStringLiteral("ring.pucker_theta"), pucker_theta_);
        QJsonObject boost = boostJson({}, refs, cadence_.sigmaRows(), false);
        last_boost_serial_count = static_cast<std::size_t>(boost.value(QStringLiteral("serial")).toArray().size());
        return boost;
    }

    const AnalysisObjectContext& context_;
    const Body& body_;
    const AnalysisCadence& cadence_;
    std::size_t ring_ = 0;
    double cutoff_ = 8.0;
    Vec3Series center_;
    Vec3Series normal_;
    ScalarSeries radius_;
    ScalarSeries flip_;
    ScalarSeries pucker_q_;
    ScalarSeries pucker_theta_;
    std::vector<std::vector<Hit>> near_center_;
};

AnalysisRing::AnalysisRing(const AnalysisObjectContext& context,
                           std::size_t ringId,
                           double nearCenterCutoffA)
    : AnalysisStructure(context, QStringLiteral("ring"), ringId),
      impl_(std::make_unique<Impl>(context, ringId, nearCenterCutoffA)) {}

void AnalysisRing::Calculate(std::size_t step) { impl_->calculate(step); }
QJsonObject AnalysisRing::Truth() const { return impl_->truth(); }

class AnalysisResidue::Impl {
public:
    Impl(const AnalysisObjectContext& context, std::size_t residue)
        : context_(context),
          body_(context.body()),
          cadence_(context.cadence()),
          residue_(residue),
          phi_(cadence_.stepCount()),
          psi_(cadence_.stepCount()),
          omega_(cadence_.stepCount()),
          chi1_(cadence_.stepCount()),
          chi2_(cadence_.stepCount()),
          chi3_(cadence_.stepCount()),
          chi4_(cadence_.stepCount()),
          ss8_(cadence_.stepCount()),
          ss3_(cadence_.stepCount()),
          ss_observed_(cadence_.stepCount()),
          frames_(cadence_.stepCount()),
          frame_present_(cadence_.stepCount(), false),
          rotamer_bins_(cadence_.stepCount(), std::array<int, 7>{-1, -1, -1, -1, -1, -1, -1}) {}

    void calculate(std::size_t step) {
        if (!body_.run.protein || residue_ >= body_.run.protein->residueCount()) return;
        const model::QtResidue& r = body_.run.protein->residue(residue_);
        foldFrame(step, r);
        const std::size_t atom = canonicalAtom(r);
        foldDihedrals(step, atom);
        const SecondaryStructureState ss = body_.idx.secondaryStructure.state(atom, step);
        ss8_.set(step, static_cast<int>(ss.ss8));
        ss3_.set(step, static_cast<int>(ss.ss3));
        ss_observed_.set(step, ss.observed);
        ++frame_folds_;
    }

    QJsonObject truth() const {
        QJsonObject root;
        root.insert(QStringLiteral("object_type"), QStringLiteral("residue"));
        root.insert(QStringLiteral("identity"), identityJson());
        root.insert(QStringLiteral("model"), modelJson());
        root.insert(QStringLiteral("series"), seriesJson());
        root.insert(QStringLiteral("boost"), buildBoostJson());
        return root;
    }

    std::size_t frameFolds() const { return frame_folds_; }
    mutable std::size_t last_boost_serial_count = 0;

private:
    struct BackboneFrame {
        Vec3 origin = Vec3::Zero();
        Vec3 ex = Vec3::Zero();
        Vec3 ey = Vec3::Zero();
        Vec3 ez = Vec3::Zero();
    };

    std::size_t canonicalAtom(const model::QtResidue& r) const {
        if (r.CA >= 0) return static_cast<std::size_t>(r.CA);
        if (!r.atomIndices.empty() && r.atomIndices.front() >= 0)
            return static_cast<std::size_t>(r.atomIndices.front());
        return 0;
    }

    void foldFrame(std::size_t step, const model::QtResidue& r) {
        if (r.N < 0 || r.CA < 0 || r.C < 0) return;
        const Vec3 n = verbs::pos(body_, static_cast<std::size_t>(r.N), step);
        const Vec3 ca = verbs::pos(body_, static_cast<std::size_t>(r.CA), step);
        const Vec3 c = verbs::pos(body_, static_cast<std::size_t>(r.C), step);
        const LocalFrame f = BuildBackboneCaFrame(ca, n, c);
        if (!f.is_valid) return;
        frames_[step] = {ca, f.x, f.y, f.z};
        frame_present_[step] = true;
    }

    void foldDihedrals(std::size_t step, std::size_t atom) {
        static constexpr std::array<DihedralKind, 7> kinds = {
            DihedralKind::Phi, DihedralKind::Psi, DihedralKind::Omega, DihedralKind::Chi1,
            DihedralKind::Chi2, DihedralKind::Chi3, DihedralKind::Chi4,
        };
        std::array<ScalarSeries*, 7> series = {&phi_, &psi_, &omega_, &chi1_, &chi2_, &chi3_, &chi4_};
        for (std::size_t i = 0; i < kinds.size(); ++i) {
            const DihedralState st = body_.idx.dihedrals.state(kinds[i], atom, step);
            if (st.present && finite(st.radians)) {
                series[i]->set(step, st.radians);
                rotamer_bins_[step][i] = st.fixed_bin;
            }
        }
    }

    QJsonObject identityJson() const {
        const model::QtResidue& r = body_.run.protein->residue(residue_);
        QJsonObject o;
        o.insert(QStringLiteral("uid"), uid(QStringLiteral("residue"), residue_));
        o.insert(QStringLiteral("object_type"), QStringLiteral("residue"));
        o.insert(QStringLiteral("model_index"), static_cast<int>(residue_));
        o.insert(QStringLiteral("residue_index"), static_cast<int>(residue_));
        o.insert(QStringLiteral("chain_id"), r.address.chainId);
        o.insert(QStringLiteral("residue_number_pdb"), r.address.residueNumber);
        o.insert(QStringLiteral("insertion_code"), r.address.insertionCode);
        o.insert(QStringLiteral("amino_acid_ord"), static_cast<int>(r.aminoAcid));
        o.insert(QStringLiteral("amino_acid"), aaName(r.aminoAcid));
        o.insert(QStringLiteral("protonation_variant_index"), r.protonationVariantIndex);
        o.insert(QStringLiteral("terminal_state_ord"), static_cast<int>(r.terminalState));
        o.insert(QStringLiteral("terminal_state"), terminalStateName(r.terminalState));
        o.insert(QStringLiteral("prev_residue_index"), r.prevResidueIndex);
        o.insert(QStringLiteral("next_residue_index"), r.nextResidueIndex);
        o.insert(QStringLiteral("is_proline"), r.isProline);
        o.insert(QStringLiteral("is_aromatic"), r.isAromatic);
        o.insert(QStringLiteral("is_titratable"), r.isTitratable);
        o.insert(QStringLiteral("has_amide_h"), r.hasAmideH);
        QJsonObject bb;
        bb.insert(QStringLiteral("N"), r.N);
        bb.insert(QStringLiteral("CA"), r.CA);
        bb.insert(QStringLiteral("C"), r.C);
        bb.insert(QStringLiteral("O"), r.O);
        bb.insert(QStringLiteral("CB"), r.CB);
        bb.insert(QStringLiteral("H"), r.H);
        bb.insert(QStringLiteral("HA"), r.HA);
        o.insert(QStringLiteral("backbone_atoms"), bb);
        return o;
    }

    QJsonObject modelJson() const {
        const model::QtResidue& r = body_.run.protein->residue(residue_);
        QJsonObject m;
        QJsonObject typing;
        typing.insert(QStringLiteral("amino_acid"), aaName(r.aminoAcid));
        typing.insert(QStringLiteral("amino_acid_ord"), static_cast<int>(r.aminoAcid));
        typing.insert(QStringLiteral("terminal_state"), terminalStateName(r.terminalState));
        m.insert(QStringLiteral("typing"), typing);
        QJsonObject category;
        category.insert(QStringLiteral("amino_acid"), aaName(r.aminoAcid));
        category.insert(QStringLiteral("terminal_state"), terminalStateName(r.terminalState));
        m.insert(QStringLiteral("category"), category);
        return m;
    }

    QJsonObject seriesJson() const {
        QJsonObject s;
        s.insert(QStringLiteral("sigma_present"), boolArrayJson(cadence_.sigmaMask()));
        QJsonArray frameArray;
        for (std::size_t i = 0; i < frames_.size(); ++i) {
            if (!frame_present_[i]) {
                frameArray.append(QJsonValue(QJsonValue::Null));
                continue;
            }
            QJsonObject f;
            f.insert(QStringLiteral("origin"), vec3Json(frames_[i].origin));
            f.insert(QStringLiteral("e_x"), vec3Json(frames_[i].ex));
            f.insert(QStringLiteral("e_y"), vec3Json(frames_[i].ey));
            f.insert(QStringLiteral("e_z"), vec3Json(frames_[i].ez));
            frameArray.append(f);
        }
        s.insert(QStringLiteral("residue.backbone_frame"), frameArray);
        s.insert(QStringLiteral("residue.phi"), phi_.json());
        s.insert(QStringLiteral("residue.psi"), psi_.json());
        s.insert(QStringLiteral("residue.omega"), omega_.json());
        s.insert(QStringLiteral("residue.chi1"), chi1_.json());
        s.insert(QStringLiteral("residue.chi2"), chi2_.json());
        s.insert(QStringLiteral("residue.chi3"), chi3_.json());
        s.insert(QStringLiteral("residue.chi4"), chi4_.json());
        QJsonArray rot;
        for (const auto& bins : rotamer_bins_) {
            QJsonArray r;
            for (int v : bins) r.append(v);
            rot.append(r);
        }
        s.insert(QStringLiteral("residue.rotamer_bins"), rot);
        s.insert(QStringLiteral("residue.ss8"), ss8_.json());
        s.insert(QStringLiteral("residue.ss3"), ss3_.json());
        s.insert(QStringLiteral("residue.ss_observed"), ss_observed_.json());
        return s;
    }

    QJsonObject buildBoostJson() const {
        std::vector<RunningSeriesRef> refs;
        auto frameComp = [&](const QString& axis, int c) {
            std::vector<double> vals(frames_.size(), kNaN);
            for (std::size_t i = 0; i < frames_.size(); ++i) {
                if (!frame_present_[i]) continue;
                if (axis == QStringLiteral("origin")) vals[i] = frames_[i].origin[c];
                if (axis == QStringLiteral("e_x")) vals[i] = frames_[i].ex[c];
                if (axis == QStringLiteral("e_y")) vals[i] = frames_[i].ey[c];
                if (axis == QStringLiteral("e_z")) vals[i] = frames_[i].ez[c];
            }
            refs.push_back({QStringLiteral("residue.backbone_frame"),
                            QStringLiteral("%1.%2").arg(axis, QString::fromLatin1("xyz").mid(c, 1)),
                            QJsonValue(QJsonValue::Null), QJsonValue(QJsonValue::Null),
                            std::move(vals)});
        };
        for (const QString& axis : {QStringLiteral("origin"), QStringLiteral("e_x"),
                                    QStringLiteral("e_y"), QStringLiteral("e_z")})
            for (int c = 0; c < 3; ++c) frameComp(axis, c);
        addScalarSeriesRef(refs, QStringLiteral("residue.phi"), phi_);
        addScalarSeriesRef(refs, QStringLiteral("residue.psi"), psi_);
        addScalarSeriesRef(refs, QStringLiteral("residue.omega"), omega_);
        addScalarSeriesRef(refs, QStringLiteral("residue.chi1"), chi1_);
        addScalarSeriesRef(refs, QStringLiteral("residue.chi2"), chi2_);
        addScalarSeriesRef(refs, QStringLiteral("residue.chi3"), chi3_);
        addScalarSeriesRef(refs, QStringLiteral("residue.chi4"), chi4_);
        QJsonObject boost = boostJson({}, refs, cadence_.sigmaRows(), false);
        last_boost_serial_count = static_cast<std::size_t>(boost.value(QStringLiteral("serial")).toArray().size());
        return boost;
    }

    const AnalysisObjectContext& context_;
    const Body& body_;
    const AnalysisCadence& cadence_;
    std::size_t residue_ = 0;
    ScalarSeries phi_;
    ScalarSeries psi_;
    ScalarSeries omega_;
    ScalarSeries chi1_;
    ScalarSeries chi2_;
    ScalarSeries chi3_;
    ScalarSeries chi4_;
    IntSeries ss8_;
    IntSeries ss3_;
    BoolTriSeries ss_observed_;
    std::vector<BackboneFrame> frames_;
    std::vector<bool> frame_present_;
    std::vector<std::array<int, 7>> rotamer_bins_;
    std::size_t frame_folds_ = 0;
};

AnalysisResidue::AnalysisResidue(const AnalysisObjectContext& context, std::size_t residueIndex)
    : AnalysisStructure(context, QStringLiteral("residue"), residueIndex),
      impl_(std::make_unique<Impl>(context, residueIndex)) {}

void AnalysisResidue::Calculate(std::size_t step) { impl_->calculate(step); }
QJsonObject AnalysisResidue::Truth() const { return impl_->truth(); }

namespace {

QString currentGitCommit() {
    QProcess git;
    git.setProgram(QStringLiteral("git"));
    git.setArguments({QStringLiteral("rev-parse"), QStringLiteral("HEAD")});
    git.setWorkingDirectory(QDir::currentPath());
    git.start();
    if (!git.waitForFinished(1000) || git.exitStatus() != QProcess::NormalExit || git.exitCode() != 0)
        return QStringLiteral("unknown");
    const QString sha = QString::fromUtf8(git.readAllStandardOutput()).trimmed();
    return sha.isEmpty() ? QStringLiteral("unknown") : sha;
}

QJsonArray seriesCatalogJson() {
    struct Item { const char* key; const char* dtype; const char* layout; bool sigma; const char* units; };
    static const std::vector<Item> items = {
        {"coord", "f64", "vec3_xyz", false, "angstrom"},
        {"molecular_frame", "object", "molecular_frame_axes", false, "dimensionless"},
        {"sigma.total", "f64", "spherical_tensor_T0_T1_T2", true, "ppm"},
        {"sigma.dia", "f64", "spherical_tensor_T0_T1_T2", true, "ppm"},
        {"sigma.para", "f64", "spherical_tensor_T0_T1_T2", true, "ppm"},
        {"sigma.total_raw", "f64", "mat3_row_major", true, "ppm"},
        {"sigma.dia_raw", "f64", "mat3_row_major", true, "ppm"},
        {"sigma.para_raw", "f64", "mat3_row_major", true, "ppm"},
        {"sigma_mol.total", "f64", "mat3_row_major_molecular_frame", true, "ppm"},
        {"sigma_mol.dia", "f64", "mat3_row_major_molecular_frame", true, "ppm"},
        {"sigma_mol.para", "f64", "mat3_row_major_molecular_frame", true, "ppm"},
        {"sigma.total_frobenius", "f64", "scalar", true, "ppm"},
        {"sigma.csa_descriptors", "object", "csa_descriptor_arrays", true, "ppm"},
        {"field.mopac_coulomb", "f64", "vec3_xyz", false, "V/angstrom"},
        {"field.ff14sb", "f64", "vec3_xyz", false, "V/angstrom"},
        {"field.apbs", "f64", "vec3_xyz", false, "V/angstrom"},
        {"field.charge_ff14sb", "f64", "vec3_xyz", false, "V/angstrom"},
        {"field_mol.mopac_coulomb", "f64", "vec3_molecular_xyz", false, "V/angstrom"},
        {"field_mol.apbs", "f64", "vec3_molecular_xyz", false, "V/angstrom"},
        {"field_mol.charge_ff14sb", "f64", "vec3_molecular_xyz", false, "V/angstrom"},
        {"field.E_z.mopac_coulomb", "f64", "scalar", false, "V/angstrom"},
        {"field.E_z.apbs", "f64", "scalar", false, "V/angstrom"},
        {"field.E_z.charge_ff14sb", "f64", "scalar", false, "V/angstrom"},
        {"field.E2.mopac_coulomb", "f64", "scalar", false, "(V/angstrom)^2"},
        {"field.E2.apbs", "f64", "scalar", false, "(V/angstrom)^2"},
        {"field.E2.charge_ff14sb", "f64", "scalar", false, "(V/angstrom)^2"},
        {"field.abs_E.mopac_coulomb", "f64", "scalar", false, "V/angstrom"},
        {"field.abs_E.apbs", "f64", "scalar", false, "V/angstrom"},
        {"field.abs_E.charge_ff14sb", "f64", "scalar", false, "V/angstrom"},
        {"field.abs_E2.mopac_coulomb", "f64", "scalar", false, "(V/angstrom)^2"},
        {"field.abs_E2.apbs", "f64", "scalar", false, "(V/angstrom)^2"},
        {"field.abs_E2.charge_ff14sb", "f64", "scalar", false, "(V/angstrom)^2"},
        {"field.charge_ff14sb_source_count", "int32", "scalar", false, "count"},
        {"field.charge_ff14sb_rejected_self", "int32", "scalar", false, "count"},
        {"field.charge_ff14sb_rejected_degenerate", "int32", "scalar", false, "count"},
        {"field.charge_ff14sb_rejected_cutoff", "int32", "scalar", false, "count"},
        {"field.charge_ff14sb_missing_charge", "int32", "scalar", false, "count"},
        {"efg.apbs", "f64", "efg_tensor_T0_T2", false, "ppm"},
        {"efg.aimnet2", "f64", "efg_tensor_T0_T2", false, "ppm"},
        {"efg_mol.apbs", "f64", "mat3_row_major_molecular_frame", false, "V/angstrom^2"},
        {"efg_mol.aimnet2", "f64", "mat3_row_major_molecular_frame", false, "V/angstrom^2"},
        {"efg.abs.apbs", "f64", "scalar", false, "V/angstrom^2"},
        {"efg.abs.aimnet2", "f64", "scalar", false, "V/angstrom^2"},
        {"efg.mopac_coulomb", "f64", "t2_only_T0_null_T2", false, "ppm"},
        {"efg.abs_T2.mopac_coulomb", "f64", "scalar", false, "ppm"},
        {"efg_mol.mopac_coulomb", "f64", "mat3_row_major_molecular_frame", false, "ppm"},
        {"bs_per_type_T2", "f64", "per_ring_type_T2_8x5_flat_40", false, "ppm_T_per_nA"},
        {"hm_per_type_T2", "f64", "per_ring_type_T2_8x5_flat_40", false, "angstrom^-1"},
        {"mc_peptide_co_fixed", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_peptide_co_bo", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_peptide_co_rhombic", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_peptide_cn_fixed", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_peptide_cn_bo", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_backbone_other_fixed", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_backbone_other_bo", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_sidechain_co_fixed", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_sidechain_co_bo", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_sidechain_other_fixed", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_sidechain_other_bo", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_disulfide_fixed", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_disulfide_bo", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_aromatic_zeroed_fixed", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_aromatic_zeroed_bo", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_nearest_co_T2", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"mc_nearest_cn_T2", "f64", "spherical_tensor_T0_T1_T2", false, "angstrom^-3"},
        {"tripeptide_bb_shielding", "f64", "spherical_tensor_T0_T1_T2", false, "ppm"},
        {"larsen_hbond_shielding", "f64", "spherical_tensor_T0_T1_T2", false, "ppm"},
        {"sasa", "f64", "scalar", false, "angstrom^2"},
        {"hbond.nearest_dir", "f64", "vec3_xyz", false, "dimensionless"},
        {"hbond.nearest_dir_mol", "f64", "vec3_molecular_xyz", false, "dimensionless"},
        {"hbond.count", "f64", "scalar", false, "count"},
        {"ff_pb_radius", "f64", "scalar", false, "angstrom"},
        {"bond.lengths", "object", "per_static_bond_length_series_with_rejection_flags", false, "angstrom"},
        {"dihedral.phi", "f64", "scalar", false, "radian"},
        {"dihedral.psi", "f64", "scalar", false, "radian"},
        {"dihedral.omega", "f64", "scalar", false, "radian"},
        {"dihedral.chi1", "f64", "scalar", false, "radian"},
        {"dihedral.chi2", "f64", "scalar", false, "radian"},
        {"dihedral.chi3", "f64", "scalar", false, "radian"},
        {"dihedral.chi4", "f64", "scalar", false, "radian"},
        {"mopac.charge", "f64", "scalar", false, "elementary_charge"},
        {"mopac.s_pop", "f64", "scalar", false, "electron"},
        {"mopac.p_pop", "f64", "scalar", false, "electron"},
        {"mopac.valency", "f64", "scalar", false, "dimensionless_wiberg"},
        {"mopac.s_character", "f64", "scalar", false, "fraction"},
        {"ring.center", "f64", "vec3_xyz", false, "angstrom"},
        {"ring.normal", "f64", "vec3_xyz", false, "angstrom"},
        {"ring.radius", "f64", "scalar", false, "angstrom"},
        {"ring.flip_state", "f64", "scalar", false, "radian"},
        {"ring.pucker_Q", "f64", "scalar", false, "angstrom"},
        {"ring.pucker_theta", "f64", "scalar", false, "degree"},
        {"ring.near_center", "object", "near_center_hits", false, "angstrom"},
        {"residue.backbone_frame", "f64", "backbone_frame", false, "angstrom"},
        {"residue.phi", "f64", "scalar", false, "radian"},
        {"residue.psi", "f64", "scalar", false, "radian"},
        {"residue.omega", "f64", "scalar", false, "radian"},
        {"residue.chi1", "f64", "scalar", false, "radian"},
        {"residue.chi2", "f64", "scalar", false, "radian"},
        {"residue.chi3", "f64", "scalar", false, "radian"},
        {"residue.chi4", "f64", "scalar", false, "radian"},
        {"residue.rotamer_bins", "int8", "rotamer_bins", false, "ordinal"},
        {"residue.ss8", "int8", "enum_scalar", false, "ordinal"},
        {"residue.ss3", "int8", "enum_scalar", false, "ordinal"},
        {"residue.ss_observed", "bool", "bool_scalar", false, "boolean"},
        {"relationships", "f64", "relationship_facets", false, "mixed"},
        {"model.category.hybridisation_series", "int8", "scalar", false, "ordinal"},
        {"model.category.hybridisation_pi_character_series", "f64", "scalar", false, "dimensionless_wiberg"},
        {"model.category.contacted_class_membership_series", "bool", "contact_class_membership", false, "boolean"},
    };
    QJsonArray a;
    for (const Item& item : items) {
        QJsonObject o;
        o.insert(QStringLiteral("key"), QString::fromLatin1(item.key));
        o.insert(QStringLiteral("dtype"), QString::fromLatin1(item.dtype));
        o.insert(QStringLiteral("layout"), QString::fromLatin1(item.layout));
        o.insert(QStringLiteral("sigma_dependent"), item.sigma);
        o.insert(QStringLiteral("units"), QString::fromLatin1(item.units));
        a.append(o);
    }
    return a;
}

QJsonObject manifestJson(const Body& body,
                         const AnalysisCadence& cadence,
                         const QJsonArray& objectInventory,
                         const AnalysisObjectPassDiagnostics& diag) {
    QJsonObject root;
    QJsonObject schema;
    schema.insert(QStringLiteral("name"), QStringLiteral("nmr-shielding-analysis-object-pass"));
    schema.insert(QStringLiteral("version"), QStringLiteral("1.3.0"));
    schema.insert(QStringLiteral("spec_source"), QStringLiteral("FORWARD_SUBSTRATE_PASS_SPEC_2026-06-18.md"));
    schema.insert(QStringLiteral("emitted_utc"), QDateTime::currentDateTimeUtc().toString(Qt::ISODate));
    schema.insert(QStringLiteral("producing_git_commit"), currentGitCommit());
    root.insert(QStringLiteral("schema"), schema);

    QJsonObject run;
    run.insert(QStringLiteral("protein_id"), body.run.manifest.protein_id);
    run.insert(QStringLiteral("calcset_path"), body.run.manifest.calcset_root_abspath);
    run.insert(QStringLiteral("step_count"), static_cast<int>(cadence.stepCount()));
    run.insert(QStringLiteral("sigma_step_count"), static_cast<int>(cadence.sigmaStepCount()));
    QJsonObject counts;
    counts.insert(QStringLiteral("atom"), static_cast<int>(diag.atom_count));
    counts.insert(QStringLiteral("ring"), static_cast<int>(diag.ring_count));
    counts.insert(QStringLiteral("residue"), static_cast<int>(diag.residue_count));
    run.insert(QStringLiteral("object_counts"), counts);
    root.insert(QStringLiteral("run"), run);

    QJsonObject boundedSigma;
    boundedSigma.insert(QStringLiteral("path"), diag.bounded_sigma_path);
    boundedSigma.insert(QStringLiteral("rows"), static_cast<qint64>(diag.bounded_sigma_rows));
    boundedSigma.insert(QStringLiteral("atoms"), static_cast<qint64>(diag.bounded_sigma_atoms));
    boundedSigma.insert(QStringLiteral("bytes"), diag.bounded_sigma_bytes);
    boundedSigma.insert(QStringLiteral("row_granularity"),
                        QStringLiteral("one row per emitted atom per sigma row"));
    boundedSigma.insert(QStringLiteral("molcomp_columns"),
                        QJsonArray{QStringLiteral("molcomp_xx"),
                                   QStringLiteral("molcomp_yy"),
                                   QStringLiteral("molcomp_xy"),
                                   QStringLiteral("molcomp_xz"),
                                   QStringLiteral("molcomp_yz"),
                                   QStringLiteral("molcomp_zz")});
    QJsonObject frameScope;
    frameScope.insert(QStringLiteral("backbone_ca"), QStringLiteral("N-CA-C molecular frame for C-alpha atoms"));
    frameScope.insert(QStringLiteral("aliphatic_carbon"), QStringLiteral("bonded-topology X/plane molecular frame for non-aromatic aliphatic carbons"));
    frameScope.insert(QStringLiteral("hydroxyl_oxygen"), QStringLiteral("O-heavy/H molecular frame for hydroxyl oxygens"));
    frameScope.insert(QStringLiteral("legacy_specific_frames"),
                      QStringLiteral("amide-H, amide-N, carbonyl, aromatic ring, MET SD, carboxylate, guanidinium, carboxamide"));
    boundedSigma.insert(QStringLiteral("molecular_frame_scope"), frameScope);
    root.insert(QStringLiteral("bounded_sigma_sidecar"), boundedSigma);

    QJsonObject classicalSource;
    classicalSource.insert(QStringLiteral("path"), diag.classical_source_path);
    classicalSource.insert(QStringLiteral("rows"), static_cast<qint64>(diag.classical_source_rows));
    classicalSource.insert(QStringLiteral("bytes"), diag.classical_source_bytes);
    classicalSource.insert(QStringLiteral("row_granularity"),
                           QStringLiteral("residue_type + IUPAC atom_name + element + backbone_role + frame_kind"));
    classicalSource.insert(QStringLiteral("leaf_audit_path"), diag.classical_source_leaf_path);
    classicalSource.insert(QStringLiteral("leaf_audit_rows"),
                           static_cast<qint64>(diag.classical_source_leaf_rows));
    classicalSource.insert(QStringLiteral("leaf_audit_bytes"), diag.classical_source_leaf_bytes);
    classicalSource.insert(QStringLiteral("fit"),
                           QStringLiteral("per IUPAC residue/atom group OLS sigma_QM ~ sigma_cl; scale_factor = sd_QM/sd_cl, not OLS slope; sigma_cl = sigma0 + buckingham_linear + buckingham_quadratic + ring + McConnell + Larsen"));
    classicalSource.insert(QStringLiteral("constant_source"),
                           QStringLiteral("src/rediscover/LiteratureConstants.h"));
    classicalSource.insert(QStringLiteral("placeholder_visible_smell"), true);
    root.insert(QStringLiteral("classical_source_terms_sidecar"), classicalSource);

    QJsonObject sourceFamilies;
    sourceFamilies.insert(QStringLiteral("path"), diag.source_family_matrix_path);
    sourceFamilies.insert(QStringLiteral("rows"), static_cast<qint64>(diag.source_family_matrix_rows));
    sourceFamilies.insert(QStringLiteral("bytes"), diag.source_family_matrix_bytes);
    sourceFamilies.insert(QStringLiteral("row_granularity"),
                          QStringLiteral("one row per emitted atom per Axis-1 registered source family"));
    sourceFamilies.insert(QStringLiteral("cross_axis_join_status"),
                          QStringLiteral("axis1_keys_ready;cross_axis_join_blocked_pending_axis2"));
    root.insert(QStringLiteral("source_family_matrices_sidecar"), sourceFamilies);

    QJsonObject subspaceOverlaps;
    subspaceOverlaps.insert(QStringLiteral("path"), diag.subspace_overlap_path);
    subspaceOverlaps.insert(QStringLiteral("rows"), static_cast<qint64>(diag.subspace_overlap_rows));
    subspaceOverlaps.insert(QStringLiteral("bytes"), diag.subspace_overlap_bytes);
    subspaceOverlaps.insert(QStringLiteral("primitive"),
                            QStringLiteral("svd_subspace_compare_v1"));
    subspaceOverlaps.insert(QStringLiteral("read_ids"),
                            QJsonArray{QStringLiteral("G2_field_sources"),
                                       QStringLiteral("G6_field_vs_efg"),
                                       QStringLiteral("G4_bs_vs_hm_tensor_components"),
                                       QStringLiteral("G6_efg_node"),
                                       QStringLiteral("G6_local_electronic_bundle"),
                                       QStringLiteral("G6_mcconnell_axial_vs_rhombic")});
    subspaceOverlaps.insert(QStringLiteral("no_relabel_policy"),
                            QStringLiteral("canonical_corr/basis_dim/explained_fraction/principal_angles only emitted by CompareSubspaces"));
    subspaceOverlaps.insert(QStringLiteral("bs_hm_independence_verdict"),
                            QStringLiteral("producer form check: BS is Biot-Savart current-loops and HM is Haigh-Mallion surface-integral under shared geometry/projection/decomposition; near-null divergence is model robustness, not a relabel"));
    root.insert(QStringLiteral("subspace_overlaps_sidecar"), subspaceOverlaps);

    QJsonObject eta2ByWell;
    eta2ByWell.insert(QStringLiteral("path"), diag.eta2_by_well_path);
    eta2ByWell.insert(QStringLiteral("rows"), static_cast<qint64>(diag.eta2_by_well_rows));
    eta2ByWell.insert(QStringLiteral("bytes"), diag.eta2_by_well_bytes);
    eta2ByWell.insert(QStringLiteral("registry"),
                      QStringLiteral("well families crossed with sigma targets; no DIHEDRAL-only gate"));
    eta2ByWell.insert(QStringLiteral("well_families"),
                      QJsonArray{QStringLiteral("dihedral_rotamer"),
                                 QStringLiteral("contacted_class"),
                                 QStringLiteral("hbond_count_tertiles"),
                                 QStringLiteral("sasa_tertiles"),
                                 QStringLiteral("mopac_charge_tertiles"),
                                 QStringLiteral("field_abs_E_mopac_tertiles")});
    root.insert(QStringLiteral("eta2_by_well_sidecar"), eta2ByWell);

    const LiteratureStatusCounts literatureCounts = CountLiteratureConstantStatuses();
    const LiteratureStatusCounts buckinghamCounts = CountBuckinghamConstantStatuses();
    QJsonObject literature;
    literature.insert(QStringLiteral("header"), QStringLiteral("src/rediscover/LiteratureConstants.h"));
    literature.insert(QStringLiteral("buckingham_header_note"),
                      QStringLiteral("Buckingham convention: sigma=sigma0-A*Eparallel-B*Eparallel^2, E in V/angstrom. Carbonyl 13C/17O signs are Augspurger CO d_sigma/dE values sign-flipped into this convention, but engine Eparallel positive-axis agreement is not verified."));
    literature.insert(QStringLiteral("buckingham_carbonyl_sign_unverified"), true);
    QJsonObject literatureStatus;
    literatureStatus.insert(QStringLiteral("cited"), literatureCounts.cited);
    literatureStatus.insert(QStringLiteral("good_enough"), literatureCounts.good_enough);
    literatureStatus.insert(QStringLiteral("placeholder"), literatureCounts.placeholder);
    literature.insert(QStringLiteral("status_counts"), literatureStatus);
    QJsonObject buckinghamStatus;
    buckinghamStatus.insert(QStringLiteral("cited"), buckinghamCounts.cited);
    buckinghamStatus.insert(QStringLiteral("good_enough"), buckinghamCounts.good_enough);
    buckinghamStatus.insert(QStringLiteral("placeholder"), buckinghamCounts.placeholder);
    literature.insert(QStringLiteral("buckingham_status_counts"), buckinghamStatus);
    literature.insert(QStringLiteral("placeholder_visible_smell"), literatureCounts.placeholder > 0);
    root.insert(QStringLiteral("literature_constants"), literature);

    QJsonObject cad;
    std::vector<std::size_t> originalByStep;
    originalByStep.reserve(cadence.stepCount());
    for (std::size_t i = 0; i < cadence.stepCount(); ++i)
        originalByStep.push_back(cadence.originalIndex(i));
    cad.insert(QStringLiteral("step_original_index"), sizeArrayJson(originalByStep));
    cad.insert(QStringLiteral("sigma_present"), boolArrayJson(cadence.sigmaMask()));
    cad.insert(QStringLiteral("sigma_rows"), sizeArrayJson(cadence.sigmaRows()));
    cad.insert(QStringLiteral("sigma_nominal_stride"),
               cadence.nominalStride() ? QJsonValue(static_cast<int>(*cadence.nominalStride()))
                                        : QJsonValue(QJsonValue::Null));
    root.insert(QStringLiteral("cadence"), cad);

    QJsonObject params;
    params.insert(QStringLiteral("rolling_window_W"), kRollingWindow);
    params.insert(QStringLiteral("null_block_B"), kNullBlock);
    params.insert(QStringLiteral("null_shifts_N"), kNullShifts);
    params.insert(QStringLiteral("null_seed"), QStringLiteral("0xC0DE5EED"));
    params.insert(QStringLiteral("ljung_box_lags_rule"), QStringLiteral("min(20, n/5)"));
    params.insert(QStringLiteral("leverage_top_k"), 1);
    root.insert(QStringLiteral("params"), params);

    QJsonObject units;
    units.insert(QStringLiteral("position"), QStringLiteral("angstrom"));
    units.insert(QStringLiteral("sigma"), QStringLiteral("ppm"));
    units.insert(QStringLiteral("field"), QStringLiteral("V/angstrom"));
    units.insert(QStringLiteral("efg_shielding"), QStringLiteral("ppm"));
    units.insert(QStringLiteral("dihedral"), QStringLiteral("radian"));
    units.insert(QStringLiteral("ring_intensity"), QStringLiteral("nA/T"));
    units.insert(QStringLiteral("ring_jb_offset"), QStringLiteral("angstrom"));
    units.insert(QStringLiteral("pucker_Q"), QStringLiteral("angstrom"));
    units.insert(QStringLiteral("pucker_theta"), QStringLiteral("degree"));
    units.insert(QStringLiteral("bond_order"), QStringLiteral("dimensionless_wiberg"));
    units.insert(QStringLiteral("charge"), QStringLiteral("elementary_charge"));
    root.insert(QStringLiteral("units"), units);

    QJsonObject tensorConventions;
    tensorConventions.insert(QStringLiteral("molcomp_order"),
                             QStringLiteral("xx,yy,xy,xz,yz,zz"));
    tensorConventions.insert(QStringLiteral("molcomp_order_asserted"), true);
    tensorConventions.insert(QStringLiteral("pas_shape_convention"),
                             QStringLiteral("haeberlen_distance_from_isotropic_v1"));
    tensorConventions.insert(QStringLiteral("signed_sorted_haeberlen_eta_refused"), true);
    tensorConventions.insert(QStringLiteral("pas_axis_continuity"),
                             QStringLiteral("trajectory_sign_continuity_used_for_csa_series"));
    root.insert(QStringLiteral("tensor_conventions"), tensorConventions);

    auto ordNameArray = [](const std::vector<QString>& names) {
        QJsonArray a;
        for (int i = 0; i < static_cast<int>(names.size()); ++i) {
            QJsonObject o;
            o.insert(QStringLiteral("ord"), i);
            o.insert(QStringLiteral("name"), names[static_cast<std::size_t>(i)]);
            a.append(o);
        }
        return a;
    };
    QJsonObject catalogs;
    catalogs.insert(QStringLiteral("mechanism"),
                    ordNameArray({QStringLiteral("field_mopac_coulomb"), QStringLiteral("charge_q_over_r3"),
                                  QStringLiteral("mc_lit_valid"), QStringLiteral("ring_jb"),
                                  QStringLiteral("hbond_larsen"), QStringLiteral("efg_coulomb_gradient")}));
    catalogs.insert(QStringLiteral("source_kind"),
                    ordNameArray({QStringLiteral("ff14sb_charge_site"), QStringLiteral("mopac_welford_charge_site"),
                                  QStringLiteral("bond_midpoint"), QStringLiteral("ring_center"),
                                  QStringLiteral("self")}));
    catalogs.insert(QStringLiteral("scope"),
                    ordNameArray({QStringLiteral("self"), QStringLiteral("bonded"),
                                  QStringLiteral("bonded_near_field"), QStringLiteral("through_space")}));
    QJsonArray sourceFamilyRegistry;
    const std::vector<std::pair<QString, QString>> registeredFamilies = {
        {QStringLiteral("field_mopac"), QStringLiteral("MOPAC Coulomb field magnitude, squared magnitude, signed/local components")},
        {QStringLiteral("field_external"), QStringLiteral("APBS and ff14SB charge-field controls")},
        {QStringLiteral("efg_node"), QStringLiteral("AIMNet2/APBS/MOPAC-Coulomb EFG tensors and eigen channels")},
        {QStringLiteral("ring_bs_tensor"), QStringLiteral("Biot-Savart ring-current molecular tensor components")},
        {QStringLiteral("ring_hm_tensor"), QStringLiteral("Haigh-Mallion/Johnson-Bovey ring-current molecular tensor components")},
        {QStringLiteral("mcconnell_axial"), QStringLiteral("McConnell tensor axial projections")},
        {QStringLiteral("mcconnell_rhombic"), QStringLiteral("McConnell tensor rhombic/off-axis projections")},
        {QStringLiteral("local_electronic_population"), QStringLiteral("MOPAC charge, s/p populations, valency")},
        {QStringLiteral("local_electronic_bonding"), QStringLiteral("hybridisation, pi character, dominant Wiberg order")},
        {QStringLiteral("hbond_larsen"), QStringLiteral("Larsen hydrogen-bond shielding and local H-bond geometry")},
    };
    int familyOrd = 0;
    for (const auto& item : registeredFamilies) {
        QJsonObject f;
        f.insert(QStringLiteral("ord"), familyOrd++);
        f.insert(QStringLiteral("family_key"), item.first);
        f.insert(QStringLiteral("description"), item.second);
        f.insert(QStringLiteral("registry_status"), QStringLiteral("registered_axis1"));
        sourceFamilyRegistry.append(f);
    }
    catalogs.insert(QStringLiteral("source_family_registry"), sourceFamilyRegistry);
    catalogs.insert(QStringLiteral("series_catalog"), seriesCatalogJson());
    root.insert(QStringLiteral("catalogs"), catalogs);

    root.insert(QStringLiteral("objects"), objectInventory);
    return root;
}

bool writeJsonFile(const QString& path, const QJsonObject& object, QString* errOut) {
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (errOut) *errOut = QStringLiteral("open failed for %1: %2").arg(path, f.errorString());
        return false;
    }
    f.write(QJsonDocument(object).toJson(QJsonDocument::Compact));
    f.write("\n");
    return true;
}

}  // namespace

bool RunAnalysisObjectPass(const Body& body,
                           const QString& outDir,
                           const AnalysisObjectPassConfig& config,
                           AnalysisObjectPassDiagnostics* diagnostics,
                           QString* errOut) {
    if (!body.run.protein) {
        if (errOut) *errOut = QStringLiteral("analysis object pass requires a loaded protein");
        return false;
    }
    QString pasErr;
    if (!AnalysisAtom::AssertPasShapeConvention(&pasErr)) {
        if (errOut) *errOut = pasErr;
        return false;
    }
    QDir out(outDir);
    if (!out.exists() && !out.mkpath(QStringLiteral("."))) {
        if (errOut) *errOut = QStringLiteral("could not create output dir: %1").arg(outDir);
        return false;
    }
    if (!out.mkpath(QStringLiteral("objects"))) {
        if (errOut) *errOut = QStringLiteral("could not create objects dir under: %1").arg(outDir);
        return false;
    }

    AnalysisCadence cadence(body.run.frameMap);
    AnalysisObjectContext context(body, cadence);
    AnalysisObjectPassDiagnostics diag;
    diag.step_count = cadence.stepCount();
    diag.sigma_step_count = cadence.sigmaStepCount();
    diag.sigma_mask_recorded = true;
    diag.field_vectors_retained = true;
    diag.full_sigma_tensors_retained = true;

    std::vector<std::unique_ptr<AnalysisElement>> objects;
    objects.reserve(body.run.protein->atomCount() + body.run.protein->topology().ringCount()
                    + body.run.protein->residueCount());
    diag.atom_count = body.run.protein->atomCount();
    diag.ring_count = body.run.protein->topology().ringCount();
    diag.residue_count = body.run.protein->residueCount();

    // Resolve the small-run emit selector (work-item 9). Each "residue:atom"
    // string -> the matching atom index. The FULL protein is still constructed
    // and walked (the source environment); only the WRITE is filtered. Empty
    // selector => all atoms emit (production byte-identical path).
    std::set<std::size_t> emitAtoms;
    const bool emitFilterActive = !config.only_atoms.isEmpty();
    if (emitFilterActive) {
        const model::QtProtein& prot = *body.run.protein;
        for (const QString& spec : config.only_atoms) {
            const int colon = spec.indexOf(QLatin1Char(':'));
            if (colon <= 0) {
                if (errOut) *errOut = QStringLiteral("invalid --only-atoms entry (want residue:atom): %1").arg(spec);
                return false;
            }
            const QString resTok = spec.left(colon).trimmed();
            const QString atomName = spec.mid(colon + 1).trimmed();
            // residue token: a PDB residue NUMBER, optionally prefixed by the
            // 3-letter name (e.g. "ASP7" or "7"). Parse the trailing integer.
            int digitStart = resTok.size();
            while (digitStart > 0 && resTok[digitStart - 1].isDigit()) --digitStart;
            bool numOk = false;
            const int resNum = resTok.mid(digitStart).toInt(&numOk);
            if (!numOk) {
                if (errOut) *errOut = QStringLiteral("invalid --only-atoms residue token: %1").arg(spec);
                return false;
            }
            bool found = false;
            for (std::size_t ai = 0; ai < prot.atomCount(); ++ai) {
                const model::QtAtom& a = prot.atom(ai);
                if (!validResidue(prot, a.residueIndex)) continue;
                if (prot.residue(static_cast<std::size_t>(a.residueIndex)).address.residueNumber != resNum)
                    continue;
                if (prot.atomLabel(ai, model::NamingConvention::Iupac) == atomName) {
                    emitAtoms.insert(ai);
                    found = true;
                    break;
                }
            }
            if (!found) {
                if (errOut) *errOut = QStringLiteral("--only-atoms entry not found in protein: %1").arg(spec);
                return false;
            }
        }
        qCInfo(cAnalysisObject).noquote()
            << "analysis_object emit-filter active | requested=" << config.only_atoms.size()
            << "| resolved_atoms=" << emitAtoms.size();
    }

    for (std::size_t atom = 0; atom < body.run.protein->atomCount(); ++atom)
        objects.push_back(std::make_unique<AnalysisAtom>(context, atom, config.per_atom));
    for (std::size_t ring = 0; ring < body.run.protein->topology().ringCount(); ++ring)
        objects.push_back(std::make_unique<AnalysisRing>(context, ring, config.per_atom.ring_cutoff_A));
    for (std::size_t residue = 0; residue < body.run.protein->residueCount(); ++residue)
        objects.push_back(std::make_unique<AnalysisResidue>(context, residue));

    qCInfo(cAnalysisObject).noquote()
        << "analysis_object objects alive | atoms=" << diag.atom_count
        << "| rings=" << diag.ring_count
        << "| residues=" << diag.residue_count
        << "| steps=" << diag.step_count
        << "| sigma_steps=" << diag.sigma_step_count;

    for (std::size_t step = 0; step < cadence.stepCount(); ++step) {
        if (step == 0 || (step % 100) == 0 || step + 1 == cadence.stepCount()) {
            qCInfo(cAnalysisObject).noquote()
                << "analysis_object step begin | step=" << step
                << "| sigma_present=" << cadence.sigmaPresent(step)
                << "| original_frame=" << cadence.originalIndex(step);
        }
        std::size_t objectIndex = 0;
        for (const std::unique_ptr<AnalysisElement>& object : objects) {
            if (step + 1 == cadence.stepCount() && (objectIndex % 100) == 0) {
                qCInfo(cAnalysisObject).noquote()
                    << "analysis_object final-step object | object_index=" << objectIndex
                    << "| type=" << object->objectType()
                    << "| model_index=" << object->modelIndex();
            }
            object->Calculate(step);
            ++diag.calculate_calls;
            ++objectIndex;
        }
        if ((step + 1) % 100 == 0 || step + 1 == cadence.stepCount()) {
            qCInfo(cAnalysisObject).noquote()
                << "analysis_object step complete | step=" << step
                << "| calculate_calls=" << static_cast<qulonglong>(diag.calculate_calls);
        }
    }

    qCInfo(cAnalysisObject).noquote()
        << "analysis_object object walk complete; writing bounded sigma sidecar";
    QString molCompErr;
    if (!AnalysisAtom::AssertMolCompOrder(&molCompErr)) {
        if (errOut) *errOut = molCompErr;
        return false;
    }
    const QString boundedSigmaPath = out.filePath(QStringLiteral("bounded_sigma.csv"));
    {
        QFile f(boundedSigmaPath);
        if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
            if (errOut) {
                *errOut = QStringLiteral("open failed for %1: %2")
                              .arg(boundedSigmaPath, f.errorString());
            }
            return false;
        }
        QTextStream ts(&f);
        AnalysisAtom::WriteBoundedSigmaHeader(ts);
        const QString datasetId = body.run.manifest.dataset_id.isEmpty()
                                      ? body.run.manifest.protein_id
                                      : body.run.manifest.dataset_id;
        const QString proteinId = body.run.manifest.protein_id.isEmpty()
                                      ? (body.run.protein ? body.run.protein->proteinId() : QString())
                                      : body.run.manifest.protein_id;
        for (const std::unique_ptr<AnalysisElement>& object : objects) {
            const auto* atom = dynamic_cast<const AnalysisAtom*>(object.get());
            if (!atom) continue;
            if (emitFilterActive && emitAtoms.count(object->modelIndex()) == 0) continue;
            diag.bounded_sigma_rows += atom->writeBoundedSigmaRows(ts, datasetId, proteinId);
            ++diag.bounded_sigma_atoms;
        }
    }
    diag.bounded_sigma_path = boundedSigmaPath;
    diag.bounded_sigma_bytes = QFileInfo(boundedSigmaPath).size();
    qCInfo(cAnalysisObject).noquote()
        << "analysis_object bounded sigma sidecar | path=" << boundedSigmaPath
        << "| rows=" << static_cast<qulonglong>(diag.bounded_sigma_rows)
        << "| atoms=" << static_cast<qulonglong>(diag.bounded_sigma_atoms)
        << "| bytes=" << diag.bounded_sigma_bytes;

    const QString classicalSourcePath = out.filePath(QStringLiteral("classical_source_terms.csv"));
    const QString classicalSourceLeafPath = out.filePath(QStringLiteral("classical_source_terms_atom_leaf.csv"));
    {
        QFile leaf(classicalSourceLeafPath);
        if (!leaf.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
            if (errOut) {
                *errOut = QStringLiteral("open failed for %1: %2")
                              .arg(classicalSourceLeafPath, leaf.errorString());
            }
            return false;
        }
        QTextStream leafTs(&leaf);
        AnalysisAtom::WriteClassicalSourceTermHeader(leafTs);
        const QString datasetId = body.run.manifest.dataset_id.isEmpty()
                                      ? body.run.manifest.protein_id
                                      : body.run.manifest.dataset_id;
        const QString proteinId = body.run.manifest.protein_id.isEmpty()
                                      ? (body.run.protein ? body.run.protein->proteinId() : QString())
                                      : body.run.manifest.protein_id;
        std::map<QString, std::vector<ClassicalSourceTermRecord>> grouped;
        for (const std::unique_ptr<AnalysisElement>& object : objects) {
            const auto* atom = dynamic_cast<const AnalysisAtom*>(object.get());
            if (!atom) continue;
            if (emitFilterActive && emitAtoms.count(object->modelIndex()) == 0) continue;
            ClassicalSourceTermRecord record = atom->classicalSourceTermRecord(datasetId, proteinId);
            if (record.atom_uid.isEmpty()) continue;
            writeClassicalSourceRecord(leafTs, record);
            ++diag.classical_source_leaf_rows;
            const QString groupKey =
                QStringLiteral("%1\x1f%2\x1f%3\x1f%4\x1f%5")
                    .arg(record.residue_type,
                         record.atom_name,
                         record.element,
                         record.backbone_role,
                         record.frame_kind);
            grouped[groupKey].push_back(std::move(record));
        }
        leaf.close();
        QFile f(classicalSourcePath);
        if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
            if (errOut) {
                *errOut = QStringLiteral("open failed for %1: %2")
                              .arg(classicalSourcePath, f.errorString());
            }
            return false;
        }
        QTextStream ts(&f);
        AnalysisAtom::WriteClassicalSourceTermHeader(ts);
        for (const auto& item : grouped) {
            if (item.second.empty()) continue;
            writeClassicalSourceRecord(ts,
                                       MergeClassicalSourceRecords(item.second.front(),
                                                                  item.second));
            ++diag.classical_source_rows;
        }
    }
    diag.classical_source_path = classicalSourcePath;
    diag.classical_source_bytes = QFileInfo(classicalSourcePath).size();
    diag.classical_source_leaf_path = classicalSourceLeafPath;
    diag.classical_source_leaf_bytes = QFileInfo(classicalSourceLeafPath).size();
    qCInfo(cAnalysisObject).noquote()
        << "analysis_object classical source-term sidecar | path=" << classicalSourcePath
        << "| rows=" << static_cast<qulonglong>(diag.classical_source_rows)
        << "| bytes=" << diag.classical_source_bytes
        << "| leaf_path=" << classicalSourceLeafPath
        << "| leaf_rows=" << static_cast<qulonglong>(diag.classical_source_leaf_rows)
        << "| leaf_bytes=" << diag.classical_source_leaf_bytes;

    const QString sourceFamilyPath = out.filePath(QStringLiteral("source_family_matrices.csv"));
    {
        QFile f(sourceFamilyPath);
        if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
            if (errOut) {
                *errOut = QStringLiteral("open failed for %1: %2")
                              .arg(sourceFamilyPath, f.errorString());
            }
            return false;
        }
        QTextStream ts(&f);
        AnalysisAtom::WriteSourceFamilyMatrixHeader(ts);
        const QString datasetId = body.run.manifest.dataset_id.isEmpty()
                                      ? body.run.manifest.protein_id
                                      : body.run.manifest.dataset_id;
        const QString proteinId = body.run.manifest.protein_id.isEmpty()
                                      ? (body.run.protein ? body.run.protein->proteinId() : QString())
                                      : body.run.manifest.protein_id;
        for (const std::unique_ptr<AnalysisElement>& object : objects) {
            const auto* atom = dynamic_cast<const AnalysisAtom*>(object.get());
            if (!atom) continue;
            if (emitFilterActive && emitAtoms.count(object->modelIndex()) == 0) continue;
            diag.source_family_matrix_rows +=
                atom->writeSourceFamilyMatrixRows(ts, datasetId, proteinId);
        }
    }
    diag.source_family_matrix_path = sourceFamilyPath;
    diag.source_family_matrix_bytes = QFileInfo(sourceFamilyPath).size();
    qCInfo(cAnalysisObject).noquote()
        << "analysis_object source-family matrix sidecar | path=" << sourceFamilyPath
        << "| rows=" << static_cast<qulonglong>(diag.source_family_matrix_rows)
        << "| bytes=" << diag.source_family_matrix_bytes;

    const QString subspaceOverlapPath = out.filePath(QStringLiteral("subspace_overlaps.csv"));
    {
        QFile f(subspaceOverlapPath);
        if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
            if (errOut) {
                *errOut = QStringLiteral("open failed for %1: %2")
                              .arg(subspaceOverlapPath, f.errorString());
            }
            return false;
        }
        QTextStream ts(&f);
        AnalysisAtom::WriteSubspaceOverlapHeader(ts);
        const QString datasetId = body.run.manifest.dataset_id.isEmpty()
                                      ? body.run.manifest.protein_id
                                      : body.run.manifest.dataset_id;
        const QString proteinId = body.run.manifest.protein_id.isEmpty()
                                      ? (body.run.protein ? body.run.protein->proteinId() : QString())
                                      : body.run.manifest.protein_id;
        for (const std::unique_ptr<AnalysisElement>& object : objects) {
            const auto* atom = dynamic_cast<const AnalysisAtom*>(object.get());
            if (!atom) continue;
            if (emitFilterActive && emitAtoms.count(object->modelIndex()) == 0) continue;
            diag.subspace_overlap_rows +=
                atom->writeSubspaceOverlapRows(ts, datasetId, proteinId);
        }
    }
    diag.subspace_overlap_path = subspaceOverlapPath;
    diag.subspace_overlap_bytes = QFileInfo(subspaceOverlapPath).size();
    qCInfo(cAnalysisObject).noquote()
        << "analysis_object subspace-overlap sidecar | path=" << subspaceOverlapPath
        << "| rows=" << static_cast<qulonglong>(diag.subspace_overlap_rows)
        << "| bytes=" << diag.subspace_overlap_bytes;

    const QString etaByWellPath = out.filePath(QStringLiteral("eta2_by_well.csv"));
    {
        QFile f(etaByWellPath);
        if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
            if (errOut) {
                *errOut = QStringLiteral("open failed for %1: %2")
                              .arg(etaByWellPath, f.errorString());
            }
            return false;
        }
        QTextStream ts(&f);
        AnalysisAtom::WriteEtaByWellHeader(ts);
        const QString datasetId = body.run.manifest.dataset_id.isEmpty()
                                      ? body.run.manifest.protein_id
                                      : body.run.manifest.dataset_id;
        const QString proteinId = body.run.manifest.protein_id.isEmpty()
                                      ? (body.run.protein ? body.run.protein->proteinId() : QString())
                                      : body.run.manifest.protein_id;
        for (const std::unique_ptr<AnalysisElement>& object : objects) {
            const auto* atom = dynamic_cast<const AnalysisAtom*>(object.get());
            if (!atom) continue;
            if (emitFilterActive && emitAtoms.count(object->modelIndex()) == 0) continue;
            diag.eta2_by_well_rows += atom->writeEtaByWellRows(ts, datasetId, proteinId);
        }
    }
    diag.eta2_by_well_path = etaByWellPath;
    diag.eta2_by_well_bytes = QFileInfo(etaByWellPath).size();
    qCInfo(cAnalysisObject).noquote()
        << "analysis_object eta2-by-well sidecar | path=" << etaByWellPath
        << "| rows=" << static_cast<qulonglong>(diag.eta2_by_well_rows)
        << "| bytes=" << diag.eta2_by_well_bytes;

    QJsonArray objectInventory;
    for (const std::unique_ptr<AnalysisElement>& object : objects) {
        QJsonObject o;
        o.insert(QStringLiteral("uid"), uid(object->objectType(), object->modelIndex()));
        o.insert(QStringLiteral("object_type"), object->objectType());
        o.insert(QStringLiteral("model_index"), static_cast<int>(object->modelIndex()));
        objectInventory.append(o);
    }

    qCInfo(cAnalysisObject).noquote() << "analysis_object writing truth";
    bool allOxygenOk = true;
    std::size_t writeIndex = 0;
    std::size_t skippedExisting = 0;
    for (std::unique_ptr<AnalysisElement>& object : objects) {
        const QString objectType = object->objectType();
        const std::size_t modelIndex = object->modelIndex();
        if ((writeIndex % 100) == 0) {
            qCInfo(cAnalysisObject).noquote()
                << "analysis_object truth begin | object_index=" << writeIndex
                << "| type=" << objectType
                << "| model_index=" << modelIndex;
        }
        const QString expectedUid = uid(objectType, modelIndex);
        if (expectedUid.isEmpty()) {
            if (errOut) *errOut = QStringLiteral("analysis object has empty expected uid");
            return false;
        }
        const QString path = out.filePath(QStringLiteral("objects/%1.json").arg(expectedUid));
        // Emit-filter (work-item 9): when active, only the selected atoms emit a
        // JSON object. Non-selected atoms + rings + residues are skipped from the
        // WRITE (they still ran the full trajectory as the source environment).
        bool skipped = false;
        if (emitFilterActive) {
            const auto* atomObj = dynamic_cast<const AnalysisAtom*>(object.get());
            skipped = !(atomObj != nullptr && emitAtoms.count(modelIndex) != 0);
        }
        QJsonObject truth;
        if (skipped) {
            ++skippedExisting;
            if ((skippedExisting % 100) == 0) {
                qCInfo(cAnalysisObject).noquote()
                    << "analysis_object truth resume skip | skipped=" << skippedExisting
                    << "| object_index=" << writeIndex
                    << "| uid=" << expectedUid;
            }
        } else {
            truth = object->Truth();
            const QString objectUid = truth.value(QStringLiteral("identity")).toObject().value(QStringLiteral("uid")).toString();
            if (objectUid.isEmpty()) {
                if (errOut) *errOut = QStringLiteral("analysis object emitted empty uid");
                return false;
            }
            if (objectUid != expectedUid) {
                if (errOut) *errOut = QStringLiteral("analysis object uid mismatch: expected %1 got %2").arg(expectedUid, objectUid);
                return false;
            }
            if (!writeJsonFile(path, truth, errOut)) return false;
        }
        if (const auto* atom = dynamic_cast<const AnalysisAtom*>(object.get())) {
            diag.atom_sigma_folds += atom->sigmaFolds();
            diag.atom_relationships += atom->relationshipCount();
            diag.mapped_bonds += atom->mappedBondCount();
            diag.mismatch_events += atom->mismatchEventCount();
            diag.accumulator_responses += atom->accumulatorResponseCount();
            diag.accumulator_contexts += atom->accumulatorContextCount();
            allOxygenOk = allOxygenOk && atom->oxygenGatePassed();
        }
        if (!skipped && objectType == QStringLiteral("ring")) {
            const QJsonArray near = truth.value(QStringLiteral("series")).toObject()
                                        .value(QStringLiteral("ring.near_center")).toArray();
            for (const QJsonValue& stepHits : near)
                diag.ring_near_center_hits += static_cast<std::size_t>(stepHits.toArray().size());
        } else if (objectType == QStringLiteral("residue")) {
            diag.residue_frame_folds += cadence.stepCount();
        }
        object.reset();
        ++writeIndex;
    }
    diag.oxygen_gate_passed = allOxygenOk;
    if (!diag.oxygen_gate_passed) {
        if (errOut) *errOut = QStringLiteral("oxygen hybridisation gate failed: expected sp2 oxygen labelled sp3 or hydroxyl/ether labelled sp2");
        return false;
    }

    const QString manifestPath = out.filePath(QStringLiteral("manifest.json"));
    diag.manifest_path = manifestPath;
    if (!writeJsonFile(manifestPath, manifestJson(body, cadence, objectInventory, diag), errOut)) return false;
    if (skippedExisting > 0) {
        qCInfo(cAnalysisObject).noquote()
            << "analysis_object truth resume skipped existing files | count=" << skippedExisting;
    }
    if (diagnostics) *diagnostics = diag;
    return true;
}

}  // namespace h5reader::rediscover
