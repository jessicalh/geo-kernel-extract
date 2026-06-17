#include "AnalysisAtom.h"

#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "SphericalBasis.h"
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
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QLoggingCategory>
#include <QTextStream>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
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
constexpr int kNullShifts = 200;
constexpr double kCoulombKe = 14.3996;

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

QString ss8Name(SecondaryStructure8 ss8) {
    switch (ss8) {
    case SecondaryStructure8::H: return QStringLiteral("H");
    case SecondaryStructure8::G: return QStringLiteral("G");
    case SecondaryStructure8::I: return QStringLiteral("I");
    case SecondaryStructure8::E: return QStringLiteral("E");
    case SecondaryStructure8::B: return QStringLiteral("B");
    case SecondaryStructure8::T: return QStringLiteral("T");
    case SecondaryStructure8::S: return QStringLiteral("S");
    case SecondaryStructure8::C: return QStringLiteral("C");
    case SecondaryStructure8::Unknown: return QStringLiteral("Unknown");
    }
    return QStringLiteral("Unknown");
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

std::vector<double> componentSeriesT1(const TensorSeries& s, int c) {
    std::vector<double> out(s.values.size(), kNaN);
    for (std::size_t i = 0; i < s.values.size(); ++i)
        if (s.present[i]) out[i] = s.values[i].T1[static_cast<std::size_t>(c)];
    return out;
}

std::vector<double> componentSeriesT2(const TensorSeries& s, int c) {
    std::vector<double> out(s.values.size(), kNaN);
    for (std::size_t i = 0; i < s.values.size(); ++i)
        if (s.present[i]) out[i] = s.values[i].T2[static_cast<std::size_t>(c)];
    return out;
}

std::vector<double> componentSeriesEfgT0(const EfgSeries& s) {
    std::vector<double> out(s.values.size(), kNaN);
    for (std::size_t i = 0; i < s.values.size(); ++i)
        if (s.present[i]) out[i] = s.values[i].t0;
    return out;
}

std::vector<double> componentSeriesEfgT2(const EfgSeries& s, int c) {
    std::vector<double> out(s.values.size(), kNaN);
    for (std::size_t i = 0; i < s.values.size(); ++i)
        if (s.present[i]) out[i] = s.values[i].t2[static_cast<std::size_t>(c)];
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

std::tuple<double, double, double> nullShiftStats(const std::vector<double>& x,
                                                  const std::vector<double>& y,
                                                  double observedSlope) {
    if (x.size() != y.size() || x.size() < 2) return {kNaN, kNaN, kNaN};
    std::mt19937 rng(kNullSeed);
    const int n = static_cast<int>(x.size());
    const int blocks = std::max(1, (n + kNullBlock - 1) / kNullBlock);
    std::uniform_int_distribution<int> dist(1, std::max(1, blocks - 1));
    std::vector<double> slopes;
    slopes.reserve(kNullShifts);
    std::vector<double> shifted(x.size(), kNaN);
    for (int i = 0; i < kNullShifts; ++i) {
        const int shift = (dist(rng) * kNullBlock) % n;
        for (int j = 0; j < n; ++j) shifted[static_cast<std::size_t>(j)] = x[static_cast<std::size_t>((j + shift) % n)];
        const OlsResult r = ols(shifted, y);
        if (finite(r.slope)) slopes.push_back(r.slope);
    }
    if (slopes.empty()) return {kNaN, kNaN, kNaN};
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
    return {mean, sd, z};
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
    auto [nullMean, nullSd, obsZ] = nullShiftStats(x, y, fit.slope);
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

void addVec3SeriesRefs(std::vector<RunningSeriesRef>& refs, const QString& key, const Vec3Series& s) {
    static constexpr std::array<const char*, 3> comps = {"x", "y", "z"};
    for (int c = 0; c < 3; ++c)
        refs.push_back({key, QString::fromLatin1(comps[static_cast<std::size_t>(c)]),
                        QJsonValue(QJsonValue::Null), QJsonValue(QJsonValue::Null),
                        componentSeries(s, c)});
}

void addScalarSeriesRef(std::vector<RunningSeriesRef>& refs, const QString& key, const ScalarSeries& s) {
    refs.push_back({key, QJsonValue(QJsonValue::Null), QJsonValue(QJsonValue::Null),
                    QJsonValue(QJsonValue::Null), componentSeries(s)});
}

void addEfgSeriesRefs(std::vector<RunningSeriesRef>& refs, const QString& key, const EfgSeries& s) {
    refs.push_back({key, QStringLiteral("T0"), QJsonValue(QJsonValue::Null),
                    QJsonValue(QJsonValue::Null), componentSeriesEfgT0(s)});
    for (int c = 0; c < 5; ++c)
        refs.push_back({key, QStringLiteral("T2.c%1").arg(c), QJsonValue(QJsonValue::Null),
                        QJsonValue(QJsonValue::Null), componentSeriesEfgT2(s, c)});
}

std::vector<std::pair<QString, std::vector<double>>> sigmaComponents(const TensorSeries& total,
                                                                      const TensorSeries& dia,
                                                                      const TensorSeries& para,
                                                                      const ScalarSeries& frob) {
    std::vector<std::pair<QString, std::vector<double>>> out;
    auto addTensor = [&](const QString& prefix, const TensorSeries& s) {
        out.push_back({prefix + QStringLiteral(".T0"), componentSeriesT0(s)});
        for (int i = 0; i < 3; ++i)
            out.push_back({prefix + QStringLiteral(".T1.c%1").arg(i), componentSeriesT1(s, i)});
        for (int i = 0; i < 5; ++i)
            out.push_back({prefix + QStringLiteral(".T2.c%1").arg(i), componentSeriesT2(s, i)});
    };
    addTensor(QStringLiteral("total"), total);
    out.push_back({QStringLiteral("total.frobenius"), componentSeries(frob)});
    addTensor(QStringLiteral("dia"), dia);
    addTensor(QStringLiteral("para"), para);
    return out;
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
    const bool carriesT0 = kind == io::FieldKind::MOPACCoulombEFG;
    if (carriesT0 && body.catalog.present(body, kind, atom, row, 0)) {
        if (const std::optional<double> t0 = body.catalog.value(body, kind, atom, row, 0))
            out.t0 = *t0;
    }
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

QJsonArray scalarAsSerialJson(const std::vector<double>& values) {
    QJsonArray a;
    for (double v : values) a.append(jd(v));
    return a;
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
    for (const MopacFrameBond& b : mopacBonds(step)) {
        if (b.atom_a == p.first && b.atom_b == p.second) return b.wiberg_order;
    }
    return std::nullopt;
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
          sigma_total_(cadence_.stepCount()),
          sigma_dia_(cadence_.stepCount()),
          sigma_para_(cadence_.stepCount()),
          sigma_frob_(cadence_.stepCount()),
          field_mopac_(cadence_.stepCount()),
          field_ff14sb_(cadence_.stepCount()),
          efg_apbs_(cadence_.stepCount()),
          efg_aimnet2_(cadence_.stepCount()),
          efg_mopac_(cadence_.stepCount()),
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
          pi_character_(cadence_.stepCount()) {
        buildStaticBonds();
        readStaticHybridisation();
    }

    void calculate(std::size_t step) {
        if (!body_.run.protein || atom_ >= body_.run.protein->atomCount()) return;
        coord_.set(step, verbs::pos(body_, atom_, step));
        foldMopacSelfState(step);
        foldTopologyJoin(step);
        foldGeometry(step);
        foldDihedrals(step);
        foldFieldsAndRelationships(step);
        foldEfg(step);
        if (cadence_.sigmaPresent(step)) foldSigma(step);
    }

    QJsonObject truth() const {
        QJsonObject root;
        qCInfo(cAnalysisObject).noquote() << "analysis_atom truth section | atom=" << atom_ << "| section=object_type";
        root.insert(QStringLiteral("object_type"), QStringLiteral("atom"));
        qCInfo(cAnalysisObject).noquote() << "analysis_atom truth section | atom=" << atom_ << "| section=identity";
        root.insert(QStringLiteral("identity"), identityJson());
        qCInfo(cAnalysisObject).noquote() << "analysis_atom truth section | atom=" << atom_ << "| section=model";
        root.insert(QStringLiteral("model"), modelJson());
        qCInfo(cAnalysisObject).noquote() << "analysis_atom truth section | atom=" << atom_ << "| section=series";
        root.insert(QStringLiteral("series"), seriesJson());
        qCInfo(cAnalysisObject).noquote() << "analysis_atom truth section | atom=" << atom_ << "| section=boost";
        root.insert(QStringLiteral("boost"), buildBoostJson());
        qCInfo(cAnalysisObject).noquote() << "analysis_atom truth section | atom=" << atom_ << "| section=done";
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
    mutable std::size_t last_boost_coupling_count = 0;
    mutable std::size_t last_boost_serial_count = 0;

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

    QJsonObject seriesJson() const {
        QJsonObject s;
        s.insert(QStringLiteral("sigma_present"), boolArrayJson(cadence_.sigmaMask()));
        s.insert(QStringLiteral("coord"), coord_.json());
        s.insert(QStringLiteral("sigma.total"), sigma_total_.json());
        s.insert(QStringLiteral("sigma.dia"), sigma_dia_.json());
        s.insert(QStringLiteral("sigma.para"), sigma_para_.json());
        s.insert(QStringLiteral("sigma.total_frobenius"), sigma_frob_.json());
        s.insert(QStringLiteral("field.mopac_coulomb"), field_mopac_.json());
        s.insert(QStringLiteral("field.ff14sb"), field_ff14sb_.json());
        s.insert(QStringLiteral("efg.apbs"), efg_apbs_.json());
        s.insert(QStringLiteral("efg.aimnet2"), efg_aimnet2_.json());
        s.insert(QStringLiteral("efg.mopac_coulomb"), efg_mopac_.json());
        s.insert(QStringLiteral("dihedral.phi"), phi_.json());
        s.insert(QStringLiteral("dihedral.psi"), psi_.json());
        s.insert(QStringLiteral("dihedral.omega"), omega_.json());
        s.insert(QStringLiteral("dihedral.chi1"), chi1_.json());
        s.insert(QStringLiteral("dihedral.chi2"), chi2_.json());
        s.insert(QStringLiteral("dihedral.chi3"), chi3_.json());
        s.insert(QStringLiteral("dihedral.chi4"), chi4_.json());
        s.insert(QStringLiteral("mopac.charge"), mopac_charge_.json());
        s.insert(QStringLiteral("mopac.s_pop"), mopac_s_pop_.json());
        s.insert(QStringLiteral("mopac.p_pop"), mopac_p_pop_.json());
        s.insert(QStringLiteral("mopac.valency"), mopac_valency_.json());
        s.insert(QStringLiteral("mopac.s_character"), mopac_s_character_.json());
        QJsonArray rels;
        int idx = 0;
        for (const auto& item : relationships_) rels.append(relationshipJson(item.first, item.second, idx++));
        s.insert(QStringLiteral("relationships"), rels);
        return s;
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
        o.insert(QStringLiteral("present"), rel.presentJson());
        QJsonObject facets;
        facets.insert(QStringLiteral("distance"), rel.distance.json(cadence_.stepCount()));
        facets.insert(QStringLiteral("cos_theta"), rel.cos_theta.json(cadence_.stepCount()));
        facets.insert(QStringLiteral("inv_r3"), rel.inv_r3.json(cadence_.stepCount()));
        facets.insert(QStringLiteral("dipolar"), rel.dipolar.json(cadence_.stepCount()));
        facets.insert(QStringLiteral("kernel_T0"), rel.kernel_T0.json(cadence_.stepCount()));
        facets.insert(QStringLiteral("kernel_T2"), rel.kernel_T2.json(cadence_.stepCount()));
        facets.insert(QStringLiteral("contribution"), rel.contribution.json(cadence_.stepCount()));
        o.insert(QStringLiteral("facets"), facets);
        return o;
    }

    QJsonObject buildBoostJson() const {
        std::vector<RunningSeriesRef> refs;
        addVec3SeriesRefs(refs, QStringLiteral("coord"), coord_);
        addVec3SeriesRefs(refs, QStringLiteral("field.mopac_coulomb"), field_mopac_);
        addVec3SeriesRefs(refs, QStringLiteral("field.ff14sb"), field_ff14sb_);
        addEfgSeriesRefs(refs, QStringLiteral("efg.apbs"), efg_apbs_);
        addEfgSeriesRefs(refs, QStringLiteral("efg.aimnet2"), efg_aimnet2_);
        addEfgSeriesRefs(refs, QStringLiteral("efg.mopac_coulomb"), efg_mopac_);
        addScalarSeriesRef(refs, QStringLiteral("dihedral.phi"), phi_);
        addScalarSeriesRef(refs, QStringLiteral("dihedral.psi"), psi_);
        addScalarSeriesRef(refs, QStringLiteral("dihedral.omega"), omega_);
        addScalarSeriesRef(refs, QStringLiteral("dihedral.chi1"), chi1_);
        addScalarSeriesRef(refs, QStringLiteral("dihedral.chi2"), chi2_);
        addScalarSeriesRef(refs, QStringLiteral("dihedral.chi3"), chi3_);
        addScalarSeriesRef(refs, QStringLiteral("dihedral.chi4"), chi4_);
        addScalarSeriesRef(refs, QStringLiteral("mopac.charge"), mopac_charge_);
        addScalarSeriesRef(refs, QStringLiteral("mopac.s_pop"), mopac_s_pop_);
        addScalarSeriesRef(refs, QStringLiteral("mopac.p_pop"), mopac_p_pop_);
        addScalarSeriesRef(refs, QStringLiteral("mopac.valency"), mopac_valency_);
        addScalarSeriesRef(refs, QStringLiteral("mopac.s_character"), mopac_s_character_);

        const auto sigmas = sigmaComponents(sigma_total_, sigma_dia_, sigma_para_, sigma_frob_);
        qCInfo(cAnalysisObject).noquote()
            << "analysis_atom boost begin | atom=" << atom_
            << "| base_refs=" << refs.size()
            << "| relationships=" << relationships_.size()
            << "| sigma_components=" << sigmas.size();
        QJsonObject out = boostJson(sigmas, refs, cadence_.sigmaRows(), true);
        QJsonArray couplingArray = out.value(QStringLiteral("coupling")).toArray();
        QJsonArray serialArray = out.value(QStringLiteral("serial")).toArray();
        qCInfo(cAnalysisObject).noquote()
            << "analysis_atom boost base complete | atom=" << atom_
            << "| coupling=" << couplingArray.size()
            << "| serial=" << serialArray.size();

        int relIndex = 0;
        for (const auto& item : relationships_) {
            const RelationshipSeries& r = item.second;
            const QJsonValue ri(relIndex);
            auto appendRel = [&](const QString& facet,
                                 QJsonValue component,
                                 std::vector<double> values) {
                RunningSeriesRef ref{QStringLiteral("relationships"), std::move(component), ri,
                                     facet, std::move(values)};
                for (const auto& s : sigmas) {
                    couplingArray.append(couplingRecord(s.first, s.second, ref, cadence_.sigmaRows()));
                }
                serialArray.append(serialRecord(ref));
            };
            appendRel(QStringLiteral("distance"), QJsonValue(QJsonValue::Null),
                      r.distance.dense(cadence_.stepCount()));
            appendRel(QStringLiteral("cos_theta"), QJsonValue(QJsonValue::Null),
                      r.cos_theta.dense(cadence_.stepCount()));
            appendRel(QStringLiteral("inv_r3"), QJsonValue(QJsonValue::Null),
                      r.inv_r3.dense(cadence_.stepCount()));
            appendRel(QStringLiteral("dipolar"), QJsonValue(QJsonValue::Null),
                      r.dipolar.dense(cadence_.stepCount()));
            appendRel(QStringLiteral("kernel_T0"), QJsonValue(QJsonValue::Null),
                      r.kernel_T0.dense(cadence_.stepCount()));
            appendRel(QStringLiteral("contribution"), QJsonValue(QJsonValue::Null),
                      r.contribution.dense(cadence_.stepCount()));
            for (int c = 0; c < 5; ++c) {
                appendRel(QStringLiteral("kernel_T2"), QStringLiteral("T2.c%1").arg(c),
                          r.kernel_T2.componentDense(cadence_.stepCount(), c));
            }
            if (relIndex > 0 && relIndex % 100 == 0) {
                qCInfo(cAnalysisObject).noquote()
                    << "analysis_atom boost relationships | atom=" << atom_
                    << "| relationship_index=" << relIndex
                    << "| coupling=" << couplingArray.size()
                    << "| serial=" << serialArray.size();
            }
            ++relIndex;
        }

        out.insert(QStringLiteral("coupling"), couplingArray);
        out.insert(QStringLiteral("serial"), serialArray);
        qCInfo(cAnalysisObject).noquote()
            << "analysis_atom boost complete | atom=" << atom_
            << "| coupling=" << couplingArray.size()
            << "| serial=" << serialArray.size();
        last_boost_coupling_count = static_cast<std::size_t>(out.value(QStringLiteral("coupling")).toArray().size());
        last_boost_serial_count = static_cast<std::size_t>(out.value(QStringLiteral("serial")).toArray().size());
        return out;
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
            static_bonds_[pairKey(a, bb)] = {b.bondIndex, b.order, b.category, other};
        }
    }

    void readStaticHybridisation() {
        if (body_.catalog.present(body_, io::FieldKind::EnrichmentHybridisation, atom_, 0, -1)) {
            const std::optional<double> h =
                body_.catalog.value(body_, io::FieldKind::EnrichmentHybridisation, atom_, 0, -1);
            if (h && finite(*h)) static_hybridisation_ord_ = static_cast<int>(*h);
        }
    }

    void foldGeometry(std::size_t step) {
        (void)step;
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
        if (body_.catalog.present(body_, ArrayId::MopacCoulombEfield, atom_, step))
            field_mopac_.set(step, body_.catalog.valueVec3(body_, ArrayId::MopacCoulombEfield, atom_, step));
    }

    void foldRelationship(std::size_t step, const PairContribution& pair) {
        const int hyb = hybridisation_.values[step] == IntSeries::kMissing
                            ? static_cast<int>(model::Hybridisation::Unassigned)
                            : hybridisation_.values[step];
        const model::QtAtom& atom = body_.run.protein->atom(atom_);
        const auto cr = contactedResidue(body_, pair);
        if (cr.first >= 0) contacted_residues_.insert({cr.first, static_cast<int>(cr.second)});
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
        r.markPresent(step, (pair.pointer_flags & PresentFlag) != 0);
        r.distance.set(step, pair.r);
        r.cos_theta.set(step, pair.cos_theta);
        r.inv_r3.set(step, pair.inv_r3);
        r.dipolar.set(step, pair.dipolar);
        r.kernel_T0.set(step, pair.kernel_T0);
        r.contribution.set(step, pair.contribution);
        r.kernel_T2.set(step, pair.kernel_T2);
    }

    void foldEfg(std::size_t step) {
        efg_apbs_.set(step, readT2Field(body_, io::FieldKind::APBSEFG, ArrayId::ApbsEfg, atom_, step));
        efg_aimnet2_.set(step, readT2Field(body_, io::FieldKind::AIMNet2EFG, ArrayId::Aimnet2Efg, atom_, step));
        efg_mopac_.set(step, readT2Field(body_, io::FieldKind::MOPACCoulombEFG,
                                         ArrayId::MopacCoulombShielding, atom_, step));
    }

    void foldSigma(std::size_t step) {
        const std::size_t original = cadence_.originalIndex(step);
        const DftTarget target = BuildTarget(body_.run, atom_, original, LocalFrame{});
        if (!target.present || !finite(target.total_decomp.T0)) return;
        sigma_total_.set(step, target.total_decomp);
        sigma_dia_.set(step, target.dia_decomp);
        sigma_para_.set(step, target.para_decomp);
        sigma_frob_.set(step, tensorFrobenius(target.total_raw));
        ++sigma_folds_;
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
    TensorSeries sigma_total_;
    TensorSeries sigma_dia_;
    TensorSeries sigma_para_;
    ScalarSeries sigma_frob_;
    Vec3Series field_mopac_;
    Vec3Series field_ff14sb_;
    EfgSeries efg_apbs_;
    EfgSeries efg_aimnet2_;
    EfgSeries efg_mopac_;
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
    ScalarSeries pi_character_;
    std::map<RelationshipKey, RelationshipSeries> relationships_;
    std::map<std::uint64_t, StaticBondInfo> static_bonds_;
    std::map<std::uint64_t, ScalarSeries> bond_series_;
    std::set<std::pair<int, int>> contacted_residues_;
    int static_hybridisation_ord_ = static_cast<int>(model::Hybridisation::Unassigned);
    bool oxygen_gate_checked_ = false;
    std::size_t sigma_folds_ = 0;
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
std::size_t AnalysisAtom::boostCouplingCount() const { return impl_->last_boost_coupling_count; }
std::size_t AnalysisAtom::boostSerialCount() const { return impl_->last_boost_serial_count; }
bool AnalysisAtom::oxygenGatePassed() const { return impl_->oxygenGatePassed(); }

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

QJsonArray seriesCatalogJson() {
    struct Item { const char* key; const char* dtype; const char* layout; bool sigma; const char* units; };
    static constexpr std::array<Item, 44> items = {{
        {"coord", "f64", "vec3_xyz", false, "angstrom"},
        {"sigma.total", "f64", "spherical_tensor_T0_T1_T2", true, "ppm"},
        {"sigma.dia", "f64", "spherical_tensor_T0_T1_T2", true, "ppm"},
        {"sigma.para", "f64", "spherical_tensor_T0_T1_T2", true, "ppm"},
        {"sigma.total_frobenius", "f64", "scalar", true, "ppm"},
        {"field.mopac_coulomb", "f64", "vec3_xyz", false, "V/angstrom"},
        {"field.ff14sb", "f64", "vec3_xyz", false, "V/angstrom"},
        {"efg.apbs", "f64", "efg_tensor_T0_T2", false, "ppm"},
        {"efg.aimnet2", "f64", "efg_tensor_T0_T2", false, "ppm"},
        {"efg.mopac_coulomb", "f64", "efg_tensor_T0_T2", false, "ppm"},
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
    }};
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
                         const std::vector<std::unique_ptr<AnalysisElement>>& objects,
                         const AnalysisObjectPassDiagnostics& diag) {
    QJsonObject root;
    QJsonObject schema;
    schema.insert(QStringLiteral("name"), QStringLiteral("nmr-shielding-analysis-object-pass"));
    schema.insert(QStringLiteral("version"), QStringLiteral("1.0.0"));
    schema.insert(QStringLiteral("spec_source"), QStringLiteral("COLLECTION_PASS_SPEC_2026-06-16.md"));
    schema.insert(QStringLiteral("emitted_utc"), QDateTime::currentDateTimeUtc().toString(Qt::ISODate));
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
    catalogs.insert(QStringLiteral("series_catalog"), seriesCatalogJson());
    root.insert(QStringLiteral("catalogs"), catalogs);

    QJsonArray inventory;
    for (const auto& obj : objects) {
        QJsonObject o;
        o.insert(QStringLiteral("uid"), uid(obj->objectType(), obj->modelIndex()));
        o.insert(QStringLiteral("object_type"), obj->objectType());
        o.insert(QStringLiteral("model_index"), static_cast<int>(obj->modelIndex()));
        inventory.append(o);
    }
    root.insert(QStringLiteral("objects"), inventory);
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

bool looksLikeCompletedJsonFile(const QString& path) {
    QFile f(path);
    if (!f.exists()) return false;
    if (!f.open(QIODevice::ReadOnly)) return false;
    if (f.size() <= 1) return false;
    if (!f.seek(f.size() - 1)) return false;
    return f.read(1) == QByteArrayLiteral("\n");
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

    qCInfo(cAnalysisObject).noquote() << "analysis_object object walk complete; writing truth";
    bool allOxygenOk = true;
    std::size_t writeIndex = 0;
    std::size_t skippedExisting = 0;
    for (const std::unique_ptr<AnalysisElement>& object : objects) {
        if ((writeIndex % 100) == 0) {
            qCInfo(cAnalysisObject).noquote()
                << "analysis_object truth begin | object_index=" << writeIndex
                << "| type=" << object->objectType()
                << "| model_index=" << object->modelIndex();
        }
        const QString expectedUid = uid(object->objectType(), object->modelIndex());
        if (expectedUid.isEmpty()) {
            if (errOut) *errOut = QStringLiteral("analysis object has empty expected uid");
            return false;
        }
        const QString path = out.filePath(QStringLiteral("objects/%1.json").arg(expectedUid));
        const bool skipped = looksLikeCompletedJsonFile(path);
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
        if (!skipped) {
            const QJsonObject boost = truth.value(QStringLiteral("boost")).toObject();
            diag.boost_coupling_results += static_cast<std::size_t>(boost.value(QStringLiteral("coupling")).toArray().size());
            diag.boost_serial_results += static_cast<std::size_t>(boost.value(QStringLiteral("serial")).toArray().size());
        }
        if (const auto* atom = dynamic_cast<const AnalysisAtom*>(object.get())) {
            diag.atom_sigma_folds += atom->sigmaFolds();
            diag.atom_relationships += atom->relationshipCount();
            diag.mapped_bonds += atom->mappedBondCount();
            diag.mismatch_events += atom->mismatchEventCount();
            allOxygenOk = allOxygenOk && atom->oxygenGatePassed();
        }
        if (!skipped && object->objectType() == QStringLiteral("ring")) {
            const QJsonArray near = truth.value(QStringLiteral("series")).toObject()
                                        .value(QStringLiteral("ring.near_center")).toArray();
            for (const QJsonValue& stepHits : near)
                diag.ring_near_center_hits += static_cast<std::size_t>(stepHits.toArray().size());
        } else if (object->objectType() == QStringLiteral("residue")) {
            diag.residue_frame_folds += cadence.stepCount();
        }
        ++writeIndex;
    }
    diag.oxygen_gate_passed = allOxygenOk;
    if (!diag.oxygen_gate_passed) {
        if (errOut) *errOut = QStringLiteral("oxygen hybridisation gate failed: expected sp2 oxygen labelled sp3 or hydroxyl/ether labelled sp2");
        return false;
    }

    const QString manifestPath = out.filePath(QStringLiteral("manifest.json"));
    diag.manifest_path = manifestPath;
    if (!writeJsonFile(manifestPath, manifestJson(body, cadence, objects, diag), errOut)) return false;
    if (skippedExisting > 0) {
        qCInfo(cAnalysisObject).noquote()
            << "analysis_object truth resume skipped existing files | count=" << skippedExisting;
    }
    if (diagnostics) *diagnostics = diag;
    return true;
}

}  // namespace h5reader::rediscover
