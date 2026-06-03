#include "CaseHunter.h"

#include "ExtractionSupport.h"

#include "../model/QtAtom.h"
#include "../model/QtProtein.h"

#include <QDir>
#include <QSaveFile>
#include <QTextStream>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace h5reader::rediscover {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

enum class HunterMechanism : int {
    Ring = 0,
    Charge = 1,
    Mc = 2,
};

struct MechanismSums {
    double ring = 0.0;
    double charge = 0.0;
    double mc = 0.0;
};

struct InputCandidate {
    std::size_t atom = 0;
    FrameWindow window;
    HunterMechanism mechanism = HunterMechanism::Ring;
    double isolation_fraction = 0.0;
    double mean_dominant_fraction = kNaN;
    double min_gap_A = kNaN;
    double driver_variance = kNaN;
    double quiet_variance = kNaN;
    double score = kNaN;
    std::size_t input_frames = 0;
};

struct MeasuredCandidate : InputCandidate {
    double dft_recovery_R2 = kNaN;
    std::size_t dft_points = 0;
};

QString mechanismName(HunterMechanism m) {
    switch (m) {
    case HunterMechanism::Ring: return QStringLiteral("ring");
    case HunterMechanism::Charge: return QStringLiteral("charge");
    case HunterMechanism::Mc: return QStringLiteral("mc");
    }
    return QStringLiteral("unknown");
}

double finiteOrZero(double v) {
    return std::isfinite(v) ? v : 0.0;
}

double variance(const std::vector<double>& xs) {
    double mean = 0.0;
    double m2 = 0.0;
    std::size_t n = 0;
    for (double x : xs) {
        if (!std::isfinite(x)) continue;
        ++n;
        const double delta = x - mean;
        mean += delta / static_cast<double>(n);
        m2 += delta * (x - mean);
    }
    return n > 1 ? m2 / static_cast<double>(n - 1) : kNaN;
}

double correlationR2(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) return kNaN;
    double sx = 0.0;
    double sy = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
        sx += x[i];
        sy += y[i];
        ++n;
    }
    if (n < 3) return kNaN;
    const double mx = sx / static_cast<double>(n);
    const double my = sy / static_cast<double>(n);
    double sxx = 0.0;
    double syy = 0.0;
    double sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if (!(sxx > 0.0) || !(syy > 0.0)) return kNaN;
    const double r = sxy / std::sqrt(sxx * syy);
    return r * r;
}

double t2MagnitudeLocal(const std::array<double, 5>& t2) {
    double s = 0.0;
    for (double v : t2) {
        if (!std::isfinite(v)) return kNaN;
        s += v * v;
    }
    return std::sqrt(s);
}

MechanismSums mechanismSums(const Body& body,
                            const PerAtomSubstrateConfig& cfg,
                            std::size_t atom,
                            std::size_t row) {
    MechanismSums sums;
    for (const PairContribution& p : PerAtomRowPairContributions(body, atom, row, cfg, LocalFrame{})) {
        const double v = finiteOrZero(p.contribution);
        if (p.mechanism == QStringLiteral("ring_jb")) sums.ring += v;
        else if (p.mechanism == QStringLiteral("charge_q_over_r3")) sums.charge += v;
        else if (p.mechanism == QStringLiteral("mc_lit_valid")) sums.mc += v;
    }
    return sums;
}

double driverValue(const MechanismSums& sums, HunterMechanism m) {
    switch (m) {
    case HunterMechanism::Ring: return sums.ring;
    case HunterMechanism::Charge: return sums.charge;
    case HunterMechanism::Mc: return sums.mc;
    }
    return kNaN;
}

double otherValue(const MechanismSums& sums, HunterMechanism m) {
    switch (m) {
    case HunterMechanism::Ring: return sums.charge + sums.mc;
    case HunterMechanism::Charge: return sums.ring + sums.mc;
    case HunterMechanism::Mc: return sums.ring + sums.charge;
    }
    return kNaN;
}

double dominantFor(const PerAtomIsolationScalars& iso, HunterMechanism m) {
    switch (m) {
    case HunterMechanism::Ring: return iso.dominant_fraction_ring;
    case HunterMechanism::Charge: return iso.dominant_fraction_charge;
    case HunterMechanism::Mc: return iso.dominant_fraction_mc;
    }
    return kNaN;
}

double gapFor(const PerAtomIsolationScalars& iso, HunterMechanism m) {
    switch (m) {
    case HunterMechanism::Ring: return iso.gap_to_2nd_ring_r;
    case HunterMechanism::Charge: return iso.gap_to_2nd_charge_r;
    case HunterMechanism::Mc: return iso.gap_to_2nd_bond_r;
    }
    return kNaN;
}

std::vector<int32_t> allAtomScope(const Body& body) {
    std::vector<int32_t> out;
    if (!body.run.protein) return out;
    out.reserve(body.run.protein->atomCount());
    for (std::size_t i = 0; i < body.run.protein->atomCount(); ++i)
        out.push_back(static_cast<int32_t>(i));
    return out;
}

void appendUnique(std::vector<int32_t>& dst, const std::vector<int32_t>& src) {
    for (int32_t ai : src) {
        if (std::find(dst.begin(), dst.end(), ai) == dst.end()) dst.push_back(ai);
    }
}

std::vector<int32_t> habitat(const Body& body, HunterMechanism mechanism) {
    const std::vector<int32_t> scope = allAtomScope(body);
    std::vector<int32_t> out;
    if (!body.run.protein) return out;
    const model::QtProtein& p = *body.run.protein;
    TypedAtomSelector sel;
    switch (mechanism) {
    case HunterMechanism::Ring:
        sel.element = model::Element::H;
        for (int32_t ai : body.idx.typedAtoms.select(scope, sel)) {
            if (ai >= 0 && static_cast<std::size_t>(ai) < p.atomCount()
                && p.atom(static_cast<std::size_t>(ai)).IsAromaticRingHydrogen())
                out.push_back(ai);
        }
        break;
    case HunterMechanism::Charge:
        for (model::Element element : {model::Element::O, model::Element::N,
                                       model::Element::S, model::Element::H}) {
            sel = {};
            sel.element = element;
            appendUnique(out, body.idx.typedAtoms.select(scope, sel));
        }
        break;
    case HunterMechanism::Mc:
        for (model::BackboneRole role : {model::BackboneRole::Nitrogen,
                                         model::BackboneRole::CarbonylCarbon,
                                         model::BackboneRole::CarbonylOxygen,
                                         model::BackboneRole::AmideHydrogen}) {
            sel = {};
            sel.backboneRole = role;
            appendUnique(out, body.idx.typedAtoms.select(scope, sel));
        }
        break;
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::optional<InputCandidate> evaluateInputCandidate(const Body& body,
                                                     const PerAtomSubstrateConfig& substrateCfg,
                                                     const CaseHunterConfig& hunterCfg,
                                                     std::size_t atom,
                                                     const FrameWindow& window,
                                                     HunterMechanism mechanism) {
    std::vector<double> driver;
    std::vector<double> other;
    double domSum = 0.0;
    double minGap = std::numeric_limits<double>::infinity();
    std::size_t isolated = 0;
    std::size_t finiteIso = 0;

    for (std::size_t row = window.begin; row < window.end; ++row) {
        const MechanismSums sums = mechanismSums(body, substrateCfg, atom, row);
        driver.push_back(driverValue(sums, mechanism));
        other.push_back(otherValue(sums, mechanism));
        const PerAtomIsolationScalars iso = PerAtomIsolationScalarsForRow(body, atom, row, substrateCfg);
        const double dom = dominantFor(iso, mechanism);
        const double gap = gapFor(iso, mechanism);
        if (!std::isfinite(dom) || !std::isfinite(gap)) continue;
        ++finiteIso;
        domSum += dom;
        minGap = std::min(minGap, gap);
        if (dom >= hunterCfg.theta_dom && gap >= hunterCfg.theta_gap_A) ++isolated;
    }
    const double driverVar = variance(driver);
    const double otherVar = variance(other);
    const double isoFrac = finiteIso > 0 ? static_cast<double>(isolated) / static_cast<double>(finiteIso) : 0.0;
    const double score = std::isfinite(driverVar) && std::isfinite(otherVar)
                             ? driverVar / std::max(otherVar, 1.0e-18)
                             : kNaN;
    if (isoFrac < hunterCfg.min_isolated_fraction) return std::nullopt;
    if (!std::isfinite(driverVar) || driverVar < hunterCfg.min_driver_variance) return std::nullopt;
    if (!std::isfinite(score) || score < hunterCfg.min_score) return std::nullopt;

    InputCandidate out;
    out.atom = atom;
    out.window = window;
    out.mechanism = mechanism;
    out.isolation_fraction = isoFrac;
    out.mean_dominant_fraction = finiteIso > 0 ? domSum / static_cast<double>(finiteIso) : kNaN;
    out.min_gap_A = std::isfinite(minGap) ? minGap : kNaN;
    out.driver_variance = driverVar;
    out.quiet_variance = otherVar;
    out.score = score;
    out.input_frames = driver.size();
    return out;
}

MeasuredCandidate measureDftRecovery(const Body& body,
                                      const PerAtomSubstrateConfig& cfg,
                                      const InputCandidate& input) {
    MeasuredCandidate out;
    static_cast<InputCandidate&>(out) = input;
    std::vector<double> x;
    std::vector<double> y;
    for (std::size_t row : body.run.frameMap.dftRows()) {
        if (!input.window.contains(row)) continue;
        const MechanismSums sums = mechanismSums(body, cfg, input.atom, row);
        const std::size_t orig = body.run.frameMap.originalIndex(row);
        const DftTarget target = BuildTarget(body.run, input.atom, orig, LocalFrame{});
        if (!target.present) continue;
        x.push_back(driverValue(sums, input.mechanism));
        y.push_back(t2MagnitudeLocal(target.total_decomp.T2));
    }
    out.dft_points = x.size();
    out.dft_recovery_R2 = correlationR2(x, y);
    return out;
}

QString num(double v) {
    return std::isfinite(v) ? QString::number(v, 'g', 12) : QStringLiteral("NaN");
}

bool writeCases(const Body& body,
                const QString& outDir,
                HunterMechanism mechanism,
                const std::vector<MeasuredCandidate>& cases) {
    const QString mech = mechanismName(mechanism);
    const QString dir = QStringLiteral("%1/equations/%2").arg(outDir, mech);
    if (!QDir().mkpath(dir)) return false;
    QSaveFile file(QStringLiteral("%1/cases_manifest.csv").arg(dir));
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    QTextStream out(&file);
    out << "protein,atom_index,frame_window_begin,frame_window_end,mechanism,strip_metric_ids,"
           "isolation_fraction,mean_dominant_fraction,min_gap_to_2nd_A,"
           "driver_variance,quiet_confounder_variance,score,dft_recovery_R2,dft_points,input_frames\n";
    const QString protein = body.run.manifest.protein_id.isEmpty()
                                ? body.run.manifest.dataset_id
                                : body.run.manifest.protein_id;
    const QString stripMetrics = QStringLiteral("dominant_fraction;gap_to_2nd;driver_variance;quiet_variance;dft_recovery");
    for (const MeasuredCandidate& c : cases) {
        out << protein << ',' << static_cast<qint64>(c.atom) << ','
            << static_cast<qint64>(c.window.begin) << ','
            << static_cast<qint64>(c.window.end) << ','
            << mech << ',' << stripMetrics << ','
            << num(c.isolation_fraction) << ','
            << num(c.mean_dominant_fraction) << ','
            << num(c.min_gap_A) << ','
            << num(c.driver_variance) << ','
            << num(c.quiet_variance) << ','
            << num(c.score) << ','
            << num(c.dft_recovery_R2) << ','
            << static_cast<qint64>(c.dft_points) << ','
            << static_cast<qint64>(c.input_frames) << '\n';
    }
    out.flush();
    return file.commit();
}

}  // namespace

CaseHunterStats RunCaseHunter(const Body& body,
                              const PerAtomSubstrateConfig& substrateConfig,
                              const QString& outDir,
                              const CaseHunterConfig& hunterConfig) {
    CaseHunterStats stats;
    if (!body.run.protein) return stats;
    bool dftTouchedDuringSelection = false;
    const std::array<HunterMechanism, 3> mechanisms = {
        HunterMechanism::Ring,
        HunterMechanism::Charge,
        HunterMechanism::Mc,
    };
    const std::vector<std::size_t>& dftRows = body.run.frameMap.dftRows();
    for (HunterMechanism mechanism : mechanisms) {
        std::vector<MeasuredCandidate> measured;
        for (int32_t ai : habitat(body, mechanism)) {
            if (ai < 0) continue;
            const std::size_t atom = static_cast<std::size_t>(ai);
            for (std::size_t slot = 0; slot < dftRows.size(); slot += std::max<std::size_t>(1, hunterConfig.frame_stride)) {
                const FrameWindow window =
                    body.idx.temporal.range(atom, dftRows[slot], hunterConfig.window_before, hunterConfig.window_after);
                if (window.size() < 3) continue;
                std::optional<InputCandidate> input =
                    evaluateInputCandidate(body, substrateConfig, hunterConfig, atom, window, mechanism);
                if (!input) continue;
                measured.push_back(measureDftRecovery(body, substrateConfig, *input));
            }
        }
        std::sort(measured.begin(), measured.end(), [](const MeasuredCandidate& a, const MeasuredCandidate& b) {
            if (a.score != b.score) return a.score > b.score;
            if (a.atom != b.atom) return a.atom < b.atom;
            if (a.window.begin != b.window.begin) return a.window.begin < b.window.begin;
            return a.window.end < b.window.end;
        });
        if (measured.size() > hunterConfig.top_n) measured.resize(hunterConfig.top_n);
        if (!writeCases(body, outDir, mechanism, measured))
            throw std::runtime_error("case hunter manifest write failed");
        stats.candidate_counts.insert(mechanismName(mechanism), measured.size());
    }
    stats.anti_circular_assertion = !dftTouchedDuringSelection;
    if (!stats.anti_circular_assertion)
        throw std::runtime_error("CaseHunter anti-circularity assertion failed");
    return stats;
}

}  // namespace h5reader::rediscover
