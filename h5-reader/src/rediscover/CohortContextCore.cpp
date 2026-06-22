#include "CohortContextAccumulator.h"

#include <QRegularExpression>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>

namespace h5reader::rediscover {
namespace {

constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

bool finite(double v) { return std::isfinite(v); }

QString cleanToken(QString v, const QString& fallback = QStringLiteral("unknown")) {
    v = v.trimmed();
    if (v.isEmpty()) v = fallback;
    v.replace(QLatin1Char('|'), QLatin1Char('_'));
    v.replace(QLatin1Char('='), QLatin1Char('_'));
    v.replace(QLatin1Char(','), QLatin1Char('_'));
    v.replace(QLatin1Char('\n'), QLatin1Char('_'));
    v.replace(QRegularExpression(QStringLiteral("\\s+")), QStringLiteral("_"));
    return v.isEmpty() ? fallback : v;
}

}  // namespace

Axis2ContextKey BuildAxis2ContextKey(const Axis2ContextKeyFields& input) {
    Axis2ContextKey out;
    out.fields = input;
    out.fields.element = cleanToken(out.fields.element);
    out.fields.residue_type = cleanToken(out.fields.residue_type);
    out.fields.atom_name = cleanToken(out.fields.atom_name);
    out.fields.hyb = cleanToken(out.fields.hyb);
    out.fields.contact_class = cleanToken(out.fields.contact_class);
    out.fields.dihedral_region = cleanToken(out.fields.dihedral_region);
    out.fields.SS = cleanToken(out.fields.SS);
    out.identity = QStringLiteral("element=%1;residue_type=%2;atom_name=%3;hyb=%4")
                       .arg(out.fields.element, out.fields.residue_type,
                            out.fields.atom_name, out.fields.hyb);
    out.context = QStringLiteral("contact_class=%1;dihedral_region=%2;SS=%3")
                      .arg(out.fields.contact_class, out.fields.dihedral_region,
                           out.fields.SS);
    out.canonical = QStringLiteral("schema=%1|identity.%2|context.%3")
                        .arg(out.fields.schema_version)
                        .arg(out.identity, out.context);
    return out;
}

SupportCredential CredentialSupport(std::size_t n_proteins,
                                    bool full_rank,
                                    SupportThresholds thresholds) {
    SupportCredential out;
    if (n_proteins >= thresholds.n_full && full_rank) {
        out.support_class = SupportClass::Full;
        out.support_name = QStringLiteral("full");
        out.may_emit_coupling = true;
        out.may_emit_full_subspace = true;
        out.underpowered_dimensions = QStringLiteral("none");
    } else if (n_proteins >= thresholds.n_min) {
        out.support_class = SupportClass::Reduced;
        out.support_name = QStringLiteral("reduced");
        out.may_emit_coupling = true;
        out.may_emit_full_subspace = false;
        out.underpowered_dimensions = full_rank ? QStringLiteral("n<N_full")
                                                : QStringLiteral("rank_deficient");
    } else {
        out.support_class = SupportClass::Insufficient;
        out.support_name = QStringLiteral("insufficient");
        out.may_emit_coupling = false;
        out.may_emit_full_subspace = false;
        out.underpowered_dimensions = QStringLiteral("n<N_min");
    }
    return out;
}

double LinearSlope(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return kNan;
    double sx = 0.0;
    double sy = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (finite(x[i]) && finite(y[i])) {
            sx += x[i];
            sy += y[i];
            ++n;
        }
    }
    if (n < 2) return kNan;
    const double mx = sx / static_cast<double>(n);
    const double my = sy / static_cast<double>(n);
    double sxx = 0.0;
    double sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!finite(x[i]) || !finite(y[i])) continue;
        const double dx = x[i] - mx;
        sxx += dx * dx;
        sxy += dx * (y[i] - my);
    }
    return sxx > 0.0 ? sxy / sxx : kNan;
}

double PearsonR(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return kNan;
    double sx = 0.0;
    double sy = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (finite(x[i]) && finite(y[i])) {
            sx += x[i];
            sy += y[i];
            ++n;
        }
    }
    if (n < 2) return kNan;
    const double mx = sx / static_cast<double>(n);
    const double my = sy / static_cast<double>(n);
    double sxx = 0.0;
    double syy = 0.0;
    double sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!finite(x[i]) || !finite(y[i])) continue;
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    const double denom = std::sqrt(sxx * syy);
    return denom > 0.0 ? sxy / denom : kNan;
}

PermutationNull ProteinLabelPermutationNull(const std::vector<double>& driver,
                                            const std::vector<double>& sigma,
                                            int k,
                                            quint32 seed) {
    PermutationNull out;
    out.permutation_K = k;
    const double obs = LinearSlope(driver, sigma);
    if (!finite(obs) || driver.size() != sigma.size() || driver.size() < 3 || k <= 0)
        return out;
    std::vector<double> shuffled = sigma;
    std::vector<double> slopes;
    slopes.reserve(static_cast<std::size_t>(k));
    std::mt19937 rng(seed);
    int extreme = 0;
    for (int i = 0; i < k; ++i) {
        std::shuffle(shuffled.begin(), shuffled.end(), rng);
        const double s = LinearSlope(driver, shuffled);
        if (finite(s)) {
            slopes.push_back(s);
            if (std::abs(s) >= std::abs(obs)) ++extreme;
        }
    }
    if (slopes.empty()) return out;
    out.null_slope_mean = std::accumulate(slopes.begin(), slopes.end(), 0.0)
                          / static_cast<double>(slopes.size());
    double ss = 0.0;
    for (double s : slopes) {
        const double d = s - out.null_slope_mean;
        ss += d * d;
    }
    out.null_slope_sd = slopes.size() > 1
                            ? std::sqrt(ss / static_cast<double>(slopes.size() - 1))
                            : 0.0;
    out.obs_slope_z = out.null_slope_sd > 0.0
                          ? (obs - out.null_slope_mean) / out.null_slope_sd
                          : 0.0;
    out.perm_p = static_cast<double>(extreme + 1) / static_cast<double>(slopes.size() + 1);
    return out;
}

double ComputeHelixDipoleField(const HelixDipoleInput& input) {
    double field = 0.0;
    for (double z : input.ca_z_A) {
        const double dz = input.target_z_A - z;
        const double r2 = dz * dz + 4.0;
        const double sign = dz >= 0.0 ? 1.0 : -1.0;
        field += input.charge_per_residue * sign / r2;
    }
    return field;
}

void CohortContextAccumulator::push(const CohortSample& sample) {
    CohortCellTruth& cell = cells_[sample.key.canonical];
    if (cell.key.canonical.isEmpty())
        cell.key = sample.key;
    cell.proteins.insert(sample.protein_id);
    cell.samples.push_back(sample);
    ++sample_count_;
}

}  // namespace h5reader::rediscover
