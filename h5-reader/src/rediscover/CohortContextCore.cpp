#include "CohortContextAccumulator.h"

#include <QRegularExpression>

#include <Eigen/Dense>

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

double percentile(const std::vector<double>& sorted, double p) {
    if (sorted.empty()) return kNan;
    if (sorted.size() == 1) return sorted.front();
    const double x = p * static_cast<double>(sorted.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(x));
    const auto hi = static_cast<std::size_t>(std::ceil(x));
    const double frac = x - static_cast<double>(lo);
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

std::array<double, 6> symmetricComponents(const model::Mat3& raw) {
    const model::Mat3 sym = 0.5 * (raw + raw.transpose());
    return {sym(0, 0), sym(1, 1), sym(0, 1), sym(0, 2), sym(1, 2), sym(2, 2)};
}

double iso(const model::Mat3& m) {
    if (!m.allFinite()) return kNan;
    const model::Mat3 sym = 0.5 * (m + m.transpose());
    return (sym(0, 0) + sym(1, 1) + sym(2, 2)) / 3.0;
}

double etaH(const model::Mat3& raw) {
    if (!raw.allFinite()) return kNan;
    const model::Mat3 sym = 0.5 * (raw + raw.transpose());
    Eigen::SelfAdjointEigenSolver<model::Mat3> es(sym);
    if (es.info() != Eigen::Success) return kNan;
    struct Eval {
        double value = 0.0;
        double distance = 0.0;
    };
    const double isotropic = es.eigenvalues().mean();
    std::array<Eval, 3> e{};
    for (int i = 0; i < 3; ++i)
        e[static_cast<std::size_t>(i)] = {es.eigenvalues()(i), std::abs(es.eigenvalues()(i) - isotropic)};
    std::sort(e.begin(), e.end(), [](const Eval& a, const Eval& b) {
        return a.distance < b.distance;
    });
    const double xx = e[0].value;
    const double yy = e[1].value;
    const double zz = e[2].value;
    const double denom = zz - isotropic;
    if (!finite(denom) || std::abs(denom) < 1.0e-12) return kNan;
    return (yy - xx) / denom;
}

}  // namespace

Axis2ContextKey BuildAxis2ContextKey(const Axis2ContextKeyFields& input) {
    Axis2ContextKey out;
    out.fields = input;
    out.fields.element = cleanToken(out.fields.element);
    out.fields.residue_type = cleanToken(out.fields.residue_type);
    out.fields.atom_name = cleanToken(out.fields.atom_name);
    out.fields.frame_kind = cleanToken(out.fields.frame_kind, QStringLiteral("none"));
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

DistantRidgeCharacterization CharacterizeDistantNonzeroRidge(std::size_t distantCount,
                                                            const QString& anySiteScope,
                                                            double slope,
                                                            const QString& channel) {
    DistantRidgeCharacterization out;
    out.flagged = distantCount > 0
                  && finite(slope)
                  && std::abs(slope) > 0.05
                  && anySiteScope == QStringLiteral("distant_from_all_sites");
    if (out.flagged) {
        out.distant_zero_check = QStringLiteral("flagged_nonzero_distant_from_all_sites");
        out.characterization = QStringLiteral("characterized_not_gated");
        out.nonzero_channel = channel;
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

void BoundedDistributionAccumulator::add(double v) {
    ++n;
    if (!finite(v)) return;
    ++finite_n;
    if (finite_n == 1) {
        mean = v;
        m2 = 0.0;
        min = v;
        max = v;
    } else {
        const double delta = v - mean;
        mean += delta / static_cast<double>(finite_n);
        m2 += delta * (v - mean);
        min = std::min(min, v);
        max = std::max(max, v);
    }
    if (reservoir.size() < kReservoirLimit) {
        reservoir.push_back(v);
    } else {
        const std::size_t idx = (finite_n * 1103515245ull + 12345ull) % kReservoirLimit;
        reservoir[idx] = v;
    }
}

DistributionSummary BoundedDistributionAccumulator::summary(std::size_t binCount) const {
    DistributionSummary s;
    s.n = n;
    s.finite_n = finite_n;
    s.finite_frac = n > 0 ? static_cast<double>(finite_n) / static_cast<double>(n) : 0.0;
    if (finite_n == 0 || reservoir.empty()) return s;

    std::vector<double> sorted = reservoir;
    std::sort(sorted.begin(), sorted.end());
    s.mean = mean;
    s.sd = finite_n > 1 ? std::sqrt(m2 / static_cast<double>(finite_n - 1)) : 0.0;
    s.min = min;
    s.p05 = percentile(sorted, 0.05);
    s.p25 = percentile(sorted, 0.25);
    s.median = percentile(sorted, 0.50);
    s.p75 = percentile(sorted, 0.75);
    s.p95 = percentile(sorted, 0.95);
    s.max = max;
    if (binCount > 0) {
        s.quantile_bins.assign(binCount, 0);
        if (max == min) {
            s.quantile_bins.front() = finite_n;
        } else {
            for (double v : sorted) {
                std::size_t b = static_cast<std::size_t>(
                    std::floor((v - min) / (max - min) * static_cast<double>(binCount)));
                if (b >= binCount) b = binCount - 1;
                ++s.quantile_bins[b];
            }
        }
    }
    return s;
}

void RunningMeanAccumulator::add(double v) {
    if (!finite(v)) return;
    ++n;
    sum += v;
}

double RunningMeanAccumulator::meanValue() const {
    return n > 0 ? sum / static_cast<double>(n) : kNan;
}

void PairAccumulator::add(double x, double y) {
    if (!finite(x) || !finite(y)) return;
    ++n;
    sx += x;
    sy += y;
    sxx += x * x;
    syy += y * y;
    sxy += x * y;
}

double PairAccumulator::slope() const {
    if (n < 2) return kNan;
    const double nn = static_cast<double>(n);
    const double den = sxx - sx * sx / nn;
    if (!(den > 0.0)) return kNan;
    return (sxy - sx * sy / nn) / den;
}

double PairAccumulator::pearson() const {
    if (n < 2) return kNan;
    const double nn = static_cast<double>(n);
    const double vx = sxx - sx * sx / nn;
    const double vy = syy - sy * sy / nn;
    const double cov = sxy - sx * sy / nn;
    const double den = std::sqrt(vx * vy);
    return den > 0.0 ? cov / den : kNan;
}

ClassicalAgreementStats ComputeClassicalAgreementForCell(const CohortCellTruth& cell,
                                                         double buckinghamA) {
    std::vector<double> cl;
    std::vector<double> qm;
    std::vector<double> residuals;
    cl.reserve(static_cast<std::size_t>(cell.protein_folds.size()));
    qm.reserve(static_cast<std::size_t>(cell.protein_folds.size()));
    residuals.reserve(static_cast<std::size_t>(cell.protein_folds.size()));
    double ss = 0.0;
    for (auto it = cell.protein_folds.begin(); it != cell.protein_folds.end(); ++it) {
        const CohortProteinFold& p = it.value();
        const double sigma = p.sigma.meanValue();
        const double apbs = p.channels.value(QStringLiteral("apbs_E_mag")).meanValue();
        const double ring = p.channels.value(QStringLiteral("ring_bs_iso")).meanValue();
        const double mc = p.channels.value(QStringLiteral("mc_lit_iso")).meanValue();
        const double classical = buckinghamA * apbs + ring + mc;
        if (!finite(sigma) || !finite(classical)) continue;
        cl.push_back(classical);
        qm.push_back(sigma);
        const double residual = sigma - classical;
        residuals.push_back(residual);
        ss += residual * residual;
    }
    ClassicalAgreementStats out;
    out.r = PearsonR(cl, qm);
    out.slope = LinearSlope(cl, qm);
    if (!residuals.empty()) {
        out.rmsd = std::sqrt(ss / static_cast<double>(residuals.size()));
        const DistributionSummary res = SummarizeDistribution(residuals);
        out.residual_mean = res.mean;
        out.residual_sd = res.sd;
    }
    return out;
}

Axis2FoldedTensor FoldAxis2TensorChannels(const model::Mat3& raw,
                                          const std::optional<model::Mat3>& molecularAxes) {
    Axis2FoldedTensor out;
    out.sigma_iso = iso(raw);
    out.sigma_eta_H = etaH(raw);
    if (molecularAxes && molecularAxes->allFinite()) {
        const model::Mat3 local = molecularAxes->transpose() * raw * (*molecularAxes);
        if (local.allFinite()) {
            out.mol_components = symmetricComponents(local);
            out.molecular_frame_projected = true;
            out.projection = QStringLiteral("molecular_frame_projected:R^T*T*R");
        }
    }
    if (!out.molecular_frame_projected) {
        out.projection = QStringLiteral("molecular_frame_unavailable_no_raw_lab_fallback");
    }
    return out;
}

Axis2FoldedTensor FoldAxis2TensorChannels(const model::Mat3& dia,
                                          const model::Mat3& para,
                                          const std::optional<model::Mat3>& molecularAxes) {
    return FoldAxis2TensorChannels(dia + para, molecularAxes);
}

std::size_t CohortCellTruth::retainedAccumulatorValueCount() const {
    std::size_t n = sigma.retainedValueCount()
                    + mol_xx.retainedValueCount()
                    + mol_yy.retainedValueCount()
                    + mol_xy.retainedValueCount()
                    + mol_xz.retainedValueCount()
                    + mol_yz.retainedValueCount()
                    + mol_zz.retainedValueCount()
                    + eta_H.retainedValueCount()
                    + helix_dipole_field.retainedValueCount();
    for (const auto& item : channel_distributions)
        n += item.retainedValueCount();
    return n;
}

void CohortContextAccumulator::push(const CohortSample& sample) {
    CohortCellTruth& cell = cells_[sample.key.canonical];
    if (cell.key.canonical.isEmpty())
        cell.key = sample.key;
    cell.proteins.insert(sample.protein_id);
    ++cell.sample_count;
    cell.sigma.add(sample.sigma_iso);
    cell.mol_xx.add(sample.mol_xx);
    cell.mol_yy.add(sample.mol_yy);
    cell.mol_xy.add(sample.mol_xy);
    cell.mol_xz.add(sample.mol_xz);
    cell.mol_yz.add(sample.mol_yz);
    cell.mol_zz.add(sample.mol_zz);
    cell.eta_H.add(sample.sigma_eta_H);
    cell.helix_dipole_field.add(sample.helix_dipole_field);
    cell.psi_iminus1_vs_sigma.add(sample.psi_iminus1, sample.sigma_iso);
    cell.psi_own_vs_sigma.add(sample.psi_own, sample.sigma_iso);

    if (sample.backbone_n && cell.psi_iminus1_region == QStringLiteral("not_backbone_N")) {
        cell.psi_iminus1_region = sample.psi_iminus1_region.isEmpty()
                                      ? QStringLiteral("unknown")
                                      : sample.psi_iminus1_region;
        cell.predecessor_identity = sample.predecessor_identity.isEmpty()
                                        ? QStringLiteral("unknown")
                                        : sample.predecessor_identity;
    }

    CohortProteinFold& protein = cell.protein_folds[sample.protein_id];
    protein.sigma.add(sample.sigma_iso);
    for (auto it = sample.channels.begin(); it != sample.channels.end(); ++it) {
        cell.channel_distributions[it.key()].add(it.value());
        cell.channel_vs_sigma[it.key()].add(it.value(), sample.sigma_iso);
        protein.channels[it.key()].add(it.value());
    }
    ++sample_count_;
}

}  // namespace h5reader::rediscover
