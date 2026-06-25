#include "SubspaceCompare.h"

#include <Eigen/Dense>

#include <QJsonArray>
#include <QJsonValue>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace h5reader::rediscover {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kScaleFloor = 1e-12;

bool finite(double v) { return std::isfinite(v); }

QJsonValue jd(double v) {
    return finite(v) ? QJsonValue(v) : QJsonValue(QJsonValue::Null);
}

QJsonArray vectorJson(const std::vector<double>& values) {
    QJsonArray out;
    for (double v : values) out.append(jd(v));
    return out;
}

QJsonArray stringListJson(const QStringList& values) {
    QJsonArray out;
    for (const QString& v : values) out.append(v);
    return out;
}

struct PreparedMatrix {
    Eigen::MatrixXd x;
    QStringList dropped;
    int input_dim = 0;
    int active_dim = 0;
};

PreparedMatrix prepare(const SubspaceFamily& f,
                       const std::vector<std::size_t>& finiteRows) {
    PreparedMatrix out;
    out.input_dim = static_cast<int>(f.channels.size());
    const int n = static_cast<int>(finiteRows.size());
    std::vector<Eigen::VectorXd> columns;
    columns.reserve(f.channels.size());

    for (const SubspaceChannel& ch : f.channels) {
        Eigen::VectorXd raw(n);
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            const double v = ch.values[finiteRows[static_cast<std::size_t>(i)]];
            raw(i) = v;
            sum += v;
        }
        const double mean = n > 0 ? sum / static_cast<double>(n) : 0.0;
        double ss = 0.0;
        double maxAbs = 0.0;
        for (int i = 0; i < n; ++i) {
            raw(i) -= mean;
            ss += raw(i) * raw(i);
            maxAbs = std::max(maxAbs, std::abs(raw(i) + mean));
        }
        const double denom = n > 1 ? std::sqrt(ss / static_cast<double>(n - 1)) : 0.0;
        const double floor = kScaleFloor * std::max(1.0, maxAbs);
        if (!(denom > floor) || !finite(denom)) {
            out.dropped.append(ch.name);
            continue;
        }
        raw /= denom;
        columns.push_back(raw);
    }

    out.active_dim = static_cast<int>(columns.size());
    out.x = Eigen::MatrixXd(n, out.active_dim);
    for (int c = 0; c < out.active_dim; ++c)
        out.x.col(c) = columns[static_cast<std::size_t>(c)];
    return out;
}

struct Basis {
    Eigen::MatrixXd q;
    std::vector<double> singular_values;
    std::vector<double> explained_spectrum;
    int basis_dim = 0;
    double explained_fraction = std::numeric_limits<double>::quiet_NaN();
    double condition_number = std::numeric_limits<double>::quiet_NaN();
};

Basis svdBasis(const Eigen::MatrixXd& x, double varianceThreshold) {
    Basis out;
    if (x.rows() < 2 || x.cols() < 1) return out;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(x, Eigen::ComputeThinU);
    const Eigen::VectorXd s = svd.singularValues();
    if (s.size() == 0) return out;

    double total = 0.0;
    for (int i = 0; i < s.size(); ++i) {
        if (!finite(s(i))) continue;
        total += s(i) * s(i);
        out.singular_values.push_back(s(i));
    }
    if (!(total > 0.0)) return out;

    double cum = 0.0;
    for (int i = 0; i < s.size(); ++i) {
        const double frac = (s(i) * s(i)) / total;
        out.explained_spectrum.push_back(frac);
        if (out.basis_dim == 0 || cum < varianceThreshold) {
            cum += frac;
            out.basis_dim = i + 1;
            out.explained_fraction = cum;
        }
    }
    if (out.basis_dim < 1) {
        out.basis_dim = 1;
        out.explained_fraction = out.explained_spectrum.empty()
                                     ? std::numeric_limits<double>::quiet_NaN()
                                     : out.explained_spectrum.front();
    }
    out.basis_dim = std::min(out.basis_dim, static_cast<int>(s.size()));
    out.q = svd.matrixU().leftCols(out.basis_dim);

    const double largest = s(0);
    double smallest = std::numeric_limits<double>::quiet_NaN();
    for (int i = s.size() - 1; i >= 0; --i) {
        if (s(i) > kScaleFloor) {
            smallest = s(i);
            break;
        }
    }
    if (finite(largest) && finite(smallest) && smallest > 0.0)
        out.condition_number = largest / smallest;
    return out;
}

QString overlapLabel(double maxCorr, double meanCorr) {
    if (!finite(maxCorr)) return QStringLiteral("not_computed");
    if (maxCorr >= 0.95 && finite(meanCorr) && meanCorr >= 0.80)
        return QStringLiteral("shared_subspace");
    if (maxCorr >= 0.80) return QStringLiteral("partial_overlap");
    return QStringLiteral("weak_overlap");
}

}  // namespace

SubspaceCompareResult CompareSubspaces(const SubspaceFamily& a,
                                       const SubspaceFamily& b,
                                       const std::vector<std::size_t>& rows,
                                       double varianceThreshold) {
    SubspaceCompareResult out;
    out.variance_threshold = varianceThreshold;
    out.input_dim_a = static_cast<int>(a.channels.size());
    out.input_dim_b = static_cast<int>(b.channels.size());
    if (a.channels.empty() || b.channels.empty()) {
        out.missing_reason = QStringLiteral("empty_family");
        return out;
    }

    auto filterPresentChannels = [&](const SubspaceFamily& f, QStringList& dropped) {
        SubspaceFamily filtered{f.name, {}};
        for (const SubspaceChannel& ch : f.channels) {
            bool any = false;
            for (std::size_t row : rows) {
                if (row < ch.values.size() && finite(ch.values[row])) {
                    any = true;
                    break;
                }
            }
            if (any) filtered.channels.push_back(ch);
            else dropped.append(QStringLiteral("%1(all_missing)").arg(ch.name));
        }
        return filtered;
    };
    QStringList allMissingA;
    QStringList allMissingB;
    const SubspaceFamily fa = filterPresentChannels(a, allMissingA);
    const SubspaceFamily fb = filterPresentChannels(b, allMissingB);
    if (fa.channels.empty() || fb.channels.empty()) {
        out.dropped_channels_a = allMissingA;
        out.dropped_channels_b = allMissingB;
        out.missing_reason = QStringLiteral("zero_present_channels");
        return out;
    }

    std::vector<std::size_t> finiteRows;
    finiteRows.reserve(rows.size());
    for (std::size_t row : rows) {
        bool ok = true;
        for (const SubspaceChannel& ch : fa.channels) {
            if (row >= ch.values.size() || !finite(ch.values[row])) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;
        for (const SubspaceChannel& ch : fb.channels) {
            if (row >= ch.values.size() || !finite(ch.values[row])) {
                ok = false;
                break;
            }
        }
        if (ok) finiteRows.push_back(row);
    }
    out.finite_n = static_cast<int>(finiteRows.size());
    if (finiteRows.size() < 3) {
        out.missing_reason = QStringLiteral("insufficient_finite_rows");
        return out;
    }

    const PreparedMatrix pa = prepare(fa, finiteRows);
    const PreparedMatrix pb = prepare(fb, finiteRows);
    out.active_dim_a = pa.active_dim;
    out.active_dim_b = pb.active_dim;
    out.dropped_channels_a = allMissingA;
    out.dropped_channels_a.append(pa.dropped);
    out.dropped_channels_b = allMissingB;
    out.dropped_channels_b.append(pb.dropped);
    if (pa.active_dim < 1 || pb.active_dim < 1) {
        out.missing_reason = QStringLiteral("zero_active_dimension_after_constant_drop");
        return out;
    }

    const Basis ba = svdBasis(pa.x, varianceThreshold);
    const Basis bb = svdBasis(pb.x, varianceThreshold);
    out.basis_dim_a = ba.basis_dim;
    out.basis_dim_b = bb.basis_dim;
    out.explained_fraction_a = ba.explained_fraction;
    out.explained_fraction_b = bb.explained_fraction;
    out.condition_number_a = ba.condition_number;
    out.condition_number_b = bb.condition_number;
    out.singular_values_a = ba.singular_values;
    out.singular_values_b = bb.singular_values;
    out.explained_spectrum_a = ba.explained_spectrum;
    out.explained_spectrum_b = bb.explained_spectrum;
    if (ba.q.cols() < 1 || bb.q.cols() < 1) {
        out.missing_reason = QStringLiteral("svd_basis_empty");
        return out;
    }

    const Eigen::MatrixXd overlap = ba.q.transpose() * bb.q;
    Eigen::JacobiSVD<Eigen::MatrixXd> overlapSvd(overlap, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::VectorXd cc = overlapSvd.singularValues();
    double sum = 0.0;
    for (int i = 0; i < cc.size(); ++i) {
        const double c = std::min(1.0, std::max(0.0, cc(i)));
        out.canonical_corrs.push_back(c);
        out.principal_angles_deg.push_back(std::acos(c) * 180.0 / kPi);
        out.max_canonical_corr = finite(out.max_canonical_corr)
                                     ? std::max(out.max_canonical_corr, c)
                                     : c;
        sum += c;
        if (c >= 0.80) ++out.n_cc_ge_0_80;
        if (c >= 0.95) ++out.n_cc_ge_0_95;
    }
    if (!out.canonical_corrs.empty()) {
        out.mean_canonical_corr = sum / static_cast<double>(out.canonical_corrs.size());
        out.min_angle_deg = *std::min_element(out.principal_angles_deg.begin(),
                                              out.principal_angles_deg.end());
    }
    out.overlap_label = overlapLabel(out.max_canonical_corr, out.mean_canonical_corr);
    out.computed = true;
    return out;
}

QJsonObject SubspaceCompareJson(const SubspaceCompareResult& r) {
    QJsonObject o;
    o.insert(QStringLiteral("provenance"), r.provenance);
    o.insert(QStringLiteral("status"), r.computed ? QStringLiteral("computed")
                                                   : QStringLiteral("missing"));
    o.insert(QStringLiteral("missing_reason"),
             r.missing_reason.isEmpty() ? QJsonValue(QJsonValue::Null)
                                        : QJsonValue(r.missing_reason));
    o.insert(QStringLiteral("variance_threshold"), jd(r.variance_threshold));
    o.insert(QStringLiteral("finite_n"), r.finite_n);
    o.insert(QStringLiteral("input_dim_a"), r.input_dim_a);
    o.insert(QStringLiteral("input_dim_b"), r.input_dim_b);
    o.insert(QStringLiteral("active_dim_a"), r.active_dim_a);
    o.insert(QStringLiteral("active_dim_b"), r.active_dim_b);
    o.insert(QStringLiteral("basis_dim_a"), r.computed ? QJsonValue(r.basis_dim_a)
                                                       : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("basis_dim_b"), r.computed ? QJsonValue(r.basis_dim_b)
                                                       : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("explained_fraction_a"), jd(r.explained_fraction_a));
    o.insert(QStringLiteral("explained_fraction_b"), jd(r.explained_fraction_b));
    o.insert(QStringLiteral("condition_number_a"), jd(r.condition_number_a));
    o.insert(QStringLiteral("condition_number_b"), jd(r.condition_number_b));
    o.insert(QStringLiteral("max_canonical_corr"), jd(r.max_canonical_corr));
    o.insert(QStringLiteral("mean_canonical_corr"), jd(r.mean_canonical_corr));
    o.insert(QStringLiteral("n_cc_ge_0_80"), r.n_cc_ge_0_80);
    o.insert(QStringLiteral("n_cc_ge_0_95"), r.n_cc_ge_0_95);
    o.insert(QStringLiteral("min_angle_deg"), jd(r.min_angle_deg));
    o.insert(QStringLiteral("canonical_corrs"), vectorJson(r.canonical_corrs));
    o.insert(QStringLiteral("principal_angles_deg"), vectorJson(r.principal_angles_deg));
    o.insert(QStringLiteral("singular_values_a"), vectorJson(r.singular_values_a));
    o.insert(QStringLiteral("singular_values_b"), vectorJson(r.singular_values_b));
    o.insert(QStringLiteral("explained_spectrum_a"), vectorJson(r.explained_spectrum_a));
    o.insert(QStringLiteral("explained_spectrum_b"), vectorJson(r.explained_spectrum_b));
    o.insert(QStringLiteral("dropped_channels_a"), stringListJson(r.dropped_channels_a));
    o.insert(QStringLiteral("dropped_channels_b"), stringListJson(r.dropped_channels_b));
    o.insert(QStringLiteral("overlap_label"),
             r.overlap_label.isEmpty() ? QJsonValue(QJsonValue::Null)
                                       : QJsonValue(r.overlap_label));
    return o;
}

}  // namespace h5reader::rediscover
