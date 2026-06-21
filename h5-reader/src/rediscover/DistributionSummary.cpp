#include "DistributionSummary.h"

#include <QJsonArray>
#include <QJsonValue>

#include <algorithm>
#include <cmath>
#include <limits>

namespace h5reader::rediscover {
namespace {

bool finite(double v) { return std::isfinite(v); }

QJsonValue jd(double v) {
    return finite(v) ? QJsonValue(v) : QJsonValue(QJsonValue::Null);
}

double percentile(const std::vector<double>& sorted, double p) {
    if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (sorted.size() == 1) return sorted.front();
    const double x = p * static_cast<double>(sorted.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(x));
    const auto hi = static_cast<std::size_t>(std::ceil(x));
    const double frac = x - static_cast<double>(lo);
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

}  // namespace

DistributionSummary SummarizeDistribution(const std::vector<double>& values,
                                          std::size_t binCount) {
    DistributionSummary s;
    s.n = values.size();

    std::vector<double> finiteValues;
    finiteValues.reserve(values.size());
    for (double v : values) {
        if (finite(v)) finiteValues.push_back(v);
    }
    s.finite_n = finiteValues.size();
    s.finite_frac = s.n > 0 ? static_cast<double>(s.finite_n) / static_cast<double>(s.n) : 0.0;
    if (finiteValues.empty()) return s;

    std::sort(finiteValues.begin(), finiteValues.end());
    s.min = finiteValues.front();
    s.p05 = percentile(finiteValues, 0.05);
    s.p25 = percentile(finiteValues, 0.25);
    s.median = percentile(finiteValues, 0.50);
    s.p75 = percentile(finiteValues, 0.75);
    s.p95 = percentile(finiteValues, 0.95);
    s.max = finiteValues.back();

    double sum = 0.0;
    for (double v : finiteValues) sum += v;
    s.mean = sum / static_cast<double>(finiteValues.size());
    if (finiteValues.size() > 1) {
        double ss = 0.0;
        for (double v : finiteValues) {
            const double d = v - s.mean;
            ss += d * d;
        }
        s.sd = std::sqrt(ss / static_cast<double>(finiteValues.size() - 1));
    }

    if (binCount > 0) {
        s.quantile_bins.assign(binCount, 0);
        if (s.max == s.min) {
            s.quantile_bins.front() = finiteValues.size();
        } else {
            for (double v : finiteValues) {
                std::size_t b = static_cast<std::size_t>(
                    std::floor((v - s.min) / (s.max - s.min) * static_cast<double>(binCount)));
                if (b >= binCount) b = binCount - 1;
                ++s.quantile_bins[b];
            }
        }
    }
    return s;
}

QJsonObject DistributionSummaryJson(const DistributionSummary& s) {
    QJsonObject o;
    o.insert(QStringLiteral("n"), static_cast<qint64>(s.n));
    o.insert(QStringLiteral("finite_n"), static_cast<qint64>(s.finite_n));
    o.insert(QStringLiteral("finite_frac"), s.finite_frac);
    o.insert(QStringLiteral("mean"), s.hasFinite() ? jd(s.mean) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("sd"), s.hasFinite() ? jd(s.sd) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("min"), s.hasFinite() ? jd(s.min) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("p05"), s.hasFinite() ? jd(s.p05) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("p25"), s.hasFinite() ? jd(s.p25) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("median"), s.hasFinite() ? jd(s.median) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("p75"), s.hasFinite() ? jd(s.p75) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("p95"), s.hasFinite() ? jd(s.p95) : QJsonValue(QJsonValue::Null));
    o.insert(QStringLiteral("max"), s.hasFinite() ? jd(s.max) : QJsonValue(QJsonValue::Null));
    QJsonArray bins;
    for (std::size_t v : s.quantile_bins) bins.append(static_cast<qint64>(v));
    o.insert(QStringLiteral("quantile_bins"), bins);
    return o;
}

}  // namespace h5reader::rediscover
