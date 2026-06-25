// DistributionSummary -- compact distribution-preserving summaries for
// rediscover emits.

#pragma once

#include <QJsonObject>
#include <QString>

#include <cstddef>
#include <vector>

namespace h5reader::rediscover {

struct DistributionSummary {
    std::size_t n = 0;
    std::size_t finite_n = 0;
    double finite_frac = 0.0;
    double mean = 0.0;
    double sd = 0.0;
    double min = 0.0;
    double p05 = 0.0;
    double p25 = 0.0;
    double median = 0.0;
    double p75 = 0.0;
    double p95 = 0.0;
    double max = 0.0;
    std::vector<std::size_t> quantile_bins;

    bool hasFinite() const { return finite_n > 0; }
};

DistributionSummary SummarizeDistribution(const std::vector<double>& values,
                                          std::size_t binCount = 5);

QJsonObject DistributionSummaryJson(const DistributionSummary& summary);

}  // namespace h5reader::rediscover
