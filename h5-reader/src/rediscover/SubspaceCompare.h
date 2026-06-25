#pragma once

#include <QJsonObject>
#include <QString>
#include <QStringList>

#include <cstddef>
#include <limits>
#include <vector>

namespace h5reader::rediscover {

struct SubspaceChannel {
    QString name;
    std::vector<double> values;
};

struct SubspaceFamily {
    QString name;
    std::vector<SubspaceChannel> channels;
};

struct SubspaceCompareResult {
    bool computed = false;
    QString provenance = QStringLiteral("svd_subspace_compare_v1");
    QString missing_reason;
    QString overlap_label;
    double variance_threshold = 0.90;

    int finite_n = 0;
    int input_dim_a = 0;
    int input_dim_b = 0;
    int active_dim_a = 0;
    int active_dim_b = 0;
    int basis_dim_a = 0;
    int basis_dim_b = 0;
    double explained_fraction_a = std::numeric_limits<double>::quiet_NaN();
    double explained_fraction_b = std::numeric_limits<double>::quiet_NaN();
    double condition_number_a = std::numeric_limits<double>::quiet_NaN();
    double condition_number_b = std::numeric_limits<double>::quiet_NaN();

    double max_canonical_corr = std::numeric_limits<double>::quiet_NaN();
    double mean_canonical_corr = std::numeric_limits<double>::quiet_NaN();
    int n_cc_ge_0_80 = 0;
    int n_cc_ge_0_95 = 0;
    double min_angle_deg = std::numeric_limits<double>::quiet_NaN();

    QStringList dropped_channels_a;
    QStringList dropped_channels_b;
    std::vector<double> canonical_corrs;
    std::vector<double> principal_angles_deg;
    std::vector<double> singular_values_a;
    std::vector<double> singular_values_b;
    std::vector<double> explained_spectrum_a;
    std::vector<double> explained_spectrum_b;
};

SubspaceCompareResult CompareSubspaces(const SubspaceFamily& a,
                                       const SubspaceFamily& b,
                                       const std::vector<std::size_t>& rows,
                                       double varianceThreshold = 0.90);

QJsonObject SubspaceCompareJson(const SubspaceCompareResult& r);

}  // namespace h5reader::rediscover
