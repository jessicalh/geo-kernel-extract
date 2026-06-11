#pragma once

#include <QString>

namespace h5reader::rediscover {

struct SpineProbeConfig {
    QString root720;
    QString run1p9j;
    QString outDir;
    QString flat720Coverage;
    QString flat1p9jCoverage;
    QString reportPath;
};

bool RunSpineReachabilityProbe(const SpineProbeConfig& cfg, QString* err_out = nullptr);

}  // namespace h5reader::rediscover
