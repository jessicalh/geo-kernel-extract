#include "ClassicalSourceMath.h"

#include <cmath>

namespace h5reader::rediscover {

namespace {

bool finite(double v) { return std::isfinite(v); }

}  // namespace

double SampleSdFinite(const std::vector<double>& values) {
    double sum = 0.0;
    std::size_t n = 0;
    for (double v : values) {
        if (!finite(v)) continue;
        sum += v;
        ++n;
    }
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    const double mean = sum / static_cast<double>(n);
    double ss = 0.0;
    for (double v : values) {
        if (!finite(v)) continue;
        const double d = v - mean;
        ss += d * d;
    }
    return std::sqrt(ss / static_cast<double>(n - 1));
}

double SdRatioScaleFactor(const std::vector<double>& sigmaQm,
                          const std::vector<double>& sigmaClassical) {
    const double sdQm = SampleSdFinite(sigmaQm);
    const double sdCl = SampleSdFinite(sigmaClassical);
    if (!finite(sdQm) || !finite(sdCl) || !(sdCl > 0.0))
        return std::numeric_limits<double>::quiet_NaN();
    return sdQm / sdCl;
}

double FoldClassicalSigma(double sigma0,
                          std::initializer_list<ClassicalTermValue> terms,
                          bool* anyPresent) {
    double folded = 0.0;
    bool any = false;
    if (finite(sigma0)) {
        folded += sigma0;
        any = true;
    }
    for (const ClassicalTermValue& term : terms) {
        if (!term.present || !finite(term.value)) continue;
        folded += term.value;
        any = true;
    }
    if (anyPresent) *anyPresent = any;
    return any ? folded : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace h5reader::rediscover
