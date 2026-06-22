#pragma once

#include <initializer_list>
#include <limits>
#include <vector>

namespace h5reader::rediscover {

struct ClassicalTermValue {
    double value = std::numeric_limits<double>::quiet_NaN();
    bool present = false;
};

double SampleSdFinite(const std::vector<double>& values);
double SdRatioScaleFactor(const std::vector<double>& sigmaQm,
                          const std::vector<double>& sigmaClassical);
double FoldClassicalSigma(double sigma0,
                          std::initializer_list<ClassicalTermValue> terms,
                          bool* anyPresent = nullptr);

}  // namespace h5reader::rediscover
