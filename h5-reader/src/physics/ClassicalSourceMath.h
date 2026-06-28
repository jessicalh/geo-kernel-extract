#pragma once

#include <initializer_list>
#include <limits>
#include <vector>

#include "../model/Types.h"

namespace h5reader::physics {

struct ClassicalTermValue {
    double value = std::numeric_limits<double>::quiet_NaN();
    bool present = false;
};

struct ClassicalSigmaInputs {
    ClassicalTermValue sigma0;
    ClassicalTermValue e_parallel_mopac;
    ClassicalTermValue buckingham_A;
    ClassicalTermValue buckingham_B;
    ClassicalTermValue ring;
    ClassicalTermValue mcconnell;
    ClassicalTermValue larsen;
};

struct ClassicalSigmaResult {
    ClassicalTermValue sigma0;
    ClassicalTermValue buckingham_linear;
    ClassicalTermValue buckingham_quadratic;
    ClassicalTermValue buckingham;
    ClassicalTermValue ring;
    ClassicalTermValue mcconnell;
    ClassicalTermValue larsen;
    ClassicalTermValue sigma_cl;
};

double SampleSdFinite(const std::vector<double>& values);
double SdRatioScaleFactor(const std::vector<double>& sigmaQm,
                          const std::vector<double>& sigmaClassical);
double FoldClassicalSigma(double sigma0,
                          std::initializer_list<ClassicalTermValue> terms,
                          bool* anyPresent = nullptr);
ClassicalSigmaResult ComputeClassicalSigma(const ClassicalSigmaInputs& inputs);

double McConnellProducerT0ToPpm(model::BondCategory category, double packedT0);

}  // namespace h5reader::physics
