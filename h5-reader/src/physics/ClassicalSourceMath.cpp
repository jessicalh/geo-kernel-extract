#include "ClassicalSourceMath.h"

#include "LiteratureAccessors.h"

#include <cmath>

namespace h5reader::physics {

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

ClassicalSigmaResult ComputeClassicalSigma(const ClassicalSigmaInputs& inputs) {
    ClassicalSigmaResult out;
    out.sigma0 = {inputs.sigma0.value,
                  inputs.sigma0.present && finite(inputs.sigma0.value)};
    out.ring = {inputs.ring.value,
                inputs.ring.present && finite(inputs.ring.value)};
    out.mcconnell = {inputs.mcconnell.value,
                     inputs.mcconnell.present && finite(inputs.mcconnell.value)};
    out.larsen = {inputs.larsen.value,
                  inputs.larsen.present && finite(inputs.larsen.value)};

    const bool ePresent = inputs.e_parallel_mopac.present
                          && finite(inputs.e_parallel_mopac.value);
    const bool aPresent = inputs.buckingham_A.present
                          && finite(inputs.buckingham_A.value);
    const bool bPresent = inputs.buckingham_B.present
                          && finite(inputs.buckingham_B.value);
    if (ePresent && aPresent) {
        out.buckingham_linear = {
            -inputs.buckingham_A.value * inputs.e_parallel_mopac.value,
            true,
        };
    }
    if (ePresent && bPresent) {
        out.buckingham_quadratic = {
            -inputs.buckingham_B.value
                * inputs.e_parallel_mopac.value
                * inputs.e_parallel_mopac.value,
            true,
        };
    }
    if (out.buckingham_linear.present || out.buckingham_quadratic.present) {
        double v = 0.0;
        if (out.buckingham_linear.present) v += out.buckingham_linear.value;
        if (out.buckingham_quadratic.present) v += out.buckingham_quadratic.value;
        out.buckingham = {v, true};
    }

    bool anyPresent = false;
    const double sigmaCl = FoldClassicalSigma(
        out.sigma0.present ? out.sigma0.value : std::numeric_limits<double>::quiet_NaN(),
        {out.buckingham_linear,
         out.buckingham_quadratic,
         out.ring,
         out.mcconnell,
         out.larsen},
        &anyPresent);
    out.sigma_cl = {sigmaCl, anyPresent && finite(sigmaCl)};
    return out;
}

double McConnellProducerT0ToPpm(model::BondCategory category, double packedT0) {
    if (!finite(packedT0)) return std::numeric_limits<double>::quiet_NaN();
    return -nmr::constants::kMcConnellMolarPrefactor.value
           * McConnellDeltaChi(category).value
           * packedT0;
}

}  // namespace h5reader::physics
