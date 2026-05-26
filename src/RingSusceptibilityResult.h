#pragma once
//
// Ring magnetic-susceptibility kernel. Uses the McConnell tensor form with
// the bond direction replaced by the ring normal and stores M_ab/r^3.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class RingSusceptibilityResult : public ConformationResult {
public:
    std::string Name() const override { return "RingSusceptibilityResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<RingSusceptibilityResult> Compute(
        ProteinConformation& conf);

    SphericalTensor SampleKernelAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
