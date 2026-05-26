#pragma once
//
// Surface-integral ring-current kernel. H_ab integrates the dipolar kernel
// over a fan triangulation of the ring; G_ab = -n_b * (H . n)_a is the
// rank-1 shielding kernel stored on RingNeighbourhood.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <array>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class HaighMallionResult : public ConformationResult {
public:
    std::string Name() const override { return "HaighMallionResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<HaighMallionResult> Compute(
        ProteinConformation& conf);

    SphericalTensor SampleKernelAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
