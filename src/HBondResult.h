#pragma once

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class HBondResult : public ConformationResult {
public:
    std::string Name() const override { return "HBondResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<HBondResult> Compute(
        ProteinConformation& conf);

    SphericalTensor SampleKernelAt(Vec3 point) const;

    size_t HBondCount() const { return hbond_midpoints_.size(); }

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;

    std::vector<Vec3> hbond_midpoints_;
    std::vector<Vec3> hbond_directions_;   // donor N to acceptor O
    std::vector<double> hbond_distances_;  // N...O distance
};

}  // namespace nmr
