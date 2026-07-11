#pragma once
//
// RingSusceptibilityResult: ring magnetic-susceptibility scalar rescue.
//
// For each accepted aromatic ring/atom pair, preserves the clean scalar
// f = (3cos²θ - 1) / r³ and shared ring-neighbour geometry. The former
// manufactured rank-2 ring-χ tensor output is intentionally absent.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

namespace ring_susceptibility_detail {

struct KernelResult {
    double scalar = 0.0;
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();
    bool valid = false;
};

// Exact production point-center susceptibility kernel.
KernelResult ComputeKernel(const Vec3& atom_pos,
                           const Vec3& ring_center,
                           const Vec3& ring_normal);

}  // namespace ring_susceptibility_detail

class RingSusceptibilityResult : public ConformationResult {
public:
    std::string Name() const override { return "RingSusceptibilityResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute ring susceptibility scalar for all atoms.
    static std::unique_ptr<RingSusceptibilityResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
