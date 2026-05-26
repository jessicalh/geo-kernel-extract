#pragma once
//
// Johnson-Bovey double-loop ring-current kernel. Positions are Angstroms
// at the public boundary; the wire-segment field is evaluated in SI and
// converted into the dimensionless shielding kernel
// G_ab = -n_b * B_a * PPM_FACTOR.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class BiotSavartResult : public ConformationResult {
public:
    std::string Name() const override { return "BiotSavartResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<BiotSavartResult> Compute(
        ProteinConformation& conf);

    // Uses ring geometries stored on the conformation by Compute().
    SphericalTensor SampleKernelAt(Vec3 point) const;

    Vec3 SampleBFieldAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
