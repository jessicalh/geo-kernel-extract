#pragma once
//
// Point axial quadrupole EFG kernel for aromatic rings. G is symmetric,
// traceless, and has units Angstrom^-5; the scalar side channel has units
// Angstrom^-4.
//
//   G_ab = 105 dn^2 d_a d_b / r^9
//        - 30 dn (n_a d_b + n_b d_a) / r^7
//        - 15 d_a d_b / r^7
//        + 6 n_a n_b / r^5
//        + delta_ab (3/r^5 - 15 dn^2/r^7)
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class PiQuadrupoleResult : public ConformationResult {
public:
    std::string Name() const override { return "PiQuadrupoleResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<PiQuadrupoleResult> Compute(
        ProteinConformation& conf);

    // Returns the summed, decomposed geometric tensor, not a final shielding.
    SphericalTensor SampleKernelAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
