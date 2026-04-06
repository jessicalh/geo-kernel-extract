#pragma once
//
// McConnellResult: Bond magnetic anisotropy shielding tensors.
//
// For each bond within 10A of each atom, computes the full McConnell
// shielding tensor (asymmetric, non-traceless, T0+T1+T2) and the
// symmetric traceless dipolar kernel K. Accumulates per-category
// totals and tracks nearest CO/CN bonds.
//
// The full tensor formula (from GEOMETRIC_KERNEL_CATALOGUE.md):
//
//   M_ab = 9 cos_theta d_hat_a b_hat_b
//        - 3 b_hat_a b_hat_b
//        - (3 d_hat_a d_hat_b - delta_ab)
//
//   Shielding contribution = M_ab / r^3  (Angstrom^-3)
//
// The dipolar kernel K_ab = (3 d_hat_a d_hat_b - delta_ab) / r^3
// is the symmetric traceless part, stored separately for features.
//
// The McConnell scalar f = (3 cos^2 theta - 1) / r^3 is the double
// contraction of K with the bond direction. It is NOT the tensor trace
// (which comes from the full M, not from K alone).
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

// Cutoff radius for bond midpoint search (Angstroms)
constexpr double MCCONNELL_CUTOFF_A = 10.0;


class McConnellResult : public ConformationResult {
public:
    std::string Name() const override { return "McConnellResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute McConnell tensors for all atoms.
    static std::unique_ptr<McConnellResult> Compute(
        ProteinConformation& conf);

    // Query methods
    double CategorySum(size_t atom_index, BondCategory cat) const;
    double NearestCOContribution(size_t atom_index) const;

    // Grid sampling: evaluate McConnell kernel at an arbitrary 3D point.
    SphericalTensor SampleShieldingAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
