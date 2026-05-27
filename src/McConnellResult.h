#pragma once
//
// McConnellResult: bond magnetic-anisotropy tensor kernels.
//
// For each bond within the configured cutoff (default 10 A), computes the
// full McConnell tensor kernel, which can contain T0, T1, and T2 components,
// and the symmetric traceless dipolar kernel K. Accumulates per-category
// totals and tracks nearest CO/CN bonds.
//
// The full tensor formula (from GEOMETRIC_KERNEL_CATALOGUE.md):
//
//   M_ab = 9 cos_theta d_hat_a b_hat_b
//        - 3 b_hat_a b_hat_b
//        - (3 d_hat_a d_hat_b - delta_ab)
//
//   Stored kernel = M_ab / r^3  (Angstrom^-3)
//
// The dipolar kernel K_ab = (3 d_hat_a d_hat_b - delta_ab) / r^3
// is symmetric and traceless, and is stored separately for features.
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

// NOTE: vestigial. The runtime cutoff comes from CalculatorConfig
// "mcconnell_bond_anisotropy_cutoff"; this constant is referenced only
// in comments/docs, never in the math path.
constexpr double MCCONNELL_CUTOFF_A = 10.0;


class McConnellResult : public ConformationResult {
public:
    std::string Name() const override { return "McConnellResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute McConnell tensors for all atoms.
    static std::unique_ptr<McConnellResult> Compute(
        ProteinConformation& conf);

    // Query methods (return the McConnell scalar f sum, not a tensor)
    double CategoryScalarSum(size_t atom_index, BondCategory cat) const;
    double NearestCOScalarContribution(size_t atom_index) const;

    // Grid sampling: evaluate McConnell kernel at an arbitrary 3D point.
    SphericalTensor SampleKernelAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
