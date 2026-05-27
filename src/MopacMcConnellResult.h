#pragma once
//
// MopacMcConnellResult: Bond anisotropy shielding weighted by MOPAC bond order.
//
// Same McConnell kernel as McConnellResult:
//   M_ab = 9 cos_theta d_hat_a b_hat_b
//        - 3 b_hat_a b_hat_b
//        - (3 d_hat_a d_hat_b - delta_ab)
//
// Each bond's contribution is weighted by its MOPAC Wiberg bond order:
//   total += bond_order_k * M_ab_k / r_k^3
//
// The bond order is a measured QM quantity (electron density sharing),
// not a parameter. It modulates which bonds dominate the angular sum.
// A C=O with order 1.8 contributes more kernel than one with 1.2.
// The model learns Delta_chi per category for this weighted kernel.
//
// Same output grouping as McConnellResult: total shielding, category
// T2 for backbone/sidechain/aromatic/nearest-CO/nearest-CN, and scalar
// sums for CO/CN/sidechain/aromatic plus nearest distances.
//
// Dependencies: MopacResult (bond orders), SpatialIndexResult, GeometryResult.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

// NOTE: vestigial. The runtime cutoff comes from CalculatorConfig
// "mopac_mcconnell_bond_anisotropy_cutoff"; this constant is referenced
// only in comments/docs, never in the math path.
constexpr double MOPAC_MCCONNELL_CUTOFF_A = 10.0;


class MopacMcConnellResult : public ConformationResult {
public:
    std::string Name() const override { return "MopacMcConnellResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute bond-order-weighted McConnell tensors for all atoms.
    static std::unique_ptr<MopacMcConnellResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
