#pragma once
//
// Bond magnetic-anisotropy kernel. ComputeBondKernel returns the full
// non-traceless tensor M_ab/r^3 plus the symmetric traceless dipolar kernel
// and scalar f for feature output.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

// Vestigial: the runtime cutoff comes from CalculatorConfig
// "mcconnell_bond_anisotropy_cutoff".
constexpr double MCCONNELL_CUTOFF_A = 10.0;


class McConnellResult : public ConformationResult {
public:
    std::string Name() const override { return "McConnellResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<McConnellResult> Compute(
        ProteinConformation& conf);

    // Return McConnell scalar f sums, not tensor components.
    double CategoryScalarSum(size_t atom_index, BondCategory cat) const;
    double NearestCOScalarContribution(size_t atom_index) const;

    SphericalTensor SampleKernelAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
