#pragma once
//
// HBondResult: hydrogen-bond clean scalar and nearest-direction features.
//
// For each atom, accumulates the angular scalar from DSSP-resolved
// backbone H-bonds that pass the filter set and records nearest
// H-bond geometry.
//
// H-bond identification comes from DsspResult, which provides
// backbone H-bond partners via the Kabsch-Sander energy criterion.
// The H-bond direction h_hat is computed from the actual atomic
// positions (donor H → acceptor O for backbone H-bonds).
//
// The former McConnell-form rank-2 H-bond kernel tensor output is
// intentionally absent.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

// Production kernel surface used by HBondResult::Compute and by the
// independent geometry forcing tests.  This stays in the result's named
// namespace (rather than a shared utility): HBondResult owns the convention
// that the source point is the explicit donor H and the axis is H→O.
namespace hbond_result_detail {

struct KernelResult {
    double f = 0.0;
    double distance = 0.0;
};

KernelResult ComputeKernel(const Vec3& atom_pos,
                           const Vec3& donor_h_pos,
                           const Vec3& h_to_o_hat);

}  // namespace hbond_result_detail

class HBondResult : public ConformationResult {
public:
    std::string Name() const override { return "HBondResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute H-bond scalar/direction values for all atoms.
    static std::unique_ptr<HBondResult> Compute(
        ProteinConformation& conf);

    // Resolved H-bond count (for diagnostics)
    size_t HBondCount() const { return hbond_count_; }

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;

    size_t hbond_count_ = 0;
};

}  // namespace nmr
