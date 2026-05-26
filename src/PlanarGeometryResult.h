#pragma once
// Per-frame planar geometry derived from positions and the LegacyAmber
// typed substrate. NaN marks not-applicable per-residue/per-ring rows;
// per-atom pyramidalization is zero for atoms with no planar group.
// Omega values are radians, pucker Q is Angstroms, and pucker theta is
// degrees.

#include "ConformationResult.h"
#include "Types.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;


class PlanarGeometryResult : public ConformationResult {
public:
    std::string Name() const override { return "PlanarGeometryResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<PlanarGeometryResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    const std::vector<double>& OmegaActual() const { return omega_actual_; }
    const std::vector<double>& OmegaDeviation() const { return omega_deviation_; }
    const std::vector<uint8_t>& OmegaIsXpro() const { return omega_is_xpro_; }

    const std::vector<double>& AromaticChi2() const { return aromatic_chi2_; }
    const std::vector<double>& PuckerQ() const { return pucker_Q_; }
    const std::vector<double>& PuckerTheta() const { return pucker_theta_; }

private:
    const ProteinConformation* conf_ = nullptr;

    // Indexed by Protein residue index. X-to-Pro peptide bonds keep
    // the actual omega value and are marked in omega_is_xpro_.
    std::vector<double>  omega_actual_;
    std::vector<double>  omega_deviation_;
    std::vector<uint8_t> omega_is_xpro_;   // 1 if i+1 is Pro, else 0; 0 at C-term

    // Indexed by LegacyAmberTopology::AromaticRingAt index.
    std::vector<double> aromatic_chi2_;

    // Indexed by LegacyAmberTopology::SaturatedAt index.
    std::vector<double> pucker_Q_;
    std::vector<double> pucker_theta_;
};

}  // namespace nmr
