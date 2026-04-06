#pragma once
//
// HBondResult: hydrogen bond dipolar shielding tensors.
//
// For each atom involved in an H-bond (identified by DSSP), computes
// the full shielding tensor from the dipolar interaction with the
// H-bond partner. Same derivation as McConnell/RingSusceptibility
// (GEOMETRIC_KERNEL_CATALOGUE.md) with b_hat replaced by the
// donor→acceptor direction h_hat.
//
// H-bond identification comes from DsspResult, which provides
// backbone H-bond partners via the Kabsch-Sander energy criterion.
// The H-bond direction h_hat is computed from the actual atomic
// positions (donor N → acceptor O for backbone H-bonds).
//
// The full tensor:
//   M_ab = 9 cosθ d̂_a h_b - 3 h_a h_b - (3 d̂_a d̂_b - δ_ab)
//   Shielding contribution = M_ab / r³  (Angstrom⁻³)
//
// where d = r_atom - r_partner, h_hat = direction from donor to acceptor.
//
// Asymmetric, non-traceless, T0+T1+T2 all non-zero.
//

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

    // Factory: compute H-bond dipolar tensors for all atoms.
    static std::unique_ptr<HBondResult> Compute(
        ProteinConformation& conf);

    // Grid sampling: evaluate H-bond kernel at arbitrary 3D point.
    SphericalTensor SampleShieldingAt(Vec3 point) const;

    // Resolved H-bond count (for diagnostics)
    size_t HBondCount() const { return hbond_midpoints_.size(); }

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;

    // Stored for SampleAt: resolved H-bond geometry
    std::vector<Vec3> hbond_midpoints_;
    std::vector<Vec3> hbond_directions_;   // donor N → acceptor O
    std::vector<double> hbond_distances_;  // N...O distance
};

}  // namespace nmr
