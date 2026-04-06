#pragma once
//
// RingSusceptibilityResult: ring magnetic susceptibility anisotropy tensors.
//
// For each aromatic ring within range of each atom, computes the full
// shielding tensor from the ring's bulk diamagnetic susceptibility
// anisotropy, treated as a point magnetic dipole at the ring center.
//
// Same derivation as McConnell (GEOMETRIC_KERNEL_CATALOGUE.md) with
// b_hat replaced by ring normal n_hat:
//
//   M_ab = 9 cosθ d̂_a n_b - 3 n_a n_b - (3 d̂_a d̂_b - δ_ab)
//
//   Shielding contribution = M_ab / r³  (Angstrom⁻³)
//
// This tensor is ASYMMETRIC (T1 ≠ 0) and NON-TRACELESS (T0 = scalar f).
// Same structure as McConnell, different source geometry.
//
// The symmetric traceless dipolar kernel K_ab is also stored separately
// for features.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class RingSusceptibilityResult : public ConformationResult {
public:
    std::string Name() const override { return "RingSusceptibilityResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute ring susceptibility tensors for all atoms.
    static std::unique_ptr<RingSusceptibilityResult> Compute(
        ProteinConformation& conf);

    // Grid sampling: evaluate ring susceptibility kernel at arbitrary 3D point.
    SphericalTensor SampleShieldingAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
