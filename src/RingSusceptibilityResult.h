#pragma once
//
// RingSusceptibilityResult: ring magnetic-susceptibility anisotropy kernels.
//
// For each aromatic ring within range of each atom, computes the tensor
// kernel for a point magnetic dipole at the ring center with its axis along
// the ring normal.
//
// Same derivation as McConnell (GEOMETRIC_KERNEL_CATALOGUE.md) with
// b_hat replaced by ring normal n_hat:
//
//   M_ab = 9 cosθ d̂_a n_b - 3 n_a n_b - (3 d̂_a d̂_b - δ_ab)
//
//   Stored kernel = M_ab / r³  (Angstrom⁻³)
//
// This tensor can be asymmetric (T1) and non-traceless (T0 = scalar f).
// Same structure as McConnell, different source geometry.
//
// The helper computes the symmetric traceless dipolar kernel K_ab, but the
// stored per-ring fields use the full M_ab / r³ tensor and scalar f.
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
    SphericalTensor SampleKernelAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
