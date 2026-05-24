#pragma once
//
// HaighMallionResult: surface integral ring current shielding tensors.
//
// For each aromatic ring within 15A of each atom, computes:
//
//   1. Raw surface integral H_ab = integral_S K_ab dS where K is the
//      dipolar kernel (3 rho_a rho_b / rho^5 - delta_ab / rho^3).
//      Fan triangulation from ring centroid, 7-point Gaussian quadrature
//      (Stroud T2:5-1 / Dunavant degree-5). Adaptive subdivision at
//      2.0 A (level 1) and 1.0 A (level 2).
//      H is SYMMETRIC and TRACELESS (integral of traceless integrand).
//      Units: Angstrom^-1. Stored as hm_H_tensor on RingNeighbourhood.
//
//   2. Full shielding kernel G_ab = n_b * (H . n)_a where n is the
//      ring normal. This is the physically correct shielding tensor
//      construction: the induced current depends on flux (n . B_0),
//      so the b-index carries n_b. The "effective B-field" is V = H.n
//      (the surface integral contracted with the normal).
//      G is RANK-1 (like Biot-Savart), with T0 != 0, T1 != 0, T2 != 0.
//      Stored as hm_G_tensor on RingNeighbourhood.
//
// This parallels BiotSavartResult: both model pi-electron circulation.
// BS uses wire segments (line integral); HM uses surface quadrature.
// Both produce rank-1 shielding tensors (G = n (x) V), but compute
// the "effective B-field" V differently. Whether they agree on T2 is
// the empirical question.
//
// Boyd & Skrynnikov JACS 2002 124:1832 eq 3 gives the tensor structure.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <array>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class HaighMallionResult : public ConformationResult {
public:
    std::string Name() const override { return "HaighMallionResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute HM surface integral tensors for all atoms.
    static std::unique_ptr<HaighMallionResult> Compute(
        ProteinConformation& conf);

    // Grid sampling: evaluate HM shielding kernel at an arbitrary 3D point.
    SphericalTensor SampleShieldingAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
