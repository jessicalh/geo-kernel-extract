#pragma once
//
// PiQuadrupoleResult: quadrupole EFG kernel from aromatic ring pi-electrons.
//
// For each aromatic ring within range of each atom, computes the EFG
// tensor from a point axial quadrupole at the ring center along the
// ring normal. The pi-electron cloud creates a quadrupole electric
// field; its EFG at nearby nuclei couples to shielding via the
// Buckingham gamma coefficient.
//
// Derivation: the potential from an axial quadrupole Theta along n is
//   Phi = Theta * (3 cos^2 theta - 1) / (2 r^3)
//
// The EFG V_ab = -d^2 Phi / dx_a dx_b is obtained via Stone's T-tensor:
//   V_ab = -(1/3) Theta_cd T_abcd
// where T_abcd = d^4(1/r) / dx_a dx_b dx_c dx_d is the rank-4 interaction
// tensor (Stone, Theory of Intermolecular Forces, OUP 2013, Ch. 3).
//
// For Theta_cd = Theta(3 n_c n_d - delta_cd)/2, the delta contraction
// vanishes by Laplace (T_abcc = nabla^2 nabla^2 (1/r) = 0), giving
// the geometric kernel:
//
//   G_ab = 105 dn^2 d_a d_b / r^9
//        - 30 dn (n_a d_b + n_b d_a) / r^7
//        - 15 d_a d_b / r^7
//        + 6 n_a n_b / r^5
//        + delta_ab (3/r^5 - 15 dn^2/r^7)
//
// where d = r_atom - r_center, r = |d|, dn = d . n.
//
// Properties:
//   - Symmetric: G_ab = G_ba
//   - Traceless: Tr(G) = 0  (Laplace equation outside sources)
//   - Pure T2: T0 = 0, T1 = 0
//   - Leading radial decay: 1/r^5
//
// Also stores the scalar (3 cos^2 theta - 1) / r^4 for the Buckingham
// A-term (E-field from quadrupole -> T0 isotropic shift). This is a
// SEPARATE physical effect from the EFG -> T2 coupling.
//
// Units: G in Angstrom^-5 (d in A, r in A). Scalar in Angstrom^-4.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class PiQuadrupoleResult : public ConformationResult {
public:
    std::string Name() const override { return "PiQuadrupoleResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute pi-quadrupole EFG tensors for all atoms.
    static std::unique_ptr<PiQuadrupoleResult> Compute(
        ProteinConformation& conf);

    // Grid sampling: evaluate PQ EFG kernel at arbitrary 3D point.
    SphericalTensor SampleShieldingAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
