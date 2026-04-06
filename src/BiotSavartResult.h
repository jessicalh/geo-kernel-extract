#pragma once
//
// BiotSavartResult: Johnson-Bovey ring current B-field and shielding tensors.
//
// For each aromatic ring within 15A of each atom, computes:
//
//   1. B-field from Johnson-Bovey double-loop model (wire segments in SI,
//      two current loops at +/- lobe_offset from ring plane, each I/2).
//
//   2. Geometric kernel G_ab = n_b * B_a * PPM_FACTOR (rank-1 outer product
//      of ring normal with B-field). Dimensionless.
//
//   3. SphericalTensor decomposition of G:
//        T0 = (n . B) * PPM_FACTOR / 3  (isotropic ring current shift)
//        T1 != 0  (asymmetric: n_b B_a != n_a B_b)
//        T2 != 0  (angular anisotropy)
//
// Boyd & Skrynnikov JACS 2002 124:1832 showed that both T0 (eq 1, known
// since Johnson-Bovey 1957) and the off-diagonal components (eq 2, new)
// are needed for the full shielding tensor. Our G_ab = n_b B_a construction
// produces exactly their eq 3 decomposition.
//
// Sign convention: G_ab = -n_b * B_a * PPM_FACTOR. The minus sign comes
// from the shielding tensor definition sigma_ab = -dB_a^sec/dB_{0,b}.
// With this convention, sigma = I * G gives the correct physical sign
// using literature ring current intensities (I < 0 for diamagnetic).
// Verified: I=-12, atom 3A above PHE -> sigma = +1.40 ppm (shielded).
// In the ring plane at 5A -> sigma = -0.16 ppm (deshielded).
//
// Units: positions in Angstroms, current in nanoamperes, B computed in
// Tesla (SI internally), G dimensionless (after PPM_FACTOR).
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class BiotSavartResult : public ConformationResult {
public:
    std::string Name() const override { return "BiotSavartResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute Biot-Savart ring current fields for all atoms.
    static std::unique_ptr<BiotSavartResult> Compute(
        ProteinConformation& conf);

    // Grid sampling: evaluate BS kernel at an arbitrary 3D point.
    // Returns the shielding tensor (sum over all rings within cutoff).
    // Requires Compute() to have run (uses stored ring geometries).
    SphericalTensor SampleShieldingAt(Vec3 point) const;

    // Grid sampling: evaluate B-field at an arbitrary 3D point.
    Vec3 SampleBFieldAt(Vec3 point) const;

    // Write BS features: shielding contribution (9), per-type T0 (8),
    // per-type T2 (8x5), ring proximity counts, B-field totals.
    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
