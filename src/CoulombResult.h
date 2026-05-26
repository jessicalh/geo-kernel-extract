#pragma once
//
// Vacuum Coulomb E-field and EFG from force-field partial charges.
// E is stored in V/A and EFG in V/A^2 after multiplying raw sums by
// COULOMB_KE.
//
// EFG kernel:
//   K_ab = (3 d_a d_b / r^5 - delta_ab / r^3)
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

// Two separate geometric kernels:
//   E_a  (rank-1)             -> T0 shielding via Buckingham A,B parameters
//   V_ab (rank-2, symmetric,  -> T2 shielding via gamma
//         traceless)
// coulomb_shielding_contribution stores Decompose(EFG) only; the T0 from E
// is applied at calibration.

class CoulombResult : public ConformationResult {
public:
    std::string Name() const override { return "CoulombResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<CoulombResult> Compute(
        ProteinConformation& conf);

    Vec3 EFieldAt(size_t atom_index) const;
    Mat3 EFGAt(size_t atom_index) const;
    SphericalTensor EFGSphericalAt(size_t atom_index) const;

    Vec3 SampleEFieldAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
