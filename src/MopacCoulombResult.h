#pragma once
//
// MopacCoulombResult: Coulomb E-field and EFG from MOPAC QM charges.
//
// Same dipolar kernel as CoulombResult:
//   E_a(i) = ke * sum_{j!=i} q_j * (r_i - r_j)_a / |r_i - r_j|^3
//   V_ab(i) = ke * sum_{j!=i} q_j * [3 (r_i-r_j)_a (r_i-r_j)_b / |r_i-r_j|^5
//                                     - delta_ab / |r_i-r_j|^3]
//
// Different charge source: reads mopac_charge (PM7 Mulliken, conformation-
// dependent) instead of partial_charge (ff14SB, fixed per atom type).
// The T2 angular pattern differs because MOPAC charges respond to the
// local electronic environment. The model learns gamma_mopac.
//
// Decomposed by source: backbone, sidechain, aromatic (same as CoulombResult).
// No APBS solvent subtraction (APBS was solved with ff14SB charges).
//
// Units: V/A (E-field), V/A^2 (EFG), same as CoulombResult.
//
// Dependencies: MopacResult (charges on ConformationAtom),
//               SpatialIndexResult (ensures geometry is resolved).
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;


class MopacCoulombResult : public ConformationResult {
public:
    std::string Name() const override { return "MopacCoulombResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute E-field and EFG from MOPAC charges for all atoms.
    static std::unique_ptr<MopacCoulombResult> Compute(
        ProteinConformation& conf);

    // Query methods
    Vec3 EFieldAt(size_t atom_index) const;
    Mat3 EFGAt(size_t atom_index) const;
    SphericalTensor EFGSphericalAt(size_t atom_index) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
