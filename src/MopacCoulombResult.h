#pragma once
//
// Coulomb E-field and EFG from conformation-dependent MOPAC charges.
// Uses the same kernel and units as CoulombResult, but reads
// ConformationAtom::mopac_charge and does not subtract APBS solvent fields.
//
//   E_a(i) = ke * sum_{j!=i} q_j * (r_i - r_j)_a / |r_i - r_j|^3
//   V_ab(i) = ke * sum_{j!=i} q_j * [3 (r_i-r_j)_a (r_i-r_j)_b / |r_i-r_j|^5
//                                     - delta_ab / |r_i-r_j|^3]
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

    static std::unique_ptr<MopacCoulombResult> Compute(
        ProteinConformation& conf);

    Vec3 EFieldAt(size_t atom_index) const;
    Mat3 EFGAt(size_t atom_index) const;
    SphericalTensor EFGSphericalAt(size_t atom_index) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
