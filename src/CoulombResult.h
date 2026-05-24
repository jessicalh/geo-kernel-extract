#pragma once
//
// CoulombResult: vacuum Coulomb E-field and EFG at each atom.
//
// For each atom, computes the electric field E and electric field gradient
// tensor V from all other partial charges. Decomposes by source:
// backbone, sidechain, aromatic. Computes solvent contribution as
// APBS - vacuum total (if ApbsFieldResult is present).
//
// The EFG tensor V_ab = sum_j q_j * K_ab(r_ij) uses the same dipolar
// kernel as McConnell:
//   K_ab = (3 d_a d_b / r^5 - delta_ab / r^3)
//
// Each term is traceless by Gauss's law (no charge at the field point).
// The sum is traceless. Traceless projection applied after accumulation
// to correct floating-point drift.
//
// Units:
//   E-field: V/A  (raw sum in e/A^2, multiplied by COULOMB_KE = 14.3996)
//   EFG:     V/A^2 (raw sum in e/A^3, multiplied by COULOMB_KE)
//   Same units as ApbsFieldResult for direct comparison.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

// Coulomb cutoff: sum over ALL atoms (no spatial cutoff).
// N^2 is cheap for N < 1000. The physics is long-range.

class CoulombResult : public ConformationResult {
public:
    std::string Name() const override { return "CoulombResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute E-field and EFG for all atoms.
    static std::unique_ptr<CoulombResult> Compute(
        ProteinConformation& conf);

    // Query methods
    Vec3 EFieldAt(size_t atom_index) const;
    Mat3 EFGAt(size_t atom_index) const;
    SphericalTensor EFGSphericalAt(size_t atom_index) const;

    // Grid sampling: evaluate Coulomb E-field at arbitrary 3D point.
    Vec3 SampleEFieldAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
