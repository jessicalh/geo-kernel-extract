#pragma once
//
// CoulombResult: vacuum Coulomb E-field and EFG at each atom.
//
// For each atom, computes the electric field E and electric field gradient
// tensor V from partial charges within the configured cutoff. Decomposes by source:
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

// Candidate source set: all atoms within coulomb_efield_cutoff (20 A default,
// TOML), followed by self and charge-noise filters.
// The 1/r^2 field is long-range but truncated at the configured radius.
//
// Two SEPARATE geometric kernels (not a single unified tensor):
//   E_a  (rank-1)             -> T0 shielding via Buckingham A,B parameters
//   V_ab (rank-2, symmetric,  -> T2 shielding via gamma
//         traceless)
// Unlike McConnell (where chi.K contraction produces an asymmetric tensor
// with non-zero T0+T1+T2 from geometry alone), there is no single "full
// tensor" that unifies E and V. coulomb_shielding_contribution stores the
// T2 (Decompose(EFG)) only; the T0 from E via Buckingham is not a pure
// geometric kernel and is applied at calibration.

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
