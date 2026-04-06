#pragma once
//
// ApbsFieldResult: electrostatic E-field and EFG at each atom.
//
// Dependencies: ChargeAssignmentResult (charges and radii).
//
// Stores per atom: apbs_efield (Vec3), apbs_efg (Mat3),
//                  apbs_efg_spherical (SphericalTensor).
//
// Linearised Poisson-Boltzmann solve via the APBS C bridge (apbs_bridge.h).
// The bridge takes in-memory arrays of positions, charges, and radii and
// returns a potential grid in kT/e. E-field and EFG are extracted by
// central-difference interpolation on the grid.
//
// No fallback. If APBS fails, Compute() returns nullptr. Solvated fields
// require a working PB solver — substituting vacuum Coulomb would silently
// produce a different physical quantity in different units.
//
// E-field from grid: E = -grad(phi) via central differences
// EFG from grid: dE_a/dx_b via central differences of E-field
// Traceless projection applied (removes self-potential artifact).
//
// Units stored on ConformationAtom: V/A for E-field, V/A^2 for EFG.
// Converted from APBS native kT/(e*A) by multiplying by kT/e = 0.025693 V
// at 298.15 K. This makes APBS fields directly comparable to CoulombResult.
//

#include "ConformationResult.h"
#include "ChargeAssignmentResult.h"
#include "ProteinConformation.h"

namespace nmr {

class ApbsFieldResult : public ConformationResult {
public:
    std::string Name() const override { return "ApbsFieldResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { std::type_index(typeid(ChargeAssignmentResult)) };
    }

    // Factory: compute E-field and EFG at each atom from charges.
    static std::unique_ptr<ApbsFieldResult> Compute(ProteinConformation& conf);

    // Query methods
    Vec3 ElectricFieldAt(size_t atom_index) const;
    Mat3 FieldGradientAt(size_t atom_index) const;
    SphericalTensor FieldGradientSphericalAt(size_t atom_index) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
