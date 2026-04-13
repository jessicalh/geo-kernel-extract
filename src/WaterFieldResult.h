#pragma once
//
// WaterFieldResult: per-atom electric field and EFG from explicit
// solvent water molecules.
//
// Same Coulomb kernel as CoulombResult, but summing over water
// charges (TIP3P: O = -0.834e, H = +0.417e) instead of protein
// partial charges.  Uses SpatialIndexResult for the protein atom
// positions; builds its own index over water oxygen positions.
//
// This is what APBS approximates with a continuum dielectric.
// The explicit field includes water orientation fluctuations,
// cavities, bridging water, structural water — none of which
// the continuum model captures.
//
// Output per protein atom:
//   water_efield        Vec3       V/A     (total from all water within cutoff)
//   water_efg           Mat3       V/A^2   (electric field gradient)
//   water_efg_spherical SphericalTensor
//   water_efield_first  Vec3       V/A     (first shell only, < 3.5A from O)
//   water_efg_first     Mat3+ST    V/A^2   (first shell EFG)
//   water_n_first       int                (water O count in first shell)
//   water_n_second      int                (water O count in second shell, 3.5-5.5A)
//
// Parameters (from TOML):
//   water_efield_cutoff       — E-field summation cutoff on O distance (default 15.0 A)
//   water_first_shell_cutoff  — first hydration shell boundary (default 3.5 A)
//   water_second_shell_cutoff — second hydration shell boundary (default 5.5 A)
//   singularity_guard_distance — MinDistanceFilter threshold (default 0.1 A)
//   efield_magnitude_sanity_clamp — E-field magnitude clamp (default 100.0 V/A)
//
// KernelFilterSet: MinDistanceFilter (Coulomb 1/r³ singularity guard).
//
// GeometryChoice: one summary record (parameters), per-atom records for
// singularity guard and E-field clamp events.
//
// Dependencies: SpatialIndexResult (for protein positions).
// Requires: SolventEnvironment (passed to Compute).
//

#include "ConformationResult.h"
#include "SpatialIndexResult.h"
#include "SolventEnvironment.h"
#include "ProteinConformation.h"

namespace nmr {

class WaterFieldResult : public ConformationResult {
public:
    std::string Name() const override { return "WaterFieldResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { std::type_index(typeid(SpatialIndexResult)) };
    }

    // Factory: compute per-atom water E-field + EFG + shell counts.
    // solvent must contain water positions for this frame.
    static std::unique_ptr<WaterFieldResult> Compute(
        ProteinConformation& conf,
        const SolventEnvironment& solvent);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
