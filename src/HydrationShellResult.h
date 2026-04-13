#pragma once
//
// HydrationShellResult: per-atom hydration shell geometry from
// explicit water positions.
//
// For each protein atom, characterises the water packing:
//   - Half-shell asymmetry (UCSB method): what fraction of first-shell
//     waters are on the solvent-exposed side vs the protein interior.
//     Computed using the solvent-accessibility vector (from atom to
//     center of mass of first-shell waters vs center of protein).
//   - Mean water dipole orientation: how ordered are the water dipoles
//     around this atom (cos angle between water dipole and atom→water vector).
//   - Nearest ion distance and charge.
//
// These quantities encode the local dielectric environment that no
// geometry-only kernel can see — you cannot know from protein
// coordinates alone whether an atom faces a water cavity or a
// hydrophobic pocket.
//
// Parameters (from TOML):
//   water_first_shell_cutoff — first hydration shell boundary (default 3.5 A)
//                              Shared with WaterFieldResult.
//   hydration_ion_cutoff     — nearest-ion search distance (default 20.0 A)
//
// No KernelFilterSet — no geometric kernel is evaluated. This calculator
// counts waters and computes cos angles, not 1/r^n fields.
//
// GeometryChoice: one summary record (parameters used).
//
// Dependencies: SpatialIndexResult.
// Requires: SolventEnvironment (passed to Compute).
//

#include "ConformationResult.h"
#include "SpatialIndexResult.h"
#include "SolventEnvironment.h"
#include "ProteinConformation.h"

namespace nmr {

class HydrationShellResult : public ConformationResult {
public:
    std::string Name() const override { return "HydrationShellResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { std::type_index(typeid(SpatialIndexResult)) };
    }

    static std::unique_ptr<HydrationShellResult> Compute(
        ProteinConformation& conf,
        const SolventEnvironment& solvent);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
