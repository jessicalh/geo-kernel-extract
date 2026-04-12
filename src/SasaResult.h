#pragma once
//
// SasaResult: per-atom solvent-accessible surface area (Shrake-Rupley).
//
// For each atom, distributes ~92 points on a sphere of radius
// (r_vdW + r_probe). Each point is checked for occlusion by nearby
// atoms via SpatialIndexResult. SASA = fraction_exposed * sphere_area.
//
// Dependencies: SpatialIndexResult.
//
// Parameters (from TOML):
//   sasa_probe_radius      — water probe radius (default 1.4 A)
//   sasa_n_points          — Fibonacci lattice points (default 92)
//
// Output: atom_sasa.npy — (N,) float64, Angstroms^2 per atom.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include <vector>

namespace nmr {

class SasaResult : public ConformationResult {
public:
    std::string Name() const override { return "SasaResult"; }
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<SasaResult> Compute(ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    double AtomSASA(size_t atom_index) const;
    const std::vector<double>& AllSASA() const { return sasa_; }

private:
    std::vector<double> sasa_;
};

}  // namespace nmr
