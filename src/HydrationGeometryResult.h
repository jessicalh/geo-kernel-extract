#pragma once
//
// HydrationGeometryResult: per-atom water polarisation features using
// the SASA-derived surface normal as reference frame.
//
// For each protein atom, characterises the first-shell water orientation
// relative to the proper surface normal (from SasaResult), replacing the
// crude protein-COM direction used by HydrationShellResult.
//
// Output: water_polarization.npy (N, 10):
//   0-2  Net water dipole vector (Vec3) — direction + coherence of local
//        electrostatic asymmetry.  Sum of individual water dipoles within
//        first shell.
//   3-5  Surface normal (Vec3, from SASA) — reference frame for interpreting
//        water orientation.  Duplicated from SasaResult intentionally.
//   6    Half-shell asymmetry (SASA-normal) — fraction of first-shell waters
//        on the solvent-exposed side.
//   7    Dipole alignment (SASA-normal) — cos(net dipole, surface normal).
//   8    Dipole coherence (|sum d_i| / n) — ordered vs random.
//   9    First-shell count — statistical weight.
//
// No KernelFilterSet — no geometric kernel is evaluated.  This calculator
// computes projections and cos angles, not 1/r^n fields.
//
// GeometryChoice: one summary record (parameters used).
//
// Dependencies: SasaResult (for surface normal).
// Requires: SolventEnvironment (passed to Compute).
//
// Design rule: existing HydrationShellResult stays untouched.  The COM-based
// values continue to exist under the same names.  This calculator writes
// distinct NPY files with the SASA-normal reference frame.
//

#include "ConformationResult.h"
#include "SasaResult.h"
#include "SolventEnvironment.h"
#include "ProteinConformation.h"

namespace nmr {

class HydrationGeometryResult : public ConformationResult {
public:
    std::string Name() const override { return "HydrationGeometryResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { std::type_index(typeid(SasaResult)) };
    }

    static std::unique_ptr<HydrationGeometryResult> Compute(
        ProteinConformation& conf,
        const SolventEnvironment& solvent);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
