#pragma once
//
// EeqCoulombResult: vacuum Coulomb E and EFG fields evaluated from the
// geometry-dependent project-local QEq/EEQ charges stored by EeqResult
// (project-local model, not a validated dftd4/multicharge port — see EeqResult.h).
//
// This is a normal ConformationResult calculator: Compute is the sole writer
// of its ConformationAtom fields and WriteFeatures is pure read-back.  It uses
// the same source classification, cutoff, sign convention, traceless
// projection, units, and E-vector clamp as CoulombResult. Unlike the legacy
// FF calculator, nonfinite input/output fails loudly instead of being
// rewritten to zero.
//
// Units: E in V/A; EFG in V/A^2.  The derivative hierarchy stops at rank 2.

#include "ConformationResult.h"
#include "Types.h"

#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;

class EeqCoulombResult : public ConformationResult {
public:
    std::string Name() const override { return "EeqCoulombResult"; }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<EeqCoulombResult> Compute(
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
