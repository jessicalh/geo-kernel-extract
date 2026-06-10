#pragma once
//
// PiQuadrupoleResult: pi-quadrupole Buckingham A-term scalar.
//
// For each accepted aromatic ring/atom pair, stores
// (3 cos^2 theta - 1) / r^4 for the Buckingham A-term. The former
// manufactured geometric EFG tensor output is intentionally absent.
// Units: scalar in Angstrom^-4.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class PiQuadrupoleResult : public ConformationResult {
public:
    std::string Name() const override { return "PiQuadrupoleResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute pi-quadrupole scalar for all atoms.
    static std::unique_ptr<PiQuadrupoleResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
