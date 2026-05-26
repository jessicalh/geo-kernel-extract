#pragma once
//
// Linearized Poisson-Boltzmann E-field and EFG from APBS. The bridge returns
// a potential grid in kT/e; this result stores central-difference E in V/A
// and EFG in V/A^2 after KT_OVER_E_298K conversion.
//
// Compute() returns nullptr on APBS failure; there is no vacuum fallback.
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

    static std::unique_ptr<ApbsFieldResult> Compute(ProteinConformation& conf);

    Vec3 ElectricFieldAt(size_t atom_index) const;
    Mat3 FieldGradientAt(size_t atom_index) const;
    SphericalTensor FieldGradientSphericalAt(size_t atom_index) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
