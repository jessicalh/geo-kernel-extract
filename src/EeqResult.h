#pragma once
//
// EeqResult: geometry-dependent partial charges from extended
// electronegativity equilibration (Caldeweyher et al. 2019).
//
// Minimises E(q) = Σ χᵢqᵢ + ½Σ ηᵢqᵢ² + ½Σ_{i≠j} qᵢqⱼγ(Rᵢⱼ)
// subject to Σqᵢ = Q_total (charge neutrality).
//
// Uses D4 parameters: element-specific electronegativity (χ), chemical
// hardness (η), CN-dependent EN shift (κ), and covalent radius (r_cov).
// Coordination number computed via error function counting (Caldeweyher 2019).
// Coulomb interaction via Ohno-Klopman kernel (Ohno 1964, Klopman 1964):
//   γ(R) = 1/√(R² + 1/(ηᵢ·ηⱼ))
//
// Pure C++ with Eigen.  No external binary, no CUDA dependency.
// One N×N linear solve per frame via Cholesky with block elimination
// for the charge neutrality constraint.
//
// Reference: Caldeweyher, Ehlert, Hansen, Neugebauer, Spicher,
// Bannwarth & Grimme, J. Chem. Phys. 150, 154122 (2019).
// DOI: 10.1063/1.5090222.
//
// Output:
//   eeq_charges.npy  (N,) float64 — partial charges (elementary charges)
//
// Parameters (from TOML):
//   eeq_total_charge — net system charge (default 0, neutral protein)
//
// No KernelFilterSet — charge calculation, not field evaluation.
// GeometryChoice: one summary record (parameters, charge statistics).
//
// Dependencies: none (reads protein geometry and element types only).
//

#include "ConformationResult.h"
#include "ProteinConformation.h"

namespace nmr {

class EeqResult : public ConformationResult {
public:
    std::string Name() const override { return "EeqResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    static std::unique_ptr<EeqResult> Compute(ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
