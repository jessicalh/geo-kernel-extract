#pragma once
//
// EeqResult: geometry-dependent partial charges from a project-local
// QEq/EEQ-style charge-equilibration model.
//
// project-local QEq/EEQ-style charge-equilibration model with error-function coordination number, CN-dependent electronegativity shift, Gaussian self term, and Ohno-Klopman off-diagonal kernel; parameters are in-repo/project-local and are not a validated dftd4/multicharge port.
//
// Solves the quadratic EEQ system
//   E(q) = Σ χ_effᵢqᵢ + ½Σ Aᵢᵢqᵢ² + Σ_{i<j} qᵢqⱼγ(Rᵢⱼ)
// subject to Σqᵢ = Q_total (net-charge constraint).
//
// Uses project-local parameters: element-specific electronegativity (χ),
// chemical hardness (η), CN-dependent EN shift (κ), Gaussian charge radius,
// and covalent radius (r_cov).
// Coordination number is computed via error-function counting.
// Coulomb interaction via Ohno-Klopman kernel (Ohno 1964, Klopman 1964):
//   γ(R) = 1/√(R² + 1/(ηᵢ·ηⱼ))
//
// Pure C++ with Eigen.  No external binary, no CUDA dependency.
// One N×N linear solve per frame via Cholesky with block elimination
// for the net-charge constraint.
//
// Output:
//   eeq_charges.npy  (N,) float64 — partial charges (elementary charges)
//   eeq_cn.npy       (N,) float64 — coordination number (intermediate,
//                                   for traceability)
//   eeq_chi_eff.npy  (N,) float64 — CN-conditioned electronegativity
//   eeq_hardness.npy (N,2) float64 — [eta, A_ii] in atomic units
//
// Net charge: Compute(conf, net_charge) receives the real system net charge
// from OperationRunner (or 0 under the charge_conditioning_neutral knob) as the
// Q_total constraint. No longer a TOML parameter.
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

    static std::unique_ptr<EeqResult> Compute(ProteinConformation& conf, int net_charge);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
