#pragma once
//
// Geometry-dependent EEQ partial charges using the D4 parameter table.
// The solve enforces sum(q_i) = eeq_total_charge.
//
// Objective:
//   E(q) = sum chi_i q_i + 1/2 sum eta_i q_i^2
//        + 1/2 sum_{i!=j} q_i q_j gamma(R_ij)
//   γ(R) = 1/√(R² + 1/(ηᵢ·ηⱼ))
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
