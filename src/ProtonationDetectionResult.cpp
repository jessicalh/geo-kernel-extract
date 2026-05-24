#include "ProtonationDetectionResult.h"
#include "Protein.h"
#include "OperationLog.h"

namespace nmr {

std::unique_ptr<ProtonationDetectionResult> ProtonationDetectionResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("ProtonationDetectionResult::Compute",
        "residues=" + std::to_string(conf.ProteinRef().ResidueCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t n_residues = protein.ResidueCount();

    auto result = std::make_unique<ProtonationDetectionResult>();
    result->conf_ = &conf;
    result->variant_names_.resize(n_residues);

    for (size_t ri = 0; ri < n_residues; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        const AminoAcidType& aatype = res.AminoAcidInfo();

        if (!aatype.is_titratable || aatype.variants.empty()) continue;

        if (!res.protonation_state_resolved) {
            result->unresolved_count_++;
            continue;
        }

        const int variant_idx = res.protonation_variant_index;
        if (variant_idx >= 0 &&
            variant_idx < static_cast<int>(aatype.variants.size())) {
            result->variant_names_[ri] = aatype.variants[variant_idx].name;
            result->assigned_count_++;

            OperationLog::Log(OperationLog::Level::Info, LogProtonation,
                "ProtonationDetectionResult",
                "residue " + std::to_string(res.sequence_number) +
                " " + ThreeLetterCodeForAminoAcid(res.type) +
                " -> " + result->variant_names_[ri]);
        }
    }

    OperationLog::Log(OperationLog::Level::Info, LogProtonation,
        "ProtonationDetectionResult::Compute",
        "assigned=" + std::to_string(result->assigned_count_) +
        " unresolved=" + std::to_string(result->unresolved_count_));

    return result;
}


std::string ProtonationDetectionResult::VariantNameAt(size_t residue_index) const {
    if (residue_index >= variant_names_.size()) return "";
    return variant_names_[residue_index];
}

}  // namespace nmr
