#pragma once
//
// Reports protonation_variant_index values already resolved during Protein
// construction. This result does not infer variants from atom names.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"

namespace nmr {

class ProtonationDetectionResult : public ConformationResult {
public:
    std::string Name() const override { return "ProtonationDetectionResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<ProtonationDetectionResult> Compute(
        ProteinConformation& conf);

    int AssignedCount() const { return assigned_count_; }

    int UnresolvedCount() const { return unresolved_count_; }

    // Empty string if the residue is not titratable or has no resolved variant.
    std::string VariantNameAt(size_t residue_index) const;

private:
    const ProteinConformation* conf_ = nullptr;
    int assigned_count_ = 0;
    int unresolved_count_ = 0;

    std::vector<std::string> variant_names_;
};

}  // namespace nmr
