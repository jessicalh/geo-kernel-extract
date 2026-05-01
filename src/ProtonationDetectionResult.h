#pragma once
//
// ProtonationDetectionResult: report prepared protonation state of titratable
// residues.
//
// Protonation/variant identity is resolved during Protein construction.
// This result remains as compatibility/reporting plumbing for the existing
// ConformationResult surface; it does not perform atom-name lookup.
//
// Prepared-state projection:
//   HIS: "HD1" present -> HID, "HE2" present -> HIE, both -> HIP, neither -> ambiguous
//   ASP: "HD2" present -> ASH (protonated), absent -> ASP (charged)
//   GLU: "HE2" present -> GLH, absent -> GLU (charged)
//   CYS: SG bonded to another SG -> CYX, "HG" present -> CYS
//   LYS: all three HZ -> LYS (charged), "HZ3" absent -> LYN
//   ARG: ARG charged default; ARN is not inferred from names in this model
//
// Dependencies: none. The construction/loading boundary performs any atom-name
// interpretation needed to prepare the residue fields this result reports.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"

namespace nmr {

class ProtonationDetectionResult : public ConformationResult {
public:
    std::string Name() const override { return "ProtonationDetectionResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Factory: project prepared protonation state into the old result surface.
    // Returns nullptr only on catastrophic failure (which should not happen).
    static std::unique_ptr<ProtonationDetectionResult> Compute(
        ProteinConformation& conf);

    // Query: how many residues had a variant assigned?
    int AssignedCount() const { return assigned_count_; }

    // Query: how many residues were titratable but unresolved (no H)?
    int UnresolvedCount() const { return unresolved_count_; }

    // Per-residue variant name (empty string if not titratable or unresolved)
    std::string VariantNameAt(size_t residue_index) const;

private:
    const ProteinConformation* conf_ = nullptr;
    int assigned_count_ = 0;
    int unresolved_count_ = 0;

    // Cache of variant name per residue for query
    std::vector<std::string> variant_names_;
};

}  // namespace nmr
