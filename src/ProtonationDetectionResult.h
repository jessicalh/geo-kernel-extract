#pragma once
//
// ProtonationDetectionResult: detect protonation state of titratable
// residues by examining which hydrogen atoms are present in the structure.
//
// PDB LOADING BOUNDARY code. Reads atom name strings ONCE to determine
// which hydrogens are present on titratable sidechains, then stores
// the result as a typed protonation_variant_index on each Residue.
//
// Detection logic (from hydrogen names in the PDB coordinate set):
//   HIS: "HD1" present -> HID, "HE2" present -> HIE, both -> HIP, neither -> ambiguous
//   ASP: "HD2" present -> ASH (protonated), absent -> ASP (charged)
//   GLU: "HE2" present -> GLH, absent -> GLU (charged)
//   CYS: SG bonded to another SG -> CYX, "HG" present -> CYS
//   LYS: all three HZ -> LYS (charged), "HZ3" absent -> LYN
//
// Dependencies: none (reads topology/atom names only).
//
// After this result, DetectAromaticRings should read from
// protonation_variant_index for HIS tautomer, not do its own string check.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"

namespace nmr {

class ProtonationDetectionResult : public ConformationResult {
public:
    std::string Name() const override { return "ProtonationDetectionResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Factory: detect protonation from hydrogen atoms present in the structure.
    // Sets protonation_variant_index on each titratable Residue in the Protein.
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
