#pragma once
//
// Per-residue protonation decisions. variant_index is an index into
// AminoAcidType::variants; -1 means the default charged form.
//

#include "Types.h"
#include <string>
#include <vector>
#include <cmath>
#include <limits>

namespace nmr {

struct ResidueProtonation {
    size_t residue_index = std::numeric_limits<size_t>::max();
    AminoAcid amino_acid = AminoAcid::Unknown;
    int variant_index = -1;       // Into AminoAcidType::variants, -1 = default (charged)
    double pKa = std::nan("");    // Predicted pKa (NaN if not computed)
    bool is_charged = false;      // True if residue carries a formal charge in this state
};


class ProtonationState {
public:
    ProtonationState() = default;

    ProtonationState(const std::string& name, double pH,
                     ProtonationTool tool, const std::string& tool_version)
        : name_(name), pH_(pH), tool_(tool), tool_version_(tool_version) {}

    bool IsEmpty() const { return residues_.empty() && tool_ == ProtonationTool::Manual; }

    const std::string& Name() const { return name_; }
    double pH() const { return pH_; }
    ProtonationTool Tool() const { return tool_; }
    const std::string& ToolVersion() const { return tool_version_; }

    void AddResidue(ResidueProtonation decision);

    size_t DecisionCount() const { return residues_.size(); }

    const std::vector<ResidueProtonation>& Decisions() const { return residues_; }

    // Lookup by residue index. Returns nullptr if no decision for this residue.
    const ResidueProtonation* ForResidue(size_t residue_index) const;

    // Only counts explicit decisions; termini are not included.
    int NetDecisionCharge() const;

    // Includes default charged states for residue_types without explicit
    // decisions; termini are not included.
    int NetChargeForProtein(const std::vector<AminoAcid>& residue_types) const;

    std::string Describe() const;

    bool operator==(const ProtonationState& other) const;
    bool operator!=(const ProtonationState& other) const { return !(*this == other); }

private:
    std::string name_;                  // e.g., "propka_pH7.0", "tleap_ff14SB_default"
    double pH_ = std::nan("");          // NaN for default protonation
    ProtonationTool tool_ = ProtonationTool::Manual;
    std::string tool_version_;          // "PROPKA 3.5.1", "KaML-CBTrees", "tleap ff14SB"
    std::vector<ResidueProtonation> residues_;
};

}  // namespace nmr
