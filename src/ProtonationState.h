#pragma once
//
// ProtonationState: per-residue protonation decisions for a protein.
//
// This is load-bearing, not metadata:
//   - Determines partial charges (ASP=-1, ASH=0)
//   - Determines ring type (HID vs HIE vs HIP → different RingTypeIndex)
//   - Determines which atoms exist (ASH has one more H than ASP)
//   - Is the input to tleap (ORCA prep) and pdb2gmx (GROMACS prep)
//
// A ProtonationState is a VALUE TYPE. It is created by a protonation tool
// (PROPKA, KaML, tleap default) and consumed by a topology builder (tleap,
// pdb2gmx) to produce a Protein with the correct hydrogen atoms and charges.
//
// The variant index into AminoAcidType::variants is the canonical form.
// The NamingRegistry translates to tool-specific names at the tool boundary.
//

#include "Types.h"
#include <string>
#include <vector>
#include <cmath>
#include <limits>

namespace nmr {

// Per-residue protonation decision.
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

    // An empty state means no protonation decisions have been made.
    // Crystal structures without hydrogens start here.
    bool IsEmpty() const { return residues_.empty() && tool_ == ProtonationTool::Manual; }

    // Identity
    const std::string& Name() const { return name_; }
    double pH() const { return pH_; }
    ProtonationTool Tool() const { return tool_; }
    const std::string& ToolVersion() const { return tool_version_; }

    // Add a per-residue decision.
    void AddResidue(ResidueProtonation decision);

    // Number of residues with protonation decisions.
    size_t DecisionCount() const { return residues_.size(); }

    // All decisions.
    const std::vector<ResidueProtonation>& Decisions() const { return residues_; }

    // Lookup by residue index. Returns nullptr if no decision for this residue.
    const ResidueProtonation* ForResidue(size_t residue_index) const;

    // Net formal charge from the decisions stored here.
    // Only counts residues with explicit decisions. For a complete protein
    // net charge, use NetChargeForProtein() which includes residues at
    // their default charged state (ARG +1, LYS +1, ASP -1, GLU -1, etc.).
    // Does NOT include backbone termini — those are from the topology builder.
    int NetDecisionCharge() const;

    // Net formal charge for a whole protein given this protonation state.
    // For each titratable residue: uses the decision if present, else the
    // default charged state from AminoAcidType::charged_formal_charge.
    // Non-titratable residues contribute 0.
    // residue_types: the AminoAcid enum for each residue in the protein.
    int NetChargeForProtein(const std::vector<AminoAcid>& residue_types) const;

    // Human-readable summary.
    std::string Describe() const;

    // Value semantics.
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
