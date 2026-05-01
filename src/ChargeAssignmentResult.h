#pragma once
//
// ChargeAssignmentResult: per-atom partial charges and PB radii.
//
// Dependencies: none.
//
// Charges come from a typed ChargeSource: ff14SB param file, CHARMM36m
// topology (.tpr), or AMBER prmtop. The ChargeSource determines the
// force field; this result projects the prepared Protein charge table
// onto ConformationAtom.
//
// For protonation variants (HID/HIE/HIP, ASH, GLH, CYX/CYM, LYN,
// ARN, TYM): the ff14SB param file path uses the residue's
// protonation_variant_index
// to select the correct charge set. Other sources (tpr, prmtop) carry
// the protonation in their topology already.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "ForceFieldChargeTable.h"
#include <string>

namespace nmr {

class ChargeSource;

class ChargeAssignmentResult : public ConformationResult {
public:
    std::string Name() const override { return "ChargeAssignmentResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Typed factory: assign charges from any ChargeSource.
    // This is the preferred path. The ChargeSource determines the
    // force field and data format.
    static std::unique_ptr<ChargeAssignmentResult> Compute(
        ProteinConformation& conf,
        const ChargeSource& source);

    // ff14SB param file factory (convenience, delegates to typed path).
    static std::unique_ptr<ChargeAssignmentResult> Compute(
        ProteinConformation& conf,
        const std::string& param_file_path);

    // No no-arg Compute: every protein enters the system with real
    // charges; there is no stub fallback.

    // Query methods
    double ChargeAt(size_t atom_index) const;
    double PbRadiusAt(size_t atom_index) const;
    const ForceFieldChargeTable& ChargeTable() const;

    // Diagnostics
    double TotalCharge() const { return total_charge_; }
    int AssignedCount() const { return assigned_count_; }
    int UnassignedCount() const { return unassigned_count_; }
    const std::string& Source() const { return source_; }

private:
    const ProteinConformation* conf_ = nullptr;
    const ForceFieldChargeTable* charge_table_ = nullptr;
    double total_charge_ = 0.0;
    int assigned_count_ = 0;
    int unassigned_count_ = 0;
    std::string source_;
};

}  // namespace nmr
