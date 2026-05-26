#pragma once
//
// Projection of Protein's prepared force-field charge table onto the
// per-conformation atom fields.
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

    // Ensures Protein has force-field charges, then projects them to atoms.
    static std::unique_ptr<ChargeAssignmentResult> Compute(
        ProteinConformation& conf,
        const ChargeSource& source);

    static std::unique_ptr<ChargeAssignmentResult> Compute(
        ProteinConformation& conf,
        const std::string& param_file_path);

    double ChargeAt(size_t atom_index) const;
    double PbRadiusAt(size_t atom_index) const;
    const ForceFieldChargeTable& ChargeTable() const;

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
