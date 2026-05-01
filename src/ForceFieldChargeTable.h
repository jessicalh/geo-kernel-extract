#pragma once
//
// ForceFieldChargeTable: loaded force-field partial charges and PB radii.
//
// This is the prepared per-atom charge/PB-radius object owned by Protein.
// ChargeAssignmentResult projects this table into the old conformation fields;
// it does not own the authoritative loaded charge model.
//

#include "ChargeSource.h"
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class Protein;
class ProteinConformation;

class ForceFieldChargeTable {
public:
    ForceFieldChargeTable(ForceField source_force_field,
                          ChargeModelKind kind,
                          std::string source_description,
                          std::vector<AtomChargeRadius> values);

    static std::unique_ptr<ForceFieldChargeTable> Build(
        const ChargeSource& source,
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out);

    static std::unique_ptr<ForceFieldChargeTable> FromValues(
        ForceField source_force_field,
        ChargeModelKind kind,
        std::string source_description,
        std::vector<AtomChargeRadius> values,
        size_t expected_atom_count,
        std::string& error_out);

    ForceField SourceForceField() const { return source_force_field_; }
    ChargeModelKind Kind() const { return kind_; }
    const std::string& SourceDescription() const { return source_description_; }

    size_t AtomCount() const { return values_.size(); }
    double PartialChargeAt(size_t atom_index) const {
        return values_[atom_index].partial_charge;
    }
    double PbRadiusAt(size_t atom_index) const {
        return values_[atom_index].pb_radius;
    }
    bool MatchedAt(size_t atom_index) const {
        return values_[atom_index].status == ChargeAssignmentStatus::Matched;
    }
    bool PbRadiusAuthoritativeAt(size_t atom_index) const {
        return values_[atom_index].status == ChargeAssignmentStatus::Matched;
    }
    const std::vector<AtomChargeRadius>& Values() const { return values_; }

    double TotalCharge() const;
    int AssignedCount() const;
    int UnassignedCount() const;
    int NonAuthoritativePbRadiusCount() const;

private:
    ForceField source_force_field_ = ForceField::Unknown;
    ChargeModelKind kind_ = ChargeModelKind::Unknown;
    std::string source_description_;
    std::vector<AtomChargeRadius> values_;
};

}  // namespace nmr
