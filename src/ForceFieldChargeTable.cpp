#include "ForceFieldChargeTable.h"
#include "ProteinConformation.h"

#include <string>
#include <utility>

namespace nmr {

ForceFieldChargeTable::ForceFieldChargeTable(
        ForceField source_force_field,
        ChargeModelKind kind,
        std::string source_description,
        std::vector<AtomChargeRadius> values)
    : source_force_field_(source_force_field)
    , kind_(kind)
    , source_description_(std::move(source_description))
    , values_(std::move(values)) {}

std::unique_ptr<ForceFieldChargeTable> ForceFieldChargeTable::Build(
        const ChargeSource& source,
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) {
    auto values = source.LoadCharges(protein, conf, error_out);
    if (values.empty()) return nullptr;

    return FromValues(source.SourceForceField(), source.Kind(),
                      source.Describe(), std::move(values),
                      conf.AtomCount(), error_out);
}

std::unique_ptr<ForceFieldChargeTable> ForceFieldChargeTable::FromValues(
        ForceField source_force_field,
        ChargeModelKind kind,
        std::string source_description,
        std::vector<AtomChargeRadius> values,
        size_t expected_atom_count,
        std::string& error_out) {
    if (values.size() != expected_atom_count) {
        error_out = "ForceFieldChargeTable: expected " +
                    std::to_string(expected_atom_count) + " atoms, have " +
                    std::to_string(values.size()) + " charge/PB-radius rows";
        return nullptr;
    }

    return std::make_unique<ForceFieldChargeTable>(
        source_force_field, kind, std::move(source_description),
        std::move(values));
}

double ForceFieldChargeTable::TotalCharge() const {
    double total = 0.0;
    for (const auto& v : values_) total += v.partial_charge;
    return total;
}

int ForceFieldChargeTable::AssignedCount() const {
    int assigned = 0;
    for (const auto& v : values_) {
        if (v.status != ChargeAssignmentStatus::Missing) ++assigned;
    }
    return assigned;
}

int ForceFieldChargeTable::UnassignedCount() const {
    int unassigned = 0;
    for (const auto& v : values_) {
        if (v.status == ChargeAssignmentStatus::Missing) ++unassigned;
    }
    return unassigned;
}

int ForceFieldChargeTable::NonAuthoritativePbRadiusCount() const {
    int count = 0;
    for (const auto& v : values_) {
        if (v.status != ChargeAssignmentStatus::Matched) ++count;
    }
    return count;
}

}  // namespace nmr
