#include "ForceFieldChargeTable.h"
#include "Bond.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "ProteinConformation.h"

#include <algorithm>
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

int ForceFieldChargeTable::MaterializeDerivedMbondi2PbRadii(
        const Protein& protein) {
    int derived = 0;
    const size_t n = std::min(values_.size(), protein.AtomCount());
    for (size_t ai = 0; ai < n; ++ai) {
        AtomChargeRadius& v = values_[ai];
        if (v.status != ChargeAssignmentStatus::PlaceholderPbRadius) continue;

        const Atom& atom = protein.AtomAt(ai);
        bool h_bonded_to_nitrogen = false;
        if (atom.element == Element::H) {
            const size_t parent = atom.parent_atom_index;
            if (parent != SIZE_MAX && parent < protein.AtomCount()) {
                h_bonded_to_nitrogen =
                    protein.AtomAt(parent).element == Element::N;
            } else {
                for (size_t bi : atom.bond_indices) {
                    const Bond& b = protein.BondAt(bi);
                    const size_t other =
                        (b.atom_index_a == ai) ? b.atom_index_b : b.atom_index_a;
                    if (other < protein.AtomCount() &&
                        protein.AtomAt(other).element == Element::N) {
                        h_bonded_to_nitrogen = true;
                        break;
                    }
                }
            }
        }

        v.pb_radius = Mbondi2Radius(atom.element, h_bonded_to_nitrogen);
        v.status = ChargeAssignmentStatus::DerivedMbondi2PbRadius;
        ++derived;
    }
    return derived;
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

int ForceFieldChargeTable::MatchedPbRadiusCount() const {
    int count = 0;
    for (const auto& v : values_) {
        if (v.status == ChargeAssignmentStatus::Matched) ++count;
    }
    return count;
}

int ForceFieldChargeTable::DerivedMbondi2PbRadiusCount() const {
    int count = 0;
    for (const auto& v : values_) {
        if (v.status == ChargeAssignmentStatus::DerivedMbondi2PbRadius) ++count;
    }
    return count;
}

int ForceFieldChargeTable::PlaceholderPbRadiusCount() const {
    int count = 0;
    for (const auto& v : values_) {
        if (v.status == ChargeAssignmentStatus::PlaceholderPbRadius) ++count;
    }
    return count;
}

int ForceFieldChargeTable::MissingCount() const {
    int count = 0;
    for (const auto& v : values_) {
        if (v.status == ChargeAssignmentStatus::Missing) ++count;
    }
    return count;
}

int ForceFieldChargeTable::NonAuthoritativePbRadiusCount() const {
    return PlaceholderPbRadiusCount() + MissingCount();
}

}  // namespace nmr
