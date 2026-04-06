#pragma once
//
// Bond: typed covalent bond between two atoms.
//
// TOPOLOGY (invariant): which atoms, bond order, category.
// GEOMETRY (conformation-dependent): length, direction, midpoint.
//

#include "Types.h"
#include <vector>

namespace nmr {

struct Bond {
    size_t atom_index_a = 0;
    size_t atom_index_b = 0;

    BondOrder    order = BondOrder::Unknown;
    BondCategory category = BondCategory::Unknown;
    bool         is_rotatable = false;

    // Geometry queries (conformation-dependent, take positions)
    Vec3 Midpoint(const std::vector<Vec3>& positions) const {
        return 0.5 * (positions[atom_index_a] + positions[atom_index_b]);
    }

    double Length(const std::vector<Vec3>& positions) const {
        return (positions[atom_index_b] - positions[atom_index_a]).norm();
    }

    Vec3 Direction(const std::vector<Vec3>& positions) const {
        Vec3 d = positions[atom_index_b] - positions[atom_index_a];
        double len = d.norm();
        return (len > 1e-15) ? Vec3(d / len) : Vec3::Zero();
    }

    // Classification helpers
    bool IsPeptideBond() const { return order == BondOrder::Peptide; }
    bool IsPeptideCO() const { return category == BondCategory::PeptideCO; }
    bool IsBackbone() const {
        return category == BondCategory::PeptideCO ||
               category == BondCategory::PeptideCN ||
               category == BondCategory::BackboneOther;
    }
    bool IsAromatic() const { return category == BondCategory::Aromatic; }
};

}  // namespace nmr
