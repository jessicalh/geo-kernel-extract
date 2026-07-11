#include "LarsenHBondGeometryCommon.h"

#include "Bond.h"
#include "LegacyAmberTopology.h"
#include "Protein.h"
#include "Residue.h"
#include "SemanticEnums.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace nmr::larsen_hbond_geometry {
namespace {

bool IsHydroxylOxygen(const Protein& protein,
                      std::size_t oxygen_idx,
                      std::size_t& bonded_hydrogen_out) {
    bonded_hydrogen_out = Residue::NONE;
    const Atom& oxygen = protein.AtomAt(oxygen_idx);
    if (oxygen.element != Element::O) return false;

    const LegacyAmberTopology& topology = protein.LegacyAmber();
    for (std::size_t bond_idx : oxygen.bond_indices) {
        const Bond& bond = protein.BondAt(bond_idx);
        const std::size_t other = bond.atom_index_a == oxygen_idx
            ? bond.atom_index_b : bond.atom_index_a;
        if (protein.AtomAt(other).element != Element::H) continue;
        const PolarHKind kind = topology.SemanticAt(other).polar_h;
        if (kind == PolarHKind::HydroxylOH_Aliphatic ||
            kind == PolarHKind::HydroxylOH_Aromatic) {
            bonded_hydrogen_out = other;
            return true;
        }
    }
    return false;
}

std::size_t FirstBondedElement(const Protein& protein,
                               std::size_t atom_idx,
                               Element wanted,
                               std::size_t excluded = SIZE_MAX) {
    for (std::size_t bond_idx : protein.AtomAt(atom_idx).bond_indices) {
        const Bond& bond = protein.BondAt(bond_idx);
        const std::size_t other = bond.atom_index_a == atom_idx
            ? bond.atom_index_b : bond.atom_index_a;
        if (other == excluded) continue;
        if (protein.AtomAt(other).element == wanted) return other;
    }
    return Residue::NONE;
}

}  // namespace

std::optional<AcceptorTriple> ClassifyAcceptor(
        const Protein& protein,
        std::size_t oxygen_atom_idx) {
    const AtomSemanticTable& sem =
        protein.LegacyAmber().SemanticAt(oxygen_atom_idx);
    if (sem.element != Element::O) return std::nullopt;

    AcceptorTriple triple;
    triple.O_idx = oxygen_atom_idx;

    // Carboxylate is checked before backbone carbonyl: a terminal backbone
    // oxygen carries CarbonylOxygen identity but carboxylate chemistry.
    if (sem.IsSidechainCarboxylateOxygen()) {
        triple.class_ = HBondAcceptorClass::CarboxylateOxygen;
        triple.C_idx = FirstBondedElement(
            protein, oxygen_atom_idx, Element::C);
        if (triple.C_idx == Residue::NONE) return std::nullopt;

        for (std::size_t bond_idx :
             protein.AtomAt(triple.C_idx).bond_indices) {
            const Bond& bond = protein.BondAt(bond_idx);
            const std::size_t other = bond.atom_index_a == triple.C_idx
                ? bond.atom_index_b : bond.atom_index_a;
            if (other == oxygen_atom_idx) continue;
            const AtomSemanticTable& other_sem =
                protein.LegacyAmber().SemanticAt(other);
            if (other_sem.element == Element::O &&
                other_sem.planar_group == PlanarGroupKind::Carboxylate) {
                triple.third_idx = other;
                break;
            }
        }
        if (triple.third_idx == Residue::NONE) return std::nullopt;
        return triple;
    }

    if (sem.IsBackboneCarbonylOxygen()) {
        triple.class_ = HBondAcceptorClass::BackboneCarbonyl;
        const std::size_t residue_idx =
            protein.AtomAt(oxygen_atom_idx).residue_index;
        const Residue& residue = protein.ResidueAt(residue_idx);
        if (residue.C == Residue::NONE) return std::nullopt;
        triple.C_idx = residue.C;
        if (const auto successor = protein.BackboneSuccessor(residue_idx)) {
            triple.i_plus_1_residue_idx = *successor;
            triple.third_idx = protein.ResidueAt(*successor).N;
        }
        return triple;
    }

    if (sem.IsSidechainAmideOxygen()) {
        triple.class_ = HBondAcceptorClass::SidechainCarbonyl;
        triple.C_idx = FirstBondedElement(
            protein, oxygen_atom_idx, Element::C);
        if (triple.C_idx == Residue::NONE) return std::nullopt;

        for (std::size_t bond_idx :
             protein.AtomAt(triple.C_idx).bond_indices) {
            const Bond& bond = protein.BondAt(bond_idx);
            const std::size_t other = bond.atom_index_a == triple.C_idx
                ? bond.atom_index_b : bond.atom_index_a;
            if (other == oxygen_atom_idx) continue;
            const AtomSemanticTable& other_sem =
                protein.LegacyAmber().SemanticAt(other);
            if (other_sem.element == Element::N &&
                other_sem.planar_group == PlanarGroupKind::SidechainAmide) {
                triple.third_idx = other;
                break;
            }
        }
        if (triple.third_idx == Residue::NONE) return std::nullopt;
        return triple;
    }

    std::size_t hydroxyl_h = Residue::NONE;
    if (IsHydroxylOxygen(protein, oxygen_atom_idx, hydroxyl_h)) {
        triple.class_ = HBondAcceptorClass::HydroxylOxygen;
        triple.C_idx = FirstBondedElement(
            protein, oxygen_atom_idx, Element::C);
        if (triple.C_idx == Residue::NONE) return std::nullopt;
        triple.third_idx = hydroxyl_h;
        return triple;
    }

    return std::nullopt;
}

std::size_t BondedHeavyAtomOfHydrogen(
        const Protein& protein,
        std::size_t hydrogen_atom_idx) {
    if (hydrogen_atom_idx >= protein.AtomCount() ||
        protein.AtomAt(hydrogen_atom_idx).element != Element::H) {
        return SIZE_MAX;
    }

    std::size_t heavy = SIZE_MAX;
    for (std::size_t bond_idx :
         protein.AtomAt(hydrogen_atom_idx).bond_indices) {
        const Bond& bond = protein.BondAt(bond_idx);
        const std::size_t other = bond.atom_index_a == hydrogen_atom_idx
            ? bond.atom_index_b : bond.atom_index_a;
        if (protein.AtomAt(other).element == Element::H) continue;
        if (heavy != SIZE_MAX && heavy != other) return SIZE_MAX;
        heavy = other;
    }
    return heavy;
}

double AngleDegrees(const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 ba = a - b;
    const Vec3 bc = c - b;
    const double denominator = ba.norm() * bc.norm();
    if (!std::isfinite(denominator) ||
        denominator <= std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double cosine = std::clamp(ba.dot(bc) / denominator, -1.0, 1.0);
    return std::acos(cosine) * (180.0 / M_PI);
}

}  // namespace nmr::larsen_hbond_geometry

