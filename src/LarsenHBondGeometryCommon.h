#pragma once
//
// Shared, typed geometry substrate for the Larsen H-bond calculators.
//
// Acceptor classification is chemistry/topology work: it reads the
// AtomSemanticTable and the explicit covalent bond graph.  Keeping it here
// gives LarsenHBondShieldingResult and LarsenSidechainDonorAuditResult one
// authoritative classification without reopening the atom-name string
// boundary.

#include "LarsenHBondGrid.h"

#include <cstddef>
#include <optional>

namespace nmr {

class Protein;

namespace larsen_hbond_geometry {

// Resolved acceptor frame (O, bonded heavy atom, dihedral third atom) plus
// the downstream residue used only by Larsen's backbone-carbonyl 2-degree
// readouts.  SIZE_MAX means the corresponding frame/routing atom is absent.
struct AcceptorTriple {
    HBondAcceptorClass class_ = HBondAcceptorClass::BackboneCarbonyl;
    std::size_t O_idx = SIZE_MAX;
    std::size_t C_idx = SIZE_MAX;
    std::size_t third_idx = SIZE_MAX;
    std::size_t i_plus_1_residue_idx = SIZE_MAX;
};

// Classify a typed oxygen acceptor and resolve its frame anchors.  Returns
// nullopt when the atom is not one of Larsen's supported acceptor classes or
// when class-defining covalent neighbours are absent.  A terminal backbone
// carbonyl is still classified with third_idx == SIZE_MAX so callers can
// emit an explicit missing-frame disposition.
std::optional<AcceptorTriple> ClassifyAcceptor(
    const Protein& protein,
    std::size_t oxygen_atom_idx);

// Resolve the unique bonded heavy atom of a hydrogen from the explicit bond
// graph.  Returns SIZE_MAX for non-hydrogens, no heavy neighbour, or an
// ambiguous multi-heavy-neighbour topology.
std::size_t BondedHeavyAtomOfHydrogen(
    const Protein& protein,
    std::size_t hydrogen_atom_idx);

// Angle a-b-c in degrees, with b as the vertex.  Returns NaN for a
// degenerate/non-finite vector rather than manufacturing a direction.
double AngleDegrees(const Vec3& a, const Vec3& b, const Vec3& c);

}  // namespace larsen_hbond_geometry
}  // namespace nmr

