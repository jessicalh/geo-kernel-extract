#pragma once
//
// Atom: identity only -- no position. Position lives in ProteinConformation.
//
// Each atom knows its element, PDB name, residue membership, and bonds.
// Element-specific properties (covalent radius, electronegativity) are
// non-virtual methods dispatching to the free functions in Types.h.
//
// Flat class: no subclass hierarchy. parent_atom_index is SIZE_MAX for
// non-hydrogen atoms and set to the nearest bonded heavy atom for hydrogens.
//
// IUPAC symbolic topology (locant, prochirality, ring position, ...) is
// carried on the typed `topology` field, populated at protein load time
// by Protein::StampAtomTopology() from the AminoAcidType table.
// Source: Markley et al. 1998, J Biomol NMR 12:1-23.
//

#include "Types.h"
#include "AtomTopology.h"
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class Atom {
public:
    Element       element = Element::Unknown;
    std::string   pdb_atom_name;
    size_t        residue_index = 0;
    std::vector<size_t> bond_indices;

    // For hydrogen atoms: index of the nearest bonded heavy atom.
    // SIZE_MAX when not hydrogen (or when parent not yet assigned).
    size_t parent_atom_index = SIZE_MAX;

    // Aromatic ring memberships: indices into Protein::Rings(). Empty if not
    // in any ring. Some atoms belong to multiple rings (Trp Cδ2, Cε2 are in
    // TRP5, TRP6, AND TRP9). Populated by Protein::StampRingMembership()
    // after DetectAromaticRings() runs.
    std::vector<size_t> ring_indices;

    // Diastereotopic partner: for atoms at a 2/3-numbered prochiral position
    // (HB2 ↔ HB3, HG2 ↔ HG3, ...), the index of the partner atom. SIZE_MAX
    // when not at such a position. Populated by Protein::StampPartnerAtomIndex()
    // from the IUPAC topology.
    size_t partner_atom_index = SIZE_MAX;

    // IUPAC symbolic topology (Markley 1998): locant, branch, diastereotopic
    // index, CIP pro-R/pro-S, planar cis/trans, pseudoatom class, ring
    // position, chi participation. Stamped at protein load; const thereafter.
    // `topology.stamped == false` flags atoms whose (residue_type, atom_name)
    // combination was not recognised — a loud diagnostic was emitted.
    AtomTopology topology;

    // Element properties -- dispatch to free functions in Types.h.
    // No virtual dispatch, no subclass override.
    double CovalentRadius() const { return CovalentRadiusForElement(element); }
    double Electronegativity() const { return ElectronegativityForElement(element); }
    bool   IsHBondDonorElement() const {
        return element == Element::N || element == Element::O;
    }
    bool   IsHBondAcceptorElement() const {
        return element == Element::N || element == Element::O;
    }

    // Factory: creates an Atom with the given element.
    static std::unique_ptr<Atom> Create(Element elem);
    static std::unique_ptr<Atom> Create(const std::string& elementSymbol);
};

}  // namespace nmr
