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

#include "Types.h"
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
