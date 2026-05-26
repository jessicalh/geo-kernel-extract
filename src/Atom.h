#pragma once
//
// Atom: identity and topology only — no position. Position lives on
// ProteinConformation.
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

    double CovalentRadius() const { return CovalentRadiusForElement(element); }
    double Electronegativity() const { return ElectronegativityForElement(element); }
    bool   IsHBondDonorElement() const {
        return element == Element::N || element == Element::O;
    }
    bool   IsHBondAcceptorElement() const {
        return element == Element::N || element == Element::O;
    }

    static std::unique_ptr<Atom> Create(Element elem);
    static std::unique_ptr<Atom> Create(const std::string& elementSymbol);
};

}  // namespace nmr
