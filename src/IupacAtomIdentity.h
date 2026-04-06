#pragma once
//
// IupacAtomIdentity: ground truth validator for IUPAC protein atom naming.
//
// This module answers: "What IS atom CA in alanine?" Not by convention,
// not by parsing -- by the IUPAC definition of the amino acid.
//
// For each atom in each amino acid, this provides:
//   - Element (redundant with name, but explicit)
//   - Whether it's backbone or sidechain
//   - Whether it's a ring member
//
// This is the VALIDATOR. When we load a PDB and assign atom names, this
// module checks: "You say this atom is CA. Is it element C? Is it backbone?"
// If not, something is wrong.
//
// Usage:
//   const auto* id = IupacAtomIdentity::Lookup("ALA", "CA");
//   id->element;       // Element::C
//   id->is_backbone;   // true
//   id->is_ring_member;// false
//
//   std::string err = IupacAtomIdentity::Validate("ALA", "CA", Element::C, true);
//   // err is empty (valid)
//

#include "Types.h"
#include <string>
#include <vector>

namespace nmr {

struct IupacAtom {
    const char* name;
    Element element;
    bool is_backbone;
    bool is_ring_member;
};

class IupacAtomIdentity {
public:
    // Look up an atom's IUPAC identity by residue code + atom name.
    // Returns nullptr if the atom is not part of this residue's canonical set.
    static const IupacAtom* Lookup(const std::string& residue_code,
                                    const std::string& atom_name);

    // Get all atoms for a residue.
    // Returns an empty vector if the residue is unknown.
    static const std::vector<IupacAtom>& AtomsForResidue(const std::string& residue_code);

    // Validate: does this atom name + element + backbone flag match IUPAC?
    // Returns empty string if OK, error message if mismatch.
    static std::string Validate(const std::string& residue_code,
                                 const std::string& atom_name,
                                 Element element,
                                 bool claimed_backbone);
};

}  // namespace nmr
