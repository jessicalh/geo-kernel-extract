// IupacAtomIdentity: IUPAC ground truth for protein atom naming.
//
// Data from IUPAC-IUB Commission on Biochemical Nomenclature (1970),
// as updated in Markley et al. (1998) Pure Appl. Chem. 70, 117-142.
//
// Built from the AminoAcidType table (the single authority for amino acid
// chemistry). Ring membership is determined from the ring definitions
// in AminoAcidType.

#include "IupacAtomIdentity.h"
#include "AminoAcidType.h"
#include <map>
#include <set>

namespace nmr {

// Build the IUPAC atom tables on first access (Meyers singleton).
static const std::map<std::string, std::vector<IupacAtom>>& GetIupacTables() {
    static const auto tables = []() {
        std::map<std::string, std::vector<IupacAtom>> t;

        for (const auto& type : AllAminoAcidTypes()) {
            if (type.index == AminoAcid::Unknown) continue;

            // Collect ring member atom names for this amino acid
            std::set<std::string> ring_atoms;
            for (const auto& ring : type.rings)
                for (const char* rn : ring.atom_names)
                    ring_atoms.insert(rn);

            std::vector<IupacAtom> atoms;
            atoms.reserve(type.atoms.size());
            for (const auto& a : type.atoms) {
                bool is_ring = ring_atoms.count(a.name) > 0;
                atoms.push_back({a.name, a.element, a.is_backbone, is_ring});
            }
            t[type.three_letter_code] = std::move(atoms);
        }

        return t;
    }();
    return tables;
}


const IupacAtom* IupacAtomIdentity::Lookup(const std::string& residue_code,
                                             const std::string& atom_name) {
    const auto& tables = GetIupacTables();
    auto it = tables.find(residue_code);
    if (it == tables.end()) return nullptr;
    for (const auto& a : it->second)
        if (atom_name == a.name) return &a;
    return nullptr;
}


const std::vector<IupacAtom>& IupacAtomIdentity::AtomsForResidue(
        const std::string& residue_code) {
    static const std::vector<IupacAtom> empty;
    const auto& tables = GetIupacTables();
    auto it = tables.find(residue_code);
    return (it != tables.end()) ? it->second : empty;
}


std::string IupacAtomIdentity::Validate(const std::string& residue_code,
                                          const std::string& atom_name,
                                          Element element,
                                          bool claimed_backbone) {
    const auto* atom = Lookup(residue_code, atom_name);
    if (!atom)
        return "Atom '" + atom_name + "' is not part of " + residue_code + " in IUPAC";

    if (atom->element != element)
        return "Atom '" + atom_name + "' in " + residue_code +
               ": element mismatch (IUPAC=" + SymbolForElement(atom->element) +
               " actual=" + SymbolForElement(element) + ")";

    if (atom->is_backbone != claimed_backbone)
        return "Atom '" + atom_name + "' in " + residue_code +
               ": backbone flag mismatch (IUPAC=" +
               std::string(atom->is_backbone ? "true" : "false") +
               " actual=" + std::string(claimed_backbone ? "true" : "false") + ")";

    return "";  // valid
}

}  // namespace nmr
