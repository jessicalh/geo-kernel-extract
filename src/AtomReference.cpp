#include "AtomReference.h"
#include "Atom.h"
#include "AminoAcidType.h"
#include "Residue.h"
#include "Protein.h"
#include "OperationLog.h"

#include <cstdio>
#include <cstdlib>

namespace nmr {

AtomReference MakeAtomReference(const Protein& protein, size_t atom_index) {
    const Atom& a = protein.AtomAt(atom_index);
    const Residue& r = protein.ResidueAt(a.residue_index);
    AtomReference ref;
    ref.residue_type     = r.type;
    ref.residue_position = r.sequence_number;
    ref.chain_id         = r.chain_id;
    ref.atom_name        = a.iupac_name;
    return ref;
}

std::unordered_map<AtomReference, size_t> BuildAtomReferenceMap(
        const Protein& protein) {
    std::unordered_map<AtomReference, size_t> map;
    map.reserve(protein.AtomCount());
    for (size_t ai = 0; ai < protein.AtomCount(); ++ai) {
        map.emplace(MakeAtomReference(protein, ai), ai);
    }
    return map;
}

std::ostream& operator<<(std::ostream& os, const AtomReference& ref) {
    os << ThreeLetterCodeForAminoAcid(ref.residue_type)
       << '-' << ref.residue_position
       << ':' << ref.atom_name;
    if (!ref.chain_id.empty()) {
        os << " [chain " << ref.chain_id << ']';
    }
    return os;
}


AtomLocator MakeAtomLocator(const Protein& protein, size_t atom_index) {
    const Atom& a = protein.AtomAt(atom_index);
    AtomLocator loc;
    loc.residue_index = a.residue_index;
    loc.atom_name     = a.iupac_name;
    return loc;
}

std::unordered_map<AtomLocator, size_t> BuildAtomLocatorMap(
        const Protein& protein) {
    std::unordered_map<AtomLocator, size_t> map;
    map.reserve(protein.AtomCount());
    for (size_t ai = 0; ai < protein.AtomCount(); ++ai) {
        AtomLocator loc = MakeAtomLocator(protein, ai);
        auto [it, inserted] = map.emplace(loc, ai);
        if (!inserted) {
            // Duplicate locator key: two atoms in the same residue share
            // the same IUPAC name. This is a topology bug — surface it
            // loud rather than letting "first wins" silently corrupt
            // cross-protein matching. Includes a useful pointer to the
            // residue's three-letter code so loaders can locate the case.
            const Residue& res = protein.ResidueAt(loc.residue_index);
            const char* res_code = res.AminoAcidInfo().three_letter_code;
            std::string detail =
                "BuildAtomLocatorMap: duplicate locator at residue " +
                std::string(res_code) + std::to_string(res.sequence_number) +
                " atom " + loc.atom_name.AsString() +
                " — atoms " + std::to_string(it->second) +
                " and " + std::to_string(ai) +
                " share (residue_index=" + std::to_string(loc.residue_index) +
                ", atom_name=" + loc.atom_name.AsString() + ")";
            OperationLog::Error("AtomReference", detail);
            fprintf(stderr,
                "FATAL: %s\n"
                "  Two atoms cannot share an AtomLocator key. Investigate\n"
                "  the loader or NamingRegistry path that produced these\n"
                "  atom names — silent first-wins matching across two\n"
                "  proteins would produce wrong DFT delta tensors.\n",
                detail.c_str());
            std::abort();
        }
    }
    return map;
}

std::ostream& operator<<(std::ostream& os, const AtomLocator& loc) {
    os << "ri=" << loc.residue_index << ':' << loc.atom_name;
    return os;
}

}  // namespace nmr
