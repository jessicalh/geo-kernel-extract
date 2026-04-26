#include "AtomReference.h"
#include "Atom.h"
#include "Residue.h"
#include "Protein.h"

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

}  // namespace nmr
