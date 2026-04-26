#pragma once
//
// AtomReference: a typed cross-protein atom identity.
//
// Identifies "the same atom" across different Protein instances — WT vs
// mutant, frame N vs frame M of an MD trajectory, classical kernel output
// vs DFT shielding output. Within a single Protein, atom_index is the
// natural identifier; AtomReference is for matching atoms BETWEEN proteins
// where atom_index is not stable (mutation sites have different atom counts;
// alternate loaders may produce different orderings).
//
// Composition:
//   residue_type      AminoAcid enum — typed, never a string
//   residue_position  int sequence number (PDB column 23-26)
//   chain_id          std::string (multi-chain proteins)
//   atom_name         IupacAtomName (typed atom-name model)
//
// Equality compares all four. Two AtomReferences compare equal iff they
// refer to the structurally-equivalent atom across protein instances.
//
// Hash and ordering are provided so AtomReference works as a
// std::unordered_map / std::map key. Typical use:
//
//   auto wt_map  = BuildAtomReferenceMap(wt_protein);   // map<AtomReference, size_t>
//   auto mut_map = BuildAtomReferenceMap(mut_protein);
//   for (size_t wi = 0; wi < wt_protein.AtomCount(); ++wi) {
//       AtomReference ref = MakeAtomReference(wt_protein, wi);
//       auto it = mut_map.find(ref);
//       if (it != mut_map.end()) {
//           // Same atom in both proteins — compare DFT shielding, etc.
//       }
//       // Mutation-site sidechain atoms naturally fall through here.
//   }
//
// No string compare in the calculator. The matching is typed object
// equality — same boundary discipline as IupacAtomName, extended to the
// (residue, atom) tuple.
//

#include "Types.h"
#include "IupacAtomName.h"

#include <cstddef>
#include <functional>
#include <ostream>
#include <string>
#include <unordered_map>

namespace nmr {

class Atom;
class Residue;
class Protein;

struct AtomReference {
    AminoAcid     residue_type     = AminoAcid::Unknown;
    int           residue_position = 0;            // PDB sequence number
    std::string   chain_id;                         // empty for single-chain
    IupacAtomName atom_name;

    bool operator==(const AtomReference& o) const {
        return residue_type     == o.residue_type
            && residue_position == o.residue_position
            && chain_id         == o.chain_id
            && atom_name        == o.atom_name;
    }
    bool operator!=(const AtomReference& o) const { return !(*this == o); }

    // Ordering for std::map keys.
    bool operator<(const AtomReference& o) const {
        if (chain_id != o.chain_id) return chain_id < o.chain_id;
        if (residue_position != o.residue_position) return residue_position < o.residue_position;
        if (residue_type != o.residue_type) return residue_type < o.residue_type;
        return atom_name < o.atom_name;
    }
};

// Build an AtomReference for a specific atom in a Protein.
// Reads residue_type / residue_position / chain_id from the atom's residue
// and atom_name from the atom itself.
AtomReference MakeAtomReference(const Protein& protein, size_t atom_index);

// Build a map from AtomReference to atom_index for a whole Protein.
// O(N) construction; O(1) average lookup. Use when comparing two proteins
// atom-by-atom (DFT compare, mutation delta, frame-to-frame identity).
std::unordered_map<AtomReference, size_t> BuildAtomReferenceMap(
    const Protein& protein);

// Stream output for diagnostics: e.g. "PHE-4:HB3" or "GLY-10:HA2 [chain A]".
std::ostream& operator<<(std::ostream& os, const AtomReference& ref);

}  // namespace nmr

namespace std {
template <>
struct hash<nmr::AtomReference> {
    size_t operator()(const nmr::AtomReference& r) const noexcept {
        // Combine field hashes via the boost-style mixer (golden ratio).
        auto mix = [](size_t seed, size_t v) {
            return seed ^ (v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
        };
        size_t h = std::hash<int>{}(static_cast<int>(r.residue_type));
        h = mix(h, std::hash<int>{}(r.residue_position));
        h = mix(h, std::hash<std::string>{}(r.chain_id));
        h = mix(h, std::hash<nmr::IupacAtomName>{}(r.atom_name));
        return h;
    }
};
}  // namespace std
