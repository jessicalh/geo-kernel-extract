#pragma once
//
// CovalentTopology: the geometry-derived connectivity of a protein.
//
// This is where geometric facts (distances between atoms) become
// symbolic facts (bonds, bond categories, H parents). This happens
// once, from one set of positions, and all conformations share it.
//
// The assumption is that covalent topology is invariant across
// conformations — true for classical MD, constant-pH excepted.
//
// Resolve() takes the output of symbolic topology resolution
// (atoms with canonical names, residues with backbone indices,
// rings with typed membership) plus one set of positions.
// It returns a self-contained object.
//
// Protein delegates BondAt/BondCount through this object and
// copies per-atom bond_indices and parent_atom_index back to
// Atom objects for convenient calculator access.
//

#include "Types.h"
#include "Bond.h"
#include <vector>
#include <memory>

namespace nmr {

class Atom;
class Ring;
class Residue;

class CovalentTopology {
public:
    // The geometry→topology boundary. Explicit.
    //
    // Requires: atoms with canonical names, residues with backbone
    // indices cached, rings detected. These are the output of
    // symbolic topology resolution (layers 1-2 of the loading pipeline).
    //
    // Positions are used ONLY for bond detection (covalent radius
    // distance check) and H parent assignment (nearest bonded heavy
    // atom). After Resolve(), the positions are not stored.
    static std::unique_ptr<CovalentTopology> Resolve(
        const std::vector<std::unique_ptr<Atom>>& atoms,
        const std::vector<std::unique_ptr<Ring>>& rings,
        const std::vector<Residue>& residues,
        const std::vector<Vec3>& positions,
        double bond_tolerance = 0.4);

    // Bond access
    size_t BondCount() const { return bonds_.size(); }
    const Bond& BondAt(size_t i) const { return bonds_[i]; }
    const std::vector<Bond>& Bonds() const { return bonds_; }

    // Per-atom connectivity (O(degree) access, same as atom.bond_indices)
    const std::vector<size_t>& BondIndicesFor(size_t atom_index) const {
        return bond_indices_[atom_index];
    }

    // H parent: nearest bonded heavy atom. SIZE_MAX if not hydrogen.
    size_t HydrogenParentOf(size_t atom_index) const {
        return h_parent_[atom_index];
    }

private:
    std::vector<Bond> bonds_;
    std::vector<std::vector<size_t>> bond_indices_;  // per atom
    std::vector<size_t> h_parent_;                    // per atom
};

}  // namespace nmr
