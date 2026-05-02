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

// Authoritative disulfide pairing fact. Recorded by the loader from
// the TPR bonded-interaction list (where pdb2gmx specbond.cpp wrote
// the SG-SG bond it decided to add) plus the readback block's CYX
// residue identity. NOT geometric inference: this is "GROMACS told us
// these two cysteines are disulfide-bonded."
//
// Lives on `LegacyAmberInvariants` (carried alongside FF-numerical
// invariants from BuildProtein → FinalizeConstruction). Applied by
// `CovalentTopology::OverrideDisulfides` after geometric bond
// detection. Empty default for non-trajectory load paths (PDB +
// ff14SB flat table, raw PDB) where no upstream authority exists.
struct DisulfidePair {
    size_t residue_a = 0;          // 0-based residue index in Protein
    size_t residue_b = 0;
    size_t atom_index_sg_a = 0;    // protein-local atom index
    size_t atom_index_sg_b = 0;
};

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

    // Apply the authoritative disulfide pairing as recorded by GROMACS
    // (TPR bonded list + topol.top rtp). Each pair is set to
    // BondCategory::Disulfide explicitly, regardless of what the
    // geometric categorization in Resolve() inferred. Any bond
    // currently tagged Disulfide that is NOT in `pairs` is reset to
    // BondCategory::SidechainOther (geometry detected an SG-SG bond
    // that chemistry says is not a disulfide — extremely rare in
    // standard MD, but the authority must win).
    //
    // If a pair's bond is missing from the geometric topology
    // (OpenBabel didn't detect it), the bond is force-added with
    // BondOrder::Single + BondCategory::Disulfide and the per-atom
    // bond_indices_ updated.
    //
    // Returns an error string on validation failure (atom indices out
    // of range, etc.); empty on success. Diagnostics about
    // geometric/chemistry disagreement are emitted via OperationLog.
    std::string OverrideDisulfides(const std::vector<DisulfidePair>& pairs);

private:
    std::vector<Bond> bonds_;
    std::vector<std::vector<size_t>> bond_indices_;  // per atom
    std::vector<size_t> h_parent_;                    // per atom
};

}  // namespace nmr
