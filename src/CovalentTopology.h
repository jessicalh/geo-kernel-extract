#pragma once
// Covalent topology is resolved once from one coordinate set and then
// shared by every conformation. That assumes the covalent graph is
// invariant across the trajectory.

#include "Types.h"
#include "Bond.h"
#include <vector>
#include <memory>

namespace nmr {

class Atom;
class Ring;
class Residue;

// Disulfide pairing recorded by a loader from GROMACS bonded terms,
// not inferred from distance. When a caller invokes OverrideDisulfides,
// an empty pair list means the source says there are zero disulfides.
struct DisulfidePair {
    size_t residue_a = 0;          // 0-based residue index in Protein
    size_t residue_b = 0;
    size_t atom_index_sg_a = 0;    // protein-local atom index
    size_t atom_index_sg_b = 0;
};

class CovalentTopology {
public:
    // Requires: atoms with canonical names, residues with backbone
    // indices cached.
    //
    // Positions are used for bond detection and nearest-heavy-atom
    // hydrogen parent assignment; they are not stored.
    static std::unique_ptr<CovalentTopology> Resolve(
        const std::vector<std::unique_ptr<Atom>>& atoms,
        const std::vector<Residue>& residues,
        const std::vector<Vec3>& positions,
        double bond_tolerance = 0.4);

    // Tags bonds whose endpoints are both in an aromatic ring.
    // Disulfide bonds keep their category even if their endpoints
    // appear in the supplied ring set.
    void TagAromaticBonds(
        const std::vector<std::unique_ptr<Ring>>& rings);

    size_t BondCount() const { return bonds_.size(); }
    const Bond& BondAt(size_t i) const { return bonds_[i]; }
    const std::vector<Bond>& Bonds() const { return bonds_; }

    const std::vector<size_t>& BondIndicesFor(size_t atom_index) const {
        return bond_indices_[atom_index];
    }

    // H parent: nearest bonded heavy atom. SIZE_MAX if not hydrogen.
    size_t HydrogenParentOf(size_t atom_index) const {
        return h_parent_[atom_index];
    }

    // Apply caller-supplied disulfide authority. Each listed pair is
    // set to BondCategory::Disulfide, regardless of Resolve() output.
    // Existing Disulfide tags not in `pairs` are reset to SidechainOther.
    //
    // An empty `pairs` vector is not a no-op once this function is
    // called. It means "authority says zero disulfides", so every
    // geometry-derived Disulfide tag is demoted.
    //
    // If a pair's bond is missing from the geometric topology
    // (OpenBabel didn't detect it), the bond is force-added with
    // BondOrder::Single + BondCategory::Disulfide and the per-atom
    // bond_indices_ updated.
    //
    // Returns an error string on validation failure; empty on success.
    // Diagnostics about geometric/chemistry disagreement are emitted
    // via OperationLog.
    std::string OverrideDisulfides(const std::vector<DisulfidePair>& pairs);

private:
    std::vector<Bond> bonds_;
    std::vector<std::vector<size_t>> bond_indices_;  // per atom
    std::vector<size_t> h_parent_;                    // per atom
};

}  // namespace nmr
