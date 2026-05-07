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
// detection. Non-trajectory load paths do not call the override; for
// those paths geometric SG-SG inference remains the source of truth.
// When an upstream authority exists, an empty pair list is meaningful:
// it means the authority says there are zero disulfides.
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
    // indices cached. These are the output of symbolic topology
    // resolution (layers 1-2 of the loading pipeline).
    //
    // Positions are used ONLY for bond detection (covalent radius
    // distance check) and H parent assignment (nearest bonded heavy
    // atom). After Resolve(), the positions are not stored.
    //
    // Bundle C / Slice B (2026-05-07): Resolve no longer takes a
    // rings parameter. The aromatic-bond categorisation overlay
    // moved to TagAromaticBonds(rings), which is called AFTER ring
    // construction (which is now substrate-driven and must run after
    // ComposeAtomSemantic). Pre-rings Resolve produces the bond
    // graph + non-aromatic categorisation; TagAromaticBonds is the
    // post-pass overlay that mirrors OverrideDisulfides's shape.
    static std::unique_ptr<CovalentTopology> Resolve(
        const std::vector<std::unique_ptr<Atom>>& atoms,
        const std::vector<Residue>& residues,
        const std::vector<Vec3>& positions,
        double bond_tolerance = 0.4);

    // Tag bonds whose both endpoints sit in any aromatic ring as
    // BondCategory::Aromatic / BondOrder::Aromatic. Lifted from
    // Resolve()'s pre-Bundle-C aromatic branch (lines 87-90 +
    // 143-146 of the pre-split CovalentTopology.cpp). Mirrors
    // OverrideDisulfides's shape: a post-Resolve overlay that
    // applies authoritative chemistry on top of geometry-derived
    // categorisation.
    //
    // Precedence note: Resolve's geometric categorisation puts S-S
    // bonds in BondCategory::Disulfide first; this overlay does NOT
    // demote disulfides to aromatic even if both S atoms were ever
    // in an aromatic ring (which is biologically impossible — no
    // aromatic ring contains an S that bonds another S). The
    // disulfide-precedence-over-aromatic behaviour is preserved
    // because the overlay only rewrites bonds whose category is
    // anything OTHER than Disulfide. See the matching note in the
    // pre-split Resolve at the original lines 137-141.
    void TagAromaticBonds(
        const std::vector<std::unique_ptr<Ring>>& rings);

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
    // (TPR bonded list + topol.top rtp). Call this only when an upstream
    // authority is present. Each pair is set to BondCategory::Disulfide
    // explicitly, regardless of what the geometric categorization in
    // Resolve() inferred. Any bond currently tagged Disulfide that is
    // NOT in `pairs` is reset to BondCategory::SidechainOther.
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
