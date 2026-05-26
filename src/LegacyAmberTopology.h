#pragma once
// LegacyAmberTopology holds the covalent graph plus invariant
// force-field facts supplied by the loader. Empty vectors and zero
// scalars mean the load path did not provide that field.

#include "ProteinTopology.h"
#include "CovalentTopology.h"
#include "RingTopology.h"
#include "SemanticEnums.h"
#include "Residue.h"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class Atom;

// Constructor payload for invariant force-field fields. Default-empty
// is the "not provided" state.
struct LegacyAmberInvariants {
    // Per-atom invariant FF-numerical fields (vectors, when populated,
    // are size = atom_count; empty otherwise).
    std::vector<double> mass;                   // Da
    std::vector<int>    ff_atom_type_index;     // index into the FF atom-type table
    std::vector<int>    ptype;                  // GROMACS ParticleType enum (cast to int)
    std::vector<std::string> atomtype_string;   // FF atom-type name (e.g. "CT", "N", "H1")

    // Exclusion lists from the FF (1-2, 1-3, 1-4 within nrexcl). Each
    // inner vector lists protein-local atom indices excluded from the
    // index's nonbonded interactions. Empty when not provided.
    std::vector<std::vector<int>> exclusions;

    // Force-field-wide constants (provenance for nonbonded calculators).
    double fudge_qq = 1.0;                      // 1-4 Coulomb scaling
    double rep_pow  = 12.0;                     // LJ repulsion exponent
    int    atnr     = 0;                        // FF atom-type table dimension

    // Number of non-perturbed interactions per InteractionFunction kind.
    std::array<int, 256> num_non_perturbed = {};

    // Disulfide pairing from GROMACS bonded terms, applied through
    // CovalentTopology::OverrideDisulfides when has_disulfide_authority
    // is true.
    //
    // Empty `disulfide_pairs` is meaningful only in conjunction with
    // `has_disulfide_authority`:
    //   - has_disulfide_authority=false (default, non-trajectory path):
    //       no upstream authority. Geometric SG-SG inference in
    //       CovalentTopology::Resolve stays as source of truth.
    //   - has_disulfide_authority=true, disulfide_pairs.empty():
    //       authority says zero disulfides. Override runs with empty
    //       authority list — any geometric Disulfide tags get demoted
    //       to SidechainOther (warned, but the authority wins).
    //   - has_disulfide_authority=true, disulfide_pairs non-empty:
    //       authority's specific pair list. Override applies and
    //       validates against geometric inference.
    std::vector<DisulfidePair> disulfide_pairs;
    bool has_disulfide_authority = false;
};


class LegacyAmberTopology final : public ProteinTopology {
public:
    LegacyAmberTopology(size_t atom_count,
                        size_t residue_count,
                        std::unique_ptr<CovalentTopology> bonds,
                        LegacyAmberInvariants invariants = {},
                        std::vector<AtomSemanticTable> atom_semantic = {},
                        std::unique_ptr<RingTopology> rings = {});

    ProteinTopologyKind Kind() const override {
        return ProteinTopologyKind::LegacyAmber;
    }
    std::string_view Name() const override { return "LegacyAmberTopology"; }
    size_t AtomCount() const override { return atom_count_; }
    size_t ResidueCount() const override { return residue_count_; }

    const CovalentTopology& Bonds() const { return *bonds_; }

    size_t BondCount() const { return bonds_->BondCount(); }
    const Bond& BondAt(size_t i) const { return bonds_->BondAt(i); }
    const std::vector<Bond>& BondList() const { return bonds_->Bonds(); }
    const std::vector<size_t>& BondIndicesFor(size_t atom_index) const {
        return bonds_->BondIndicesFor(atom_index);
    }
    size_t HydrogenParentOf(size_t atom_index) const {
        return bonds_->HydrogenParentOf(atom_index);
    }

    // rings_ is always non-null. It may contain zero rings when the
    // atom-semantic substrate is absent.
    const RingTopology& Rings() const { return *rings_; }

    size_t AromaticRingCount() const { return rings_->AromaticCount(); }
    const Ring& AromaticRingAt(size_t i) const {
        return rings_->AromaticAt(i);
    }
    const std::vector<std::unique_ptr<Ring>>& AromaticRingList() const {
        return rings_->Aromatic();
    }

    size_t SaturatedRingCount() const { return rings_->SaturatedCount(); }
    const Ring& SaturatedRingAt(size_t i) const {
        return rings_->SaturatedAt(i);
    }
    const std::vector<std::unique_ptr<Ring>>& SaturatedRingList() const {
        return rings_->Saturated();
    }

    // Empty / zero when the load path didn't carry FF data. Calculators
    // gate on `!Mass().empty()` (etc.) when needed.

    const std::vector<double>& Mass() const { return mass_; }
    const std::vector<int>& FfAtomTypeIndex() const { return ff_atom_type_index_; }
    const std::vector<int>& Ptype() const { return ptype_; }
    const std::vector<std::string>& AtomtypeString() const { return atomtype_string_; }
    const std::vector<std::vector<int>>& Exclusions() const { return exclusions_; }

    double FudgeQq() const { return fudge_qq_; }
    double RepPow() const { return rep_pow_; }
    int    Atnr() const { return atnr_; }
    const std::array<int, 256>& NumNonPerturbed() const { return num_non_perturbed_; }

    // Empty atom_semantic_ means no substrate was populated; when
    // populated, size equals AtomCount().

    bool HasAtomSemantic() const { return !atom_semantic_.empty(); }

    // FATAL+abort if atom_semantic_ is empty or atom_index is out of range.
    const AtomSemanticTable& SemanticAt(size_t atom_index) const;

    const std::vector<AtomSemanticTable>& AtomSemantic() const {
        return atom_semantic_;
    }

    AtomMechanicalIdentity IdentityAt(size_t atom_index) const {
        const AtomSemanticTable& sem = SemanticAt(atom_index);
        return AtomMechanicalIdentity{
            sem.element, sem.locant, sem.branch, sem.di_index,
            sem.backbone_role
        };
    }

    // `residues` is the protein's residue vector (the topology does
    // not own residues; the caller passes them in for atom-index
    // lookup). Empty result means no match. FATAL+abort if
    // residue_index is out of range.
    std::vector<size_t> ResidueAtomsWithIdentity(
        size_t residue_index,
        const AtomMechanicalIdentity& identity,
        const std::vector<Residue>& residues) const;

    // Returns Residue::NONE if no atom in the residue has the role.
    // FATAL+abort if residue_index is out of range.
    size_t AtomWithRole(size_t residue_index,
                        BackboneRole role,
                        const std::vector<Residue>& residues) const;

private:
    size_t atom_count_ = 0;
    size_t residue_count_ = 0;
    std::unique_ptr<CovalentTopology> bonds_;

    std::vector<double> mass_;
    std::vector<int>    ff_atom_type_index_;
    std::vector<int>    ptype_;
    std::vector<std::string> atomtype_string_;
    std::vector<std::vector<int>> exclusions_;

    double fudge_qq_ = 1.0;
    double rep_pow_  = 12.0;
    int    atnr_     = 0;
    std::array<int, 256> num_non_perturbed_ = {};

    // Per-atom typed chemistry-substrate, composed at FinalizeConstruction
    // via ComposeAtomSemantic. Empty when the load path delivers stub
    // atoms with no names; non-empty has size == atom_count_.
    std::vector<AtomSemanticTable> atom_semantic_;

    // Aromatic + saturated rings. Built by RingTopology::ConstructFromSubstrate
    // at FinalizeConstruction. Always non-null after construction; empty
    // (zero rings) is the legitimate stub-fixture signal.
    std::unique_ptr<RingTopology> rings_;
};


// Compose the per-atom AtomSemanticTable vector from canonical atom
// names, heavy-atom parent names, and generated substrate tables.
//
// If no atom has a PDB atom name, returns an empty vector.
//
// Fails fatally on:
//   - chain atom whose (residue, variant, identity) has no LookupBy match.
//   - cap-only atom whose (terminal_state, identity) has no LookupCap match.
//   - named AminoAcid::Unknown residue.
std::vector<AtomSemanticTable> ComposeAtomSemantic(
    const std::vector<std::unique_ptr<Atom>>& atoms,
    const std::vector<Residue>& residues,
    const CovalentTopology& bonds);

}  // namespace nmr
