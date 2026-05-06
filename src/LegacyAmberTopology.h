#pragma once
//
// LegacyAmberTopology: explicit topology contract for the existing
// Amber-derived calculator language.
//
// THE topology object. Holds the covalent bond graph plus invariant
// FF-numerical facts the source (TPR, PRMTOP, ...) gives us. All
// fields are plain const data populated at construction. Loaders
// pass a `LegacyAmberInvariants` value-pack when they have FF
// numerical data; loaders that don't (PDB + ff14SB flat table, raw
// PDB) pass nothing and the corresponding fields are empty.
//
// Empty vectors / zero scalars are the legitimate "this load path
// didn't have it" signal; calculators that need a specific field
// check for emptiness directly. There is NO `Optional` slot, NO
// `Attach*Data()` lifecycle, NO `Has*Data()` gate, NO post-hoc
// mutable accessor on the parent. See memory
// `feedback_no_attach_lifecycle_for_invariant_data` and
// `feedback_capture_at_the_boundary`.
//
// Water/ion FF data (SETTLE per water moltype, virtual-site geometry
// for TIP4P-class M-sites) is NOT held here. Waters and ions are
// not in the typed Protein/Residue substrate per the 2026-04-30
// walk-back; calculators that need water FF facts read them at
// per-frame side-channel scope, not from the protein topology.
//

#include "ProteinTopology.h"
#include "CovalentTopology.h"
#include "SemanticEnums.h"
#include "Residue.h"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class Atom;

// Transient bundle of invariant FF-numerical fields passed to the
// `LegacyAmberTopology` constructor. Default-empty: loaders without
// FF data pass `{}`. Fields are moved into the topology and the
// bundle goes out of scope. NOT an object that lives on the topology;
// just constructor-arg packaging.
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

    // Number of non-perturbed interactions per InteractionFunction kind,
    // for FEP-perturbation diagnostics. Surfaced for completeness; we
    // do not currently consume FEP-perturbed TPRs.
    std::array<int, 256> num_non_perturbed = {};

    // Authoritative disulfide pairing fact recorded by pdb2gmx
    // (specbond.cpp + rtp comment line). Applied at FinalizeConstruction
    // time via CovalentTopology::OverrideDisulfides — geometric SG-SG
    // inference stops being the source of truth on the consume side.
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
                        std::vector<AtomSemanticTable> atom_semantic = {});

    ProteinTopologyKind Kind() const override {
        return ProteinTopologyKind::LegacyAmber;
    }
    std::string_view Name() const override { return "LegacyAmberTopology"; }
    size_t AtomCount() const override { return atom_count_; }
    size_t ResidueCount() const override { return residue_count_; }

    // ── Covalent topology ───────────────────────────────────────────
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

    // ── Invariant FF-numerical fields ───────────────────────────────
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

    // ── Per-atom chemistry-substrate (AtomSemanticTable) ────────────
    //
    // `atom_semantic_` is composed at FinalizeConstruction time by
    // ComposeAtomSemantic (free function below). Empty vector means
    // "not populated" — a legitimate signal for stub calculator-physics
    // fixtures whose atoms carry no PDB names. When populated, size
    // equals AtomCount().
    //
    // Whole-row primary surface: calculators read fields off the row
    // returned by SemanticAt(atom_index) directly. Predicate methods
    // live on AtomSemanticTable itself, not here, per the
    // Ring::Intensity() / Ring::IsFused() pattern (objects answer
    // questions about themselves). Per-field shortcuts on this class
    // were considered and rejected because they channel calculator
    // authors toward narrow consumption that hides substrate richness
    // (e.g. IsAromatic(ai) collapses 13 distinct RingPositionLabel
    // values into one boolean).

    bool HasAtomSemantic() const { return !atom_semantic_.empty(); }

    // Returns the AtomSemanticTable for the given atom. FATAL+abort if
    // atom_semantic_ is empty (caller must gate via HasAtomSemantic())
    // or if atom_index is out of range. Mirrors the bonds_ FATAL pattern
    // in the constructor.
    const AtomSemanticTable& SemanticAt(size_t atom_index) const;

    const std::vector<AtomSemanticTable>& AtomSemantic() const {
        return atom_semantic_;
    }

    // Project the identity-tuple subset from the semantic row. Cheap
    // convenience for typed-graph queries below; calculators that need
    // identity directly call this rather than building one by hand.
    AtomMechanicalIdentity IdentityAt(size_t atom_index) const {
        const AtomSemanticTable& sem = SemanticAt(atom_index);
        return AtomMechanicalIdentity{
            sem.element, sem.locant, sem.branch, sem.di_index,
            sem.backbone_role
        };
    }

    // Typed-graph query: return atom indices in `residue_index`'s
    // atom_indices whose mechanical identity equals `identity`.
    // Empty result means no match. Used by Bundle B's typed
    // CacheResidueBackboneIndices for Glycine HA disambiguation
    // (Locant::Alpha + DiastereotopicIndex::Position2) and CB cache
    // (Locant::Beta). Future bundles expand callers.
    //
    // `residues` is the protein's residue vector (the topology does
    // not own residues; the caller passes them in for atom-index
    // lookup). FATAL+abort if residue_index out of range.
    std::vector<size_t> ResidueAtomsWithIdentity(
        size_t residue_index,
        const AtomMechanicalIdentity& identity,
        const std::vector<Residue>& residues) const;

    // Typed backbone-role lookup within a residue. Returns
    // Residue::NONE if no atom has the role. Used by Bundle B's typed
    // CacheResidueBackboneIndices (one call per backbone slot per
    // residue). FATAL+abort if residue_index out of range.
    size_t AtomWithRole(size_t residue_index,
                        BackboneRole role,
                        const std::vector<Residue>& residues) const;

private:
    size_t atom_count_ = 0;
    size_t residue_count_ = 0;
    std::unique_ptr<CovalentTopology> bonds_;

    // Plain const-style fields populated at construction. No Optional,
    // no Attach lifecycle, no Has gate. Empty / zero is the signal.
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
};


// Compose the per-atom AtomSemanticTable vector by walking each atom in
// each residue, computing typed mechanical identity from the canonical
// atom name + heavy-atom parent's canonical name (via the lifted
// parser in src/generated/LegacyAmberSemanticTables.h), and looking up
// the substrate row via LookupBy / LookupCap / ApplyCapDelta per
// spec/plan/topology-encoding-dependencies-2026-05-05.md §H.5.
//
// Stub-fixture guard: if no residue carries non-empty atom names
// across all of its atoms (calculator-physics tests with raw element-
// and-position fixtures), returns an empty vector. Calculators that
// need substrate gate on LegacyAmberTopology::HasAtomSemantic().
//
// Fails fatally on:
//   - chain atom whose (residue, variant, identity) has no LookupBy match.
//   - cap-only atom whose (terminal_state, identity) has no LookupCap match.
// All these cases indicate substrate gap or naming variance not caught
// by post-protonation re-canonicalisation; they are bugs to fix, not
// silent skips. Per the Bundle B brief: "Fail-loudly on substrate
// misses." (Predecessor stash@{0} chose degrade-and-skip; rejected.)
std::vector<AtomSemanticTable> ComposeAtomSemantic(
    const std::vector<std::unique_ptr<Atom>>& atoms,
    const std::vector<Residue>& residues,
    const CovalentTopology& bonds);

}  // namespace nmr
