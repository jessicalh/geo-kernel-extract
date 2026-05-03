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
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace nmr {

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
                        LegacyAmberInvariants invariants = {});

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
};

}  // namespace nmr
