#pragma once
//
// Protein: sequence and science data. What the molecule IS, independent
// of geometry. Does NOT hold positions. Does NOT hold computed properties.
//
// Owns conformations. Factory methods create typed conformations.
// Owns residues, atoms, rings (topology), bonds (topology), build context.
//

#include "Types.h"
#include "Atom.h"
#include "Residue.h"
#include "Bond.h"
#include "Ring.h"

#include <optional>
#include "CovalentTopology.h"
#include "ProteinTopology.h"
#include "LegacyAmberTopology.h"
#include "ForceFieldChargeTable.h"
#include "ProteinBuildContext.h"
#include "ProteinConformation.h"
#include <vector>
#include <memory>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <cstdio>
#include <cstdlib>

namespace nmr {

class ChargeSource;
class LegacyAmberTopology;

class Protein {
public:
    Protein() = default;
    ~Protein() = default;

    // Proteins live on the heap via unique_ptr. They do not move.
    // Conformations hold a raw Protein* back-pointer that is valid
    // for the Protein's lifetime. No move means no dangling pointers.
    Protein(Protein&&) = delete;
    Protein& operator=(Protein&&) = delete;
    Protein(const Protein&) = delete;
    Protein& operator=(const Protein&) = delete;

    // ================================================================
    // Atom access
    // ================================================================

    size_t AtomCount() const { return atoms_.size(); }
    const Atom& AtomAt(size_t i) const { return *atoms_[i]; }
    Atom& MutableAtomAt(size_t i) { return *atoms_[i]; }
    const std::vector<std::unique_ptr<Atom>>& Atoms() const { return atoms_; }
    size_t AddAtom(std::unique_ptr<Atom> atom);

    // ================================================================
    // Residue access
    // ================================================================

    size_t ResidueCount() const { return residues_.size(); }
    const Residue& ResidueAt(size_t i) const { return residues_[i]; }
    Residue& MutableResidueAt(size_t i) { return residues_[i]; }
    size_t AddResidue(Residue residue);

    // ================================================================
    // Backbone connectivity (canonical, bond-graph-driven)
    //
    // Returns true iff a `BondOrder::Peptide` (alias
    // `BondCategory::PeptideCN`) bond connects res(a).C to res(b).N
    // in the cifpp-derived bond graph (i.e. a precedes b along the
    // chain). The peptide-bond classification is set by
    // CovalentTopology.cpp's loader-time labelling: a bond IS tagged
    // PeptideCN exactly when one atom is `res_a.C` and the other is
    // `res_b.N` of a different residue. The walk and the loader tag
    // are mutually-reinforcing — we ask the bond for its category
    // rather than re-deriving it from atom slots.
    //
    // Geometry-native substrate: bypasses every classification quirk
    // that chain_id / sequence_number / terminal_state / insertion_code
    // can produce (antibody insertion codes, engineered chimeras with
    // non-monotonic numbering, residue numbering gaps with intact
    // covalent bonds, cyclic peptides where chain_id wraps, loader
    // leaving terminal_state Unknown).
    //
    // If C(a) or N(b) is missing from the residue's backbone-cache
    // (incomplete structure) the answer is false. The query is
    // O(degree_of_C) ≈ small constant.
    //
    // ALL calc-side residue-adjacency walks (omega, phi, psi, Tripeptide
    // BB/Neighbor DFT lookups, Larsen H-bond donor frames) MUST go
    // through these methods. Chain_id / sequence_number / terminal_state
    // / insertion_code-based adjacency is a banned anti-pattern; see
    // PATTERNS.md and OBJECT_MODEL.md "Backbone connectivity discipline
    // (2026-05-19)".
    bool BackboneConnected(size_t residue_a_idx, size_t residue_b_idx) const;

    // Predecessor / successor on the backbone (PeptideCN bond).
    //
    // BackbonePredecessor(ri): returns the residue index whose C is
    //   bonded to res(ri).N via a `BondOrder::Peptide` bond. Returns
    //   nullopt at chain-N-term, loader gaps, or when N is missing
    //   from the backbone-cache.
    // BackboneSuccessor(ri): symmetric -- residue whose N is bonded to
    //   res(ri).C via a peptide bond.
    //
    // Same `Bond::IsPeptideBond()` filter as `BackboneConnected` above;
    // see that method's docstring for the loader-tag-invariant argument.
    //
    // These are the wrap-correct primitives for "walk to the next/prev
    // residue on the backbone." For cyclic peptides where C(last) is
    // bonded to N(0), Predecessor(0) returns last and Successor(last)
    // returns 0 -- the bond graph carries the answer directly. For
    // antibody insertion-coded structures (100 -> 100A -> 100B -> 101),
    // the walk follows the actual covalent topology, not the seq_number
    // labels. For ACE/NME caps, the cap residue's C/N participate in the
    // bond graph and are returned as the neighbour just as any internal
    // residue would be.
    //
    // **Test coverage status (2026-05-19):** Only linear single-chain
    // proteins with monotonic numbering and no insertion codes
    // (1P9J / 1Z9B / current OF3 fleet) are NUMERICALLY exercised by
    // the per-calc tests. Multi-chain, ACE/NME-capped, insertion-coded,
    // and cyclic structures are not exercised by any current fixture.
    // The methods walk the cifpp bond graph and are correct on those
    // shapes by-substrate-construction (the bond graph is authoritative),
    // but the numerical paths through them have not been validated
    // against real data of those classes. A multi-chain / cyclic /
    // insertion-coded fixture would close this.
    //
    // **Cyclic-peptide caveat:** `Residue.terminal_state` is loader-
    // assigned chain-order metadata, NOT a validity signal. For cyclic
    // peptides Predecessor(0) / Successor(N-1) return wrap edges and
    // dihedrals at "termini" are finite; the residue_terminal_state
    // field may still report NTerminus/CTerminus for the same rows.
    // Consumers must use isfinite(phi/psi/omega), not terminal_state,
    // for validity. The DihedralTimeSeries H5 group documents this in
    // its residue_terminal_state_legend attr.
    std::optional<size_t> BackbonePredecessor(size_t residue_idx) const;
    std::optional<size_t> BackboneSuccessor(size_t residue_idx) const;

    // ================================================================
    // Ring access (delegated through RingTopology on LegacyAmberTopology)
    //
    // Bundle C / Slice B (2026-05-07): ring storage moved to
    // LegacyAmberTopology::Rings(), parallel to how bonds delegate
    // through CovalentTopology. The public Protein API is unchanged
    // for aromatic rings — `RingCount()` / `RingAt(i)` / `Rings()`
    // still mean aromatic-only — and gains symmetric accessors for
    // saturated rings (Pro pyrrolidine). Calculator code consumers
    // see the same calling shape.
    //
    // Stub-fixture path: when LegacyAmberTopology has no ring topology
    // attached (no FinalizeConstruction yet), these accessors return
    // 0 / abort respectively, mirroring BondCount() / BondAt().
    // ================================================================

    size_t RingCount() const;
    const Ring& RingAt(size_t i) const;
    const std::vector<std::unique_ptr<Ring>>& Rings() const;

    size_t SaturatedRingCount() const;
    const Ring& SaturatedRingAt(size_t i) const;
    const std::vector<std::unique_ptr<Ring>>& SaturatedRings() const;

    // ================================================================
    // Bond access (delegated through CovalentTopology)
    // ================================================================

    size_t BondCount() const;
    const Bond& BondAt(size_t i) const;
    const std::vector<Bond>& Bonds() const;
    const CovalentTopology& BondTopology() const;

    // ================================================================
    // Explicit topology / loaded charge contracts
    // ================================================================

    bool HasTopology() const { return protein_topology_ != nullptr; }
    const ProteinTopology& TopologyBase() const;

    template<class TopologyT>
    const TopologyT& TopologyAs() const {
        static_assert(std::is_base_of_v<ProteinTopology, TopologyT>,
                      "TopologyAs<T> requires a ProteinTopology subtype");
        const auto* typed = dynamic_cast<const TopologyT*>(protein_topology_.get());
        if (!typed) {
            fprintf(stderr, "FATAL: requested topology %s is not attached.\n",
                    typeid(TopologyT).name());
            std::abort();
        }
        return *typed;
    }

    const LegacyAmberTopology& LegacyAmber() const;

    bool HasForceFieldCharges() const { return force_field_charges_ != nullptr; }
    const ForceFieldChargeTable& ForceFieldCharges() const;
    void SetForceFieldCharges(std::unique_ptr<ForceFieldChargeTable> charges);
    bool PrepareForceFieldCharges(const ChargeSource& source,
                                  const ProteinConformation& conf,
                                  std::string& error_out);

    // ================================================================
    // Build context
    // ================================================================

    const ProteinBuildContext& BuildContext() const { return *build_context_; }
    void SetBuildContext(std::unique_ptr<ProteinBuildContext> ctx) {
        build_context_ = std::move(ctx);
    }

    // ================================================================
    // Conformation factory methods
    // ================================================================

    // Base conformation — provenance unknown or not yet classified.
    // Use when the source doesn't warrant a typed subclass.
    ProteinConformation& AddConformation(
        std::vector<Vec3> positions,
        std::string description = "");

    CrystalConformation& AddCrystalConformation(
        std::vector<Vec3> positions,
        double resolution, double r_factor,
        double temperature, std::string pdb_id);

    PredictionConformation& AddPrediction(
        std::vector<Vec3> positions,
        std::string method,
        double confidence = std::nan(""));

    MDFrameConformation& AddMDFrame(
        std::vector<Vec3> positions,
        int walker, double time_ps, double weight,
        double rmsd_nm, double rg_nm);

    DerivedConformation& AddDerived(
        const ProteinConformation& parent,
        std::string description,
        std::vector<Vec3> positions);

    // Access crystal conformation (exactly one)
    CrystalConformation& CrystalConf();
    const CrystalConformation& CrystalConf() const;
    bool HasCrystalConformation() const { return crystal_index_ != SIZE_MAX; }

    // Access predictions
    size_t PredictionCount() const { return prediction_indices_.size(); }
    PredictionConformation& PredictionAt(size_t i);
    const PredictionConformation& PredictionAt(size_t i) const;

    // Access MD frames
    size_t MDFrameCount() const { return md_frame_indices_.size(); }
    MDFrameConformation& MDFrameAt(size_t i);
    const MDFrameConformation& MDFrameAt(size_t i) const;

    // Primary conformation: the one every loader creates. This is what
    // calculators and tests use. Returns ProteinConformation& regardless
    // of the subtype. Aborts if no conformations exist.
    ProteinConformation& Conformation();
    const ProteinConformation& Conformation() const;

    // All conformations (by index)
    size_t ConformationCount() const { return conformations_.size(); }
    ProteinConformation& ConformationAt(size_t i) { return *conformations_[i]; }
    const ProteinConformation& ConformationAt(size_t i) const { return *conformations_[i]; }

    // ================================================================
    // Construction helpers
    // ================================================================

    // FinalizeConstruction: must be called by every loader after adding
    // all atoms and residues. Caches backbone indices, detects bonds
    // (via OpenBabel), detects aromatic rings from residue types, and
    // constructs the LegacyAmberTopology — moving any source-provided
    // invariant FF data into its plain fields. Positions are needed
    // for bond detection (covalent radius check). Call BEFORE creating
    // any ProteinConformation.
    //
    // `invariants` is an optional value-pack of FF-numerical fields the
    // loader extracted from the source (TPR, PRMTOP). Loaders without
    // FF data pass `{}` (default); the topology's corresponding fields
    // remain empty and that is the legitimate signal. There is no
    // post-construction attach path; the value-pack is moved into the
    // topology and then goes out of scope.
    void FinalizeConstruction(const std::vector<Vec3>& positions,
                              LegacyAmberInvariants invariants = {},
                              double bond_tolerance = 0.4);

    // Individual steps (public for testing, prefer FinalizeConstruction).
    // CacheResidueBackboneIndices() runs the string-matched pass; the
    // typed pass at the end of FinalizeConstruction is private and
    // overwrites the cache from the substrate.
    //
    // Bundle C / Slice B (2026-05-07): DetectAromaticRings was deleted.
    // Ring construction is now substrate-driven via
    // RingTopology::ConstructFromSubstrate, called inside
    // FinalizeConstruction after ComposeAtomSemantic. Tests that
    // previously called DetectAromaticRings() directly must call
    // FinalizeConstruction() with positions instead.
    void CacheResidueBackboneIndices();

private:
    void ResolveResidueTerminalStates();
    // Two-pass discipline:
    //   bonds == nullptr: first pass (no covalent topology yet). Detects
    //     HIS / ASP / GLU / LYS / TYR / ARG variants from explicit
    //     hydrogen presence. CYS variant requires bonds; deferred.
    //   bonds != nullptr: second pass. Adds CYS -> CYX detection from
    //     the resolved disulfide bond list.
    // The second pass runs BEFORE the LegacyAmberTopology that carries
    // atom_semantic is constructed, so substrate composition has the
    // final variant_index for every residue.
    void ResolveProtonationStates(const CovalentTopology* bonds);

    // Typed CacheResidueBackboneIndices: reads BackboneRole + Locant
    // from `LegacyAmber().AtomSemantic()` and overwrites res.{N, CA,
    // C, O, H, HA, CB} with substrate-driven indices. Glycine HA
    // resolves via Locant::Alpha + DiastereotopicIndex::Position2;
    // CB via Locant::Beta. Pro res.H stays NONE because PRO chain
    // table drops backbone H per substrate dependencies §H.10. Chi-
    // angle resolver inside this function is intentionally NOT typed
    // — kept string-matched per Audit Hotspot 2; separate slice.
    void CacheResidueBackboneIndices_Typed();

    std::vector<std::unique_ptr<Atom>> atoms_;
    std::vector<Residue> residues_;
    // Rings live on LegacyAmberTopology::Rings() post-Bundle-C / Slice B.
    // Protein keeps no ring storage of its own; ring access goes through
    // delegation methods declared above.
    std::unique_ptr<ProteinTopology> protein_topology_;
    std::unique_ptr<ForceFieldChargeTable> force_field_charges_;
    std::unique_ptr<ProteinBuildContext> build_context_ =
        std::make_unique<ProteinBuildContext>();
    std::vector<std::unique_ptr<ProteinConformation>> conformations_;
    size_t crystal_index_ = SIZE_MAX;
    std::vector<size_t> prediction_indices_;
    std::vector<size_t> md_frame_indices_;
};

}  // namespace nmr
