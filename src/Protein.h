#pragma once
//
// Protein: chemical identity and topology — what the molecule IS.
// Holds no coordinates and no computed properties; those live on
// ProteinConformation and ConformationAtom.
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

    // Non-copyable and non-movable: conformations hold a raw Protein*
    // back-pointer valid for the Protein's lifetime, so no move means no
    // dangling pointers.
    Protein(Protein&&) = delete;
    Protein& operator=(Protein&&) = delete;
    Protein(const Protein&) = delete;
    Protein& operator=(const Protein&) = delete;

    size_t AtomCount() const { return atoms_.size(); }
    const Atom& AtomAt(size_t i) const { return *atoms_[i]; }
    Atom& MutableAtomAt(size_t i) { return *atoms_[i]; }
    const std::vector<std::unique_ptr<Atom>>& Atoms() const { return atoms_; }
    size_t AddAtom(std::unique_ptr<Atom> atom);

    size_t ResidueCount() const { return residues_.size(); }
    const Residue& ResidueAt(size_t i) const { return residues_[i]; }
    Residue& MutableResidueAt(size_t i) { return residues_[i]; }
    size_t AddResidue(Residue residue);

    // Returns true iff a BondCategory::PeptideCN bond connects res(a).C
    // to res(b).N — i.e. a precedes b on the chain. False if C(a) or
    // N(b) is absent. Uses the bond graph, not residue numbering / chain
    // labels. Requires attached topology (aborts otherwise).
    bool BackboneConnected(size_t residue_a_idx, size_t residue_b_idx) const;

    // Predecessor / successor on the backbone (PeptideCN bond).
    //
    // BackbonePredecessor(ri): residue whose C is bonded to res(ri).N;
    //   nullopt at the N-terminus, a loader gap, or when N is absent.
    // BackboneSuccessor(ri): symmetric, off res(ri).C.
    //
    // Bond-graph-driven and wrap-correct on cyclic / multi-chain /
    // insertion-coded structures by construction; they do not consult
    // Residue.terminal_state.
    std::optional<size_t> BackbonePredecessor(size_t residue_idx) const;
    std::optional<size_t> BackboneSuccessor(size_t residue_idx) const;

    // Ring access — delegated to LegacyAmberTopology::Rings().
    // RingCount() / RingAt() / Rings() are aromatic-only; the
    // SaturatedRing* accessors cover saturated rings (Pro pyrrolidine).
    // Before FinalizeConstruction (no ring topology attached) the count
    // methods return 0 and the element/list accessors abort.

    size_t RingCount() const;
    const Ring& RingAt(size_t i) const;
    const std::vector<std::unique_ptr<Ring>>& Rings() const;

    size_t SaturatedRingCount() const;
    const Ring& SaturatedRingAt(size_t i) const;
    const std::vector<std::unique_ptr<Ring>>& SaturatedRings() const;

    size_t BondCount() const;
    const Bond& BondAt(size_t i) const;
    const std::vector<Bond>& Bonds() const;
    const CovalentTopology& BondTopology() const;

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

    const ProteinBuildContext& BuildContext() const { return *build_context_; }
    void SetBuildContext(std::unique_ptr<ProteinBuildContext> ctx) {
        build_context_ = std::move(ctx);
    }

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

    // The tracked crystal conformation (aborts if none).
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

    // The primary conformation, as ProteinConformation& regardless of
    // subtype. Aborts if there are none.
    ProteinConformation& Conformation();
    const ProteinConformation& Conformation() const;

    // All conformations (by index)
    size_t ConformationCount() const { return conformations_.size(); }
    ProteinConformation& ConformationAt(size_t i) { return *conformations_[i]; }
    const ProteinConformation& ConformationAt(size_t i) const { return *conformations_[i]; }

    // FinalizeConstruction: every loader must call this after adding all
    // atoms and residues, and BEFORE creating any ProteinConformation.
    // Caches backbone indices, resolves the covalent topology and rings,
    // and constructs the LegacyAmberTopology — moving any source-provided
    // FF invariants into its plain fields. Positions are needed for bond
    // detection (covalent-radius check).
    //
    // `invariants` is an optional value-pack of FF-numerical fields the
    // loader extracted from the source (TPR, PRMTOP). Loaders without FF
    // data pass `{}`; the topology's corresponding fields stay empty, and
    // that emptiness is the legitimate signal. The pack is moved into the
    // topology and then goes out of scope — there is no post-construction
    // attach path.
    void FinalizeConstruction(const std::vector<Vec3>& positions,
                              LegacyAmberInvariants invariants = {},
                              double bond_tolerance = 0.4);

    // CacheResidueBackboneIndices() (public for testing) runs the
    // string-matched backbone-index pass.
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
    // table drops backbone H per substrate dependencies §H.10.
    void CacheResidueBackboneIndices_Typed();

    std::vector<std::unique_ptr<Atom>> atoms_;
    std::vector<Residue> residues_;
    // Protein keeps no ring storage of its own; rings live on
    // LegacyAmberTopology::Rings() (delegation methods above).
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
