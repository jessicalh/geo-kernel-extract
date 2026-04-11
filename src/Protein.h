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
#include "CovalentTopology.h"
#include "ProteinBuildContext.h"
#include "ProteinConformation.h"
#include <vector>
#include <memory>
#include <string>

namespace nmr {

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
    // Ring access
    // ================================================================

    size_t RingCount() const { return rings_.size(); }
    const Ring& RingAt(size_t i) const { return *rings_[i]; }
    const std::vector<std::unique_ptr<Ring>>& Rings() const { return rings_; }

    // ================================================================
    // Bond access (delegated through CovalentTopology)
    // ================================================================

    size_t BondCount() const { return topology_ ? topology_->BondCount() : 0; }
    const Bond& BondAt(size_t i) const { return topology_->BondAt(i); }
    const std::vector<Bond>& Bonds() const { return topology_->Bonds(); }
    const CovalentTopology& Topology() const { return *topology_; }

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
    // (via OpenBabel), and detects aromatic rings from residue types.
    // Positions are needed for bond detection (covalent radius check).
    // Call BEFORE creating any ProteinConformation.
    void FinalizeConstruction(const std::vector<Vec3>& positions,
                              double bond_tolerance = 0.4);

    // Individual steps (public for testing, prefer FinalizeConstruction)
    void DetectAromaticRings();
    void CacheResidueBackboneIndices();

private:
    std::vector<std::unique_ptr<Atom>> atoms_;
    std::vector<Residue> residues_;
    std::vector<std::unique_ptr<Ring>> rings_;
    std::unique_ptr<CovalentTopology> topology_;
    std::unique_ptr<ProteinBuildContext> build_context_ =
        std::make_unique<ProteinBuildContext>();
    std::vector<std::unique_ptr<ProteinConformation>> conformations_;
    size_t crystal_index_ = SIZE_MAX;
    std::vector<size_t> prediction_indices_;
    std::vector<size_t> md_frame_indices_;
};

}  // namespace nmr
