#pragma once
//
// TrajectoryProtein: model of a protein in trajectory context.
//
// Replaces GromacsProtein (moved to learn/bones/). Role is unchanged
// in kind — the running buffer at trajectory scope — but cleaned up
// in content per spec/WIP_OBJECT_MODEL.md §3:
//
//   - Per-atom state is TrajectoryAtom (typed output fields only, no
//     Welford instances). Accumulator state lives inside the relevant
//     TrajectoryResult subclass, not on the per-atom struct.
//
//   - Per-frame work dispatches polymorphically through attached
//     TrajectoryResults (`DispatchCompute` iterates ResultsInAttachOrder
//     and calls Compute on each). No monolithic AccumulateFrame that
//     knows every subsystem; the orchestration is named and lives on
//     TrajectoryProtein where the attached Results live.
//
//   - Serialisation traverses attached TrajectoryResults, each one
//     writing its own group (WriteH5Group). No AllWelfords() central
//     enumeration.
//
// Identity delegation: topology (atoms, bonds, rings, residues, build
// context) lives on the wrapped Protein. TrajectoryProtein does not
// duplicate it.
//
// Ownership: TrajectoryProtein OWNS the Protein (via unique_ptr) and
// the ChargeSource (also unique_ptr). Its lifetime exceeds the
// wrapped Protein's and the attached TrajectoryResults'.
//

#include "Protein.h"
#include "ChargeSource.h"
#include "BondedEnergyResult.h"       // BondedParameters (topology-scope, set at build)
#include "FullSystemReader.h"
#include "TrajectoryAtom.h"
#include "TrajectoryResult.h"
#include "DenseBuffer.h"

#include <memory>
#include <string>
#include <typeindex>
#include <unordered_map>
#include <vector>

namespace HighFive { class File; }

namespace nmr {

class ProteinConformation;
class Trajectory;

class TrajectoryProtein {
public:
    TrajectoryProtein();
    ~TrajectoryProtein();

    TrajectoryProtein(const TrajectoryProtein&) = delete;
    TrajectoryProtein& operator=(const TrajectoryProtein&) = delete;

    // ── Build ────────────────────────────────────────────────────

    // Build from a full-system trajectory directory. dir_path must
    // contain md.tpr, md.xtc, md.edr. Reads TPR once for topology
    // (atom ranges, bonded parameters, charges). Builds the wrapped
    // Protein (no conformations yet — FinalizeProtein() completes
    // construction from first frame positions).
    //
    // Returns false on error; call Error() for message.
    bool BuildFromTrajectory(const std::string& dir_path);

    // Seed the trajectory with the canonical conformation (conf0).
    //
    // conf0 is the protein's static-analysis anchor within the
    // trajectory — the same ProteinConformation that static paths
    // work with via Protein::Conformation(). It lives permanently in
    // Protein.conformations_; its ConformationResults persist for the
    // whole run so any consumer that wants topology-invariant fields
    // (enrichment flags, partial_charge, graph distances, ring
    // memberships) can read them directly from the canonical
    // conformation via CanonicalConformation().
    //
    // Does three things:
    //   1. Protein::FinalizeConstruction(positions) — bond + ring
    //      detection from frame-0 geometry.
    //   2. protein.AddMDFrame(...) — creates conf0 as an
    //      MDFrameConformation and seats it in Protein.conformations_.
    //   3. InitTrajectoryAtoms — allocates one TrajectoryAtom per
    //      protein atom.
    //
    // Does NOT run per-frame calculators on conf0 — that's
    // Trajectory::Run's job after Seed returns.
    void Seed(std::vector<Vec3> first_frame_positions,
              double time_ps);

    // Canonical conformation (conf0). Same object Protein::Conformation()
    // returns in static paths. Accessor-only — permanence is provided
    // by Protein, not by TrajectoryProtein.
    ProteinConformation& CanonicalConformation();
    const ProteinConformation& CanonicalConformation() const;

    // Create an ephemeral per-frame ProteinConformation for frame i>0.
    // The returned conformation points at the wrapped Protein for
    // topology but is NOT added to Protein.conformations_ — its
    // lifetime is the caller's per-frame iteration. Per-frame
    // ConformationResults attach to it, populate its ConformationAtoms,
    // and die with it.
    std::unique_ptr<ProteinConformation> TickConformation(
            std::vector<Vec3> positions) const;

    // ── Identity delegation ──────────────────────────────────────

    Protein& ProteinRef() { return *protein_; }
    const Protein& ProteinRef() const { return *protein_; }

    size_t AtomCount() const { return protein_ ? protein_->AtomCount() : 0; }

    const std::string& ProteinId() const { return protein_id_; }
    ChargeSource* Charges() { return charges_.get(); }
    const ChargeSource* Charges() const { return charges_.get(); }
    int NetCharge() const { return net_charge_; }

    // FullSystemReader: topology + frame-splitting helper. Borrowed
    // by GromacsFrameHandler for ExtractFrame() on each frame.
    const FullSystemReader& SysReader() const { return sys_reader_; }

    // Bonded force-field parameters from TPR (topology-scope, set
    // once at BuildFromTrajectory). Lives here because it describes
    // protein topology, not a per-frame quantity.
    const BondedParameters& BondedParams() const { return bonded_params_; }
    bool HasBondedParams() const { return !bonded_params_.interactions.empty(); }

    // ── Per-atom trajectory data store ───────────────────────────

    const TrajectoryAtom& AtomAt(size_t i) const { return atoms_[i]; }
    TrajectoryAtom& MutableAtomAt(size_t i) { return atoms_[i]; }
    const std::vector<TrajectoryAtom>& Atoms() const { return atoms_; }

    // ── Attached TrajectoryResults (singleton-per-class) ─────────

    // Attach a TrajectoryResult. Checks singleton (one per type) and
    // dependencies (every declared type_index must already be attached
    // as another TrajectoryResult OR must be a ConformationResult type
    // — the latter is validated by Trajectory::Run against the
    // RunConfiguration's per-frame set).
    //
    // Returns false if rejected; reason is logged via OperationLog.
    bool AttachResult(std::unique_ptr<TrajectoryResult> result);

    template <typename T>
    T& Result();

    template <typename T>
    const T& Result() const;

    template <typename T>
    bool HasResult() const;

    const std::unordered_map<std::type_index,
                             std::unique_ptr<TrajectoryResult>>&
    AllResults() const { return results_; }

    // Attach-order iteration. Per-frame dispatch and serialisation
    // iterate in this order (respects declared dependencies by
    // construction: dependencies must be attached first).
    const std::vector<TrajectoryResult*>& ResultsInAttachOrder() const {
        return results_attach_order_;
    }

    // ── Per-frame dispatch ───────────────────────────────────────

    // The trajectory protein sends this frame's state to every
    // attached TrajectoryResult, in attach order. Called by
    // Trajectory::Run once per frame after the per-frame
    // ConformationResult pipeline has populated `conf`. Never
    // called by the frame reader — orchestration lives here and in
    // Trajectory::Run, not in the format-specific reader.
    void DispatchCompute(const ProteinConformation& conf,
                         Trajectory& traj,
                         std::size_t frame_idx,
                         double time_ps);

    // ── End-of-stream ────────────────────────────────────────────

    // Called by Trajectory::Run after the last frame. Iterates
    // attached results and calls Finalize on each, then sets
    // finalized_ = true.
    void FinalizeAllResults(Trajectory& traj);

    bool IsFinalized() const { return finalized_; }

    // ── Dense buffers ────────────────────────────────────────────

    // TrajectoryResults with per-atom × (frame|lag|frequency) output
    // transfer ownership of a DenseBuffer to TrajectoryProtein at
    // Finalize. The buffer is keyed by the owning TrajectoryResult's
    // type_index; query methods on that Result dereference into it.
    template <typename T>
    void AdoptDenseBuffer(std::unique_ptr<DenseBuffer<T>> buffer,
                          std::type_index owner) {
        dense_buffers_[owner] = std::move(buffer);
    }

    // Typed retrieval. Returns nullptr if no buffer is registered
    // under the given owner type or if the stored buffer's element
    // type differs.
    template <typename T>
    DenseBuffer<T>* GetDenseBuffer(std::type_index owner) {
        auto it = dense_buffers_.find(owner);
        if (it == dense_buffers_.end()) return nullptr;
        return dynamic_cast<DenseBuffer<T>*>(it->second.get());
    }

    // ── Serialisation ────────────────────────────────────────────

    // Each attached TrajectoryResult writes its own group/arrays.
    // TrajectoryProtein does not own a central schema; it traverses
    // the model. The /atoms/ group is emitted here as a single
    // per-atom collection of primary identifiers (residue index,
    // element) sourced from the wrapped Protein.
    void WriteH5(HighFive::File& file) const;

    // NPY per-result output into output_dir.
    int WriteFeatures(const std::string& output_dir) const;

    // ── Error ────────────────────────────────────────────────────

    const std::string& Error() const { return error_; }

private:
    void InitTrajectoryAtoms();

    std::unique_ptr<Protein> protein_;
    std::unique_ptr<ChargeSource> charges_;
    int net_charge_ = 0;
    std::string protein_id_;
    std::string error_;

    FullSystemReader sys_reader_;
    BondedParameters bonded_params_;

    std::vector<TrajectoryAtom> atoms_;

    std::unordered_map<std::type_index, std::unique_ptr<TrajectoryResult>>
        results_;
    std::vector<TrajectoryResult*> results_attach_order_;

    std::unordered_map<std::type_index, std::unique_ptr<DenseBufferBase>>
        dense_buffers_;

    bool finalized_ = false;
};


// ── Template method bodies (inline in header) ────────────────────

template <typename T>
T& TrajectoryProtein::Result() {
    auto it = results_.find(std::type_index(typeid(T)));
    if (it == results_.end()) {
        fprintf(stderr,
                "FATAL: TrajectoryProtein::Result<%s> not attached\n",
                typeid(T).name());
        std::abort();
    }
    return static_cast<T&>(*it->second);
}

template <typename T>
const T& TrajectoryProtein::Result() const {
    auto it = results_.find(std::type_index(typeid(T)));
    if (it == results_.end()) {
        fprintf(stderr,
                "FATAL: TrajectoryProtein::Result<%s> not attached\n",
                typeid(T).name());
        std::abort();
    }
    return static_cast<const T&>(*it->second);
}

template <typename T>
bool TrajectoryProtein::HasResult() const {
    return results_.find(std::type_index(typeid(T))) != results_.end();
}

}  // namespace nmr
