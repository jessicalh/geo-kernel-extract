#pragma once
//
// TrajectoryProtein: the protein model in trajectory context.
// Canonical description in OBJECT_MODEL.md (trajectory-scope entities).
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

    // Parse md.tpr from dir_path; build Protein + charges + bonded
    // params. Does not seat conformation 0 (Seed does that on first
    // frame). Returns false on error; call Error() for the message.
    bool BuildFromTrajectory(const std::string& dir_path);

    // Finalize the wrapped Protein from first-frame geometry (bond +
    // ring detection), seat conformation 0 as an MDFrameConformation,
    // and allocate the TrajectoryAtom vector. Per-frame calculators
    // are NOT run here — that is Trajectory::Run's job.
    void Seed(std::vector<Vec3> first_frame_positions,
              double time_ps);

    // Conformation 0, permanent on the wrapped Protein after Seed.
    ProteinConformation& CanonicalConformation();
    const ProteinConformation& CanonicalConformation() const;

    // Ephemeral per-frame conformation pointing at the wrapped
    // Protein for topology. Caller owns; lifetime is the iteration.
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

    // Singleton-per-type check only. Returns false if another TR of
    // the same type is already attached. Dependency validation is in
    // Trajectory::Run Phase 4, not here.
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

    // Attach order is dispatch order.
    const std::vector<TrajectoryResult*>& ResultsInAttachOrder() const {
        return results_attach_order_;
    }

    // ── Per-frame dispatch ───────────────────────────────────────

    // Called once per frame by Trajectory::Run after the per-frame
    // ConformationResult pipeline populates `conf`. Iterates TRs in
    // attach order.
    void DispatchCompute(const ProteinConformation& conf,
                         Trajectory& traj,
                         std::size_t frame_idx,
                         double time_ps);

    // ── End-of-stream ────────────────────────────────────────────

    // Called by Trajectory::Run Phase 8 after the last frame.
    void FinalizeAllResults(Trajectory& traj);

    bool IsFinalized() const { return finalized_; }

    // ── Dense buffers ────────────────────────────────────────────

    // Keyed by the owning TR's type_index.
    template <typename T>
    void AdoptDenseBuffer(std::unique_ptr<DenseBuffer<T>> buffer,
                          std::type_index owner) {
        dense_buffers_[owner] = std::move(buffer);
    }

    // Returns nullptr if no buffer under the given owner type or if
    // the stored element type differs.
    template <typename T>
    DenseBuffer<T>* GetDenseBuffer(std::type_index owner) {
        auto it = dense_buffers_.find(owner);
        if (it == dense_buffers_.end()) return nullptr;
        return dynamic_cast<DenseBuffer<T>*>(it->second.get());
    }

    // ── Serialisation ────────────────────────────────────────────

    // Emits /atoms/ then delegates to each attached TR's WriteH5Group.
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
