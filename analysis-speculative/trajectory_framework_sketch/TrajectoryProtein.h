// TrajectoryProtein.h
//
// The running buffer at trajectory scope. Per WIP_OBJECT_MODEL.md §3.
//
// Wraps a Protein (identity + topology) and adds:
//   - a vector<TrajectoryAtom>, parallel to protein.atoms_ in indexing;
//   - attached TrajectoryResults, typed singleton-per-class;
//   - adopted dense buffers transferred from TrajectoryResults at
//     their Finalize.
//
// This class replaces GromacsProtein's role — same location in the
// system, corrected content per §3's anti-pattern warnings.

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_PROTEIN_H
#define TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_PROTEIN_H

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeindex>
#include <unordered_map>
#include <vector>

#include "DenseBuffer.h"
#include "Stubs.h"
#include "TrajectoryAtom.h"
#include "TrajectoryResult.h"

class TrajectoryProtein {
public:
    explicit TrajectoryProtein(std::unique_ptr<::Protein> protein);

    // Identity delegation — TrajectoryProtein does not duplicate
    // topology. Ask the wrapped Protein.
    const ::Protein& Protein() const;
    std::size_t AtomCount() const;

    // Per-atom trajectory data store.
    const TrajectoryAtom& AtomAt(std::size_t atom_idx) const;
    TrajectoryAtom& MutableAtomAt(std::size_t atom_idx);
    const std::vector<TrajectoryAtom>& Atoms() const;

    // Attach a TrajectoryResult. Typed singleton-per-class. Returns
    // false on rejection (duplicate attach or failed dependency check);
    // returns true on success. Per §3 code block + §4 attach discipline.
    bool AttachResult(std::unique_ptr<TrajectoryResult> result);

    // Typed result access, mirroring ProteinConformation::Result<T>().
    // Throws if missing. This is the ONE template pattern, per
    // PATTERNS.md "Template metaprogramming beyond Result<T>()".
    template <typename T>
    T& Result() const;

    template <typename T>
    bool HasResult() const;

    const std::unordered_map<std::type_index,
          std::unique_ptr<TrajectoryResult>>& AllResults() const;

    // Iteration order matches attach order. Per §5 "Ordering guarantees"
    // this is load-bearing for per-frame dispatch + Finalize dispatch.
    const std::vector<TrajectoryResult*>& ResultsInAttachOrder() const;

    // Dense buffers transferred from TrajectoryResults at their
    // Finalize. Per §4 "Dense buffers on TrajectoryProtein". The
    // owner type_index keys the storage; query methods on the owning
    // TrajectoryResult dereference.
    template <typename T>
    void AdoptDenseBuffer(std::unique_ptr<DenseBuffer<T>> buffer,
                          std::type_index owner);

    // State for gated reads. Set to true at end of Trajectory::Run
    // Finalize phase.
    bool IsFinalized() const { return finalized_; }
    void SetFinalized(bool v) { finalized_ = v; }

private:
    std::unique_ptr<::Protein> protein_;
    std::vector<TrajectoryAtom> atoms_;
    std::unordered_map<std::type_index,
          std::unique_ptr<TrajectoryResult>> results_;
    std::vector<TrajectoryResult*> results_attach_order_;
    std::unordered_map<std::type_index,
          std::unique_ptr<DenseBufferBase>> dense_buffers_;
    bool finalized_ = false;
};

// Template definitions must live in the header for instantiation.

template <typename T>
T& TrajectoryProtein::Result() const {
    const auto it = results_.find(std::type_index(typeid(T)));
    if (it == results_.end()) {
        throw std::runtime_error(
            std::string("TrajectoryProtein::Result<T>(): ") +
            typeid(T).name() + " not attached");
    }
    return *static_cast<T*>(it->second.get());
}

template <typename T>
bool TrajectoryProtein::HasResult() const {
    return results_.find(std::type_index(typeid(T))) != results_.end();
}

template <typename T>
void TrajectoryProtein::AdoptDenseBuffer(std::unique_ptr<DenseBuffer<T>> buffer,
                                         std::type_index owner) {
    dense_buffers_[owner] = std::move(buffer);
}

#endif // TRAJECTORY_FRAMEWORK_SKETCH_TRAJECTORY_PROTEIN_H
