// TrajectoryProtein.cpp
//
// Per WIP_OBJECT_MODEL.md §3 and §4's attach discipline.

#include "TrajectoryProtein.h"

#include <stdexcept>
#include <utility>

TrajectoryProtein::TrajectoryProtein(std::unique_ptr<::Protein> protein)
    : protein_(std::move(protein)) {
    // Build atoms_ once, never resize. Per WIP §3: "the vector is built
    // once in the ProteinConformation constructor from the Protein's
    // atom count...and is never resized" (same discipline at trajectory
    // scope). Each TrajectoryAtom has a private ctor; this class is its
    // friend, so emplace_back works.
    const std::size_t n = protein_->AtomCount();
    atoms_.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        atoms_.emplace_back(TrajectoryAtom{});
    }
}

const ::Protein& TrajectoryProtein::Protein() const {
    return *protein_;
}

std::size_t TrajectoryProtein::AtomCount() const {
    return atoms_.size();
}

const TrajectoryAtom& TrajectoryProtein::AtomAt(std::size_t atom_idx) const {
    return atoms_[atom_idx];
}

TrajectoryAtom& TrajectoryProtein::MutableAtomAt(std::size_t atom_idx) {
    return atoms_[atom_idx];
}

const std::vector<TrajectoryAtom>& TrajectoryProtein::Atoms() const {
    return atoms_;
}

bool TrajectoryProtein::AttachResult(std::unique_ptr<TrajectoryResult> result) {
    if (!result) {
        return false;
    }
    const std::type_index tid(typeid(*result));

    // Singleton check. Per WIP §4 code block: "rejected: already attached".
    if (results_.find(tid) != results_.end()) {
        return false;
    }

    // Dependency validation happens in Trajectory::Run Phase 2, not here.
    // WIP §12 item 4 RESOLVED says Phase 2 is the one place with both
    // attached results and run config in scope; attach does not need to
    // know about run configuration.

    results_[tid] = std::move(result);
    results_attach_order_.push_back(results_[tid].get());
    return true;
}

const std::unordered_map<std::type_index,
      std::unique_ptr<TrajectoryResult>>&
TrajectoryProtein::AllResults() const {
    return results_;
}

const std::vector<TrajectoryResult*>&
TrajectoryProtein::ResultsInAttachOrder() const {
    return results_attach_order_;
}
