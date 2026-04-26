#include "TrajectoryProtein.h"
#include "BuildResult.h"
#include "OperationLog.h"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <filesystem>

namespace fs = std::filesystem;

namespace nmr {


TrajectoryProtein::TrajectoryProtein() = default;
TrajectoryProtein::~TrajectoryProtein() = default;


// ── BuildFromTrajectory ──────────────────────────────────────────
//
// Single TPR parse: atom ranges, bonded parameters, stored topology
// for BuildProtein. Charges come from the TPR. The Protein is not
// finalized here — Seed does that on the first frame (bond + ring
// detection needs geometry).

bool TrajectoryProtein::BuildFromTrajectory(const std::string& dir_path) {
    protein_id_ = fs::path(dir_path).filename().string();

    const std::string tpr_path = dir_path + "/md.tpr";

    if (!sys_reader_.ReadTopology(tpr_path)) {
        error_ = "TPR topology: " + sys_reader_.error();
        return false;
    }

    // Bonded parameters from the same TPR parse (no re-read).
    bonded_params_ = sys_reader_.BondedParams();

    // Build Protein + ChargeSource from stored TPR data.
    auto build = sys_reader_.BuildProtein(protein_id_);
    if (!build.Ok()) {
        error_ = "TPR protein: " + build.error;
        return false;
    }

    protein_ = std::move(build.protein);
    charges_ = std::move(build.charges);
    net_charge_ = build.net_charge;

    OperationLog::Info(LogCalcOther, "TrajectoryProtein::BuildFromTrajectory",
        protein_id_ + ": " +
        std::to_string(protein_->AtomCount()) + " protein atoms, " +
        std::to_string(sys_reader_.Topology().water_count) + " water, " +
        std::to_string(sys_reader_.Topology().ion_count) + " ions" +
        (HasBondedParams() ?
            ", " + std::to_string(bonded_params_.interactions.size()) +
            " bonded interactions" : ""));

    return true;
}


// ── Seed ─────────────────────────────────────────────────────────
//
// Creates the canonical conformation (conf0) from first-frame
// geometry. See TrajectoryProtein.h for full contract.

void TrajectoryProtein::Seed(
        std::vector<Vec3> first_frame_positions, double time_ps) {
    // Bond detection needs first-frame geometry.
    protein_->FinalizeConstruction(first_frame_positions);

    // Canonical conformation: permanent, lives in Protein.conformations_.
    // Tick conformations created by TickConformation are free-standing
    // and not added to this vector.
    protein_->AddMDFrame(std::move(first_frame_positions),
                         /*walker=*/0, time_ps, /*weight=*/1.0,
                         /*rmsd_nm=*/0.0, /*rg_nm=*/0.0);

    InitTrajectoryAtoms();

    OperationLog::Info(LogCalcOther, "TrajectoryProtein::Seed",
        protein_id_ + ": " +
        std::to_string(protein_->AtomCount()) + " atoms, " +
        std::to_string(protein_->BondCount()) + " bonds, " +
        std::to_string(protein_->RingCount()) + " rings, " +
        "canonical conformation at t=" + std::to_string(time_ps) + "ps");
}


// ── CanonicalConformation ────────────────────────────────────────

const ProteinConformation& TrajectoryProtein::CanonicalConformation() const {
    return protein_->ConformationAt(0);
}

// Non-const path, private + friend-only. See TrajectoryProtein.h.
ProteinConformation& TrajectoryProtein::MutableCanonicalConformation_() {
    return protein_->ConformationAt(0);
}


// ── TickConformation ─────────────────────────────────────────────

std::unique_ptr<ProteinConformation>
TrajectoryProtein::TickConformation(std::vector<Vec3> positions) const {
    return std::make_unique<ProteinConformation>(
        protein_.get(), std::move(positions));
}


// ── InitTrajectoryAtoms ──────────────────────────────────────────
//
// TrajectoryAtom has a private constructor; friend access lets
// TrajectoryProtein construct them.

void TrajectoryProtein::InitTrajectoryAtoms() {
    atoms_.clear();
    atoms_.reserve(protein_->AtomCount());
    for (size_t i = 0; i < protein_->AtomCount(); ++i) {
        (void)i;
        atoms_.emplace_back(TrajectoryAtom{});
    }
}


// ── AttachResult ─────────────────────────────────────────────────
//
// Singleton-per-type check only. Dependency validation is in
// Trajectory::Run Phase 4 (this method cannot see RunConfiguration).

bool TrajectoryProtein::AttachResult(
        std::unique_ptr<TrajectoryResult> result) {
    if (!result) return false;

    const std::string name = result->Name();
    const std::type_index tid(typeid(*result));

    // Singleton check.
    if (results_.find(tid) != results_.end()) {
        OperationLog::Warn("TrajectoryProtein::AttachResult",
            "rejected " + name + ": already attached");
        return false;
    }

    results_[tid] = std::move(result);
    results_attach_order_.push_back(results_[tid].get());

    OperationLog::Info(LogResultAttach,
        "TrajectoryProtein::AttachResult",
        "attached " + name);
    return true;
}


// ── DispatchCompute ──────────────────────────────────────────────
//
// Attach order is dispatch order: a TR reading another TR's per-atom
// fields sees values already written earlier in the same frame.

void TrajectoryProtein::DispatchCompute(
        const ProteinConformation& conf,
        Trajectory& traj,
        std::size_t frame_idx, double time_ps) {
    for (TrajectoryResult* r : results_attach_order_) {
        r->Compute(conf, *this, traj, frame_idx, time_ps);
    }
}


// ── FinalizeAllResults ───────────────────────────────────────────

void TrajectoryProtein::FinalizeAllResults(Trajectory& traj) {
    for (TrajectoryResult* r : results_attach_order_) {
        r->Finalize(*this, traj);
    }
    finalized_ = true;
}


// ── WriteH5 ──────────────────────────────────────────────────────
//
// File-root attributes + /atoms/ passthrough; each attached TR
// writes its own group.

void TrajectoryProtein::WriteH5(HighFive::File& file) const {
    if (!protein_) return;

    const size_t N = atoms_.size();

    file.createAttribute("protein_id", protein_id_);
    file.createAttribute("n_atoms", N);
    file.createAttribute("finalized", finalized_);

    // /atoms/ — minimal per-atom identity passthrough: element,
    // residue_index, iupac_name. Richer typed identity
    // (NmrAtomIdentity) is deferred.
    std::vector<int> elements(N);
    std::vector<size_t> residue_indices(N);
    std::vector<std::string> atom_names(N);
    for (size_t i = 0; i < N; ++i) {
        const auto& a = protein_->AtomAt(i);
        elements[i] = static_cast<int>(a.element);
        residue_indices[i] = a.residue_index;
        atom_names[i] = a.iupac_name.AsString();
    }
    auto atoms = file.createGroup("/atoms");
    atoms.createDataSet("element", elements);
    atoms.createDataSet("residue_index", residue_indices);
    atoms.createDataSet("iupac_name", atom_names);

    // Each attached TrajectoryResult writes its own group.
    for (const TrajectoryResult* r : results_attach_order_) {
        r->WriteH5Group(*this, file);
    }
}


// ── WriteFeatures ────────────────────────────────────────────────

int TrajectoryProtein::WriteFeatures(const std::string& output_dir) const {
    int total = 0;
    for (const TrajectoryResult* r : results_attach_order_) {
        total += r->WriteFeatures(*this, output_dir);
    }
    return total;
}

}  // namespace nmr
