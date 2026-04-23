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
// for BuildProtein. Charges come from the TPR (CHARMM36m convention
// for the fleet). The protein is not finalized here — FinalizeProtein
// must be called by the frame handler once first-frame positions are
// available (bond detection needs geometry).

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

ProteinConformation& TrajectoryProtein::CanonicalConformation() {
    return protein_->ConformationAt(0);
}

const ProteinConformation& TrajectoryProtein::CanonicalConformation() const {
    return protein_->ConformationAt(0);
}


// ── TickConformation ─────────────────────────────────────────────
//
// Ephemeral per-frame conformation. Caller owns the returned
// unique_ptr and releases it at the end of its iteration; the
// ProteinConformation's ConformationResults die with it.

std::unique_ptr<ProteinConformation>
TrajectoryProtein::TickConformation(std::vector<Vec3> positions) const {
    return std::make_unique<ProteinConformation>(
        protein_.get(), std::move(positions));
}


// ── InitTrajectoryAtoms ──────────────────────────────────────────
//
// One TrajectoryAtom per wrapped Protein atom, constructed once and
// never resized. TrajectoryAtom has a private constructor; friend
// access lets TrajectoryProtein construct them.

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
// Same shape as ProteinConformation::AttachResult: singleton-per-type
// check; TrajectoryResult-dependency check. ConformationResult-type
// dependencies (those referring to per-frame calculators that must
// run on each ProteinConformation) are validated at a different layer
// — by Trajectory::Run against the RunConfiguration's per-frame
// calculator set. This method cannot see the RunConfiguration; it
// only reports which declared dependencies are not met by already-
// attached TrajectoryResults, leaving the rest for Trajectory::Run.

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
// The trajectory protein's named per-frame operation: send this
// frame's state to every attached TrajectoryResult in attach order.
// Attach order is dispatch order so that a Result which reads
// another Result's per-atom output (via tp.AtomAt(i) fields) sees
// values already written earlier in the same frame's iteration.
//
// Thread-through of Trajectory is explicit: Results that need
// run-scope access (env, selections, EDR) receive it as a parameter
// during the only call where it's meaningful.

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
// TrajectoryProtein emits its own metadata (protein identity, atom
// count) and delegates the field-level output to each attached
// TrajectoryResult via WriteH5Group. No central column list.

void TrajectoryProtein::WriteH5(HighFive::File& file) const {
    if (!protein_) return;

    const size_t N = atoms_.size();

    file.createAttribute("protein_id", protein_id_);
    file.createAttribute("n_atoms", N);
    file.createAttribute("finalized", finalized_);

    // /atoms/identity/ — minimal per-atom identity passthrough from
    // the wrapped Protein. Per-atom typed identity (NmrAtomIdentity)
    // is deferred per WIP Appendix H — not yet computed, so we emit
    // only element + residue_index here. Downstream consumers that
    // want the richer identity read from /topology/ if present.
    std::vector<int> elements(N);
    std::vector<size_t> residue_indices(N);
    std::vector<std::string> atom_names(N);
    for (size_t i = 0; i < N; ++i) {
        const auto& a = protein_->AtomAt(i);
        elements[i] = static_cast<int>(a.element);
        residue_indices[i] = a.residue_index;
        atom_names[i] = a.pdb_atom_name;
    }
    auto atoms = file.createGroup("/atoms");
    atoms.createDataSet("element", elements);
    atoms.createDataSet("residue_index", residue_indices);
    atoms.createDataSet("pdb_atom_name", atom_names);

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
