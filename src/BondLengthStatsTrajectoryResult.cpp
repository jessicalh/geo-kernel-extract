#include "BondLengthStatsTrajectoryResult.h"
#include "TrajectoryProtein.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Bond.h"
#include "TrajectoryMoments.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>

namespace nmr {


std::unique_ptr<BondLengthStatsTrajectoryResult>
BondLengthStatsTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<BondLengthStatsTrajectoryResult>();
    r->per_bond_.assign(tp.ProteinRef().BondCount(), PerBondWelford{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Per frame, per bond: compute current length, update the Welford on
// per_bond_[i] + the frame-to-frame delta tracker. Same pattern as
// BsWelford, just bond-indexed instead of atom-indexed.
//
// Protein's Bond objects carry atom_index_a/b as topology — the same
// pair across all frames. Length varies with frame geometry.

void BondLengthStatsTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const Protein& protein = tp.ProteinRef();
    const std::size_t B = protein.BondCount();

    for (std::size_t bi = 0; bi < B; ++bi) {
        const Bond& bond = protein.BondAt(bi);
        const Vec3 a = conf.PositionAt(bond.atom_index_a);
        const Vec3 b = conf.PositionAt(bond.atom_index_b);
        const double length = (b - a).norm();

        PerBondWelford& pb = per_bond_[bi];
        const std::size_t n_new = pb.n_frames + 1;

        MomentsUpdate(pb.length_mean, pb.length_m2, length, n_new);
        MinMaxUpdate(pb.length_min, pb.length_max, length);

        if (pb.has_prev) {
            const double delta = length - pb.prev_length;
            const std::size_t dn_new = pb.delta_n + 1;
            MomentsUpdate(pb.delta_mean, pb.delta_m2, delta, dn_new);
            pb.delta_n = dn_new;
        }
        pb.prev_length = length;
        pb.has_prev = true;
        pb.n_frames = n_new;
    }

    ++n_frames_;
    (void)frame_idx;
}


// ── Finalize ─────────────────────────────────────────────────────

void BondLengthStatsTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                               Trajectory& traj) {
    (void)traj;
    for (auto& pb : per_bond_) {
        pb.length_std = MomentsStd(pb.length_m2, pb.n_frames);
        pb.delta_std  = MomentsStd(pb.delta_m2,  pb.delta_n);
    }
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "BondLengthStatsTrajectoryResult::Finalize",
        "finalized " + std::to_string(per_bond_.size()) +
        " bonds across " + std::to_string(n_frames_) + " frames");
    (void)tp;
}


// ── WriteH5Group ─────────────────────────────────────────────────

void BondLengthStatsTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const Protein& protein = tp.ProteinRef();
    const std::size_t B = per_bond_.size();

    std::vector<double> mean(B), std_(B), minv(B), maxv(B);
    std::vector<double> dmean(B), dstd(B);
    std::vector<std::uint64_t> aa(B), ab(B);
    std::vector<std::int8_t>   ord(B), cat(B);

    for (std::size_t i = 0; i < B; ++i) {
        const auto& pb = per_bond_[i];
        const Bond& bond = protein.BondAt(i);
        mean[i]  = pb.length_mean;
        std_[i]  = pb.length_std;
        minv[i]  = pb.length_min;
        maxv[i]  = pb.length_max;
        dmean[i] = pb.delta_mean;
        dstd[i]  = pb.delta_std;
        aa[i]    = static_cast<std::uint64_t>(bond.atom_index_a);
        ab[i]    = static_cast<std::uint64_t>(bond.atom_index_b);
        ord[i]   = static_cast<std::int8_t>(bond.order);
        cat[i]   = static_cast<std::int8_t>(bond.category);
    }

    auto grp = file.createGroup("/trajectory/bond_length_stats");
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_bonds",  B);
    grp.createAttribute("n_frames", n_frames_);
    grp.createAttribute("finalized", finalized_);
    grp.createAttribute("units",    std::string("Angstrom"));

    grp.createDataSet("length_mean",       mean);
    grp.createDataSet("length_std",        std_);
    grp.createDataSet("length_min",        minv);
    grp.createDataSet("length_max",        maxv);
    grp.createDataSet("length_delta_mean", dmean);
    grp.createDataSet("length_delta_std",  dstd);
    grp.createDataSet("atom_a",   aa);
    grp.createDataSet("atom_b",   ab);
    grp.createDataSet("order",    ord);
    grp.createDataSet("category", cat);
}

}  // namespace nmr
