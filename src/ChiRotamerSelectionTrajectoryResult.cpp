#include "ChiRotamerSelectionTrajectoryResult.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "AminoAcidType.h"
#include "SelectionRecord.h"
#include "OperationLog.h"

#include <cmath>
#include <typeinfo>

namespace nmr {


namespace {

// Dihedral from four positions. Returns angle in radians in [-pi, pi].
// Robust to degenerate cases (returns 0 if any cross product is
// near-zero).
double ChiRadians(const Vec3& p0, const Vec3& p1,
                  const Vec3& p2, const Vec3& p3) {
    Vec3 b1 = p1 - p0;
    Vec3 b2 = p2 - p1;
    Vec3 b3 = p3 - p2;
    Vec3 n1 = b1.cross(b2);
    Vec3 n2 = b2.cross(b3);
    const double n1n = n1.norm();
    const double n2n = n2.norm();
    if (n1n < 1e-10 || n2n < 1e-10) return 0.0;
    n1 /= n1n;
    n2 /= n2n;
    const double b2n = b2.norm();
    if (b2n < 1e-10) return 0.0;
    const Vec3 m1 = n1.cross(b2 / b2n);
    const double x = n1.dot(n2);
    const double y = m1.dot(n2);
    return std::atan2(y, x);
}

// Three-bin classification of chi in radians. Bins are 120° segments
// centered on +60° (gauche+), 180° (trans), and -60° (gauche-).
ChiRotamerSelectionTrajectoryResult::RotamerBin
BinForChi(double chi_rad) {
    using Bin = ChiRotamerSelectionTrajectoryResult::RotamerBin;
    const double third = 2.0 * M_PI / 3.0;  // 120°
    if (chi_rad >  third) return Bin::Trans;
    if (chi_rad < -third) return Bin::Trans;
    if (chi_rad >  0.0)   return Bin::Gplus;
    return Bin::Gminus;
}

const char* BinName(ChiRotamerSelectionTrajectoryResult::RotamerBin b) {
    using Bin = ChiRotamerSelectionTrajectoryResult::RotamerBin;
    switch (b) {
        case Bin::Gplus:  return "gplus";
        case Bin::Trans:  return "trans";
        case Bin::Gminus: return "gminus";
    }
    return "unknown";
}

}  // namespace


std::unique_ptr<ChiRotamerSelectionTrajectoryResult>
ChiRotamerSelectionTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<ChiRotamerSelectionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    r->prev_bin_.assign(R, {RotamerBin::Gplus, RotamerBin::Gplus,
                            RotamerBin::Gplus, RotamerBin::Gplus});
    r->bin_valid_.assign(R, {false, false, false, false});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Per residue, per chi index k=0..3:
//   1. Skip if this residue doesn't have chi[k] (Residue.chi[k].Valid()
//      reports false when the AminoAcidType has fewer chi angles).
//   2. Compute chi_k from four position lookups through the cached
//      atom indices on Residue.
//   3. Classify into rotamer bin.
//   4. If bin_valid_[r][k] && new_bin != prev_bin_[r][k], push a
//      SelectionRecord onto traj.MutableSelections() describing the
//      transition.
//   5. Update prev_bin_ / bin_valid_ for next frame.

void ChiRotamerSelectionTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp;
    const Protein& protein = conf.ProteinRef();
    const std::size_t R = protein.ResidueCount();

    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);

        for (int k = 0; k < 4; ++k) {
            if (!res.chi[k].Valid()) break;

            const Vec3 p0 = conf.PositionAt(res.chi[k].a[0]);
            const Vec3 p1 = conf.PositionAt(res.chi[k].a[1]);
            const Vec3 p2 = conf.PositionAt(res.chi[k].a[2]);
            const Vec3 p3 = conf.PositionAt(res.chi[k].a[3]);
            const double chi_rad = ChiRadians(p0, p1, p2, p3);
            const RotamerBin new_bin = BinForChi(chi_rad);

            if (bin_valid_[ri][k] && new_bin != prev_bin_[ri][k]) {
                const char* res_name =
                    GetAminoAcidType(res.type).three_letter_code;
                const char* prev_name = BinName(prev_bin_[ri][k]);
                const char* new_name  = BinName(new_bin);

                std::string reason =
                    std::string("chi") + std::to_string(k + 1) +
                    "_rotamer_transition_" + res_name + "_" +
                    std::to_string(res.sequence_number) + "_" +
                    prev_name + "_to_" + new_name;

                traj.MutableSelections().Push(SelectionRecord(
                    std::type_index(typeid(
                        ChiRotamerSelectionTrajectoryResult)),
                    frame_idx, time_ps,
                    std::move(reason),
                    {
                        {"residue_index",
                            std::to_string(ri)},
                        {"residue_sequence",
                            std::to_string(res.sequence_number)},
                        {"chi_index",
                            std::to_string(k)},
                        {"bin_before", prev_name},
                        {"bin_after",  new_name},
                    }));
                ++n_transitions_;
            }

            prev_bin_[ri][k] = new_bin;
            bin_valid_[ri][k] = true;
        }
    }
}

}  // namespace nmr
