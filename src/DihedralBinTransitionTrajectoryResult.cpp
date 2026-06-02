#include "DihedralBinTransitionTrajectoryResult.h"

#include "AminoAcidType.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <typeinfo>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Dihedral from four positions. Same atan2-based formulation as
// PlanarGeometryResult, ChiRotamerSelectionTrajectoryResult, and
// DihedralTimeSeriesTrajectoryResult. Returns signed angle in radians in
// [-π, π]; NaN at degenerate inputs
// (collinear b2, zero-norm n1 or n2) so consumers can distinguish
// "indeterminate" from a real 0 rad measurement.
//
//   D(p1,p2,p3,p4) = atan2( (n1×b̂2)·n2, n1·n2 )
//     b1 = p2−p1, b2 = p3−p2, b3 = p4−p3
//     n1 = b1×b2, n2 = b2×b3
double Dihedral(const Vec3& p1, const Vec3& p2,
                const Vec3& p3, const Vec3& p4) {
    const Vec3 b1 = p2 - p1;
    const Vec3 b2 = p3 - p2;
    const Vec3 b3 = p4 - p3;
    const double b2n = b2.norm();
    if (b2n < 1e-10) return kNaN;
    const Vec3 n1 = b1.cross(b2);
    const Vec3 n2 = b2.cross(b3);
    if (n1.norm() < 1e-10 || n2.norm() < 1e-10) return kNaN;
    const Vec3 m1 = n1.cross(b2 / b2n);
    return std::atan2(m1.dot(n2), n1.dot(n2));
}

// Ramachandran-region binning. Lovell-Richardson 2003-aligned grid,
// kept bit-identical to DihedralTimeSeriesTrajectoryResult.cpp's copy
// (science-review HIGH 1-4, 2026-05-19). Downstream consumers reading
// both groups see identical bin labels and boundaries. Returns
// kBinUnassigned for NaN inputs.
//
// Boundaries (degrees, inclusive):
//   αR  : phi ∈ [-180, -30], psi ∈ [-90,  30]
//   β   : phi ∈ [-180, -45], psi ∈ [60, 180] ∪ [-180,-150]
//   αL  : phi ∈ [  30, 100], psi ∈ [-10,  80]
//   PPII: phi ∈ [-75,  -50], psi ∈ [140, 165]  (Berkholz/Adzhubei cone)
// Resolution: αR → αL → PPII → β → other (first match wins).
// References: Lovell 2003 Proteins 50:437; Berkholz 2010 Structure
// 18:1257; Adzhubei 2013 J. Mol. Biol. 425:2100.
std::uint8_t RamachandranBin(double phi_rad, double psi_rad) {
    using R = DihedralBinTransitionTrajectoryResult;
    if (!std::isfinite(phi_rad) || !std::isfinite(psi_rad))
        return R::kBinUnassigned;

    const double deg_per_rad = 180.0 / M_PI;
    const double phi = phi_rad * deg_per_rad;
    const double psi = psi_rad * deg_per_rad;

    if (phi >= -180.0 && phi <= -30.0 &&
        psi >=  -90.0 && psi <=  30.0)
        return R::kBinAlphaR;

    if (phi >=  30.0 && phi <= 100.0 &&
        psi >= -10.0 && psi <=  80.0)
        return R::kBinAlphaL;

    if (phi >= -75.0 && phi <= -50.0 &&
        psi >= 140.0 && psi <= 165.0)
        return R::kBinPPII;

    if (phi >= -180.0 && phi <= -45.0 &&
        ((psi >=  60.0 && psi <= 180.0) ||
         (psi >= -180.0 && psi <= -150.0)))
        return R::kBinBeta;

    return R::kBinOther;
}

// Chi rotamer 120° three-bin classification, Lovell-Richardson 2003
// convention (math-review MED-4, 2026-05-19): the boundaries between
// rotamer wells centered at ±180° and ±60° are at ±120°, and those
// boundary points themselves are conventionally assigned to trans
// (not gauche). This differs from ChiRotamerSelectionTrajectoryResult
// .cpp:42-52 which uses strict `>` boundaries (places ±120° in
// gauche). Both calculators are live in PerFrameExtractionSet; the two
// conventions are deliberately documented-divergent rather than aligned
// (they serve different downstream consumers — bin-transition counting
// here vs DFT-pose rotamer selection there).
//   trans = |chi| ≥ 120°    (chi ≈ ±180°)
//   g+    = 0° < chi < 120° (chi ≈ +60°)
//   g-    = -120° < chi ≤ 0° (chi ≈ -60°)
// chi == 0° exactly lands in g- via the else branch (rare in MD).
// Returns kChiBinUnassigned for NaN inputs.
std::uint8_t ChiRotamerBin(double chi_rad) {
    using R = DihedralBinTransitionTrajectoryResult;
    if (!std::isfinite(chi_rad)) return R::kChiBinUnassigned;
    const double third = 2.0 * M_PI / 3.0;  // 120°
    if (chi_rad >=  third) return R::kChiBinTrans;
    if (chi_rad <= -third) return R::kChiBinTrans;
    if (chi_rad >   0.0)   return R::kChiBinGplus;
    return R::kChiBinGminus;
}

// Argmax over a fixed-size bin-occupancy array. Returns the index of
// the bin with the highest frame count, or `unassigned` when all
// counts are zero (no observed frames for this residue).
template <typename Arr, typename Idx>
Idx ArgmaxOrUnassigned(const Arr& arr, Idx unassigned) {
    std::uint32_t best_count = 0;
    Idx best_idx = unassigned;
    for (std::size_t i = 0; i < arr.size(); ++i) {
        if (arr[i] > best_count) {
            best_count = arr[i];
            best_idx = static_cast<Idx>(i);
        }
    }
    return best_idx;
}

}  // anonymous namespace


std::unique_ptr<DihedralBinTransitionTrajectoryResult>
DihedralBinTransitionTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<DihedralBinTransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();

    r->prev_backbone_bin_.assign(R, kBinUnassigned);
    r->prev_chi_bin_.assign(R, {kChiBinUnassigned, kChiBinUnassigned,
                                kChiBinUnassigned, kChiBinUnassigned});

    r->backbone_transition_count_.assign(R, 0);
    r->backbone_bin_occupancy_.assign(R, {});
    r->n_frames_observed_.assign(R, 0);

    r->chi_transition_count_.assign(R, {});
    r->chi_n_frames_observed_.assign(R, {});
    r->chi_rotamer_occupancy_.assign(R, {});

    r->backbone_dominant_region_.assign(R, kBinUnassigned);
    r->chi_dominant_rotamer_.assign(R, {kChiBinUnassigned, kChiBinUnassigned,
                                        kChiBinUnassigned, kChiBinUnassigned});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────────
//
// Per residue: compute phi (via Predecessor), psi (via Successor), and
// chi[k] from cached Residue.chi[k] atom indices. Classify into bins;
// update running counters in place. Consecutive-frame transition gate:
// previous bin must be observed (non-unassigned) AND current bin must
// be observed; if either is unassigned, no transition is counted but
// the current bin is recorded for the next frame's comparison.

void DihedralBinTransitionTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const Protein& protein = conf.ProteinRef();
    const std::size_t R = protein.ResidueCount();

    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);

        // Backbone: phi/psi via the canonical bond-graph queries.
        const auto prev_idx_opt = protein.BackbonePredecessor(ri);
        const auto next_idx_opt = protein.BackboneSuccessor(ri);

        double phi_val = kNaN;
        if (prev_idx_opt &&
            res.CA != Residue::NONE && res.C != Residue::NONE) {
            const Residue& res_prev = protein.ResidueAt(*prev_idx_opt);
            phi_val = Dihedral(conf.PositionAt(res_prev.C),
                                conf.PositionAt(res.N),
                                conf.PositionAt(res.CA),
                                conf.PositionAt(res.C));
        }
        double psi_val = kNaN;
        if (next_idx_opt &&
            res.N != Residue::NONE && res.CA != Residue::NONE) {
            const Residue& res_next = protein.ResidueAt(*next_idx_opt);
            psi_val = Dihedral(conf.PositionAt(res.N),
                                conf.PositionAt(res.CA),
                                conf.PositionAt(res.C),
                                conf.PositionAt(res_next.N));
        }

        const std::uint8_t cur_bin = RamachandranBin(phi_val, psi_val);
        // Bin 0 (kBinUnassigned) IS populated here (codex review 2026-05-19):
        // backbone_bin_occupancy[:, 0] previously stayed zero forever, but
        // the H5 legend names bin 0 as "unassigned" — so consumers reading
        // occupancy[:, 0] expected an unassigned count and got zero. Now
        // every frame contributes to exactly one bin, and
        // sum(backbone_bin_occupancy[ri, :]) == T for all residues.
        // n_frames_observed[ri] = sum(backbone_bin_occupancy[ri, 1:]).
        ++backbone_bin_occupancy_[ri][cur_bin];
        if (cur_bin != kBinUnassigned) {
            ++n_frames_observed_[ri];

            const std::uint8_t prev_bin = prev_backbone_bin_[ri];
            if (prev_bin != kBinUnassigned && prev_bin != cur_bin) {
                ++backbone_transition_count_[ri];
            }
        }
        // Always update prev (even unassigned-to-observed transition
        // does NOT count, but next frame's compare needs the current
        // bin). When cur_bin is unassigned, prev becomes unassigned so
        // the next consecutive-frame check requires re-observation.
        prev_backbone_bin_[ri] = cur_bin;

        // Chi: bins per defined chi index k.
        for (int k = 0; k < 4; ++k) {
            double chi_k = kNaN;
            if (res.chi[k].Valid()) {
                chi_k = Dihedral(conf.PositionAt(res.chi[k].a[0]),
                                  conf.PositionAt(res.chi[k].a[1]),
                                  conf.PositionAt(res.chi[k].a[2]),
                                  conf.PositionAt(res.chi[k].a[3]));
            }
            const std::uint8_t cur = ChiRotamerBin(chi_k);
            if (cur != kChiBinUnassigned) {
                ++chi_rotamer_occupancy_[ri][k][cur];
                ++chi_n_frames_observed_[ri][k];
                const std::uint8_t prev = prev_chi_bin_[ri][k];
                if (prev != kChiBinUnassigned && prev != cur) {
                    ++chi_transition_count_[ri][k];
                }
            }
            prev_chi_bin_[ri][k] = cur;
        }
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(1u);  // positions always present
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────────
//
// Compute dominant bins per residue / per chi from the accumulated
// occupancy counts. Idempotent: a second call recomputes the same
// answer from the same accumulators.

void DihedralBinTransitionTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                    Trajectory& traj) {
    (void)tp; (void)traj;
    const std::size_t R = backbone_transition_count_.size();
    for (std::size_t ri = 0; ri < R; ++ri) {
        backbone_dominant_region_[ri] = ArgmaxOrUnassigned<
            std::array<std::uint32_t, kBackboneBinCount>, std::uint8_t>(
            backbone_bin_occupancy_[ri], kBinUnassigned);

        for (int k = 0; k < 4; ++k) {
            chi_dominant_rotamer_[ri][k] = ArgmaxOrUnassigned<
                std::array<std::uint32_t, kChiBinCount>, std::uint8_t>(
                chi_rotamer_occupancy_[ri][k], kChiBinUnassigned);
        }
    }
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "DihedralBinTransitionTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(R) + " residues.");
}


// ── WriteH5Group ─────────────────────────────────────────────────────
//
// Per-residue stats + per-residue per-bin occupancy + per-residue per-
// chi per-rotamer occupancy. Frame metadata for cadence + provenance.
// Convention pins (legend, boundaries, resolution order, chi rotamer
// endpoints) match DihedralTS verbatim; consumers see identical
// vocabulary across both groups.

void DihedralBinTransitionTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t R = backbone_transition_count_.size();
    const std::size_t T = n_frames_;
    const std::size_t N = tp.AtomCount();

    auto grp = file.createGroup("/trajectory/dihedral_bin_transition");

    grp.createAttribute("result_name",  Name());
    grp.createAttribute("n_residues",   R);
    grp.createAttribute("n_atoms",      N);
    grp.createAttribute("n_frames",     T);
    grp.createAttribute("finalized",    finalized_);

    grp.createAttribute("backbone_bin_count",     kBackboneBinCount);
    grp.createAttribute("chi_count",              kChiCount);
    grp.createAttribute("chi_bin_count",          kChiBinCount);

    grp.createAttribute("backbone_bin_legend", std::string(
        "0=unassigned (phi or psi NaN at termini / non-bonded gaps), "
        "1=alphaR, 2=beta, 3=alphaL, 4=PPII, 5=other. "
        "Identical labelling to DihedralTimeSeries.rama_region. "
        "Codex-review-2026-05-19 fix: bin 0 IS populated when phi/psi "
        "are NaN, so sum(backbone_bin_occupancy[ri, :]) == n_frames for "
        "all residues. n_frames_observed[ri] = sum(occupancy[ri, 1:])."));
    grp.createAttribute("backbone_bin_boundaries", std::string(
        "alphaR: phi[-180,-30], psi[-90,30]; "
        "beta: phi[-180,-45], psi[60,180]U[-180,-150]; "
        "alphaL: phi[30,100], psi[-10,80]; "
        "PPII: phi[-75,-50], psi[140,165] (tight Berkholz/Adzhubei cone); "
        "boundaries in degrees, inclusive both ends. "
        "Resolution order: alphaR -> alphaL -> PPII -> beta -> other "
        "(first match wins). Matches DihedralTimeSeries verbatim. "
        "References: Lovell 2003 Proteins 50:437; Berkholz 2010 "
        "Structure 18:1257; Adzhubei 2013 J. Mol. Biol. 425:2100."));
    grp.createAttribute("chi_rotamer_legend", std::string(
        "0=gplus (0 < chi < 120 deg), 1=trans (|chi| >= 120 deg), "
        "2=gminus (-120 < chi <= 0 deg), 255=unassigned (chi not "
        "defined for this AA or per-frame geometry degenerate). "
        "Lovell-Richardson 2003 convention: ±120 deg lands in trans "
        "(boundary between wells at ±180 and ±60). chi == 0 exactly "
        "lands in g- via the else branch (rare in MD)."));
    grp.createAttribute("cadence_caveat", std::string(
        "Transition counts are stride-sensitive. At <= 100 ps cadence, "
        "surface-exposed sidechain rotamer transitions are captured; "
        "buried-residue transitions (dwell ~ns-us) require finer cadence "
        "to resolve. Cadence is set by --stride (default 1 = every TRR "
        "frame); coarser strides miss progressively faster transitions "
        "(e.g. stride=300 / 3 ns/frame misses most rotamer transitions)."));
    grp.createAttribute("transition_gate", std::string(
        "Both prev and curr frame must have an observed bin (non-"
        "unassigned) for a transition to count. Wrap-tolerant via bin "
        "labels, not raw angle deltas. Consecutive-frame gate uses "
        "frame-index order; if intermediate frames are unobserved "
        "(NaN bin), the transition chain breaks and re-starts on the "
        "next observation."));
    grp.createAttribute("angle_convention", std::string(
        "IUPAC signed dihedral atan2(y,x). phi = C(i-1)-N(i)-CA(i)-C(i); "
        "psi = N(i)-CA(i)-C(i)-N(i+1); chi from AminoAcidType.chi_angles. "
        "Connectivity via Protein::BackbonePredecessor / "
        "BackboneSuccessor (bond-graph walk, PeptideCN filter); see "
        "DihedralTimeSeries chain_break_policy for the canonical "
        "discipline."));
    grp.createAttribute("residue_axis", std::string("protein_residue_index"));
    grp.createAttribute("source", std::string(
        "positions + AminoAcidType.chi_angles (Residue.chi[k]) + "
        "Protein backbone bond graph. No source ConformationResult "
        "dependency; independent compute (does not cross-read "
        "DihedralTimeSeries)."));
    grp.createAttribute("source_attached_policy", std::string(
        "always_attached -- positions always present at tp.Seed time. "
        "source_attached_per_frame trivially all-1 for SDK uniformity "
        "(OBJECT_MODEL.md Conditional-attach TR discipline)."));

    // ── Per-residue stats (R,) ───────────────────────────────────────
    grp.createDataSet("backbone_transition_count", backbone_transition_count_)
       .createAttribute("units", std::string("transitions"));
    grp.createDataSet("backbone_dominant_region", backbone_dominant_region_)
       .createAttribute("units", std::string("category"));
    grp.createDataSet("n_frames_observed", n_frames_observed_)
       .createAttribute("units", std::string("frame_count"));

    // ── Per-residue per backbone bin (R, 6) ──────────────────────────
    {
        std::vector<std::uint32_t> flat(R * kBackboneBinCount);
        for (std::size_t ri = 0; ri < R; ++ri) {
            for (std::size_t b = 0; b < kBackboneBinCount; ++b) {
                flat[ri * kBackboneBinCount + b] =
                    backbone_bin_occupancy_[ri][b];
            }
        }
        const std::vector<std::size_t> dims = {R, kBackboneBinCount};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<std::uint32_t>(
            "backbone_bin_occupancy", space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("frame_count"));
    }

    // ── Per-residue per chi (R, 4) ───────────────────────────────────
    {
        std::vector<std::uint32_t> flat_t(R * kChiCount);
        std::vector<std::uint8_t>  flat_d(R * kChiCount);
        std::vector<std::uint32_t> flat_n(R * kChiCount);
        for (std::size_t ri = 0; ri < R; ++ri) {
            for (std::size_t k = 0; k < kChiCount; ++k) {
                flat_t[ri * kChiCount + k] = chi_transition_count_[ri][k];
                flat_d[ri * kChiCount + k] = chi_dominant_rotamer_[ri][k];
                flat_n[ri * kChiCount + k] = chi_n_frames_observed_[ri][k];
            }
        }
        const std::vector<std::size_t> dims = {R, kChiCount};
        HighFive::DataSpace s(dims);
        auto dt = grp.createDataSet<std::uint32_t>(
            "chi_transition_count", s);
        dt.write_raw(flat_t.data());
        dt.createAttribute("units", std::string("transitions"));

        auto dd = grp.createDataSet<std::uint8_t>(
            "chi_dominant_rotamer", s);
        dd.write_raw(flat_d.data());
        dd.createAttribute("units", std::string("category"));
        dd.createAttribute("note", std::string(
            "255 = no observation for this chi (chi not defined for "
            "this residue's AA, or no finite chi observed in any frame)."));

        auto dn = grp.createDataSet<std::uint32_t>(
            "chi_n_frames_observed", s);
        dn.write_raw(flat_n.data());
        dn.createAttribute("units", std::string("frame_count"));
    }

    // ── Per-residue per chi per bin (R, 4, 3) ────────────────────────
    {
        std::vector<std::uint32_t> flat(R * kChiCount * kChiBinCount);
        for (std::size_t ri = 0; ri < R; ++ri) {
            for (std::size_t k = 0; k < kChiCount; ++k) {
                for (std::size_t b = 0; b < kChiBinCount; ++b) {
                    flat[(ri * kChiCount + k) * kChiBinCount + b] =
                        chi_rotamer_occupancy_[ri][k][b];
                }
            }
        }
        const std::vector<std::size_t> dims = {R, kChiCount, kChiBinCount};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<std::uint32_t>(
            "chi_rotamer_occupancy", space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("frame_count"));
    }

    // ── Per-frame metadata (T,) ──────────────────────────────────────
    grp.createDataSet("frame_indices", frame_indices_)
       .createAttribute("units", std::string("frame_index"));
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));
    auto ds_attached = grp.createDataSet(
        "source_attached_per_frame", source_attached_per_frame_);
    ds_attached.createAttribute("units", std::string("dimensionless"));

    OperationLog::Info(LogCalcOther,
        "DihedralBinTransitionTrajectoryResult::WriteH5Group",
        "wrote /trajectory/dihedral_bin_transition with " +
        std::to_string(R) + " residues x " + std::to_string(T) + " frames");
}


}  // namespace nmr
