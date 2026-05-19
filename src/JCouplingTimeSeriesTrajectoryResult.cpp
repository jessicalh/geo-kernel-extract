#include "JCouplingTimeSeriesTrajectoryResult.h"

#include "Atom.h"
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

#include <cmath>
#include <limits>
#include <string>
#include <typeinfo>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Karplus coefficients (A, B, C) per channel. Karplus form:
//   3J(θ) = A·cos²(θ) + B·cos(θ) + C
// θ in radians via IUPAC signed dihedral atan2. References pinned in
// the header and the H5 attrs.

// Wang & Bax 1996 JACS 118:2483 — 3J(HN, Hα).
constexpr double kJHNHAlpha_A =  6.51;
constexpr double kJHNHAlpha_B = -1.76;
constexpr double kJHNHAlpha_C =  1.60;

// Pérez et al. 2001 JACS 123:7081 — 3J(N, Cγ).
constexpr double kJNCgamma_A =  1.29;
constexpr double kJNCgamma_B = -0.49;
constexpr double kJNCgamma_C =  0.37;

// Pérez et al. 2001 — 3J(C', Cγ).
constexpr double kJCprimeCgamma_A =  1.74;
constexpr double kJCprimeCgamma_B = -0.57;
constexpr double kJCprimeCgamma_C =  0.25;


// Dihedral helper — same atan2-based formulation as DihedralTimeSeries
// .cpp:45, ChiRotamerSelection, PlanarGeometry, etc. NaN-on-degenerate
// (collinear b2 or zero-norm n1/n2). Per PATTERNS Lesson 10
// (equation in comment, no utility namespace).
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


// Karplus: 3J(θ) = A·cos²(θ) + B·cos(θ) + C. θ in radians.
double Karplus(double A, double B, double C, double theta_rad) {
    if (!std::isfinite(theta_rad)) return kNaN;
    const double cos_theta = std::cos(theta_rad);
    return A * cos_theta * cos_theta + B * cos_theta + C;
}

}  // anonymous namespace


std::unique_ptr<JCouplingTimeSeriesTrajectoryResult>
JCouplingTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<JCouplingTimeSeriesTrajectoryResult>();
    const Protein& protein = tp.ProteinRef();
    const std::size_t R = protein.ResidueCount();
    const std::size_t N = tp.AtomCount();

    r->j_hn_halpha_.assign(R, {});
    r->j_n_cgamma_.assign(R, {});
    r->j_cprime_cgamma_.assign(R, {});

    // Static masks: which channels can structurally exist per residue.
    r->j_hn_halpha_exists_.assign(R, 0);
    r->j_chi1_exists_.assign(R, 0);
    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.H  != Residue::NONE && res.N  != Residue::NONE &&
            res.CA != Residue::NONE && res.HA != Residue::NONE) {
            r->j_hn_halpha_exists_[ri] = 1u;
        }
        if (res.chi[0].Valid()) {
            r->j_chi1_exists_[ri] = 1u;
        }
    }

    // Atom-axis broadcast lookup.
    r->residue_index_per_atom_.assign(N, -1);
    for (std::size_t ai = 0; ai < N; ++ai) {
        r->residue_index_per_atom_[ai] =
            static_cast<std::int32_t>(protein.AtomAt(ai).residue_index);
    }
    return r;
}


// ── Compute ──────────────────────────────────────────────────────────
//
// Per residue, per frame: compute the three Karplus 3J observables.
// Each channel NaN-fills where the structural atoms are missing
// (PRO: no H; GLY/ALA: no chi1).

void JCouplingTimeSeriesTrajectoryResult::Compute(
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

        // ── J(HN, Hα) via H-N-CA-HA ──────────────────────────────────
        double j_hn = kNaN;
        if (j_hn_halpha_exists_[ri]) {
            const double theta = Dihedral(
                conf.PositionAt(res.H),
                conf.PositionAt(res.N),
                conf.PositionAt(res.CA),
                conf.PositionAt(res.HA));
            j_hn = Karplus(kJHNHAlpha_A, kJHNHAlpha_B, kJHNHAlpha_C, theta);
        }
        j_hn_halpha_[ri].push_back(j_hn);

        // ── J(N, Cγ) via N-CA-CB-CG = chi1 ───────────────────────────
        // chi[0].a = (N, CA, CB, CG_first); use directly.
        double j_n_cg = kNaN;
        if (j_chi1_exists_[ri]) {
            const auto& c1 = res.chi[0];
            const double theta = Dihedral(
                conf.PositionAt(c1.a[0]),
                conf.PositionAt(c1.a[1]),
                conf.PositionAt(c1.a[2]),
                conf.PositionAt(c1.a[3]));
            j_n_cg = Karplus(kJNCgamma_A, kJNCgamma_B, kJNCgamma_C, theta);
        }
        j_n_cgamma_[ri].push_back(j_n_cg);

        // ── J(C', Cγ) via C-CA-CB-CG ──────────────────────────────────
        // Replace chi1's first atom (N) with res.C; preserve CB, CG.
        double j_cp_cg = kNaN;
        if (j_chi1_exists_[ri] && res.C != Residue::NONE &&
            res.CA != Residue::NONE) {
            const auto& c1 = res.chi[0];
            const double theta = Dihedral(
                conf.PositionAt(res.C),
                conf.PositionAt(res.CA),
                conf.PositionAt(c1.a[2]),
                conf.PositionAt(c1.a[3]));
            j_cp_cg = Karplus(kJCprimeCgamma_A, kJCprimeCgamma_B,
                              kJCprimeCgamma_C, theta);
        }
        j_cprime_cgamma_[ri].push_back(j_cp_cg);
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(1u);  // positions always present
    ++n_frames_;
}


void JCouplingTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                    Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "JCouplingTimeSeriesTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(j_hn_halpha_.size()) +
        " residues.");
}


void JCouplingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t R = j_hn_halpha_.size();
    const std::size_t T = n_frames_;
    const std::size_t N = tp.AtomCount();

    auto grp = file.createGroup("/trajectory/j_coupling_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_residues",  R);
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("karplus_form", std::string(
        "3J(theta) = A * cos^2(theta) + B * cos(theta) + C; "
        "theta in radians via IUPAC signed dihedral atan2."));
    grp.createAttribute("J_HN_Halpha_coefficients", std::string(
        "A=6.51, B=-1.76, C=1.60 (Wang & Bax 1996 JACS 118:2483); "
        "dihedral H-N-CA-HA."));
    grp.createAttribute("J_N_Cgamma_coefficients", std::string(
        "A=1.29, B=-0.49, C=0.37 (Perez et al. 2001 JACS 123:7081); "
        "dihedral N-CA-CB-CG (= chi1)."));
    grp.createAttribute("J_Cprime_Cgamma_coefficients", std::string(
        "A=1.74, B=-0.57, C=0.25 (Perez et al. 2001 JACS 123:7081); "
        "dihedral C-CA-CB-CG (chi1 with C' as the leading atom; 120 "
        "degrees offset from chi1 on the same Cbeta-Calpha axis)."));
    grp.createAttribute("dihedral_convention", std::string(
        "IUPAC signed dihedral atan2(y,x). H-N-CA-HA for J(HN,Halpha); "
        "N-CA-CB-CG for J(N,Cgamma); C-CA-CB-CG for J(C',Cgamma)."));
    grp.createAttribute("GLY_caveat", std::string(
        "GLY Halpha uses Residue.HA which is HA2 by Residue.h cache "
        "convention. HA3 is NOT separately measured here; consumers "
        "needing pro-R/pro-S resolution should compute directly via "
        "the residue's two Halpha atom indices. The slight averaging "
        "error vs GLY's two Halpha is documented; calibration absorbs."));
    grp.createAttribute("units",       std::string("Hz"));
    grp.createAttribute("absent_sentinel", std::string("NaN"));
    grp.createAttribute("residue_axis", std::string("protein_residue_index"));
    grp.createAttribute("atom_axis",    std::string("protein_atom_index"));
    grp.createAttribute("source", std::string(
        "positions + Residue backbone-cache (H, N, CA, HA, C) + "
        "Residue.chi[0] (chi1 atom indices from AminoAcidType.chi_"
        "angles). No source ConformationResult dependency; positions "
        "always present at tp.Seed time."));
    grp.createAttribute("source_attached_policy", std::string(
        "always_attached -- source_attached_per_frame trivially all-1 "
        "for SDK uniformity (OBJECT_MODEL.md Conditional-attach TR "
        "discipline)."));

    // ── Per-frame (T,) ───────────────────────────────────────────────
    grp.createDataSet("frame_indices", frame_indices_)
       .createAttribute("units", std::string("frame_index"));
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_)
       .createAttribute("units", std::string("dimensionless"));

    // ── Helper: emit (R, T) flat float64 dataset, NaN-fill ───────────
    auto emit_2d_f64 = [&](const std::string& name,
                            const std::vector<std::vector<double>>& src) {
        std::vector<double> flat(R * T, kNaN);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto& row = src[ri];
            for (std::size_t f = 0; f < T && f < row.size(); ++f) {
                flat[ri * T + f] = row[f];
            }
        }
        const std::vector<std::size_t> dims = {R, T};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("Hz"));
    };

    emit_2d_f64("J_HN_Halpha",     j_hn_halpha_);
    emit_2d_f64("J_N_Cgamma",      j_n_cgamma_);
    emit_2d_f64("J_Cprime_Cgamma", j_cprime_cgamma_);

    // ── Static per-residue masks (R,) ────────────────────────────────
    grp.createDataSet("J_HN_Halpha_exists", j_hn_halpha_exists_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("J_chi1_exists", j_chi1_exists_)
       .createAttribute("units", std::string("dimensionless"));

    // ── Per-atom lookup (N,) ─────────────────────────────────────────
    grp.createDataSet("residue_index_per_atom", residue_index_per_atom_)
       .createAttribute("units", std::string("residue_index"));

    OperationLog::Info(LogCalcOther,
        "JCouplingTimeSeriesTrajectoryResult::WriteH5Group",
        "wrote /trajectory/j_coupling_time_series with " +
        std::to_string(R) + " residues x " + std::to_string(T) + " frames");
}


}  // namespace nmr
