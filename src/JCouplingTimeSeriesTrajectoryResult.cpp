#include "JCouplingTimeSeriesTrajectoryResult.h"

#include "Atom.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
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

// Karplus coefficients live in PhysicalConstants.h with full literature
// citations + reference PDFs (cited per-symbol). Karplus form for the
// nine emitted datasets across eight channel families:
//   3J(theta) = A * cos^2(theta) + B * cos(theta) + C   [Hz]
// Backbone channels use project_phi + project_theta_offset; chi1
// channels use the actual atomic 4-atom dihedral. See
// PhysicalConstants.h "Karplus 3J-coupling parameters" for the
// per-channel convention and arithmetic bounds.


// Dihedral helper — same atan2-based formulation as DihedralTimeSeries,
// ChiRotamerSelection, PlanarGeometry, etc. NaN-on-degenerate
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
    r->j_hn_halpha_vogeli_.assign(R, {});
    r->j_hn_cbeta_.assign(R, {});
    r->j_hn_cprime_.assign(R, {});
    r->j_halpha_cprime_.assign(R, {});
    r->j_n_cgamma_.assign(R, {});
    r->j_cprime_cgamma_.assign(R, {});
    r->j_halpha_hbeta2_.assign(R, {});
    r->j_halpha_hbeta3_.assign(R, {});

    // Static masks: which channels can structurally exist per residue.
    r->j_hn_halpha_exists_.assign(R, 0);
    r->j_hn_cbeta_exists_.assign(R, 0);
    r->j_hn_cprime_exists_.assign(R, 0);
    r->j_halpha_cprime_exists_.assign(R, 0);
    r->j_chi1_exists_.assign(R, 0);
    r->j_n_cgamma_exists_.assign(R, 0);
    r->j_cprime_cgamma_exists_.assign(R, 0);
    r->j_halpha_hbeta_exists_.assign(R, 0);
    r->hb2_index_.assign(R, Residue::NONE);
    r->hb3_index_.assign(R, Residue::NONE);
    r->c_prev_index_.assign(R, Residue::NONE);

    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);

        // C(prev): canonical bond-graph predecessor lookup. The
        // backbone phi-derived Karplus channels need a defined
        // C(prev)-N-CA-C phi. Halpha-C' also needs C(prev) for the
        // actual 3-bond path across the peptide bond.
        if (const auto prev_idx_opt = protein.BackbonePredecessor(ri)) {
            const Residue& prev = protein.ResidueAt(*prev_idx_opt);
            if (prev.C != Residue::NONE) {
                r->c_prev_index_[ri] = prev.C;
            }
        }

        if (res.H  != Residue::NONE && res.N  != Residue::NONE &&
            res.CA != Residue::NONE && res.C  != Residue::NONE &&
            res.HA != Residue::NONE &&
            r->c_prev_index_[ri] != Residue::NONE) {
            r->j_hn_halpha_exists_[ri] = 1u;
        }
        if (res.H  != Residue::NONE && res.N  != Residue::NONE &&
            res.CA != Residue::NONE && res.C  != Residue::NONE &&
            res.CB != Residue::NONE &&
            r->c_prev_index_[ri] != Residue::NONE) {
            r->j_hn_cbeta_exists_[ri] = 1u;
        }
        if (res.H  != Residue::NONE && res.N  != Residue::NONE &&
            res.CA != Residue::NONE && res.C  != Residue::NONE &&
            r->c_prev_index_[ri] != Residue::NONE) {
            r->j_hn_cprime_exists_[ri] = 1u;
        }
        if (res.chi[0].Valid()) {
            r->j_chi1_exists_[ri] = 1u;
            // chi1-existence is necessary but NOT sufficient for the
            // J(N,Cγ) / J(C',Cγ) channels: chi[0].a[3] is the chi1
            // terminal which is Cγ only for residues with a carbon
            // chi1 terminal. SER (OG), CYS (SG), THR (OG1) carry
            // non-carbon chi1 terminals -- the Pérez 2001 J(N',Cγ)
            // table technically tabulates per-residue coefficients
            // for those (page 7086 SER / CYS / THR rows) but the
            // consensus row used here is intended for a true carbon
            // Cγ. Gate on Element::C at chi[0].a[3]. THR specifically:
            // chi1 terminal is OG1, not the methyl CG2; the per-channel
            // CG2 J value is not measured by this 4-atom path.
            if (protein.AtomAt(res.chi[0].a[3]).element == Element::C) {
                r->j_n_cgamma_exists_[ri] = 1u;
                if (res.C != Residue::NONE && res.CA != Residue::NONE) {
                    r->j_cprime_cgamma_exists_[ri] = 1u;
                }
            }
        }
        if (res.HA != Residue::NONE && res.CA != Residue::NONE &&
            res.N  != Residue::NONE && res.C  != Residue::NONE &&
            r->c_prev_index_[ri] != Residue::NONE) {
            r->j_halpha_cprime_exists_[ri] = 1u;
        }

        // Hbeta lookup: walk this residue's atoms by pdb_atom_name.
        // Most residues: prochiral methylene pair HB2 + HB3.
        // Ile/Val/Thr: single Hbeta methine ("HB"); emit it in BOTH
        // slots so the mask + downstream consumer is uniform.
        // Ala: methyl HB1/HB2/HB3; deliberately leave HB2 + HB3 caches
        // as NONE so both channels NaN-fill (the methyl 3J is not the
        // same chemical observable as the methylene 3J).
        // Gly: no Cbeta; leave NONE.
        if (res.CB != Residue::NONE && res.HA != Residue::NONE &&
            res.CA != Residue::NONE) {
            std::size_t single_hb = Residue::NONE;
            std::size_t hb2 = Residue::NONE;
            std::size_t hb3 = Residue::NONE;
            std::size_t ala_methyl_count = 0;
            for (std::size_t ai : res.atom_indices) {
                const std::string& name = protein.AtomAt(ai).pdb_atom_name;
                if      (name == "HB")  single_hb = ai;
                else if (name == "HB2") hb2 = ai;
                else if (name == "HB3") hb3 = ai;
                else if (name == "HB1") ++ala_methyl_count;
            }
            // Ile/Val/Thr methine path: a single "HB" atom and no
            // HB1/HB2/HB3 methylene names. Mirror to both slots so
            // mask is uniform.
            if (single_hb != Residue::NONE &&
                hb2 == Residue::NONE && hb3 == Residue::NONE) {
                r->hb2_index_[ri] = single_hb;
                r->hb3_index_[ri] = single_hb;
            } else if (ala_methyl_count == 0 &&
                       hb2 != Residue::NONE && hb3 != Residue::NONE) {
                // Prochiral methylene Hbeta pair (Ser/Cys/Asp/Asn/...
                // also Phe/Tyr/Trp/His/Leu/Met/Glu/Gln/Arg/Lys/Pro).
                r->hb2_index_[ri] = hb2;
                r->hb3_index_[ri] = hb3;
            }
            // Ala (3 methyl HBs): leave both slots NONE so the channel
            // emits NaN; documented in the H5 attr.
        }
        if (r->hb2_index_[ri] != Residue::NONE ||
            r->hb3_index_[ri] != Residue::NONE) {
            r->j_halpha_hbeta_exists_[ri] = 1u;
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
// Per residue, per frame: compute the eight Karplus 3J channel
// families (nine datasets, including the Hbeta methylene pair).
// Each channel NaN-fills where the structural atoms are missing
// (PRO: no H; GLY: no Cβ; GLY/ALA: no chi1; N-terminus: no C(prev);
// SER/CYS/THR: chi1 terminal is non-carbon, so J(N,Cγ) / J(C',Cγ)
// NaN; ALA Hβ deliberately excluded — methyl Cβ is not the methylene
// observable).
//
// Convention -- per codex F6 plus project-sign repair (2026-05-20):
// The four BACKBONE channels (HN-Hα, HN-Cβ, HN-C', Hα-C') feed
// `(phi + theta_offset)` into the Karplus equation, where phi is the
// canonical Ramachandran C(prev)-N-CA-C dihedral computed via
// IUPAC signed atan2 directly from positions, and theta_offset is
// the per-channel KARPLUS_*_THETA constant from PhysicalConstants.h.
// The project phi sign is opposite to the Wang-Bax / Vogeli plotting
// convention, so PhysicalConstants.h stores theta_offset in the
// PROJECT convention (HN-Hα=+pi/3, HN-Cβ=-pi/3, HN-C'=0,
// Hα-C'=+pi/3).
//
// The four CHI1 channels (J(N,Cγ), J(C',Cγ), J(Hα,Hβ{2,3})) feed
// the actual atomic 4-atom dihedral computed from positions (their
// theta_offset = 0). Pérez 2001 Table 2 footnote c states the per-
// coupling (A, B, C) internalize the substituent rotation offset
// around Cα, so feeding the atomic dihedral matches the Table 2
// consensus row directly.

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

        // Compute canonical Ramachandran phi = C(prev)-N-CA-C once
        // per residue; used by all four backbone channels.
        // BackbonePredecessor handles the typed bond-graph query
        // (correct on chain start, ACE/NME caps, insertion codes,
        // cyclic peptides). NaN at chain N-terminus.
        double phi = kNaN;
        if (res.N != Residue::NONE && res.CA != Residue::NONE &&
            res.C != Residue::NONE && c_prev_index_[ri] != Residue::NONE) {
            phi = Dihedral(
                conf.PositionAt(c_prev_index_[ri]),
                conf.PositionAt(res.N),
                conf.PositionAt(res.CA),
                conf.PositionAt(res.C));
        }

        // ── J(HN, Hα) via Karplus(phi + theta) — Vuister & Bax 1993 ──
        // theta_project=+pi/3 (Wang-Bax theta_pub=-pi/3; project phi
        // sign is opposite).
        double j_hn = kNaN;
        double j_hn_vogeli = kNaN;
        if (j_hn_halpha_exists_[ri] && std::isfinite(phi)) {
            const double th = phi + KARPLUS_HN_HA_THETA;
            j_hn = Karplus(KARPLUS_HN_HA_A, KARPLUS_HN_HA_B,
                           KARPLUS_HN_HA_C, th);
            // Same project phi + theta; alternate parametrization
            // (methods-accumulate).
            const double th_v = phi + KARPLUS_HN_HA_VOGELI_THETA;
            j_hn_vogeli = Karplus(KARPLUS_HN_HA_VOGELI_A,
                                  KARPLUS_HN_HA_VOGELI_B,
                                  KARPLUS_HN_HA_VOGELI_C, th_v);
        }
        j_hn_halpha_[ri].push_back(j_hn);
        j_hn_halpha_vogeli_[ri].push_back(j_hn_vogeli);

        // ── J(HN, Cβ) via Karplus(phi + theta) — Wang & Bax 1996 ─────
        // theta_project=-pi/3 (Wang-Bax theta_pub=+pi/3).
        double j_hn_cb = kNaN;
        if (j_hn_cbeta_exists_[ri] && std::isfinite(phi)) {
            const double th = phi + KARPLUS_HN_CB_THETA;
            j_hn_cb = Karplus(KARPLUS_HN_CB_A, KARPLUS_HN_CB_B,
                              KARPLUS_HN_CB_C, th);
        }
        j_hn_cbeta_[ri].push_back(j_hn_cb);

        // ── J(HN, C') via Karplus(phi + theta) — Wang-Bax 1996 row 4 ─
        // theta = 0 (Wang-Bax row 4); equivalent Vogeli form would use
        // eta = pi and B with the opposite sign, but Wang-Bax's
        // theta=0 / B>0 form is the one shipped.
        double j_hn_cp = kNaN;
        if (j_hn_cprime_exists_[ri] && std::isfinite(phi)) {
            const double th = phi + KARPLUS_HN_CP_THETA;
            j_hn_cp = Karplus(KARPLUS_HN_CP_A, KARPLUS_HN_CP_B,
                              KARPLUS_HN_CP_C, th);
        }
        j_hn_cprime_[ri].push_back(j_hn_cp);

        // ── J(Hα, C') via Karplus(phi + theta) — Wang-Bax 1996 row 2 ─
        // theta_project=+pi/3 (Wang-Bax theta_pub=-pi/3). The 3-bond path is
        // Halpha-CA-N-C'(prev), rotation around N-CA (phi axis),
        // per Vuister teaching lecture sect 6.1 + Vogeli 2007 page
        // 9384. We evaluate against project phi directly with the
        // project-sign offset; no need to compute the atomic dihedral.
        // NaN at N-terminus (phi NaN due to no C(prev)).
        double j_ha_cp = kNaN;
        if (j_halpha_cprime_exists_[ri] && std::isfinite(phi)) {
            const double th = phi + KARPLUS_HA_CP_THETA;
            j_ha_cp = Karplus(KARPLUS_HA_CP_A, KARPLUS_HA_CP_B,
                              KARPLUS_HA_CP_C, th);
        }
        j_halpha_cprime_[ri].push_back(j_ha_cp);

        // ── J(N, Cγ) via N-CA-CB-CG = chi1 — Pérez 2001 ──────────────
        // chi[0].a = (N, CA, CB, chi1-terminal). For ARG/ASN/ASP/GLN/
        // GLU/HIS/LEU/LYS/MET/PHE/PRO/TRP/TYR the chi1 terminal IS Cγ.
        // For ILE/VAL the chi1 terminal is CG1 (also carbon). For SER
        // the chi1 terminal is OG (O), for CYS SG (S), for THR OG1 (O)
        // -- those residues' j_n_cgamma_exists_ is 0 by element gate
        // (chi[0].a[3] element != C). Pérez 2001 page 7086 Table 2
        // does tabulate per-residue rows for SER/CYS/THR with their
        // own coefficients (treating the heteroatomic terminal as
        // "Cγ" in the table's residue-keyed parametrization). We
        // omit those channels here -- they are a different chemical
        // observable than the carbon-Cγ consensus row, and shipping a
        // separate "chi1 terminal" channel for them is a methods-
        // accumulate decision pending a calibration trigger (memory
        // entry feedback_methods_accumulate).
        double j_n_cg = kNaN;
        if (j_n_cgamma_exists_[ri]) {
            const auto& c1 = res.chi[0];
            const double theta = Dihedral(
                conf.PositionAt(c1.a[0]),
                conf.PositionAt(c1.a[1]),
                conf.PositionAt(c1.a[2]),
                conf.PositionAt(c1.a[3]));
            j_n_cg = Karplus(KARPLUS_N_CG_A, KARPLUS_N_CG_B,
                             KARPLUS_N_CG_C, theta);
        }
        j_n_cgamma_[ri].push_back(j_n_cg);

        // ── J(C', Cγ) via C-CA-CB-CG — Pérez 2001 ────────────────────
        // Replace chi1's first atom (N) with res.C; preserve CB, CG.
        // Gated by j_cprime_cgamma_exists_ which already requires
        // chi[0].a[3] = Element::C (NaN's SER/CYS/THR for the same
        // reason as J(N, Cγ)).
        double j_cp_cg = kNaN;
        if (j_cprime_cgamma_exists_[ri]) {
            const auto& c1 = res.chi[0];
            const double theta = Dihedral(
                conf.PositionAt(res.C),
                conf.PositionAt(res.CA),
                conf.PositionAt(c1.a[2]),
                conf.PositionAt(c1.a[3]));
            j_cp_cg = Karplus(KARPLUS_CP_CG_A, KARPLUS_CP_CG_B,
                              KARPLUS_CP_CG_C, theta);
        }
        j_cprime_cgamma_[ri].push_back(j_cp_cg);

        // ── J(Hα, Hβ2 / Hβ3) via HA-CA-CB-HB{2,3} — Pérez 2001 ──────
        double j_ha_hb2 = kNaN;
        double j_ha_hb3 = kNaN;
        if (j_halpha_hbeta_exists_[ri]) {
            if (hb2_index_[ri] != Residue::NONE) {
                const double theta = Dihedral(
                    conf.PositionAt(res.HA),
                    conf.PositionAt(res.CA),
                    conf.PositionAt(res.CB),
                    conf.PositionAt(hb2_index_[ri]));
                j_ha_hb2 = Karplus(KARPLUS_HA_HB_A, KARPLUS_HA_HB_B,
                                   KARPLUS_HA_HB_C, theta);
            }
            if (hb3_index_[ri] != Residue::NONE) {
                const double theta = Dihedral(
                    conf.PositionAt(res.HA),
                    conf.PositionAt(res.CA),
                    conf.PositionAt(res.CB),
                    conf.PositionAt(hb3_index_[ri]));
                j_ha_hb3 = Karplus(KARPLUS_HA_HB_A, KARPLUS_HA_HB_B,
                                   KARPLUS_HA_HB_C, theta);
            }
        }
        j_halpha_hbeta2_[ri].push_back(j_ha_hb2);
        j_halpha_hbeta3_[ri].push_back(j_ha_hb3);
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
        "Backbone channels (HN-Halpha, HN-Cbeta, HN-C', Halpha-C'): "
        "3J(phi + theta_offset) = A * cos^2(phi + theta_offset) + B * "
        "cos(phi + theta_offset) + C; phi = canonical Ramachandran "
        "C(prev)-N-CA-C dihedral (radians, IUPAC signed atan2 directly "
        "from positions); theta_offset = per-channel project-convention "
        "constant (KARPLUS_*_THETA in src/PhysicalConstants.h: "
        "HN-Halpha=+pi/3, HN-Cbeta=-pi/3, HN-C'=0, Halpha-C'=+pi/3; "
        "these are the negatives of the Wang-Bax/Vogeli printed offsets "
        "where the project phi sign is opposite). "
        "Chi1 channels (J(N,Cgamma), J(C',Cgamma), J(Halpha,Hbeta)): "
        "3J(alpha) = A * cos^2(alpha) + B * cos(alpha) + C; alpha = "
        "the actual 4-atom atomic dihedral (radians, IUPAC signed "
        "atan2 directly from positions); Perez 2001 Table 2 footnote c "
        "internalizes the Cα-substituent offset in the per-coupling "
        "(A, B, C), so feeding the atomic dihedral matches the "
        "consensus row directly. Backbone form repaired on 2026-05-20 "
        "(codex F6 + project-sign fix) and checked on 1UBQ by "
        "LiteratureAnchoredProbeOn1UBQ. Coefficients defined in "
        "src/PhysicalConstants.h."));

    // Coefficient strings are formatted from the PhysicalConstants.h
    // values at write time so the H5 attrs stay in lockstep with the
    // compiled-in numerics. If a value here ever disagrees with what
    // PhysicalConstants.h declares, the linker is the adversary
    // (feedback_compiler_check).
    auto coeff_str = [](double A, double B, double C) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "A=%.4f, B=%.4f, C=%.4f",
                      A, B, C);
        return std::string(buf);
    };
    grp.createAttribute("J_HN_Halpha_coefficients",
        coeff_str(KARPLUS_HN_HA_A, KARPLUS_HN_HA_B, KARPLUS_HN_HA_C) +
        ", theta_offset_project=+pi/3 (Vuister & Bax 1993 JACS 115:7772, "
        "DOI 10.1021/ja00070a024); evaluated as J = A*cos^2(phi + "
        "theta_offset) + B*cos(phi + theta_offset) + C. Reference PDF "
        "references/vuister-lecture-j-couplings.pdf.");
    grp.createAttribute("J_HN_Halpha_Vogeli_coefficients",
        coeff_str(KARPLUS_HN_HA_VOGELI_A, KARPLUS_HN_HA_VOGELI_B,
                  KARPLUS_HN_HA_VOGELI_C) +
        ", theta_offset_project=+pi/3 (Vogeli, Ying, Grishaev & Bax 2007 "
        "JACS 129:9377, DOI 10.1021/ja070324o, Table 1 'rigid' row); "
        "same project theta_offset as J_HN_Halpha (Vogeli eq 5 prints "
        "eta_ik=-pi/3 in the opposite phi convention). Methods-accumulate "
        "alternate parametrization. "
        "Reference PDF references/vogeli-2007-limits-backbone-dynamics-"
        "3j-couplings-gb3.pdf (byte-verified 2026-05-19).");
    grp.createAttribute("J_HN_Cbeta_coefficients",
        coeff_str(KARPLUS_HN_CB_A, KARPLUS_HN_CB_B, KARPLUS_HN_CB_C) +
        ", theta_offset_project=-pi/3 (Wang & Bax 1996 JACS 118:2483, DOI "
        "10.1021/ja9535524, Table 1 NMR/X-ray refined fit row 3, "
        "theta=+60 deg). Reference PDF "
        "references/wang-bax-1996-karplus-phi-ubiquitin.pdf "
        "(byte-verified 2026-05-19).");
    grp.createAttribute("J_HN_Cprime_coefficients",
        coeff_str(KARPLUS_HN_CP_A, KARPLUS_HN_CP_B, KARPLUS_HN_CP_C) +
        ", theta_offset=0 (Wang & Bax 1996 JACS 118:2483, DOI "
        "10.1021/ja9535524, Table 1 NMR/X-ray refined fit ROW 4, "
        "theta=0 deg). Note B is POSITIVE; J can be slightly negative "
        "at the vertex (~-0.04 Hz; physical). Reference PDF "
        "references/wang-bax-1996-karplus-phi-ubiquitin.pdf "
        "(byte-verified 2026-05-19; row-mapping fixed 2026-05-20 per "
        "codex F1 -- prior bundle attributed Wang-Bax row 2 (Halpha-C') "
        "values (3.75, +2.19, 1.28) to this channel).");
    grp.createAttribute("J_Halpha_Cprime_coefficients",
        coeff_str(KARPLUS_HA_CP_A, KARPLUS_HA_CP_B, KARPLUS_HA_CP_C) +
        ", theta_offset_project=+pi/3 (Wang & Bax 1996 JACS 118:2483, DOI "
        "10.1021/ja9535524, Table 1 NMR/X-ray refined fit ROW 2, "
        "theta=-60 deg). 3-bond path Halpha-CA-N-C'(prev) crosses the "
        "peptide bond at the previous-residue's C'; rotation around "
        "N-CA (phi axis). Note B is POSITIVE -- arithmetic max is "
        "f(+1)=A+B+C, vertex (MIN) at u*=-B/(2A). Reference PDFs "
        "references/wang-bax-1996-karplus-phi-ubiquitin.pdf "
        "(byte-verified 2026-05-19) + Vuister teaching lecture "
        "sect 6.1 + Vogeli 2007 page 9384 (3J(C'(i-1), Halpha) "
        "is one of the six phi-related couplings). Row-mapping + "
        "atom-path fixed 2026-05-20 per codex F1 + F2; prior bundle "
        "attributed row 4 values + used HA-CA-C-N(next) (psi axis, "
        "wrong observable). Eval form switched to (phi + theta_offset) "
        "per codex F6 same date.");
    grp.createAttribute("J_N_Cgamma_coefficients",
        coeff_str(KARPLUS_N_CG_A, KARPLUS_N_CG_B, KARPLUS_N_CG_C) +
        " (Perez, Lohr, Ruterjans & Schmidt 2001 JACS 123:7081, "
        "DOI 10.1021/ja003724j, Table 2 consensus row); dihedral "
        "N-CA-CB-CG (= chi1). Reference PDF "
        "references/perez-2001-self-consistent-karplus-3j-chi1.pdf "
        "(byte-verified 2026-05-19, page 7086).");
    grp.createAttribute("J_Cprime_Cgamma_coefficients",
        coeff_str(KARPLUS_CP_CG_A, KARPLUS_CP_CG_B, KARPLUS_CP_CG_C) +
        " (Perez, Lohr, Ruterjans & Schmidt 2001 JACS 123:7081, "
        "DOI 10.1021/ja003724j, Table 2 consensus row); dihedral "
        "C-CA-CB-CG. (Prior commit carried (1.74, -0.57, 0.25); "
        "corrected to Table 2 consensus values at byte-verification "
        "2026-05-19.) Reference PDF "
        "references/perez-2001-self-consistent-karplus-3j-chi1.pdf "
        "(page 7086).");
    grp.createAttribute("J_Halpha_Hbeta_coefficients",
        coeff_str(KARPLUS_HA_HB_A, KARPLUS_HA_HB_B, KARPLUS_HA_HB_C) +
        " (Perez, Lohr, Ruterjans & Schmidt 2001 JACS 123:7081, "
        "DOI 10.1021/ja003724j, Table 2 consensus row); dihedral "
        "HA-CA-CB-HB{2 or 3} -- two separate channels J_Halpha_Hbeta2 "
        "and J_Halpha_Hbeta3 for the prochiral methylene pair. "
        "Ile/Val/Thr (single methine Hbeta): same atom in both "
        "channels. ALA / GLY: NaN both channels. Reference PDF "
        "references/perez-2001-self-consistent-karplus-3j-chi1.pdf "
        "(byte-verified 2026-05-19, page 7086).");
    grp.createAttribute("dihedral_convention", std::string(
        "Backbone channels: phi = C(prev)-N-CA-C atomic dihedral "
        "(IUPAC signed atan2, radians); plug into Karplus as "
        "(phi + theta_offset_project). Each backbone channel's offset "
        "uses the project phi sign -- HN-Halpha=+pi/3, HN-Cbeta=-pi/3, "
        "HN-C'=0, Halpha-C'=+pi/3. The 4-atom atomic dihedrals "
        "(H-N-CA-HA, H-N-CA-CB, H-N-CA-C, Halpha-CA-N-C'(prev)) "
        "are documented but NOT plugged in directly -- the "
        "(phi + theta_offset) form matches the Wang-Bax / Vogeli "
        "published fit independent of any ideal-tetrahedral-geometry "
        "assumption about the atomic 4-atom path. Halpha-C' uses "
        "C'(prev) (resolved via Protein::BackbonePredecessor) and is "
        "rotation around N-CA = phi axis per Vuister teaching lecture "
        "sect 6.1 + Vogeli 2007 page 9384. "
        "Chi1 channels: feed the actual atomic 4-atom dihedral "
        "directly (Perez 2001 Table 2 footnote c form); paths are "
        "N-CA-CB-CG for J(N,Cgamma); C-CA-CB-CG for J(C',Cgamma) -- "
        "carbon-Cgamma only (SER OG / CYS SG / THR OG1 NaN by "
        "element gate on chi[0].a[3]); HA-CA-CB-HB{2,3} for "
        "J(Halpha,Hbeta2) / J(Halpha,Hbeta3)."));
    grp.createAttribute("GLY_caveat", std::string(
        "GLY: Halpha uses Residue.HA which is HA2 by Residue.h cache "
        "convention; HA3 is NOT separately measured. Vuister & Bax "
        "1993 fit DID include glycine 3J values, so the published "
        "(A, B, C) absorb the pro-R/pro-S averaging error. Consumers "
        "needing strict pro-R/pro-S resolution should compute directly "
        "from the two Halpha atom indices. J_HN_Cbeta is NaN on GLY "
        "(no Cbeta). J_Halpha_Hbeta2/3 is NaN on GLY (no Cbeta) and "
        "on ALA (methyl Cbeta; the methyl 3J is not the same chemical "
        "observable as the methylene 3J). "
        "SER/CYS/THR: J_N_Cgamma and J_Cprime_Cgamma are NaN "
        "(chi1 terminal is OG/SG/OG1, not Cγ; the consensus Pérez "
        "Table 2 row is for a carbon Cγ -- per-residue heteroatom-"
        "terminal coefficients exist in Pérez Table 2 SER/CYS/THR "
        "rows but are a different chemical observable; not shipped, "
        "see PhysicalConstants.h header)."));
    grp.createAttribute("units",       std::string("Hz"));
    grp.createAttribute("absent_sentinel", std::string("NaN"));
    grp.createAttribute("residue_axis", std::string("protein_residue_index"));
    grp.createAttribute("atom_axis",    std::string("protein_atom_index"));
    grp.createAttribute("source", std::string(
        "positions + Residue backbone-cache (H, N, CA, HA, CB, C) + "
        "C(prev) from Protein::BackbonePredecessor + Residue.chi[0] "
        "(chi1 atom indices from AminoAcidType.chi_angles) + Hbeta "
        "atom-name lookup. No source ConformationResult dependency; "
        "positions always present at tp.Seed time."));
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

    emit_2d_f64("J_HN_Halpha",        j_hn_halpha_);
    emit_2d_f64("J_HN_Halpha_Vogeli", j_hn_halpha_vogeli_);
    emit_2d_f64("J_HN_Cbeta",         j_hn_cbeta_);
    emit_2d_f64("J_HN_Cprime",        j_hn_cprime_);
    emit_2d_f64("J_Halpha_Cprime",    j_halpha_cprime_);
    emit_2d_f64("J_N_Cgamma",         j_n_cgamma_);
    emit_2d_f64("J_Cprime_Cgamma",    j_cprime_cgamma_);
    emit_2d_f64("J_Halpha_Hbeta2",    j_halpha_hbeta2_);
    emit_2d_f64("J_Halpha_Hbeta3",    j_halpha_hbeta3_);

    // ── Static per-residue masks (R,) ────────────────────────────────
    grp.createDataSet("J_HN_Halpha_exists", j_hn_halpha_exists_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("J_HN_Cbeta_exists", j_hn_cbeta_exists_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("J_HN_Cprime_exists", j_hn_cprime_exists_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("J_Halpha_Cprime_exists", j_halpha_cprime_exists_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("J_chi1_exists", j_chi1_exists_)
       .createAttribute("units", std::string("dimensionless"));
    // J_N_Cgamma_exists / J_Cprime_Cgamma_exists are STRICTER than
    // J_chi1_exists -- they additionally require chi[0].a[3] to be
    // Element::C (i.e. exclude SER/CYS/THR whose chi1 terminal is
    // OG/SG/OG1).
    grp.createDataSet("J_N_Cgamma_exists", j_n_cgamma_exists_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("J_Cprime_Cgamma_exists", j_cprime_cgamma_exists_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("J_Halpha_Hbeta_exists", j_halpha_hbeta_exists_)
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
