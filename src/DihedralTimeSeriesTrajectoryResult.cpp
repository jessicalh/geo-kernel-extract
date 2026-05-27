#include "DihedralTimeSeriesTrajectoryResult.h"

#include "AminoAcidType.h"
#include "Atom.h"
#include "Bond.h"
#include "LegacyAmberTopology.h"
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
#include <highfive/H5PropertyList.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <typeinfo>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr std::uint8_t kRamaUnassigned = 0;
constexpr std::uint8_t kRamaAlphaR     = 1;
constexpr std::uint8_t kRamaBeta       = 2;
constexpr std::uint8_t kRamaAlphaL     = 3;
constexpr std::uint8_t kRamaPPII       = 4;
constexpr std::uint8_t kRamaOther      = 5;

// Dihedral from four positions. Returns signed angle in radians in
// [-π, π] (closed range — atan2 can return both -π and +π exactly).
// Same atan2(y, x) formulation as PlanarGeometryResult and
// ChiRotamerSelectionTrajectoryResult, but with strict NaN guards at every
// degeneracy site — the other two sites differ in degenerate behaviour
// (PlanarGeometry NaN-propagates implicitly via zero-norm normalize();
// ChiRotamer returns 0.0). NaN is the honest signal for "indeterminate
// angle" and lets consumers distinguish from a real 0 rad measurement.
// Per PATTERNS.md the utility-namespace anti-pattern blocks extracting a
// shared helper; equation in comment per PATTERNS Lesson 10 ("Equations in
// comments, Eigen in code"). The drift between the three sites is tracked
// as a follow-up audit per feedback_audit_the_formula_family.
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
    const double x = n1.dot(n2);
    const double y = m1.dot(n2);
    return std::atan2(y, x);
}

// Wrap an angle to [-π, π]. Uses std::remainder which is single-call,
// bounded, and correct for arbitrary inputs — vs the historical while-
// loop pattern which spins forever on pathologically-large inputs.
double WrapPi(double a) {
    if (!std::isfinite(a)) return a;
    return std::remainder(a, 2.0 * M_PI);
}

// Ramachandran-region binning, Lovell-Richardson 2003-aligned grid
// (science-review HIGH 1-4, 2026-05-19). Returns the bin code from
// {kRamaUnassigned, kRamaAlphaR, kRamaBeta, kRamaAlphaL, kRamaPPII,
// kRamaOther}. Inputs in radians; NaN yields kRamaUnassigned.
//
// Boundaries (degrees, inclusive both ends):
//   αR  : phi ∈ [-180,  -30], psi ∈ [-90,  30]   widened from earlier
//                                                 [-100,-30] × [-65,-15]
//                                                 to catch upper-psi
//                                                 helix edge per Lovell
//                                                 2003 favored region.
//   β   : phi ∈ [-180,  -45], psi ∈ [60, 180]    psi lower edge extended
//                                ∪ [-180,-150]   from 90 to 60 to catch
//                                                 the parallel-β / β-
//                                                 PPII overlap zone.
//   αL  : phi ∈ [  30,  100], psi ∈ [-10,  80]   upper edge extended
//                                                 from 50 to 80 per the
//                                                 Lovell αL favored
//                                                 region.
//   PPII: phi ∈ [-75,   -50], psi ∈ [140, 165]   narrowed from earlier
//                                                 broad box to the
//                                                 Berkholz/Adzhubei
//                                                 tight PPII cone, so
//                                                 antiparallel β residues
//                                                 at (-60, +145) no
//                                                 longer get mislabelled
//                                                 PPII.
// Anything else with finite phi+psi → other.
//
// PPII still resolved before β (and αL before PPII, αR before αL).
// Boundary points land in the FIRST matching bin.
//
// References:
//   Lovell, S.C., Davis, I.W., et al. (2003). Structure validation by
//     Calpha geometry: phi, psi and Cbeta deviation. Proteins 50:437.
//   Berkholz, D.S., et al. (2010). Conformation dependence of backbone
//     geometry in proteins. Structure 18:1257.
//   Adzhubei, A.A., et al. (2013). Polyproline-II helix in proteins:
//     structure and function. J. Mol. Biol. 425:2100.
//
// Gly + Pro + pre-Pro have their own static masks (`is_glycine`,
// `is_proline`, `is_pre_proline`); downstream re-bins with type-aware
// Rama maps (or Lovell penultimate-rotamer-library variants) as needed.
std::uint8_t RamachandranBin(double phi_rad, double psi_rad) {
    if (!std::isfinite(phi_rad) || !std::isfinite(psi_rad))
        return kRamaUnassigned;

    const double deg_per_rad = 180.0 / M_PI;
    const double phi = phi_rad * deg_per_rad;
    const double psi = psi_rad * deg_per_rad;

    if (phi >= -180.0 && phi <= -30.0 &&
        psi >=  -90.0 && psi <=  30.0)
        return kRamaAlphaR;

    if (phi >=  30.0 && phi <= 100.0 &&
        psi >= -10.0 && psi <=  80.0)
        return kRamaAlphaL;

    if (phi >= -75.0 && phi <= -50.0 &&
        psi >= 140.0 && psi <= 165.0)
        return kRamaPPII;

    if (phi >= -180.0 && phi <= -45.0 &&
        ((psi >= 60.0 && psi <= 180.0) || (psi >= -180.0 && psi <= -150.0)))
        return kRamaBeta;

    return kRamaOther;
}

// Backbone connectivity is delegated to Protein::BackboneConnected and its
// predecessor/successor helpers; this file's residue-adjacency walks route
// through that bond-graph discipline.

}  // anonymous namespace


std::unique_ptr<DihedralTimeSeriesTrajectoryResult>
DihedralTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<DihedralTimeSeriesTrajectoryResult>();

    const Protein& protein = tp.ProteinRef();
    const std::size_t R = protein.ResidueCount();
    const std::size_t N = tp.AtomCount();

    r->phi_.assign(R, {});
    r->psi_.assign(R, {});
    r->omega_.assign(R, {});
    r->omega_deviation_.assign(R, {});
    r->chi_.assign(R, {});
    r->rama_region_.assign(R, {});

    // ── Static per-residue masks ─────────────────────────────────────
    r->chi_exists_.assign(R * 4, 0);
    r->omega_is_xpro_.assign(R, 0);
    r->is_glycine_.assign(R, 0);
    r->is_proline_.assign(R, 0);
    r->is_pre_proline_.assign(R, 0);
    r->residue_terminal_state_.assign(R, 0);
    r->chain_id_per_residue_.assign(R, std::string{});

    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        for (int k = 0; k < 4; ++k) {
            if (res.chi[k].Valid())
                r->chi_exists_[ri * 4 + k] = 1u;
        }
        if (res.type == AminoAcid::GLY) r->is_glycine_[ri] = 1u;
        if (res.type == AminoAcid::PRO) r->is_proline_[ri] = 1u;
        r->residue_terminal_state_[ri] =
            static_cast<std::uint8_t>(res.terminal_state);
        r->chain_id_per_residue_[ri] = res.chain_id;
    }

    // is_pre_proline + omega_is_xpro: both set on residue i if i's
    // backbone successor (per the canonical Protein::BackboneSuccessor
    // bond-graph query) is a Pro. Wrap-correct: in a cyclic peptide
    // where the successor of res(N-1) is res(0), the mask follows.
    // PlanarGeometryResult sets omega_is_xpro_ only on the i row (the
    // row carrying that peptide bond's omega); we follow that convention.
    for (std::size_t ri = 0; ri < R; ++ri) {
        auto next_idx = protein.BackboneSuccessor(ri);
        if (!next_idx) continue;
        const Residue& res_next = protein.ResidueAt(*next_idx);
        if (res_next.type == AminoAcid::PRO) {
            r->omega_is_xpro_[ri]  = 1u;
            r->is_pre_proline_[ri] = 1u;
        }
    }

    // ── Per-atom lookup ──────────────────────────────────────────────
    r->residue_index_per_atom_.assign(N, -1);
    for (std::size_t ai = 0; ai < N; ++ai) {
        r->residue_index_per_atom_[ai] =
            static_cast<std::int32_t>(protein.AtomAt(ai).residue_index);
    }

    return r;
}


// ── Compute ──────────────────────────────────────────────────────────
//
// Per residue per frame: compute phi, psi, omega, chi[k], rama_region.
// All in radians, NaN where the dihedral is undefined (terminal residue,
// chain break, missing backbone-cache atom, or chi[k] not defined for
// the AA). The per-frame growth buffers are flattened at Finalize.

void DihedralTimeSeriesTrajectoryResult::Compute(
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

        // Resolve backbone predecessor / successor once per residue
        // via the canonical bond-graph queries. Wrap-correct for cyclic
        // peptides; correct on ACE/NME caps; correct on antibody
        // insertion-coded structures.
        const auto prev_idx_opt = protein.BackbonePredecessor(ri);
        const auto next_idx_opt = protein.BackboneSuccessor(ri);

        // ── phi: C(i−1)-N(i)-Cα(i)-C(i) ──────────────────────────────
        // Predecessor query guarantees res.N != NONE and res(prev).C
        // exists; still need res.CA / res.C non-NONE for the dihedral.
        double phi_val = kNaN;
        if (prev_idx_opt &&
            res.CA != Residue::NONE && res.C != Residue::NONE) {
            const Residue& res_prev = protein.ResidueAt(*prev_idx_opt);
            phi_val = Dihedral(
                conf.PositionAt(res_prev.C),
                conf.PositionAt(res.N),
                conf.PositionAt(res.CA),
                conf.PositionAt(res.C));
        }
        phi_[ri].push_back(phi_val);

        // ── psi: N(i)-Cα(i)-C(i)-N(i+1) ──────────────────────────────
        // Successor query guarantees res.C != NONE and res(next).N
        // exists; still need res.N / res.CA non-NONE.
        double psi_val = kNaN;
        if (next_idx_opt &&
            res.N != Residue::NONE && res.CA != Residue::NONE) {
            const Residue& res_next = protein.ResidueAt(*next_idx_opt);
            psi_val = Dihedral(
                conf.PositionAt(res.N),
                conf.PositionAt(res.CA),
                conf.PositionAt(res.C),
                conf.PositionAt(res_next.N));
        }
        psi_[ri].push_back(psi_val);

        // ── omega: Cα(i)-C(i)-N(i+1)-Cα(i+1) ─────────────────────────
        // Successor query guarantees res.C != NONE and res(next).N
        // exists; still need res.CA / res_next.CA non-NONE.
        double omega_val = kNaN;
        if (next_idx_opt && res.CA != Residue::NONE) {
            const Residue& res_next = protein.ResidueAt(*next_idx_opt);
            if (res_next.CA != Residue::NONE) {
                omega_val = Dihedral(
                    conf.PositionAt(res.CA),
                    conf.PositionAt(res.C),
                    conf.PositionAt(res_next.N),
                    conf.PositionAt(res_next.CA));
            }
        }
        omega_[ri].push_back(omega_val);

        // omega_deviation = WrapPi(omega − π) ∈ [-π, π]. Emitted for
        // every well-defined peptide bond INCLUDING X→Pro (cis/trans
        // isomerism is real signal, not a deviation — use the
        // omega_is_xpro static mask to flag those rows). Matches the
        // PlanarGeometryResult production implementation.
        const double omega_dev = std::isfinite(omega_val)
            ? WrapPi(omega_val - M_PI) : kNaN;
        omega_deviation_[ri].push_back(omega_dev);

        // ── chi[k] from Residue.chi[k] pre-cached atom indices ───────
        std::array<double, 4> chi_row{kNaN, kNaN, kNaN, kNaN};
        for (int k = 0; k < 4; ++k) {
            if (!res.chi[k].Valid()) continue;
            chi_row[k] = Dihedral(
                conf.PositionAt(res.chi[k].a[0]),
                conf.PositionAt(res.chi[k].a[1]),
                conf.PositionAt(res.chi[k].a[2]),
                conf.PositionAt(res.chi[k].a[3]));
        }
        chi_[ri].push_back(chi_row);

        rama_region_[ri].push_back(RamachandranBin(phi_val, psi_val));
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(1u);  // positions always present
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────────
//
// Idempotent in state: sets finalized_ true and leaves the per-residue
// growth buffers populated for WriteH5Group.

void DihedralTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                  Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "DihedralTimeSeriesTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(phi_.size()) + " residues.");
}


// ── WriteH5Group ─────────────────────────────────────────────────────
//
// Flat-array emission for each per-residue dataset. Residue-major
// layout matches the natural reader: res-0-frame-0..T-1, res-1-frame-0..,
// ..., res-R-1-frame-T-1.

void DihedralTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t R = phi_.size();
    const std::size_t T = n_frames_;
    const std::size_t N = tp.AtomCount();

    auto grp = file.createGroup("/trajectory/dihedral_time_series");

    grp.createAttribute("result_name",  Name());
    grp.createAttribute("n_residues",   R);
    grp.createAttribute("n_atoms",      N);
    grp.createAttribute("n_frames",     T);
    grp.createAttribute("finalized",    finalized_);

    grp.createAttribute("angle_units",  std::string("radians"));
    grp.createAttribute("periodicity",  std::string("2pi"));
    grp.createAttribute("value_range",  std::string(
        "[-pi, pi] (atan2 closed-range; both -pi and +pi can occur "
        "exactly). Only omega_deviation is explicitly wrapped post-hoc. "
        "Consumers comparing across frames must handle the +-pi "
        "discontinuity (use circular differences)."));
    grp.createAttribute("angle_convention", std::string(
        "IUPAC signed dihedral atan2(y,x). "
        "phi   = C(i-1)-N(i)-CA(i)-C(i); "
        "psi   = N(i)-CA(i)-C(i)-N(i+1); "
        "omega = CA(i)-C(i)-N(i+1)-CA(i+1); "
        "chi_k from AminoAcidType.chi_angles (Residue.chi[k] pre-cached "
        "atom indices, IUPAC sidechain order). "
        "NOTE: DsspResult.Phi/Psi forward libdssp/libcifpp values which "
        "use the NEGATED IUPAC sign convention "
        "(phi_DSSP = -phi_IUPAC). This TR and PlanarGeometryResult use "
        "IUPAC directly; downstream code comparing DsspResult to this "
        "TR must negate one side. Verified by "
        "test_dihedral_time_series.cpp CrossResultConsistency*."));
    grp.createAttribute("chain_break_policy", std::string(
        "NaN at any dihedral spanning a non-bonded residue boundary. "
        "Connectivity queried via Protein::BackbonePredecessor / "
        "Protein::BackboneSuccessor, which walk the cifpp bond graph "
        "and filter on Bond::IsPeptideBond() (BondOrder::Peptide / "
        "BondCategory::PeptideCN, assigned at loader time by "
        "CovalentTopology.cpp when atoms are res_a.C / res_b.N of "
        "different residues). Geometry-native substrate; cleanly "
        "handles ACE/NME-capped termini (cap C/N participate in the "
        "bond graph and get the PeptideCN tag), antibody insertion-"
        "coded structures, engineered chimeras with non-monotonic "
        "numbering, and residue numbering gaps with intact peptide "
        "bonds. Cyclic peptides: API returns the wrap edge "
        "(Successor(N-1) = 0) when the loader tagged the head-to-tail "
        "bond as PeptideCN; cyclic case not numerically validated on "
        "the current linear-only test fleet -- see Protein.h."));
    grp.createAttribute("omega_deviation_policy", std::string(
        "WrapPi(omega - pi) in [-pi, pi]. Emitted for EVERY well-defined "
        "peptide bond INCLUDING X->Pro bonds (cis/trans isomerism at "
        "X-Pro is real signal, not a deviation -- use the omega_is_xpro "
        "static mask to flag those rows for consumer-side interpretation). "
        "Matches the PlanarGeometryResult.cpp:302-303 production impl. "
        "PG's own header doc claims NaN-fill at X->Pro but the impl "
        "emits the actual value; this TR aligns with PG impl. "
        "Recommended consumer pattern for aggregate stats (X-Pro cis "
        "lobe is bimodal and would otherwise leak into 'normal' tails): "
        "  omega_dev_normal = omega_deviation[~omega_is_xpro[:, None]] "
        "  omega_dev_xpro   = omega_deviation[omega_is_xpro[:, None]] "
        "  # analyse the X-Pro distribution separately (~5% cis vs trans). "
        "WrapPi uses std::remainder for bit-identical agreement with "
        "PlanarGeometryResult.cpp::WrapPi (math-review MED-2 fix, "
        "2026-05-19)."));
    grp.createAttribute("rama_region_legend", std::string(
        "0=unassigned, 1=alphaR, 2=beta, 3=alphaL, 4=PPII, 5=other"));
    grp.createAttribute("rama_region_boundaries", std::string(
        "alphaR: phi[-180,-30], psi[-90,30]; "
        "beta: phi[-180,-45], psi[60,180]U[-180,-150]; "
        "alphaL: phi[30,100], psi[-10,80]; "
        "PPII: phi[-75,-50], psi[140,165] (tight Berkholz/Adzhubei cone); "
        "boundaries in degrees, inclusive both ends. "
        "Resolution order: alphaR -> alphaL -> PPII -> beta -> other "
        "(first match wins). PPII narrowed (2026-05-19) so antiparallel "
        "beta residues near (-60,+145) are no longer mislabeled PPII. "
        "References: Lovell et al. 2003 Proteins 50:437; Berkholz et al. "
        "2010 Structure 18:1257; Adzhubei et al. 2013 J. Mol. Biol. "
        "425:2100. Downstream re-binners can use raw phi+psi (also "
        "emitted) for Lovell penultimate-rotamer-library variants."));
    grp.createAttribute("chi_symmetry_caveats", std::string(
        "chi mod-pi (or near-mod-pi) symmetries that consumers must apply "
        "themselves -- raw chi here is the IUPAC signed value: "
        "PHE chi2 (CD1<->CD2 ring flip), TYR chi2 (CD1<->CD2 ring flip), "
        "ASP chi2 (OD1<->OD2 carboxylate flip), GLU chi3 (OE1<->OE2 "
        "carboxylate flip). Near-mod-pi at equilibrium: ARG chi-terminal "
        "(guanidinium NH1<->NH2). Mod-(2pi/3): LYS chi-terminal (NH3+ "
        "3-fold). NOT SYMMETRIC: TRP chi2 (CD1 / CD2 chemically distinct "
        "across 5/6 ring junction), HIS chi2 (ND1 / CD2 chemically "
        "distinct). HIS chi2 atom convention per Markley et al. 1998 "
        "(IUPAC nomenclature, Pure Appl. Chem. 70:117): CA-CB-CG-ND1; "
        "older sources sometimes use CA-CB-CG-CD2 which differs by "
        "~120 degrees (HIS asymmetric). Rotamer counters that need "
        "modular reduction must apply it residue-by-residue."));
    grp.createAttribute("residue_terminal_state_legend", std::string(
        "0=internal, 1=n_terminus, 2=c_terminus, 3=n_and_c_terminus, 4=unknown. "
        "NOTE: residue_terminal_state is loader-assigned CHAIN-ORDER "
        "METADATA, not a validity signal for dihedrals. For cyclic "
        "peptides the loader may still report NTerminus/CTerminus on "
        "the wrap residues even though phi/psi/omega are finite (the "
        "bond graph carries the wrap edge). Consumers must use "
        "isfinite(phi/psi/omega), not terminal_state, to test whether "
        "a dihedral was actually computed."));
    grp.createAttribute("residue_axis", std::string("protein_residue_index"));
    grp.createAttribute("atom_axis",    std::string("protein_atom_index"));
    grp.createAttribute("source",       std::string(
        "positions + Residue.chi[k] (AminoAcidType.chi_angles) + Protein "
        "chain structure; no source ConformationResult dependency."));
    grp.createAttribute("chunking_policy", std::string(
        "Per-residue datasets chunked as {R, min(T, 64)} -- frame slabs "
        "(`phi[:, t]`) fit in one chunk read for the movie-target viewer."));
    grp.createAttribute("source_attached_policy", std::string(
        "always_attached -- positions are always present at tp.Seed time, "
        "so this TR has no conditional source. source_attached_per_frame "
        "is emitted as all-1 for SDK uniformity with conditionally-"
        "attached-source TRs (see OBJECT_MODEL.md 'Conditional-attach TR "
        "discipline, 2026-05-15')."));

    // ── Per-frame (T,) ───────────────────────────────────────────────
    grp.createDataSet("frame_indices", frame_indices_)
       .createAttribute("units", std::string("frame_index"));
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));
    auto ds_attached = grp.createDataSet(
        "source_attached_per_frame", source_attached_per_frame_);
    ds_attached.createAttribute("units", std::string("dimensionless"));
    ds_attached.createAttribute("note", std::string(
        "Trivially always 1 for this TR. Positions are always present at "
        "tp.Seed time; this TR has no conditional source ConformationResult. "
        "Dataset emitted to match the Conditional-attach TR discipline "
        "(OBJECT_MODEL.md 'Conditional-attach TR discipline, 2026-05-15') "
        "for SDK uniformity -- consumers reading either flavour of TR "
        "see the same dataset shape."));

    // Movie-target chunking: chunk shape {R, min(T, 64)} so per-frame
    // slabs (`phi[:, t]`) are at most one chunk read. Larger chunks
    // amortise the HDF5 chunk-cache overhead — 64 frames × R × 8 bytes
    // is ~50 KB at R=100, in the recommended 16 KB-1 MB band. Smaller
    // chunks would page-thrash; larger would over-read for single-frame
    // movie playback.
    const std::size_t frame_chunk = std::min<std::size_t>(T, 64);

    auto make_chunk_props_2d = [&]() {
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking(std::vector<hsize_t>{
            static_cast<hsize_t>(R),
            static_cast<hsize_t>(std::max<std::size_t>(frame_chunk, 1u))
        }));
        return props;
    };

    // ── Helpers: emit a residue-major flat 2D dataset (R, T). ────────
    auto emit_2d_f64 = [&](const std::string& name,
                            const std::vector<std::vector<double>>& src,
                            const std::string& units) {
        std::vector<double> flat(R * T);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto& row = src[ri];
            for (std::size_t f = 0; f < T; ++f) {
                flat[ri * T + f] = (f < row.size()) ? row[f] : kNaN;
            }
        }
        const std::vector<std::size_t> dims = {R, T};
        HighFive::DataSpace space(dims);
        auto props = make_chunk_props_2d();
        auto ds = grp.createDataSet<double>(name, space, props);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };

    auto emit_2d_u8 = [&](const std::string& name,
                           const std::vector<std::vector<std::uint8_t>>& src,
                           const std::string& units) {
        std::vector<std::uint8_t> flat(R * T, 0);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto& row = src[ri];
            for (std::size_t f = 0; f < T && f < row.size(); ++f) {
                flat[ri * T + f] = row[f];
            }
        }
        const std::vector<std::size_t> dims = {R, T};
        HighFive::DataSpace space(dims);
        auto props = make_chunk_props_2d();
        auto ds = grp.createDataSet<std::uint8_t>(name, space, props);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };

    // ── Per-residue per-frame (R, T) ─────────────────────────────────
    emit_2d_f64("phi",             phi_,             "radians");
    emit_2d_f64("psi",             psi_,             "radians");
    emit_2d_f64("omega",           omega_,           "radians");
    emit_2d_f64("omega_deviation", omega_deviation_, "radians");
    emit_2d_u8("rama_region",      rama_region_,     "category");

    // ── chi (R, T, 4) ────────────────────────────────────────────────
    {
        std::vector<double> flat(R * T * 4, kNaN);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto& row = chi_[ri];
            for (std::size_t f = 0; f < T && f < row.size(); ++f) {
                for (int k = 0; k < 4; ++k) {
                    flat[(ri * T + f) * 4 + k] = row[f][k];
                }
            }
        }
        const std::vector<std::size_t> dims = {R, T, std::size_t(4)};
        HighFive::DataSpace space(dims);
        HighFive::DataSetCreateProps chi_props;
        chi_props.add(HighFive::Chunking(std::vector<hsize_t>{
            static_cast<hsize_t>(R),
            static_cast<hsize_t>(std::max<std::size_t>(frame_chunk, 1u)),
            hsize_t(4)
        }));
        auto ds = grp.createDataSet<double>("chi", space, chi_props);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("radians"));
        ds.createAttribute("axis_3", std::string("chi_index_0_to_3"));
    }

    // ── Per-residue static (R,) and (R, 4) ───────────────────────────
    {
        const std::vector<std::size_t> dims4 = {R, std::size_t(4)};
        HighFive::DataSpace s4(dims4);
        auto ds = grp.createDataSet<std::uint8_t>("chi_exists", s4);
        ds.write_raw(chi_exists_.data());
        ds.createAttribute("units", std::string("dimensionless"));
        ds.createAttribute("note", std::string(
            "1 = chi[k] is STRUCTURALLY CACHEABLE at this residue (all "
            "4 atom indices in Residue.chi[k] are non-NONE); 0 = chi[k] "
            "not defined by the residue's AminoAcidType.chi_angles OR "
            "the residue is structurally broken (e.g. ARG missing CZ). "
            "Note: chi_exists==1 does NOT guarantee a finite chi value "
            "at runtime -- per-frame geometry can be degenerate. Use "
            "`isfinite(chi[ri, t, k])` for per-frame validity."));
    }
    auto put_u8 = [&](const std::string& name,
                       const std::vector<std::uint8_t>& v,
                       const std::string& note) {
        auto ds = grp.createDataSet(name, v);
        ds.createAttribute("units", std::string("dimensionless"));
        ds.createAttribute("note", note);
    };
    put_u8("omega_is_xpro", omega_is_xpro_,
        "1 if residue i+1 is PRO; cis/trans isomerism at X-Pro is real "
        "conformational signal (X-Pro cis is the canonical case). Flag "
        "is set on residue i (whose peptide bond INTO Pro is the one to "
        "interpret). omega_deviation at this row is still the actual "
        "WrapPi(omega - pi) value -- consumer must apply the X-Pro "
        "interpretation themselves.");
    put_u8("is_glycine", is_glycine_,
        "1 if residue is GLY; Rama allowed region is much wider for Gly.");
    put_u8("is_proline", is_proline_,
        "1 if residue is PRO; phi constrained to ~[-90, -30] by the ring.");
    put_u8("is_pre_proline", is_pre_proline_,
        "1 if residue i+1 is PRO; flag is set on residue i (whose psi is "
        "constrained by the next-residue Pro side chain). i has its own "
        "constrained Rama region (separate plot from internal residues).");
    put_u8("residue_terminal_state", residue_terminal_state_,
        "0=internal, 1=n_terminus, 2=c_terminus, 3=n_and_c_terminus, 4=unknown");

    // chain_id_per_residue (R,) variable-length strings
    grp.createDataSet("chain_id_per_residue", chain_id_per_residue_)
       .createAttribute("units", std::string("chain_label"));

    // residue_index_per_atom (N,) int32 atom→residue lookup
    grp.createDataSet("residue_index_per_atom", residue_index_per_atom_)
       .createAttribute("units", std::string("residue_index"));

    OperationLog::Info(LogCalcOther,
        "DihedralTimeSeriesTrajectoryResult::WriteH5Group",
        "wrote /trajectory/dihedral_time_series with " +
        std::to_string(R) + " residues x " + std::to_string(T) + " frames");
}


}  // namespace nmr
