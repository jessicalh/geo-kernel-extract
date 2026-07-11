#include "LarsenHBondShieldingResult.h"

#include "Bond.h"
#include "ConformationAtom.h"
#include "GeometryChoice.h"
#include "LarsenHBondGeometryCommon.h"
#include "LarsenHBondGrid.h"
#include "LegacyAmberTopology.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "SemanticEnums.h"
#include "SpatialIndexResult.h"
#include "Types.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <set>
#include <string>
#include <typeindex>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace nmr {

namespace {

// Larsen 2015 Δσ_w: amide H atoms with no H-bond pair (geometric)
// get this isotropic offset (NMA + water complex, OPBE/6-31G(d,p)).
constexpr double kWaterTerm_ppm = 2.07;

// Spatial enumeration cutoff. Larsen scan ranges:
//   rOH in [1.5, 3.0] Å for NMA donor, [1.8, 4.0] for ALA donor.
// We pad the cutoff to 4.2 Å (slightly past the ALA donor max so
// boundary candidates are evaluated by the grid's own r-range gate
// rather than silently dropped by the spatial sweep).
constexpr double kSpatialCutoff_A = 4.2;
constexpr double kThetaMinDeg     = 90.0;
constexpr double kThetaMaxDeg     = 180.0;


// Apply a per-class contribution to a target atom's per-class Mat3 field.
void AccumulateContribution(
    ConformationAtom& atom,
    LarsenContribDispatch::Term term,
    const Mat3& contribution_lab) {

    using Term = LarsenContribDispatch::Term;
    switch (term) {
        case Term::Primary_HB:    atom.larsen_hbond_1pHB_tensor  += contribution_lab; break;
        case Term::Secondary_HB:  atom.larsen_hbond_2pHB_tensor  += contribution_lab; break;
        case Term::Primary_HaB:   atom.larsen_hbond_1pHaB_tensor += contribution_lab; break;
        case Term::Secondary_HaB: atom.larsen_hbond_2pHaB_tensor += contribution_lab; break;
        default:                  break;  // RC, w handled elsewhere
    }
}


// Donor-side readout-name → TargetAtom. ALA grid has 6 readouts;
// NMA grid has 5 (no CB).
struct DonorReadout {
    LarsenContribDispatch::TargetAtom target;
    Mat3                              canonical_tensor;
    bool                              present;
};
std::array<DonorReadout, 6> ExtractDonorReadouts(const LarsenHBondRecord& rec) {
    using TA = LarsenContribDispatch::TargetAtom;
    return {{
        {TA::N,  rec.donor_N,  true},
        {TA::CA, rec.donor_CA, true},
        {TA::CB, rec.donor_CB, rec.has_donor_CB},
        {TA::C,  rec.donor_C,  true},
        {TA::HA, rec.donor_HA, true},
        {TA::HN, rec.donor_HN, true},
    }};
}


// Acceptor-side readout-name → TargetAtom. Populated only for
// BackboneCarbonyl acceptor (NMA acceptor archives). The NMA acceptor's
// amide group physically straddles two protein residues j (acceptor's
// own) and j+1, so readouts route to two DIFFERENT target residues:
//
//   • N, HN, CA, HA → residue j+1 (the residue downstream in the chain)
//   • C             → residue j (the acceptor's OWN residue's carbonyl C)
//
// Larsen 2015 Table 2 says Δσ_2°HB applies to {N, Cα, C', Hα, HN}; this
// is the full readout set. The C readout in particular carries the
// largest 2°HB contribution per Larsen Table 1 (~2.1 ppm RMSD on Ub).
struct AcceptorReadout {
    LarsenContribDispatch::TargetAtom target;
    Mat3                              canonical_tensor;
    bool                              present;
    bool                              target_is_acceptor_residue;
};
std::array<AcceptorReadout, 5> ExtractAcceptorReadouts(const LarsenHBondRecord& rec) {
    using TA = LarsenContribDispatch::TargetAtom;
    return {{
        {TA::N,  rec.acceptor_N,  rec.has_acceptor_readouts, false},
        {TA::CA, rec.acceptor_CA, rec.has_acceptor_readouts, false},
        {TA::C,  rec.acceptor_C,  rec.has_acceptor_readouts, true},
        {TA::HA, rec.acceptor_HA, rec.has_acceptor_readouts, false},
        {TA::HN, rec.acceptor_HN, rec.has_acceptor_readouts, false},
    }};
}


// Map LarsenContribDispatch::TargetAtom to the protein-side atom
// indices on a specific Residue. Returns a list because TA::HA can
// resolve to MULTIPLE atoms on glycine (HA2 + HA3 — both prochiral
// α-hydrogens). For non-GLY residues the returned list has one
// element (the single Hα). For other target slots the list is at most
// one element (or empty if absent — e.g. CB for GLY).
//
// Implementation: HA enumeration is substrate-driven via
// AtomSemanticTable::IsAnyAlphaHydrogen() rather than hard-coded
// (matches `feedback_identity_from_chemistry_not_position`). The
// predicate fires on both `BackboneRole::AlphaHydrogen` (non-GLY HA)
// and `(Element::H, Locant::Alpha, BackboneRole::None)` (GLY HA2/HA3
// per Markley convention; see SemanticEnums.h:85-98).
//
// The same enumeration handles the Hα donor case: each GLY HA atom
// independently enumerates as a donor candidate (via the spatial
// sweep over IsAnyAlphaHydrogen) AND each receives readout fan-out
// when its residue is the target of someone else's H-bond.
std::vector<std::size_t> TargetAtomIndices(
        const Protein& protein,
        const Residue& res,
        LarsenContribDispatch::TargetAtom t) {
    using TA = LarsenContribDispatch::TargetAtom;
    std::vector<std::size_t> out;
    out.reserve(2);  // GLY HA is the only multi-atom case (2 entries).
    auto add_if_set = [&](std::size_t idx) {
        if (idx != Residue::NONE) out.push_back(idx);
    };
    switch (t) {
        case TA::N:  add_if_set(res.N);  break;
        case TA::CA: add_if_set(res.CA); break;
        case TA::CB: add_if_set(res.CB); break;
        case TA::C:  add_if_set(res.C);  break;
        case TA::HN: add_if_set(res.H);  break;
        case TA::HA: {
            const auto& topo = protein.LegacyAmber();
            for (std::size_t ai : res.atom_indices) {
                const AtomSemanticTable& sem = topo.SemanticAt(ai);
                if (sem.IsAnyAlphaHydrogen()) out.push_back(ai);
            }
            break;
        }
        default: break;
    }
    return out;
}


}  // namespace


void larsen_hbond_shielding_detail::RecordTable2Contribution(
        LarsenHBondShieldingResult::PairRecord& pair,
        LarsenContribDispatch::TargetAtom target,
        LarsenContribDispatch::Term term,
        const Mat3& contribution_lab) {
    const std::uint32_t target_bit =
        std::uint32_t{1} << static_cast<unsigned>(target);
    const double isotropic = contribution_lab.trace() / 3.0;

    double* term_iso = nullptr;
    std::uint32_t* term_mask = nullptr;
    using Term = LarsenContribDispatch::Term;
    switch (term) {
        case Term::Primary_HB:
            term_iso = &pair.iso_1pHB;
            term_mask = &pair.target_mask_1pHB;
            break;
        case Term::Secondary_HB:
            term_iso = &pair.iso_2pHB;
            term_mask = &pair.target_mask_2pHB;
            break;
        case Term::Primary_HaB:
            term_iso = &pair.iso_1pHaB;
            term_mask = &pair.target_mask_1pHaB;
            break;
        case Term::Secondary_HaB:
            term_iso = &pair.iso_2pHaB;
            term_mask = &pair.target_mask_2pHaB;
            break;
        case Term::RingCurrent:
        case Term::Water:
        case Term::Count:
            return;
    }

    if (!std::isfinite(*term_iso)) *term_iso = 0.0;
    *term_iso += isotropic;
    *term_mask |= target_bit;
}


void larsen_hbond_shielding_detail::RecordDiagnosticCbContribution(
        LarsenHBondShieldingResult::PairRecord& pair,
        const Mat3& contribution_lab) {
    if (!std::isfinite(pair.iso_diagnostic_CB)) {
        pair.iso_diagnostic_CB = 0.0;
    }
    pair.iso_diagnostic_CB += contribution_lab.trace() / 3.0;
    pair.target_mask_diagnostic_CB |=
        std::uint32_t{1} << static_cast<unsigned>(
            LarsenContribDispatch::TargetAtom::CB);
}


// ============================================================================
// Dependencies
// ============================================================================

std::vector<std::type_index> LarsenHBondShieldingResult::Dependencies() const {
    return { std::type_index(typeid(SpatialIndexResult)) };
}


// ============================================================================
// Compute — spatial enumeration over both donor classes (amide H + Hα).
// ============================================================================
//
// Each donor atom is queried against SpatialIndexResult for candidate
// acceptor O atoms within kSpatialCutoff_A. Each candidate is classified
// (BackboneCarbonyl / SidechainCarbonyl / HydroxylOxygen /
// CarboxylateOxygen) via ClassifyAcceptor. Geometry computed; grid
// queried with (donor_class, acceptor_class); tensors rotated to lab
// frame and distributed per Larsen 2015 Table 2 dispatch. Water term
// Δσ_w applies to amide H atoms with zero geometric H-bond candidates
// in the range. Larsen's framework is geometric, not DSSP-energy-based;
// the spatial sweep IS the H-bond finder.

std::unique_ptr<LarsenHBondShieldingResult> LarsenHBondShieldingResult::Compute(
        ProteinConformation& conf,
        const LarsenHBondGrid& grid) {

    OperationLog::Scope scope("LarsenHBondShieldingResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    if (!grid.IsLoaded()) {
        OperationLog::Warn("LarsenHBondShieldingResult::Compute",
            "grid not loaded; returning empty result.");
        return nullptr;
    }
    if (conf.AtomCount() == 0) return nullptr;

    const Protein&    protein = conf.ProteinRef();
    const auto&       spatial = conf.Result<SpatialIndexResult>();
    const auto&       topo    = protein.LegacyAmber();
    const std::size_t n_atoms     = conf.AtomCount();
    const std::size_t n_residues  = protein.ResidueCount();

    auto result_ptr = std::make_unique<LarsenHBondShieldingResult>();
    result_ptr->conf_ = &conf;
    result_ptr->imputed_pair_count_by_atom_.assign(n_atoms, 0);
    result_ptr->sidechain_carbonyl_pair_count_by_atom_.assign(n_atoms, 0);

    GeometryChoiceBuilder choices(conf);
    std::size_t resolution_key = 0;

    // amide_h_geometric_paired: true if ANY candidate acceptor passed
    // the geometric H-bond criterion (θ ∈ [90°, 180°]) for this amide H.
    // Larsen's framework defines an H-bond geometrically; the spatial
    // enumeration we run IS the H-bond finder, with θ as the lone gate.
    // Water term Δσ_w applies to amide Hs with NO such candidate.
    //
    // Critical: the flag is set ONLY when θ passes — spatial proximity
    // alone is not enough. A backbone HN(i) frequently sits within 4.2 Å
    // of O(i-1) but at θ ≈ 20° (the i,i geometry of the chain), which is
    // anti-bonded, not an H-bond. Gating on classification-only would
    // suppress the water term for every amide H regardless of solvent
    // exposure (codex finding F2, 2026-05-12).
    std::vector<bool> amide_h_geometric_paired(n_atoms, false);

    // n_pairs_grid_skipped counts every processed candidate that the
    // grid path could not turn into contributions. Every non-success
    // disposition increments it exactly once:
    //   • MissingFrameAtoms  — classification ok, but a frame anchor
    //                           is None (e.g. chain boundary).
    //   • ThetaOutOfRange    — θ < 90° or > 180°. NOT an H-bond.
    //   • GridMiss           — θ in range, but r outside the grid's
    //                           r-axis bounds. H-bond confirmed; just
    //                           outside Larsen's scan range.
    //   • InvalidFrame       — coordinate/frame geometry is malformed.
    //                           This is NOT a confirmed H-bond.
    //   • CarboxylateSymmetryFiltered — classified farther sibling;
    //                           retained for audit but not evaluated.
    // Mass conservation: processed_candidates == pairs_found
    //                                      + n_pairs_grid_skipped.
    int n_pairs_grid_skipped = 0;

    // backbone_prev_of: canonical Protein::BackbonePredecessor query.
    // Returns the residue index of the predecessor (whose C is bonded
    // to res(d_ri).N) or nullopt. Wrap-correct for cyclic peptides;
    // correct on ACE-capped N-termini (ACE.C is a real atom in the
    // bond graph); correct on antibody insertion-coded structures.
    // Replaces the retired chain_id-based has_same_chain_prev / the
    // intermediate (d_ri-1)-arithmetic gate (2026-05-19).
    auto backbone_prev_of =
        [&](std::size_t d_ri) -> std::optional<std::size_t> {
        return protein.BackbonePredecessor(d_ri);
    };

    // PairResult — disposition for one donor/acceptor candidate.
    //
    // The outer sweep uses these to (a) increment n_pairs_grid_skipped
    // exactly once per non-Success case and (b) decide whether the amide
    // H counts as geometrically paired for the water-term gate (only
    // Success and GridMiss confirm θ ≥ 90°, so only those gate water).
    enum class PairResult : std::uint8_t {
        MissingFrameAtoms = 0,  // Structural; θ never computed.
        ThetaOutOfRange = 1,    // Geometry computed; θ < 90° or > 180°.
        GridMiss = 2,           // θ ok; r outside grid's r-axis bounds.
        Success = 3,            // Tensors accumulated.
        InvalidFrame = 4,       // Non-finite/degenerate geometry or donor frame.
        CarboxylateSymmetryFiltered = 5,  // Diagnostic-only farther COO sibling.
    };
    static_assert(
        static_cast<int>(PairResult::CarboxylateSymmetryFiltered) ==
        static_cast<int>(PairDisposition::CarboxylateSymmetryFiltered));

    using AcceptorTriple = larsen_hbond_geometry::AcceptorTriple;

    // Process one donor → acceptor pair. Resolves donor frame (per
    // donor class), classifies acceptor, computes geometry, queries
    // grid, dispatches per Larsen 2015 Table 2.
    auto process_pair =
        [&](HBondDonorClass donor_class,
            std::size_t donor_H_idx,
            std::size_t donor_anchor_idx,
            std::size_t donor_third_idx,
            std::size_t donor_residue_idx,
            const AcceptorTriple& acc,
            bool selected_carboxylate_representative) -> PairResult {

        PairRecord pair;
        pair.donor_atom_idx = donor_H_idx;
        pair.acceptor_atom_idx = acc.O_idx;
        pair.donor_residue_idx = donor_residue_idx;
        pair.acceptor_residue_idx =
            protein.AtomAt(acc.O_idx).residue_index;
        pair.donor_class = donor_class;
        pair.acceptor_class = acc.class_;
        pair.donor_anchor_atom_idx = donor_anchor_idx;
        pair.donor_third_atom_idx = donor_third_idx;
        pair.acceptor_C_atom_idx = acc.C_idx;
        pair.acceptor_third_atom_idx = acc.third_idx;

        auto finish = [&](PairResult disposition) {
            pair.disposition = static_cast<PairDisposition>(disposition);
            result_ptr->pairs_.push_back(pair);
            if (disposition == PairResult::Success) {
                ++result_ptr->successful_pairs_;
            }
            return disposition;
        };

        // Skip degenerate (donor H == any frame anchor) configurations.
        if (donor_anchor_idx == Residue::NONE ||
            donor_third_idx  == Residue::NONE ||
            acc.C_idx        == Residue::NONE ||
            acc.third_idx    == Residue::NONE) {
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "missing frame anchor",
                [donor_class, ac_class = acc.class_,
                 donor_anchor_idx, donor_third_idx,
                 acc_C = acc.C_idx, acc_thd = acc.third_idx]
                (GeometryChoice& gc) {
                    AddNumber(gc, "donor_class",
                        static_cast<double>(donor_class), "enum");
                    AddNumber(gc, "acceptor_class",
                        static_cast<double>(ac_class), "enum");
                    AddNumber(gc, "donor_anchor_idx",
                        static_cast<double>(donor_anchor_idx), "atom_idx");
                    AddNumber(gc, "donor_third_idx",
                        static_cast<double>(donor_third_idx), "atom_idx");
                    AddNumber(gc, "acceptor_C_idx",
                        static_cast<double>(acc_C), "atom_idx");
                    AddNumber(gc, "acceptor_third_idx",
                        static_cast<double>(acc_thd), "atom_idx");
                    AddNumber(gc, "rejection", 1.0, "missing_frame_anchor");
                });
            return finish(PairResult::MissingFrameAtoms);
        }

        Vec3 donor_H_pos    = conf.PositionAt(donor_H_idx);
        Vec3 donor_anchor   = conf.PositionAt(donor_anchor_idx);
        Vec3 donor_third    = conf.PositionAt(donor_third_idx);
        Vec3 accept_O_pos   = conf.PositionAt(acc.O_idx);
        Vec3 accept_C_pos   = conf.PositionAt(acc.C_idx);
        Vec3 accept_thd_pos = conf.PositionAt(acc.third_idx);

        LarsenHBondGeometry geom = ComputeLarsenHBondGeometry(
            donor_H_pos, accept_O_pos, accept_C_pos, accept_thd_pos);
        pair.r_angstrom = geom.r_angstrom;
        pair.theta_deg = geom.theta_deg;
        pair.rho_deg = geom.rho_deg;

        const LarsenDonorFrame donor_frame = ComputeLarsenDonorFrame(
            donor_H_pos, donor_anchor, donor_third);
        pair.frame_valid = geom.IsFinite() && donor_frame.valid;
        if (!pair.frame_valid) {
            OperationLog::Warn(
                "LarsenHBondShieldingResult::Compute",
                "invalid Larsen pair frame: donor_atom=" +
                    std::to_string(donor_H_idx) +
                    " acceptor_atom=" + std::to_string(acc.O_idx) +
                    " geometry_finite=" +
                    std::to_string(geom.IsFinite() ? 1 : 0) +
                    " donor_frame_failure=" + std::to_string(
                        static_cast<int>(donor_frame.failure)));
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "invalid Larsen pair frame",
                [geom, donor_class, ac_class = acc.class_, donor_frame]
                (GeometryChoice& gc) {
                    AddNumber(gc, "donor_class",
                        static_cast<double>(donor_class), "enum");
                    AddNumber(gc, "acceptor_class",
                        static_cast<double>(ac_class), "enum");
                    AddNumber(gc, "r_angstrom", geom.r_angstrom, "A");
                    AddNumber(gc, "theta_deg", geom.theta_deg, "degrees");
                    AddNumber(gc, "rho_deg", geom.rho_deg, "degrees");
                    AddNumber(gc, "geometry_finite",
                        geom.IsFinite() ? 1.0 : 0.0, "bool");
                    AddNumber(gc, "donor_frame_failure",
                        static_cast<double>(donor_frame.failure), "enum");
                    AddNumber(gc, "rejection", 1.0, "invalid_frame");
                });
            return finish(PairResult::InvalidFrame);
        }

        // Preserve the physical one-oxygen acetate selection, but retain
        // the classified farther sibling as an explicit diagnostic row.
        // No grid query, tensor rotation, count, or water-gate decision is
        // made from this row.
        if (!selected_carboxylate_representative) {
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "carboxylate symmetry sibling filtered",
                [donor_H_idx, acceptor_O_idx = acc.O_idx, geom]
                (GeometryChoice& gc) {
                    AddNumber(gc, "donor_atom",
                        static_cast<double>(donor_H_idx), "atom_idx");
                    AddNumber(gc, "acceptor_O",
                        static_cast<double>(acceptor_O_idx), "atom_idx");
                    AddNumber(gc, "r_angstrom", geom.r_angstrom, "A");
                    AddNumber(gc, "theta_deg", geom.theta_deg, "degrees");
                    AddNumber(gc, "rho_deg", geom.rho_deg, "degrees");
                    AddNumber(gc, "rejection", 1.0,
                        "carboxylate_symmetry_filtered");
                });
            return finish(PairResult::CarboxylateSymmetryFiltered);
        }

        if (geom.theta_deg < kThetaMinDeg || geom.theta_deg > kThetaMaxDeg) {
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "theta out of range",
                [geom, donor_class](GeometryChoice& gc) {
                    AddNumber(gc, "donor_class",
                        static_cast<double>(donor_class), "enum");
                    AddNumber(gc, "theta_deg", geom.theta_deg, "degrees");
                    AddNumber(gc, "rejection", 1.0, "theta_out_of_range");
                });
            return finish(PairResult::ThetaOutOfRange);
        }

        LarsenHBondRecord rec = grid.QueryNearest(
            donor_class, acc.class_, geom);
        if (!rec.IsHit()) {
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "grid query miss",
                [geom, donor_class, ac_class = acc.class_](GeometryChoice& gc) {
                    AddNumber(gc, "donor_class",
                        static_cast<double>(donor_class), "enum");
                    AddNumber(gc, "acceptor_class",
                        static_cast<double>(ac_class), "enum");
                    AddNumber(gc, "r_angstrom", geom.r_angstrom, "A");
                    AddNumber(gc, "theta_deg",  geom.theta_deg, "degrees");
                    AddNumber(gc, "rho_deg",    geom.rho_deg, "degrees");
                    AddNumber(gc, "rejection",  1.0, "grid_miss");
                });
            return finish(PairResult::GridMiss);
        }

        // A successful row has finite, explicit zeroes for terms that do
        // not apply. Each actual tensor write below adds its own trace/3
        // through the production audit helper.
        // Pair geometry is the raw protein geometry. QueryNearest wraps rho
        // for its periodic axis internally; do not replace the audit value
        // with that compatibility-normalized coordinate.
        pair.iso_1pHB = 0.0;
        pair.iso_2pHB = 0.0;
        pair.iso_1pHaB = 0.0;
        pair.iso_2pHaB = 0.0;
        pair.iso_diagnostic_CB = 0.0;
        pair.iso_table2_total = 0.0;
        pair.isotropic_total = 0.0;
        pair.any_corner_imputed = rec.any_corner_imputed;
        pair.imputed_corner_count = rec.imputed_corner_count;

        const Mat3& R_protein = donor_frame.rotation;

        // Pick the Table 2 Term names per donor class (primary +
        // secondary). HB terms apply when donor is amide H; HαB terms
        // when donor is Hα. The substrate of the calculation is the
        // same — only the term labels differ.
        using Term = LarsenContribDispatch::Term;
        using TA   = LarsenContribDispatch::TargetAtom;
        const Term primary_term =
            (donor_class == HBondDonorClass::AmideHydrogen)
                ? Term::Primary_HB : Term::Primary_HaB;
        const Term secondary_term =
            (donor_class == HBondDonorClass::AmideHydrogen)
                ? Term::Secondary_HB : Term::Secondary_HaB;

        const Residue& don_res = protein.ResidueAt(donor_residue_idx);

        auto donor_readouts    = ExtractDonorReadouts(rec);
        auto acceptor_readouts = ExtractAcceptorReadouts(rec);

        std::set<std::size_t> table2_contributors;
        std::set<std::size_t> diagnostic_contributors;

        // Cβ DIAGNOSTIC — Table 2 says Cβ gets NO contribution under
        // any term. We still rotate and emit Cβ so a downstream test
        // can verify the parser → loader → frame-rotation pipeline
        // produces near-zero (the physics statement). Emitted BEFORE
        // the Table 2 dispatch because the dispatch would short-
        // circuit on Cβ otherwise.
        for (const auto& dr_readout : donor_readouts) {
            if (!dr_readout.present)         continue;
            if (dr_readout.target != TA::CB) continue;
            auto targets = TargetAtomIndices(protein, don_res, TA::CB);
            if (targets.empty())             continue;
            Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                dr_readout.canonical_tensor, R_protein);
            for (std::size_t target_ai : targets) {
                conf.MutableAtomAt(target_ai).larsen_hbond_diagnostic_CB +=
                    sigma_lab;
                larsen_hbond_shielding_detail::
                    RecordDiagnosticCbContribution(pair, sigma_lab);
                diagnostic_contributors.insert(target_ai);
            }
        }

        // 1° readouts → donor residue i atoms. TargetAtomIndices may
        // return MULTIPLE atoms for HA on GLY (HA2 + HA3); same tensor
        // is applied to each (each Hα is an independent atom subject
        // to the H-bond geometry Larsen's grid encodes).
        for (const auto& dr_readout : donor_readouts) {
            if (!dr_readout.present)         continue;
            if (dr_readout.target == TA::CB) continue;  // diagnostic above
            if (!LarsenContribDispatch::Applies(dr_readout.target, primary_term))
                continue;
            auto targets = TargetAtomIndices(protein, don_res, dr_readout.target);
            if (targets.empty())             continue;
            Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                dr_readout.canonical_tensor, R_protein);
            for (std::size_t target_ai : targets) {
                AccumulateContribution(
                    conf.MutableAtomAt(target_ai), primary_term, sigma_lab);
                larsen_hbond_shielding_detail::RecordTable2Contribution(
                    pair, dr_readout.target, primary_term, sigma_lab);
                table2_contributors.insert(target_ai);
            }
        }

        // 2° readouts → routed per acceptor_readout.target_is_acceptor_residue:
        //   • N, CA, HA, HN → residue j+1 (the residue downstream of the
        //                     H-bond acceptor in the chain).
        //   • C            → residue j (the acceptor's OWN residue's
        //                     carbonyl C — the C atom whose C=O accepts
        //                     the H-bond). Per Larsen 2015 Table 2 this
        //                     is the largest 2°HB term (~2.1 ppm RMSD).
        //
        // Gated on acc.i_plus_1_residue_idx != SIZE_MAX, which by
        // construction in ClassifyAcceptor is equivalent to "the acceptor
        // is BackboneCarbonyl with a valid chain context." When this is
        // false, the geometry's `third` anchor is also None, so
        // process_pair would have returned MissingFrameAtoms before
        // reaching this point — the gate is belt-and-suspenders.
        //
        // HOMe / Acetate / SidechainCarbonyl never carry acceptor
        // readouts (has_acceptor_readouts=false for those grid files);
        // the per-readout `present` flag filters them.
        if (acc.i_plus_1_residue_idx != SIZE_MAX) {
            const std::size_t acc_j_residue_idx =
                protein.AtomAt(acc.O_idx).residue_index;
            const Residue& res_j        = protein.ResidueAt(acc_j_residue_idx);
            const Residue& res_j_plus_1 =
                protein.ResidueAt(acc.i_plus_1_residue_idx);
            for (const auto& ac_readout : acceptor_readouts) {
                if (!ac_readout.present) continue;
                if (!LarsenContribDispatch::Applies(ac_readout.target, secondary_term))
                    continue;
                const Residue& target_res =
                    ac_readout.target_is_acceptor_residue ? res_j : res_j_plus_1;
                auto targets = TargetAtomIndices(protein, target_res, ac_readout.target);
                if (targets.empty()) continue;
                Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                    ac_readout.canonical_tensor, R_protein);
                for (std::size_t target_ai : targets) {
                    AccumulateContribution(
                        conf.MutableAtomAt(target_ai), secondary_term, sigma_lab);
                    larsen_hbond_shielding_detail::RecordTable2Contribution(
                        pair, ac_readout.target, secondary_term, sigma_lab);
                    table2_contributors.insert(target_ai);
                }
            }
        }

        // Per-pair bookkeeping. n_pairs counts only Table 2 classes;
        // the Cβ diagnostic does NOT inflate n_pairs (per
        // ConformationAtom field doc). any_corner_imputed propagates
        // to every target atom that received any write (Table 2 or
        // diagnostic) so downstream introspection sees the imputation
        // flag wherever a contribution landed.
        for (std::size_t ai : table2_contributors) {
            ConformationAtom& a = conf.MutableAtomAt(ai);
            a.larsen_hbond_n_pairs += 1;
            if (rec.any_corner_imputed) {
                a.larsen_hbond_any_corner_imputed = true;
                ++result_ptr->imputed_pair_count_by_atom_[ai];
            }
            if (acc.class_ == HBondAcceptorClass::SidechainCarbonyl) {
                ++result_ptr->sidechain_carbonyl_pair_count_by_atom_[ai];
            }
        }
        for (std::size_t ai : diagnostic_contributors) {
            if (table2_contributors.count(ai)) continue;
            ConformationAtom& a = conf.MutableAtomAt(ai);
            if (rec.any_corner_imputed) a.larsen_hbond_any_corner_imputed = true;
        }

        choices.Record(CalculatorId::LarsenHBond, resolution_key++,
            "pair included",
            [donor_class, ac_class = acc.class_,
             donor_residue_idx, acc_O_idx = acc.O_idx,
             rec,
             n_table2 = table2_contributors.size(),
             n_diag   = diagnostic_contributors.size()]
            (GeometryChoice& gc) {
                AddNumber(gc, "donor_class",
                    static_cast<double>(donor_class), "enum");
                AddNumber(gc, "acceptor_class",
                    static_cast<double>(ac_class), "enum");
                AddNumber(gc, "donor_residue",
                    static_cast<double>(donor_residue_idx), "index");
                AddNumber(gc, "acceptor_O",
                    static_cast<double>(acc_O_idx), "atom_idx");
                AddNumber(gc, "r_angstrom", rec.r_angstrom, "A");
                AddNumber(gc, "theta_deg",  rec.theta_deg,  "degrees");
                AddNumber(gc, "rho_deg",    rec.rho_deg,    "degrees");
                AddNumber(gc, "n_table2_contributors",
                    static_cast<double>(n_table2), "atoms");
                AddNumber(gc, "n_diagnostic_contributors",
                    static_cast<double>(n_diag), "atoms");
                AddNumber(gc, "any_corner_imputed",
                    rec.any_corner_imputed ? 1.0 : 0.0, "bool");
            });

        pair.iso_table2_total = pair.iso_1pHB + pair.iso_2pHB
                              + pair.iso_1pHaB + pair.iso_2pHaB;
        pair.isotropic_total = pair.iso_table2_total;
        return finish(PairResult::Success);
    };

    // ------------------------------------------------------------------
    // Donor sweep — one pass over all atoms in the protein, dispatching
    // amide H and α-hydrogen donors to the spatial enumeration. The
    // donor frame anchors per donor class:
    //   AmideHydrogen: anchor = res.N, third = prev_res.C
    //                  (prev_res via Protein::BackbonePredecessor —
    //                  bond-graph backbone predecessor of d_ri).
    //   AlphaHydrogen: anchor = res.CA, third = res.N (this residue's
    //                  own N — no preceding-residue dependency).
    // ------------------------------------------------------------------
    for (std::size_t ai = 0; ai < n_atoms; ++ai) {
        const auto& sem = topo.SemanticAt(ai);
        const auto& atom = protein.AtomAt(ai);
        const Residue& res = protein.ResidueAt(atom.residue_index);

        HBondDonorClass donor_class;
        std::size_t donor_anchor_idx = Residue::NONE;
        std::size_t donor_third_idx  = Residue::NONE;

        if (sem.IsBackboneAmideHydrogen()) {
            donor_class = HBondDonorClass::AmideHydrogen;
            donor_anchor_idx = res.N;
            // C'(prev) third atom -- bond-graph predecessor of this
            // residue. No assumption of consecutive numbering: handles
            // cyclic peptide wrap, ACE-capped N-termini, and antibody
            // insertion-coded structures uniformly.
            auto prev_idx = backbone_prev_of(atom.residue_index);
            if (!prev_idx) {
                choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                    "amide donor at chain N-term (no backbone predecessor)",
                    [ri = atom.residue_index](GeometryChoice& gc) {
                        AddNumber(gc, "residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "rejection", 1.0,
                            "no_preceding_C_in_bond_graph");
                    });
            } else {
                const Residue& prev_res = protein.ResidueAt(*prev_idx);
                donor_third_idx = prev_res.C;
            }
        } else if (sem.IsAnyAlphaHydrogen()) {
            donor_class = HBondDonorClass::AlphaHydrogen;
            donor_anchor_idx = res.CA;
            donor_third_idx  = res.N;
        } else {
            continue;  // not a donor candidate
        }
        Vec3 donor_pos = conf.PositionAt(ai);
        auto candidate_atoms =
            spatial.AtomsWithinRadius(donor_pos, kSpatialCutoff_A);

        // For each candidate atom: skip non-O, skip same-residue
        // self-contacts (donor H bonded to its own residue's atoms
        // can't H-bond to them), skip the donor's own anchor and
        // third atom positions. Classify; process.
        bool found_geometric_h_bond = false;  // θ ≥ 90° confirmed
        for (std::size_t o_idx : candidate_atoms) {
            if (o_idx == ai) continue;
            if (protein.AtomAt(o_idx).element != Element::O) continue;
            // Same-residue exclusion: a donor H and an O in the SAME
            // residue are bonded (e.g. an amide H and its own carbonyl
            // O three bonds away — possible for some sidechains).
            // Sequence separation check protects amide-H to backbone-O
            // of the same residue (geometrically infeasible H-bond).
            if (protein.AtomAt(o_idx).residue_index == atom.residue_index) {
                continue;
            }

            auto classified =
                larsen_hbond_geometry::ClassifyAcceptor(protein, o_idx);
            if (!classified.has_value()) continue;

            // Carboxylate symmetry: ASP/GLU/C-term carboxylate Os come in
            // pairs (OD1/OD2, OE1/OE2, O/OXT) bonded to a common C. Both
            // Os iterate as independent spatial candidates. Larsen
            // scanned the acetate grid with one O fixed as the
            // H-bond acceptor — only the closer-to-donor O is the real
            // H-bond. The further sibling remains a diagnostic attempt,
            // but it must not enter the grid/tensor/water-gate physics.
            // Tie-break by atom index when equidistant so exactly one O
            // is selected. The closer-to-donor O is also the one whose
            // `third` (the OTHER O on the same C) points consistently
            // with Larsen's reference geometry.
            bool selected_carboxylate_representative = true;
            if (classified->class_ == HBondAcceptorClass::CarboxylateOxygen) {
                const double d_self =
                    (conf.PositionAt(o_idx) - donor_pos).norm();
                const double d_sib  =
                    (conf.PositionAt(classified->third_idx) - donor_pos).norm();
                if (d_sib < d_self ||
                    (d_sib == d_self && classified->third_idx < o_idx)) {
                    selected_carboxylate_representative = false;
                }
            }

            PairResult r = process_pair(
                donor_class, ai, donor_anchor_idx, donor_third_idx,
                atom.residue_index, *classified,
                selected_carboxylate_representative);

            // Mass-conservation increments + water-term gate. Both
            // Success and GridMiss confirm θ ≥ 90° — only those count
            // toward "this amide H is geometrically H-bonded" so the
            // water term is correctly suppressed only for real H-bonds.
            switch (r) {
                case PairResult::Success:
                    found_geometric_h_bond = true;
                    break;
                case PairResult::GridMiss:
                    found_geometric_h_bond = true;
                    ++n_pairs_grid_skipped;
                    break;
                case PairResult::ThetaOutOfRange:
                case PairResult::MissingFrameAtoms:
                case PairResult::InvalidFrame:
                case PairResult::CarboxylateSymmetryFiltered:
                    ++n_pairs_grid_skipped;
                    break;
            }
        }

        if (donor_class == HBondDonorClass::AmideHydrogen &&
            found_geometric_h_bond) {
            amide_h_geometric_paired[ai] = true;
        }
    }

    // ------------------------------------------------------------------
    // Water term sweep: Δσ_w = 2.07 ppm applies to amide H atoms with
    // ZERO geometric H-bond candidates found in the spatial sweep.
    // Replaces the DSSP-based gate from the earlier draft — Larsen's
    // framework is geometric, and the spatial enumeration we run IS
    // the H-bond finder for this calculator.
    // ------------------------------------------------------------------
    for (std::size_t ri = 0; ri < n_residues; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.H == Residue::NONE) continue;  // PRO etc.
        if (amide_h_geometric_paired[res.H]) continue;
        conf.MutableAtomAt(res.H).larsen_hbond_water_term = kWaterTerm_ppm;
        ++result_ptr->amide_hs_unbound_;
    }

    // ------------------------------------------------------------------
    // Step 4: Per-atom totals + SphericalTensor decomposition for every
    // tensor field (Pattern 11 — both representations always present).
    // ------------------------------------------------------------------
    for (std::size_t ai = 0; ai < n_atoms; ++ai) {
        ConformationAtom& a = conf.MutableAtomAt(ai);
        Mat3 shielding = a.larsen_hbond_1pHB_tensor
                       + a.larsen_hbond_2pHB_tensor
                       + a.larsen_hbond_1pHaB_tensor
                       + a.larsen_hbond_2pHaB_tensor;
        a.larsen_hbond_shielding_tensor    = shielding;
        a.larsen_hbond_shielding_spherical = SphericalTensor::Decompose(shielding);
        a.larsen_hbond_1pHB_spherical      = SphericalTensor::Decompose(a.larsen_hbond_1pHB_tensor);
        a.larsen_hbond_2pHB_spherical      = SphericalTensor::Decompose(a.larsen_hbond_2pHB_tensor);
        a.larsen_hbond_1pHaB_spherical     = SphericalTensor::Decompose(a.larsen_hbond_1pHaB_tensor);
        a.larsen_hbond_2pHaB_spherical     = SphericalTensor::Decompose(a.larsen_hbond_2pHaB_tensor);
        a.larsen_hbond_diagnostic_CB_spherical =
            SphericalTensor::Decompose(a.larsen_hbond_diagnostic_CB);
        // Count: non-zero contribution = at least one class's tensor
        // is non-zero. Use a tiny threshold to avoid FP-noise counting.
        const double kThresh = 1e-9;
        bool has_any = a.larsen_hbond_1pHB_tensor.norm()  > kThresh
                    || a.larsen_hbond_2pHB_tensor.norm()  > kThresh
                    || a.larsen_hbond_1pHaB_tensor.norm() > kThresh
                    || a.larsen_hbond_2pHaB_tensor.norm() > kThresh;
        if (has_any) ++result_ptr->atoms_with_contribution_;
    }

    result_ptr->pairs_grid_skipped_ = n_pairs_grid_skipped;

    OperationLog::Info(LogCalcOther, "LarsenHBondShieldingResult::Compute",
        "pairs_grid_included=" +
            std::to_string(result_ptr->successful_pairs_)
        + " pairs_grid_skipped=" + std::to_string(n_pairs_grid_skipped)
        + " pair_diagnostic_rows=" +
            std::to_string(result_ptr->pairs_.size())
        + " atoms_with_contribution=" + std::to_string(result_ptr->atoms_with_contribution_)
        + " amide_hs_with_water_term=" + std::to_string(result_ptr->amide_hs_unbound_));

    return result_ptr;
}


// ============================================================================
// WriteFeatures — emits SphericalTensor-packed per-class shieldings plus
// the raw per-pair audit tables and the two approximation/imputation counts.
// Packing is T0 (1) + T1 (3) + T2 (5) = 9
// columns per atom, matching the HBondResult convention. The Mat3
// runtime fields stay on ConformationAtom (Pattern 11) but NPY emission
// uses the spherical decomposition — same information, different layout
// downstream consumers key on.
// ============================================================================

int LarsenHBondShieldingResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const std::size_t N = conf.AtomCount();
    if (N == 0) return 0;

    std::vector<double> shielding(N * 9, std::nan(""));
    std::vector<double> sh_1pHB  (N * 9, std::nan(""));
    std::vector<double> sh_2pHB  (N * 9, std::nan(""));
    std::vector<double> sh_1pHaB (N * 9, std::nan(""));
    std::vector<double> sh_2pHaB (N * 9, std::nan(""));
    std::vector<double> sh_CB    (N * 9, std::nan(""));
    std::vector<double> water    (N,     std::nan(""));
    std::vector<std::int32_t> n_pairs(N, 0);
    std::vector<std::int32_t> imputed_pair_count(N, 0);
    std::vector<std::int32_t> sidechain_carbonyl_pair_count(N, 0);
    std::vector<std::int8_t> corner_imputed(N, 0);

    for (std::size_t i = 0; i < N; ++i) {
        const ConformationAtom& a = conf.AtomAt(i);
        a.larsen_hbond_shielding_spherical.PackFull9(    &shielding[i * 9]);
        a.larsen_hbond_1pHB_spherical.PackFull9(         &sh_1pHB  [i * 9]);
        a.larsen_hbond_2pHB_spherical.PackFull9(         &sh_2pHB  [i * 9]);
        a.larsen_hbond_1pHaB_spherical.PackFull9(        &sh_1pHaB [i * 9]);
        a.larsen_hbond_2pHaB_spherical.PackFull9(        &sh_2pHaB [i * 9]);
        a.larsen_hbond_diagnostic_CB_spherical.PackFull9(&sh_CB    [i * 9]);
        water[i]   = a.larsen_hbond_water_term;
        n_pairs[i] = a.larsen_hbond_n_pairs;
        if (i < imputed_pair_count_by_atom_.size()) {
            imputed_pair_count[i] = imputed_pair_count_by_atom_[i];
        }
        if (i < sidechain_carbonyl_pair_count_by_atom_.size()) {
            sidechain_carbonyl_pair_count[i] =
                sidechain_carbonyl_pair_count_by_atom_[i];
        }
        corner_imputed[i] = a.larsen_hbond_any_corner_imputed ? 1 : 0;
    }

    auto npy_index = [](std::size_t index) -> std::int32_t {
        if (index == SIZE_MAX ||
            index > static_cast<std::size_t>(
                std::numeric_limits<std::int32_t>::max())) {
            return -1;
        }
        return static_cast<std::int32_t>(index);
    };

    const std::size_t P = pairs_.size();
    std::vector<std::int32_t> pair_index(P * 16, 0);
    std::vector<double> pair_geometry(P * 6,
        std::numeric_limits<double>::quiet_NaN());
    std::vector<double> pair_isotropic(P * 6,
        std::numeric_limits<double>::quiet_NaN());
    std::vector<double> pair_compat(P * 28,
        std::numeric_limits<double>::quiet_NaN());

    for (std::size_t row = 0; row < P; ++row) {
        const PairRecord& pair = pairs_[row];
        std::int32_t* idx = &pair_index[row * 16];
        idx[0] = npy_index(pair.donor_atom_idx);
        idx[1] = npy_index(pair.acceptor_atom_idx);
        idx[2] = npy_index(pair.donor_residue_idx);
        idx[3] = npy_index(pair.acceptor_residue_idx);
        idx[4] = static_cast<std::int32_t>(pair.donor_class);
        idx[5] = static_cast<std::int32_t>(pair.acceptor_class);
        idx[6] = static_cast<std::int32_t>(pair.disposition);
        idx[7] = npy_index(pair.donor_anchor_atom_idx);
        idx[8] = npy_index(pair.donor_third_atom_idx);
        idx[9] = npy_index(pair.acceptor_C_atom_idx);
        idx[10] = npy_index(pair.acceptor_third_atom_idx);
        idx[11] = static_cast<std::int32_t>(pair.target_mask_1pHB);
        idx[12] = static_cast<std::int32_t>(pair.target_mask_2pHB);
        idx[13] = static_cast<std::int32_t>(pair.target_mask_1pHaB);
        idx[14] = static_cast<std::int32_t>(pair.target_mask_2pHaB);
        idx[15] =
            static_cast<std::int32_t>(pair.target_mask_diagnostic_CB);

        double* geometry = &pair_geometry[row * 6];
        geometry[0] = pair.r_angstrom;
        geometry[1] = pair.theta_deg;
        geometry[2] = pair.rho_deg;
        geometry[3] = pair.any_corner_imputed ? 1.0 : 0.0;
        geometry[4] = static_cast<double>(pair.imputed_corner_count);
        geometry[5] = pair.frame_valid ? 1.0 : 0.0;

        double* isotropic = &pair_isotropic[row * 6];
        isotropic[0] = pair.iso_1pHB;
        isotropic[1] = pair.iso_2pHB;
        isotropic[2] = pair.iso_1pHaB;
        isotropic[3] = pair.iso_2pHaB;
        isotropic[4] = pair.iso_diagnostic_CB;
        isotropic[5] = pair.iso_table2_total;

        double* compat = &pair_compat[row * 28];
        for (int col = 0; col < 16; ++col) {
            compat[col] = static_cast<double>(idx[col]);
        }
        for (int col = 0; col < 6; ++col) {
            compat[16 + col] = geometry[col];
            compat[22 + col] = isotropic[col];
        }
    }

    fs::path dir(output_dir);
    int n_written = 0;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_shielding.npy").string(),
                            shielding.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_1pHB_shielding.npy").string(),
                            sh_1pHB.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_2pHB_shielding.npy").string(),
                            sh_2pHB.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_1pHaB_shielding.npy").string(),
                            sh_1pHaB.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_2pHaB_shielding.npy").string(),
                            sh_2pHaB.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_diagnostic_CB_shielding.npy").string(),
                            sh_CB.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_water_term.npy").string(),
                            water.data(), N, 1); ++n_written;
    NpyWriter::WriteInt32((dir / "larsen_hbond_count.npy").string(),
                          n_pairs.data(), N); ++n_written;
    NpyWriter::WriteInt32(
        (dir / "larsen_imputed_pair_count.npy").string(),
        imputed_pair_count.data(), N); ++n_written;
    NpyWriter::WriteInt32(
        (dir / "larsen_sidechain_carbonyl_pair_count.npy").string(),
        sidechain_carbonyl_pair_count.data(), N); ++n_written;
    NpyWriter::WriteInt8((dir / "larsen_corner_imputed.npy").string(),
                         corner_imputed.data(), N); ++n_written;
    NpyWriter::WriteInt32(
        (dir / "larsen_hbond_pairs_index.npy").string(),
        pair_index.data(), P, 16); ++n_written;
    NpyWriter::WriteFloat64(
        (dir / "larsen_hbond_pairs_geometry.npy").string(),
        pair_geometry.data(), P, 6); ++n_written;
    NpyWriter::WriteFloat64(
        (dir / "larsen_hbond_pairs_isotropic.npy").string(),
        pair_isotropic.data(), P, 6); ++n_written;
    NpyWriter::WriteFloat64(
        (dir / "larsen_hbond_pairs.npy").string(),
        pair_compat.data(), P, 28); ++n_written;
    return n_written;
}


}  // namespace nmr
