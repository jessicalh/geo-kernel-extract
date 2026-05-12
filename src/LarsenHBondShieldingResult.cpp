#include "LarsenHBondShieldingResult.h"

#include "ConformationAtom.h"
#include "DsspResult.h"
#include "GeometryChoice.h"
#include "LarsenHBondGrid.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "SpatialIndexResult.h"
#include "Types.h"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <set>
#include <string>
#include <typeindex>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace nmr {

namespace {

// Larsen 2015 Δσ_w: amide H atoms that form NO H-bond pair get this
// isotropic offset (NMA + water complex, OPBE/6-31G(d,p), CPCM).
constexpr double kWaterTerm_ppm = 2.07;

// Geometric criteria (Larsen 2015 §H-bond scans):
//   r_OH in [1.5, 3.0] Å for NMA donor, [1.8, 4.0] for ALA donor.
//   θ in [90°, 180°].
// We let the grid itself reject out-of-range; the only thing we gate
// here at calculator side is the max r used to limit pair enumeration
// (DSSP already filters to plausible H-bond range, so DSSP-resolved
// pairs land within the grid bounds with rare exceptions).

// Map AtomSemanticTable BackboneRole to LarsenContribDispatch::TargetAtom.
// Returns nullopt for atoms outside the 6 ProCS-relevant target roles.
struct TargetSlot {
    LarsenContribDispatch::TargetAtom target;
    bool valid;
};
TargetSlot ResolveTargetSlot(const AtomSemanticTable& sem) {
    using TA = LarsenContribDispatch::TargetAtom;
    switch (sem.backbone_role) {
        case BackboneRole::Nitrogen:        return {TA::N,  true};
        case BackboneRole::AlphaCarbon:     return {TA::CA, true};
        case BackboneRole::CarbonylCarbon:  return {TA::C,  true};
        case BackboneRole::AlphaHydrogen:   return {TA::HA, true};
        case BackboneRole::AmideHydrogen:   return {TA::HN, true};
        default:                            break;
    }
    // CB: not a backbone role; check Element+Locant.
    if (sem.element == Element::C && sem.locant == Locant::Beta
        && sem.backbone_role == BackboneRole::None) {
        return {TA::CB, true};
    }
    // GLY HA2/HA3: Locant::Alpha + BackboneRole::None + Element::H.
    if (sem.element == Element::H && sem.locant == Locant::Alpha
        && sem.backbone_role == BackboneRole::None) {
        return {TA::HA, true};
    }
    return {TA::N, false};  // invalid; caller checks .valid
}


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


// Acceptor-side readout-name → TargetAtom (in residue j+1's coord space).
// Populated only for BackboneCarbonyl acceptor (NMA acceptor archives).
struct AcceptorReadout {
    LarsenContribDispatch::TargetAtom target;
    Mat3                              canonical_tensor;
    bool                              present;
};
std::array<AcceptorReadout, 3> ExtractAcceptorReadouts(const LarsenHBondRecord& rec) {
    using TA = LarsenContribDispatch::TargetAtom;
    return {{
        {TA::N,  rec.acceptor_N,  rec.has_acceptor_readouts},
        {TA::HN, rec.acceptor_HN, rec.has_acceptor_readouts},
        {TA::HA, rec.acceptor_HA, rec.has_acceptor_readouts},
    }};
}


// Map LarsenContribDispatch::TargetAtom to the protein-side atom
// index on a specific Residue. Returns Residue::NONE if not present.
std::size_t TargetAtomIndex(const Residue& res, LarsenContribDispatch::TargetAtom t) {
    using TA = LarsenContribDispatch::TargetAtom;
    switch (t) {
        case TA::N:  return res.N;
        case TA::CA: return res.CA;
        case TA::CB: return res.CB;
        case TA::C:  return res.C;
        case TA::HA: return res.HA;
        case TA::HN: return res.H;
        default:     return Residue::NONE;
    }
}


}  // namespace


// ============================================================================
// Dependencies
// ============================================================================

std::vector<std::type_index> LarsenHBondShieldingResult::Dependencies() const {
    return {
        std::type_index(typeid(DsspResult)),
        std::type_index(typeid(SpatialIndexResult)),
    };
}


// ============================================================================
// Compute
// ============================================================================

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
    if (conf.AtomCount() == 0) {
        return nullptr;
    }

    const Protein&    protein = conf.ProteinRef();
    const auto&       dssp    = conf.Result<DsspResult>();
    // SpatialIndexResult is a declared dep for Phase 2 (Hα donors via
    // spatial search); Phase 1 uses DSSP only.
    (void)conf.Result<SpatialIndexResult>();

    const std::size_t n_atoms     = conf.AtomCount();
    const std::size_t n_residues  = protein.ResidueCount();

    auto result_ptr = std::make_unique<LarsenHBondShieldingResult>();
    result_ptr->conf_ = &conf;

    GeometryChoiceBuilder choices(conf);

    // Track which amide H atoms received any contribution (for Δσ_w
    // sweep at the end).
    std::vector<bool> amide_h_has_pair(n_atoms, false);

    // ------------------------------------------------------------------
    // Step 1: DSSP-driven pair resolution. Same iteration shape as
    // HBondResult — both directions per residue, deduped via set keyed
    // on (donor_N, acceptor_O).
    //
    // Phase 1: AmideHydrogen donor + BackboneCarbonyl acceptor only.
    // Phase 2 adds spatial-search-driven Hα donors and sidechain-O
    // acceptors.
    // ------------------------------------------------------------------

    struct ResolvedPair {
        std::size_t donor_residue_idx;
        std::size_t acceptor_residue_idx;
        std::size_t donor_N_idx;
        std::size_t donor_H_idx;
        std::size_t acceptor_O_idx;
        std::size_t acceptor_C_idx;
    };
    std::vector<ResolvedPair> resolved_pairs;
    std::set<std::pair<std::size_t, std::size_t>> seen;
    std::size_t resolution_key = 0;

    for (std::size_t ri = 0; ri < n_residues; ++ri) {
        const auto&     dr  = dssp.AllResidues()[ri];
        const Residue&  res = protein.ResidueAt(ri);

        // This residue's N-H donates to acceptor residues.
        for (int bi = 0; bi < 2; ++bi) {
            std::size_t acc_ri = dr.acceptors[bi].residue_index;
            if (acc_ri == SIZE_MAX || acc_ri >= n_residues) continue;
            const Residue& acc_res = protein.ResidueAt(acc_ri);

            // Donor frame requires: HN (this res), N (this res),
            // C'(i−1) (preceding residue). Skip if any is missing.
            if (res.H == Residue::NONE || res.N == Residue::NONE) continue;
            if (acc_res.O == Residue::NONE || acc_res.C == Residue::NONE) continue;
            if (ri == 0) {
                choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                    "amide donor at N-terminus",
                    [ri](GeometryChoice& gc) {
                        AddNumber(gc, "residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "rejection", 1.0, "no_preceding_C_for_frame");
                    });
                continue;
            }
            auto key = std::make_pair(res.N, acc_res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            ResolvedPair rp;
            rp.donor_residue_idx    = ri;
            rp.acceptor_residue_idx = acc_ri;
            rp.donor_N_idx          = res.N;
            rp.donor_H_idx          = res.H;
            rp.acceptor_O_idx       = acc_res.O;
            rp.acceptor_C_idx       = acc_res.C;
            resolved_pairs.push_back(rp);
        }

        // This residue's C=O accepts from donor residues — same
        // discipline, opposite direction. (HBondResult does this too;
        // DSSP reports both edges so the set-dedup eliminates dupes.)
        for (int bi = 0; bi < 2; ++bi) {
            std::size_t don_ri = dr.donors[bi].residue_index;
            if (don_ri == SIZE_MAX || don_ri >= n_residues) continue;
            const Residue& don_res = protein.ResidueAt(don_ri);
            if (don_res.H == Residue::NONE || don_res.N == Residue::NONE) continue;
            if (res.O == Residue::NONE || res.C == Residue::NONE) continue;
            if (don_ri == 0) continue;  // donor at N-terminus, can't build frame

            auto key = std::make_pair(don_res.N, res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            ResolvedPair rp;
            rp.donor_residue_idx    = don_ri;
            rp.acceptor_residue_idx = ri;
            rp.donor_N_idx          = don_res.N;
            rp.donor_H_idx          = don_res.H;
            rp.acceptor_O_idx       = res.O;
            rp.acceptor_C_idx       = res.C;
            resolved_pairs.push_back(rp);
        }
    }

    // ------------------------------------------------------------------
    // Step 2: For each pair: geometry, grid lookup, Kabsch rotation,
    // Table 2 dispatch, accumulation onto target atoms.
    // ------------------------------------------------------------------

    for (const auto& pair : resolved_pairs) {
        const Residue& don_res = protein.ResidueAt(pair.donor_residue_idx);
        const Residue& acc_res = protein.ResidueAt(pair.acceptor_residue_idx);

        // donor_third = C'(i-1).
        const Residue& prev_res = protein.ResidueAt(pair.donor_residue_idx - 1);
        if (prev_res.C == Residue::NONE) continue;

        // acceptor_third = N(j+1). For C-terminus acceptor, no j+1 →
        // 2° terms cannot be assigned; we still apply 1° terms.
        bool acceptor_has_next = (pair.acceptor_residue_idx + 1 < n_residues);
        std::size_t accthird_idx = Residue::NONE;
        std::size_t i_plus_1 = SIZE_MAX;
        if (acceptor_has_next) {
            i_plus_1 = pair.acceptor_residue_idx + 1;
            const Residue& next_res = protein.ResidueAt(i_plus_1);
            accthird_idx = next_res.N;
            if (accthird_idx == Residue::NONE) {
                // i+1 exists but its N is missing — degenerate frame.
                // Treat as no-i+1.
                acceptor_has_next = false;
                accthird_idx      = Residue::NONE;
                i_plus_1          = SIZE_MAX;
            }
        }

        Vec3 donor_H_pos = conf.PositionAt(pair.donor_H_idx);
        Vec3 donor_N_pos = conf.PositionAt(pair.donor_N_idx);
        Vec3 donor_thd_pos = conf.PositionAt(prev_res.C);
        Vec3 accept_O_pos = conf.PositionAt(pair.acceptor_O_idx);
        Vec3 accept_C_pos = conf.PositionAt(pair.acceptor_C_idx);

        // Need acceptor_third for dihedral. If no i+1, use a synthetic
        // axis from the acceptor's local backbone (e.g. a vector
        // perpendicular to C=O in the residue plane). Simplest: skip
        // pair entirely if no third — Larsen's grid was built with
        // the (N j+1) atom in the dihedral and we can't lookup
        // faithfully without it.
        if (!acceptor_has_next) {
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "C-terminus acceptor",
                [pair](GeometryChoice& gc) {
                    AddNumber(gc, "acceptor_residue",
                        static_cast<double>(pair.acceptor_residue_idx), "index");
                    AddNumber(gc, "rejection", 1.0, "no_next_residue_for_dihedral");
                });
            continue;
        }
        Vec3 accept_thd_pos = conf.PositionAt(accthird_idx);

        // Compute geometry per the contract in LarsenHBondGrid.h.
        LarsenHBondGeometry geom = ComputeLarsenHBondGeometry(
            donor_H_pos, accept_O_pos, accept_C_pos, accept_thd_pos);

        // Out-of-range guard. Grid will also reject; we record the
        // rejection event here.
        if (geom.theta_deg < 90.0 || geom.theta_deg > 180.0) {
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "theta out of range",
                [geom](GeometryChoice& gc) {
                    AddNumber(gc, "theta_deg", geom.theta_deg, "degrees");
                    AddNumber(gc, "rejection", 1.0, "theta_out_of_range");
                });
            continue;
        }

        // Query the grid for amide-H donor × backbone-O acceptor.
        LarsenHBondRecord rec = grid.QueryNearest(
            HBondDonorClass::AmideHydrogen,
            HBondAcceptorClass::BackboneCarbonyl,
            geom);
        if (!rec.IsHit()) {
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "grid query miss",
                [geom](GeometryChoice& gc) {
                    AddNumber(gc, "r_angstrom", geom.r_angstrom, "A");
                    AddNumber(gc, "theta_deg",  geom.theta_deg, "degrees");
                    AddNumber(gc, "rho_deg",    geom.rho_deg, "degrees");
                    AddNumber(gc, "rejection",  1.0, "grid_miss");
                });
            continue;
        }

        // Build donor frame rotation from PROTEIN-side atom positions.
        Mat3 R_protein = ComputeLarsenDonorFrame(
            donor_H_pos, donor_N_pos, donor_thd_pos);

        // Extract donor- and acceptor-side readouts (in canonical frame).
        auto donor_readouts    = ExtractDonorReadouts(rec);
        auto acceptor_readouts = ExtractAcceptorReadouts(rec);

        // Donor-side dispatch: 1°HB readouts apply to donor residue i.
        // For each readout (N, CA, CB, C, HA, HN): rotate to lab, look
        // up Table 2 cell for (target, 1°HB), accumulate.
        using Term = LarsenContribDispatch::Term;
        using TA = LarsenContribDispatch::TargetAtom;
        for (const auto& dr_readout : donor_readouts) {
            if (!dr_readout.present) continue;
            if (!LarsenContribDispatch::Applies(dr_readout.target, Term::Primary_HB))
                continue;
            std::size_t target_ai = TargetAtomIndex(don_res, dr_readout.target);
            if (target_ai == Residue::NONE) continue;
            Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                dr_readout.canonical_tensor, R_protein);
            AccumulateContribution(
                conf.MutableAtomAt(target_ai), Term::Primary_HB, sigma_lab);
            // Cβ diagnostic: even though Table 2 says CB gets no
            // contribution, we still rotate and record the canonical
            // CB tensor (it should produce a near-zero contribution).
            if (dr_readout.target == TA::CB) {
                conf.MutableAtomAt(target_ai).larsen_hbond_diagnostic_CB += sigma_lab;
            }
        }

        // Acceptor-side dispatch: 2°HB readouts apply to acceptor's i+1
        // residue. Same pattern with the 2°HB Table 2 cell.
        const Residue& next_res = protein.ResidueAt(i_plus_1);
        for (const auto& ac_readout : acceptor_readouts) {
            if (!ac_readout.present) continue;
            if (!LarsenContribDispatch::Applies(ac_readout.target, Term::Secondary_HB))
                continue;
            std::size_t target_ai = TargetAtomIndex(next_res, ac_readout.target);
            if (target_ai == Residue::NONE) continue;
            Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                ac_readout.canonical_tensor, R_protein);
            AccumulateContribution(
                conf.MutableAtomAt(target_ai), Term::Secondary_HB, sigma_lab);
        }

        // Bookkeeping: mark this amide H as having a pair (for water-term sweep).
        amide_h_has_pair[pair.donor_H_idx] = true;

        // Per-pair record for diagnostic NPY emission (Phase 2).
        PairRecord pr;
        pr.donor_atom_idx       = pair.donor_H_idx;
        pr.acceptor_atom_idx    = pair.acceptor_O_idx;
        pr.donor_residue_idx    = pair.donor_residue_idx;
        pr.acceptor_residue_idx = pair.acceptor_residue_idx;
        pr.donor_class          = HBondDonorClass::AmideHydrogen;
        pr.acceptor_class       = HBondAcceptorClass::BackboneCarbonyl;
        pr.r_angstrom           = rec.r_angstrom;
        pr.theta_deg            = rec.theta_deg;
        pr.rho_deg              = rec.rho_deg;
        pr.isotropic_total      = (rec.donor_HN.trace() / 3.0);
        pr.any_corner_imputed   = rec.any_corner_imputed;
        result_ptr->pairs_.push_back(pr);

        // Mark any-corner-imputed on every target atom involved in
        // this pair (downstream calculator can downweight).
        if (rec.any_corner_imputed) {
            conf.MutableAtomAt(pair.donor_H_idx).larsen_hbond_any_corner_imputed = true;
            conf.MutableAtomAt(pair.acceptor_O_idx).larsen_hbond_any_corner_imputed = true;
        }
    }

    // ------------------------------------------------------------------
    // Step 3: Water term sweep. For each amide H atom, if zero pairs
    // involved it as donor, set larsen_hbond_water_term = 2.07 ppm.
    // ------------------------------------------------------------------
    for (std::size_t ri = 0; ri < n_residues; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.H == Residue::NONE) continue;  // PRO etc.
        if (amide_h_has_pair[res.H]) continue;
        conf.MutableAtomAt(res.H).larsen_hbond_water_term = kWaterTerm_ppm;
        ++result_ptr->amide_hs_unbound_;
    }

    // ------------------------------------------------------------------
    // Step 4: Compute totals + spherical decomposition + counts.
    // ------------------------------------------------------------------
    for (std::size_t ai = 0; ai < n_atoms; ++ai) {
        ConformationAtom& a = conf.MutableAtomAt(ai);
        Mat3 total = a.larsen_hbond_1pHB_tensor
                   + a.larsen_hbond_2pHB_tensor
                   + a.larsen_hbond_1pHaB_tensor
                   + a.larsen_hbond_2pHaB_tensor;
        a.larsen_hbond_total_tensor = total;
        a.larsen_hbond_total_spherical = SphericalTensor::Decompose(total);
        // Count: non-zero contribution = at least one class's tensor
        // is non-zero. Use a tiny threshold to avoid FP-noise counting.
        const double kThresh = 1e-9;
        bool has_any = a.larsen_hbond_1pHB_tensor.norm()  > kThresh
                    || a.larsen_hbond_2pHB_tensor.norm()  > kThresh
                    || a.larsen_hbond_1pHaB_tensor.norm() > kThresh
                    || a.larsen_hbond_2pHaB_tensor.norm() > kThresh;
        if (has_any) ++result_ptr->atoms_with_contribution_;
    }

    OperationLog::Info(LogCalcOther, "LarsenHBondShieldingResult::Compute",
        "pairs=" + std::to_string(result_ptr->pairs_.size())
        + " atoms_with_contribution=" + std::to_string(result_ptr->atoms_with_contribution_)
        + " amide_hs_with_water_term=" + std::to_string(result_ptr->amide_hs_unbound_));

    return result_ptr;
}


// ============================================================================
// WriteFeatures — Phase 1 emits 6 NPYs
// ============================================================================

int LarsenHBondShieldingResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const std::size_t N = conf.AtomCount();
    if (N == 0) return 0;

    auto pack_mat3 = [&](const Mat3& m, std::vector<double>& out, std::size_t i_row) {
        // Row-major flatten of 3x3, stored as 9 columns per row.
        std::size_t base = i_row * 9;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                out[base + r * 3 + c] = m(r, c);
    };

    // larsen_hbond_total_tensor (N, 9)
    std::vector<double> total(N * 9, std::nan(""));
    // larsen_hbond_1pHB_tensor  (N, 9)
    std::vector<double> ten_1pHB(N * 9, std::nan(""));
    // larsen_hbond_2pHB_tensor  (N, 9)
    std::vector<double> ten_2pHB(N * 9, std::nan(""));
    // larsen_hbond_diagnostic_CB (N, 9)
    std::vector<double> diag_CB(N * 9, std::nan(""));
    // larsen_hbond_water_term (N,)
    std::vector<double> water(N, std::nan(""));
    // larsen_hbond_count (N,) int32
    std::vector<std::int32_t> n_pairs(N, 0);

    for (std::size_t i = 0; i < N; ++i) {
        const ConformationAtom& a = conf.AtomAt(i);
        pack_mat3(a.larsen_hbond_total_tensor, total, i);
        pack_mat3(a.larsen_hbond_1pHB_tensor,  ten_1pHB, i);
        pack_mat3(a.larsen_hbond_2pHB_tensor,  ten_2pHB, i);
        pack_mat3(a.larsen_hbond_diagnostic_CB, diag_CB, i);
        water[i]    = a.larsen_hbond_water_term;
        n_pairs[i]  = a.larsen_hbond_n_pairs;
    }

    fs::path dir(output_dir);
    int n_written = 0;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_total_tensor.npy").string(),
                            total.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_1pHB_tensor.npy").string(),
                            ten_1pHB.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_2pHB_tensor.npy").string(),
                            ten_2pHB.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_diagnostic_CB.npy").string(),
                            diag_CB.data(), N, 9); ++n_written;
    NpyWriter::WriteFloat64((dir / "larsen_hbond_water_term.npy").string(),
                            water.data(), N, 1); ++n_written;
    NpyWriter::WriteInt32((dir / "larsen_hbond_count.npy").string(),
                          n_pairs.data(), N); ++n_written;
    return n_written;
}


}  // namespace nmr
