#include "LarsenHBondShieldingResult.h"

#include "ConformationAtom.h"
#include "DsspResult.h"
#include "GeometryChoice.h"
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
// This same enumeration carries over to Phase 2 (Hα donor work): each
// GLY HA atom enumerates separately as a donor candidate, but at the
// readout-target side the fan-out is captured here.
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

    // Three separate bookkeeping buckets for amide H atoms:
    //   amide_h_dssp_paired : DSSP detected ANY H-bond involving this
    //     amide H (donor side), regardless of whether our grid lookup
    //     could be performed. Used to suppress the Δσ_w water term —
    //     if DSSP saw the bond, it's not solvent-exposed even if our
    //     grid couldn't compute its contribution (e.g. C-terminus
    //     acceptor, out-of-range θ).
    //   amide_h_grid_paired : grid lookup succeeded and a tensor
    //     contribution landed somewhere. Used for atoms-with-
    //     contribution counting.
    //   The two diverge whenever a DSSP pair is detected but the
    //     grid path skips (C-term acceptor, NaN geometry, etc.).
    //     Without this split, those donor Hs were spuriously assigned
    //     the water term despite being H-bonded.
    std::vector<bool> amide_h_dssp_paired(n_atoms, false);
    std::vector<bool> amide_h_grid_paired(n_atoms, false);
    int n_pairs_dssp_only = 0;  // DSSP-detected, grid-skipped — for log

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

    // Helper: given a candidate donor residue index `d_ri`, return true
    // iff a preceding residue d_ri-1 exists in the SAME chain (needed
    // for donor_third = C'(i-1)). The chain-boundary check is essential
    // for multi-chain proteins where the index-adjacent residue may
    // belong to a different chain (no peptide bond).
    auto has_same_chain_prev = [&](std::size_t d_ri) -> bool {
        if (d_ri == 0) return false;
        const Residue& cur  = protein.ResidueAt(d_ri);
        const Residue& prev = protein.ResidueAt(d_ri - 1);
        return cur.chain_id == prev.chain_id;
    };

    // Helper: mark `amide_h_dssp_paired` for a donor amide H that
    // DSSP detected as forming an H-bond, regardless of whether our
    // grid path can later process it. Suppresses the spurious water
    // term on grid-omitted-but-DSSP-detected pairs.
    auto mark_dssp_pair = [&](std::size_t donor_h_idx) {
        if (donor_h_idx != Residue::NONE) {
            amide_h_dssp_paired[donor_h_idx] = true;
        }
    };

    for (std::size_t ri = 0; ri < n_residues; ++ri) {
        const auto&     dr  = dssp.AllResidues()[ri];
        const Residue&  res = protein.ResidueAt(ri);

        // This residue's N-H donates to acceptor residues.
        for (int bi = 0; bi < 2; ++bi) {
            std::size_t acc_ri = dr.acceptors[bi].residue_index;
            if (acc_ri == SIZE_MAX || acc_ri >= n_residues) continue;
            const Residue& acc_res = protein.ResidueAt(acc_ri);

            // Donor amide H must exist for any DSSP H-bond here; mark
            // dssp-paired before any frame-related skip.
            if (res.H == Residue::NONE || res.N == Residue::NONE) continue;
            if (acc_res.O == Residue::NONE || acc_res.C == Residue::NONE) continue;
            mark_dssp_pair(res.H);

            // Frame requires C'(i−1) in the same chain.
            if (!has_same_chain_prev(ri)) {
                choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                    "amide donor at N-term or chain boundary",
                    [ri](GeometryChoice& gc) {
                        AddNumber(gc, "residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "rejection", 1.0, "no_preceding_C_in_same_chain");
                    });
                ++n_pairs_dssp_only;  // DSSP saw it, grid can't process
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
            mark_dssp_pair(don_res.H);

            if (!has_same_chain_prev(don_ri)) {
                ++n_pairs_dssp_only;
                continue;  // can't build frame
            }

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

        // donor_third = C'(i-1) in the SAME chain (already filtered at
        // resolution time via has_same_chain_prev).
        const Residue& prev_res = protein.ResidueAt(pair.donor_residue_idx - 1);
        if (prev_res.C == Residue::NONE) continue;

        // acceptor_third = N(j+1) in the SAME chain. For C-terminus
        // acceptor (or chain-boundary), no j+1 → 2° terms cannot be
        // assigned; we still apply 1° terms.
        bool acceptor_has_next = (pair.acceptor_residue_idx + 1 < n_residues);
        std::size_t accthird_idx = Residue::NONE;
        std::size_t i_plus_1 = SIZE_MAX;
        if (acceptor_has_next) {
            i_plus_1 = pair.acceptor_residue_idx + 1;
            const Residue& next_res = protein.ResidueAt(i_plus_1);
            // Chain-boundary check: if next residue is in a different
            // chain, no peptide bond, no i+1 readout target.
            if (next_res.chain_id != acc_res.chain_id) {
                acceptor_has_next = false;
                accthird_idx      = Residue::NONE;
                i_plus_1          = SIZE_MAX;
            } else {
                accthird_idx = next_res.N;
                if (accthird_idx == Residue::NONE) {
                    acceptor_has_next = false;
                    accthird_idx      = Residue::NONE;
                    i_plus_1          = SIZE_MAX;
                }
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
            ++n_pairs_dssp_only;
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
            ++n_pairs_dssp_only;
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
            ++n_pairs_dssp_only;
            continue;
        }

        // Build donor frame rotation from PROTEIN-side atom positions.
        Mat3 R_protein = ComputeLarsenDonorFrame(
            donor_H_pos, donor_N_pos, donor_thd_pos);

        // Extract donor- and acceptor-side readouts (in canonical frame).
        auto donor_readouts    = ExtractDonorReadouts(rec);
        auto acceptor_readouts = ExtractAcceptorReadouts(rec);

        // Track contributors in TWO sets so the diagnostic CB doesn't
        // contaminate the Table-2-only larsen_hbond_n_pairs counter:
        //
        //   table2_contributors: atoms receiving a contribution under
        //     one of the four real Table 2 classes (1°HB, 2°HB,
        //     1°HαB, 2°HαB). Increment n_pairs + propagate
        //     any_corner_imputed.
        //   diagnostic_contributors: atoms receiving only the Cβ
        //     diagnostic write. Propagate any_corner_imputed only —
        //     n_pairs is the "Table-2-contribution count" per
        //     ConformationAtom doc.
        std::set<std::size_t> table2_contributors;
        std::set<std::size_t> diagnostic_contributors;

        using Term = LarsenContribDispatch::Term;
        using TA   = LarsenContribDispatch::TargetAtom;

        // Cβ DIAGNOSTIC — emitted BEFORE the Table 2 dispatch loop
        // because the Cβ row in Table 2 is intentionally all-false
        // (Larsen 2015 Cβ contribution = zero by construction). The
        // diagnostic exists to verify the parser → loader → frame-
        // rotation pipeline produces a near-zero result where the
        // physics expects it; non-zero is a methodology signal
        // (feedback_methods_accumulate). In Phase 1 the NMA donor
        // archive has no CB readout (NMA has no sidechain), so the
        // diagnostic is trivially zero. Phase 2 ALA donor adds CB and
        // the diagnostic becomes load-bearing.
        for (const auto& dr_readout : donor_readouts) {
            if (!dr_readout.present)         continue;
            if (dr_readout.target != TA::CB) continue;
            auto targets = TargetAtomIndices(protein, don_res, TA::CB);
            if (targets.empty())             continue;
            Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                dr_readout.canonical_tensor, R_protein);
            for (std::size_t target_ai : targets) {
                conf.MutableAtomAt(target_ai).larsen_hbond_diagnostic_CB += sigma_lab;
                diagnostic_contributors.insert(target_ai);
            }
        }

        // Donor-side dispatch: 1°HB readouts apply to donor residue i.
        // CB is handled above (Cβ-diagnostic branch). All other
        // readouts go through the Table 2 dispatch table.
        for (const auto& dr_readout : donor_readouts) {
            if (!dr_readout.present)         continue;
            if (dr_readout.target == TA::CB) continue;  // handled above
            if (!LarsenContribDispatch::Applies(dr_readout.target, Term::Primary_HB))
                continue;
            // TargetAtomIndices returns a LIST — for GLY HA it's two
            // atoms (HA2 + HA3); for everything else, one. The same
            // tensor is applied to each atom in the list (each is an
            // α-hydrogen at its own position, subject to the same
            // H-bond perturbation Larsen's grid encodes).
            auto targets = TargetAtomIndices(protein, don_res, dr_readout.target);
            if (targets.empty())             continue;
            Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                dr_readout.canonical_tensor, R_protein);
            for (std::size_t target_ai : targets) {
                AccumulateContribution(
                    conf.MutableAtomAt(target_ai), Term::Primary_HB, sigma_lab);
                table2_contributors.insert(target_ai);
            }
        }

        // Acceptor-side dispatch: 2°HB readouts apply to acceptor's i+1
        // residue (same chain, already verified above).
        if (acceptor_has_next) {
            const Residue& next_res = protein.ResidueAt(i_plus_1);
            for (const auto& ac_readout : acceptor_readouts) {
                if (!ac_readout.present) continue;
                if (!LarsenContribDispatch::Applies(ac_readout.target, Term::Secondary_HB))
                    continue;
                auto targets = TargetAtomIndices(protein, next_res, ac_readout.target);
                if (targets.empty())     continue;
                Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                    ac_readout.canonical_tensor, R_protein);
                for (std::size_t target_ai : targets) {
                    AccumulateContribution(
                        conf.MutableAtomAt(target_ai), Term::Secondary_HB, sigma_lab);
                    table2_contributors.insert(target_ai);
                }
            }
        }

        // Per-pair bookkeeping:
        //   Table 2 contributors: n_pairs += 1, any_corner_imputed |= imp.
        //   Diagnostic-only contributors (CB): any_corner_imputed |= imp.
        //   (n_pairs counts ONLY real Table 2 contributions per the
        //    ConformationAtom field doc.)
        for (std::size_t ai : table2_contributors) {
            ConformationAtom& a = conf.MutableAtomAt(ai);
            a.larsen_hbond_n_pairs += 1;
            if (rec.any_corner_imputed) a.larsen_hbond_any_corner_imputed = true;
        }
        for (std::size_t ai : diagnostic_contributors) {
            if (table2_contributors.count(ai)) continue;  // already handled
            ConformationAtom& a = conf.MutableAtomAt(ai);
            if (rec.any_corner_imputed) a.larsen_hbond_any_corner_imputed = true;
        }

        // Inclusion GeometryChoice record (per PATTERNS.md "every
        // inclusion, exclusion, and triggered event gets a Record()").
        choices.Record(CalculatorId::LarsenHBond, resolution_key++,
            "pair included",
            [pair, rec,
             n_table2 = table2_contributors.size(),
             n_diag   = diagnostic_contributors.size()]
            (GeometryChoice& gc) {
                AddNumber(gc, "donor_residue",
                    static_cast<double>(pair.donor_residue_idx), "index");
                AddNumber(gc, "acceptor_residue",
                    static_cast<double>(pair.acceptor_residue_idx), "index");
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

        // Bookkeeping: grid lookup succeeded and contributed.
        amide_h_grid_paired[pair.donor_H_idx] = true;

        // Per-pair record stored on the result (for Phase 2 per-pair NPY).
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
        // any_corner_imputed propagated to each target atom above via
        // the contributors_this_pair loop (covers the actual Table 2
        // target atoms, not the donor H / acceptor O which are not
        // readout targets).
    }

    // ------------------------------------------------------------------
    // Step 3: Water term sweep. Δσ_w = 2.07 ppm applies to amide Hs
    // that are SOLVENT-EXPOSED (no H-bond at all). We use
    // amide_h_dssp_paired (NOT amide_h_grid_paired) — if DSSP detected
    // any H-bond involving this amide H, it's not solvent-exposed even
    // if our grid path couldn't compute the contribution (C-terminus
    // acceptor, out-of-range θ, etc.). Spurious water assignment was
    // the codex M2 finding.
    // ------------------------------------------------------------------
    for (std::size_t ri = 0; ri < n_residues; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.H == Residue::NONE) continue;  // PRO etc.
        if (amide_h_dssp_paired[res.H]) continue;
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

    result_ptr->pairs_dssp_only_ = n_pairs_dssp_only;

    OperationLog::Info(LogCalcOther, "LarsenHBondShieldingResult::Compute",
        "pairs_grid_included=" + std::to_string(result_ptr->pairs_.size())
        + " pairs_dssp_only_grid_skipped=" + std::to_string(n_pairs_dssp_only)
        + " atoms_with_contribution=" + std::to_string(result_ptr->atoms_with_contribution_)
        + " amide_hs_with_water_term=" + std::to_string(result_ptr->amide_hs_unbound_));

    return result_ptr;
}


// ============================================================================
// WriteFeatures — emits SphericalTensor-packed per-class shieldings.
// Phase 1 emits 8 NPYs: total + 4 per-class shieldings + CB diagnostic
// + water term + pair count. Packing is T0 (1) + T1 (3) + T2 (5) = 9
// columns per atom, matching the HBondResult convention. The Mat3
// runtime fields stay on ConformationAtom (Pattern 11) but NPY emission
// uses the spherical decomposition — same information, different layout
// downstream consumers key on.
// ============================================================================

static void PackSphericalST(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1 + i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4 + i] = st.T2[i];
}

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

    for (std::size_t i = 0; i < N; ++i) {
        const ConformationAtom& a = conf.AtomAt(i);
        PackSphericalST(a.larsen_hbond_shielding_spherical,    &shielding[i * 9]);
        PackSphericalST(a.larsen_hbond_1pHB_spherical,         &sh_1pHB  [i * 9]);
        PackSphericalST(a.larsen_hbond_2pHB_spherical,         &sh_2pHB  [i * 9]);
        PackSphericalST(a.larsen_hbond_1pHaB_spherical,        &sh_1pHaB [i * 9]);
        PackSphericalST(a.larsen_hbond_2pHaB_spherical,        &sh_2pHaB [i * 9]);
        PackSphericalST(a.larsen_hbond_diagnostic_CB_spherical,&sh_CB    [i * 9]);
        water[i]   = a.larsen_hbond_water_term;
        n_pairs[i] = a.larsen_hbond_n_pairs;
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
    return n_written;
}


}  // namespace nmr
