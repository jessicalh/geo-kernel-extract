#include "LarsenHBondShieldingResult.h"

#include "Bond.h"
#include "ConformationAtom.h"
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


// IsHydroxylOxygen: walk to bonded H atoms; return true if any is
// labelled HydroxylOH_Aliphatic (SER OG, THR OG1) or HydroxylOH_Aromatic
// (TYR OH). Identifies hydroxyl acceptor Os in a substrate-typed way
// (no string traversal on atom names). Also returns the bonded H atom
// index via out-param (becomes the "third atom" for the H-bond dihedral
// when the O is acted on as an acceptor — Larsen's HOMe acceptor model
// uses Hα..O-H as the dihedral reference).
bool IsHydroxylOxygen(const Protein& protein,
                      std::size_t O_idx,
                      std::size_t& bonded_H_out) {
    bonded_H_out = Residue::NONE;
    const auto& o_atom = protein.AtomAt(O_idx);
    if (o_atom.element != Element::O) return false;
    const auto& topo = protein.LegacyAmber();
    for (std::size_t bi : o_atom.bond_indices) {
        const Bond& b = protein.BondAt(bi);
        std::size_t other = (b.atom_index_a == O_idx) ? b.atom_index_b : b.atom_index_a;
        if (protein.AtomAt(other).element != Element::H) continue;
        PolarHKind k = topo.SemanticAt(other).polar_h;
        if (k == PolarHKind::HydroxylOH_Aliphatic ||
            k == PolarHKind::HydroxylOH_Aromatic) {
            bonded_H_out = other;
            return true;
        }
    }
    return false;
}


// AcceptorTriple — resolved (O, C, third) for the dihedral plus
// optional i+1 residue for the Larsen 2° term routing. Built by
// ClassifyAcceptor below; one shape across all 4 acceptor classes.
struct AcceptorTriple {
    HBondAcceptorClass class_;
    std::size_t O_idx;
    std::size_t C_idx;
    std::size_t third_idx;
    std::size_t i_plus_1_residue_idx;  // SIZE_MAX = not applicable
};


// ClassifyAcceptor: given an acceptor candidate O atom, decide its
// chemistry class and resolve the frame anchors. Returns nullopt if
// the atom does not qualify as any Larsen acceptor class.
//
// Class dispatch (substrate-typed):
//   BackboneCarbonyl       — res.O carbonyl O; C = res.C; third = N(j+1)
//                            (no third if at C-terminus → 2° term skipped
//                            but 1° still applies).
//   SidechainCarbonyl      — Asn OD1 / Gln OE1 (PlanarGroupKind::SidechainAmide
//                            + Element::O); C = the sidechain carbonyl C;
//                            third = the sidechain amide N.
//   HydroxylOxygen         — Ser OG / Thr OG1 / Tyr OH (bond-walk for a
//                            HydroxylOH_* polar H); C = the bonded sp3
//                            heavy atom; third = the bonded hydroxyl H.
//   CarboxylateOxygen      — Asp OD1/OD2 / Glu OE1/OE2 / C-term carboxylate
//                            O (PlanarGroupKind::Carboxylate + Element::O);
//                            C = the carboxylate C; third = the OTHER
//                            carboxylate O on the same C.
//
// SidechainCarbonyl uses the NMA acceptor grid as a documented
// approximation (Larsen did not separately scan sidechain primary amide
// acceptors). 2° term is NOT routed for SidechainCarbonyl because the
// i+1 residue mapping doesn't exist outside backbone.
std::optional<AcceptorTriple> ClassifyAcceptor(const Protein& protein,
                                                std::size_t O_idx) {
    const auto& sem = protein.LegacyAmber().SemanticAt(O_idx);
    if (sem.element != Element::O) return std::nullopt;

    AcceptorTriple t{};
    t.O_idx = O_idx;
    t.i_plus_1_residue_idx = SIZE_MAX;

    // (1) BackboneCarbonyl: substrate flag, residue-cached C, i+1 from
    //     next residue's N (same-chain).
    if (sem.IsBackboneCarbonylOxygen()) {
        t.class_ = HBondAcceptorClass::BackboneCarbonyl;
        const auto& o_atom = protein.AtomAt(O_idx);
        const Residue& res_j = protein.ResidueAt(o_atom.residue_index);
        if (res_j.C == Residue::NONE) return std::nullopt;
        t.C_idx = res_j.C;
        // i+1 third: next residue's N, only if same chain.
        std::size_t j_plus_1 = o_atom.residue_index + 1;
        if (j_plus_1 < protein.ResidueCount()) {
            const Residue& next_res = protein.ResidueAt(j_plus_1);
            if (next_res.chain_id == res_j.chain_id &&
                next_res.N != Residue::NONE) {
                t.third_idx = next_res.N;
                t.i_plus_1_residue_idx = j_plus_1;
            } else {
                // C-terminus or chain boundary: no i+1, 2° term skipped.
                t.third_idx = Residue::NONE;
            }
        } else {
            t.third_idx = Residue::NONE;
        }
        return t;
    }

    // (2) SidechainCarbonyl: substrate flag plus locate the bonded C
    //     (carbonyl) and the sidechain amide N.
    if (sem.IsSidechainAmideOxygen()) {
        t.class_ = HBondAcceptorClass::SidechainCarbonyl;
        const auto& o_atom = protein.AtomAt(O_idx);
        std::size_t bonded_C = Residue::NONE;
        for (std::size_t bi : o_atom.bond_indices) {
            const Bond& b = protein.BondAt(bi);
            std::size_t other = (b.atom_index_a == O_idx) ?
                                 b.atom_index_b : b.atom_index_a;
            if (protein.AtomAt(other).element == Element::C) {
                bonded_C = other;
                break;
            }
        }
        if (bonded_C == Residue::NONE) return std::nullopt;
        t.C_idx = bonded_C;
        // third = the sidechain amide N bonded to the carbonyl C
        // (same SidechainAmide planar group).
        std::size_t amide_N = Residue::NONE;
        for (std::size_t bi : protein.AtomAt(bonded_C).bond_indices) {
            const Bond& b = protein.BondAt(bi);
            std::size_t other = (b.atom_index_a == bonded_C) ?
                                 b.atom_index_b : b.atom_index_a;
            if (other == O_idx) continue;
            const auto& other_sem = protein.LegacyAmber().SemanticAt(other);
            if (other_sem.element == Element::N &&
                other_sem.planar_group == PlanarGroupKind::SidechainAmide) {
                amide_N = other;
                break;
            }
        }
        if (amide_N == Residue::NONE) return std::nullopt;
        t.third_idx = amide_N;
        return t;
    }

    // (3) CarboxylateOxygen: substrate flag plus the symmetric partner
    //     pick. The "third" is the OTHER carboxylate O on the same C.
    if (sem.IsSidechainCarboxylateOxygen()) {
        t.class_ = HBondAcceptorClass::CarboxylateOxygen;
        const auto& o_atom = protein.AtomAt(O_idx);
        std::size_t bonded_C = Residue::NONE;
        for (std::size_t bi : o_atom.bond_indices) {
            const Bond& b = protein.BondAt(bi);
            std::size_t other = (b.atom_index_a == O_idx) ?
                                 b.atom_index_b : b.atom_index_a;
            if (protein.AtomAt(other).element == Element::C) {
                bonded_C = other;
                break;
            }
        }
        if (bonded_C == Residue::NONE) return std::nullopt;
        t.C_idx = bonded_C;
        // third = the OTHER carboxylate O on the same C.
        std::size_t other_O = Residue::NONE;
        for (std::size_t bi : protein.AtomAt(bonded_C).bond_indices) {
            const Bond& b = protein.BondAt(bi);
            std::size_t other = (b.atom_index_a == bonded_C) ?
                                 b.atom_index_b : b.atom_index_a;
            if (other == O_idx) continue;
            const auto& other_sem = protein.LegacyAmber().SemanticAt(other);
            if (other_sem.element == Element::O &&
                other_sem.planar_group == PlanarGroupKind::Carboxylate) {
                other_O = other;
                break;
            }
        }
        if (other_O == Residue::NONE) return std::nullopt;
        t.third_idx = other_O;
        return t;
    }

    // (4) HydroxylOxygen: bond-walk discovers the hydroxyl H; the
    //     bonded sp3 C is the "acceptor C"; the hydroxyl H is "third".
    {
        std::size_t hydroxyl_H = Residue::NONE;
        if (IsHydroxylOxygen(protein, O_idx, hydroxyl_H)) {
            t.class_ = HBondAcceptorClass::HydroxylOxygen;
            const auto& o_atom = protein.AtomAt(O_idx);
            std::size_t bonded_C = Residue::NONE;
            for (std::size_t bi : o_atom.bond_indices) {
                const Bond& b = protein.BondAt(bi);
                std::size_t other = (b.atom_index_a == O_idx) ?
                                     b.atom_index_b : b.atom_index_a;
                if (protein.AtomAt(other).element == Element::C) {
                    bonded_C = other;
                    break;
                }
            }
            if (bonded_C == Residue::NONE) return std::nullopt;
            t.C_idx = bonded_C;
            t.third_idx = hydroxyl_H;
            return t;
        }
    }

    return std::nullopt;
}

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

    GeometryChoiceBuilder choices(conf);
    std::size_t resolution_key = 0;

    // amide_h_geometric_paired: true if ANY candidate acceptor was
    // found in the geometric H-bond range for this amide H (any class,
    // pre-grid). Water term applies to amide Hs that have NO candidate
    // — the geometric criterion replaces the DSSP-driven check.
    std::vector<bool> amide_h_geometric_paired(n_atoms, false);
    int n_pairs_grid_skipped = 0;  // geometric pair found, grid path skipped

    auto has_same_chain_prev = [&](std::size_t d_ri) -> bool {
        if (d_ri == 0) return false;
        const Residue& cur  = protein.ResidueAt(d_ri);
        const Residue& prev = protein.ResidueAt(d_ri - 1);
        return cur.chain_id == prev.chain_id;
    };

    // Process one donor → acceptor pair. Resolves donor frame (per
    // donor class), classifies acceptor, computes geometry, queries
    // grid, dispatches per Larsen 2015 Table 2. Returns true if the
    // grid lookup succeeded and contributions were applied.
    auto process_pair =
        [&](HBondDonorClass donor_class,
            std::size_t donor_H_idx,
            std::size_t donor_anchor_idx,
            std::size_t donor_third_idx,
            std::size_t donor_residue_idx,
            const AcceptorTriple& acc) -> bool {

        // Skip degenerate (donor H == any frame anchor) configurations.
        if (donor_anchor_idx == Residue::NONE ||
            donor_third_idx  == Residue::NONE ||
            acc.C_idx        == Residue::NONE ||
            acc.third_idx    == Residue::NONE) {
            return false;
        }

        Vec3 donor_H_pos    = conf.PositionAt(donor_H_idx);
        Vec3 donor_anchor   = conf.PositionAt(donor_anchor_idx);
        Vec3 donor_third    = conf.PositionAt(donor_third_idx);
        Vec3 accept_O_pos   = conf.PositionAt(acc.O_idx);
        Vec3 accept_C_pos   = conf.PositionAt(acc.C_idx);
        Vec3 accept_thd_pos = conf.PositionAt(acc.third_idx);

        LarsenHBondGeometry geom = ComputeLarsenHBondGeometry(
            donor_H_pos, accept_O_pos, accept_C_pos, accept_thd_pos);

        if (geom.theta_deg < kThetaMinDeg || geom.theta_deg > kThetaMaxDeg) {
            choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                "theta out of range",
                [geom, donor_class](GeometryChoice& gc) {
                    AddNumber(gc, "donor_class",
                        static_cast<double>(donor_class), "enum");
                    AddNumber(gc, "theta_deg", geom.theta_deg, "degrees");
                    AddNumber(gc, "rejection", 1.0, "theta_out_of_range");
                });
            ++n_pairs_grid_skipped;
            return false;
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
            ++n_pairs_grid_skipped;
            return false;
        }

        Mat3 R_protein = ComputeLarsenDonorFrame(
            donor_H_pos, donor_anchor, donor_third);

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
                table2_contributors.insert(target_ai);
            }
        }

        // 2° readouts → acceptor residue i+1 atoms, ONLY when the
        // acceptor class carries an i+1 mapping (BackboneCarbonyl).
        // HOMe / Acetate / SidechainCarbonyl (NMA-grid-approximated)
        // have no defined i+1, so 2° terms are skipped for those.
        if (acc.i_plus_1_residue_idx != SIZE_MAX) {
            const Residue& next_res = protein.ResidueAt(acc.i_plus_1_residue_idx);
            for (const auto& ac_readout : acceptor_readouts) {
                if (!ac_readout.present) continue;
                if (!LarsenContribDispatch::Applies(ac_readout.target, secondary_term))
                    continue;
                auto targets = TargetAtomIndices(protein, next_res, ac_readout.target);
                if (targets.empty())     continue;
                Mat3 sigma_lab = RotateTensorToProteinLabFrame(
                    ac_readout.canonical_tensor, R_protein);
                for (std::size_t target_ai : targets) {
                    AccumulateContribution(
                        conf.MutableAtomAt(target_ai), secondary_term, sigma_lab);
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
            if (rec.any_corner_imputed) a.larsen_hbond_any_corner_imputed = true;
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

        PairRecord pr;
        pr.donor_atom_idx       = donor_H_idx;
        pr.acceptor_atom_idx    = acc.O_idx;
        pr.donor_residue_idx    = donor_residue_idx;
        pr.acceptor_residue_idx = protein.AtomAt(acc.O_idx).residue_index;
        pr.donor_class          = donor_class;
        pr.acceptor_class       = acc.class_;
        pr.r_angstrom           = rec.r_angstrom;
        pr.theta_deg            = rec.theta_deg;
        pr.rho_deg              = rec.rho_deg;
        pr.isotropic_total      = (rec.donor_HN.trace() / 3.0);
        pr.any_corner_imputed   = rec.any_corner_imputed;
        result_ptr->pairs_.push_back(pr);
        return true;
    };

    // ------------------------------------------------------------------
    // Donor sweep — one pass over all atoms in the protein, dispatching
    // amide H and α-hydrogen donors to the spatial enumeration. The
    // donor frame anchors per donor class:
    //   AmideHydrogen: anchor = res.N, third = prev_res.C  (same chain;
    //                  prev_res guaranteed via has_same_chain_prev).
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
            // C'(i-1) — chain-boundary check is essential here.
            if (!has_same_chain_prev(atom.residue_index)) {
                choices.Record(CalculatorId::LarsenHBond, resolution_key++,
                    "amide donor at N-term or chain boundary",
                    [ri = atom.residue_index](GeometryChoice& gc) {
                        AddNumber(gc, "residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "rejection", 1.0,
                            "no_preceding_C_in_same_chain");
                    });
                continue;
            }
            const Residue& prev_res =
                protein.ResidueAt(atom.residue_index - 1);
            donor_third_idx = prev_res.C;
        } else if (sem.IsAnyAlphaHydrogen()) {
            donor_class = HBondDonorClass::AlphaHydrogen;
            donor_anchor_idx = res.CA;
            donor_third_idx  = res.N;
        } else {
            continue;  // not a donor candidate
        }
        if (donor_anchor_idx == Residue::NONE ||
            donor_third_idx  == Residue::NONE) continue;

        Vec3 donor_pos = conf.PositionAt(ai);
        auto candidate_atoms =
            spatial.AtomsWithinRadius(donor_pos, kSpatialCutoff_A);

        // For each candidate atom: skip non-O, skip same-residue
        // self-contacts (donor H bonded to its own residue's atoms
        // can't H-bond to them), skip the donor's own anchor and
        // third atom positions. Classify; process.
        bool found_any_geometric_pair = false;
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

            auto classified = ClassifyAcceptor(protein, o_idx);
            if (!classified.has_value()) continue;

            // The geometric pair count is independent of whether the
            // grid path can process it (the same pair-found semantics
            // that DSSP previously gave us for amide-H water-term gate;
            // here generalised to any donor / acceptor).
            found_any_geometric_pair = true;

            (void)process_pair(donor_class, ai, donor_anchor_idx,
                                donor_third_idx, atom.residue_index,
                                *classified);
        }

        if (donor_class == HBondDonorClass::AmideHydrogen &&
            found_any_geometric_pair) {
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
        "pairs_grid_included=" + std::to_string(result_ptr->pairs_.size())
        + " pairs_grid_skipped=" + std::to_string(n_pairs_grid_skipped)
        + " atoms_with_contribution=" + std::to_string(result_ptr->atoms_with_contribution_)
        + " amide_hs_with_water_term=" + std::to_string(result_ptr->amide_hs_unbound_));

    return result_ptr;
}


// ============================================================================
// WriteFeatures — emits SphericalTensor-packed per-class shieldings.
// Emits 8 NPYs: total + 4 per-class shieldings + CB diagnostic
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
