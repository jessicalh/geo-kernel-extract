#include "TripeptidePoseAssembler.h"

#include "AminoAcidType.h"
#include "Atom.h"
#include "ConformationAtom.h"
#include "LegacyAmberTopology.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "SemanticEnums.h"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <set>

namespace nmr {

namespace {


// ─────────────────────────────────────────────────────────────────
// Kabsch
// ─────────────────────────────────────────────────────────────────

struct KabschResult {
    Mat3   rotation        = Mat3::Identity();
    Vec3   source_centroid = Vec3::Zero();
    Vec3   target_centroid = Vec3::Zero();
    double rmsd            = 0.0;
};


KabschResult KabschAlign(const Vec3 src[3], const Vec3 dst[3]) {
    KabschResult result;
    // centroids
    result.source_centroid = (src[0] + src[1] + src[2]) / 3.0;
    result.target_centroid = (dst[0] + dst[1] + dst[2]) / 3.0;
    Eigen::Matrix<double, 3, 3> src_centered, dst_centered;
    for (int i = 0; i < 3; ++i) {
        src_centered.col(i) = src[i] - result.source_centroid;
        dst_centered.col(i) = dst[i] - result.target_centroid;
    }
    // cross-covariance (src_centered · dst_centeredᵀ)
    const Mat3 cross_covariance = src_centered * dst_centered.transpose();
    Eigen::JacobiSVD<Mat3> svd(cross_covariance,
        Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Mat3& U = svd.matrixU();
    const Mat3& V = svd.matrixV();
    // reflection guard
    // Canonical Kabsch: sign(det(V·Uᵀ)) for the reflection guard, not
    // the raw determinant — the product is orthogonal (|det|=1), so
    // only the sign matters.
    const double det = (V * U.transpose()).determinant();
    Eigen::DiagonalMatrix<double, 3> D(1.0, 1.0,
        (det < 0.0) ? -1.0 : 1.0);
    result.rotation = V * D * U.transpose();
    // fitted RMSD
    double sumSq = 0.0;
    for (int i = 0; i < 3; ++i) {
        const Vec3 aligned =
            result.rotation * src_centered.col(i) + result.target_centroid;
        sumSq += (aligned - dst[i]).squaredNorm();
    }
    result.rmsd = std::sqrt(sumSq / 3.0);
    return result;
}


inline Vec3 ApplyKabsch(const KabschResult& K, const Vec3& p) {
    return K.rotation * (p - K.source_centroid) + K.target_centroid;
}


// R σ Rᵀ: rotate a rank-2 Cartesian tensor by R (source→target frame;
// here DFT→protein via K.rotation).
inline Mat3 RotateTensor(const Mat3& sigma, const Mat3& R) {
    return R * sigma * R.transpose();
}


// ─────────────────────────────────────────────────────────────────
// Substrate cross-check: for a given canonical role slot in the cap,
// what's the expected BackboneRole / Locant on the protein side?
// Each slot's expected role is a constant of the canonical ordering;
// the protein atom we mapped to is looked up via SemanticAt and the
// two are compared.
// ─────────────────────────────────────────────────────────────────

enum class SlotRole {
    BackboneN, BackboneCA, BackboneC, BackboneO,
    BackboneAmideH, BackboneHA, SidechainCB,
};


bool SubstrateRoleMatches(const Protein& protein,
                           size_t protein_atom_idx,
                           SlotRole slot) {
    if (!protein.LegacyAmber().HasAtomSemantic()) {
        // Substrate not populated (stub fixture): permit; the residual
        // gate catches mismappings.
        return true;
    }
    const AtomSemanticTable& sem =
        protein.LegacyAmber().SemanticAt(protein_atom_idx);
    switch (slot) {
        case SlotRole::BackboneN:
            return sem.backbone_role == BackboneRole::Nitrogen;
        case SlotRole::BackboneCA:
            return sem.backbone_role == BackboneRole::AlphaCarbon;
        case SlotRole::BackboneC:
            return sem.backbone_role == BackboneRole::CarbonylCarbon;
        case SlotRole::BackboneO:
            return sem.backbone_role == BackboneRole::CarbonylOxygen;
        case SlotRole::BackboneAmideH:
            return sem.backbone_role == BackboneRole::AmideHydrogen;
        case SlotRole::BackboneHA:
            // The locant==Alpha clause is load-bearing for GLY HA2/HA3,
            // stamped (Locant::Alpha, BackboneRole::None) per Markley; the
            // element==H guard closes the CA over-match (unreachable today
            // since this only validates the assigned res.HA).
            return sem.backbone_role == BackboneRole::AlphaHydrogen ||
                   (sem.locant == Locant::Alpha &&
                    sem.element == Element::H);
        case SlotRole::SidechainCB:
            return sem.locant == Locant::Beta;
    }
    return false;
}


// ─────────────────────────────────────────────────────────────────
// Add one aligned cap atom to the result, with substrate and residual
// validation.
// ─────────────────────────────────────────────────────────────────

void EmitAlignedAtom(
        AssembledTripeptide& out,
        const Protein& protein,
        const ProteinConformation& conf,
        const KabschResult& K,
        const TripeptideDftAtom& src_atom,
        int dft_atom_idx,
        std::size_t protein_atom_idx,
        SlotRole expected_slot,
        double validation_threshold_A,
        bool substrate_check_strict) {

    AlignedDftAtom out_atom;
    out_atom.dft_atom_idx     = dft_atom_idx;
    out_atom.protein_atom_idx = protein_atom_idx;
    out_atom.element          = src_atom.element;

    out_atom.aligned_position  = ApplyKabsch(K, src_atom.position);
    out_atom.residual_vec      =
        out_atom.aligned_position - conf.PositionAt(protein_atom_idx);
    out_atom.residual_distance = out_atom.residual_vec.norm();

    out_atom.substrate_role_agrees =
        SubstrateRoleMatches(protein, protein_atom_idx, expected_slot);

    if (!out_atom.substrate_role_agrees) {
        ++out.n_substrate_disagreements;
        OperationLog::Warn(
            "TripeptidePoseAssembler::EmitAlignedAtom",
            "substrate typology disagreement: protein atom " +
                std::to_string(protein_atom_idx) +
                " (element " +
                std::to_string(static_cast<int>(out_atom.element)) +
                ") does not carry the expected role for this cap "
                "slot");
        if (substrate_check_strict) return;
    }

    if (out_atom.residual_distance > validation_threshold_A) {
        ++out.n_above_threshold;
        return;
    }

    // rotate shielding tensor
    out_atom.shielding_tensor_aligned =
        RotateTensor(src_atom.shielding_tensor, K.rotation);
    out_atom.shielding_spherical_aligned =
        SphericalTensor::Decompose(out_atom.shielding_tensor_aligned);

    out.aligned_atoms.push_back(std::move(out_atom));
}


// ─────────────────────────────────────────────────────────────────
// Cap assembly (NTerm / CTerm — always ALA).
// ─────────────────────────────────────────────────────────────────

// Resolve a LarsenResidue local-atom index to its 0-indexed position in
// rec.atoms by matching dft_atom_idx. Returns -1 on miss.
static int LarsenLocalAtomToRecordIndex(const LarsenResidue& piece, int local_idx,
                                const TripeptideDftRecord& rec) {
    if (local_idx < 0 ||
        local_idx >= static_cast<int>(piece.atoms.size())) return -1;
    const int target_dft = piece.atoms[local_idx].dft_atom_idx;
    for (std::size_t k = 0; k < rec.atoms.size(); ++k) {
        if (rec.atoms[k].atom_idx == target_dft) return static_cast<int>(k);
    }
    return -1;
}

bool AssembleAlaCap(
        const Protein& protein,
        const ProteinConformation& conf,
        std::size_t residue_idx,
        const TripeptideDftRecord& rec,
        TripeptidePoseSide side,
        double validation_threshold_A,
        bool substrate_check_strict,
        AssembledTripeptide& out) {

    const Residue& res = protein.ResidueAt(residue_idx);
    if (res.N  == Residue::NONE ||
        res.CA == Residue::NONE ||
        res.C  == Residue::NONE) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleAlaCap",
            "protein residue " +
                std::to_string(res.sequence_number) +
                " has incomplete N/CA/C backbone cache");
        return false;
    }

    // Cross-substrate matching needs typed AtomSemantic; without it
    // SubstrateRoleMatches is a no-op, so decline rather than emit
    // untyped data.
    if (!protein.LegacyAmber().HasAtomSemantic()) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleAlaCap",
            "residue " + std::to_string(res.sequence_number) +
            ": LegacyAmberTopology has no atom semantic — typed cap "
            "match is impossible; calling residue is unprotected");
        return false;
    }

    // Cap slots come from typed LarsenResidue perception; no heuristic
    // fallback (perception is the single source of DFT-side identity).
    if (!rec.larsen.has_value()) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleAlaCap",
            "calc_id=" + std::to_string(rec.calc_id) +
            ": no perceived LarsenTripeptide — declining cap assembly. "
            "See prior PerceiveLarsenTripeptide warning for the failure "
            "reason.");
        return false;
    }
    const LarsenResidue& cap = (side == TripeptidePoseSide::NTerm)
                                ? rec.larsen->n_cap
                                : rec.larsen->c_cap;
    const int cap_n  = LarsenLocalAtomToRecordIndex(cap, cap.N_idx,  rec);
    const int cap_h  = LarsenLocalAtomToRecordIndex(cap, cap.H_idx,  rec);
    const int cap_ca = LarsenLocalAtomToRecordIndex(cap, cap.CA_idx, rec);
    const int cap_ha = LarsenLocalAtomToRecordIndex(cap, cap.HA_idx, rec);
    const int cap_cb = LarsenLocalAtomToRecordIndex(cap, cap.CB_idx, rec);
    const int cap_c  = LarsenLocalAtomToRecordIndex(cap, cap.C_idx,  rec);
    const int cap_o  = LarsenLocalAtomToRecordIndex(cap, cap.O_idx,  rec);
    if (cap_n < 0 || cap_ca < 0 || cap_c < 0 || cap_o < 0) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleAlaCap",
            "calc_id=" + std::to_string(rec.calc_id) +
            ": perceived cap missing BB slots N=" +
            std::to_string(cap.N_idx) + " CA=" +
            std::to_string(cap.CA_idx) + " C=" +
            std::to_string(cap.C_idx) + " O=" +
            std::to_string(cap.O_idx));
        return false;
    }

    // Kabsch from cap N/CA/C to protein residue i N/CA/C.
    const Vec3 src[3] = {
        rec.atoms[cap_n].position,
        rec.atoms[cap_ca].position,
        rec.atoms[cap_c].position
    };
    const Vec3 dst[3] = {
        conf.PositionAt(res.N),
        conf.PositionAt(res.CA),
        conf.PositionAt(res.C)
    };
    const KabschResult K = KabschAlign(src, dst);
    out.rotation              = K.rotation;
    out.source_centroid       = K.source_centroid;
    out.target_centroid       = K.target_centroid;
    out.backbone_kabsch_rmsd  = K.rmsd;

    // Emit each cap-slot atom with the protein-side counterpart.
    struct Slot {
        int           dft_idx;
        std::size_t   protein_idx;
        SlotRole      role;
    };
    const Slot slots[] = {
        {cap_n,  res.N,  SlotRole::BackboneN},
        {cap_h,  res.H,  SlotRole::BackboneAmideH},
        {cap_ca, res.CA, SlotRole::BackboneCA},
        {cap_ha, res.HA, SlotRole::BackboneHA},
        {cap_cb, res.CB, SlotRole::SidechainCB},
        {cap_c,  res.C,  SlotRole::BackboneC},
        {cap_o,  res.O,  SlotRole::BackboneO},
    };
    for (const Slot& s : slots) {
        // Defensive guard against rec.atoms[-1]: optional cap atoms can
        // be absent from either the perceived DFT cap or the protein.
        if (s.dft_idx < 0) continue;
        if (s.protein_idx == Residue::NONE) continue;
        EmitAlignedAtom(out, protein, conf, K,
                        rec.atoms[s.dft_idx], s.dft_idx,
                        s.protein_idx, s.role,
                        validation_threshold_A,
                        substrate_check_strict);
    }
    return true;
}


// ─────────────────────────────────────────────────────────────────
// Typed-identity central assembly. Used when the DFT record carries a
// perceived LarsenTripeptide. Matching is typed-identity equality with
// nearest-spatial tiebreak within equivalence classes (methyl Hs,
// prochiral methylenes whose perception didn't disambiguate diastereo
// position, aromatic CD/CE pairs whose graph signatures are symmetric).
// No element-walks; no nearest-element heuristic. The SER OG ↔ O swap
// is structurally impossible because both atoms carry distinct typed
// identities (Locant::Gamma vs BackboneRole::CarbonylOxygen).
// ─────────────────────────────────────────────────────────────────

bool IdentityCompatible(const AtomMechanicalIdentity& a,
                         const AtomMechanicalIdentity& b,
                         bool relaxed) {
    if (a.element != b.element) return false;
    if (a.locant != b.locant) return false;
    if (a.backbone_role != b.backbone_role) return false;
    if (relaxed) return true;
    if (!(a.branch == b.branch)) return false;
    if (a.di_index != b.di_index) return false;
    return true;
}


AtomMechanicalIdentity ProteinIdentityAt(const Protein& protein,
                                          std::size_t ai) {
    const AtomSemanticTable& sem = protein.LegacyAmber().SemanticAt(ai);
    AtomMechanicalIdentity id;
    id.element       = sem.element;
    id.locant        = sem.locant;
    id.branch        = sem.branch;
    id.di_index      = sem.di_index;
    id.backbone_role = sem.backbone_role;
    return id;
}


bool AssembleCentralTyped(
        const Protein& protein,
        const ProteinConformation& conf,
        std::size_t residue_idx,
        const TripeptideDftRecord& rec,
        double validation_threshold_A,
        bool substrate_check_strict,
        AssembledTripeptide& out) {

    const Residue& res = protein.ResidueAt(residue_idx);
    if (res.N == Residue::NONE || res.CA == Residue::NONE ||
        res.C == Residue::NONE) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleCentralTyped",
            "residue " + std::to_string(res.sequence_number) +
            " has incomplete BB cache (N/CA/C)");
        return false;
    }
    if (!protein.LegacyAmber().HasAtomSemantic()) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleCentralTyped",
            "residue " + std::to_string(res.sequence_number) +
            ": LegacyAmberTopology has no atom semantic — typed match "
            "is impossible; calling residue is unprotected");
        return false;
    }

    const LarsenResidue& larsen = rec.larsen->central;
    if (larsen.N_idx < 0 || larsen.CA_idx < 0 || larsen.C_idx < 0) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleCentralTyped",
            "calc_id=" + std::to_string(rec.calc_id) +
            " perception did not populate central N/CA/C slots: " +
            "N=" + std::to_string(larsen.N_idx) +
            " CA=" + std::to_string(larsen.CA_idx) +
            " C=" + std::to_string(larsen.C_idx));
        return false;
    }

    // Kabsch from perceived N/CA/C to protein N/CA/C.
    const Vec3 src[3] = {
        larsen.atoms[larsen.N_idx].position,
        larsen.atoms[larsen.CA_idx].position,
        larsen.atoms[larsen.C_idx].position,
    };
    const Vec3 dst[3] = {
        conf.PositionAt(res.N),
        conf.PositionAt(res.CA),
        conf.PositionAt(res.C),
    };
    const KabschResult K = KabschAlign(src, dst);
    out.rotation             = K.rotation;
    out.source_centroid      = K.source_centroid;
    out.target_centroid      = K.target_centroid;
    out.backbone_kabsch_rmsd = K.rmsd;

    // Aligned positions via the single BB Kabsch K — no sidechain
    // re-rotation. χ1-grid coarseness is left in residual_vec as an ML
    // feature (feedback_residual_as_ml_feature) rather than rotated away.
    std::vector<Vec3> aligned(larsen.atoms.size());
    for (std::size_t i = 0; i < larsen.atoms.size(); ++i) {
        aligned[i] = ApplyKabsch(K, larsen.atoms[i].position);
    }

    // Per-atom typed-identity match. For each perceived atom, find the
    // protein atom whose typed identity equals (or is compatible with,
    // for diastereotopic / aromatic equivalence classes) the perceived
    // identity. Within an equivalence class, pick by nearest aligned
    // position.
    std::set<std::size_t> used;
    auto candidate_protein_atoms = [&](const AtomMechanicalIdentity& perceived_identity,
                                        bool relaxed) {
        std::vector<std::size_t> matches;
        for (std::size_t ai : res.atom_indices) {
            if (used.count(ai)) continue;
            if (IdentityCompatible(perceived_identity,
                                    ProteinIdentityAt(protein, ai),
                                    relaxed)) {
                matches.push_back(ai);
            }
        }
        return matches;
    };

    for (std::size_t i = 0; i < larsen.atoms.size(); ++i) {
        const auto& perceived = larsen.atoms[i];

        // identity candidates
        // Strict identity for chemistry-distinct branches (e.g. ILE
        // CG1/CG2); relaxed (drop BranchAddress/DiastereotopicIndex,
        // nearest-spatial tiebreak) for graph-automorphic pairs (PHE
        // CD1/CD2, ARG NH1/NH2, methyl Hs) that no WL round can split.
        // See larsen-residue-design-2026-05-11.md.
        const bool relaxed = perceived.canonical_assignment_ambiguous;
        std::vector<std::size_t> candidates =
            candidate_protein_atoms(perceived.identity, relaxed);
        if (candidates.empty()) {
            ++out.n_substrate_disagreements;
            continue;
        }

        // nearest aligned pose
        // The identity equality IS the validation; we do not reject on
        // residual distance because chi-grid coarseness puts deep-
        // sidechain atoms (Arg / Lys terminal-N region) at 2-4 Å residual
        // routinely. Residual is captured per-atom as an ML feature.
        std::size_t best_atom = SIZE_MAX;
        double      best_dist = std::numeric_limits<double>::infinity();
        for (std::size_t ai : candidates) {
            const double d = (aligned[i] - conf.PositionAt(ai)).norm();
            if (d < best_dist) { best_atom = ai; best_dist = d; }
        }
        if (best_atom == SIZE_MAX) continue;
        used.insert(best_atom);

        // Outlier stat (for diagnostics; does NOT reject).
        if (best_dist > validation_threshold_A) {
            ++out.n_above_threshold;
        }

        // emit (mirrors EmitAlignedAtom; central path keeps outliers)
        AlignedDftAtom aligned_atom;
        aligned_atom.dft_atom_idx     = perceived.dft_atom_idx;
        aligned_atom.protein_atom_idx = best_atom;
        aligned_atom.element          = perceived.element;
        aligned_atom.aligned_position = aligned[i];
        aligned_atom.residual_vec     =
            aligned_atom.aligned_position - conf.PositionAt(best_atom);
        aligned_atom.residual_distance = best_dist;
        // substrate_role_agrees: relaxed identity match with a non-
        // empty candidate set means BOTH sides agree on chemistry at
        // the (element, locant, backbone_role) level; the within-class
        // assignment (BranchAddress / DiastereotopicIndex) is resolved
        // by nearest-spatial. That's a principled match, not a
        // disagreement.
        aligned_atom.substrate_role_agrees = true;
        // rotate shielding tensor
        aligned_atom.shielding_tensor_aligned =
            RotateTensor(perceived.shielding_tensor, K.rotation);
        aligned_atom.shielding_spherical_aligned =
            SphericalTensor::Decompose(aligned_atom.shielding_tensor_aligned);
        out.aligned_atoms.push_back(std::move(aligned_atom));
    }

    // Residue-level diagnostic: if any perceived atoms failed identity
    // match, log a single summary (don't spam per-atom). Common cause
    // on production trajectories: disulfide-bonded CYS (state CYX) —
    // tensorcs15 CYS rows were computed against reduced free Cys-SH so
    // the DB row carries an HG hydrogen the disulfide-bonded protein
    // side correctly does not. Unmatched atoms remain unassigned/NaN
    // until CYX DB rows land; matched atoms can still be emitted.
    if (out.n_substrate_disagreements > 0) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleCentralTyped",
            "residue " + std::to_string(res.sequence_number) + " " +
            std::string(GetAminoAcidType(res.type).three_letter_code) +
            " calc_id=" + std::to_string(rec.calc_id) +
            ": " + std::to_string(out.n_substrate_disagreements) +
            " perceived atoms had no protein-side identity match — "
            "likely protonation-variant mismatch (e.g. CYX disulfide "
            "where DB row carries HG that the bonded protein lacks) "
            "or non-standard residue chemistry");
    }

    return true;
}


bool AssembleCentral(
        const Protein& protein,
        const ProteinConformation& conf,
        std::size_t residue_idx,
        const TripeptideDftRecord& rec,
        double validation_threshold_A,
        bool substrate_check_strict,
        AssembledTripeptide& out) {

    // Typed path: when perception succeeded, route through the typed-
    // identity matcher. No mixed-state of typed-BB + heuristic-sidechain.
    if (rec.larsen.has_value()) {
        return AssembleCentralTyped(protein, conf, residue_idx, rec,
                                     validation_threshold_A,
                                     substrate_check_strict, out);
    }

    // Perception or nothing. There is no longer a heuristic fallback;
    // an absent LarsenTripeptide means we decline the residue and the
    // calculator continues with one fewer assignment. The specific
    // perception-failure reason is already in the OperationLog via
    // PerceiveLarsenTripeptide's structured warning.
    OperationLog::Warn("TripeptidePoseAssembler::AssembleCentral",
        "residue " +
        std::to_string(protein.ResidueAt(residue_idx).sequence_number) +
        " calc_id=" + std::to_string(rec.calc_id) +
        ": no perceived LarsenTripeptide; declining central assembly.");
    return false;
}


}  // anonymous namespace


// ============================================================================
// Public API
// ============================================================================

AssembledTripeptide AssembleTripeptide(
        const Protein& protein,
        const ProteinConformation& conf,
        std::size_t protein_residue_idx,
        const TripeptideDftRecord& rec,
        TripeptidePoseSide side,
        double validation_threshold_A,
        bool substrate_check_strict) {

    AssembledTripeptide out;
    out.calc_id    = rec.calc_id;
    out.frame_type = rec.frame_type;
    out.side       = side;
    if (!rec.IsHit()) return out;

    bool ok = false;
    if (side == TripeptidePoseSide::Central) {
        ok = AssembleCentral(protein, conf, protein_residue_idx, rec,
                              validation_threshold_A,
                              substrate_check_strict, out);
    } else {
        ok = AssembleAlaCap(protein, conf, protein_residue_idx, rec,
                             side, validation_threshold_A,
                             substrate_check_strict, out);
    }
    out.ok = ok && !out.aligned_atoms.empty();

    // Aggregate residual stats.
    if (!out.aligned_atoms.empty()) {
        double sum = 0.0;
        for (const auto& a : out.aligned_atoms) {
            sum += a.residual_distance;
            if (a.residual_distance > out.max_residual_A) {
                out.max_residual_A = a.residual_distance;
            }
        }
        out.mean_residual_A = sum / (double)out.aligned_atoms.size();
    }
    return out;
}


}  // namespace nmr
