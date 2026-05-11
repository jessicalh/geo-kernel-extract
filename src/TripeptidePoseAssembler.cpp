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
    KabschResult r;
    r.source_centroid = (src[0] + src[1] + src[2]) / 3.0;
    r.target_centroid = (dst[0] + dst[1] + dst[2]) / 3.0;
    Eigen::Matrix<double, 3, 3> P, Q;
    for (int i = 0; i < 3; ++i) {
        P.col(i) = src[i] - r.source_centroid;
        Q.col(i) = dst[i] - r.target_centroid;
    }
    const Mat3 H = P * Q.transpose();
    Eigen::JacobiSVD<Mat3> svd(H,
        Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Mat3& U = svd.matrixU();
    const Mat3& V = svd.matrixV();
    const double det = (V * U.transpose()).determinant();
    Eigen::DiagonalMatrix<double, 3> D(1.0, 1.0, det);
    r.rotation = V * D * U.transpose();
    double sumSq = 0.0;
    for (int i = 0; i < 3; ++i) {
        const Vec3 aligned = r.rotation * P.col(i) + r.target_centroid;
        sumSq += (aligned - dst[i]).squaredNorm();
    }
    r.rmsd = std::sqrt(sumSq / 3.0);
    return r;
}


inline Vec3 ApplyKabsch(const KabschResult& K, const Vec3& p) {
    return K.rotation * (p - K.source_centroid) + K.target_centroid;
}


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
        // Substrate not populated (stub-fixture path). Cannot
        // cross-check; permit and let the residual gate catch
        // mismappings.
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
            return sem.backbone_role == BackboneRole::AlphaHydrogen ||
                   (sem.locant == Locant::Alpha);
        case SlotRole::SidechainCB:
            return sem.locant == Locant::Beta;
    }
    return false;
}


// ─────────────────────────────────────────────────────────────────
// Add one aligned cap atom to the result, with both-path validation.
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

    AlignedDftAtom rec;
    rec.dft_atom_idx     = dft_atom_idx;
    rec.protein_atom_idx = protein_atom_idx;
    rec.element          = src_atom.element;

    rec.aligned_position  = ApplyKabsch(K, src_atom.position);
    rec.residual_vec      =
        rec.aligned_position - conf.PositionAt(protein_atom_idx);
    rec.residual_distance = rec.residual_vec.norm();

    rec.substrate_role_agrees =
        SubstrateRoleMatches(protein, protein_atom_idx, expected_slot);

    if (!rec.substrate_role_agrees) {
        ++out.n_substrate_disagreements;
        OperationLog::Warn(
            "TripeptidePoseAssembler::EmitAlignedAtom",
            "substrate typology disagreement: protein atom " +
                std::to_string(protein_atom_idx) +
                " (element " +
                std::to_string(static_cast<int>(rec.element)) +
                ") does not carry the expected role for this cap "
                "slot");
        if (substrate_check_strict) return;
    }

    if (rec.residual_distance > validation_threshold_A) {
        ++out.n_above_threshold;
        return;
    }

    rec.shielding_tensor_aligned =
        RotateTensor(src_atom.shielding_tensor, K.rotation);
    rec.shielding_spherical_aligned =
        SphericalTensor::Decompose(rec.shielding_tensor_aligned);

    out.aligned_atoms.push_back(std::move(rec));
}


// ─────────────────────────────────────────────────────────────────
// Cap assembly (NTerm / CTerm — always ALA).
// ─────────────────────────────────────────────────────────────────

// Resolve a LarsenResidue local-atom index to its 0-indexed position in
// rec.atoms by matching dft_atom_idx. Returns -1 on miss.
static int LarsenLocalToRecIdx(const LarsenResidue& piece, int local_idx,
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

    // Fail-loud parity with AssembleCentralTyped: cross-substrate
    // matching requires typed AtomSemantic on the protein side. Without
    // it, SubstrateRoleMatches returns true unconditionally and the
    // cap's substrate cross-check degenerates to no-op. Reject early so
    // a missing-substrate build doesn't silently emit untyped data.
    if (!protein.LegacyAmber().HasAtomSemantic()) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleAlaCap",
            "residue " + std::to_string(res.sequence_number) +
            ": LegacyAmberTopology has no atom semantic — typed cap "
            "match is impossible; calling residue is unprotected");
        return false;
    }

    // Cap atom indices come from the typed LarsenResidue slots
    // populated by perception. No heuristic fallback: if perception
    // did not produce a LarsenTripeptide for this record, we decline
    // the residue rather than silently degrade to positional-ordering
    // matching (per the project discipline that perception is the
    // single source of truth for DFT-side identity).
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
    const int cap_n  = LarsenLocalToRecIdx(cap, cap.N_idx,  rec);
    const int cap_h  = LarsenLocalToRecIdx(cap, cap.H_idx,  rec);
    const int cap_ca = LarsenLocalToRecIdx(cap, cap.CA_idx, rec);
    const int cap_ha = LarsenLocalToRecIdx(cap, cap.HA_idx, rec);
    const int cap_cb = LarsenLocalToRecIdx(cap, cap.CB_idx, rec);
    const int cap_c  = LarsenLocalToRecIdx(cap, cap.C_idx,  rec);
    const int cap_o  = LarsenLocalToRecIdx(cap, cap.O_idx,  rec);
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
        // s.dft_idx is -1 when the corresponding LarsenResidue slot is
        // absent (HasAllRequiredSlots guarantees N/CA/C/O/H/HA/CB are
        // populated for ALA caps, so this is currently unreachable —
        // but a defensive guard against future invariant relaxation
        // closes the rec.atoms[-1] UB).
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
// Central assembly. Maps the central residue's full atom set to the
// protein residue's atoms via element+nearest-distance (with the
// gotham sidechain re-rotation around CA-CB to absorb chi-grid
// coarseness). Substrate cross-check applies to the BB+Cβ slots only;
// sidechain-beyond-Cβ atoms are trusted to the element+distance
// matcher (which is what Larsen does in their reference impl).
// ─────────────────────────────────────────────────────────────────

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

    // Aligned positions for every perceived atom — single rotation,
    // the BB Kabsch K. NO sidechain re-rotation: a Rsc on positions
    // without the same Rsc applied to the shielding tensor would put
    // the tensor's principal axes in a different frame than the atom's
    // chemical environment. Either both rotate together (consistent
    // frame but discards an ML signal) or neither does (position
    // residual carries χ1-grid coarseness as a Vec3 feature the
    // upstream model can attend to). We choose neither — per
    // feedback_residual_as_ml_feature, direction + magnitude of the
    // residual are exactly what the calibration wants. The earlier
    // sidechain re-rotation Rsc was a holdover from the validation-
    // threshold-rejection design that no longer exists.
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
    auto candidate_protein_atoms = [&](const AtomMechanicalIdentity& pid,
                                        bool relaxed) {
        std::vector<std::size_t> matches;
        for (std::size_t ai : res.atom_indices) {
            if (used.count(ai)) continue;
            if (IdentityCompatible(pid, ProteinIdentityAt(protein, ai),
                                    relaxed)) {
                matches.push_back(ai);
            }
        }
        return matches;
    };

    for (std::size_t i = 0; i < larsen.atoms.size(); ++i) {
        const auto& perc = larsen.atoms[i];

        // ALWAYS use relaxed identity matching (drop BranchAddress +
        // DiastereotopicIndex from the equality check) and let nearest-
        // spatial within the candidate set resolve the within-class
        // assignment. This closes the within-pair scramble for
        // graph-automorphic atom pairs that K=3 Weisfeiler-Lehman in
        // MatchPiece cannot split (CD1↔CD2, CE1↔CE2 on PHE/TYR,
        // NH1↔NH2 on ARG, HD21↔HD22 on ASN, etc.). Strict match would
        // lock in MatchPiece's arbitrary within-class assignment and
        // silently swap tensors across symmetric positions.
        std::vector<std::size_t> cand =
            candidate_protein_atoms(perc.identity, /*relaxed=*/true);
        if (cand.empty()) {
            ++out.n_substrate_disagreements;
            continue;
        }

        // Nearest aligned within the candidate set. With typed-identity
        // matching, the identity equality IS the validation — we do not
        // reject on residual distance because chi-grid coarseness puts
        // deep-sidechain atoms (Arg / Lys terminal-N region) at 2-4 Å
        // residual routinely. Residual is captured per-atom as an ML
        // feature; downstream analysis filters on it.
        std::size_t best_atom = SIZE_MAX;
        double      best_dist = std::numeric_limits<double>::infinity();
        for (std::size_t ai : cand) {
            const double d = (aligned[i] - conf.PositionAt(ai)).norm();
            if (d < best_dist) { best_atom = ai; best_dist = d; }
        }
        if (best_atom == SIZE_MAX) continue;
        used.insert(best_atom);

        // Outlier stat (for diagnostics; does NOT reject).
        if (best_dist > validation_threshold_A) {
            ++out.n_above_threshold;
        }

        AlignedDftAtom adat;
        adat.dft_atom_idx     = perc.dft_atom_idx;
        adat.protein_atom_idx = best_atom;
        adat.element          = perc.element;
        adat.aligned_position = aligned[i];
        adat.residual_vec     = adat.aligned_position - conf.PositionAt(best_atom);
        adat.residual_distance = best_dist;
        // substrate_role_agrees: relaxed identity match with a non-
        // empty candidate set means BOTH sides agree on chemistry at
        // the (element, locant, backbone_role) level; the within-class
        // assignment (BranchAddress / DiastereotopicIndex) is resolved
        // by nearest-spatial. That's a principled match, not a
        // disagreement.
        adat.substrate_role_agrees = true;
        adat.shielding_tensor_aligned =
            RotateTensor(perc.shielding_tensor, K.rotation);
        adat.shielding_spherical_aligned =
            SphericalTensor::Decompose(adat.shielding_tensor_aligned);
        out.aligned_atoms.push_back(std::move(adat));
    }

    // Residue-level diagnostic: if any perceived atoms failed identity
    // match, log a single summary (don't spam per-atom).
    if (out.n_substrate_disagreements > 0) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleCentralTyped",
            "residue " + std::to_string(res.sequence_number) + " " +
            std::string(GetAminoAcidType(res.type).three_letter_code) +
            " calc_id=" + std::to_string(rec.calc_id) +
            ": " + std::to_string(out.n_substrate_disagreements) +
            " perceived atoms had no protein-side identity match — "
            "likely protonation-variant mismatch or non-standard "
            "residue chemistry");
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
