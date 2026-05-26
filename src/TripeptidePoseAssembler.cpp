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


struct KabschResult {
    Mat3   rotation        = Mat3::Identity();
    Vec3   source_centroid = Vec3::Zero();
    Vec3   target_centroid = Vec3::Zero();
    double rmsd            = 0.0;
};


KabschResult KabschAlign(const Vec3 src[3], const Vec3 dst[3]) {
    KabschResult result;
    result.source_centroid = (src[0] + src[1] + src[2]) / 3.0;
    result.target_centroid = (dst[0] + dst[1] + dst[2]) / 3.0;
    Eigen::Matrix<double, 3, 3> src_centered, dst_centered;
    for (int i = 0; i < 3; ++i) {
        src_centered.col(i) = src[i] - result.source_centroid;
        dst_centered.col(i) = dst[i] - result.target_centroid;
    }
    const Mat3 cross_covariance = src_centered * dst_centered.transpose();
    Eigen::JacobiSVD<Mat3> svd(cross_covariance,
        Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Mat3& U = svd.matrixU();
    const Mat3& V = svd.matrixV();
    // Canonical Kabsch: sign(det(V*U^T)) for the reflection guard, not
    // the raw determinant — the product is orthogonal (|det|=1), so
    // only the sign matters.
    const double det = (V * U.transpose()).determinant();
    Eigen::DiagonalMatrix<double, 3> D(1.0, 1.0,
        (det < 0.0) ? -1.0 : 1.0);
    result.rotation = V * D * U.transpose();
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


// Rotate a rank-2 Cartesian tensor by R from DFT to protein frame.
inline Mat3 RotateTensor(const Mat3& sigma, const Mat3& R) {
    return R * sigma * R.transpose();
}


enum class SlotRole {
    BackboneN, BackboneCA, BackboneC, BackboneO,
    BackboneAmideH, BackboneHA, SidechainCB,
};


bool SubstrateRoleMatches(const Protein& protein,
                           size_t protein_atom_idx,
                           SlotRole slot) {
    if (!protein.LegacyAmber().HasAtomSemantic()) {
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
            // GLY HA2/HA3 can carry Locant::Alpha without
            // BackboneRole::AlphaHydrogen; require H to avoid CA matches.
            return sem.backbone_role == BackboneRole::AlphaHydrogen ||
                   (sem.locant == Locant::Alpha &&
                    sem.element == Element::H);
        case SlotRole::SidechainCB:
            return sem.locant == Locant::Beta;
    }
    return false;
}


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

    out_atom.shielding_tensor_aligned =
        RotateTensor(src_atom.shielding_tensor, K.rotation);
    out_atom.shielding_spherical_aligned =
        SphericalTensor::Decompose(out_atom.shielding_tensor_aligned);

    out.aligned_atoms.push_back(std::move(out_atom));
}


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

    if (!protein.LegacyAmber().HasAtomSemantic()) {
        OperationLog::Warn("TripeptidePoseAssembler::AssembleAlaCap",
            "residue " + std::to_string(res.sequence_number) +
            ": LegacyAmberTopology has no atom semantic — typed cap "
            "match is impossible; calling residue is unprotected");
        return false;
    }

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

    // Residuals keep sidechain chi-grid mismatch; no sidechain
    // re-rotation is applied after backbone Kabsch alignment.
    std::vector<Vec3> aligned(larsen.atoms.size());
    for (std::size_t i = 0; i < larsen.atoms.size(); ++i) {
        aligned[i] = ApplyKabsch(K, larsen.atoms[i].position);
    }

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

        // Relaxed identity drops branch and diastereotopic labels only
        // for graph-ambiguous assignments marked by perception.
        const bool relaxed = perceived.canonical_assignment_ambiguous;
        std::vector<std::size_t> candidates =
            candidate_protein_atoms(perceived.identity, relaxed);
        if (candidates.empty()) {
            ++out.n_substrate_disagreements;
            continue;
        }

        std::size_t best_atom = SIZE_MAX;
        double      best_dist = std::numeric_limits<double>::infinity();
        for (std::size_t ai : candidates) {
            const double d = (aligned[i] - conf.PositionAt(ai)).norm();
            if (d < best_dist) { best_atom = ai; best_dist = d; }
        }
        if (best_atom == SIZE_MAX) continue;
        used.insert(best_atom);

        if (best_dist > validation_threshold_A) {
            ++out.n_above_threshold;
        }

        AlignedDftAtom aligned_atom;
        aligned_atom.dft_atom_idx     = perceived.dft_atom_idx;
        aligned_atom.protein_atom_idx = best_atom;
        aligned_atom.element          = perceived.element;
        aligned_atom.aligned_position = aligned[i];
        aligned_atom.residual_vec     =
            aligned_atom.aligned_position - conf.PositionAt(best_atom);
        aligned_atom.residual_distance = best_dist;
        aligned_atom.substrate_role_agrees = true;
        aligned_atom.shielding_tensor_aligned =
            RotateTensor(perceived.shielding_tensor, K.rotation);
        aligned_atom.shielding_spherical_aligned =
            SphericalTensor::Decompose(aligned_atom.shielding_tensor_aligned);
        out.aligned_atoms.push_back(std::move(aligned_atom));
    }

    // Log one residue-level summary rather than one warning per atom.
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

    if (rec.larsen.has_value()) {
        return AssembleCentralTyped(protein, conf, residue_idx, rec,
                                     validation_threshold_A,
                                     substrate_check_strict, out);
    }

    OperationLog::Warn("TripeptidePoseAssembler::AssembleCentral",
        "residue " +
        std::to_string(protein.ResidueAt(residue_idx).sequence_number) +
        " calc_id=" + std::to_string(rec.calc_id) +
        ": no perceived LarsenTripeptide; declining central assembly.");
    return false;
}


}  // anonymous namespace


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
