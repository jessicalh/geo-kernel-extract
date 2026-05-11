#include "TripeptideBackboneShieldingResult.h"

#include "AminoAcidType.h"
#include "Atom.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "TripeptidePoseAssembler.h"

#include <cmath>
#include <limits>

namespace nmr {

namespace {

// Dihedral in DEGREES, [-180, 180).
double DihedralDegrees(const Vec3& a, const Vec3& b,
                        const Vec3& c, const Vec3& d) {
    const Vec3 b1 = b - a;
    const Vec3 b2 = c - b;
    const Vec3 b3 = d - c;
    const Vec3 n1 = b1.cross(b2);
    const Vec3 n2 = b2.cross(b3);
    const double x = n1.dot(n2);
    const double y = n1.cross(n2).dot(b2.normalized());
    return std::atan2(y, x) * 180.0 / M_PI;
}


char ResidueOneLetterCode(AminoAcid type) {
    if (type == AminoAcid::Unknown) return 'X';
    return GetAminoAcidType(type).one_letter_code;
}


bool SameChain(const Protein& protein, std::size_t a, std::size_t b) {
    return protein.ResidueAt(a).chain_id ==
           protein.ResidueAt(b).chain_id;
}


// frame_type → method tag enum (0 unknown, 1 OPBE, 2 ORCA-PBE).
uint8_t MethodTagFromFrameType(const std::string& ft) {
    if (ft == "gaussian_standard_orientation") return 1;
    if (ft == "orca_input_orientation")        return 2;
    return 0;
}


}  // anonymous namespace


// ============================================================================
// Dependencies
// ============================================================================

std::vector<std::type_index>
TripeptideBackboneShieldingResult::Dependencies() const {
    return {
        std::type_index(typeid(GeometryResult)),
        std::type_index(typeid(EnrichmentResult))
    };
}


// ============================================================================
// Compute
// ============================================================================

std::unique_ptr<TripeptideBackboneShieldingResult>
TripeptideBackboneShieldingResult::Compute(
        ProteinConformation& conf,
        const TripeptideDftTable& table) {

    OperationLog::Scope scope(
        "TripeptideBackboneShieldingResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " residues=" +
            std::to_string(conf.ProteinRef().ResidueCount()));

    if (conf.AtomCount() == 0) {
        OperationLog::Error(
            "TripeptideBackboneShieldingResult::Compute",
            "Zero atoms.");
        return nullptr;
    }
    if (!table.IsConnected()) {
        OperationLog::Error(
            "TripeptideBackboneShieldingResult::Compute",
            "TripeptideDftTable not connected.");
        return nullptr;
    }

    auto result = std::make_unique<TripeptideBackboneShieldingResult>();
    result->conf_ = &conf;

    const Protein& protein = conf.ProteinRef();
    const std::size_t N_res = protein.ResidueCount();
    result->residue_matches_.assign(N_res, ResidueMatch{});

    double rmsd_sum = 0.0;

    for (std::size_t ri = 0; ri < N_res; ++ri) {
        ++result->residues_attempted_;
        const Residue& res = protein.ResidueAt(ri);
        ResidueMatch&  rm  = result->residue_matches_[ri];

        const char letter = ResidueOneLetterCode(res.type);
        if (letter == 'X') continue;

        if (res.N  == Residue::NONE ||
            res.CA == Residue::NONE ||
            res.C  == Residue::NONE) {
            ++result->residues_failed_;
            continue;
        }

        // φ/ψ from same-chain neighbours.
        bool has_phi = false, has_psi = false;
        double phi = 0.0, psi = 0.0;
        if (ri > 0 && SameChain(protein, ri - 1, ri)) {
            const std::size_t prev_C = protein.ResidueAt(ri - 1).C;
            if (prev_C != Residue::NONE) {
                phi = DihedralDegrees(
                    conf.PositionAt(prev_C),
                    conf.PositionAt(res.N),
                    conf.PositionAt(res.CA),
                    conf.PositionAt(res.C));
                has_phi = true;
            }
        }
        if (ri + 1 < N_res && SameChain(protein, ri, ri + 1)) {
            const std::size_t next_N = protein.ResidueAt(ri + 1).N;
            if (next_N != Residue::NONE) {
                psi = DihedralDegrees(
                    conf.PositionAt(res.N),
                    conf.PositionAt(res.CA),
                    conf.PositionAt(res.C),
                    conf.PositionAt(next_N));
                has_psi = true;
            }
        }
        if (!has_phi || !has_psi) continue;

        // χ angles from cached chi atom indices.
        double chis[4] = {0.0, 0.0, 0.0, 0.0};
        int n_chi_set = 0;
        for (int k = 0; k < 4; ++k) {
            if (!res.chi[k].Valid()) break;
            chis[k] = DihedralDegrees(
                conf.PositionAt(res.chi[k].a[0]),
                conf.PositionAt(res.chi[k].a[1]),
                conf.PositionAt(res.chi[k].a[2]),
                conf.PositionAt(res.chi[k].a[3]));
            n_chi_set = k + 1;
        }

        rm.phi_actual = phi;
        rm.psi_actual = psi;
        for (int k = 0; k < 4; ++k) rm.chi_actual[k] = chis[k];

        // HIS variant hint: pass the protein residue's protonation
        // variant index so perception locks onto the matching HID/HIE/
        // HIP. -1 (no hint) for non-HIS residues.
        const int his_hint = (res.type == AminoAcid::HIS)
            ? res.protonation_variant_index : -1;

        // The tensorcs15 DB schema as inspected 2026-05-11 carries the
        // HIP variant exclusively for HIS (AHA rows uniformly have 18
        // central atoms = HID + HE2 + HD1). A protein-side HID
        // (variant_idx=0) or HIE (variant_idx=1) residue will fail
        // perception (piece_n=17 vs canonical HIP=18 atom-count
        // mismatch in the variant try-loop) and yield no σ_BB^i for
        // that residue. Surface this loudly rather than letting the
        // residue silently disappear from the output; the user can
        // then either re-protonate the input to HIP or accept the
        // missing-σ_BB^i for that residue. Trigger to revisit:
        // tensorcs15 ingest gains HID/HIE rows for HIS.
        if (res.type == AminoAcid::HIS &&
            (his_hint == 0 || his_hint == 1)) {
            OperationLog::Warn(
                "TripeptideBackboneShieldingResult::Compute",
                "residue " + std::to_string(res.sequence_number) +
                " HIS variant_idx=" + std::to_string(his_hint) +
                " (" + (his_hint == 0 ? "HID" : "HIE") +
                ") is not in tensorcs15 (which carries HIP only); "
                "perception will fail and σ_BB^i will be absent for "
                "this residue. Re-protonate the input to HIP if "
                "σ_BB^i coverage at this residue is required.");
        }

        // Query DB with chi-fallback.
        TripeptideDftRecord entry;
        int n_chi_used = n_chi_set;
        for (; n_chi_used >= 0; --n_chi_used) {
            try {
                entry = table.QueryNearest(
                    letter, phi, psi,
                    chis[0], chis[1], chis[2], chis[3],
                    n_chi_used, his_hint);
            } catch (const std::exception& e) {
                OperationLog::Warn(
                    "TripeptideBackboneShieldingResult::Compute",
                    "DB query failed at residue " +
                        std::to_string(res.sequence_number) +
                        " (" + std::string(1, letter) + "): " +
                        e.what());
                break;
            }
            // A row was found AND perception produced a typed model.
            // Without the larsen check, a perception failure at chi
            // depth N would lock the residue out instead of falling
            // back to shallower chi where a clean row exists.
            if (entry.IsHit() && entry.larsen.has_value()) break;
        }
        const bool entry_usable = entry.IsHit() && entry.larsen.has_value();
        rm.n_chi_axes_used = (entry_usable ? n_chi_used : 0);
        if (!entry_usable) {
            ++result->residues_failed_;
            continue;
        }

        // Two-path validation assembly.
        AssembledTripeptide asm_ = AssembleTripeptide(
            protein, conf, ri, entry,
            TripeptidePoseSide::Central);
        if (!asm_.ok) {
            ++result->residues_failed_;
            continue;
        }

        rm.calc_id        = asm_.calc_id;
        rm.frame_type     = asm_.frame_type;
        rm.backbone_rmsd  = asm_.backbone_kabsch_rmsd;
        rm.phi_db         = entry.phi;
        rm.psi_db         = entry.psi;
        rm.chi_db[0]      = entry.chi1;
        rm.chi_db[1]      = entry.chi2;
        rm.chi_db[2]      = entry.chi3;
        rm.chi_db[3]      = entry.chi4;

        const uint8_t method_tag =
            MethodTagFromFrameType(asm_.frame_type);

        for (const auto& a : asm_.aligned_atoms) {
            ConformationAtom& ca = conf.MutableAtomAt(a.protein_atom_idx);
            ca.tripeptide_bb_shielding_tensor    = a.shielding_tensor_aligned;
            ca.tripeptide_bb_shielding_spherical = a.shielding_spherical_aligned;
            ca.tripeptide_bb_match_distance      = a.residual_distance;
            ca.tripeptide_bb_residual_vec        = a.residual_vec;
            ca.tripeptide_bb_has_match           = true;
            ca.tripeptide_bb_method_tag          = method_tag;

            if (a.protein_atom_idx == res.CA) {
                rm.ca_match_dist = a.residual_distance;
            }
            ++result->atoms_assigned_;
        }
        rm.n_atoms_matched = (int)asm_.aligned_atoms.size();

        ++result->residues_matched_;
        rmsd_sum += asm_.backbone_kabsch_rmsd;
        if (asm_.backbone_kabsch_rmsd > result->max_backbone_rmsd_) {
            result->max_backbone_rmsd_ = asm_.backbone_kabsch_rmsd;
        }
    }

    if (result->residues_matched_ > 0) {
        result->mean_backbone_rmsd_ =
            rmsd_sum / (double)result->residues_matched_;
    }

    OperationLog::Info(LogCalcOther,
        "TripeptideBackboneShieldingResult::Compute",
        std::to_string(result->residues_matched_) + "/" +
        std::to_string(result->residues_attempted_) +
        " residues matched, " +
        std::to_string(result->residues_failed_) +
        " failed, " +
        std::to_string(result->atoms_assigned_) + "/" +
        std::to_string(conf.AtomCount()) +
        " atoms assigned, mean BB RMSD=" +
        std::to_string(result->mean_backbone_rmsd_) +
        " Å, max BB RMSD=" +
        std::to_string(result->max_backbone_rmsd_) + " Å");

    return result;
}


// ============================================================================
// WriteFeatures
// ============================================================================

namespace {

void PackSphericalTensor(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1 + i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4 + i] = st.T2[i];
}

}  // anonymous namespace


int TripeptideBackboneShieldingResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const std::size_t N = conf.AtomCount();
    int written = 0;
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // tripeptide_bb_shielding.npy (N, 9) float64
    {
        std::vector<double> data(N * 9, kNaN);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            if (ca.tripeptide_bb_has_match) {
                PackSphericalTensor(
                    ca.tripeptide_bb_shielding_spherical,
                    &data[i * 9]);
            }
        }
        NpyWriter::WriteFloat64(
            output_dir + "/tripeptide_bb_shielding.npy",
            data.data(), N, 9);
        ++written;
    }

    // tripeptide_bb_residual_vec.npy (N, 3) float64 — displacement Vec3.
    // The ML model consumes this alongside the tensor; direction +
    // magnitude both load-bearing. NaN where !has_match.
    {
        std::vector<double> data(N * 3, kNaN);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            if (ca.tripeptide_bb_has_match) {
                data[i * 3 + 0] = ca.tripeptide_bb_residual_vec.x();
                data[i * 3 + 1] = ca.tripeptide_bb_residual_vec.y();
                data[i * 3 + 2] = ca.tripeptide_bb_residual_vec.z();
            }
        }
        NpyWriter::WriteFloat64(
            output_dir + "/tripeptide_bb_residual_vec.npy",
            data.data(), N, 3);
        ++written;
    }

    // tripeptide_bb_match_distance.npy (N,) float64 — magnitude, Å.
    // Redundant with residual_vec but cheap and useful for diagnostics.
    {
        std::vector<double> data(N, kNaN);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            if (ca.tripeptide_bb_has_match) {
                data[i] = ca.tripeptide_bb_match_distance;
            }
        }
        NpyWriter::WriteFloat64(
            output_dir + "/tripeptide_bb_match_distance.npy",
            data.data(), N);
        ++written;
    }

    // tripeptide_bb_method_tag.npy (N,) int8.
    {
        std::vector<int8_t> data(N, 0);
        for (std::size_t i = 0; i < N; ++i) {
            data[i] = static_cast<int8_t>(
                conf.AtomAt(i).tripeptide_bb_method_tag);
        }
        NpyWriter::WriteInt8(
            output_dir + "/tripeptide_bb_method_tag.npy",
            data.data(), N);
        ++written;
    }

    return written;
}


}  // namespace nmr
