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

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace nmr {

// Exposed in a named (external-linkage) per-file namespace SOLELY so the
// fixed-coordinate ±60° forcing test can pin THIS production helper directly
// rather than a copy of it (vet finding 2026-07). This is NOT a shared
// utility: each tripeptide result keeps its own file-local helper (the
// neighbor result uses a distinct namespace); nothing is extracted to a
// common header — the trajectory-code reason extraction is blocked still holds.
namespace tripeptide_backbone_dihedral {

// Dihedral in degrees. atan2 returns [-180, 180].
double DihedralDegrees(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const Vec3 b1 = b - a;
    const Vec3 b2 = c - b;
    const Vec3 b3 = d - c;
    const double b2n = b2.norm();
    if (b2n < 1e-10) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const Vec3 n1 = b1.cross(b2);
    const Vec3 n2 = b2.cross(b3);
    if (n1.norm() < 1e-10 || n2.norm() < 1e-10) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const Vec3 m1 = n1.cross(b2 / b2n);
    const double x = n1.dot(n2);
    const double y = m1.dot(n2);
    return std::atan2(y, x) * 180.0 / M_PI;
}

}  // namespace tripeptide_backbone_dihedral

// Production-owned M14 mapping. External linkage in this per-file namespace
// gives the non-skipping forcing test the exact branch used by Compute while
// keeping the policy local to this calculator.
namespace tripeptide_backbone_status {

int HisVariantDiagnosticCode(AminoAcid type, int variant_index) {
    if (type != AminoAcid::HIS)
        return 0;
    if (variant_index == 2)
        return 1;  // HIP
    if (variant_index == 0)
        return 2;  // HID
    if (variant_index == 1)
        return 3;  // HIE
    return 0;
}

TripeptideMatchStatus UnsupportedHisStatus(AminoAcid type, int variant_index) {
    if (type == AminoAcid::HIS && variant_index == 0)
        return TripeptideMatchStatus::UnsupportedHid;
    if (type == AminoAcid::HIS && variant_index == 1)
        return TripeptideMatchStatus::UnsupportedHie;
    return TripeptideMatchStatus::Miss;
}

}  // namespace tripeptide_backbone_status

namespace {

using tripeptide_backbone_dihedral::DihedralDegrees;


char ResidueOneLetterCode(AminoAcid type) {
    if (type == AminoAcid::Unknown)
        return 'X';
    return GetAminoAcidType(type).one_letter_code;
}


// Local SameChain helper retired 2026-05-19: replaced by the canonical
// Protein::BackboneConnected query (covalent C-N bond graph). chain_id
// matching was an ad-hoc inference that missed within-chain numbering
// gaps with intact bonds, antibody insertion-coded structures, and
// cyclic peptides. See PATTERNS.md + OBJECT_MODEL.md "Backbone
// connectivity discipline" + commit log for the substrate-correction
// sweep.


// frame_type → method tag enum (0 unknown, 1 OPBE, 2 ORCA-PBE).
uint8_t MethodTagFromFrameType(const std::string& ft) {
    if (ft == "gaussian_standard_orientation")
        return 1;
    if (ft == "orca_input_orientation")
        return 2;
    return 0;
}


void CopyRecordDiagnostics(TripeptideBackboneShieldingResult::ResidueMatch& rm, const TripeptideDftRecord& entry) {
    if (!entry.IsHit())
        return;
    rm.calc_id = entry.calc_id;
    rm.frame_type = entry.frame_type;
    rm.natural_chi_axes = entry.natural_chi_axes;
    rm.n_chi_axes_used = entry.n_chi_axes_used;
    rm.dropped_chi_distance_deg = entry.dropped_chi_distance_deg;
    rm.target_phi_grid = entry.target_phi_grid_deg;
    rm.target_psi_grid = entry.target_psi_grid_deg;
    rm.phi_db = entry.phi;
    rm.psi_db = entry.psi;
    for (int k = 0; k < 4; ++k) {
        rm.target_chi_grid[k] = entry.target_chi_grid_deg[k];
    }
    rm.chi_db[0] = entry.chi1;
    rm.chi_db[1] = entry.chi2;
    rm.chi_db[2] = entry.chi3;
    rm.chi_db[3] = entry.chi4;
}


}  // anonymous namespace


// Production diagnostic row packer. It has external linkage in this
// per-file named namespace so the non-skipping schema/status forcing test
// executes the exact code used by WriteFeatures rather than a copied layout.
namespace tripeptide_backbone_diagnostics {

std::array<double, 28> PackDiagnosticRow(const TripeptideBackboneShieldingResult::ResidueMatch& rm) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::array<double, 28> row;
    row.fill(nan);

    row[0] = static_cast<double>(rm.calc_id);
    row[1] = static_cast<double>(rm.status);
    row[2] = rm.status == TripeptideMatchStatus::Ok ? 1.0 : 0.0;
    row[3] = rm.backbone_rmsd;
    row[4] = rm.ca_match_dist;
    row[5] = static_cast<double>(MethodTagFromFrameType(rm.frame_type));
    row[6] = static_cast<double>(rm.n_atoms_matched);
    row[7] = static_cast<double>(rm.natural_chi_axes);
    row[8] = static_cast<double>(rm.n_chi_axes_used);
    row[9] = rm.dropped_chi_distance_deg;
    row[10] = rm.phi_actual;
    row[11] = rm.psi_actual;

    if (rm.calc_id != 0) {
        row[12] = static_cast<double>(rm.phi_db);
        row[13] = static_cast<double>(rm.psi_db);
        for (int k = 0; k < rm.natural_chi_axes && k < 4; ++k) {
            row[18 + k] = static_cast<double>(rm.target_chi_grid[k]);
            row[22 + k] = static_cast<double>(rm.chi_db[k]);
        }
    }
    for (int k = 0; k < 4; ++k)
        row[14 + k] = rm.chi_actual[k];

    row[26] = static_cast<double>(rm.residue_type_code);
    row[27] = static_cast<double>(rm.his_variant_hint);
    return row;
}

}  // namespace tripeptide_backbone_diagnostics


// ============================================================================
// Dependencies
// ============================================================================

std::vector<std::type_index> TripeptideBackboneShieldingResult::Dependencies() const {
    return {std::type_index(typeid(GeometryResult)), std::type_index(typeid(EnrichmentResult))};
}


// ============================================================================
// Compute
// ============================================================================

std::unique_ptr<TripeptideBackboneShieldingResult>
TripeptideBackboneShieldingResult::Compute(ProteinConformation& conf, const TripeptideDftTable& table) {
    OperationLog::Scope scope(
        "TripeptideBackboneShieldingResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) + " residues=" + std::to_string(conf.ProteinRef().ResidueCount()));

    if (conf.AtomCount() == 0) {
        OperationLog::Error("TripeptideBackboneShieldingResult::Compute", "Zero atoms.");
        return nullptr;
    }
    if (!table.IsConnected()) {
        OperationLog::Error("TripeptideBackboneShieldingResult::Compute", "TripeptideDftTable not connected.");
        return nullptr;
    }

    auto result = std::make_unique<TripeptideBackboneShieldingResult>();
    result->conf_ = &conf;

    const Protein& protein = conf.ProteinRef();
    const std::size_t N_res = protein.ResidueCount();
    result->residue_matches_.assign(N_res, ResidueMatch{});
    result->matched_dft_atom_idx_by_atom_.assign(conf.AtomCount(), 0);

    double rmsd_sum = 0.0;

    for (std::size_t ri = 0; ri < N_res; ++ri) {
        ++result->residues_attempted_;
        const Residue& res = protein.ResidueAt(ri);
        ResidueMatch&  rm  = result->residue_matches_[ri];

        rm.residue_type_code = static_cast<int>(res.type);
        rm.his_variant_hint =
            tripeptide_backbone_status::HisVariantDiagnosticCode(res.type, res.protonation_variant_index);

        const char letter = ResidueOneLetterCode(res.type);
        if (letter == 'X')
            continue;

        if (res.N == Residue::NONE || res.CA == Residue::NONE || res.C == Residue::NONE) {
            ++result->residues_failed_;
            continue;
        }

        // φ/ψ from backbone-connected neighbours via the canonical
        // Protein::BackbonePredecessor / BackboneSuccessor bond-graph
        // walk. Wrap-correct for cyclic peptides; covers ACE/NME caps
        // and antibody insertion-coded structures by walking the bond
        // graph rather than chain_id labels.
        bool has_phi = false, has_psi = false;
        double phi = 0.0, psi = 0.0;
        if (auto prev_idx = protein.BackbonePredecessor(ri); prev_idx) {
            // Predecessor guarantees prev.C != NONE.
            const std::size_t prev_C = protein.ResidueAt(*prev_idx).C;
            phi = DihedralDegrees(conf.PositionAt(prev_C),
                conf.PositionAt(res.N),
                conf.PositionAt(res.CA),
                conf.PositionAt(res.C));
            has_phi = true;
        }
        if (auto next_idx = protein.BackboneSuccessor(ri); next_idx) {
            // Successor guarantees next.N != NONE.
            const std::size_t next_N = protein.ResidueAt(*next_idx).N;
            psi = DihedralDegrees(conf.PositionAt(res.N),
                conf.PositionAt(res.CA),
                conf.PositionAt(res.C),
                conf.PositionAt(next_N));
            has_psi = true;
        }
        if (has_phi)
            rm.phi_actual = phi;
        if (has_psi)
            rm.psi_actual = psi;
        if (!has_phi || !has_psi)
            continue;

        // χ angles from cached chi atom indices.
        double chis[4] = {0.0, 0.0, 0.0, 0.0};
        int n_chi_set = 0;
        for (int k = 0; k < 4; ++k) {
            if (!res.chi[k].Valid())
                break;
            chis[k] = DihedralDegrees(conf.PositionAt(res.chi[k].a[0]),
                conf.PositionAt(res.chi[k].a[1]),
                conf.PositionAt(res.chi[k].a[2]),
                conf.PositionAt(res.chi[k].a[3]));
            n_chi_set = k + 1;
        }

        for (int k = 0; k < n_chi_set; ++k)
            rm.chi_actual[k] = chis[k];

        // HIS variant hint: pass the protein residue's protonation
        // variant index so perception locks onto the matching HID/HIE/
        // HIP. -1 (no hint) for non-HIS residues.
        const int his_hint = (res.type == AminoAcid::HIS) ? res.protonation_variant_index : -1;

        // The tensorcs15 DB schema as inspected 2026-05-11 carries the
        // HIP variant exclusively for HIS (AHA rows uniformly have 18
        // central atoms = HID + HE2 + HD1). A protein-side HID
        // (variant_idx=0) or HIE (variant_idx=1) residue is explicitly
        // censored before query and yields no σ_BB^i for that residue.
        // Surface this loudly rather than letting the
        // residue silently disappear from the output; the user can
        // then either re-protonate the input to HIP or accept the
        // missing-σ_BB^i for that residue. Trigger to revisit:
        // tensorcs15 ingest gains HID/HIE rows for HIS.
        const TripeptideMatchStatus his_censor =
            tripeptide_backbone_status::UnsupportedHisStatus(res.type, his_hint);
        if (his_censor != TripeptideMatchStatus::Miss) {
            OperationLog::Warn("TripeptideBackboneShieldingResult::Compute",
                               "residue " + std::to_string(res.sequence_number) + " HIS variant_idx=" + std::to_string(his_hint)
                                   + " (" + (his_hint == 0 ? "HID" : "HIE")
                                   + ") is not in tensorcs15 (which carries HIP only); "
                                     "the lookup is censored and σ_BB^i will be absent for "
                "this residue. Re-protonate the input to HIP if "
                "σ_BB^i coverage at this residue is required.");
            rm.status = his_censor;
            ++result->residues_failed_;
            continue;
        }

        // Query DB with chi-fallback.
        TripeptideDftRecord entry;
        TripeptideDftRecord perception_failed_entry;
        bool saw_perception_failure = false;
        bool entry_usable = false;
        for (int n_chi_used = n_chi_set; n_chi_used >= 0; --n_chi_used) {
            TripeptideDftRecord candidate;
            try {
                candidate = table.QueryNearest(letter, phi, psi, chis[0], chis[1], chis[2], chis[3], n_chi_used, his_hint);
            } catch (const std::exception& e) {
                OperationLog::Warn("TripeptideBackboneShieldingResult::Compute",
                                   "DB query failed at residue " + std::to_string(res.sequence_number) + " ("
                                       + std::string(1, letter) + "): " + e.what());
                break;
            }
            // A row was found AND perception produced a typed model.
            // Without the larsen check, a perception failure at chi
            // depth N would lock the residue out instead of falling
            // back to shallower chi where a clean row exists.
            if (candidate.IsHit() && candidate.larsen.has_value()) {
                entry = std::move(candidate);
                entry_usable = true;
                break;
            }
            if (candidate.IsHit() && !candidate.larsen.has_value() && !saw_perception_failure) {
                perception_failed_entry = std::move(candidate);
                saw_perception_failure = true;
            }
        }
        if (!entry_usable) {
            if (saw_perception_failure) {
                CopyRecordDiagnostics(rm, perception_failed_entry);
                rm.status = TripeptideMatchStatus::PerceptionFailed;
            }
            ++result->residues_failed_;
            continue;
        }

        // Preserve selected-row fallback/query provenance even if the later
        // pose assembly fails. Assembly failure remains status=miss by the
        // frozen status definition, but the selected calc_id and grids are
        // still valuable diagnostics rather than being silently discarded.
        CopyRecordDiagnostics(rm, entry);

        // Two-path validation assembly.
        AssembledTripeptide asm_ = AssembleTripeptide(protein, conf, ri, entry, TripeptidePoseSide::Central);
        if (!asm_.ok) {
            ++result->residues_failed_;
            continue;
        }

        rm.backbone_rmsd  = asm_.backbone_kabsch_rmsd;

        const uint8_t method_tag = MethodTagFromFrameType(asm_.frame_type);

        for (const auto& a : asm_.aligned_atoms) {
            ConformationAtom& ca = conf.MutableAtomAt(a.protein_atom_idx);
            ca.tripeptide_bb_shielding_tensor    = a.shielding_tensor_aligned;
            ca.tripeptide_bb_shielding_spherical = a.shielding_spherical_aligned;
            ca.tripeptide_bb_match_distance      = a.residual_distance;
            ca.tripeptide_bb_residual_vec        = a.residual_vec;
            ca.tripeptide_bb_has_match           = true;
            ca.tripeptide_bb_method_tag          = method_tag;
            result->matched_dft_atom_idx_by_atom_[a.protein_atom_idx] = a.dft_atom_idx;

            if (a.protein_atom_idx == res.CA) {
                rm.ca_match_dist = a.residual_distance;
            }
            ++result->atoms_assigned_;
        }
        rm.n_atoms_matched = (int)asm_.aligned_atoms.size();
        rm.status = TripeptideMatchStatus::Ok;

        ++result->residues_matched_;
        rmsd_sum += asm_.backbone_kabsch_rmsd;
        if (asm_.backbone_kabsch_rmsd > result->max_backbone_rmsd_) {
            result->max_backbone_rmsd_ = asm_.backbone_kabsch_rmsd;
        }
    }

    if (result->residues_matched_ > 0) {
        result->mean_backbone_rmsd_ = rmsd_sum / (double)result->residues_matched_;
    }

    OperationLog::Info(LogCalcOther,
        "TripeptideBackboneShieldingResult::Compute",
                       std::to_string(result->residues_matched_) + "/" + std::to_string(result->residues_attempted_)
                           + " residues matched, " + std::to_string(result->residues_failed_) + " failed, "
                           + std::to_string(result->atoms_assigned_) + "/" + std::to_string(conf.AtomCount())
                           + " atoms assigned, mean BB RMSD=" + std::to_string(result->mean_backbone_rmsd_)
                           + " Å, max BB RMSD=" + std::to_string(result->max_backbone_rmsd_) + " Å");

    return result;
}


// ============================================================================
// WriteFeatures
// ============================================================================

int TripeptideBackboneShieldingResult::WriteFeatures(const ProteinConformation& conf, const std::string& output_dir) const {
    const std::size_t N = conf.AtomCount();
    int written = 0;
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // tripeptide_bb_shielding.npy (N, 9) float64
    {
        std::vector<double> data(N * 9, kNaN);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            if (ca.tripeptide_bb_has_match) {
                ca.tripeptide_bb_shielding_spherical.PackFull9(&data[i * 9]);
            }
        }
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_bb_shielding.npy", data.data(), N, 9);
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
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_bb_residual_vec.npy", data.data(), N, 3);
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
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_bb_match_distance.npy", data.data(), N);
        ++written;
    }

    // tripeptide_bb_method_tag.npy (N,) int8.
    {
        std::vector<int8_t> data(N, 0);
        for (std::size_t i = 0; i < N; ++i) {
            data[i] = static_cast<int8_t>(conf.AtomAt(i).tripeptide_bb_method_tag);
        }
        NpyWriter::WriteInt8(output_dir + "/tripeptide_bb_method_tag.npy", data.data(), N);
        ++written;
    }

    // tripeptide_bb_match_atoms.npy (N, 5) float64.
    // Columns: residue_index, has_match, matched_dft_atom_idx,
    // match_distance_A, method_tag.
    {
        std::vector<double> data(N * 5, kNaN);
        const Protein& protein = conf.ProteinRef();
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            data[i * 5 + 0] = static_cast<double>(protein.AtomAt(i).residue_index);
            data[i * 5 + 1] = ca.tripeptide_bb_has_match ? 1.0 : 0.0;
            data[i * 5 + 2] = (i < matched_dft_atom_idx_by_atom_.size()) ? static_cast<double>(matched_dft_atom_idx_by_atom_[i])
                    : 0.0;
            if (ca.tripeptide_bb_has_match) {
                data[i * 5 + 3] = ca.tripeptide_bb_match_distance;
            }
            data[i * 5 + 4] = static_cast<double>(ca.tripeptide_bb_method_tag);
        }
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_bb_match_atoms.npy", data.data(), N, 5);
        ++written;
        }

    // tripeptide_bb_diagnostics.npy (R, 28) float64. This is the single
    // owner for C3 fallback provenance and M14 HIS censor status; there are
    // deliberately no separate status sidecars.
    {
        constexpr std::size_t C = 28;
        const std::size_t R = residue_matches_.size();
        std::vector<double> data(R * C, kNaN);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto row = tripeptide_backbone_diagnostics::PackDiagnosticRow(residue_matches_[ri]);
            std::copy(row.begin(), row.end(), data.begin() + ri * C);
        }
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_bb_diagnostics.npy", data.data(), R, C);
        ++written;
    }

    return written;
}


}  // namespace nmr
