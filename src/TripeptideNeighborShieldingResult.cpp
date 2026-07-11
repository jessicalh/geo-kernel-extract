#include "TripeptideNeighborShieldingResult.h"

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
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>

namespace nmr {

// Exposed in a named (external-linkage) per-file namespace SOLELY so the
// fixed-coordinate ±60° forcing test can pin THIS production helper directly
// rather than a copy of it (vet finding 2026-07). This is NOT a shared
// utility: each tripeptide result keeps its own file-local helper (the
// backbone result uses a distinct namespace); nothing is extracted to a
// common header — the trajectory-code reason extraction is blocked still holds.
namespace tripeptide_neighbor_dihedral {

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

}  // namespace tripeptide_neighbor_dihedral

// Production-owned M14 mapping for neighbour directions. This remains a
// per-file policy seam and is called by Compute directly.
namespace tripeptide_neighbor_status {

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

}  // namespace tripeptide_neighbor_status

namespace {

using tripeptide_neighbor_dihedral::DihedralDegrees;

// AAA reference standard angles per Larsen 2015 Eq 3.
constexpr int kPhiStd = -120;
constexpr int kPsiStd =  140;


char ResidueOneLetterCode(AminoAcid type) {
    if (type == AminoAcid::Unknown)
        return 'X';
    return GetAminoAcidType(type).one_letter_code;
}

// frame_type → method tag enum (0 unknown, 1 OPBE, 2 ORCA-PBE).
std::uint8_t MethodTagFromFrameType(const std::string& ft) {
    if (ft == "gaussian_standard_orientation")
        return 1;
    if (ft == "orca_input_orientation")
        return 2;
    return 0;
}


void CopyRecordDiagnostics(TripeptideNeighborShieldingResult::SideMatch& sm, const TripeptideDftRecord& entry) {
    if (!entry.IsHit())
        return;
    sm.calc_id = entry.calc_id;
    sm.frame_type = entry.frame_type;
    sm.natural_chi_axes = entry.natural_chi_axes;
    sm.n_chi_axes_used = entry.n_chi_axes_used;
    sm.dropped_chi_distance_deg = entry.dropped_chi_distance_deg;
    sm.target_phi_grid = entry.target_phi_grid_deg;
    sm.target_psi_grid = entry.target_psi_grid_deg;
    sm.phi_db = entry.phi;
    sm.psi_db = entry.psi;
    for (int k = 0; k < 4; ++k) {
        sm.target_chi_grid[k] = entry.target_chi_grid_deg[k];
    }
    sm.chi_db[0] = entry.chi1;
    sm.chi_db[1] = entry.chi2;
    sm.chi_db[2] = entry.chi3;
    sm.chi_db[3] = entry.chi4;
}


// Local SameChain helper retired 2026-05-19: replaced by the canonical
// Protein::BackboneConnected query (covalent C-N bond graph). chain_id
// matching missed within-chain numbering gaps with intact bonds,
// antibody insertion-coded structures, and cyclic peptides. See
// PATTERNS.md + OBJECT_MODEL.md "Backbone connectivity discipline."


// One per-direction outcome for the per-residue match record.
struct DirOutcome {
    TripeptideNeighborShieldingResult::SideMatch match;
    // Per-protein-atom map of residual vec for this direction. Used
    // by the caller to project residual onto the corresponding atoms.
    std::unordered_map<std::size_t, Vec3> per_atom_residual;
    // Per-protein-atom rotated tensor for this direction.
    std::unordered_map<std::size_t, Mat3> per_atom_tensor;
};


}  // anonymous namespace


// Production diagnostic row packer. External linkage in this per-file
// namespace lets the non-skipping forcing test pin the exact (R,59) layout
// and status codes used by WriteFeatures, without duplicating the writer.
namespace tripeptide_neighbor_diagnostics {

std::array<double, 59>
PackDiagnosticRow(const TripeptideNeighborShieldingResult::ResidueMatch& rm, std::uint8_t aaa_frame_type_code) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::array<double, 59> row;
    row.fill(nan);

    const auto pack_side = [&](const TripeptideNeighborShieldingResult::SideMatch& sm, std::size_t base) {
        row[base + 0] = static_cast<double>(sm.status);
        row[base + 1] = static_cast<double>(sm.calc_id);
        row[base + 2] = sm.status == TripeptideMatchStatus::Ok ? 1.0 : 0.0;
        row[base + 3] = static_cast<double>(sm.neighbor_residue_index);
        row[base + 4] = static_cast<double>(sm.neighbor_residue_type_code);
        row[base + 5] = sm.backbone_rmsd;
        row[base + 6] = static_cast<double>(MethodTagFromFrameType(sm.frame_type));
        row[base + 7] = static_cast<double>(sm.n_atoms_matched);
        row[base + 8] = static_cast<double>(sm.natural_chi_axes);
        row[base + 9] = static_cast<double>(sm.n_chi_axes_used);
        row[base + 10] = sm.dropped_chi_distance_deg;
        row[base + 11] = sm.phi_actual;
        row[base + 12] = sm.psi_actual;
        if (sm.calc_id != 0) {
            row[base + 13] = static_cast<double>(sm.phi_db);
            row[base + 14] = static_cast<double>(sm.psi_db);
            for (int k = 0; k < sm.natural_chi_axes && k < 4; ++k) {
                row[base + 19 + k] = static_cast<double>(sm.target_chi_grid[k]);
                row[base + 23 + k] = static_cast<double>(sm.chi_db[k]);
            }
        }
        for (int k = 0; k < 4; ++k) {
            row[base + 15 + k] = sm.chi_actual[k];
        }
        row[base + 27] = static_cast<double>(sm.his_variant_hint);
    };

    pack_side(rm.prev, 0);
    pack_side(rm.next, 28);

    const std::uint8_t prev_code = MethodTagFromFrameType(rm.prev.frame_type);
    const std::uint8_t next_code = MethodTagFromFrameType(rm.next.frame_type);
    const bool mixed = (rm.prev.status == TripeptideMatchStatus::Ok && prev_code != aaa_frame_type_code)
                       || (rm.next.status == TripeptideMatchStatus::Ok && next_code != aaa_frame_type_code);
    const bool any_match = rm.prev.status == TripeptideMatchStatus::Ok || rm.next.status == TripeptideMatchStatus::Ok;
    row[56] = static_cast<double>(aaa_frame_type_code);
    row[57] = mixed ? 1.0 : 0.0;
    row[58] = any_match ? 1.0 : 0.0;
    return row;
}

}  // namespace tripeptide_neighbor_diagnostics


// ============================================================================
// Dependencies
// ============================================================================

std::vector<std::type_index> TripeptideNeighborShieldingResult::Dependencies() const {
    return {std::type_index(typeid(GeometryResult)), std::type_index(typeid(EnrichmentResult))};
}


// ============================================================================
// Compute
// ============================================================================

std::unique_ptr<TripeptideNeighborShieldingResult>
TripeptideNeighborShieldingResult::Compute(ProteinConformation& conf, const TripeptideDftTable& table) {
    OperationLog::Scope scope(
        "TripeptideNeighborShieldingResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) + " residues=" + std::to_string(conf.ProteinRef().ResidueCount()));

    if (conf.AtomCount() == 0) {
        OperationLog::Error("TripeptideNeighborShieldingResult::Compute", "Zero atoms.");
        return nullptr;
    }
    if (!table.IsConnected()) {
        OperationLog::Error("TripeptideNeighborShieldingResult::Compute", "TripeptideDftTable not connected.");
        return nullptr;
    }

    // AAA reference at standard angles. Cached once for all residues.
    TripeptideDftRecord aaa_ref;
    try {
        aaa_ref = table.QueryNearest('A', (double)kPhiStd, (double)kPsiStd);
    } catch (const std::exception& e) {
        OperationLog::Error("TripeptideNeighborShieldingResult::Compute",
            std::string("AAA reference query failed: ") + e.what());
        return nullptr;
    }
    if (!aaa_ref.IsHit()) {
        OperationLog::Error("TripeptideNeighborShieldingResult::Compute", "AAA reference at (-120, 140) missing.");
        return nullptr;
    }
    // AAA reference must also perceive cleanly — every per-residue
    // assembly subtracts this reference, and AssembleAlaCap declines
    // records with no LarsenTripeptide. A nullopt here would silently
    // zero every residue's neighbor contribution. Fail loud at module
    // entry instead.
    if (!aaa_ref.larsen.has_value()) {
        OperationLog::Error("TripeptideNeighborShieldingResult::Compute",
                            "AAA reference perception failed (calc_id=" + std::to_string(aaa_ref.calc_id)
                                + "); cannot compute Δσ_BB^{i±1} for any residue. "
                "See prior PerceiveLarsenTripeptide warning.");
        return nullptr;
    }

    auto result = std::make_unique<TripeptideNeighborShieldingResult>();
    result->conf_ = &conf;
    result->aaa_reference_calc_id_ = aaa_ref.calc_id;
    result->aaa_reference_frame_type_ = aaa_ref.frame_type;
    result->aaa_reference_method_tag_ = MethodTagFromFrameType(aaa_ref.frame_type);

    const Protein& protein = conf.ProteinRef();
    const std::size_t N_res = protein.ResidueCount();
    result->residue_matches_.assign(N_res, ResidueMatch{});
    result->prev_shielding_.assign(conf.AtomCount(), SphericalTensor{});
    result->next_shielding_.assign(conf.AtomCount(), SphericalTensor{});
    result->prev_has_tensor_.assign(conf.AtomCount(), 0);
    result->next_has_tensor_.assign(conf.AtomCount(), 0);

    // Reset per-atom accumulators.
    //
    // Per-direction residual vectors are NaN-initialised (not zero) so
    // that downstream ML can distinguish "no contribution from this
    // direction" from a coincidental zero-magnitude residual. The
    // conditional setters below overwrite only the directions that
    // contributed; absent directions stay NaN and propagate through
    // WriteFeatures to the NPY. The summed shielding tensor stays
    // Zero-initialised because there's no ambiguity: it's the sum, and
    // tripeptide_neighbor_has_match flags whether ANY direction
    // contributed.
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    const Vec3 kNanVec3{kNaN, kNaN, kNaN};
    for (std::size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        ca.tripeptide_neighbor_shielding_tensor    = Mat3::Zero();
        ca.tripeptide_neighbor_shielding_spherical = SphericalTensor{};
        ca.tripeptide_neighbor_residual_vec_prev   = kNanVec3;
        ca.tripeptide_neighbor_residual_vec_next   = kNanVec3;
        ca.tripeptide_neighbor_has_match           = false;
    }

    for (std::size_t ri = 0; ri < N_res; ++ri) {
        ResidueMatch& rm = result->residue_matches_[ri];
        bool any_neighbor_contribution = false;

        auto do_side = [&](int delta_sign, DirOutcome& out) {
            SideMatch& sm = out.match;
            // Walk the backbone via the canonical bond-graph predecessor/
            // successor query. Wrap-correct for cyclic peptides; correct
            // on ACE/NME caps; correct on antibody insertion-coded
            // structures. Replaces the retired delta-arithmetic + chain_id
            // gate (2026-05-19).
            const auto ni_opt = (delta_sign < 0) ? protein.BackbonePredecessor(ri) : protein.BackboneSuccessor(ri);
            if (!ni_opt)
                return;
            const std::size_t ni = *ni_opt;

            const Residue& neigh = protein.ResidueAt(ni);
            sm.neighbor_residue_index = static_cast<int>(ni);
            sm.neighbor_residue_type_code = static_cast<int>(neigh.type);
            sm.his_variant_hint =
                tripeptide_neighbor_status::HisVariantDiagnosticCode(neigh.type, neigh.protonation_variant_index);
            const char letter = ResidueOneLetterCode(neigh.type);
            if (letter == 'X')
                return;
            if (neigh.N == Residue::NONE || neigh.CA == Residue::NONE || neigh.C == Residue::NONE)
                return;

            // Neighbor's actual φ/ψ via the bond graph.
            bool has_phi = false, has_psi = false;
            double phi = 0.0, psi = 0.0;
            if (auto prev = protein.BackbonePredecessor(ni); prev) {
                const std::size_t prev_C = protein.ResidueAt(*prev).C;
                phi = DihedralDegrees(conf.PositionAt(prev_C),
                    conf.PositionAt(neigh.N),
                    conf.PositionAt(neigh.CA),
                    conf.PositionAt(neigh.C));
                has_phi = true;
            }
            if (auto next = protein.BackboneSuccessor(ni); next) {
                const std::size_t next_N = protein.ResidueAt(*next).N;
                psi = DihedralDegrees(conf.PositionAt(neigh.N),
                    conf.PositionAt(neigh.CA),
                    conf.PositionAt(neigh.C),
                    conf.PositionAt(next_N));
                has_psi = true;
            }
            if (has_phi)
                sm.phi_actual = phi;
            if (has_psi)
                sm.psi_actual = psi;
            if (!has_phi || !has_psi)
                return;

            double chis[4] = {0.0, 0.0, 0.0, 0.0};
            int n_chi_set = 0;
            for (int k = 0; k < 4; ++k) {
                if (!neigh.chi[k].Valid())
                    break;
                chis[k] = DihedralDegrees(conf.PositionAt(neigh.chi[k].a[0]),
                    conf.PositionAt(neigh.chi[k].a[1]),
                    conf.PositionAt(neigh.chi[k].a[2]),
                    conf.PositionAt(neigh.chi[k].a[3]));
                n_chi_set = k + 1;
            }

            for (int k = 0; k < n_chi_set; ++k)
                sm.chi_actual[k] = chis[k];

            // HIS variant hint: the NEIGHBOR's protonation variant
            // (perception is for the neighbor's tripeptide, where the
            // neighbor sits as central).
            const int neigh_his_hint = (neigh.type == AminoAcid::HIS) ? neigh.protonation_variant_index : -1;

            // tensorcs15 carries HIP only (see comment in
            // TripeptideBackboneShieldingResult::Compute). A HID/HIE
            // neighbour is explicitly censored before query, so Δσ_BB
            // from THAT direction never contributes to residue `ri`'s sum. The
            // accumulation continues from the other direction if its
            // neighbour is present. Surface the gap per-direction.
            const TripeptideMatchStatus his_censor =
                tripeptide_neighbor_status::UnsupportedHisStatus(neigh.type, neigh_his_hint);
            if (his_censor != TripeptideMatchStatus::Miss) {
                OperationLog::Warn("TripeptideNeighborShieldingResult::Compute",
                                   "neighbour residue " + std::to_string(neigh.sequence_number) + " (HIS variant_idx="
                                       + std::to_string(neigh_his_hint) + ", " + (neigh_his_hint == 0 ? "HID" : "HIE")
                                       + ") is not in tensorcs15 (HIP only); this direction "
                                         "will contribute no Δσ_BB to residue "
                                       + std::to_string(protein.ResidueAt(ri).sequence_number));
                sm.status = his_censor;
                return;
            }

            // Query with chi-fallback.
            TripeptideDftRecord rec_axa;
            TripeptideDftRecord perception_failed_record;
            bool saw_perception_failure = false;
            bool record_usable = false;
            for (int n_chi_used = n_chi_set; n_chi_used >= 0; --n_chi_used) {
                TripeptideDftRecord candidate;
                try {
                    candidate =
                        table.QueryNearest(letter, phi, psi, chis[0], chis[1], chis[2], chis[3], n_chi_used, neigh_his_hint);
                } catch (const std::exception& e) {
                    OperationLog::Warn("TripeptideNeighborShieldingResult::Compute",
                                       "DB query failed at neighbor residue " + std::to_string(neigh.sequence_number) + " ("
                                           + std::string(1, letter) + "): " + e.what());
                    break;
                }
                // Found a row AND perception produced a typed model.
                // See the matching invariant in
                // TripeptideBackboneShieldingResult::Compute.
                if (candidate.IsHit() && candidate.larsen.has_value()) {
                    rec_axa = std::move(candidate);
                    record_usable = true;
                    break;
                }
                if (candidate.IsHit() && !candidate.larsen.has_value() && !saw_perception_failure) {
                    perception_failed_record = std::move(candidate);
                    saw_perception_failure = true;
                }
            }
            if (!record_usable) {
                if (saw_perception_failure) {
                    CopyRecordDiagnostics(sm, perception_failed_record);
                    sm.status = TripeptideMatchStatus::PerceptionFailed;
            }
                return;
            }

            // Keep selected-row query/fallback fields even if cap assembly
            // later fails. Assembly failure itself remains status=miss.
            CopyRecordDiagnostics(sm, rec_axa);

            // For Δσ_BB^{i-1} (delta_sign=-1): read at C-terminal ALA
            // cap of the (i-1)-centered tripeptide.
            // For Δσ_BB^{i+1} (delta_sign=+1): read at N-terminal ALA
            // cap of the (i+1)-centered tripeptide.
            const TripeptidePoseSide side = (delta_sign == -1) ? TripeptidePoseSide::CTerm : TripeptidePoseSide::NTerm;

            AssembledTripeptide asm_axa = AssembleTripeptide(protein, conf, ri, rec_axa, side);
            AssembledTripeptide asm_aaa = AssembleTripeptide(protein, conf, ri, aaa_ref, side);
            if (!asm_axa.ok || !asm_aaa.ok)
                return;

            // Build a protein-atom keyed lookup of AAA tensors for
            // pairing with AXA atoms.
            std::unordered_map<std::size_t, Mat3> aaa_tensor;
            std::unordered_map<std::size_t, Vec3> aaa_residual;
            aaa_tensor.reserve(asm_aaa.aligned_atoms.size());
            aaa_residual.reserve(asm_aaa.aligned_atoms.size());
            for (const auto& a : asm_aaa.aligned_atoms) {
                aaa_tensor[a.protein_atom_idx] = a.shielding_tensor_aligned;
                aaa_residual[a.protein_atom_idx] = a.residual_vec;
            }

            sm.backbone_rmsd = asm_axa.backbone_kabsch_rmsd;

            for (const auto& a : asm_axa.aligned_atoms) {
                auto it = aaa_tensor.find(a.protein_atom_idx);
                if (it == aaa_tensor.end())
                    continue;
                const Mat3 delta = a.shielding_tensor_aligned - it->second;

                // Track the AXA-side residual (the AAA residual is
                // generally similar in magnitude but we store the AXA
                // one — that's the one the calibration responds to).
                out.per_atom_residual[a.protein_atom_idx] = a.residual_vec;
                out.per_atom_tensor[a.protein_atom_idx] = delta;
                ++sm.n_atoms_matched;
            }
            if (sm.n_atoms_matched > 0) {
                sm.status = TripeptideMatchStatus::Ok;
                any_neighbor_contribution = true;
            }
        };

        DirOutcome prev_out, next_out;
        do_side(-1, prev_out);
        do_side(+1, next_out);

        // Accumulate Δσ_{i-1} + Δσ_{i+1} at each protein atom in res i.
        const Residue& res = protein.ResidueAt(ri);
        for (std::size_t ai : res.atom_indices) {
            const auto pit = prev_out.per_atom_tensor.find(ai);
            const auto nit = next_out.per_atom_tensor.find(ai);
            const bool have_prev = (pit != prev_out.per_atom_tensor.end());
            const bool have_next = (nit != next_out.per_atom_tensor.end());
            if (!have_prev && !have_next)
                continue;

            Mat3 sum = Mat3::Zero();
            if (have_prev) {
                sum += pit->second;
                result->prev_shielding_[ai] = SphericalTensor::Decompose(pit->second);
                result->prev_has_tensor_[ai] = 1;
            }
            if (have_next) {
                sum += nit->second;
                result->next_shielding_[ai] = SphericalTensor::Decompose(nit->second);
                result->next_has_tensor_[ai] = 1;
            }

            ConformationAtom& ca = conf.MutableAtomAt(ai);
            ca.tripeptide_neighbor_shielding_tensor    = sum;
            ca.tripeptide_neighbor_shielding_spherical = SphericalTensor::Decompose(sum);
            ca.tripeptide_neighbor_has_match           = true;
            if (have_prev) {
                ca.tripeptide_neighbor_residual_vec_prev = prev_out.per_atom_residual.at(ai);
            }
            if (have_next) {
                ca.tripeptide_neighbor_residual_vec_next = next_out.per_atom_residual.at(ai);
            }
            ++result->atoms_accumulated_;
        }

        rm.prev = std::move(prev_out.match);
        rm.next = std::move(next_out.match);
        rm.prev_calc_id = rm.prev.calc_id;
        rm.prev_backbone_rmsd = rm.prev.backbone_rmsd;
        rm.prev_frame_type = rm.prev.frame_type;
        rm.prev_n_atoms_matched = rm.prev.n_atoms_matched;
        rm.next_calc_id = rm.next.calc_id;
        rm.next_backbone_rmsd = rm.next.backbone_rmsd;
        rm.next_frame_type = rm.next.frame_type;
        rm.next_n_atoms_matched = rm.next.n_atoms_matched;
        if (any_neighbor_contribution)
            ++result->residues_any_;
    }

    OperationLog::Info(LogCalcOther,
        "TripeptideNeighborShieldingResult::Compute",
                       std::to_string(result->residues_any_) + "/" + std::to_string(N_res)
                           + " residues received >=1 neighbor contribution; " + std::to_string(result->atoms_accumulated_)
                           + " per-atom Δσ accumulations applied");

    return result;
}


// ============================================================================
// WriteFeatures
// ============================================================================

int TripeptideNeighborShieldingResult::WriteFeatures(const ProteinConformation& conf, const std::string& output_dir) const {
    const std::size_t N = conf.AtomCount();
    int written = 0;
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // (N, 9) packed irreps for the summed Δσ_{i-1} + Δσ_{i+1} tensor.
    {
        std::vector<double> data(N * 9, kNaN);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            if (ca.tripeptide_neighbor_has_match) {
                ca.tripeptide_neighbor_shielding_spherical.PackFull9(&data[i * 9]);
            }
        }
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_neighbor_shielding.npy", data.data(), N, 9);
        ++written;
    }

    auto emit_direction_shielding =
        [&](const std::string& filename, const std::vector<SphericalTensor>& src, const std::vector<std::uint8_t>& has_tensor) {
            std::vector<double> data(N * 9, kNaN);
            for (std::size_t i = 0; i < N && i < src.size() && i < has_tensor.size(); ++i) {
                if (has_tensor[i]) {
                    src[i].PackFull9(&data[i * 9]);
                }
            }
            NpyWriter::WriteFloat64(output_dir + "/" + filename, data.data(), N, 9);
            ++written;
        };
    emit_direction_shielding("tripeptide_neighbor_shielding_prev.npy", prev_shielding_, prev_has_tensor_);
    emit_direction_shielding("tripeptide_neighbor_shielding_next.npy", next_shielding_, next_has_tensor_);

    // (N, 3) residual vec from i-1 direction. NaN where that direction
    // had no contribution (so downstream ML distinguishes absent
    // contribution from a coincidentally-zero residual). The
    // ConformationAtom field is NaN-initialised in Compute's reset loop
    // and only overwritten when the i-1 direction contributed; the
    // copy here is unconditional and NaN naturally propagates.
    {
        std::vector<double> data(N * 3, kNaN);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            data[i * 3 + 0] = ca.tripeptide_neighbor_residual_vec_prev.x();
            data[i * 3 + 1] = ca.tripeptide_neighbor_residual_vec_prev.y();
            data[i * 3 + 2] = ca.tripeptide_neighbor_residual_vec_prev.z();
        }
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_neighbor_residual_vec_prev.npy", data.data(), N, 3);
        ++written;
    }

    // (N, 3) residual vec from i+1 direction. NaN where that direction
    // had no contribution; see prev-direction comment for rationale.
    {
        std::vector<double> data(N * 3, kNaN);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            data[i * 3 + 0] = ca.tripeptide_neighbor_residual_vec_next.x();
            data[i * 3 + 1] = ca.tripeptide_neighbor_residual_vec_next.y();
            data[i * 3 + 2] = ca.tripeptide_neighbor_residual_vec_next.z();
        }
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_neighbor_residual_vec_next.npy", data.data(), N, 3);
        ++written;
    }

    // tripeptide_neighbor_reference.npy (1, 5) float64.
    // Columns: aaa_calc_id, aaa_frame_type_code, aaa_phi_db_deg,
    // aaa_psi_db_deg, any_mixed_method_flag.
    {
        double any_mixed_method = 0.0;
        for (const auto& rm : residue_matches_) {
            const bool prev_matched = rm.prev_calc_id != 0 && rm.prev_n_atoms_matched > 0;
            const bool next_matched = rm.next_calc_id != 0 && rm.next_n_atoms_matched > 0;
            const std::uint8_t prev_tag = MethodTagFromFrameType(rm.prev_frame_type);
            const std::uint8_t next_tag = MethodTagFromFrameType(rm.next_frame_type);
            if ((prev_matched && prev_tag != 0 && prev_tag != aaa_reference_method_tag_)
                || (next_matched && next_tag != 0 && next_tag != aaa_reference_method_tag_)) {
                any_mixed_method = 1.0;
                break;
            }
        }

        const double data[5] = {
            static_cast<double>(aaa_reference_calc_id_),
            static_cast<double>(aaa_reference_method_tag_),
            static_cast<double>(kPhiStd),
            static_cast<double>(kPsiStd),
            any_mixed_method,
        };
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_neighbor_reference.npy", data, 1, 5);
        ++written;
    }

    // tripeptide_neighbor_diagnostics.npy (R, 59) float64. Columns 0-27
    // describe i-1, 28-55 describe i+1, and 56-58 carry the AAA method,
    // per-row mixed-method flag, and any-match flag. This single table owns
    // both C3 fallback provenance and M14 censor status.
    {
        constexpr std::size_t C = 59;
        const std::size_t R = residue_matches_.size();
        std::vector<double> data(R * C, kNaN);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto row =
                tripeptide_neighbor_diagnostics::PackDiagnosticRow(residue_matches_[ri], aaa_reference_method_tag_);
            std::copy(row.begin(), row.end(), data.begin() + ri * C);
        }
        NpyWriter::WriteFloat64(output_dir + "/tripeptide_neighbor_diagnostics.npy", data.data(), R, C);
        ++written;
    }

    return written;
}


}  // namespace nmr
