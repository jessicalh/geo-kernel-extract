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
#include <limits>
#include <unordered_map>

namespace nmr {

namespace {

// AAA reference standard angles per Larsen 2015 Eq 3.
constexpr int kPhiStd = -120;
constexpr int kPsiStd =  140;


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


// One per-direction outcome for the per-residue match record.
struct DirOutcome {
    int    calc_id = 0;
    double rmsd    = 0.0;
    std::string frame_type;
    int    n_atoms_matched = 0;
    // Per-protein-atom map of residual vec for this direction. Used
    // by the caller to project residual onto the corresponding atoms.
    std::unordered_map<std::size_t, Vec3> per_atom_residual;
    // Per-protein-atom rotated tensor for this direction.
    std::unordered_map<std::size_t, Mat3> per_atom_tensor;
};


}  // anonymous namespace


// ============================================================================
// Dependencies
// ============================================================================

std::vector<std::type_index>
TripeptideNeighborShieldingResult::Dependencies() const {
    return {
        std::type_index(typeid(GeometryResult)),
        std::type_index(typeid(EnrichmentResult))
    };
}


// ============================================================================
// Compute
// ============================================================================

std::unique_ptr<TripeptideNeighborShieldingResult>
TripeptideNeighborShieldingResult::Compute(
        ProteinConformation& conf,
        const TripeptideDftTable& table) {

    OperationLog::Scope scope(
        "TripeptideNeighborShieldingResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " residues=" +
            std::to_string(conf.ProteinRef().ResidueCount()));

    if (conf.AtomCount() == 0) {
        OperationLog::Error(
            "TripeptideNeighborShieldingResult::Compute",
            "Zero atoms.");
        return nullptr;
    }
    if (!table.IsConnected()) {
        OperationLog::Error(
            "TripeptideNeighborShieldingResult::Compute",
            "TripeptideDftTable not connected.");
        return nullptr;
    }

    // AAA reference at standard angles. Cached once for all residues.
    TripeptideDftRecord aaa_ref;
    try {
        aaa_ref = table.QueryNearest(
            'A', (double)kPhiStd, (double)kPsiStd);
    } catch (const std::exception& e) {
        OperationLog::Error(
            "TripeptideNeighborShieldingResult::Compute",
            std::string("AAA reference query failed: ") + e.what());
        return nullptr;
    }
    if (!aaa_ref.IsHit()) {
        OperationLog::Error(
            "TripeptideNeighborShieldingResult::Compute",
            "AAA reference at (-120, 140) missing.");
        return nullptr;
    }
    // AAA reference must also perceive cleanly — every per-residue
    // assembly subtracts this reference, and AssembleAlaCap declines
    // records with no LarsenTripeptide. A nullopt here would silently
    // zero every residue's neighbor contribution. Fail loud at module
    // entry instead.
    if (!aaa_ref.larsen.has_value()) {
        OperationLog::Error(
            "TripeptideNeighborShieldingResult::Compute",
            "AAA reference perception failed (calc_id=" +
                std::to_string(aaa_ref.calc_id) +
                "); cannot compute Δσ_BB^{i±1} for any residue. "
                "See prior PerceiveLarsenTripeptide warning.");
        return nullptr;
    }

    auto result = std::make_unique<TripeptideNeighborShieldingResult>();
    result->conf_ = &conf;

    const Protein& protein = conf.ProteinRef();
    const std::size_t N_res = protein.ResidueCount();
    result->residue_matches_.assign(N_res, ResidueMatch{});

    // Reset per-atom accumulators.
    for (std::size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        ca.tripeptide_neighbor_shielding_tensor    = Mat3::Zero();
        ca.tripeptide_neighbor_shielding_spherical = SphericalTensor{};
        ca.tripeptide_neighbor_residual_vec_prev   = Vec3::Zero();
        ca.tripeptide_neighbor_residual_vec_next   = Vec3::Zero();
        ca.tripeptide_neighbor_has_match           = false;
    }

    for (std::size_t ri = 0; ri < N_res; ++ri) {
        ResidueMatch& rm = result->residue_matches_[ri];
        bool any_neighbor_contribution = false;

        auto do_side = [&](int delta_sign,
                            DirOutcome& out) {
            const long long ni_signed =
                (long long)ri + (long long)delta_sign;
            if (ni_signed < 0 || ni_signed >= (long long)N_res) return;
            const std::size_t ni = (std::size_t)ni_signed;
            if (!SameChain(protein, ri, ni)) return;

            const Residue& neigh = protein.ResidueAt(ni);
            const char letter = ResidueOneLetterCode(neigh.type);
            if (letter == 'X') return;
            if (neigh.N  == Residue::NONE ||
                neigh.CA == Residue::NONE ||
                neigh.C  == Residue::NONE) return;

            // Neighbor's actual φ/ψ.
            bool has_phi = false, has_psi = false;
            double phi = 0.0, psi = 0.0;
            if (ni > 0 && SameChain(protein, ni - 1, ni)) {
                const std::size_t prev_C =
                    protein.ResidueAt(ni - 1).C;
                if (prev_C != Residue::NONE) {
                    phi = DihedralDegrees(
                        conf.PositionAt(prev_C),
                        conf.PositionAt(neigh.N),
                        conf.PositionAt(neigh.CA),
                        conf.PositionAt(neigh.C));
                    has_phi = true;
                }
            }
            if (ni + 1 < N_res && SameChain(protein, ni, ni + 1)) {
                const std::size_t next_N = protein.ResidueAt(ni + 1).N;
                if (next_N != Residue::NONE) {
                    psi = DihedralDegrees(
                        conf.PositionAt(neigh.N),
                        conf.PositionAt(neigh.CA),
                        conf.PositionAt(neigh.C),
                        conf.PositionAt(next_N));
                    has_psi = true;
                }
            }
            if (!has_phi || !has_psi) return;

            double chis[4] = {0.0, 0.0, 0.0, 0.0};
            int n_chi_set = 0;
            for (int k = 0; k < 4; ++k) {
                if (!neigh.chi[k].Valid()) break;
                chis[k] = DihedralDegrees(
                    conf.PositionAt(neigh.chi[k].a[0]),
                    conf.PositionAt(neigh.chi[k].a[1]),
                    conf.PositionAt(neigh.chi[k].a[2]),
                    conf.PositionAt(neigh.chi[k].a[3]));
                n_chi_set = k + 1;
            }

            // HIS variant hint: the NEIGHBOR's protonation variant
            // (perception is for the neighbor's tripeptide, where the
            // neighbor sits as central).
            const int neigh_his_hint = (neigh.type == AminoAcid::HIS)
                ? neigh.protonation_variant_index : -1;

            // Query with chi-fallback.
            TripeptideDftRecord rec_axa;
            for (int n_chi_used = n_chi_set;
                 n_chi_used >= 0; --n_chi_used) {
                try {
                    rec_axa = table.QueryNearest(
                        letter, phi, psi,
                        chis[0], chis[1], chis[2], chis[3],
                        n_chi_used, neigh_his_hint);
                } catch (const std::exception& e) {
                    OperationLog::Warn(
                        "TripeptideNeighborShieldingResult::Compute",
                        "DB query failed at neighbor residue " +
                            std::to_string(neigh.sequence_number) +
                            " (" + std::string(1, letter) + "): " +
                            e.what());
                    return;
                }
                // Found a row AND perception produced a typed model.
                // See the matching invariant in
                // TripeptideBackboneShieldingResult::Compute.
                if (rec_axa.IsHit() && rec_axa.larsen.has_value()) break;
            }
            if (!rec_axa.IsHit() || !rec_axa.larsen.has_value()) return;

            // For Δσ_BB^{i-1} (delta_sign=-1): read at C-terminal ALA
            // cap of the (i-1)-centered tripeptide.
            // For Δσ_BB^{i+1} (delta_sign=+1): read at N-terminal ALA
            // cap of the (i+1)-centered tripeptide.
            const TripeptidePoseSide side = (delta_sign == -1)
                ? TripeptidePoseSide::CTerm
                : TripeptidePoseSide::NTerm;

            AssembledTripeptide asm_axa = AssembleTripeptide(
                protein, conf, ri, rec_axa, side);
            AssembledTripeptide asm_aaa = AssembleTripeptide(
                protein, conf, ri, aaa_ref, side);
            if (!asm_axa.ok || !asm_aaa.ok) return;

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

            out.calc_id    = asm_axa.calc_id;
            out.frame_type = asm_axa.frame_type;
            out.rmsd       = asm_axa.backbone_kabsch_rmsd;

            for (const auto& a : asm_axa.aligned_atoms) {
                auto it = aaa_tensor.find(a.protein_atom_idx);
                if (it == aaa_tensor.end()) continue;
                const Mat3 delta =
                    a.shielding_tensor_aligned - it->second;

                // Track the AXA-side residual (the AAA residual is
                // generally similar in magnitude but we store the AXA
                // one — that's the one the calibration responds to).
                out.per_atom_residual[a.protein_atom_idx] =
                    a.residual_vec;
                out.per_atom_tensor[a.protein_atom_idx] = delta;
                ++out.n_atoms_matched;
            }
            if (out.n_atoms_matched > 0) any_neighbor_contribution = true;
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
            if (!have_prev && !have_next) continue;

            Mat3 sum = Mat3::Zero();
            if (have_prev) sum += pit->second;
            if (have_next) sum += nit->second;

            ConformationAtom& ca = conf.MutableAtomAt(ai);
            ca.tripeptide_neighbor_shielding_tensor    = sum;
            ca.tripeptide_neighbor_shielding_spherical = SphericalTensor::Decompose(sum);
            ca.tripeptide_neighbor_has_match           = true;
            if (have_prev) {
                ca.tripeptide_neighbor_residual_vec_prev =
                    prev_out.per_atom_residual.at(ai);
            }
            if (have_next) {
                ca.tripeptide_neighbor_residual_vec_next =
                    next_out.per_atom_residual.at(ai);
            }
            ++result->atoms_accumulated_;
        }

        rm.prev_calc_id          = prev_out.calc_id;
        rm.prev_backbone_rmsd    = prev_out.rmsd;
        rm.prev_frame_type       = prev_out.frame_type;
        rm.prev_n_atoms_matched  = prev_out.n_atoms_matched;
        rm.next_calc_id          = next_out.calc_id;
        rm.next_backbone_rmsd    = next_out.rmsd;
        rm.next_frame_type       = next_out.frame_type;
        rm.next_n_atoms_matched  = next_out.n_atoms_matched;
        if (any_neighbor_contribution) ++result->residues_any_;
    }

    OperationLog::Info(LogCalcOther,
        "TripeptideNeighborShieldingResult::Compute",
        std::to_string(result->residues_any_) + "/" +
        std::to_string(N_res) +
        " residues received >=1 neighbor contribution; " +
        std::to_string(result->atoms_accumulated_) +
        " per-atom Δσ accumulations applied");

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


int TripeptideNeighborShieldingResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const std::size_t N = conf.AtomCount();
    int written = 0;
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // (N, 9) packed irreps for the summed Δσ_{i-1} + Δσ_{i+1} tensor.
    {
        std::vector<double> data(N * 9, kNaN);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            if (ca.tripeptide_neighbor_has_match) {
                PackSphericalTensor(
                    ca.tripeptide_neighbor_shielding_spherical,
                    &data[i * 9]);
            }
        }
        NpyWriter::WriteFloat64(
            output_dir + "/tripeptide_neighbor_shielding.npy",
            data.data(), N, 9);
        ++written;
    }

    // (N, 3) residual vec from i-1 direction. Zero where that
    // direction had no contribution. Fed to ML alongside the tensor.
    {
        std::vector<double> data(N * 3, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            data[i * 3 + 0] = ca.tripeptide_neighbor_residual_vec_prev.x();
            data[i * 3 + 1] = ca.tripeptide_neighbor_residual_vec_prev.y();
            data[i * 3 + 2] = ca.tripeptide_neighbor_residual_vec_prev.z();
        }
        NpyWriter::WriteFloat64(
            output_dir + "/tripeptide_neighbor_residual_vec_prev.npy",
            data.data(), N, 3);
        ++written;
    }

    // (N, 3) residual vec from i+1 direction.
    {
        std::vector<double> data(N * 3, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            data[i * 3 + 0] = ca.tripeptide_neighbor_residual_vec_next.x();
            data[i * 3 + 1] = ca.tripeptide_neighbor_residual_vec_next.y();
            data[i * 3 + 2] = ca.tripeptide_neighbor_residual_vec_next.z();
        }
        NpyWriter::WriteFloat64(
            output_dir + "/tripeptide_neighbor_residual_vec_next.npy",
            data.data(), N, 3);
        ++written;
    }

    return written;
}


}  // namespace nmr
