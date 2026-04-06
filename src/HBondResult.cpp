#include "HBondResult.h"
#include "Protein.h"
#include "DsspResult.h"
#include "SpatialIndexResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <set>

namespace nmr {


std::vector<std::type_index> HBondResult::Dependencies() const {
    return {
        std::type_index(typeid(DsspResult)),
        std::type_index(typeid(SpatialIndexResult))
    };
}


// ============================================================================
// An identified H-bond from DSSP, resolved to atom positions.
//
// DSSP identifies backbone H-bonds by the Kabsch-Sander energy criterion
// between residue pairs. Each H-bond has:
//   - donor residue (N-H donates)
//   - acceptor residue (C=O accepts)
//
// We resolve this to atoms:
//   - donor_N: the backbone N of the donor residue
//   - acceptor_O: the backbone O of the acceptor residue
//   - h_hat: unit direction from donor N to acceptor O
//   - midpoint: midpoint of N...O (the source point for the dipolar field)
//   - distance: |N...O| distance
//   - is_backbone: true (all DSSP H-bonds are backbone)
// ============================================================================

struct ResolvedHBond {
    size_t donor_N = SIZE_MAX;
    size_t acceptor_O = SIZE_MAX;
    size_t donor_residue = SIZE_MAX;
    size_t acceptor_residue = SIZE_MAX;
    Vec3 midpoint = Vec3::Zero();
    Vec3 h_hat = Vec3::Zero();        // donor N → acceptor O direction
    double distance = 0.0;            // N...O distance
    int sequence_separation = 0;
};


// ============================================================================
// The full H-bond dipolar tensor from one H-bond at one atom.
//
// Same derivation as McConnell with b_hat → h_hat:
//
//   M_ab = 9 cosθ d̂_a h_b  -  3 h_a h_b  -  (3 d̂_a d̂_b - δ_ab)
//
// Returns M_ab / r³ (Angstrom⁻³).
// ============================================================================

struct HBondKernelResult {
    Mat3 M_over_r3 = Mat3::Zero();
    double f = 0.0;
    double distance = 0.0;
};


static HBondKernelResult ComputeHBondKernel(
        const Vec3& atom_pos,
        const Vec3& hbond_midpoint,
        const Vec3& h_hat) {

    HBondKernelResult result;

    Vec3 d = atom_pos - hbond_midpoint;
    double r = d.norm();

    if (r < MIN_DISTANCE) return result;

    result.distance = r;
    double r3 = r * r * r;
    Vec3 d_hat = d / r;
    double cos_theta = d_hat.dot(h_hat);

    result.f = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.M_over_r3(a, b) =
                (9.0 * cos_theta * d_hat(a) * h_hat(b)
                 - 3.0 * h_hat(a) * h_hat(b)
                 - (3.0 * d_hat(a) * d_hat(b) - (a == b ? 1.0 : 0.0)))
                / r3;
        }
    }

    return result;
}


// ============================================================================
// HBondResult::Compute
// ============================================================================

std::unique_ptr<HBondResult> HBondResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("HBondResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& dssp = conf.Result<DsspResult>();
    // SpatialIndexResult is a declared dependency (ensures it's computed
    // before us) but we iterate H-bonds directly, not via spatial search.
    (void)conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_residues = protein.ResidueCount();

    auto result_ptr = std::make_unique<HBondResult>();
    result_ptr->conf_ = &conf;

    // ------------------------------------------------------------------
    // Step 1: Resolve DSSP H-bond partners to atom positions.
    //
    // DSSP provides up to 2 acceptor partners and 2 donor partners per
    // residue. Each partner is a residue index. We resolve to backbone
    // N (donor) and O (acceptor) atoms.
    //
    // Skip H-bonds where:
    //   - partner residue index is SIZE_MAX (no partner)
    //   - backbone N or O atoms are not present
    //   - sequence separation < SEQUENTIAL_EXCLUSION_THRESHOLD (too close)
    //   - N...O distance > HBOND_MAX_DIST
    // ------------------------------------------------------------------

    std::vector<ResolvedHBond> hbonds;

    // Use a set to deduplicate: (donor_N, acceptor_O) pairs
    std::set<std::pair<size_t, size_t>> seen;

    for (size_t ri = 0; ri < n_residues; ++ri) {
        const auto& dr = dssp.AllResidues()[ri];
        const Residue& res = protein.ResidueAt(ri);

        // This residue's N-H donates to acceptor residues
        for (int bi = 0; bi < 2; ++bi) {
            size_t acc_ri = dr.acceptors[bi].residue_index;
            if (acc_ri == SIZE_MAX || acc_ri >= n_residues) continue;

            const Residue& acc_res = protein.ResidueAt(acc_ri);
            if (res.N == Residue::NONE || acc_res.O == Residue::NONE) continue;

            int seq_sep = std::abs(static_cast<int>(ri) - static_cast<int>(acc_ri));
            if (seq_sep < SEQUENTIAL_EXCLUSION_THRESHOLD) continue;

            auto key = std::make_pair(res.N, acc_res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            Vec3 N_pos = conf.PositionAt(res.N);
            Vec3 O_pos = conf.PositionAt(acc_res.O);
            Vec3 d = O_pos - N_pos;
            double dist = d.norm();

            if (dist < MIN_DISTANCE || dist > HBOND_MAX_DIST) continue;

            ResolvedHBond hb;
            hb.donor_N = res.N;
            hb.acceptor_O = acc_res.O;
            hb.donor_residue = ri;
            hb.acceptor_residue = acc_ri;
            hb.midpoint = 0.5 * (N_pos + O_pos);
            hb.h_hat = d / dist;
            hb.distance = dist;
            hb.sequence_separation = seq_sep;
            hbonds.push_back(hb);
        }

        // This residue's C=O accepts from donor residues
        for (int bi = 0; bi < 2; ++bi) {
            size_t don_ri = dr.donors[bi].residue_index;
            if (don_ri == SIZE_MAX || don_ri >= n_residues) continue;

            const Residue& don_res = protein.ResidueAt(don_ri);
            if (don_res.N == Residue::NONE || res.O == Residue::NONE) continue;

            int seq_sep = std::abs(static_cast<int>(ri) - static_cast<int>(don_ri));
            if (seq_sep < SEQUENTIAL_EXCLUSION_THRESHOLD) continue;

            auto key = std::make_pair(don_res.N, res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            Vec3 N_pos = conf.PositionAt(don_res.N);
            Vec3 O_pos = conf.PositionAt(res.O);
            Vec3 d = O_pos - N_pos;
            double dist = d.norm();

            if (dist < MIN_DISTANCE || dist > HBOND_MAX_DIST) continue;

            ResolvedHBond hb;
            hb.donor_N = don_res.N;
            hb.acceptor_O = res.O;
            hb.donor_residue = don_ri;
            hb.acceptor_residue = ri;
            hb.midpoint = 0.5 * (N_pos + O_pos);
            hb.h_hat = d / dist;
            hb.distance = dist;
            hb.sequence_separation = seq_sep;
            hbonds.push_back(hb);
        }
    }

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "resolved " + std::to_string(hbonds.size()) +
        " unique backbone H-bonds from DSSP");

    // Store resolved H-bond geometry for SampleAt grid queries
    for (const auto& hb : hbonds) {
        result_ptr->hbond_midpoints_.push_back(hb.midpoint);
        result_ptr->hbond_directions_.push_back(hb.h_hat);
        result_ptr->hbond_distances_.push_back(hb.distance);
    }

    if (hbonds.empty()) {
        return result_ptr;
    }

    // ------------------------------------------------------------------
    // Step 2: Build the filter set for H-bond kernel evaluations.
    //
    // SelfSourceFilter: atom cannot be a field point for an H-bond
    //   where it is the donor N or acceptor O.
    // DipolarNearFieldFilter: point-source model invalid when field
    //   point is inside the N...O source distribution.
    // ------------------------------------------------------------------

    KernelFilterSet filters;
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "filter set: " + filters.Describe());

    // ------------------------------------------------------------------
    // Step 3: For each atom, compute the dipolar tensor from all H-bonds
    // that pass the filter set. 1/r³ decay handles range naturally.
    // ------------------------------------------------------------------

    int total_pairs = 0;
    int filtered_out = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        // Residue index of this atom (for sequence separation)
        size_t ai_res = protein.AtomAt(ai).residue_index;

        Mat3 M_total = Mat3::Zero();
        double nearest_dist = NO_DATA_SENTINEL;
        size_t nearest_hb_idx = SIZE_MAX;
        int count_3_5 = 0;

        for (size_t hi = 0; hi < hbonds.size(); ++hi) {
            const auto& hb = hbonds[hi];

            HBondKernelResult kernel = ComputeHBondKernel(
                atom_pos, hb.midpoint, hb.h_hat);

            if (kernel.distance < MIN_DISTANCE) continue;

            // Build evaluation context from already-computed geometry
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = hb.distance;  // N...O distance
            ctx.atom_index = ai;
            ctx.source_atom_a = hb.donor_N;
            ctx.source_atom_b = hb.acceptor_O;

            // Sequence separation: min distance to either endpoint residue
            int sep_don = std::abs(static_cast<int>(ai_res)
                                 - static_cast<int>(hb.donor_residue));
            int sep_acc = std::abs(static_cast<int>(ai_res)
                                 - static_cast<int>(hb.acceptor_residue));
            ctx.sequence_separation = std::min(sep_don, sep_acc);

            if (!filters.AcceptAll(ctx)) {
                filtered_out++;
                continue;
            }

            // Count H-bonds within 3.5A of this atom
            if (kernel.distance < HBOND_COUNT_RADIUS) count_3_5++;

            // Track nearest (among accepted evaluations only)
            if (kernel.distance < nearest_dist) {
                nearest_dist = kernel.distance;
                nearest_hb_idx = hi;
            }

            // Accumulate tensor from all H-bonds (1/r³ decay handles range)
            M_total += kernel.M_over_r3;
            total_pairs++;
        }

        ca.hbond_count_within_3_5A = count_3_5;

        if (nearest_hb_idx != SIZE_MAX) {
            const auto& nearest_hb = hbonds[nearest_hb_idx];
            ca.hbond_nearest_dist = nearest_dist;
            ca.hbond_nearest_dir = (atom_pos - nearest_hb.midpoint).normalized();
            ca.hbond_is_backbone = true;  // all DSSP H-bonds are backbone

            HBondKernelResult nearest_kernel = ComputeHBondKernel(
                atom_pos, nearest_hb.midpoint, nearest_hb.h_hat);

            ca.hbond_nearest_tensor = nearest_kernel.M_over_r3;
            ca.hbond_nearest_spherical = SphericalTensor::Decompose(
                nearest_kernel.M_over_r3);
            ca.hbond_inv_d3 = 1.0 / (nearest_dist * nearest_dist * nearest_dist);
        }

        // Check if this atom is a donor or acceptor
        for (const auto& hb : hbonds) {
            if (ai == hb.donor_N) ca.hbond_is_donor = true;
            if (ai == hb.acceptor_O) ca.hbond_is_acceptor = true;
        }

        ca.hbond_shielding_contribution = SphericalTensor::Decompose(M_total);
    }

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "atom_hbond_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " hbonds=" + std::to_string(hbonds.size()) +
        " atoms=" + std::to_string(n_atoms));

    return result_ptr;
}


SphericalTensor HBondResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_ || hbond_midpoints_.empty()) return SphericalTensor{};

    Mat3 M_total = Mat3::Zero();

    for (size_t hi = 0; hi < hbond_midpoints_.size(); ++hi) {
        auto kernel = ComputeHBondKernel(
            point, hbond_midpoints_[hi], hbond_directions_[hi]);
        if (kernel.distance < MIN_DISTANCE) continue;

        // DipolarNearFieldFilter: skip if inside the N...O distribution
        if (kernel.distance < 0.5 * hbond_distances_[hi]) continue;

        M_total += kernel.M_over_r3;
    }

    return SphericalTensor::Decompose(M_total);
}


static void PackST_HB(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int HBondResult::WriteFeatures(const ProteinConformation& conf,
                                const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> scalars(N * 3);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_HB(ca.hbond_shielding_contribution, &shielding[i*9]);
        scalars[i*3+0] = ca.hbond_nearest_dist;
        scalars[i*3+1] = ca.hbond_inv_d3;
        scalars[i*3+2] = static_cast<double>(ca.hbond_count_within_3_5A);
    }

    NpyWriter::WriteFloat64(output_dir + "/hbond_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/hbond_scalars.npy", scalars.data(), N, 3);
    return 2;
}

}  // namespace nmr
