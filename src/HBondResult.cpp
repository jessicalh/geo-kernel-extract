#include "HBondResult.h"
#include "Protein.h"
#include "DsspResult.h"
#include "SpatialIndexResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
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


struct ResolvedHBond {
    size_t donor_N = SIZE_MAX;
    size_t acceptor_O = SIZE_MAX;
    size_t donor_residue = SIZE_MAX;
    size_t acceptor_residue = SIZE_MAX;
    Vec3 midpoint = Vec3::Zero();
    Vec3 h_hat = Vec3::Zero();        // donor N to acceptor O direction
    double distance = 0.0;            // N...O distance
    int sequence_separation = 0;
};


// H-bond dipolar tensor from one H-bond at one atom:
//
//   M_ab = 9 cos(theta) d_hat_a h_b - 3 h_a h_b
//          - (3 d_hat_a d_hat_b - delta_ab)
//
// Returns M_ab / r^3 in A^-3.
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

    if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;

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

    GeometryChoiceBuilder choices(conf);

    std::vector<ResolvedHBond> hbonds;

    std::set<std::pair<size_t, size_t>> seen;
    size_t choice_record_seq = 0;

    for (size_t ri = 0; ri < n_residues; ++ri) {
        const auto& dr = dssp.AllResidues()[ri];
        const Residue& res = protein.ResidueAt(ri);

        for (int partner_slot = 0; partner_slot < 2; ++partner_slot) {
            size_t acceptor_residue_idx = dr.acceptors[partner_slot].residue_index;
            if (acceptor_residue_idx == SIZE_MAX || acceptor_residue_idx >= n_residues) continue;

            const Residue& acc_res = protein.ResidueAt(acceptor_residue_idx);
            if (res.N == Residue::NONE || acc_res.O == Residue::NONE) continue;

            int seq_sep = std::abs(static_cast<int>(ri) - static_cast<int>(acceptor_residue_idx));
            if (seq_sep < static_cast<int>(CalculatorConfig::Get("hbond_sequential_exclusion_residues"))) {
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [ri, acceptor_residue_idx, seq_sep](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(acceptor_residue_idx), "index");
                        AddNumber(gc, "sequence_separation", static_cast<double>(seq_sep), "residues");
                        AddNumber(gc, "rejection", 1.0, "seq_sep_too_small");
                    });
                continue;
            }

            auto key = std::make_pair(res.N, acc_res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            Vec3 N_pos = conf.PositionAt(res.N);
            Vec3 O_pos = conf.PositionAt(acc_res.O);
            Vec3 d = O_pos - N_pos;
            double dist = d.norm();

            if (dist < CalculatorConfig::Get("singularity_guard_distance") || dist > CalculatorConfig::Get("hbond_dipolar_max_distance")) {
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [ri, acceptor_residue_idx, dist](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(acceptor_residue_idx), "index");
                        AddNumber(gc, "distance", dist, "A");
                        AddNumber(gc, "rejection", 1.0, "distance_out_of_range");
                    });
                continue;
            }

            ResolvedHBond hb;
            hb.donor_N = res.N;
            hb.acceptor_O = acc_res.O;
            hb.donor_residue = ri;
            hb.acceptor_residue = acceptor_residue_idx;
            hb.midpoint = 0.5 * (N_pos + O_pos);
            hb.h_hat = d / dist;
            hb.distance = dist;
            hb.sequence_separation = seq_sep;
            hbonds.push_back(hb);
        }

        for (int partner_slot = 0; partner_slot < 2; ++partner_slot) {
            size_t donor_residue_idx = dr.donors[partner_slot].residue_index;
            if (donor_residue_idx == SIZE_MAX || donor_residue_idx >= n_residues) continue;

            const Residue& don_res = protein.ResidueAt(donor_residue_idx);
            if (don_res.N == Residue::NONE || res.O == Residue::NONE) continue;

            int seq_sep = std::abs(static_cast<int>(ri) - static_cast<int>(donor_residue_idx));
            if (seq_sep < static_cast<int>(CalculatorConfig::Get("hbond_sequential_exclusion_residues"))) {
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [ri, donor_residue_idx, seq_sep](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(donor_residue_idx), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "sequence_separation", static_cast<double>(seq_sep), "residues");
                        AddNumber(gc, "rejection", 1.0, "seq_sep_too_small");
                    });
                continue;
            }

            auto key = std::make_pair(don_res.N, res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            Vec3 N_pos = conf.PositionAt(don_res.N);
            Vec3 O_pos = conf.PositionAt(res.O);
            Vec3 d = O_pos - N_pos;
            double dist = d.norm();

            if (dist < CalculatorConfig::Get("singularity_guard_distance") || dist > CalculatorConfig::Get("hbond_dipolar_max_distance")) {
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [donor_residue_idx, ri, dist](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(donor_residue_idx), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "distance", dist, "A");
                        AddNumber(gc, "rejection", 1.0, "distance_out_of_range");
                    });
                continue;
            }

            ResolvedHBond hb;
            hb.donor_N = don_res.N;
            hb.acceptor_O = res.O;
            hb.donor_residue = donor_residue_idx;
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

    for (const auto& hb : hbonds) {
        result_ptr->hbond_midpoints_.push_back(hb.midpoint);
        result_ptr->hbond_directions_.push_back(hb.h_hat);
        result_ptr->hbond_distances_.push_back(hb.distance);
    }

    if (hbonds.empty()) {
        return result_ptr;
    }

    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "filter set: " + filters.Describe());

    int total_pairs = 0;
    int filtered_out = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        size_t ai_res = protein.AtomAt(ai).residue_index;

        Mat3 M_total = Mat3::Zero();
        double f_sum = 0.0;
        double nearest_dist = NO_DATA_SENTINEL;
        size_t nearest_hb_idx = SIZE_MAX;
        int nearby_hbond_count = 0;

        for (size_t hi = 0; hi < hbonds.size(); ++hi) {
            const auto& hb = hbonds[hi];

            HBondKernelResult kernel = ComputeHBondKernel(
                atom_pos, hb.midpoint, hb.h_hat);

            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = hb.distance;  // N...O distance
            ctx.atom_index = ai;
            ctx.source_atom_a = hb.donor_N;
            ctx.source_atom_b = hb.acceptor_O;

            int sep_don = std::abs(static_cast<int>(ai_res)
                                 - static_cast<int>(hb.donor_residue));
            int sep_acc = std::abs(static_cast<int>(ai_res)
                                 - static_cast<int>(hb.acceptor_residue));
            ctx.sequence_separation = std::min(sep_don, sep_acc);

            if (!filters.AcceptAll(ctx)) {
                choices.Record(CalculatorId::HBond, hi, "filter exclusion",
                    [&](GeometryChoice& gc) {
                        AddAtom(gc, &conf.AtomAt(hb.donor_N), hb.donor_N,
                                EntityRole::Context, EntityOutcome::Included);
                        AddAtom(gc, &conf.AtomAt(hb.acceptor_O), hb.acceptor_O,
                                EntityRole::Context, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                filtered_out++;
                continue;
            }

            if (kernel.distance < CalculatorConfig::Get("hbond_counting_radius")) nearby_hbond_count++;

            if (kernel.distance < nearest_dist) {
                nearest_dist = kernel.distance;
                nearest_hb_idx = hi;
            }

            M_total += kernel.M_over_r3;
            f_sum += kernel.f;
            total_pairs++;
        }

        ca.hbond_count_within_3_5A = nearby_hbond_count;

        choices.Record(CalculatorId::HBond, ai, "hbond neighbourhood",
            [&ca, ai, nearby_hbond_count](GeometryChoice& gc) {
                AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Included);
                AddNumber(gc, "hbond_count_within_3_5A", static_cast<double>(nearby_hbond_count), "count");
            });

        if (nearest_hb_idx != SIZE_MAX) {
            const auto& nearest_hb = hbonds[nearest_hb_idx];
            ca.hbond_nearest_dist = nearest_dist;
            ca.hbond_nearest_dir = (atom_pos - nearest_hb.midpoint).normalized();
            ca.hbond_is_backbone = true;

            // recompute kernel for the nearest H-bond only (inner loop did not
            // retain per-bond kernels)
            HBondKernelResult nearest_kernel = ComputeHBondKernel(
                atom_pos, nearest_hb.midpoint, nearest_hb.h_hat);

            ca.hbond_nearest_tensor = nearest_kernel.M_over_r3;
            ca.hbond_nearest_spherical = SphericalTensor::Decompose(
                nearest_kernel.M_over_r3);
            // Atom-to-midpoint r, distinct from the N...O source extent.
            ca.hbond_inv_d3 = 1.0 / (nearest_dist * nearest_dist * nearest_dist);
        }

        for (const auto& hb : hbonds) {
            if (ai == hb.donor_N) ca.hbond_is_donor = true;
            if (ai == hb.acceptor_O) ca.hbond_is_acceptor = true;
        }

        ca.hbond_shielding_contribution = SphericalTensor::Decompose(M_total);
        ca.hbond_mcconnell_scalar = f_sum;
    }

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "atom_hbond_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " hbonds=" + std::to_string(hbonds.size()) +
        " atoms=" + std::to_string(n_atoms));

    return result_ptr;
}


SphericalTensor HBondResult::SampleKernelAt(Vec3 point) const {
    if (!conf_ || hbond_midpoints_.empty()) return SphericalTensor{};

    Mat3 M_total = Mat3::Zero();

    for (size_t hi = 0; hi < hbond_midpoints_.size(); ++hi) {
        auto kernel = ComputeHBondKernel(
            point, hbond_midpoints_[hi], hbond_directions_[hi]);
        // The singularity guard leaves distance at zero.
        if (kernel.distance < CalculatorConfig::Get("singularity_guard_distance")) continue;

        if (kernel.distance < CalculatorConfig::Get("near_field_exclusion_ratio") * hbond_distances_[hi]) continue;

        M_total += kernel.M_over_r3;
    }

    return SphericalTensor::Decompose(M_total);
}


int HBondResult::WriteFeatures(const ProteinConformation& conf,
                                const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> scalars(N * 4);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        ca.hbond_shielding_contribution.PackFull9(&shielding[i*9]);
        scalars[i*4+0] = ca.hbond_nearest_dist;
        scalars[i*4+1] = ca.hbond_inv_d3;
        scalars[i*4+2] = static_cast<double>(ca.hbond_count_within_3_5A);
        scalars[i*4+3] = ca.hbond_mcconnell_scalar;
    }

    NpyWriter::WriteFloat64(output_dir + "/hbond_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/hbond_scalars.npy", scalars.data(), N, 4);
    return 2;
}

}  // namespace nmr
