#include "MopacMcConnellResult.h"
#include "Protein.h"
#include "MopacResult.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <algorithm>

namespace nmr {


std::vector<std::type_index> MopacMcConnellResult::Dependencies() const {
    return {
        std::type_index(typeid(MopacResult)),
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// Bond kernel computation (same physics as McConnellResult).
//
// M_ab = 9 cos_theta d_hat_a b_hat_b
//      - 3 b_hat_a b_hat_b
//      - (3 d_hat_a d_hat_b - delta_ab)
//
// Returns M_ab / r^3 (Angstrom^-3), plus symmetric traceless K and scalar f.
// ============================================================================

struct MopacBondKernelResult {
    Mat3 M_over_r3 = Mat3::Zero();
    Mat3 K = Mat3::Zero();
    double f = 0.0;
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();
};


static MopacBondKernelResult ComputeBondKernel(
        const Vec3& atom_pos,
        const Vec3& bond_midpoint,
        const Vec3& bond_direction) {

    MopacBondKernelResult result;

    Vec3 d = atom_pos - bond_midpoint;
    double r = d.norm();

    if (r < MIN_DISTANCE) return result;

    result.distance = r;

    double r3 = r * r * r;
    Vec3 d_hat = d / r;
    result.direction = d_hat;

    double cos_theta = d_hat.dot(bond_direction);

    // McConnell scalar: (3 cos^2 theta - 1) / r^3
    result.f = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    // Symmetric traceless dipolar kernel K_ab
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.K(a, b) = (3.0 * d_hat(a) * d_hat(b)
                              - (a == b ? 1.0 : 0.0)) / r3;

    // Full McConnell tensor M_ab / r^3
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            result.M_over_r3(a, b) =
                (9.0 * cos_theta * d_hat(a) * bond_direction(b)
                 - 3.0 * bond_direction(a) * bond_direction(b)
                 - (3.0 * d_hat(a) * d_hat(b) - (a == b ? 1.0 : 0.0)))
                / r3;
        }
    }

    return result;
}


// ============================================================================
// MopacMcConnellResult::Compute
//
// Same loop as McConnellResult but each bond's contribution is weighted
// by its MOPAC Wiberg bond order. Bonds with no MOPAC order (0.0)
// contribute nothing — they are electronically insignificant.
// ============================================================================

std::unique_ptr<MopacMcConnellResult> MopacMcConnellResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("MopacMcConnellResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " bonds=" + std::to_string(conf.ProteinRef().BondCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const auto& mopac = conf.Result<MopacResult>();
    const size_t n_atoms = conf.AtomCount();

    auto result_ptr = std::make_unique<MopacMcConnellResult>();
    result_ptr->conf_ = &conf;

    KernelFilterSet filters;
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    int total_pairs = 0;
    int filtered_out = 0;
    int zero_bo_skipped = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_bonds = spatial.BondsWithinRadius(atom_pos, MOPAC_MCCONNELL_CUTOFF_A);

        // Per-category accumulators (bond-order-weighted)
        double co_sum = 0.0, cn_sum = 0.0, sidechain_sum = 0.0, aromatic_sum = 0.0;
        Mat3 M_backbone_total = Mat3::Zero();
        Mat3 M_sidechain_total = Mat3::Zero();
        Mat3 M_aromatic_total = Mat3::Zero();
        Mat3 M_total = Mat3::Zero();

        // Nearest CO and CN tracking (bond-order-weighted scalar)
        double best_co_dist = NO_DATA_SENTINEL;
        double best_cn_dist = NO_DATA_SENTINEL;
        double best_co_f_weighted = 0.0;
        MopacBondKernelResult best_co_kernel;
        MopacBondKernelResult best_cn_kernel;
        double best_co_bo = 0.0;
        double best_cn_bo = 0.0;

        for (size_t bi : nearby_bonds) {
            const Bond& bond = protein.BondAt(bi);

            // MOPAC Wiberg bond order for this topology bond
            double bo = mopac.TopologyBondOrder(bi);
            if (bo < 1e-6) { zero_bo_skipped++; continue; }

            Vec3 midpoint = conf.bond_midpoints[bi];
            Vec3 direction = conf.bond_directions[bi];

            MopacBondKernelResult kernel = ComputeBondKernel(atom_pos, midpoint, direction);
            if (kernel.distance < MIN_DISTANCE) continue;

            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = conf.bond_lengths[bi];
            ctx.atom_index = ai;
            ctx.source_atom_a = bond.atom_index_a;
            ctx.source_atom_b = bond.atom_index_b;
            if (!filters.AcceptAll(ctx)) { filtered_out++; continue; }

            // Bond-order-weighted accumulation
            Mat3 weighted_M = bo * kernel.M_over_r3;
            double weighted_f = bo * kernel.f;

            M_total += weighted_M;

            switch (bond.category) {
                case BondCategory::PeptideCO:
                    co_sum += weighted_f;
                    M_backbone_total += weighted_M;
                    if (kernel.distance < best_co_dist) {
                        best_co_dist = kernel.distance;
                        best_co_f_weighted = weighted_f;
                        best_co_kernel = kernel;
                        best_co_bo = bo;
                    }
                    break;

                case BondCategory::PeptideCN:
                    cn_sum += weighted_f;
                    M_backbone_total += weighted_M;
                    if (kernel.distance < best_cn_dist) {
                        best_cn_dist = kernel.distance;
                        best_cn_kernel = kernel;
                        best_cn_bo = bo;
                    }
                    break;

                case BondCategory::BackboneOther:
                    M_backbone_total += weighted_M;
                    break;

                case BondCategory::SidechainCO:
                    sidechain_sum += weighted_f;
                    M_sidechain_total += weighted_M;
                    break;

                case BondCategory::Aromatic:
                    aromatic_sum += weighted_f;
                    M_aromatic_total += weighted_M;
                    break;

                case BondCategory::SidechainOther:
                    M_sidechain_total += weighted_M;
                    break;

                default:
                    break;
            }

            total_pairs++;
        }

        // Store per-atom totals
        ca.mopac_mc_co_sum = co_sum;
        ca.mopac_mc_cn_sum = cn_sum;
        ca.mopac_mc_sidechain_sum = sidechain_sum;
        ca.mopac_mc_aromatic_sum = aromatic_sum;

        ca.mopac_mc_co_nearest = best_co_f_weighted;
        ca.mopac_mc_nearest_CO_dist = best_co_dist;
        ca.mopac_mc_nearest_CN_dist = best_cn_dist;

        if (best_co_dist < NO_DATA_SENTINEL) {
            ca.mopac_mc_T2_CO_nearest =
                SphericalTensor::Decompose(best_co_bo * best_co_kernel.K);
        }
        if (best_cn_dist < NO_DATA_SENTINEL) {
            ca.mopac_mc_T2_CN_nearest =
                SphericalTensor::Decompose(best_cn_bo * best_cn_kernel.K);
        }

        // Category T2 totals — extract symmetric traceless part
        Mat3 K_backbone = 0.5 * (M_backbone_total + M_backbone_total.transpose());
        K_backbone -= (K_backbone.trace() / 3.0) * Mat3::Identity();
        ca.mopac_mc_T2_backbone_total = SphericalTensor::Decompose(K_backbone);

        Mat3 K_sidechain = 0.5 * (M_sidechain_total + M_sidechain_total.transpose());
        K_sidechain -= (K_sidechain.trace() / 3.0) * Mat3::Identity();
        ca.mopac_mc_T2_sidechain_total = SphericalTensor::Decompose(K_sidechain);

        Mat3 K_aromatic = 0.5 * (M_aromatic_total + M_aromatic_total.transpose());
        K_aromatic -= (K_aromatic.trace() / 3.0) * Mat3::Identity();
        ca.mopac_mc_T2_aromatic_total = SphericalTensor::Decompose(K_aromatic);

        // Full shielding contribution
        ca.mopac_mc_shielding_contribution = SphericalTensor::Decompose(M_total);
    }

    OperationLog::Info(LogCalcMcConnell, "MopacMcConnellResult::Compute",
        "atom_bond_pairs=" + std::to_string(total_pairs) +
        " zero_bo_skipped=" + std::to_string(zero_bo_skipped) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms));

    return result_ptr;
}


// ============================================================================
// WriteFeatures: mopac_mc_shielding (9),
// mopac_mc_category_T2 (5 categories × 5 T2),
// mopac_mc_scalars (weighted CO/CN/sidechain/aromatic sums, nearest dists).
// ============================================================================

static void PackST_MMC(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int MopacMcConnellResult::WriteFeatures(const ProteinConformation& conf,
                                         const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> cat_T2(N * 25);
    std::vector<double> scalars(N * 6);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_MMC(ca.mopac_mc_shielding_contribution, &shielding[i*9]);

        const SphericalTensor* cats[5] = {
            &ca.mopac_mc_T2_backbone_total, &ca.mopac_mc_T2_sidechain_total,
            &ca.mopac_mc_T2_aromatic_total, &ca.mopac_mc_T2_CO_nearest,
            &ca.mopac_mc_T2_CN_nearest
        };
        for (int c = 0; c < 5; ++c)
            for (int m = 0; m < 5; ++m)
                cat_T2[i*25 + c*5 + m] = cats[c]->T2[m];

        scalars[i*6+0] = ca.mopac_mc_co_sum;
        scalars[i*6+1] = ca.mopac_mc_cn_sum;
        scalars[i*6+2] = ca.mopac_mc_sidechain_sum;
        scalars[i*6+3] = ca.mopac_mc_aromatic_sum;
        scalars[i*6+4] = ca.mopac_mc_nearest_CO_dist;
        scalars[i*6+5] = ca.mopac_mc_nearest_CN_dist;
    }

    NpyWriter::WriteFloat64(output_dir + "/mopac_mc_shielding.npy",
                            shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mopac_mc_category_T2.npy",
                            cat_T2.data(), N, 25);
    NpyWriter::WriteFloat64(output_dir + "/mopac_mc_scalars.npy",
                            scalars.data(), N, 6);
    return 3;
}

}  // namespace nmr
