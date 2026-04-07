#include "McConnellResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <algorithm>

namespace nmr {


std::vector<std::type_index> McConnellResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// The full McConnell shielding tensor from one bond at one atom.
//
// M_ab = 9 cos_theta d_hat_a b_hat_b
//      - 3 b_hat_a b_hat_b
//      - (3 d_hat_a d_hat_b - delta_ab)
//
// Returns M_ab / r^3 (Angstrom^-3).
//
// Also computes:
//   K_ab = (3 d_hat_a d_hat_b - delta_ab) / r^3  (symmetric traceless)
//   f    = (3 cos^2 theta - 1) / r^3              (McConnell scalar)
//
// Three terms in M:
//   Term 1: 9 cos_theta d_hat ⊗ b_hat     — asymmetric, gives T1
//   Term 2: -3 b_hat ⊗ b_hat              — symmetric, gives T0
//   Term 3: -(3 d_hat ⊗ d_hat - I)        — symmetric traceless, gives T2
// ============================================================================

struct BondKernelResult {
    Mat3 M_over_r3 = Mat3::Zero();   // full McConnell tensor (asymmetric)
    Mat3 K = Mat3::Zero();           // symmetric traceless dipolar kernel
    double f = 0.0;                  // McConnell scalar
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();   // unit vector from bond midpoint to atom
};


static BondKernelResult ComputeBondKernel(
        const Vec3& atom_pos,
        const Vec3& bond_midpoint,
        const Vec3& bond_direction) {

    BondKernelResult result;

    Vec3 d = atom_pos - bond_midpoint;
    double r = d.norm();

    if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;

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
    //   = [9 cos_theta d_hat_a b_hat_b - 3 b_hat_a b_hat_b
    //      - (3 d_hat_a d_hat_b - delta_ab)] / r^3
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
// McConnellResult::Compute
// ============================================================================

std::unique_ptr<McConnellResult> McConnellResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("McConnellResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " bonds=" + std::to_string(conf.ProteinRef().BondCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_bonds = protein.BondCount();

    auto result_ptr = std::make_unique<McConnellResult>();
    result_ptr->conf_ = &conf;

    // Filter set: SelfSourceFilter (atom is bond endpoint) +
    // DipolarNearFieldFilter (source extent = bond length).
    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcMcConnell, "McConnellResult::Compute",
        "filter set: " + filters.Describe());

    GeometryChoiceBuilder choices(conf);

    int total_pairs = 0;
    int filtered_out = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        // Find nearby bonds via spatial index
        auto nearby_bonds = spatial.BondsWithinRadius(atom_pos, CalculatorConfig::Get("mcconnell_bond_anisotropy_cutoff"));

        // Per-category accumulators
        double co_sum = 0.0, cn_sum = 0.0, sidechain_sum = 0.0, aromatic_sum = 0.0;
        Mat3 M_backbone_total = Mat3::Zero();
        Mat3 M_sidechain_total = Mat3::Zero();
        Mat3 M_aromatic_total = Mat3::Zero();
        Mat3 M_total = Mat3::Zero();

        // Nearest CO and CN tracking
        double best_co_dist = NO_DATA_SENTINEL;
        double best_cn_dist = NO_DATA_SENTINEL;
        double best_co_f = 0.0;
        Vec3 best_co_midpoint = Vec3::Zero();
        Vec3 best_co_direction = Vec3::Zero();
        BondKernelResult best_co_kernel;
        BondKernelResult best_cn_kernel;

        for (size_t bi : nearby_bonds) {
            const Bond& bond = protein.BondAt(bi);
            Vec3 midpoint = conf.bond_midpoints[bi];
            Vec3 direction = conf.bond_directions[bi];

            BondKernelResult kernel = ComputeBondKernel(atom_pos, midpoint, direction);

            // Build evaluation context and apply filter set
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = conf.bond_lengths[bi];
            ctx.atom_index = ai;
            ctx.source_atom_a = bond.atom_index_a;
            ctx.source_atom_b = bond.atom_index_b;
            if (!filters.AcceptAll(ctx)) {
                // ---- GeometryChoice: filter exclusion ----
                choices.Record(CalculatorId::McConnell, bi, "filter exclusion",
                    [&](GeometryChoice& gc) {
                        AddBond(gc, &bond, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                filtered_out++; continue;
            }

            // Store in BondNeighbourhood
            BondNeighbourhood bn;
            bn.bond_index = bi;
            bn.bond_category = bond.category;
            bn.distance_to_midpoint = kernel.distance;
            bn.direction_to_midpoint = kernel.direction;
            bn.dipolar_tensor = kernel.K;
            bn.dipolar_spherical = SphericalTensor::Decompose(kernel.K);
            bn.mcconnell_scalar = kernel.f;
            ca.bond_neighbours.push_back(bn);

            // Accumulate full McConnell tensor by category
            M_total += kernel.M_over_r3;

            switch (bond.category) {
                case BondCategory::PeptideCO:
                    co_sum += kernel.f;
                    M_backbone_total += kernel.M_over_r3;
                    if (kernel.distance < best_co_dist) {
                        best_co_dist = kernel.distance;
                        best_co_f = kernel.f;
                        best_co_midpoint = midpoint;
                        best_co_direction = kernel.direction;
                        best_co_kernel = kernel;
                    }
                    break;

                case BondCategory::PeptideCN:
                    cn_sum += kernel.f;
                    M_backbone_total += kernel.M_over_r3;
                    if (kernel.distance < best_cn_dist) {
                        best_cn_dist = kernel.distance;
                        best_cn_kernel = kernel;
                    }
                    break;

                case BondCategory::BackboneOther:
                    M_backbone_total += kernel.M_over_r3;
                    break;

                case BondCategory::SidechainCO:
                    sidechain_sum += kernel.f;
                    M_sidechain_total += kernel.M_over_r3;
                    break;

                case BondCategory::Aromatic:
                    aromatic_sum += kernel.f;
                    M_aromatic_total += kernel.M_over_r3;
                    break;

                case BondCategory::SidechainOther:
                    M_sidechain_total += kernel.M_over_r3;
                    break;

                default:
                    break;
            }

            total_pairs++;
        }

        // Store per-atom totals
        ca.mcconnell_co_sum = co_sum;
        ca.mcconnell_cn_sum = cn_sum;
        ca.mcconnell_sidechain_sum = sidechain_sum;
        ca.mcconnell_aromatic_sum = aromatic_sum;

        // Nearest CO
        ca.mcconnell_co_nearest = best_co_f;
        ca.nearest_CO_midpoint = best_co_midpoint;
        ca.nearest_CO_dist = best_co_dist;
        ca.nearest_CN_dist = best_cn_dist;

        if (best_co_dist < NO_DATA_SENTINEL) {
            double dir_norm = best_co_direction.norm();
            ca.dir_nearest_CO = (dir_norm > CalculatorConfig::Get("near_zero_vector_norm_threshold"))
                ? Vec3(best_co_direction / dir_norm) : Vec3::Zero();
            ca.T2_CO_nearest = SphericalTensor::Decompose(best_co_kernel.K);
        }
        if (best_cn_dist < NO_DATA_SENTINEL) {
            ca.T2_CN_nearest = SphericalTensor::Decompose(best_cn_kernel.K);
        }

        // Category T2 totals (from symmetric dipolar kernel sums)
        // Apply traceless projection to fix floating-point drift
        auto project_traceless = [](Mat3& m) {
            double trace = m.trace();
            m -= (trace / 3.0) * Mat3::Identity();
        };

        // Extract symmetric part for T2 features
        Mat3 K_backbone = 0.5 * (M_backbone_total + M_backbone_total.transpose());
        K_backbone -= (K_backbone.trace() / 3.0) * Mat3::Identity();
        ca.T2_backbone_total = SphericalTensor::Decompose(K_backbone);

        Mat3 K_sidechain = 0.5 * (M_sidechain_total + M_sidechain_total.transpose());
        K_sidechain -= (K_sidechain.trace() / 3.0) * Mat3::Identity();
        ca.T2_sidechain_total = SphericalTensor::Decompose(K_sidechain);

        Mat3 K_aromatic = 0.5 * (M_aromatic_total + M_aromatic_total.transpose());
        K_aromatic -= (K_aromatic.trace() / 3.0) * Mat3::Identity();
        ca.T2_aromatic_total = SphericalTensor::Decompose(K_aromatic);

        // Full McConnell shielding contribution (from the full M tensor sum)
        ca.mc_shielding_contribution = SphericalTensor::Decompose(M_total);

        // ---- GeometryChoice: bond anisotropy ----
        choices.Record(CalculatorId::McConnell, ai, "bond anisotropy",
            [&ca, ai, best_co_dist, best_cn_dist](GeometryChoice& gc) {
                AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Included);
                AddNumber(gc, "nearest_CO_dist", best_co_dist, "A");
                AddNumber(gc, "nearest_CN_dist", best_cn_dist, "A");
            });
    }

    OperationLog::Info(LogCalcMcConnell, "McConnellResult::Compute",
        "atom_bond_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " bonds=" + std::to_string(n_bonds));

    return result_ptr;
}


// ============================================================================
// Query methods
// ============================================================================

double McConnellResult::CategorySum(size_t atom_index, BondCategory cat) const {
    const auto& ca = conf_->AtomAt(atom_index);
    switch (cat) {
        case BondCategory::PeptideCO: return ca.mcconnell_co_sum;
        case BondCategory::PeptideCN: return ca.mcconnell_cn_sum;
        case BondCategory::SidechainCO: return ca.mcconnell_sidechain_sum;
        case BondCategory::Aromatic: return ca.mcconnell_aromatic_sum;
        default: return 0.0;
    }
}

double McConnellResult::NearestCOContribution(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mcconnell_co_nearest;
}


// ============================================================================
// SampleShieldingAt: evaluate McConnell kernel at arbitrary 3D point.
// Sums over all bonds within MCCONNELL_CUTOFF_A. No atom-specific filters.
// ============================================================================

SphericalTensor McConnellResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 M_total = Mat3::Zero();

    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        auto kernel = ComputeBondKernel(
            point, conf_->bond_midpoints[bi], conf_->bond_directions[bi]);
        if (kernel.distance < CalculatorConfig::Get("singularity_guard_distance")) continue;
        if (kernel.distance > CalculatorConfig::Get("mcconnell_bond_anisotropy_cutoff")) continue;

        // DipolarNearFieldFilter: skip if inside the bond
        double bond_len = conf_->bond_lengths[bi];
        if (kernel.distance < CalculatorConfig::Get("near_field_exclusion_ratio") * bond_len) continue;

        M_total += kernel.M_over_r3;
    }

    return SphericalTensor::Decompose(M_total);
}


// ============================================================================
// WriteFeatures: mc_shielding (9), mc_category_T2 (5 categories × 5 T2),
// mc_scalars (CO/CN/sidechain/aromatic sums, nearest distances).
// ============================================================================

static void PackST_MC(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int McConnellResult::WriteFeatures(const ProteinConformation& conf,
                                    const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> cat_T2(N * 25);
    std::vector<double> scalars(N * 6);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_MC(ca.mc_shielding_contribution, &shielding[i*9]);

        const SphericalTensor* cats[5] = {
            &ca.T2_backbone_total, &ca.T2_sidechain_total,
            &ca.T2_aromatic_total, &ca.T2_CO_nearest, &ca.T2_CN_nearest
        };
        for (int c = 0; c < 5; ++c)
            for (int m = 0; m < 5; ++m)
                cat_T2[i*25 + c*5 + m] = cats[c]->T2[m];

        scalars[i*6+0] = ca.mcconnell_co_sum;
        scalars[i*6+1] = ca.mcconnell_cn_sum;
        scalars[i*6+2] = ca.mcconnell_sidechain_sum;
        scalars[i*6+3] = ca.mcconnell_aromatic_sum;
        scalars[i*6+4] = ca.nearest_CO_dist;
        scalars[i*6+5] = ca.nearest_CN_dist;
    }

    NpyWriter::WriteFloat64(output_dir + "/mc_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mc_category_T2.npy", cat_T2.data(), N, 25);
    NpyWriter::WriteFloat64(output_dir + "/mc_scalars.npy", scalars.data(), N, 6);
    return 3;
}

}  // namespace nmr
