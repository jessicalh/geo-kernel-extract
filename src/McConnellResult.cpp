#include "McConnellResult.h"

#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "MopacResult.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "SpatialIndexResult.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

namespace fs = std::filesystem;

namespace nmr {

namespace {

constexpr double kNearFieldAuditDistanceA = 3.0;

using McAccum = std::array<std::array<Mat3, kMcConnellChannelCount>,
                           kMcConnellSourceCategoryCount>;

std::size_t CatIndex(McConnellSourceCategory cat) {
    return static_cast<std::size_t>(cat);
}

std::size_t ChannelIndex(McConnellChannel channel) {
    return static_cast<std::size_t>(channel);
}

void ZeroAccum(McAccum& accum) {
    for (auto& by_channel : accum)
        for (auto& m : by_channel)
            m = Mat3::Zero();
}

bool IsStrictBackboneAtom(const Protein& protein, std::size_t atom_index) {
    if (atom_index >= protein.AtomCount()) return false;
    const Atom& atom = protein.AtomAt(atom_index);
    if (atom.residue_index >= protein.ResidueCount()) return false;
    const Residue& residue = protein.ResidueAt(atom.residue_index);
    return atom_index == residue.N ||
           atom_index == residue.CA ||
           atom_index == residue.C ||
           atom_index == residue.O;
}

bool IsXHBond(const Protein& protein, const Bond& bond) {
    const Atom& a = protein.AtomAt(bond.atom_index_a);
    const Atom& b = protein.AtomAt(bond.atom_index_b);
    auto is_x = [](Element e) {
        return e == Element::C || e == Element::N || e == Element::O;
    };
    return (a.element == Element::H && is_x(b.element)) ||
           (b.element == Element::H && is_x(a.element));
}

McConnellSourceCategory SourceCategoryForMcConnell(
        const Protein& protein,
        const Bond& bond) {
    if (bond.category == BondCategory::Aromatic)
        return McConnellSourceCategory::AromaticZeroed;
    if (bond.category == BondCategory::Disulfide)
        return McConnellSourceCategory::Disulfide;
    if (bond.category == BondCategory::PeptideCO)
        return McConnellSourceCategory::PeptideCO;
    if (bond.category == BondCategory::PeptideCN)
        return McConnellSourceCategory::PeptideCN;
    if (bond.category == BondCategory::SidechainCO)
        return McConnellSourceCategory::SidechainCO;

    const bool strict_backbone =
        IsStrictBackboneAtom(protein, bond.atom_index_a) &&
        IsStrictBackboneAtom(protein, bond.atom_index_b);

    if (strict_backbone)
        return McConnellSourceCategory::BackboneOther;
    return McConnellSourceCategory::SidechainOther;
}

SphericalTensor DecomposeAccum(const McAccum& accum,
                               McConnellSourceCategory cat,
                               McConnellChannel channel) {
    return SphericalTensor::Decompose(
        accum[CatIndex(cat)][ChannelIndex(channel)]);
}

Mat3 SumFixedCategories(const McAccum& accum,
                        std::initializer_list<McConnellSourceCategory> cats) {
    Mat3 total = Mat3::Zero();
    for (McConnellSourceCategory cat : cats)
        total += accum[CatIndex(cat)][ChannelIndex(McConnellChannel::Fixed)];
    return total;
}

Mat3 SumChannel(const McAccum& accum, McConnellChannel channel) {
    Mat3 total = Mat3::Zero();
    const std::size_t ch = ChannelIndex(channel);
    for (const auto& by_channel : accum)
        total += by_channel[ch];
    return total;
}

void MergeMcConnellManifest(const std::string& output_dir,
                            bool include_xh_sources) {
    const fs::path path = fs::path(output_dir) / "extraction_manifest.json";
    nlohmann::ordered_json j = nlohmann::ordered_json::object();

    if (fs::exists(path)) {
        std::ifstream in(path);
        if (in) {
            try {
                j = nlohmann::ordered_json::parse(in);
            } catch (const std::exception& e) {
                OperationLog::Warn(
                    "McConnellResult::WriteManifest",
                    "could not parse existing manifest " + path.string() +
                    ": " + e.what() + "; rewriting minimal manifest");
                j = nlohmann::ordered_json::object();
            }
        }
    }

    if (!j.is_object()) j = nlohmann::ordered_json::object();
    if (!j.contains("schema_version")) j["schema_version"] = "1.0";
    if (!j.contains("extractor")) j["extractor"] = "nmr_extract";
    j["feature_metadata"]["mcconnell"] =
        McConnellResult::FeatureMetadata(include_xh_sources);

    std::ofstream out(path);
    if (!out) {
        OperationLog::Error("McConnellResult::WriteManifest",
            "could not open " + path.string() + " for write");
        return;
    }
    out << j.dump(2, ' ', false,
                  nlohmann::ordered_json::error_handler_t::replace) << "\n";
}

}  // namespace


// static
nlohmann::ordered_json McConnellResult::FeatureMetadata(
        bool include_xh_sources) {
    nlohmann::ordered_json arrays = nlohmann::ordered_json::array();
    nlohmann::ordered_json channels = nlohmann::ordered_json::array();
    nlohmann::ordered_json categories = nlohmann::ordered_json::array();

    for (std::size_t h = 0; h < kMcConnellChannelCount; ++h) {
        channels.push_back(McConnellChannelStem(
            static_cast<McConnellChannel>(h)));
    }
    for (std::size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
        const auto cat = static_cast<McConnellSourceCategory>(c);
        categories.push_back(McConnellSourceCategoryStem(cat));
        for (std::size_t h = 0; h < kMcConnellChannelCount; ++h) {
            const auto channel = static_cast<McConnellChannel>(h);
            arrays.push_back(
                std::string("mc_") + McConnellSourceCategoryStem(cat) +
                "_" + McConnellChannelStem(channel) + ".npy");
        }
    }

    return nlohmann::ordered_json{
        {"source_model", "unit susceptibility shape; scale learned"},
        {"bo_source", "MOPAC Wiberg bond order"},
        {"aromatic_zeroed_when_ring_active", true},
        {"aromatic_zeroed_reason",
         "BS/HM always compute the aromatic ring-current; McConnell zeros aromatic to avoid the double-count."},
        {"irrep_layout", kMcConnellPackFull9IrrepLayout},
        {"units", "Angstrom^-3"},
        {"rhombic_status", "absent_no_primary_table"},
        {"channels", channels},
        {"source_categories", categories},
        {"arrays", arrays},
        {"xh_sources", nlohmann::ordered_json{
            {"separable", true},
            {"toggle", "mcconnell_include_xh_sources"},
            {"included", include_xh_sources},
            {"definition", "C-H/N-H/O-H source bonds from ALLBONDS; scaffold for Part-1 ablation"}}},
    };
}


std::vector<std::type_index> McConnellResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


McConnellPairKernel McConnellResult::ComputePairKernel(
        const Vec3& atom_pos,
        const Vec3& source_center,
        const Vec3& source_axis) {
    McConnellPairKernel result;

    Vec3 disp = atom_pos - source_center;
    const double r = disp.norm();
    result.distance = r;
    if (r < CalculatorConfig::Get("singularity_guard_distance"))
        return result;

    const double axis_norm = source_axis.norm();
    if (axis_norm < 1e-15)
        return result;

    const double r3 = r * r * r;
    const Vec3 n = disp / r;
    const Vec3 u = source_axis / axis_norm;
    result.direction = n;

    result.dipolar = (3.0 * (n * n.transpose()) - Mat3::Identity()) / r3;
    result.source_shape = (u * u.transpose()) - (Mat3::Identity() / 3.0);
    result.response = result.dipolar * result.source_shape;
    result.pcs_scalar = result.response.trace() / 3.0;
    return result;
}


std::unique_ptr<McConnellResult> McConnellResult::Compute(
        ProteinConformation& conf) {
    OperationLog::Scope scope("McConnellResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " bonds=" + std::to_string(conf.ProteinRef().BondCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const MopacResult* mopac =
        conf.HasResult<MopacResult>() ? &conf.Result<MopacResult>() : nullptr;
    const bool include_xh_sources =
        CalculatorConfig::Get("mcconnell_include_xh_sources") != 0.0;
    const double cutoff =
        CalculatorConfig::Get("mcconnell_bond_anisotropy_cutoff");
    const double bo_floor =
        CalculatorConfig::Get("mopac_bond_order_noise_floor");
    const size_t n_atoms = conf.AtomCount();

    auto result_ptr = std::make_unique<McConnellResult>();
    result_ptr->conf_ = &conf;

    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcMcConnell, "McConnellResult::Compute",
        "source model: A=D(r)Qhat; categories=7 channels={fixed,bo}; "
        "aromatic_zeroed=unconditional; xh_included=" +
        std::string(include_xh_sources ? "true" : "false") +
        "; optional_mopac=" + std::string(mopac ? "present" : "absent") +
        "; filter set: " + filters.Describe());

    GeometryChoiceBuilder choices(conf);

    std::size_t total_pairs = 0;
    std::size_t filtered_out = 0;
    std::size_t xh_skipped = 0;
    std::size_t aromatic_zeroed_pairs = 0;
    std::size_t bo_zeroed_pairs = 0;
    std::size_t near_accepted_lt3 = 0;
    std::size_t near_rejected_lt3 = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        const Vec3 atom_pos = conf.PositionAt(ai);

        McAccum accum;
        ZeroAccum(accum);

        double best_co_dist = NO_DATA_SENTINEL;
        double best_cn_dist = NO_DATA_SENTINEL;
        double best_co_scalar = 0.0;
        Vec3 best_co_midpoint = Vec3::Zero();
        Vec3 best_co_direction = Vec3::Zero();
        McConnellPairKernel best_co_kernel;
        McConnellPairKernel best_cn_kernel;
        double best_co_bo = 0.0;
        double best_cn_bo = 0.0;
        std::size_t atom_near_accepted_lt3 = 0;
        std::size_t atom_near_rejected_lt3 = 0;

        const auto nearby_bonds = spatial.BondsWithinRadius(atom_pos, cutoff);
        for (size_t bi : nearby_bonds) {
            const Bond& bond = protein.BondAt(bi);
            if (!include_xh_sources && IsXHBond(protein, bond)) {
                ++xh_skipped;
                continue;
            }

            const Vec3 midpoint = conf.bond_midpoints[bi];
            const Vec3 direction = conf.bond_directions[bi];
            const McConnellPairKernel kernel =
                ComputePairKernel(atom_pos, midpoint, direction);

            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = conf.bond_lengths[bi];
            ctx.atom_index = ai;
            ctx.source_atom_a = bond.atom_index_a;
            ctx.source_atom_b = bond.atom_index_b;
            ctx.conf = &conf;

            if (!filters.AcceptAll(ctx)) {
                if (kernel.distance < kNearFieldAuditDistanceA) {
                    ++near_rejected_lt3;
                    ++atom_near_rejected_lt3;
                }
                choices.Record(CalculatorId::McConnell, bi,
                    "filter exclusion",
                    [&](GeometryChoice& gc) {
                        AddBond(gc, &bond, EntityRole::Source,
                                EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target,
                                EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                ++filtered_out;
                continue;
            }

            if (kernel.distance < kNearFieldAuditDistanceA) {
                ++near_accepted_lt3;
                ++atom_near_accepted_lt3;
            }

            BondNeighbourhood bn;
            bn.bond_index = bi;
            bn.bond_category = bond.category;
            bn.distance_to_midpoint = kernel.distance;
            bn.direction_to_midpoint = kernel.direction;
            bn.dipolar_tensor = kernel.dipolar;
            bn.dipolar_spherical = SphericalTensor::Decompose(kernel.dipolar);
            bn.mcconnell_scalar = kernel.pcs_scalar;
            ca.bond_neighbours.push_back(bn);

            const McConnellSourceCategory cat =
                SourceCategoryForMcConnell(protein, bond);

            const double raw_bo = mopac ? mopac->TopologyBondOrder(bi) : 0.0;
            const double bo = (raw_bo >= bo_floor) ? raw_bo : 0.0;
            if (bo == 0.0) ++bo_zeroed_pairs;

            if (cat == McConnellSourceCategory::AromaticZeroed) {
                ++aromatic_zeroed_pairs;
            } else {
                accum[CatIndex(cat)][ChannelIndex(McConnellChannel::Fixed)]
                    += kernel.response;
                accum[CatIndex(cat)][ChannelIndex(McConnellChannel::BondOrder)]
                    += bo * kernel.response;
            }

            if (cat == McConnellSourceCategory::PeptideCO &&
                kernel.distance < best_co_dist) {
                best_co_dist = kernel.distance;
                best_co_scalar = kernel.pcs_scalar;
                best_co_midpoint = midpoint;
                best_co_direction = kernel.direction;
                best_co_kernel = kernel;
                best_co_bo = bo;
            } else if (cat == McConnellSourceCategory::PeptideCN &&
                       kernel.distance < best_cn_dist) {
                best_cn_dist = kernel.distance;
                best_cn_kernel = kernel;
                best_cn_bo = bo;
            }

            ++total_pairs;
        }

        for (std::size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
            for (std::size_t h = 0; h < kMcConnellChannelCount; ++h) {
                ca.mcconnell_source_tensors[c][h] =
                    SphericalTensor::Decompose(accum[c][h]);
            }
        }

        ca.mcconnell_near_field_accepted_lt3A =
            static_cast<int>(atom_near_accepted_lt3);
        ca.mcconnell_near_field_rejected_lt3A =
            static_cast<int>(atom_near_rejected_lt3);

        const SphericalTensor fixed_co = DecomposeAccum(
            accum, McConnellSourceCategory::PeptideCO,
            McConnellChannel::Fixed);
        const SphericalTensor fixed_cn = DecomposeAccum(
            accum, McConnellSourceCategory::PeptideCN,
            McConnellChannel::Fixed);
        const SphericalTensor fixed_sidechain_co = DecomposeAccum(
            accum, McConnellSourceCategory::SidechainCO,
            McConnellChannel::Fixed);
        const SphericalTensor fixed_sidechain_other = DecomposeAccum(
            accum, McConnellSourceCategory::SidechainOther,
            McConnellChannel::Fixed);
        const SphericalTensor fixed_disulfide = DecomposeAccum(
            accum, McConnellSourceCategory::Disulfide,
            McConnellChannel::Fixed);

        ca.mcconnell_co_sum = fixed_co.T0;
        ca.mcconnell_cn_sum = fixed_cn.T0;
        ca.mcconnell_sidechain_sum =
            fixed_sidechain_co.T0 + fixed_sidechain_other.T0 + fixed_disulfide.T0;
        ca.mcconnell_aromatic_sum = 0.0;

        ca.mcconnell_co_nearest = best_co_scalar;
        ca.nearest_CO_midpoint = best_co_midpoint;
        ca.nearest_CO_dist = best_co_dist;
        ca.nearest_CN_dist = best_cn_dist;

        if (best_co_dist < NO_DATA_SENTINEL) {
            const double dir_norm = best_co_direction.norm();
            ca.dir_nearest_CO =
                (dir_norm > CalculatorConfig::Get("near_zero_vector_norm_threshold"))
                ? Vec3(best_co_direction / dir_norm) : Vec3::Zero();
            ca.T2_CO_nearest =
                SphericalTensor::Decompose(best_co_kernel.response);
            ca.mopac_mc_T2_CO_nearest =
                SphericalTensor::Decompose(best_co_bo * best_co_kernel.response);
        }
        if (best_cn_dist < NO_DATA_SENTINEL) {
            ca.T2_CN_nearest =
                SphericalTensor::Decompose(best_cn_kernel.response);
            ca.mopac_mc_T2_CN_nearest =
                SphericalTensor::Decompose(best_cn_bo * best_cn_kernel.response);
        }

        ca.T2_backbone_total = SphericalTensor::Decompose(SumFixedCategories(
            accum, {McConnellSourceCategory::PeptideCO,
                    McConnellSourceCategory::PeptideCN,
                    McConnellSourceCategory::BackboneOther}));
        ca.T2_sidechain_total = SphericalTensor::Decompose(SumFixedCategories(
            accum, {McConnellSourceCategory::SidechainCO,
                    McConnellSourceCategory::SidechainOther,
                    McConnellSourceCategory::Disulfide}));
        ca.T2_aromatic_total = SphericalTensor{};
        ca.mc_shielding_contribution =
            SphericalTensor::Decompose(SumChannel(accum, McConnellChannel::Fixed));

        const SphericalTensor bo_co = DecomposeAccum(
            accum, McConnellSourceCategory::PeptideCO,
            McConnellChannel::BondOrder);
        const SphericalTensor bo_cn = DecomposeAccum(
            accum, McConnellSourceCategory::PeptideCN,
            McConnellChannel::BondOrder);
        const SphericalTensor bo_sidechain_co = DecomposeAccum(
            accum, McConnellSourceCategory::SidechainCO,
            McConnellChannel::BondOrder);
        const SphericalTensor bo_sidechain_other = DecomposeAccum(
            accum, McConnellSourceCategory::SidechainOther,
            McConnellChannel::BondOrder);
        const SphericalTensor bo_disulfide = DecomposeAccum(
            accum, McConnellSourceCategory::Disulfide,
            McConnellChannel::BondOrder);

        ca.mopac_mc_co_sum = bo_co.T0;
        ca.mopac_mc_cn_sum = bo_cn.T0;
        ca.mopac_mc_sidechain_sum =
            bo_sidechain_co.T0 + bo_sidechain_other.T0 + bo_disulfide.T0;
        ca.mopac_mc_aromatic_sum = 0.0;
        ca.mopac_mc_co_nearest = best_co_bo * best_co_scalar;
        ca.mopac_mc_nearest_CO_dist = best_co_dist;
        ca.mopac_mc_nearest_CN_dist = best_cn_dist;
        ca.mopac_mc_T2_backbone_total = SphericalTensor::Decompose(
            accum[CatIndex(McConnellSourceCategory::PeptideCO)]
                 [ChannelIndex(McConnellChannel::BondOrder)] +
            accum[CatIndex(McConnellSourceCategory::PeptideCN)]
                 [ChannelIndex(McConnellChannel::BondOrder)] +
            accum[CatIndex(McConnellSourceCategory::BackboneOther)]
                 [ChannelIndex(McConnellChannel::BondOrder)]);
        ca.mopac_mc_T2_sidechain_total = SphericalTensor::Decompose(
            accum[CatIndex(McConnellSourceCategory::SidechainCO)]
                 [ChannelIndex(McConnellChannel::BondOrder)] +
            accum[CatIndex(McConnellSourceCategory::SidechainOther)]
                 [ChannelIndex(McConnellChannel::BondOrder)] +
            accum[CatIndex(McConnellSourceCategory::Disulfide)]
                 [ChannelIndex(McConnellChannel::BondOrder)]);
        ca.mopac_mc_T2_aromatic_total = SphericalTensor{};
        ca.mopac_mc_shielding_contribution =
            SphericalTensor::Decompose(SumChannel(
                accum, McConnellChannel::BondOrder));

        choices.Record(CalculatorId::McConnell, ai,
            "bond susceptibility anisotropy",
            [&ca, ai, best_co_dist, best_cn_dist](GeometryChoice& gc) {
                AddAtom(gc, &ca, ai, EntityRole::Target,
                        EntityOutcome::Included);
                AddNumber(gc, "nearest_CO_dist", best_co_dist, "A");
                AddNumber(gc, "nearest_CN_dist", best_cn_dist, "A");
                AddNumber(gc, "near_field_accepted_lt3A",
                          ca.mcconnell_near_field_accepted_lt3A, "count");
                AddNumber(gc, "near_field_rejected_lt3A",
                          ca.mcconnell_near_field_rejected_lt3A, "count");
            });
    }

    OperationLog::Info(LogCalcMcConnell, "McConnellResult::Compute",
        "atom_bond_pairs=" + std::to_string(total_pairs) +
        " xh_skipped=" + std::to_string(xh_skipped) +
        " aromatic_zeroed_pairs=" + std::to_string(aromatic_zeroed_pairs) +
        " bo_zeroed_pairs=" + std::to_string(bo_zeroed_pairs) +
        " near_lt3A={accepted=" + std::to_string(near_accepted_lt3) +
        ", rejected=" + std::to_string(near_rejected_lt3) + "}" +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " bonds=" + std::to_string(protein.BondCount()));

    return result_ptr;
}


double McConnellResult::CategoryScalarSum(size_t atom_index,
                                          BondCategory cat) const {
    const auto& ca = conf_->AtomAt(atom_index);
    switch (cat) {
        case BondCategory::PeptideCO: return ca.mcconnell_co_sum;
        case BondCategory::PeptideCN: return ca.mcconnell_cn_sum;
        case BondCategory::SidechainCO:
        case BondCategory::SidechainOther:
        case BondCategory::Disulfide:
            return ca.mcconnell_sidechain_sum;
        case BondCategory::Aromatic:
            return ca.mcconnell_aromatic_sum;
        default:
            return 0.0;
    }
}

double McConnellResult::NearestCOScalarContribution(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mcconnell_co_nearest;
}


SphericalTensor McConnellResult::SampleKernelAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    const bool include_xh_sources =
        CalculatorConfig::Get("mcconnell_include_xh_sources") != 0.0;
    const double cutoff =
        CalculatorConfig::Get("mcconnell_bond_anisotropy_cutoff");

    Mat3 total = Mat3::Zero();
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        if (!include_xh_sources && IsXHBond(protein, bond)) continue;
        const McConnellSourceCategory cat =
            SourceCategoryForMcConnell(protein, bond);
        if (cat == McConnellSourceCategory::AromaticZeroed) continue;

        const auto kernel = ComputePairKernel(
            point, conf_->bond_midpoints[bi], conf_->bond_directions[bi]);
        if (kernel.distance < CalculatorConfig::Get("singularity_guard_distance"))
            continue;
        if (kernel.distance > cutoff)
            continue;

        const double threshold =
            CalculatorConfig::Get("near_field_exclusion_ratio") *
            conf_->bond_lengths[bi];
        if (kernel.distance < threshold)
            continue;

        total += kernel.response;
    }

    return SphericalTensor::Decompose(total);
}


int McConnellResult::WriteFeatures(const ProteinConformation& conf,
                                   const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    constexpr int kCols = 9;

    int written = 0;
    for (std::size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
        const auto cat = static_cast<McConnellSourceCategory>(c);
        for (std::size_t h = 0; h < kMcConnellChannelCount; ++h) {
            const auto channel = static_cast<McConnellChannel>(h);
            std::vector<double> data(N * kCols, 0.0);
            for (size_t i = 0; i < N; ++i) {
                conf.AtomAt(i).mcconnell_source_tensors[c][h]
                    .PackFull9(&data[i * kCols]);
            }

            const std::string stem =
                std::string("mc_") + McConnellSourceCategoryStem(cat) +
                "_" + McConnellChannelStem(channel);
            if (NpyWriter::WriteFloat64(output_dir + "/" + stem + ".npy",
                                        data.data(), N, kCols)) {
                ++written;
            }
        }
    }

    {
        constexpr int kCountCols = 2;
        std::vector<int32_t> counts(N * kCountCols, 0);
        for (size_t i = 0; i < N; ++i) {
            counts[i * kCountCols + 0] = static_cast<int32_t>(
                conf.AtomAt(i).mcconnell_near_field_accepted_lt3A);
            counts[i * kCountCols + 1] = static_cast<int32_t>(
                conf.AtomAt(i).mcconnell_near_field_rejected_lt3A);
        }
        if (NpyWriter::WriteInt32(output_dir + "/mc_nearfield_counts.npy",
                                  counts.data(), N, kCountCols)) {
            ++written;
        }
    }

    MergeMcConnellManifest(
        output_dir,
        CalculatorConfig::Get("mcconnell_include_xh_sources") != 0.0);
    return written;
}

}  // namespace nmr
