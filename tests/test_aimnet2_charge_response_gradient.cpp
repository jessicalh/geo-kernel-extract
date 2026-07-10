// Smoke test for AIMNet2ChargeResponseGradientResult.
//
// Mirrors the test_apbs_ff14sb.cpp pattern: load 1UBQ via
// BuildFromProtonatedPdb, attach the dependency chain, run
// AIMNet2ChargeResponseGradientResult and verify gradients + WriteFeatures.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AIMNet2ChargeResponseGradientResult.h"
#include "AIMNet2Result.h"
#include "CalculatorConfig.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "ConformationAtom.h"
#include "EeqResult.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "PdbFileReader.h"
#include "PhysicalConstants.h"
#include "OperationRunner.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "SpatialIndexResult.h"
#include "generated/AIMNet2AimProjection.h"

#include <torch/cuda.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

using namespace nmr;
namespace fs = std::filesystem;

namespace {

struct NpyArray {
    std::string descr;
    std::vector<std::size_t> shape;
    std::vector<char> bytes;
};

std::string Trim(std::string value) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(),
                             [&](char c) { return !is_space(c); }));
    value.erase(std::find_if(value.rbegin(), value.rend(),
                             [&](char c) { return !is_space(c); }).base(),
                value.end());
    return value;
}

NpyArray ReadNpy(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    EXPECT_TRUE(input.is_open()) << path;
    NpyArray array;
    if (!input.is_open()) return array;

    char magic[6] = {};
    input.read(magic, sizeof(magic));
    EXPECT_EQ(std::string(magic, sizeof(magic)),
              std::string("\x93NUMPY", sizeof(magic))) << path;
    std::uint8_t major = 0;
    std::uint8_t minor = 0;
    input.read(reinterpret_cast<char*>(&major), sizeof(major));
    input.read(reinterpret_cast<char*>(&minor), sizeof(minor));
    EXPECT_EQ(major, 1u) << path;
    EXPECT_EQ(minor, 0u) << path;

    std::uint16_t header_length = 0;
    input.read(reinterpret_cast<char*>(&header_length),
               sizeof(header_length));
    std::string header(header_length, '\0');
    input.read(header.data(), static_cast<std::streamsize>(header.size()));

    const std::string descr_key = "'descr': '";
    auto descr_begin = header.find(descr_key);
    if (descr_begin == std::string::npos) {
        ADD_FAILURE() << header;
        return array;
    }
    descr_begin += descr_key.size();
    const auto descr_end = header.find('\'', descr_begin);
    if (descr_end == std::string::npos) {
        ADD_FAILURE() << header;
        return array;
    }
    array.descr = header.substr(descr_begin, descr_end - descr_begin);

    const auto shape_key = header.find("'shape': (");
    if (shape_key == std::string::npos) {
        ADD_FAILURE() << header;
        return array;
    }
    const auto shape_begin = header.find('(', shape_key);
    const auto shape_end = header.find(')', shape_begin);
    if (shape_begin == std::string::npos || shape_end == std::string::npos) {
        ADD_FAILURE() << header;
        return array;
    }
    std::stringstream shape_stream(
        header.substr(shape_begin + 1, shape_end - shape_begin - 1));
    std::string token;
    while (std::getline(shape_stream, token, ',')) {
        token = Trim(token);
        if (!token.empty()) {
            array.shape.push_back(
                static_cast<std::size_t>(std::stoull(token)));
        }
    }
    array.bytes.assign(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
    return array;
}

template <typename T>
const T* DataAs(const NpyArray& array) {
    return reinterpret_cast<const T*>(array.bytes.data());
}

void RemoveDirectoryContents(const fs::path& directory) {
    std::error_code ec;
    if (!fs::exists(directory, ec)) return;
    for (const auto& entry : fs::directory_iterator(directory, ec)) {
        fs::remove(entry.path(), ec);
    }
    fs::remove(directory, ec);
}

// Independent closed form of the frozen model asset's `LRCoulomb` head.
// The committed model reports method="simple", rc=4.6 A and the symmetric
// neighbour-list factor 7.199822675975274 eV A/e^2.  Summing i<j therefore
// uses twice that factor.  This test helper never calls a Torch output head.
double RecomputeSimpleLongRangeEnergy(
        const ProteinConformation& conf,
        const std::vector<double>& charges) {
    constexpr double kLongRangeCutoffA = 15.0;
    constexpr double kSwitchRadiusA = 4.6;
    constexpr double kPairFactor = 2.0 * 7.199822675975274;
    constexpr double kMaxSwitchCoordinate = 0.999999;
    const double exp_minus_one = std::exp(-1.0);
    double energy = 0.0;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        for (std::size_t j = i + 1; j < conf.AtomCount(); ++j) {
            const double distance =
                (conf.PositionAt(i) - conf.PositionAt(j)).norm();
            if (distance <= 0.0 || distance > kLongRangeCutoffA) continue;
            const double x = std::clamp(
                distance / kSwitchRadiusA, 0.0, kMaxSwitchCoordinate);
            const double exp_cutoff =
                std::exp(-1.0 / (1.0 - x*x)) / exp_minus_one;
            energy += kPairFactor * charges[i] * charges[j] *
                      (1.0 - exp_cutoff) / distance;
        }
    }
    return energy;
}

class ScopedNeutralConditioning {
public:
    ScopedNeutralConditioning() {
        const std::string suffix = std::to_string(::getpid());
        enable_path_ = fs::temp_directory_path() /
            ("aimnet2_neutral_on_" + suffix + ".toml");
        restore_path_ = fs::temp_directory_path() /
            ("aimnet2_neutral_off_" + suffix + ".toml");
        WriteConfig(enable_path_, 1.0);
        WriteConfig(restore_path_, 0.0);
        CalculatorConfig::Load(enable_path_.string());
    }

    ScopedNeutralConditioning(const ScopedNeutralConditioning&) = delete;
    ScopedNeutralConditioning& operator=(
        const ScopedNeutralConditioning&) = delete;

    ~ScopedNeutralConditioning() {
        // CalculatorConfig::Load intentionally retains keys absent from a
        // later file.  Restore the knob explicitly before reloading the
        // normal test config; merely loading calculator_params.toml would
        // leave the neutral=1 override live for later tests.
        CalculatorConfig::Load(restore_path_.string());
        test::TestEnvironment::LoadCalculatorConfig();
        std::error_code ec;
        fs::remove(enable_path_, ec);
        fs::remove(restore_path_, ec);
    }

private:
    static void WriteConfig(const fs::path& path, double value) {
        std::ofstream output(path);
        if (!output.is_open()) {
            throw std::runtime_error("cannot create test config " +
                                     path.string());
        }
        output << "charge_conditioning_neutral = " << value << "\n";
    }

    fs::path enable_path_;
    fs::path restore_path_;
};

float ExpectedProjectionWeight(std::uint64_t element_slot,
                               std::uint64_t component,
                               std::uint64_t input_dim) {
    std::uint64_t x = 0xA17E20260708ULL ^ (element_slot << 40) ^
                      (component << 24) ^ input_dim;
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    const std::uint64_t u = x ^ (x >> 31);
    constexpr float scale = static_cast<float>(
        0.306186217847897235151236163473);  // sqrt(3/32)
    if (u % 6u == 0u) return scale;
    if (u % 6u == 1u) return -scale;
    return 0.0f;
}

}  // namespace


TEST(AIMNet2AimProjectionBasis, LiteralTableMatchesPinnedProvenance) {
    static_assert(kAimnet2AimProjectionDims == 32);
    static_assert(AIMNET2_AIM_DIMS == 256);

    constexpr float kScale = 0x1.3988e2p-2f;
    std::size_t positive = 0;
    std::size_t negative = 0;
    std::size_t zero = 0;
    for (std::size_t e = 0; e < kAimnet2AimProjectionElementSlots; ++e) {
        for (std::size_t k = 0; k < kAimnet2AimProjectionDims; ++k) {
            for (std::size_t d = 0; d < AIMNET2_AIM_DIMS; ++d) {
                const float value = kAimnet2AimProjectionBasis[e][k][d];
                EXPECT_FLOAT_EQ(value, ExpectedProjectionWeight(e, k, d))
                    << "basis mismatch at [" << e << "][" << k
                    << "][" << d << "]";
                if (value > 0.0f) {
                    ++positive;
                    EXPECT_FLOAT_EQ(value, kScale);
                } else if (value < 0.0f) {
                    ++negative;
                    EXPECT_FLOAT_EQ(value, -kScale);
                } else {
                    ++zero;
                }
            }
        }
    }
    EXPECT_EQ(positive, 6813u);
    EXPECT_EQ(negative, 6791u);
    EXPECT_EQ(zero, 27356u);

    // Pinned cells were produced independently from the binding offline
    // SplitMix64 recipe; they catch axis transposition as well as sign drift.
    EXPECT_FLOAT_EQ(kAimnet2AimProjectionBasis[0][0][0], kScale);
    EXPECT_FLOAT_EQ(kAimnet2AimProjectionBasis[0][0][1], 0.0f);
    EXPECT_FLOAT_EQ(kAimnet2AimProjectionBasis[2][17][93], -kScale);
    EXPECT_FLOAT_EQ(kAimnet2AimProjectionBasis[4][31][255], -kScale);
    EXPECT_EQ(kAimnet2AimProjectionBasisId,
              "splitmix64_0xA17E20260708_achlioptas_32x256_element_HCNOS");
}


TEST(AIMNet2FrozenOutputs,
     IndependentOffGpuEfgT2AndEnergyRecomputePinsBlessedPayloads) {
#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif
    test::TestEnvironment::LoadCalculatorConfig();
    auto build = BuildFromProtonatedPdb(
        test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto& conf = build.protein->Conformation();
    const Protein& protein = conf.ProteinRef();
    const std::size_t N = conf.AtomCount();
    ASSERT_GT(N, 0u);

    const fs::path golden_dir = fs::path(NMR_TEST_DATA_DIR).parent_path() /
        "golden" / "blessed" / "nodft";
    const auto charges_array = ReadNpy(golden_dir / "aimnet2_charges.npy");
    const auto efg_total_array = ReadNpy(golden_dir / "aimnet2_efg.npy");
    const auto efg_backbone_array =
        ReadNpy(golden_dir / "aimnet2_efg_backbone.npy");
    const auto efg_sidechain_array =
        ReadNpy(golden_dir / "aimnet2_efg_sidechain.npy");
    const auto efg_aromatic_array =
        ReadNpy(golden_dir / "aimnet2_efg_aromatic.npy");

    ASSERT_EQ(charges_array.descr, "<f8");
    ASSERT_EQ(charges_array.shape, (std::vector<std::size_t>{N}));
    for (const auto* array : {&efg_total_array, &efg_backbone_array,
                              &efg_sidechain_array, &efg_aromatic_array}) {
        ASSERT_EQ(array->descr, "<f8");
        ASSERT_EQ(array->shape,
                  (std::vector<std::size_t>{N, std::size_t(5)}));
        ASSERT_EQ(array->bytes.size(), N * 5 * sizeof(double));
    }
    ASSERT_EQ(charges_array.bytes.size(), N * sizeof(double));

    const double* charge_payload = DataAs<double>(charges_array);
    std::vector<double> charges(charge_payload, charge_payload + N);
    const std::array<const double*, 4> frozen_t2 = {
        DataAs<double>(efg_total_array),
        DataAs<double>(efg_backbone_array),
        DataAs<double>(efg_sidechain_array),
        DataAs<double>(efg_aromatic_array),
    };

    std::vector<bool> is_backbone(N, false);
    std::vector<bool> is_aromatic(N, false);
    for (std::size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& residue = protein.ResidueAt(ri);
        for (std::size_t index : {residue.N, residue.CA, residue.C,
                                  residue.O, residue.H, residue.HA,
                                  residue.CB}) {
            if (index != Residue::NONE && index < N) {
                is_backbone[index] = true;
            }
        }
    }
    for (std::size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (std::size_t index : protein.RingAt(ri).atom_indices) {
            ASSERT_LT(index, N);
            is_aromatic[index] = true;
        }
    }
    ASSERT_GT(std::count(is_backbone.begin(), is_backbone.end(), true), 0);
    ASSERT_GT(std::count(is_aromatic.begin(), is_aromatic.end(), true), 0);

    constexpr double kCutoffA = 20.0;
    constexpr double kChargeFloor = 1e-15;
    constexpr double kSingularityGuardA = 0.1;
    constexpr double kSqrt2 = 1.4142135623730950488;
    constexpr double kSqrtThreeHalves = 1.2247448713915890491;
    constexpr double kT2Tolerance = 5e-4;
    std::array<double, 4> max_signal = {};
    std::array<double, 4> max_error = {};

    for (std::size_t target = 0; target < N; ++target) {
        std::array<Mat3, 4> expected = {
            Mat3::Zero(), Mat3::Zero(), Mat3::Zero(), Mat3::Zero()
        };
        for (std::size_t source = 0; source < N; ++source) {
            if (source == target ||
                std::abs(charges[source]) < kChargeFloor) {
                continue;
            }
            const Vec3 r = conf.PositionAt(target) - conf.PositionAt(source);
            const double distance = r.norm();
            if (distance > kCutoffA || distance < kSingularityGuardA) {
                continue;
            }
            const double r3 = distance * distance * distance;
            const double r5 = r3 * distance * distance;
            const Mat3 contribution = COULOMB_KE * charges[source] *
                (3.0 * r * r.transpose() / r5 - Mat3::Identity() / r3);
            expected[0] += contribution;
            // Aromatic precedence is part of the frozen source partition.
            if (is_aromatic[source]) {
                expected[3] += contribution;
            } else if (is_backbone[source]) {
                expected[1] += contribution;
            } else {
                expected[2] += contribution;
            }
        }

        for (std::size_t role = 0; role < expected.size(); ++role) {
            Mat3& matrix = expected[role];
            matrix -= (matrix.trace() / 3.0) * Mat3::Identity();

            // Manual real-T2 basis: deliberately do not call
            // SphericalTensor::Decompose or PackT2.  A sign or output-index
            // change in the production basis must fail this off-GPU test.
            const std::array<double, 5> manual_t2 = {
                kSqrt2 * matrix(0, 1),
                kSqrt2 * matrix(1, 2),
                kSqrtThreeHalves * matrix(2, 2),
                kSqrt2 * matrix(0, 2),
                (matrix(0, 0) - matrix(1, 1)) / kSqrt2,
            };
            for (std::size_t component = 0; component < 5; ++component) {
                const double actual =
                    frozen_t2[role][target * 5 + component];
                const double error = std::abs(actual - manual_t2[component]);
                max_error[role] = std::max(max_error[role], error);
                max_signal[role] = std::max(max_signal[role],
                                            std::abs(actual));
                EXPECT_NEAR(actual, manual_t2[component], kT2Tolerance)
                    << "target=" << target << " role=" << role
                    << " T2_component=" << component;
            }
        }
    }
    for (std::size_t role = 0; role < max_signal.size(); ++role) {
        EXPECT_GT(max_signal[role], 1.0)
            << "frozen fixture must exercise role " << role;
        EXPECT_LT(max_error[role], kT2Tolerance);
    }

    // Non-telescoping energy check on the committed NPY contract.  Local and
    // D3 terms come from independent per-atom payload sums; long-range
    // Coulomb is recomputed directly from charges/geometry using the frozen
    // model head's documented closed form.  None is obtained by subtracting
    // adjacent total-energy stages.
    const auto energy_mlp_array =
        ReadNpy(golden_dir / "aimnet2_energy_mlp.npy");
    const auto shifted_array =
        ReadNpy(golden_dir / "aimnet2_energy_shifted_local.npy");
    const auto d3_atom_array =
        ReadNpy(golden_dir / "aimnet2_d3_e_disp_atom.npy");
    const auto energy_terms_array =
        ReadNpy(golden_dir / "aimnet2_energy_terms.npy");
    const auto d3_cn_array = ReadNpy(golden_dir / "aimnet2_d3_cn.npy");
    const auto d3_c6_array =
        ReadNpy(golden_dir / "aimnet2_d3_c6_stats.npy");
    for (const auto* array : {&energy_mlp_array, &shifted_array,
                              &d3_atom_array, &d3_cn_array}) {
        ASSERT_EQ(array->descr, "<f8");
        ASSERT_EQ(array->shape, (std::vector<std::size_t>{N}));
        ASSERT_EQ(array->bytes.size(), N * sizeof(double));
    }
    ASSERT_EQ(energy_terms_array.descr, "<f8");
    ASSERT_EQ(energy_terms_array.shape,
              (std::vector<std::size_t>{1u, 6u}));
    ASSERT_EQ(d3_c6_array.descr, "<f8");
    ASSERT_EQ(d3_c6_array.shape,
              (std::vector<std::size_t>{N, 3u}));
    ASSERT_EQ(d3_c6_array.bytes.size(), N * 3 * sizeof(double));

    const double* energy_mlp = DataAs<double>(energy_mlp_array);
    const double* shifted = DataAs<double>(shifted_array);
    const double* d3_atom = DataAs<double>(d3_atom_array);
    const double* terms = DataAs<double>(energy_terms_array);
    const double* d3_cn = DataAs<double>(d3_cn_array);
    const double* d3_c6 = DataAs<double>(d3_c6_array);
    double local_sum = 0.0;
    double d3_sum = 0.0;
    double charge_sum = 0.0;
    int positive_c6_rows = 0;
    for (std::size_t i = 0; i < N; ++i) {
        EXPECT_TRUE(std::isfinite(energy_mlp[i]));
        EXPECT_TRUE(std::isfinite(shifted[i]));
        EXPECT_TRUE(std::isfinite(d3_atom[i]));
        EXPECT_TRUE(std::isfinite(d3_cn[i]));
        EXPECT_GE(d3_cn[i], 0.0);
        local_sum += shifted[i];
        d3_sum += d3_atom[i];
        charge_sum += charges[i];

        const double c6_sum = d3_c6[i*3];
        const double c6_mean = d3_c6[i*3 + 1];
        const double c6_max = d3_c6[i*3 + 2];
        EXPECT_TRUE(std::isfinite(c6_sum));
        EXPECT_TRUE(std::isfinite(c6_mean));
        EXPECT_TRUE(std::isfinite(c6_max));
        EXPECT_GE(c6_sum, 0.0);
        EXPECT_GE(c6_mean, 0.0);
        EXPECT_GE(c6_max, 0.0);
        EXPECT_GE(c6_sum + 1e-10, c6_max);
        EXPECT_GE(c6_max + 1e-10, c6_mean);
        if (c6_sum > 0.0) ++positive_c6_rows;
    }
    EXPECT_EQ(positive_c6_rows, static_cast<int>(N));

    // CN/C6 values are private-model diagnostics: recomputing them here
    // would duplicate the model internals and be circular.  The checks above
    // pin their frozen numeric shape/finiteness and [sum,mean,max] column
    // semantics.  D3 energy itself has the independent per-atom-to-term sum.
    const double lrc_recomputed =
        RecomputeSimpleLongRangeEnergy(conf, charges);
    EXPECT_NEAR(terms[0], local_sum, 5e-8);
    EXPECT_NEAR(terms[1], lrc_recomputed, 5e-5);
    EXPECT_NEAR(terms[2], d3_sum, 1e-5);
    EXPECT_NEAR(terms[3], local_sum + lrc_recomputed + d3_sum, 1e-4);
    EXPECT_NEAR(terms[4], charge_sum, 1e-5);
    EXPECT_DOUBLE_EQ(terms[5], 0.0);
}


class AIMNet2ChargeResponseGradientTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!torch::cuda::is_available()) {
            GTEST_SKIP() << "CUDA not available; AIMNet2 requires GPU.";
        }
        if (!fs::exists(test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found at "
                         << test::TestEnvironment::UbqProtonated();
        }
        if (!fs::exists(test::TestEnvironment::Ff14sbParams())) {
            GTEST_SKIP() << "ff14sb_params.dat not found at "
                         << test::TestEnvironment::Ff14sbParams();
        }
        if (!fs::exists(test::TestEnvironment::Aimnet2Model())) {
            GTEST_SKIP() << "AIMNet2 .jpt not found at "
                         << test::TestEnvironment::Aimnet2Model();
        }

        auto r = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);

        model = AIMNet2Model::Load(test::TestEnvironment::Aimnet2Model());
        if (!model) {
            GTEST_SKIP() << "AIMNet2Model::Load returned null";
        }
    }

    std::unique_ptr<Protein> protein;
    std::unique_ptr<AIMNet2Model> model;
};


TEST_F(AIMNet2ChargeResponseGradientTest, PipelineProducesNonZeroChargeResponseGradient) {
    auto& conf = protein->Conformation();
    constexpr int kInputNetCharge = -2;
    ParamFileChargeSource charge_source(
        test::TestEnvironment::Ff14sbParams());
    RunOptions options;
    options.charge_source = &charge_source;
    options.net_charge = kInputNetCharge;
    options.skip_dssp = true;
    options.skip_mopac = true;
    options.skip_apbs = true;
    options.skip_coulomb = true;
    options.aimnet2_model = model.get();

    // A5 forcing: enter through the real producer order so the
    // charge_conditioning_neutral policy, not a direct calculator argument,
    // determines what both AIMNet2 calculators receive.
    const RunResult run = OperationRunner::Run(conf, options);
    ASSERT_TRUE(run.Ok()) << run.error;
    ASSERT_TRUE(conf.HasResult<ChargeAssignmentResult>());
    ASSERT_TRUE(conf.HasResult<AIMNet2Result>());
    ASSERT_TRUE(conf.HasResult<AIMNet2ChargeResponseGradientResult>());
    const int expected_conditioned_charge =
        CalculatorConfig::Get("charge_conditioning_neutral") != 0.0
            ? 0 : kInputNetCharge;

    // Independent forcing functions for the new AIMNet2 products.
    const auto& aim_result = conf.Result<AIMNet2Result>();
    EXPECT_NEAR(aim_result.EnergyTotal(),
                aim_result.EnergyLocalSum() +
                    aim_result.EnergyLongRangeCoulomb() +
                    aim_result.EnergyDftD3(),
                1.0e-8);
    EXPECT_DOUBLE_EQ(aim_result.ConditionedNetCharge(),
                     static_cast<double>(expected_conditioned_charge));
    EXPECT_DOUBLE_EQ(
        aim_result.NeutralConditioningFlag(),
        CalculatorConfig::Get("charge_conditioning_neutral") != 0.0
            ? 1.0 : 0.0);

    double shifted_local_sum = 0.0;
    double d3_atom_sum = 0.0;
    double aim_charge_sum = 0.0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        shifted_local_sum += ca.aimnet2_energy_shifted_local;
        d3_atom_sum += ca.aimnet2_d3_e_disp_atom;
        aim_charge_sum += ca.aimnet2_charge;
        EXPECT_TRUE(std::isfinite(ca.aimnet2_energy_mlp));
        EXPECT_TRUE(std::isfinite(ca.aimnet2_energy_shifted_local));
        EXPECT_TRUE(std::isfinite(ca.aimnet2_d3_cn));
        for (double value : ca.aimnet2_d3_c6_stats) {
            EXPECT_TRUE(std::isfinite(value));
        }

        EXPECT_LT((ca.aimnet2_E_total -
                   (ca.aimnet2_E_backbone + ca.aimnet2_E_sidechain +
                    ca.aimnet2_E_aromatic)).norm(),
                  1.0e-9);
        EXPECT_LT((ca.aimnet2_EFG_total -
                   (ca.aimnet2_EFG_backbone + ca.aimnet2_EFG_sidechain +
                    ca.aimnet2_EFG_aromatic)).norm(),
                  1.0e-9);

        const Element element = conf.ProteinRef().AtomAt(i).element;
        std::size_t slot = 0;
        switch (element) {
            case Element::H: slot = 0; break;
            case Element::C: slot = 1; break;
            case Element::N: slot = 2; break;
            case Element::O: slot = 3; break;
            case Element::S: slot = 4; break;
            case Element::Unknown: FAIL() << "guarded AIMNet2 element"; break;
        }
        for (std::size_t k = 0; k < kAimnet2AimProjectionDims; ++k) {
            double expected = 0.0;
            for (std::size_t d = 0; d < AIMNET2_AIM_DIMS; ++d) {
                expected += static_cast<double>(
                    kAimnet2AimProjectionBasis[slot][k][d]) *
                    static_cast<double>(ca.aimnet2_aim[d]);
            }
            EXPECT_FLOAT_EQ(ca.aimnet2_aim_projection[k],
                            static_cast<float>(expected));
        }
    }
    const double d3_tolerance = std::max(
        1.0e-5,
        1.0e-6 * std::max(1.0, std::abs(aim_result.EnergyDftD3())));
    EXPECT_NEAR(shifted_local_sum, aim_result.EnergyLocalSum(), 1.0e-7);
    EXPECT_NEAR(d3_atom_sum, aim_result.EnergyDftD3(), d3_tolerance);
    EXPECT_NEAR(aim_charge_sum,
                static_cast<double>(expected_conditioned_charge), 1.0e-5);

    // Brute-force Coulomb source of truth (independent of SpatialIndexResult)
    // pins the rank-1 sign convention r = target - source.
    const double cutoff =
        CalculatorConfig::Get("aimnet2_coulomb_efg_cutoff");
    const double charge_floor =
        CalculatorConfig::Get("coulomb_charge_noise_floor");
    const double singularity_guard =
        CalculatorConfig::Get("singularity_guard_distance");
    for (std::size_t target :
         {std::size_t(0), conf.AtomCount() / 2, conf.AtomCount() - 1}) {
        Vec3 expected_E = Vec3::Zero();
        for (std::size_t source = 0; source < conf.AtomCount(); ++source) {
            if (source == target) continue;
            const double q = conf.AtomAt(source).aimnet2_charge;
            if (std::abs(q) < charge_floor) continue;
            const Vec3 r = conf.PositionAt(target) - conf.PositionAt(source);
            const double distance = r.norm();
            if (distance > cutoff || distance < singularity_guard) continue;
            expected_E += COULOMB_KE * q * r /
                (distance * distance * distance);
        }
        const Vec3 actual_E = conf.AtomAt(target).aimnet2_E_total;
        EXPECT_LT((actual_E - expected_E).norm(),
                  1.0e-9 * std::max(1.0, expected_E.norm()));
    }

    const size_t N = conf.AtomCount();
    ASSERT_GT(N, 0u);

    int finite_count = 0;
    int nonzero_count = 0;
    double max_scalar = 0.0;
    double sum_scalar = 0.0;
    double max_consistency_diff = 0.0;

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        const Vec3& v = ca.aimnet2_charge_response_gradient_vector;
        const double s = ca.aimnet2_charge_response_gradient_scalar;

        if (std::isfinite(v.x()) && std::isfinite(v.y()) &&
            std::isfinite(v.z()) && std::isfinite(s)) {
            finite_count++;
        }

        const double v_norm = v.norm();
        const double diff = std::abs(v_norm - s);
        if (diff > max_consistency_diff) max_consistency_diff = diff;

        if (s > 1e-10) nonzero_count++;
        if (s > max_scalar) max_scalar = s;
        sum_scalar += s;
    }
    const double mean_scalar = sum_scalar / static_cast<double>(N);

    fprintf(stderr,
        "\n=== AIMNet2 charge-response gradient summary (1UBQ, %zu atoms) ===\n"
        "  finite values:         %d / %zu\n"
        "  non-zero scalar:       %d / %zu\n"
        "  max scalar:            %.6e\n"
        "  mean scalar:           %.6e\n"
        "  max |scalar - ||v|||:  %.6e\n"
        "==========================================================\n",
        N, finite_count, N, nonzero_count, N,
        max_scalar, mean_scalar, max_consistency_diff);

    EXPECT_EQ(finite_count, static_cast<int>(N))
        << "Some charge-response gradient values are NaN or Inf";
    EXPECT_GT(nonzero_count, static_cast<int>(N) / 4)
        << "Too few atoms with non-zero charge-response gradient.";
    EXPECT_GT(max_scalar, 1e-8)
        << "Maximum gradient norm is suspiciously small";
    EXPECT_LT(max_consistency_diff, 1e-9)
        << "Scalar field does not match L2 norm of vector field";
}


TEST_F(AIMNet2ChargeResponseGradientTest,
       TrpCageEnergyTermsUseNeutralConditioningAndIndependentRecompute) {
#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif
    const fs::path trp_cage = fs::path(NMR_TEST_DATA_DIR) /
        "illustrative_peptides" / "trp_cage_1l2y_model1.pdb";
    ASSERT_TRUE(fs::exists(trp_cage)) << trp_cage;
    auto build = BuildFromProtonatedPdb(trp_cage.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto& conf = build.protein->Conformation();

    // Trp-cage TC5b has a nonzero physical charge (+1).  Turn on the actual
    // producer ablation knob, enter through OperationRunner, and require both
    // charge models to receive 0 rather than merely testing a direct
    // AIMNet2Result::Compute argument.
    ScopedNeutralConditioning neutral_scope;
    ASSERT_DOUBLE_EQ(CalculatorConfig::Get("charge_conditioning_neutral"),
                     1.0);
    ParamFileChargeSource charge_source(
        test::TestEnvironment::Ff14sbParams());
    RunOptions options;
    options.charge_source = &charge_source;
    options.net_charge = +1;
    options.skip_dssp = true;
    options.skip_mopac = true;
    options.skip_apbs = true;
    options.skip_coulomb = true;
    options.aimnet2_model = model.get();

    const RunResult run = OperationRunner::Run(conf, options);
    ASSERT_TRUE(run.Ok()) << run.error;
    ASSERT_TRUE(conf.HasResult<EeqResult>());
    ASSERT_TRUE(conf.HasResult<AIMNet2Result>());
    ASSERT_TRUE(conf.HasResult<AIMNet2ChargeResponseGradientResult>());

    const auto& aim = conf.Result<AIMNet2Result>();
    EXPECT_DOUBLE_EQ(aim.ConditionedNetCharge(), 0.0);
    EXPECT_DOUBLE_EQ(aim.NeutralConditioningFlag(), 1.0);
    std::vector<double> charges(conf.AtomCount());
    double aim_charge_sum = 0.0;
    double eeq_charge_sum = 0.0;
    double local_sum = 0.0;
    double d3_sum = 0.0;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& atom = conf.AtomAt(i);
        charges[i] = atom.aimnet2_charge;
        aim_charge_sum += atom.aimnet2_charge;
        eeq_charge_sum += atom.eeq_charge;
        local_sum += atom.aimnet2_energy_shifted_local;
        d3_sum += atom.aimnet2_d3_e_disp_atom;
    }
    EXPECT_NEAR(aim_charge_sum, 0.0, 1e-5);
    EXPECT_NEAR(eeq_charge_sum, 0.0, 1e-10);

    const double lrc_recomputed =
        RecomputeSimpleLongRangeEnergy(conf, charges);
    const double d3_tolerance = std::max(
        1e-5, 1e-6 * std::max(1.0, std::abs(aim.EnergyDftD3())));
    EXPECT_NEAR(aim.EnergyLocalSum(), local_sum, 1e-6);
    EXPECT_NEAR(aim.EnergyLongRangeCoulomb(), lrc_recomputed, 5e-5);
    EXPECT_NEAR(aim.EnergyDftD3(), d3_sum, d3_tolerance);
    EXPECT_NEAR(aim.EnergyTotal(),
                local_sum + lrc_recomputed + d3_sum,
                2e-4);
}


TEST_F(AIMNet2ChargeResponseGradientTest, WriteFeaturesEmitsBothNpys) {
    auto& conf = protein->Conformation();

    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    auto spatial = SpatialIndexResult::Compute(conf);
    ASSERT_NE(spatial, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(spatial)));

    auto enrich = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrich, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(enrich)));

    auto charges = ChargeAssignmentResult::Compute(
        conf, test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(charges, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    auto aim = AIMNet2Result::Compute(conf, *model, /*net_charge=*/0);
    ASSERT_NE(aim, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(aim)));

    auto pol = AIMNet2ChargeResponseGradientResult::Compute(conf, *model, /*net_charge=*/0);
    ASSERT_NE(pol, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(pol)));

    const fs::path output_dir = fs::temp_directory_path() /
        ("aimnet2_charge_response_gradient_test_writefeatures_" +
         std::to_string(::getpid()));
    RemoveDirectoryContents(output_dir);
    ASSERT_FALSE(fs::exists(output_dir));
    ASSERT_TRUE(fs::create_directories(output_dir));

    // Hardened calculator-language invariant: projection is born in Compute
    // and WriteFeatures reads it back. Destroy the raw embedding after
    // Compute; the emitted projection must retain the stored value.
    const float stored_projection_0 =
        conf.AtomAt(0).aimnet2_aim_projection[0];
    conf.MutableAtomAt(0).aimnet2_aim.fill(12345.0f);

    const auto& aim_result = conf.Result<AIMNet2Result>();
    const int aim_written = aim_result.WriteFeatures(conf, output_dir.string());
    ASSERT_EQ(aim_written, 17);

    const size_t N = conf.AtomCount();
    const auto projection =
        ReadNpy(output_dir / "aimnet2_aim_projection.npy");
    ASSERT_EQ(projection.descr, "<f4");
    ASSERT_EQ(projection.shape,
              (std::vector<size_t>{N, kAimnet2AimProjectionDims}));
    ASSERT_EQ(projection.bytes.size(),
              N * kAimnet2AimProjectionDims * sizeof(float));
    const float* projection_data = DataAs<float>(projection);
    EXPECT_FLOAT_EQ(projection_data[0], stored_projection_0);
    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < kAimnet2AimProjectionDims; ++k) {
            EXPECT_FLOAT_EQ(
                projection_data[i*kAimnet2AimProjectionDims + k],
                conf.AtomAt(i).aimnet2_aim_projection[k]);
        }
    }

    struct VectorOutput {
        const char* filename;
        Vec3 ConformationAtom::* member;
    };
    const VectorOutput vector_outputs[] = {
        {"aimnet2_E.npy", &ConformationAtom::aimnet2_E_total},
        {"aimnet2_E_backbone.npy", &ConformationAtom::aimnet2_E_backbone},
        {"aimnet2_E_sidechain.npy", &ConformationAtom::aimnet2_E_sidechain},
        {"aimnet2_E_aromatic.npy", &ConformationAtom::aimnet2_E_aromatic},
    };
    for (const auto& output : vector_outputs) {
        const auto array = ReadNpy(output_dir / output.filename);
        ASSERT_EQ(array.descr, "<f8") << output.filename;
        ASSERT_EQ(array.shape, (std::vector<size_t>{N, 3u}))
            << output.filename;
        ASSERT_EQ(array.bytes.size(), N * 3 * sizeof(double));
        const double* data = DataAs<double>(array);
        for (size_t i = 0; i < N; ++i) {
            const Vec3& expected = conf.AtomAt(i).*(output.member);
            for (int d = 0; d < 3; ++d) {
                EXPECT_DOUBLE_EQ(data[i*3 + static_cast<size_t>(d)],
                                 expected(d))
                    << output.filename << " atom=" << i << " d=" << d;
            }
        }
    }

    const auto efg_sidechain =
        ReadNpy(output_dir / "aimnet2_efg_sidechain.npy");
    ASSERT_EQ(efg_sidechain.descr, "<f8");
    ASSERT_EQ(efg_sidechain.shape, (std::vector<size_t>{N, 5u}));
    ASSERT_EQ(efg_sidechain.bytes.size(), N * 5 * sizeof(double));
    const double* efg_sidechain_data = DataAs<double>(efg_sidechain);
    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < 5; ++k) {
            EXPECT_DOUBLE_EQ(efg_sidechain_data[i*5 + k],
                             conf.AtomAt(i)
                                 .aimnet2_EFG_sidechain_spherical.T2[k]);
        }
    }

    struct ScalarOutput {
        const char* filename;
        double ConformationAtom::* member;
    };
    const ScalarOutput scalar_outputs[] = {
        {"aimnet2_energy_mlp.npy",
         &ConformationAtom::aimnet2_energy_mlp},
        {"aimnet2_energy_shifted_local.npy",
         &ConformationAtom::aimnet2_energy_shifted_local},
        {"aimnet2_d3_e_disp_atom.npy",
         &ConformationAtom::aimnet2_d3_e_disp_atom},
        {"aimnet2_d3_cn.npy", &ConformationAtom::aimnet2_d3_cn},
    };
    for (const auto& output : scalar_outputs) {
        const auto array = ReadNpy(output_dir / output.filename);
        ASSERT_EQ(array.descr, "<f8") << output.filename;
        ASSERT_EQ(array.shape, (std::vector<size_t>{N}))
            << output.filename;
        ASSERT_EQ(array.bytes.size(), N * sizeof(double));
        const double* data = DataAs<double>(array);
        for (size_t i = 0; i < N; ++i) {
            EXPECT_DOUBLE_EQ(data[i], conf.AtomAt(i).*(output.member))
                << output.filename << " atom=" << i;
        }
    }

    const auto c6_stats =
        ReadNpy(output_dir / "aimnet2_d3_c6_stats.npy");
    ASSERT_EQ(c6_stats.descr, "<f8");
    ASSERT_EQ(c6_stats.shape, (std::vector<size_t>{N, 3u}));
    ASSERT_EQ(c6_stats.bytes.size(), N * 3 * sizeof(double));
    const double* c6_data = DataAs<double>(c6_stats);
    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < 3; ++k) {
            EXPECT_DOUBLE_EQ(c6_data[i*3 + k],
                             conf.AtomAt(i).aimnet2_d3_c6_stats[k]);
        }
    }

    const auto energy_terms =
        ReadNpy(output_dir / "aimnet2_energy_terms.npy");
    ASSERT_EQ(energy_terms.descr, "<f8");
    ASSERT_EQ(energy_terms.shape, (std::vector<size_t>{1u, 6u}));
    ASSERT_EQ(energy_terms.bytes.size(), 6 * sizeof(double));
    const double* energy_data = DataAs<double>(energy_terms);
    EXPECT_DOUBLE_EQ(energy_data[0], aim_result.EnergyLocalSum());
    EXPECT_DOUBLE_EQ(energy_data[1], aim_result.EnergyLongRangeCoulomb());
    EXPECT_DOUBLE_EQ(energy_data[2], aim_result.EnergyDftD3());
    EXPECT_DOUBLE_EQ(energy_data[3], aim_result.EnergyTotal());
    EXPECT_DOUBLE_EQ(energy_data[4], aim_result.ConditionedNetCharge());
    EXPECT_DOUBLE_EQ(energy_data[5], aim_result.NeutralConditioningFlag());

    const auto& result = conf.Result<AIMNet2ChargeResponseGradientResult>();
    const int written = result.WriteFeatures(conf, output_dir.string());
    ASSERT_EQ(written, 2);
    const auto gradient = ReadNpy(
        output_dir / "aimnet2_charge_response_gradient.npy");
    const auto gradient_scalar = ReadNpy(
        output_dir / "aimnet2_charge_response_gradient_scalar.npy");
    ASSERT_EQ(gradient.descr, "<f8");
    ASSERT_EQ(gradient.shape, (std::vector<size_t>{N, 3u}));
    ASSERT_EQ(gradient_scalar.descr, "<f8");
    ASSERT_EQ(gradient_scalar.shape, (std::vector<size_t>{N}));
    ASSERT_EQ(gradient.bytes.size(), N * 3 * sizeof(double));
    ASSERT_EQ(gradient_scalar.bytes.size(), N * sizeof(double));
    const double* gradient_data = DataAs<double>(gradient);
    const double* gradient_scalar_data = DataAs<double>(gradient_scalar);
    for (size_t i = 0; i < N; ++i) {
        const auto& atom = conf.AtomAt(i);
        for (int d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(
                gradient_data[i*3 + static_cast<size_t>(d)],
                atom.aimnet2_charge_response_gradient_vector(d));
        }
        EXPECT_DOUBLE_EQ(
            gradient_scalar_data[i],
            atom.aimnet2_charge_response_gradient_scalar);
    }

    // Clean every file without std::filesystem::remove_all: the libtorch
    // package's interposed remove_all symbol is broken in this environment.
    // A unique, pre-cleaned directory plus numeric payload reads ensures a
    // failed writer can no longer pass using stale files from a prior run.
    RemoveDirectoryContents(output_dir);
    EXPECT_FALSE(fs::exists(output_dir));
}
