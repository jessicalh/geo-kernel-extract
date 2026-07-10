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

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>

using namespace nmr;
namespace fs = std::filesystem;

namespace {

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
        "aimnet2_charge_response_gradient_test_writefeatures";
    fs::create_directories(output_dir);

    // Hardened calculator-language invariant: projection is born in Compute
    // and WriteFeatures reads it back. Destroy the raw embedding after
    // Compute; the emitted projection must retain the stored value.
    const float stored_projection_0 =
        conf.AtomAt(0).aimnet2_aim_projection[0];
    conf.MutableAtomAt(0).aimnet2_aim.fill(12345.0f);

    const auto& aim_result = conf.Result<AIMNet2Result>();
    const int aim_written = aim_result.WriteFeatures(conf, output_dir.string());
    EXPECT_EQ(aim_written, 17);
    const fs::path projection_path =
        output_dir / "aimnet2_aim_projection.npy";
    ASSERT_TRUE(fs::exists(projection_path));
    {
        std::ifstream input(projection_path, std::ios::binary);
        ASSERT_TRUE(input.good());
        char prefix[8] = {};
        input.read(prefix, sizeof(prefix));
        ASSERT_EQ(input.gcount(), static_cast<std::streamsize>(sizeof(prefix)));
        std::uint16_t header_length = 0;
        input.read(reinterpret_cast<char*>(&header_length),
                   sizeof(header_length));
        input.seekg(static_cast<std::streamoff>(header_length),
                    std::ios::cur);
        float emitted_projection_0 = 0.0f;
        input.read(reinterpret_cast<char*>(&emitted_projection_0),
                   sizeof(emitted_projection_0));
        ASSERT_TRUE(input.good());
        EXPECT_FLOAT_EQ(emitted_projection_0, stored_projection_0);
    }

    const auto& result = conf.Result<AIMNet2ChargeResponseGradientResult>();
    int written = result.WriteFeatures(conf, output_dir.string());
    EXPECT_EQ(written, 2);

    const fs::path vec_path = output_dir / "aimnet2_charge_response_gradient.npy";
    const fs::path scalar_path =
        output_dir / "aimnet2_charge_response_gradient_scalar.npy";
    EXPECT_TRUE(fs::exists(vec_path)) << "missing " << vec_path;
    EXPECT_TRUE(fs::exists(scalar_path)) << "missing " << scalar_path;
    EXPECT_GT(fs::file_size(vec_path), 0u);
    EXPECT_GT(fs::file_size(scalar_path), 0u);

    // No fs::remove_all here. The PyTorch-shipped libtorch.so exports a
    // broken `std::filesystem::remove_all` symbol that the dynamic
    // linker resolves before libstdc++ (because LD_LIBRARY_PATH lists
    // torch's lib dir first); it dispatches through a NULL pointer and
    // crashes. Tests that clean their own /tmp artifacts must avoid
    // `fs::remove_all` and use unlink() / rmdir() / `std::system("rm
    // -rf …")` instead. The artifacts are kilobytes; /tmp handles
    // them. This is a libtorch packaging defect, not a test bug.
}
