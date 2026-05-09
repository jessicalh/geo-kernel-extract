// Smoke test for AIMNet2PolarisabilityResult.
//
// Mirrors the test_apbs_ff14sb.cpp pattern: load 1UBQ via
// BuildFromProtonatedPdb, attach the dependency chain, run
// AIMNet2PolarisabilityResult and verify gradients + WriteFeatures.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AIMNet2PolarisabilityResult.h"
#include "AIMNet2Result.h"
#include "ChargeAssignmentResult.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "SpatialIndexResult.h"

#include <torch/cuda.h>

#include <cmath>
#include <cstdio>
#include <filesystem>

using namespace nmr;
namespace fs = std::filesystem;


class AIMNet2PolarisabilityTest : public ::testing::Test {
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


TEST_F(AIMNet2PolarisabilityTest, PipelineProducesNonZeroPolarisability) {
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
    ASSERT_NE(charges, nullptr) << "ff14SB charge loading failed";
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    auto aim = AIMNet2Result::Compute(conf, *model);
    ASSERT_NE(aim, nullptr) << "AIMNet2Result::Compute returned nullptr";
    ASSERT_TRUE(conf.AttachResult(std::move(aim)));

    auto pol = AIMNet2PolarisabilityResult::Compute(conf, *model);
    ASSERT_NE(pol, nullptr) << "AIMNet2PolarisabilityResult returned nullptr";
    ASSERT_TRUE(conf.AttachResult(std::move(pol)));

    const size_t N = conf.AtomCount();
    ASSERT_GT(N, 0u);

    int finite_count = 0;
    int nonzero_count = 0;
    double max_scalar = 0.0;
    double sum_scalar = 0.0;
    double max_consistency_diff = 0.0;

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        const Vec3& v = ca.aimnet2_polarisability_vector;
        const double s = ca.aimnet2_polarisability_scalar;

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
        "\n=== AIMNet2 polarisability summary (1UBQ, %zu atoms) ===\n"
        "  finite values:         %d / %zu\n"
        "  non-zero scalar:       %d / %zu\n"
        "  max scalar:            %.6e\n"
        "  mean scalar:           %.6e\n"
        "  max |scalar - ||v|||:  %.6e\n"
        "==========================================================\n",
        N, finite_count, N, nonzero_count, N,
        max_scalar, mean_scalar, max_consistency_diff);

    EXPECT_EQ(finite_count, static_cast<int>(N))
        << "Some polarisability values are NaN or Inf";
    EXPECT_GT(nonzero_count, static_cast<int>(N) / 4)
        << "Too few atoms with non-zero polarisability gradient.";
    EXPECT_GT(max_scalar, 1e-8)
        << "Maximum gradient norm is suspiciously small";
    EXPECT_LT(max_consistency_diff, 1e-9)
        << "Scalar field does not match L2 norm of vector field";
}


TEST_F(AIMNet2PolarisabilityTest, WriteFeaturesEmitsBothNpys) {
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

    auto aim = AIMNet2Result::Compute(conf, *model);
    ASSERT_NE(aim, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(aim)));

    auto pol = AIMNet2PolarisabilityResult::Compute(conf, *model);
    ASSERT_NE(pol, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(pol)));

    const fs::path output_dir = fs::temp_directory_path() /
        "aimnet2_polarisability_test_writefeatures";
    fs::create_directories(output_dir);

    const auto& result = conf.Result<AIMNet2PolarisabilityResult>();
    int written = result.WriteFeatures(conf, output_dir.string());
    EXPECT_EQ(written, 2);

    const fs::path vec_path = output_dir / "aimnet2_polarisability.npy";
    const fs::path scalar_path = output_dir / "aimnet2_polarisability_scalar.npy";
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
