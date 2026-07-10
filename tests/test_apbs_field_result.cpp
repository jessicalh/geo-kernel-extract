#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "ChargeAssignmentResult.h"
#include "ApbsFieldResult.h"
#include "CalculatorConfig.h"
#include "PhysicalConstants.h"
#include <algorithm>
#include <filesystem>
#include <cmath>
#include <fstream>
#include <iterator>
#include <unistd.h>

using namespace nmr;

namespace {

void WriteApbsTestConfig(const std::filesystem::path& path,
                         double pdie,
                         double sdie,
                         double ionic_strength,
                         double clamp) {
    std::ofstream out(path);
    // 161 is the production grid and avoids APBS's coarse-grid direct-solver
    // fallback, whose legacy teardown is not safe at 65^3 for this protein.
    out << "apbs_grid_dim = 161\n"
        << "apbs_manual_grid_padding_A = 70.0\n"
        << "apbs_manual_grid_min_dim_A = 70.0\n"
        << "apbs_protein_dielectric = " << pdie << "\n"
        << "apbs_solvent_dielectric = " << sdie << "\n"
        << "apbs_temperature_K = 298.15\n"
        << "apbs_ionic_strength_M = " << ionic_strength << "\n"
        << "efield_magnitude_sanity_clamp = " << clamp << "\n";
}

std::unique_ptr<Protein> BuildApbsChildProtein() {
    auto build = BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    if (!build.Ok()) return nullptr;
    return std::move(build.protein);
}

void RunHomogeneousReferenceForcingFunctionInChild() {
    const auto config_path = std::filesystem::temp_directory_path() /
        ("apbs_homogeneous_reference_" + std::to_string(::getpid()) +
         ".toml");
    WriteApbsTestConfig(config_path, 1.0, 1.0, 0.0, 100.0);
    CalculatorConfig::Load(config_path.string());

    auto protein = BuildApbsChildProtein();
    bool ok = protein != nullptr;
    if (ok) {
        auto& conf = protein->Conformation();
        auto charges = ChargeAssignmentResult::Compute(
            conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
        ok = charges != nullptr && conf.AttachResult(std::move(charges));
        auto apbs = ok ? ApbsFieldResult::Compute(conf) : nullptr;
        ok = ok && apbs != nullptr;

        double max_reaction = 0.0;
        double max_total_diagnostic = 0.0;
        if (ok) {
            for (size_t i = 0; i < conf.AtomCount(); ++i) {
                const auto& atom = conf.AtomAt(i);
                max_reaction = std::max(max_reaction,
                    atom.apbs_efield.norm() + atom.apbs_efg.norm() +
                    std::abs(atom.apbs_phi));
                max_total_diagnostic = std::max(max_total_diagnostic,
                    atom.apbs_efield_total_diagnostic.norm() +
                    atom.apbs_efg_total_diagnostic.norm());
                ok = ok && atom.apbs_efield_clamp_mask == 0u
                        && atom.apbs_efield_clamp_scale == 1.0;
            }
            // With production physics set equal to the independently defined
            // homogeneous reference, the two deterministic grids must cancel.
            ok = ok && max_reaction < 1e-8
                    && max_total_diagnostic > 1e-8;
        }
    }

    std::filesystem::remove(config_path);
    _exit(ok ? 0 : 1);
}

void RunPostConversionClampForcingFunctionInChild() {
    constexpr double kClampVPerA = 1e-8;
    const auto config_path = std::filesystem::temp_directory_path() /
        ("apbs_post_conversion_clamp_" + std::to_string(::getpid()) +
         ".toml");
    WriteApbsTestConfig(config_path, 4.0, 78.54, 0.15, kClampVPerA);
    CalculatorConfig::Load(config_path.string());

    auto protein = BuildApbsChildProtein();
    bool ok = protein != nullptr;
    std::size_t n_clamped = 0;
    bool total_diagnostic_exceeds_clamp = false;
    if (ok) {
        auto& conf = protein->Conformation();
        auto charges = ChargeAssignmentResult::Compute(
            conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
        ok = charges != nullptr && conf.AttachResult(std::move(charges));
        auto apbs = ok ? ApbsFieldResult::Compute(conf) : nullptr;
        ok = ok && apbs != nullptr;
        if (ok) {
            for (size_t i = 0; i < conf.AtomCount(); ++i) {
                const auto& atom = conf.AtomAt(i);
                total_diagnostic_exceeds_clamp =
                    total_diagnostic_exceeds_clamp ||
                    atom.apbs_efield_total_diagnostic.norm() > kClampVPerA;
                if (atom.apbs_efield_clamp_mask != 0u) {
                    ++n_clamped;
                    // This distinguishes a V/A clamp after conversion from
                    // the rejected native-kT/e clamp (which would be smaller
                    // by the thermal-voltage factor after conversion).
                    ok = ok && std::abs(atom.apbs_efield.norm() -
                                         kClampVPerA) < 1e-12
                            && atom.apbs_efield_clamp_scale > 0.0
                            && atom.apbs_efield_clamp_scale < 1.0;
                } else {
                    ok = ok && atom.apbs_efield_clamp_scale == 1.0
                            && atom.apbs_efield.norm() <= kClampVPerA;
                }
            }
            ok = ok && n_clamped > 0 && total_diagnostic_exceeds_clamp;
        }
    }

    std::filesystem::remove(config_path);
    _exit(ok ? 0 : 1);
}

}  // namespace


TEST(ApbsFieldResultStandalone, ThermalVoltageUsesConfiguredTemperature) {
    EXPECT_DOUBLE_EQ(ApbsThermalVoltage(298.15), KT_OVER_E_298K);
    EXPECT_DOUBLE_EQ(ApbsThermalVoltage(310.0),
                     KT_OVER_E_298K * 310.0 / 298.15);
}

TEST(ApbsFieldResultStandalone, HomogeneousReferenceCancelsCanonicalField) {
    if (!std::filesystem::exists(
            nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    EXPECT_EXIT({ RunHomogeneousReferenceForcingFunctionInChild(); },
                ::testing::ExitedWithCode(0), "");
}

TEST(ApbsFieldResultStandalone, ClampIsAppliedInVoltsPerAngstrom) {
    if (!std::filesystem::exists(
            nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    EXPECT_EXIT({ RunPostConversionClampForcingFunctionInChild(); },
                ::testing::ExitedWithCode(0), "");
}


class ApbsFieldResultTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};

TEST_F(ApbsFieldResultTest, ComputeSucceeds) {
    auto& conf = protein->Conformation();

    // Attach ChargeAssignmentResult (dependency)
    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    ASSERT_NE(charges, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    // Compute and attach ApbsFieldResult
    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(apbs)));

    const auto& result = conf.Result<ApbsFieldResult>();
    const auto& dims = result.GridDims();
    const auto& lengths = result.GridLengthsA();
    const auto& origin = result.GridOriginA();
    const auto& spacing = result.GridSpacingA();
    Vec3 bbox_min = conf.PositionAt(0);
    Vec3 bbox_max = bbox_min;
    for (size_t i = 1; i < conf.AtomCount(); ++i) {
        bbox_min = bbox_min.cwiseMin(conf.PositionAt(i));
        bbox_max = bbox_max.cwiseMax(conf.PositionAt(i));
    }
    const Vec3 extent = bbox_max - bbox_min;
    const auto configured_dim = static_cast<std::uint64_t>(
        CalculatorConfig::Get("apbs_grid_dim"));
    for (size_t d = 0; d < 3; ++d) {
        EXPECT_EQ(dims[d], configured_dim);
        const double expected_length = std::max(
            extent(static_cast<Eigen::Index>(d)) +
                CalculatorConfig::Get("apbs_manual_grid_padding_A"),
            CalculatorConfig::Get("apbs_manual_grid_min_dim_A"));
        EXPECT_DOUBLE_EQ(lengths[d], expected_length);
        EXPECT_TRUE(std::isfinite(origin[d]));
        EXPECT_NEAR(spacing[d] * static_cast<double>(dims[d] - 1),
                    lengths[d], 1e-12);
    }
    EXPECT_DOUBLE_EQ(result.ManualGridPaddingA(),
                     CalculatorConfig::Get("apbs_manual_grid_padding_A"));
    EXPECT_DOUBLE_EQ(result.ManualGridMinDimA(),
                     CalculatorConfig::Get("apbs_manual_grid_min_dim_A"));
    EXPECT_DOUBLE_EQ(result.TemperatureK(),
                     CalculatorConfig::Get("apbs_temperature_K"));
    EXPECT_DOUBLE_EQ(result.ThermalVoltageV(),
                     ApbsThermalVoltage(result.TemperatureK()));
    EXPECT_DOUBLE_EQ(result.EfieldClampThreshold(),
                     CalculatorConfig::Get("efield_magnitude_sanity_clamp"));

    const auto output_dir = std::filesystem::temp_directory_path() /
        ("apbs_static_schema_" + std::to_string(::getpid()));
    std::filesystem::create_directories(output_dir);
    EXPECT_EQ(result.WriteFeatures(conf, output_dir.string()), 7);
    for (const char* filename : {
            "apbs_E.npy", "apbs_efg.npy", "apbs_phi.npy",
            "apbs_E_clamp_mask.npy", "apbs_E_clamp_scale.npy",
            "apbs_E_total_diagnostic.npy",
            "apbs_efg_total_diagnostic.npy"}) {
        EXPECT_TRUE(std::filesystem::exists(output_dir / filename))
            << filename;
    }
    std::ifstream mask_file(output_dir / "apbs_E_clamp_mask.npy",
                            std::ios::binary);
    const std::string mask_bytes(
        (std::istreambuf_iterator<char>(mask_file)),
        std::istreambuf_iterator<char>());
    EXPECT_NE(mask_bytes.find("|u1"), std::string::npos);
    mask_file.close();
    for (const char* filename : {
            "apbs_E.npy", "apbs_efg.npy", "apbs_phi.npy",
            "apbs_E_clamp_mask.npy", "apbs_E_clamp_scale.npy",
            "apbs_E_total_diagnostic.npy",
            "apbs_efg_total_diagnostic.npy"}) {
        std::filesystem::remove(output_dir / filename);
    }
    std::filesystem::remove(output_dir);
}

TEST_F(ApbsFieldResultTest, EFieldNonZeroForSomeAtoms) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    int nonzero = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Vec3 E = result.ElectricFieldAt(ai);
        if (E.norm() > 1e-10) nonzero++;
    }
    // A normal heterogeneous PB solve must retain a nonzero reaction field
    // for a substantial fraction of atoms.
    EXPECT_GT(nonzero, static_cast<int>(conf.AtomCount()) / 2);
}

TEST_F(ApbsFieldResultTest, NoNanOrInf) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Vec3 E = result.ElectricFieldAt(ai);
        Mat3 EFG = result.FieldGradientAt(ai);
        const auto& atom = conf.AtomAt(ai);

        for (int d = 0; d < 3; ++d) {
            EXPECT_FALSE(std::isnan(E(d))) << "NaN in E-field, atom " << ai;
            EXPECT_FALSE(std::isinf(E(d))) << "Inf in E-field, atom " << ai;
        }

        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                EXPECT_FALSE(std::isnan(EFG(a,b)))
                    << "NaN in EFG, atom " << ai << " (" << a << "," << b << ")";
                EXPECT_FALSE(std::isinf(EFG(a,b)))
                    << "Inf in EFG, atom " << ai << " (" << a << "," << b << ")";
            }
        }
        EXPECT_TRUE(std::isfinite(atom.apbs_phi));
        EXPECT_TRUE(atom.apbs_efield_total_diagnostic.allFinite());
        EXPECT_TRUE(atom.apbs_efg_total_diagnostic.allFinite());
        EXPECT_TRUE(atom.apbs_efield_clamp_mask == 0u ||
                    atom.apbs_efield_clamp_mask == 1u);
        if (atom.apbs_efield_clamp_mask == 0u)
            EXPECT_DOUBLE_EQ(atom.apbs_efield_clamp_scale, 1.0);
        else {
            EXPECT_GT(atom.apbs_efield_clamp_scale, 0.0);
            EXPECT_LT(atom.apbs_efield_clamp_scale, 1.0);
        }
    }
}

TEST_F(ApbsFieldResultTest, EFGDecomposesCorrectly) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    // SphericalTensor roundtrip: decompose and reconstruct should match
    // Test on a few atoms
    for (size_t ai = 0; ai < std::min(conf.AtomCount(), size_t(20)); ++ai) {
        Mat3 EFG = result.FieldGradientAt(ai);
        SphericalTensor st = result.FieldGradientSphericalAt(ai);

        // Reconstruct from spherical
        Mat3 reconstructed = st.Reconstruct();

        // Should match within numerical tolerance
        double diff = (EFG - reconstructed).norm();
        EXPECT_LT(diff, 1e-10)
            << "SphericalTensor roundtrip failed at atom " << ai
            << " diff=" << diff;
    }
}

TEST_F(ApbsFieldResultTest, EFGIsTraceless) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    // EFG tensor should be traceless (trace projected out)
    for (size_t ai = 0; ai < std::min(conf.AtomCount(), size_t(50)); ++ai) {
        Mat3 EFG = result.FieldGradientAt(ai);
        double trace = EFG.trace();
        EXPECT_NEAR(trace, 0.0, 1e-10)
            << "Non-zero trace at atom " << ai << " trace=" << trace;
    }
}

TEST_F(ApbsFieldResultTest, FullTensorStored) {
    // Both Mat3 AND SphericalTensor must be stored on ConformationAtom
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    // Check that ConformationAtom fields are populated
    for (size_t ai = 0; ai < std::min(conf.AtomCount(), size_t(10)); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        EXPECT_EQ(ca.apbs_nonfinite_sanitizer_mask, 0u);

        // At least one of these should be nonzero in the normal PB solve.
        bool efield_ok = ca.apbs_efield.norm() > 1e-15;

        if (efield_ok) {
            // SphericalTensor should also be set
            SphericalTensor st = ca.apbs_efg_spherical;
            Mat3 roundtrip = st.Reconstruct();
            double diff = (ca.apbs_efg - roundtrip).norm();
            EXPECT_LT(diff, 1e-10);

            const Mat3 total_roundtrip =
                ca.apbs_efg_total_diagnostic_spherical.Reconstruct();
            EXPECT_LT((ca.apbs_efg_total_diagnostic - total_roundtrip).norm(),
                      1e-10);
        }
    }
}

TEST_F(ApbsFieldResultTest, DependencyEnforced) {
    auto& conf = protein->Conformation();

    // Try to attach ApbsFieldResult WITHOUT ChargeAssignmentResult
    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);

    // AttachResult should fail because ChargeAssignmentResult is missing
    EXPECT_FALSE(conf.AttachResult(std::move(apbs)));
}
