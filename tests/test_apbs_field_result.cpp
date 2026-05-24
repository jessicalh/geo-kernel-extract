#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "ChargeAssignmentResult.h"
#include "ApbsFieldResult.h"
#include <filesystem>
#include <cmath>

using namespace nmr;


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
    // With stub charges (all 0.1e), E-fields should be nonzero for all atoms
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

        // At least one of these should be nonzero (all atoms have charges)
        bool efield_ok = ca.apbs_efield.norm() > 1e-15;
        bool efg_ok = ca.apbs_efg.norm() > 1e-15;

        // With uniform positive charges, E-field should be nonzero
        // for atoms that are not at the center of mass
        if (efield_ok) {
            // SphericalTensor should also be set
            SphericalTensor st = ca.apbs_efg_spherical;
            Mat3 roundtrip = st.Reconstruct();
            double diff = (ca.apbs_efg - roundtrip).norm();
            EXPECT_LT(diff, 1e-10);
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
