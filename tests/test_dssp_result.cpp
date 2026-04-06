#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "DsspResult.h"
#include "PhysicalConstants.h"
#include <filesystem>
#include <cmath>

using namespace nmr;


class DsspResultTest : public ::testing::Test {
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

TEST_F(DsspResultTest, ComputeSucceeds) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));
}

TEST_F(DsspResultTest, HasSecondaryStructure) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    // Check that we have per-residue data
    EXPECT_EQ(d.AllResidues().size(), protein->ResidueCount());
}

TEST_F(DsspResultTest, HelixResiduesDetected) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    // 1UBQ has an alpha helix from residue 23 to 34 (PDB numbering).
    // Find the residue indices corresponding to these sequence numbers.
    int helix_count = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        int seq = protein->ResidueAt(ri).sequence_number;
        if (seq >= 23 && seq <= 34) {
            char ss = d.SecondaryStructure(ri);
            // H = alpha helix, G = 3-10 helix, I = pi helix
            if (ss == 'H' || ss == 'G' || ss == 'I')
                helix_count++;
        }
    }
    // At least some residues in this range should be helical
    EXPECT_GT(helix_count, 5) << "Expected helix residues in range 23-34";
}

TEST_F(DsspResultTest, SheetResiduesDetected) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    // 1UBQ has beta strands. Check for at least some E/B assignments.
    int sheet_count = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        char ss = d.SecondaryStructure(ri);
        if (ss == 'E' || ss == 'B')
            sheet_count++;
    }
    EXPECT_GT(sheet_count, 5) << "Expected some beta strand residues";
}

TEST_F(DsspResultTest, PhiPsiInRange) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    // Phi and Psi should be in [-pi, pi] for internal residues.
    // Terminal residues may have special values (360 from DSSP = 2*pi).
    int valid_count = 0;
    for (size_t ri = 1; ri < protein->ResidueCount() - 1; ++ri) {
        double phi = d.Phi(ri);
        double psi = d.Psi(ri);
        // Allow some tolerance for DSSP conversion: [-pi-0.1, pi+0.1]
        // DSSP returns 360 for undefined angles, which converts to ~6.28 rad
        if (std::abs(phi) <= PI + 0.1 && std::abs(psi) <= PI + 0.1)
            valid_count++;
    }
    // Most internal residues should have valid phi/psi
    EXPECT_GT(valid_count, static_cast<int>(protein->ResidueCount()) / 2);
}

TEST_F(DsspResultTest, SASANonNegative) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        EXPECT_GE(d.SASA(ri), 0.0) << "Negative SASA at residue " << ri;
    }
}

TEST_F(DsspResultTest, SomeSASANonZero) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    int nonzero = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        if (d.SASA(ri) > 0.0) nonzero++;
    }
    // Surface residues should have nonzero SASA
    EXPECT_GT(nonzero, 10);
}
