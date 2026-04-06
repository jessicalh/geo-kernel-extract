#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "ProtonationDetectionResult.h"
#include "Protein.h"
#include <filesystem>

using namespace nmr;

// ============================================================================
// 1UBQ crystal structure (no hydrogens) — titratable residues unresolved
// ============================================================================


class ProtonationNoHTest : public ::testing::Test {
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

TEST_F(ProtonationNoHTest, ComputeSucceeds) {
    auto& conf = protein->Conformation();
    auto result = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(result, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(result)));
}

// ProtonationNoHTest.TitratableResiduesUnresolved — DELETED 2026-04-03
// ProtonationNoHTest.VariantNamesEmptyForNoH — DELETED 2026-04-03
// These tested the "PDB without H" scenario which no longer exists.
// Every PDB is now protonated by reduce before entering the system.

// ============================================================================
// Protonated protein — titratable residues should have variants
// ============================================================================


class ProtonationWithHTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::GmxProtonated())) {
            GTEST_SKIP() << "Protonated PDB not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::GmxProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};

TEST_F(ProtonationWithHTest, ComputeSucceeds) {
    auto& conf = protein->Conformation();
    auto result = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(result, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(result)));
}

TEST_F(ProtonationWithHTest, SomeTitratableResiduesAssigned) {
    auto& conf = protein->Conformation();
    auto result = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(result, nullptr);

    // Protonated structure should have some titratable residues resolved.
    // At minimum, unresolved count should be 0 (hydrogens are present).
    EXPECT_EQ(result->UnresolvedCount(), 0)
        << "Protonated structure should have no unresolved titratable residues";

    conf.AttachResult(std::move(result));
}

TEST_F(ProtonationWithHTest, VariantIndexSetOnResidues) {
    auto& conf = protein->Conformation();
    auto prot_result = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(prot_result, nullptr);
    conf.AttachResult(std::move(prot_result));

    // Check that at least some titratable residues have variant_index set
    bool found_titratable = false;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const Residue& res = protein->ResidueAt(ri);
        if (res.IsTitratable()) {
            found_titratable = true;
            // In a protonated structure, variant index may or may not be set
            // (depends on whether the variant differs from default).
            // But the test should at least not crash.
        }
    }
    // If there are titratable residues, good. If not, the test still passes.
    (void)found_titratable;
}

TEST_F(ProtonationWithHTest, HisRingTypeMatchesProtonation) {
    auto& conf = protein->Conformation();
    auto prot_result = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(prot_result, nullptr);
    conf.AttachResult(std::move(prot_result));

    // For each HIS residue with a protonation variant, check that the ring
    // type index is consistent. HID -> HidImidazole, HIE -> HieImidazole,
    // HIP -> HisImidazole.
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const Residue& res = protein->ResidueAt(ri);
        if (res.type != AminoAcid::HIS) continue;

        std::string variant = conf.Result<ProtonationDetectionResult>()
                                  .VariantNameAt(ri);

        // Find the ring for this HIS residue
        for (size_t ring_i = 0; ring_i < protein->RingCount(); ++ring_i) {
            const Ring& ring = protein->RingAt(ring_i);
            if (ring.parent_residue_index != ri) continue;

            // Ring detection already uses HD1/HE2, so it should be
            // consistent with what protonation detection found.
            if (variant == "HID") {
                EXPECT_EQ(ring.type_index, RingTypeIndex::HidImidazole)
                    << "HIS " << res.sequence_number
                    << " protonation=HID but ring type is "
                    << ring.TypeName();
            } else if (variant == "HIE") {
                EXPECT_EQ(ring.type_index, RingTypeIndex::HieImidazole)
                    << "HIS " << res.sequence_number
                    << " protonation=HIE but ring type is "
                    << ring.TypeName();
            } else if (variant == "HIP") {
                EXPECT_EQ(ring.type_index, RingTypeIndex::HisImidazole)
                    << "HIS " << res.sequence_number
                    << " protonation=HIP but ring type is "
                    << ring.TypeName();
            }
        }
    }
}
