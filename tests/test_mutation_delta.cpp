#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <cmath>
#include <iostream>

#include "MutationDeltaResult.h"
#include "OrcaRunLoader.h"
#include "OrcaShieldingResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "GeometryResult.h"
#include "EnrichmentResult.h"
#include "SpatialIndexResult.h"
#include "DsspResult.h"
#include "ApbsFieldResult.h"
#include "MopacResult.h"
#include "MolecularGraphResult.h"
#include "ProtonationDetectionResult.h"

namespace fs = std::filesystem;
using namespace nmr;



// ============================================================================
// Helper: load an ORCA protein and run the full Layer 0 pipeline
// ============================================================================

struct LoadedProtein {
    std::unique_ptr<Protein> protein;
    bool ok = false;
};

static LoadedProtein LoadAndPrepare(const std::string& prefix) {
    OrcaRunFiles files;
    files.pdb_path    = std::string(nmr::test::TestEnvironment::OrcaDir()) + prefix + ".pdb";
    files.xyz_path    = std::string(nmr::test::TestEnvironment::OrcaDir()) + prefix + ".xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + prefix + ".prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        return {nullptr, false};

    auto load = BuildFromOrca(files);
    if (!load.Ok()) return {nullptr, false};

    auto& conf = load.protein->Conformation();

    // Full Layer 0 pipeline
    conf.AttachResult(ProtonationDetectionResult::Compute(conf));
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(EnrichmentResult::Compute(conf));

    PrmtopChargeSource charge_source(files.prmtop_path);
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));

    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto dssp = DsspResult::Compute(conf);
    if (dssp) conf.AttachResult(std::move(dssp));

    conf.AttachResult(ApbsFieldResult::Compute(conf));
    conf.AttachResult(MolecularGraphResult::Compute(conf));

    // ORCA shielding tensors
    std::string nmr_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + prefix + "_nmr.out";
    if (fs::exists(nmr_path)) {
        auto orca = OrcaShieldingResult::Compute(conf, nmr_path);
        if (orca) conf.AttachResult(std::move(orca));
    }

    return {std::move(load.protein), true};
}


// ============================================================================
// Test fixture
// ============================================================================

class MutationDeltaTest : public ::testing::Test {
protected:
    void SetUp() override {
        wt_ = LoadAndPrepare("A0A7C5FAR6_WT");
        ala_ = LoadAndPrepare("A0A7C5FAR6_ALA");

        if (!wt_.ok || !ala_.ok)
            GTEST_SKIP() << "ORCA test data not available";

        auto& wt_conf = wt_.protein->Conformation();
        auto& ala_conf = ala_.protein->Conformation();

        if (!wt_conf.HasResult<OrcaShieldingResult>() ||
            !ala_conf.HasResult<OrcaShieldingResult>())
            GTEST_SKIP() << "ORCA shielding tensors not loaded";
    }

    LoadedProtein wt_;
    LoadedProtein ala_;
};


// ============================================================================
// Tests
// ============================================================================

TEST_F(MutationDeltaTest, ComputeSucceeds) {
    auto& wt_conf = wt_.protein->Conformation();
    const auto& ala_conf = ala_.protein->Conformation();

    auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
    ASSERT_NE(delta, nullptr);
    ASSERT_TRUE(wt_conf.AttachResult(std::move(delta)));
}


TEST_F(MutationDeltaTest, FourMutationSitesDetected) {
    auto& wt_conf = wt_.protein->Conformation();
    const auto& ala_conf = ala_.protein->Conformation();

    auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
    ASSERT_NE(delta, nullptr);

    const auto& sites = delta->MutationSites();
    EXPECT_EQ(sites.size(), 4u);

    for (const auto& site : sites) {
        EXPECT_EQ(site.mut_type, AminoAcid::ALA);
        std::cout << "  Mutation site " << site.residue_index << ": "
                  << ThreeLetterCodeForAminoAcid(site.wt_type) << " -> "
                  << ThreeLetterCodeForAminoAcid(site.mut_type)
                  << " (" << site.wt_ring_indices.size() << " rings)\n";
    }

    // Check all aromatic types present
    std::set<AminoAcid> wt_types;
    for (const auto& s : sites) wt_types.insert(s.wt_type);
    EXPECT_TRUE(wt_types.count(AminoAcid::TRP) > 0);
    EXPECT_TRUE(wt_types.count(AminoAcid::TYR) > 0);
    EXPECT_TRUE(wt_types.count(AminoAcid::HIS) > 0);
    EXPECT_TRUE(wt_types.count(AminoAcid::PHE) > 0);

    // TRP should have 3 rings (benzene + pyrrole + perimeter), others 1
    for (const auto& s : sites) {
        if (s.wt_type == AminoAcid::TRP)
            EXPECT_EQ(s.wt_ring_indices.size(), 3u);
        else
            EXPECT_EQ(s.wt_ring_indices.size(), 1u);
    }

    wt_conf.AttachResult(std::move(delta));
}


TEST_F(MutationDeltaTest, AtomMatchingReasonable) {
    auto& wt_conf = wt_.protein->Conformation();
    const auto& ala_conf = ala_.protein->Conformation();

    auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
    ASSERT_NE(delta, nullptr);

    std::cout << "  WT=" << wt_conf.AtomCount()
              << " ALA=" << ala_conf.AtomCount()
              << " matched=" << delta->MatchedAtomCount()
              << " unmatched=" << delta->UnmatchedWtAtomCount() << "\n";

    EXPECT_GT(delta->MatchedAtomCount(), 400u);
    EXPECT_GT(delta->UnmatchedWtAtomCount(), 0u);
    EXPECT_EQ(delta->MatchedAtomCount() + delta->UnmatchedWtAtomCount(),
              wt_conf.AtomCount());
}


TEST_F(MutationDeltaTest, BackboneAtomsMatch) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    const Protein& p = wt_conf.ProteinRef();
    size_t res0_N = p.ResidueAt(0).N;
    if (res0_N != Residue::NONE) {
        EXPECT_TRUE(delta->HasMatch(res0_N));
        EXPECT_LT(delta->MatchDistanceAt(res0_N), 0.1);
    }
}


TEST_F(MutationDeltaTest, DeltaShieldingFinite) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    int checked = 0;
    for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
        if (!delta->HasMatch(ai)) continue;
        const auto& st = delta->DeltaShieldingSphericalAt(ai);
        EXPECT_FALSE(std::isnan(st.T0)) << "NaN T0 at atom " << ai;
        for (int i = 0; i < 5; ++i)
            EXPECT_FALSE(std::isnan(st.T2[i])) << "NaN T2[" << i << "] at " << ai;
        checked++;
    }
    EXPECT_GT(checked, 400);
}


TEST_F(MutationDeltaTest, RingProximityComputed) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    // Count total removed rings
    size_t total_rings = 0;
    for (const auto& s : delta->MutationSites())
        total_rings += s.wt_ring_indices.size();

    int has_proximity = 0;
    for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
        if (!delta->HasMatch(ai)) continue;
        const auto& m = delta->MatchedDataAt(ai);
        if (!m.removed_ring_proximity.empty()) {
            // Each matched atom should have proximity to ALL removed rings
            EXPECT_EQ(m.removed_ring_proximity.size(), total_rings);

            // Nearest removed ring distance should be positive
            EXPECT_GT(m.nearest_removed_ring_dist, 0.0);
            EXPECT_LT(m.nearest_removed_ring_dist, 50.0);

            // Cylindrical coords should be consistent
            for (const auto& rp : m.removed_ring_proximity) {
                EXPECT_GT(rp.distance, 0.0);
                double r_from_cyl = std::sqrt(rp.z * rp.z + rp.rho * rp.rho);
                EXPECT_NEAR(r_from_cyl, rp.distance, 0.01)
                    << "Cylindrical coords inconsistent at atom " << ai;
            }
            has_proximity++;
        }
    }
    EXPECT_GT(has_proximity, 400) << "Most matched atoms should have ring proximity";
}


TEST_F(MutationDeltaTest, DistanceDecayCurve) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    const auto& summary = delta->Summary();

    std::cout << "  Distance decay of |delta T0|:\n";
    for (const auto& bin : summary.by_distance) {
        if (bin.count > 0) {
            std::cout << "    " << (int)bin.bin_start << "-"
                      << (int)bin.bin_end << " A: n="
                      << bin.count << " mean_|dT0|="
                      << bin.mean_abs_delta_t0
                      << " mean_|T2|=" << bin.mean_t2_magnitude << "\n";
        }
    }

    // Signal should be stronger near rings than far away
    // Find the 3-4A bin and the 10-11A bin
    double near_signal = 0, far_signal = 0;
    for (const auto& bin : summary.by_distance) {
        if (bin.bin_start == 3.0 && bin.count > 0) near_signal = bin.mean_abs_delta_t0;
        if (bin.bin_start == 10.0 && bin.count > 0) far_signal = bin.mean_abs_delta_t0;
    }
    if (near_signal > 0 && far_signal > 0) {
        EXPECT_GT(near_signal, far_signal)
            << "Signal should decay with distance from removed rings";
    }
}


TEST_F(MutationDeltaTest, ElementStratification) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    const auto& summary = delta->Summary();

    std::cout << "  Delta by element:\n";
    for (const auto& bin : summary.by_element) {
        std::cout << "    " << SymbolForElement(bin.element)
                  << ": n=" << bin.count
                  << " mean_dT0=" << bin.mean_delta_t0
                  << " mean_|dT0|=" << bin.mean_abs_delta_t0
                  << " max_|dT0|=" << bin.max_abs_delta_t0
                  << " mean_|T2|=" << bin.mean_t2_magnitude << "\n";
    }

    EXPECT_FALSE(summary.by_element.empty());

    // Backbone vs sidechain
    std::cout << "  Backbone: n=" << summary.backbone_count
              << " mean_|dT0|=" << summary.backbone_mean_abs_t0 << "\n";
    std::cout << "  Sidechain: n=" << summary.sidechain_count
              << " mean_|dT0|=" << summary.sidechain_mean_abs_t0 << "\n";
}


TEST_F(MutationDeltaTest, DsspDeltaAvailable) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    EXPECT_TRUE(delta->HasDsspDelta());

    if (delta->HasDsspDelta()) {
        // SASA should change at mutation sites (exposed surface changes)
        double max_sasa_delta = 0;
        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            if (!delta->HasMatch(ai)) continue;
            max_sasa_delta = std::max(max_sasa_delta,
                std::abs(delta->MatchedDataAt(ai).delta_sasa));
        }
        std::cout << "  Max |SASA delta|: " << max_sasa_delta << " A^2\n";
        EXPECT_GT(max_sasa_delta, 0.0) << "Some SASA should change";
    }
}


TEST_F(MutationDeltaTest, GraphDeltaAvailable) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    EXPECT_TRUE(delta->HasGraphDelta());

    if (delta->HasGraphDelta()) {
        // Ring atoms in WT had graph_dist_ring=0; in mutant those rings
        // don't exist, so graph_dist_ring is large. Delta should be negative
        // (WT was closer to rings).
        int large_delta = 0;
        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            if (!delta->HasMatch(ai)) continue;
            int dd = delta->MatchedDataAt(ai).delta_graph_dist_ring;
            if (std::abs(dd) > 3) large_delta++;
        }
        std::cout << "  Atoms with |graph delta| > 3: " << large_delta << "\n";
    }
}


TEST_F(MutationDeltaTest, AccessViaTemplateAfterAttach) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);
    ASSERT_TRUE(wt_conf.AttachResult(std::move(delta)));

    ASSERT_TRUE(wt_conf.HasResult<MutationDeltaResult>());
    const auto& d = wt_conf.Result<MutationDeltaResult>();
    EXPECT_GT(d.MatchedAtomCount(), 0u);
    EXPECT_EQ(d.MutationSites().size(), 4u);
    EXPECT_FALSE(d.Summary().by_element.empty());
}
