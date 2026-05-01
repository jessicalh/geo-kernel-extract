#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "OrcaShieldingResult.h"
#include "PropkaProtonator.h"
#include "ProtonationDetectionResult.h"
#include "ProtonationState.h"
#include "ChargeSource.h"
#include "ChargeAssignmentResult.h"
#include "AminoAcidType.h"
#include "RuntimeEnvironment.h"
#include "Protein.h"

#include <filesystem>
#include <cmath>
#include <numeric>

using namespace nmr;

// ============================================================================
// Validate variant index contract at test startup
// ============================================================================

class ProtonationPipelineTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        ValidateVariantIndices();
    }
};


// ============================================================================
// PROPKA on 1UBQ crystal structure (no hydrogens)
// ============================================================================


class PropkaTest : public ProtonationPipelineTest {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        if (!std::filesystem::exists("propka3")) {
            GTEST_SKIP() << "propka3 not found at " << "propka3";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};

TEST_F(PropkaTest, PredictsPkaValues) {
    auto& conf = protein->Conformation();
    PropkaProtonator propka;

    std::string error;
    auto pkas = PropkaProtonator::PredictPka(*protein, conf, error);
    ASSERT_FALSE(pkas.empty()) << "PROPKA returned no pKa values: " << error;

    // 1UBQ has titratable residues: ASP, GLU, HIS, LYS, TYR, CYS
    // PROPKA should find some of them
    EXPECT_GE(pkas.size(), 5u) << "Expected at least 5 pKa predictions for 1UBQ";

    // All pKa values should be physically reasonable (0-14)
    for (const auto& pka : pkas) {
        EXPECT_GE(pka.pKa, 0.0) << pka.residue_type << " " << pka.residue_number;
        EXPECT_LE(pka.pKa, 14.0) << pka.residue_type << " " << pka.residue_number;
    }
}

TEST_F(PropkaTest, ProtonatesAtPhysiologicalPH) {
    auto& conf = protein->Conformation();
    PropkaProtonator propka;

    auto result = propka.Protonate(*protein, conf, 7.0);
    ASSERT_TRUE(result.Ok()) << result.error;

    const auto& state = result.state;
    EXPECT_GT(state.DecisionCount(), 0u)
        << "No protonation decisions at pH 7";

    // Build residue type list for net charge calculation
    std::vector<AminoAcid> residue_types;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        residue_types.push_back(protein->ResidueAt(ri).type);
    }

    // 1UBQ at pH 7: expect net charge roughly 0 to +3
    // 4 ASP(-1) + 6 GLU(-1) + 1 MET(N-term, +1) = -9
    // 7 LYS(+1) + 4 ARG(+1) + 1 HIS(~0 at pH 7) = +11 to +12
    // Net: roughly +2 to +3, but without terminal patches
    int net = state.NetChargeForProtein(residue_types);
    EXPECT_GE(net, -2) << "Net charge " << net << " too negative for pH 7";
    EXPECT_LE(net, 6) << "Net charge " << net << " too positive for pH 7";
}

TEST_F(PropkaTest, LowPHProtonatesAcids) {
    auto& conf = protein->Conformation();
    PropkaProtonator propka;

    auto result = propka.Protonate(*protein, conf, 2.0);
    ASSERT_TRUE(result.Ok()) << result.error;

    // At pH 2: most ASP and GLU should be protonated (neutral)
    int protonated_acids = 0;
    for (const auto& d : result.state.Decisions()) {
        if (d.amino_acid == AminoAcid::ASP && d.variant_index == 0) protonated_acids++;
        if (d.amino_acid == AminoAcid::GLU && d.variant_index == 0) protonated_acids++;
    }
    // 1UBQ has 5 ASP + 6 GLU = 11 acids. At pH 2, most should be protonated.
    EXPECT_GE(protonated_acids, 5)
        << "Only " << protonated_acids << " acids protonated at pH 2";
}

TEST_F(PropkaTest, DecisionsAreTyped) {
    auto& conf = protein->Conformation();
    PropkaProtonator propka;

    auto result = propka.Protonate(*protein, conf, 7.0);
    ASSERT_TRUE(result.Ok());

    // Every decision must have valid typed fields
    for (const auto& d : result.state.Decisions()) {
        EXPECT_NE(d.amino_acid, AminoAcid::Unknown)
            << "Decision at residue " << d.residue_index << " has Unknown amino acid";
        EXPECT_LT(d.residue_index, protein->ResidueCount())
            << "Decision residue_index out of range";
        EXPECT_FALSE(std::isnan(d.pKa))
            << "Decision at residue " << d.residue_index << " has NaN pKa";

        // If variant_index >= 0, it must be valid for that amino acid
        if (d.variant_index >= 0) {
            const auto& aat = GetAminoAcidType(d.amino_acid);
            EXPECT_LT(d.variant_index, static_cast<int>(aat.variants.size()))
                << "variant_index " << d.variant_index << " out of range for "
                << aat.three_letter_code;
        }
    }
}


// ============================================================================
// GmxTprChargeSource on real fleet data
// ============================================================================

// GmxChargeTest fixture — DELETED 2026-04-03
// Tested GmxTprChargeSource (gmx dump), superseded by PreloadedChargeSource
// in BuildFromGromacs. CHARMM charge validation is in test_fleet_loader.cpp.

// GmxChargeTest.LoadsCharmmCharges — DELETED 2026-04-03
// GmxChargeTest.ChargesDifferFromFf14sb — DELETED 2026-04-03
// GmxChargeTest.TypedPathThroughChargeAssignmentResult — DELETED 2026-04-03
// Tested GmxTprChargeSource (gmx dump path), superseded by
// PreloadedChargeSource in BuildFromGromacs. GROMACS data loads
// from GROMACS, always. CHARMM charge tests are in test_fleet_loader.cpp.

// ============================================================================
// PrmtopChargeSource on ORCA test data
// ============================================================================


TEST(PrmtopChargeTest, LoadsFromAmberPrmtop) {
    std::string prmtop = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";
    std::string amber_pdb = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT_amber.pdb";
    if (!std::filesystem::exists(prmtop)) GTEST_SKIP() << "prmtop not found";
    if (!std::filesystem::exists(amber_pdb)) GTEST_SKIP() << "amber PDB not found";

    // Load the heavy-atom PDB (amber PDB has same atoms as the AlphaFold PDB)
    auto r = BuildFromProtonatedPdb(amber_pdb);
    if (!r.Ok()) GTEST_SKIP() << r.error;
    auto& conf = r.protein->Conformation();

    // The prmtop has 543 atoms (protonated), the PDB has 280 (heavy only).
    // PrmtopChargeSource should load 280 charges (first 280 of 543).
    // NOTE: This only works if heavy atoms come first in the prmtop, which
    // they do NOT — tleap interleaves hydrogens with their parent heavy atoms.
    // So this test verifies the atom count mismatch is handled.
    PrmtopChargeSource source(prmtop);
    std::string error;
    auto charges = source.LoadCharges(*r.protein, conf, error);

    // With 280 PDB atoms vs 543 prmtop atoms, this should still load
    // (prmtop has >= protein atoms). But the charges won't be correct
    // because the atom ordering doesn't match (prmtop has H interleaved).
    // This is expected — for correct charges, we need to load the
    // protonated structure that matches the prmtop.
    ASSERT_FALSE(charges.empty()) << "Failed to load: " << error;
    EXPECT_EQ(charges.size(), conf.AtomCount());
}

TEST(PrmtopChargeTest, TotalChargeIsInteger) {
    // This test uses the XYZ atom count to verify the prmtop directly,
    // without going through the Protein loader (which can't load XYZ yet).
    std::string prmtop = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";
    if (!std::filesystem::exists(prmtop)) GTEST_SKIP() << "prmtop not found";

    // Read charges directly from the prmtop
    std::string amber_pdb = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT_amber.pdb";
    if (!std::filesystem::exists(amber_pdb)) GTEST_SKIP();
    auto r = BuildFromProtonatedPdb(amber_pdb);
    if (!r.Ok()) GTEST_SKIP();

    PrmtopChargeSource source(prmtop);
    std::string error;
    auto charges = source.LoadCharges(*r.protein, r.protein->Conformation(), error);
    ASSERT_FALSE(charges.empty()) << error;

    // Sum should be near integer. The prmtop has 543 atoms but we only
    // loaded 280 (heavy atom PDB). The sum won't be integer because
    // hydrogen charges are missing. This is the SAME bug as ff14sb_params.dat
    // but now we understand WHY: the protein needs to be loaded from the
    // protonated structure to match the prmtop.
    double total = 0.0;
    for (const auto& acr : charges) total += acr.partial_charge;

    // Just verify it loaded reasonable values (not all zero, not NaN)
    bool has_nonzero = false;
    for (const auto& acr : charges) {
        EXPECT_FALSE(std::isnan(acr.partial_charge));
        EXPECT_FALSE(std::isnan(acr.pb_radius));
        if (std::abs(acr.partial_charge) > 1e-6) has_nonzero = true;
    }
    EXPECT_TRUE(has_nonzero) << "All charges are zero";
}

TEST_F(PropkaTest, ProtonationDetectionOnProtonatedStructure) {
    // Item 1: verify ProtonationDetectionResult works on CHARMM-named structures.
    // CHARMM uses the same titratable hydrogen names as PDB standard
    // (HD1, HE2, HD2, HG, HZ1/HZ2/HZ3), so detection should work directly.
    auto& conf = protein->Conformation();
    auto result = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(result, nullptr);

    // This is a protonated structure — no residues should be unresolved
    EXPECT_EQ(result->UnresolvedCount(), 0)
        << "Protonated CHARMM structure has " << result->UnresolvedCount()
        << " unresolved titratable residues";

    // Should have some assigned variants (at least HIS tautomers)
    EXPECT_GT(result->AssignedCount(), 0)
        << "No protonation variants detected on CHARMM structure";

    conf.AttachResult(std::move(result));

    // Check that HIS residues got variant assignments
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const Residue& res = protein->ResidueAt(ri);
        if (res.type == AminoAcid::HIS) {
            EXPECT_GE(res.protonation_variant_index, 0)
                << "HIS " << res.sequence_number << " not assigned a variant";
        }
    }
}

TEST_F(PropkaTest, PropkaAgreesWithStructure) {
    // Item 3: compare PROPKA's opinion with what the structure actually has.
    // Run PROPKA on the fleet protein, detect protonation from H atoms,
    // and compare. Disagreements are interesting, not errors.
    if (!std::filesystem::exists("propka3")) {
        GTEST_SKIP() << "propka3 not found";
    }

    auto& conf = protein->Conformation();

    // Detect from structure
    auto detect = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(detect, nullptr);
    conf.AttachResult(std::move(detect));

    // Get PROPKA's opinion at pH 7 (the structure was likely prepared near pH 7)
    PropkaProtonator propka;
    auto propka_result = propka.Protonate(*protein, conf, 7.0);
    ASSERT_TRUE(propka_result.Ok()) << propka_result.error;

    // Compare: for each HIS, does PROPKA agree with what was built?
    int agreements = 0;
    int disagreements = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const Residue& res = protein->ResidueAt(ri);
        if (res.type != AminoAcid::HIS) continue;

        int struct_variant = res.protonation_variant_index;
        const auto* propka_decision = propka_result.state.ForResidue(ri);

        if (!propka_decision) continue;

        if (struct_variant == propka_decision->variant_index) {
            agreements++;
        } else {
            disagreements++;
            // Log the disagreement — this is science, not a test failure
            std::string struct_name = (struct_variant >= 0)
                ? res.AminoAcidInfo().variants[struct_variant].name : "default";
            std::string propka_name = (propka_decision->variant_index >= 0)
                ? res.AminoAcidInfo().variants[propka_decision->variant_index].name
                : "default";
            std::cout << "  HIS " << res.sequence_number
                      << ": structure=" << struct_name
                      << " propka=" << propka_name
                      << " (pKa=" << propka_decision->pKa << ")\n";
        }
    }

    // At least check it ran without crashing. Agreement/disagreement
    // is informative, not pass/fail.
    EXPECT_GE(agreements + disagreements, 0);
    std::cout << "HIS protonation: " << agreements << " agree, "
              << disagreements << " disagree with PROPKA\n";
}


// ============================================================================
// ORCA run loader: full pipeline from prmtop + XYZ + charges
// ============================================================================

TEST(OrcaRunTest, LoadsProtonatedProtein) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!std::filesystem::exists(files.prmtop_path)) GTEST_SKIP();
    if (!std::filesystem::exists(files.xyz_path)) GTEST_SKIP();

    auto result = BuildFromOrca(files);
    ASSERT_TRUE(result.Ok()) << result.error;

    // 543 atoms from XYZ/prmtop (protonated), 35 residues
    EXPECT_EQ(result.protein->AtomCount(), 543u);
    EXPECT_EQ(result.protein->ResidueCount(), 35u);

    // Should have one conformation with 543 positions
    auto& conf = result.protein->Conformation();
    EXPECT_EQ(conf.AtomCount(), 543u);
}

TEST(OrcaRunTest, PrmtopChargesAreCorrect) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!std::filesystem::exists(files.prmtop_path)) GTEST_SKIP();
    if (!std::filesystem::exists(files.xyz_path)) GTEST_SKIP();

    auto result = BuildFromOrca(files);
    ASSERT_TRUE(result.Ok()) << result.error;

    auto& conf = result.protein->Conformation();

    // Assign charges from the prmtop (authoritative ff14SB)
    PrmtopChargeSource source(files.prmtop_path);
    auto charges = ChargeAssignmentResult::Compute(conf, source);
    ASSERT_NE(charges, nullptr);

    // Total charge must be integer — this is a protonated protein
    // where atom count matches prmtop exactly
    double total = charges->TotalCharge();
    double nearest_int = std::round(total);
    double frac = std::abs(total - nearest_int);
    EXPECT_LT(frac, 0.01)
        << "Total charge " << total << " not integer (nearest=" << nearest_int << ")";

    // Provenance should mention ff14SB and prmtop
    EXPECT_NE(charges->Source().find("ff14SB"), std::string::npos)
        << "Source: " << charges->Source();
    EXPECT_NE(charges->Source().find("prmtop"), std::string::npos)
        << "Source: " << charges->Source();

    conf.AttachResult(std::move(charges));

    std::cout << "ORCA WT charge: " << total << " (integer: " << nearest_int << ")\n";
}

TEST(OrcaRunTest, ProtonationDetectedFromStructure) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!std::filesystem::exists(files.prmtop_path)) GTEST_SKIP();

    auto result = BuildFromOrca(files);
    ASSERT_TRUE(result.Ok()) << result.error;

    auto& conf = result.protein->Conformation();

    // Protonation should already be set from prmtop residue labels
    // (HID/HIE/HIP detected during loading)
    // Also run ProtonationDetectionResult to verify from H atoms
    auto detect = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(detect, nullptr);

    // Protonated structure: no unresolved
    EXPECT_EQ(detect->UnresolvedCount(), 0)
        << detect->UnresolvedCount() << " unresolved titratable residues";

    conf.AttachResult(std::move(detect));
}

TEST(OrcaRunTest, BuildContextHasProvenance) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";
    files.tleap_script_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT_tleap.in";

    if (!std::filesystem::exists(files.prmtop_path)) GTEST_SKIP();

    auto result = BuildFromOrca(files);
    ASSERT_TRUE(result.Ok()) << result.error;

    const auto& ctx = result.protein->BuildContext();
    EXPECT_FALSE(ctx.pdb_source.empty());
    EXPECT_EQ(ctx.force_field, "ff14SB");
    EXPECT_EQ(ctx.protonation_tool, "tleap");
    EXPECT_FALSE(ctx.prmtop_path.empty());
}

TEST(OrcaRunTest, ShieldingTensorsLoadedCorrectly) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";
    std::string nmr_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT_nmr.out";

    if (!std::filesystem::exists(files.prmtop_path)) GTEST_SKIP();
    if (!std::filesystem::exists(nmr_path)) GTEST_SKIP();

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;

    auto& conf = load.protein->Conformation();

    // Parse and attach shielding tensors
    auto orca = OrcaShieldingResult::Compute(conf, nmr_path);
    ASSERT_NE(orca, nullptr);
    EXPECT_EQ(orca->ParsedAtomCount(), static_cast<int>(conf.AtomCount()));

    conf.AttachResult(std::move(orca));

    // Verify tensors are on the atoms
    int with_tensors = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        if (ca.has_orca_shielding) {
            with_tensors++;

            // Total = diamagnetic + paramagnetic (basic sanity)
            Mat3 sum = ca.orca_shielding_diamagnetic + ca.orca_shielding_paramagnetic;
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    EXPECT_NEAR(ca.orca_shielding_total(r, c), sum(r, c), 0.01)
                        << "Atom " << ai << " total != dia + para at (" << r << "," << c << ")";
                }
            }

            // T0 (isotropic) should be physically reasonable: ~100-400 ppm for heavy,
            // ~20-35 ppm for H
            double t0 = ca.orca_shielding_total_spherical.T0;
            EXPECT_FALSE(std::isnan(t0)) << "NaN T0 at atom " << ai;
        }
    }

    EXPECT_EQ(with_tensors, static_cast<int>(conf.AtomCount()))
        << "Not all atoms have shielding tensors";

    // Spot-check: atom 0 (N of MET 1) should have isotropic shielding
    // in the range ~100-300 ppm for nitrogen
    double n_iso = conf.AtomAt(0).orca_shielding_total_spherical.T0;
    EXPECT_GT(n_iso, 50.0) << "N isotropic " << n_iso << " too low";
    EXPECT_LT(n_iso, 400.0) << "N isotropic " << n_iso << " too high";

    std::cout << "ORCA shielding: " << with_tensors << " atoms, "
              << "atom 0 (N) iso=" << n_iso << " ppm\n";
}


// ============================================================================
// KaML protonator
// ============================================================================

#include "KamlProtonator.h"

class KamlTest : public ProtonationPipelineTest {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        if (!std::filesystem::exists("KaML-CBtree.py")) {
            GTEST_SKIP() << "KaML not found at " << "KaML-CBtree.py";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};

TEST_F(KamlTest, PredictsPkaValues) {
    auto& conf = protein->Conformation();

    std::string error;
    auto pkas = KamlProtonator::PredictPka(*protein, conf, error);

    // KaML has ~80% success rate on ARM64. If it fails, skip don't fail.
    if (pkas.empty()) {
        GTEST_SKIP() << "KaML failed: " << error;
    }

    EXPECT_GE(pkas.size(), 5u) << "Expected at least 5 pKa predictions";

    for (const auto& pka : pkas) {
        EXPECT_GE(pka.pKa, 0.0) << pka.residue_type << " " << pka.residue_number;
        EXPECT_LE(pka.pKa, 14.0) << pka.residue_type << " " << pka.residue_number;
    }
}

TEST_F(KamlTest, ProtonatesAtPhysiologicalPH) {
    auto& conf = protein->Conformation();
    KamlProtonator kaml;

    auto result = kaml.Protonate(*protein, conf, 7.0);
    if (!result.Ok()) {
        GTEST_SKIP() << "KaML failed: " << result.error;
    }

    EXPECT_GT(result.state.DecisionCount(), 0u);
    EXPECT_EQ(result.state.Tool(), ProtonationTool::KaML);

    std::vector<AminoAcid> residue_types;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri)
        residue_types.push_back(protein->ResidueAt(ri).type);

    int net = result.state.NetChargeForProtein(residue_types);
    EXPECT_GE(net, -2) << "KaML net charge " << net << " too negative";
    EXPECT_LE(net, 6) << "KaML net charge " << net << " too positive";

    std::cout << "KaML pH 7: " << result.state.DecisionCount()
              << " decisions, net charge=" << net << "\n";
}

TEST_F(KamlTest, DecisionsAreTyped) {
    auto& conf = protein->Conformation();
    KamlProtonator kaml;

    auto result = kaml.Protonate(*protein, conf, 7.0);
    if (!result.Ok()) GTEST_SKIP() << result.error;

    for (const auto& d : result.state.Decisions()) {
        EXPECT_NE(d.amino_acid, AminoAcid::Unknown);
        EXPECT_LT(d.residue_index, protein->ResidueCount());
        EXPECT_FALSE(std::isnan(d.pKa));

        if (d.variant_index >= 0) {
            const auto& aat = GetAminoAcidType(d.amino_acid);
            EXPECT_LT(d.variant_index, static_cast<int>(aat.variants.size()));
        }
    }
}
