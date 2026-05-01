#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "Protein.h"
#include "ProteinBuildContext.h"
#include "ChargeSource.h"
#include "AmberChargeResolver.h"
#include "AminoAcidType.h"
#include "Atom.h"
#include "Residue.h"
#include <filesystem>
#include <fstream>

using namespace nmr;

// ============================================================================
// Synthetic Protein builders, parallel to BuildSingleAtomVariantProtein in
// test_foundation_results.cpp but exposing the full atom set so the verdict
// can walk every (terminal, resname, atom_name) triple.
// ============================================================================

namespace {

// Build a one-residue protein with all the atoms named in the AminoAcidType
// canonical_atoms list. terminal_state defaults to Internal.
std::unique_ptr<Protein> BuildOneResidueProtein(
        AminoAcid aa,
        int variant_index,
        ResidueTerminalState terminal_state = ResidueTerminalState::Internal) {
    auto protein = std::make_unique<Protein>();

    Residue res;
    res.type = aa;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.protonation_variant_index = variant_index;
    res.protonation_state_resolved = true;
    res.terminal_state = terminal_state;
    size_t ri = protein->AddResidue(res);

    const AminoAcidType& aa_type = GetAminoAcidType(aa);
    for (const auto& templ : aa_type.atoms) {
        auto atom = Atom::Create(templ.element);
        atom->pdb_atom_name = templ.name;
        atom->residue_index = ri;
        size_t ai = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(ai);
    }

    std::vector<Vec3> positions(protein->AtomCount(), Vec3::Zero());
    protein->AddConformation(std::move(positions), "verdict-test");
    return protein;
}

// A custom dat file that includes a terminal template for a residue but
// drops one specific atom. Caller cleans up.
std::filesystem::path WriteCustomDatWithMissingAtom(
        const std::string& terminal_token,
        const std::string& ff_resname,
        const std::string& dropped_atom_name) {
    auto path = std::filesystem::temp_directory_path() /
        "nmr_amber_resolver_missing_atom.dat";
    std::ofstream out(path);
    EXPECT_TRUE(out.is_open());

    // Write minimal rows: every canonical atom for the residue except one.
    AminoAcid aa = AminoAcidFromThreeLetterCode(ff_resname);
    const AminoAcidType& aa_type = GetAminoAcidType(
        aa == AminoAcid::Unknown ? AminoAcid::TYR : aa);
    for (const auto& templ : aa_type.atoms) {
        if (std::string(templ.name) == dropped_atom_name) continue;
        out << terminal_token << " " << ff_resname << " " << templ.name
            << " 0.0 1.5\n";
    }
    return path;
}

}  // namespace


// ============================================================================
// Verdict structure invariants
// ============================================================================

TEST(AmberFlatTableCoverageVerdict, OkVerdictHasNoFailuresAndDescribesItself) {
    AmberFlatTableCoverageVerdict v;
    EXPECT_TRUE(v.Ok());
    EXPECT_TRUE(v.failures.empty());
    EXPECT_EQ(v.Detail(), "Satisfiable");
}

TEST(AmberFlatTableCoverageVerdict, NonOkDetailContainsCanonicalFallbackPhrase) {
    AmberFlatTableCoverageVerdict v;
    v.ok = false;
    AmberFlatTableCoverageFailure f;
    f.kind = AmberFlatTableCoverageKind::UnsupportedTerminalVariant;
    f.terminal_token = "NTERM";
    f.ff_residue_name = "ASH";
    f.residue_sequence_number = 1;
    v.failures.push_back(f);

    const std::string detail = v.Detail();
    EXPECT_NE(detail.find("NTERM"), std::string::npos) << detail;
    EXPECT_NE(detail.find("ASH"), std::string::npos) << detail;
    EXPECT_NE(detail.find("no canonical fallback"), std::string::npos) << detail;
}

TEST(AmberPreparationPolicyName, AllValuesNamed) {
    EXPECT_STREQ(
        "UseStockTermini",
        AmberPreparationPolicyName(AmberPreparationPolicy::UseStockTermini));
    EXPECT_STREQ(
        "UseCappedFragmentsForUnsupportedTerminalVariants",
        AmberPreparationPolicyName(
            AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants));
    EXPECT_STREQ(
        "FailOnUnsupportedTerminalVariants",
        AmberPreparationPolicyName(
            AmberPreparationPolicy::FailOnUnsupportedTerminalVariants));
}


// ============================================================================
// AnalyzeFlatTableCoverage on real and synthetic proteins
// ============================================================================

class AmberChargeResolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) {
            GTEST_SKIP() << "ff14sb_params.dat not found";
        }
    }
};

TEST_F(AmberChargeResolverTest, SatisfiableOnUbq) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(r.Ok()) << r.error;

    auto verdict = AnalyzeFlatTableCoverage(
        *r.protein, nmr::test::TestEnvironment::Ff14sbParams());
    EXPECT_TRUE(verdict.Ok()) << verdict.Detail();
    EXPECT_TRUE(verdict.failures.empty());
}

TEST_F(AmberChargeResolverTest, UnsupportedTerminalVariantOnNTermAsh) {
    // ASP variant 0 = ASH (per the AminoAcidType variants contract).
    // ff14SB has INTERNAL ASH but no NTERM ASH.
    auto protein = BuildOneResidueProtein(
        AminoAcid::ASP, 0, ResidueTerminalState::NTerminus);

    auto verdict = AnalyzeFlatTableCoverage(
        *protein, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_FALSE(verdict.Ok()) << "expected NTERM ASH to be unsupported";
    ASSERT_FALSE(verdict.failures.empty());
    EXPECT_EQ(verdict.failures.front().kind,
              AmberFlatTableCoverageKind::UnsupportedTerminalVariant);
    EXPECT_EQ(verdict.failures.front().terminal_token, "NTERM");
    EXPECT_EQ(verdict.failures.front().ff_residue_name, "ASH");
    EXPECT_EQ(verdict.failures.front().residue_sequence_number, 1);
    EXPECT_EQ(verdict.failures.front().chain_id, "A");
}

TEST_F(AmberChargeResolverTest, UnsupportedTerminalVariantOnCTermLyn) {
    // LYS variant 0 = LYN. ff14SB has INTERNAL LYN but no CTERM LYN.
    auto protein = BuildOneResidueProtein(
        AminoAcid::LYS, 0, ResidueTerminalState::CTerminus);

    auto verdict = AnalyzeFlatTableCoverage(
        *protein, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_FALSE(verdict.Ok());
    EXPECT_EQ(verdict.failures.front().kind,
              AmberFlatTableCoverageKind::UnsupportedTerminalVariant);
    EXPECT_EQ(verdict.failures.front().terminal_token, "CTERM");
    EXPECT_EQ(verdict.failures.front().ff_residue_name, "LYN");
}

TEST_F(AmberChargeResolverTest, UnsupportedResidueOnTymAtAnyTerminal) {
    // TYR variant 0 = TYM. ff14SB has no TYM at any terminal state in our
    // regenerated table — so this is UnsupportedResidue, not just an
    // unsupported terminal.
    auto protein = BuildOneResidueProtein(
        AminoAcid::TYR, 0, ResidueTerminalState::Internal);

    auto verdict = AnalyzeFlatTableCoverage(
        *protein, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_FALSE(verdict.Ok());
    EXPECT_EQ(verdict.failures.front().kind,
              AmberFlatTableCoverageKind::UnsupportedResidue);
    EXPECT_EQ(verdict.failures.front().ff_residue_name, "TYM");
}

TEST_F(AmberChargeResolverTest, MissingAtomNameWhenTemplateExistsButLacksAtom) {
    // Build a custom dat that has every TYR atom EXCEPT OH. A protein
    // containing OH must verdict-fail with MissingAtomName.
    auto custom_dat = WriteCustomDatWithMissingAtom("INTERNAL", "TYR", "OH");
    ASSERT_TRUE(std::filesystem::exists(custom_dat));

    auto protein = BuildOneResidueProtein(
        AminoAcid::TYR, -1, ResidueTerminalState::Internal);

    auto verdict = AnalyzeFlatTableCoverage(*protein, custom_dat.string());
    ASSERT_FALSE(verdict.Ok()) << verdict.Detail();
    EXPECT_EQ(verdict.failures.front().kind,
              AmberFlatTableCoverageKind::MissingAtomName);
    EXPECT_EQ(verdict.failures.front().ff_residue_name, "TYR");
    EXPECT_EQ(verdict.failures.front().atom_name, "OH");

    std::error_code ec;
    std::filesystem::remove(custom_dat, ec);
}

TEST_F(AmberChargeResolverTest, FirstFailureMatchesEmittedDetailSubstrings) {
    // The Detail() string must always carry the first-failure terminal
    // token, the ff_residue_name, and the canonical-fallback phrase.
    // Existing tests (test_foundation_results.cpp) rely on these
    // substrings appearing in error_out.
    auto protein = BuildOneResidueProtein(
        AminoAcid::ASP, 0, ResidueTerminalState::NTerminus);

    auto verdict = AnalyzeFlatTableCoverage(
        *protein, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_FALSE(verdict.Ok());

    const std::string detail = verdict.Detail();
    EXPECT_NE(detail.find(verdict.failures.front().terminal_token),
              std::string::npos) << detail;
    EXPECT_NE(detail.find(verdict.failures.front().ff_residue_name),
              std::string::npos) << detail;
    EXPECT_NE(detail.find("no canonical fallback"), std::string::npos)
        << detail;
}

TEST_F(AmberChargeResolverTest, ResolverPicksPrmtopWhenBuildContextHasIt) {
    // BuildFromOrca always sets prmtop_path; the resolver's branch 1
    // returns a PrmtopChargeSource regardless of flat-table coverage.
    const std::string orca_dir = nmr::test::TestEnvironment::OrcaDir();
    if (orca_dir.empty()) {
        GTEST_SKIP() << "ORCA test data not available";
    }
    OrcaRunFiles files;
    files.pdb_path     = orca_dir + "A0A7C5FAR6_WT.pdb";
    files.xyz_path     = orca_dir + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path  = orca_dir + "A0A7C5FAR6_WT.prmtop";
    if (!std::filesystem::exists(files.prmtop_path)) {
        GTEST_SKIP() << "ORCA test prmtop not found at " << files.prmtop_path;
    }

    auto r = BuildFromOrca(files);
    ASSERT_TRUE(r.Ok()) << r.error;
    ASSERT_NE(r.charges, nullptr);
    EXPECT_EQ(r.charges->Kind(), ChargeModelKind::AmberPrmtop);
    EXPECT_EQ(r.charges->SourceForceField(), ForceField::Amber_ff14SB);
}

TEST_F(AmberChargeResolverTest, ResolverPicksFlatTableForStockProtein) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(r.Ok()) << r.error;
    ASSERT_NE(r.charges, nullptr);
    EXPECT_EQ(r.charges->Kind(), ChargeModelKind::Ff14SBParamFile);
}

TEST_F(AmberChargeResolverTest, ResolverFailsForMissingPrmtopFile) {
    auto protein = std::make_unique<Protein>();
    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.terminal_state = ResidueTerminalState::Internal;
    protein->AddResidue(res);

    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->prmtop_path = "/tmp/this_file_does_not_exist_for_resolver_test.prmtop";
    protein->SetBuildContext(std::move(ctx));

    AmberSourceConfig cfg;
    cfg.flat_table_path = nmr::test::TestEnvironment::Ff14sbParams();
    cfg.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;

    std::string err;
    auto src = ResolveAmberChargeSource(
        *protein, protein->BuildContext(), cfg, err);
    EXPECT_EQ(src, nullptr);
    EXPECT_NE(err.find("does not exist"), std::string::npos) << err;
}

TEST_F(AmberChargeResolverTest, ResolverFailsForUnsupportedTerminalVariantUnderFailPolicy) {
    auto protein = BuildOneResidueProtein(
        AminoAcid::ASP, 0, ResidueTerminalState::NTerminus);
    auto ctx = std::make_unique<ProteinBuildContext>();
    protein->SetBuildContext(std::move(ctx));

    AmberSourceConfig cfg;
    cfg.flat_table_path = nmr::test::TestEnvironment::Ff14sbParams();
    cfg.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;

    std::string err;
    auto src = ResolveAmberChargeSource(
        *protein, protein->BuildContext(), cfg, err);
    EXPECT_EQ(src, nullptr);
    EXPECT_NE(err.find("FailOnUnsupportedTerminalVariants"),
              std::string::npos) << err;
    EXPECT_NE(err.find("ASH"), std::string::npos) << err;
    EXPECT_NE(err.find("NTERM"), std::string::npos) << err;
}

TEST_F(AmberChargeResolverTest, ResolverProducesAmberPreparedSourceUnderCappingPolicy) {
    // Step-4 behaviour: the resolver constructs a real
    // AmberPreparedChargeSource under
    // UseCappedFragmentsForUnsupportedTerminalVariants. The source
    // identifies itself via Kind() == AmberPreparedPrmtop and carries
    // the verdict reason for downstream provenance.
    auto protein = BuildOneResidueProtein(
        AminoAcid::ASP, 0, ResidueTerminalState::NTerminus);
    auto ctx = std::make_unique<ProteinBuildContext>();
    protein->SetBuildContext(std::move(ctx));

    AmberSourceConfig cfg;
    cfg.flat_table_path = nmr::test::TestEnvironment::Ff14sbParams();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;
    cfg.tleap_path = "/dev/null/placeholder-tleap-not-actually-invoked-step-4";

    std::string err;
    auto src = ResolveAmberChargeSource(
        *protein, protein->BuildContext(), cfg, err);
    ASSERT_NE(src, nullptr) << err;
    EXPECT_EQ(src->Kind(), ChargeModelKind::AmberPreparedPrmtop);
    EXPECT_EQ(src->SourceForceField(), ForceField::Amber_ff14SB);

    const std::string desc = src->Describe();
    EXPECT_NE(desc.find("AmberPreparedPrmtop"), std::string::npos) << desc;
    EXPECT_NE(desc.find("UseCappedFragmentsForUnsupportedTerminalVariants"),
              std::string::npos) << desc;
    EXPECT_NE(desc.find("ASH"), std::string::npos) << desc;
}

TEST_F(AmberChargeResolverTest, ResolverFailsWhenFlatTablePathEmpty) {
    auto protein = BuildOneResidueProtein(
        AminoAcid::ALA, -1, ResidueTerminalState::Internal);
    auto ctx = std::make_unique<ProteinBuildContext>();
    protein->SetBuildContext(std::move(ctx));

    AmberSourceConfig cfg;  // flat_table_path empty
    cfg.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;

    std::string err;
    auto src = ResolveAmberChargeSource(
        *protein, protein->BuildContext(), cfg, err);
    EXPECT_EQ(src, nullptr);
    EXPECT_NE(err.find("flat_table_path"), std::string::npos) << err;
}

TEST_F(AmberChargeResolverTest, MultipleFailuresCollectedInOneWalk) {
    // Two-residue chain: NTERM ASH at residue 1, CTERM LYN at residue 2.
    // Both are unsupported terminal variants. The verdict must collect
    // both failures (capping policy needs each end's diagnostic
    // separately in step 6).
    auto protein = std::make_unique<Protein>();

    auto add_residue = [&](AminoAcid aa, int variant_index,
                           ResidueTerminalState ts, int seq) {
        Residue res;
        res.type = aa;
        res.sequence_number = seq;
        res.chain_id = "A";
        res.protonation_variant_index = variant_index;
        res.protonation_state_resolved = true;
        res.terminal_state = ts;
        size_t ri = protein->AddResidue(res);

        const AminoAcidType& aa_type = GetAminoAcidType(aa);
        for (const auto& templ : aa_type.atoms) {
            auto atom = Atom::Create(templ.element);
            atom->pdb_atom_name = templ.name;
            atom->residue_index = ri;
            size_t ai = protein->AddAtom(std::move(atom));
            protein->MutableResidueAt(ri).atom_indices.push_back(ai);
        }
    };

    add_residue(AminoAcid::ASP, 0, ResidueTerminalState::NTerminus, 1);
    add_residue(AminoAcid::LYS, 0, ResidueTerminalState::CTerminus, 2);
    std::vector<Vec3> positions(protein->AtomCount(), Vec3::Zero());
    protein->AddConformation(std::move(positions), "two-residue-test");

    auto verdict = AnalyzeFlatTableCoverage(
        *protein, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_FALSE(verdict.Ok());
    ASSERT_GE(verdict.failures.size(), 2u);

    bool saw_nterm_ash = false;
    bool saw_cterm_lyn = false;
    for (const auto& f : verdict.failures) {
        if (f.terminal_token == "NTERM" && f.ff_residue_name == "ASH")
            saw_nterm_ash = true;
        if (f.terminal_token == "CTERM" && f.ff_residue_name == "LYN")
            saw_cterm_lyn = true;
    }
    EXPECT_TRUE(saw_nterm_ash);
    EXPECT_TRUE(saw_cterm_lyn);
}
