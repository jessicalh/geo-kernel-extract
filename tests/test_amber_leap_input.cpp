#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "AmberChargeResolver.h"
#include "AmberLeapInput.h"
#include "AmberPreparedChargeSource.h"
#include "AminoAcidType.h"
#include "Atom.h"
#include "ChargeSource.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinBuildContext.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RuntimeEnvironment.h"

#include <filesystem>
#include <sstream>

using namespace nmr;

// ============================================================================
// Synthetic-Protein helpers used to drive deterministic generator output
// without requiring an on-disk fixture.
// ============================================================================

namespace {

void AppendResidueAtomsFromAaType(Protein& protein,
                                   AminoAcid aa,
                                   size_t residue_index,
                                   std::vector<Vec3>& positions_out) {
    const AminoAcidType& aa_type = GetAminoAcidType(aa);
    for (const auto& templ : aa_type.atoms) {
        auto atom = Atom::Create(templ.element);
        atom->pdb_atom_name = templ.name;
        atom->residue_index = residue_index;
        size_t ai = protein.AddAtom(std::move(atom));
        protein.MutableResidueAt(residue_index).atom_indices.push_back(ai);
        positions_out.push_back(Vec3(0.0, 0.0, 0.0));
    }
}

std::unique_ptr<Protein> BuildOneResidueProtein(
        AminoAcid aa,
        int variant_index,
        ResidueTerminalState terminal_state = ResidueTerminalState::Internal,
        const std::string& chain_id = "A",
        int sequence_number = 1) {
    auto protein = std::make_unique<Protein>();

    Residue res;
    res.type = aa;
    res.sequence_number = sequence_number;
    res.chain_id = chain_id;
    res.protonation_variant_index = variant_index;
    res.protonation_state_resolved = true;
    res.terminal_state = terminal_state;
    size_t ri = protein->AddResidue(res);

    std::vector<Vec3> positions;
    AppendResidueAtomsFromAaType(*protein, aa, ri, positions);
    protein->AddConformation(std::move(positions), "leap-input-test");
    return protein;
}

// Build a minimal two-CYS protein with the SG atoms placed close enough
// to register as a disulfide via DetectDisulfides (typed SG-SG distance).
// No CovalentTopology / OpenBabel involvement — that is the whole point.
std::unique_ptr<Protein> BuildTwoCysWithSgDistance(double sg_distance) {
    auto protein = std::make_unique<Protein>();

    auto add_cys = [&](int seq, const std::string& chain_id) {
        Residue res;
        res.type = AminoAcid::CYS;
        res.sequence_number = seq;
        res.chain_id = chain_id;
        res.protonation_variant_index = -1;
        res.terminal_state = ResidueTerminalState::Internal;
        return protein->AddResidue(res);
    };

    size_t r0 = add_cys(1, "A");
    size_t r1 = add_cys(2, "A");

    // Place atoms with arbitrary coords; the only ones that matter for
    // DetectDisulfides are the two SG atoms.
    std::vector<Vec3> positions;
    auto add_atoms = [&](size_t residue_index,
                          const Vec3& offset) {
        const AminoAcidType& aa_type = GetAminoAcidType(AminoAcid::CYS);
        for (const auto& templ : aa_type.atoms) {
            auto atom = Atom::Create(templ.element);
            atom->pdb_atom_name = templ.name;
            atom->residue_index = residue_index;
            size_t ai = protein->AddAtom(std::move(atom));
            protein->MutableResidueAt(residue_index).atom_indices.push_back(ai);
            // Position SG atoms specifically; everything else far away
            // so it doesn't accidentally match anything.
            if (std::string(templ.name) == "SG") {
                positions.push_back(offset);
            } else {
                positions.push_back(offset + Vec3(100.0, 0.0, 0.0));
            }
        }
    };

    add_atoms(r0, Vec3(0.0, 0.0, 0.0));
    add_atoms(r1, Vec3(sg_distance, 0.0, 0.0));
    protein->AddConformation(std::move(positions), "two-cys-typed-sg-test");
    return protein;
}

// Convenience: count substring occurrences.
size_t CountSubstring(const std::string& haystack, const std::string& needle) {
    size_t count = 0;
    for (size_t pos = haystack.find(needle); pos != std::string::npos;
         pos = haystack.find(needle, pos + needle.size())) {
        ++count;
    }
    return count;
}

}  // namespace


// ============================================================================
// GenerateAmberPdb tests
// ============================================================================

class AmberLeapInputTest : public ::testing::Test {
protected:
    AmberSourceConfig MakeCfg() {
        AmberSourceConfig cfg;
        cfg.flat_table_path = nmr::test::TestEnvironment::Ff14sbParams();
        cfg.preparation_policy =
            AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;
        cfg.tleap_path =
            "/dev/null/placeholder-tleap-not-actually-invoked-test";
        return cfg;
    }

    AmberFlatTableCoverageVerdict MakeUnsupportedVerdict() {
        AmberFlatTableCoverageVerdict v;
        v.ok = false;
        AmberFlatTableCoverageFailure f;
        f.kind = AmberFlatTableCoverageKind::UnsupportedTerminalVariant;
        f.terminal_token = "NTERM";
        f.ff_residue_name = "ASH";
        f.residue_sequence_number = 1;
        v.failures.push_back(f);
        return v;
    }
};

TEST_F(AmberLeapInputTest, GeneratedPdbHisVariantUsesHidHieHipNames) {
    struct Case { int variant_index; const char* expected; };
    const Case cases[] = {
        {0, "HID"},
        {1, "HIE"},
        {2, "HIP"},
    };
    for (const auto& c : cases) {
        auto protein = BuildOneResidueProtein(AminoAcid::HIS, c.variant_index);
        AmberPreparedChargeSource src(
            *protein, MakeCfg().preparation_policy,
            MakeUnsupportedVerdict(), MakeCfg());
        const std::string pdb = src.GeneratedPdb(protein->Conformation());
        EXPECT_NE(pdb.find(c.expected), std::string::npos)
            << "variant " << c.variant_index << " expected " << c.expected
            << " in PDB:\n" << pdb;
    }
}

TEST_F(AmberLeapInputTest, GeneratedPdbAshFromAspVariantZero) {
    auto protein = BuildOneResidueProtein(
        AminoAcid::ASP, 0, ResidueTerminalState::Internal);
    AmberPreparedChargeSource src(
        *protein, MakeCfg().preparation_policy,
        MakeUnsupportedVerdict(), MakeCfg());
    const std::string pdb = src.GeneratedPdb(protein->Conformation());
    EXPECT_NE(pdb.find("ASH"), std::string::npos) << pdb;
    EXPECT_EQ(pdb.find("ASP"), std::string::npos) << pdb;
}

TEST_F(AmberLeapInputTest, DetectDisulfidesFindsPairAtBondingDistance) {
    auto protein = BuildTwoCysWithSgDistance(2.05);  // typical S-S bond
    auto pairs = amber_leap::DetectDisulfides(*protein, protein->Conformation());
    ASSERT_EQ(pairs.size(), 1u);
    EXPECT_EQ(pairs[0].first, 0u);
    EXPECT_EQ(pairs[0].second, 1u);
}

TEST_F(AmberLeapInputTest, DetectDisulfidesIgnoresPairAtVdwContact) {
    auto protein = BuildTwoCysWithSgDistance(3.6);  // non-bonded vdW contact
    auto pairs = amber_leap::DetectDisulfides(*protein, protein->Conformation());
    EXPECT_EQ(pairs.size(), 0u);
}

TEST_F(AmberLeapInputTest, GeneratedPdbDisulfidesEmitCYX) {
    auto protein = BuildTwoCysWithSgDistance(2.05);
    AmberPreparedChargeSource src(
        *protein, MakeCfg().preparation_policy,
        MakeUnsupportedVerdict(), MakeCfg());
    const std::string pdb = src.GeneratedPdb(protein->Conformation());

    // Both CYS residues should be emitted as CYX. Each CYS has multiple
    // atoms so CYX appears once per atom line of each residue.
    EXPECT_GT(CountSubstring(pdb, "CYX"), 2u) << pdb;
    EXPECT_EQ(CountSubstring(pdb, " CYS "), 0u) << pdb;
}

TEST_F(AmberLeapInputTest, GeneratedPdbAtomSerialIsSequential) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(r.Ok()) << r.error;

    AmberPreparedChargeSource src(
        *r.protein, MakeCfg().preparation_policy,
        MakeUnsupportedVerdict(), MakeCfg());
    const std::string pdb = src.GeneratedPdb(r.protein->Conformation());

    // First ATOM serial must be 1, last must equal AtomCount.
    auto first = pdb.find("ATOM  ");
    ASSERT_NE(first, std::string::npos);
    int first_serial = std::stoi(pdb.substr(first + 6, 5));
    EXPECT_EQ(first_serial, 1);

    auto last = pdb.rfind("ATOM  ");
    ASSERT_NE(last, std::string::npos);
    int last_serial = std::stoi(pdb.substr(last + 6, 5));
    EXPECT_EQ(static_cast<size_t>(last_serial), r.protein->AtomCount());
}

TEST_F(AmberLeapInputTest, GeneratedPdbTerRecordsBetweenChains) {
    // Two single-residue chains. Expect one TER between them and a
    // closing TER + END at the end.
    auto protein = std::make_unique<Protein>();
    auto add = [&](int seq, const std::string& chain_id) {
        Residue res;
        res.type = AminoAcid::ALA;
        res.sequence_number = seq;
        res.chain_id = chain_id;
        res.terminal_state = ResidueTerminalState::Internal;
        return protein->AddResidue(res);
    };
    size_t r0 = add(1, "A");
    size_t r1 = add(1, "B");
    std::vector<Vec3> positions;
    AppendResidueAtomsFromAaType(*protein, AminoAcid::ALA, r0, positions);
    AppendResidueAtomsFromAaType(*protein, AminoAcid::ALA, r1, positions);
    protein->AddConformation(std::move(positions), "two-chain-test");

    AmberPreparedChargeSource src(
        *protein, MakeCfg().preparation_policy,
        MakeUnsupportedVerdict(), MakeCfg());
    const std::string pdb = src.GeneratedPdb(protein->Conformation());

    EXPECT_GE(CountSubstring(pdb, "TER\n"), 2u) << pdb;
    EXPECT_NE(pdb.find("END\n"), std::string::npos) << pdb;
}

// ----------------------------------------------------------------------------
// Step-6 capping policy tests
// ----------------------------------------------------------------------------

namespace {

// Build a single-residue protein with full atom set (so backbone N and C
// indices are populated by FinalizeConstruction-equivalent logic) and
// the requested terminal_state/variant.
std::unique_ptr<Protein> BuildSingleResidueProteinWithBackbone(
        AminoAcid aa,
        int variant_index,
        ResidueTerminalState terminal_state,
        int sequence_number = 1,
        const std::string& chain_id = "A") {
    auto protein = std::make_unique<Protein>();

    Residue res;
    res.type = aa;
    res.sequence_number = sequence_number;
    res.chain_id = chain_id;
    res.protonation_variant_index = variant_index;
    res.protonation_state_resolved = true;
    res.terminal_state = terminal_state;
    size_t ri = protein->AddResidue(res);

    std::vector<Vec3> positions;
    const AminoAcidType& aa_type = GetAminoAcidType(aa);
    for (const auto& templ : aa_type.atoms) {
        auto atom = Atom::Create(templ.element);
        atom->pdb_atom_name = templ.name;
        atom->residue_index = ri;
        size_t ai = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(ai);
        positions.push_back(Vec3(0.0, 0.0, 0.0));

        // Cache backbone slot indices the same way Protein::CacheResidueBackboneIndices
        // does — the cap geometry generator reads res.N and res.C.
        std::string aname(templ.name);
        if (aname == "N")  protein->MutableResidueAt(ri).N  = ai;
        if (aname == "CA") protein->MutableResidueAt(ri).CA = ai;
        if (aname == "C")  protein->MutableResidueAt(ri).C  = ai;
        if (aname == "O")  protein->MutableResidueAt(ri).O  = ai;
    }
    protein->AddConformation(std::move(positions), "single-res-cap-test");
    return protein;
}

AmberFlatTableCoverageVerdict MakeNTermVerdict(
        const std::string& ff_resname,
        int seq_num = 1,
        const std::string& chain_id = "A") {
    AmberFlatTableCoverageVerdict v;
    v.ok = false;
    AmberFlatTableCoverageFailure f;
    f.kind = AmberFlatTableCoverageKind::UnsupportedTerminalVariant;
    f.terminal_token = "NTERM";
    f.ff_residue_name = ff_resname;
    f.residue_sequence_number = seq_num;
    f.chain_id = chain_id;
    v.failures.push_back(f);
    return v;
}

AmberFlatTableCoverageVerdict MakeCTermVerdict(
        const std::string& ff_resname,
        int seq_num = 1,
        const std::string& chain_id = "A") {
    AmberFlatTableCoverageVerdict v;
    v.ok = false;
    AmberFlatTableCoverageFailure f;
    f.kind = AmberFlatTableCoverageKind::UnsupportedTerminalVariant;
    f.terminal_token = "CTERM";
    f.ff_residue_name = ff_resname;
    f.residue_sequence_number = seq_num;
    f.chain_id = chain_id;
    v.failures.push_back(f);
    return v;
}

AmberFlatTableCoverageVerdict MakeNCTermVerdict(
        const std::string& ff_resname,
        int seq_num = 1,
        const std::string& chain_id = "A") {
    AmberFlatTableCoverageVerdict v;
    v.ok = false;
    AmberFlatTableCoverageFailure f;
    f.kind = AmberFlatTableCoverageKind::UnsupportedTerminalVariant;
    f.terminal_token = "NCTERM";
    f.ff_residue_name = ff_resname;
    f.residue_sequence_number = seq_num;
    f.chain_id = chain_id;
    v.failures.push_back(f);
    return v;
}

}  // namespace

TEST_F(AmberLeapInputTest, GeneratedPdbNTermAshUnderCappingHasAce) {
    auto protein = BuildSingleResidueProteinWithBackbone(
        AminoAcid::ASP, 0, ResidueTerminalState::NTerminus);
    AmberSourceConfig cfg = MakeCfg();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  MakeNTermVerdict("ASH"), cfg);
    const std::string pdb = src.GeneratedPdb(protein->Conformation());

    EXPECT_NE(pdb.find(" ACE "), std::string::npos) << pdb;
    EXPECT_NE(pdb.find(" ASH "), std::string::npos) << pdb;
    EXPECT_GT(CountSubstring(pdb, "HH31"), 0u);
    // ACE precedes ASH in the file.
    EXPECT_LT(pdb.find(" ACE "), pdb.find(" ASH ")) << pdb;
}

TEST_F(AmberLeapInputTest, GeneratedPdbCTermGlhUnderCappingHasNme) {
    auto protein = BuildSingleResidueProteinWithBackbone(
        AminoAcid::GLU, 0, ResidueTerminalState::CTerminus);
    AmberSourceConfig cfg = MakeCfg();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  MakeCTermVerdict("GLH"), cfg);
    const std::string pdb = src.GeneratedPdb(protein->Conformation());

    EXPECT_NE(pdb.find(" GLH "), std::string::npos) << pdb;
    EXPECT_NE(pdb.find(" NME "), std::string::npos) << pdb;
    // GLH precedes NME in the file.
    EXPECT_LT(pdb.find(" GLH "), pdb.find(" NME ")) << pdb;
}

TEST_F(AmberLeapInputTest, GeneratedPdbNCTermLynUnderCappingHasBoth) {
    auto protein = BuildSingleResidueProteinWithBackbone(
        AminoAcid::LYS, 0, ResidueTerminalState::NAndCTerminus);
    AmberSourceConfig cfg = MakeCfg();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  MakeNCTermVerdict("LYN"), cfg);
    const std::string pdb = src.GeneratedPdb(protein->Conformation());

    EXPECT_NE(pdb.find(" ACE "), std::string::npos) << pdb;
    EXPECT_NE(pdb.find(" LYN "), std::string::npos) << pdb;
    EXPECT_NE(pdb.find(" NME "), std::string::npos) << pdb;
    // Order: ACE < LYN < NME.
    EXPECT_LT(pdb.find(" ACE "), pdb.find(" LYN ")) << pdb;
    EXPECT_LT(pdb.find(" LYN "), pdb.find(" NME ")) << pdb;
}

TEST_F(AmberLeapInputTest, GeneratedPdbFailPolicyDoesNotCap) {
    auto protein = BuildSingleResidueProteinWithBackbone(
        AminoAcid::ASP, 0, ResidueTerminalState::NTerminus);
    AmberSourceConfig cfg = MakeCfg();
    cfg.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  MakeNTermVerdict("ASH"), cfg);
    const std::string pdb = src.GeneratedPdb(protein->Conformation());
    // No caps under Fail policy. The PDB still emits "ASH" but no ACE.
    EXPECT_EQ(pdb.find(" ACE "), std::string::npos) << pdb;
    EXPECT_EQ(pdb.find(" NME "), std::string::npos) << pdb;
    EXPECT_NE(pdb.find(" ASH "), std::string::npos) << pdb;
}

TEST_F(AmberLeapInputTest, GeneratedPdbHidNTermDoesNotCapEvenUnderCappingPolicy) {
    // HID is supported as a NTERM template in stock ff14SB, so the
    // verdict should not include it as an unsupported variant. But even
    // if some upstream code accidentally fed a "HID" UnsupportedTerminalVariant
    // failure into the source, capping policy must NOT cap HID — only
    // the cappable set (ASH, CYM, GLH, LYN) gets caps.
    auto protein = BuildSingleResidueProteinWithBackbone(
        AminoAcid::HIS, 0 /* HID */, ResidueTerminalState::NTerminus);
    AmberSourceConfig cfg = MakeCfg();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  MakeNTermVerdict("HID"), cfg);
    const std::string pdb = src.GeneratedPdb(protein->Conformation());
    EXPECT_EQ(pdb.find(" ACE "), std::string::npos) << pdb;
}

TEST_F(AmberLeapInputTest, ResidueMappingMarksCapsAsNoneForCap) {
    auto protein = BuildSingleResidueProteinWithBackbone(
        AminoAcid::ASP, 0, ResidueTerminalState::NTerminus);
    AmberSourceConfig cfg = MakeCfg();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  MakeNTermVerdict("ASH"), cfg);
    (void)src.GeneratedPdb(protein->Conformation());

    const auto& mapping = src.ResidueMapping();
    ASSERT_EQ(mapping.extractor_index_for_prmtop_residue.size(), 2u);
    EXPECT_EQ(mapping.extractor_index_for_prmtop_residue[0],
              amber_leap::ResidueAmberMapping::NONE_FOR_CAP);
    EXPECT_EQ(mapping.extractor_index_for_prmtop_residue[1], 0u);
}

TEST_F(AmberLeapInputTest, DescribeRecordsCappingDecisions) {
    auto protein = BuildSingleResidueProteinWithBackbone(
        AminoAcid::GLU, 0, ResidueTerminalState::CTerminus);
    AmberSourceConfig cfg = MakeCfg();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  MakeCTermVerdict("GLH"), cfg);
    const std::string desc = src.Describe();
    EXPECT_NE(desc.find("UseCappedFragmentsForUnsupportedTerminalVariants"),
              std::string::npos) << desc;
    EXPECT_NE(desc.find("GLH"), std::string::npos) << desc;
    EXPECT_NE(desc.find("NME"), std::string::npos) << desc;
}

TEST_F(AmberLeapInputTest, ResidueAmberMappingNoCapsIsIdentity) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(r.Ok()) << r.error;

    AmberPreparedChargeSource src(
        *r.protein, MakeCfg().preparation_policy,
        MakeUnsupportedVerdict(), MakeCfg());
    (void)src.GeneratedPdb(r.protein->Conformation());
    const auto& mapping = src.ResidueMapping();

    ASSERT_EQ(mapping.extractor_index_for_prmtop_residue.size(),
              r.protein->ResidueCount());
    for (size_t i = 0; i < mapping.extractor_index_for_prmtop_residue.size(); ++i) {
        EXPECT_EQ(mapping.extractor_index_for_prmtop_residue[i], i)
            << "step-4 mapping must be identity (no caps)";
    }
}


// ============================================================================
// GenerateLeapScript tests
// ============================================================================

TEST_F(AmberLeapInputTest, GeneratedLeapScriptNoDisulfidesHasFiveLines) {
    auto protein = BuildOneResidueProtein(AminoAcid::ALA, -1);
    AmberPreparedChargeSource src(
        *protein, MakeCfg().preparation_policy,
        MakeUnsupportedVerdict(), MakeCfg());
    (void)src.GeneratedPdb(protein->Conformation());  // populate mapping
    const std::string script = src.GeneratedLeapScript(
        "/tmp/in.pdb", "/tmp/out.prmtop", "/tmp/out.inpcrd");

    EXPECT_NE(script.find("source leaprc.protein.ff14SB\n"),
              std::string::npos);
    EXPECT_NE(script.find("set default PBRadii mbondi2\n"),
              std::string::npos);
    EXPECT_NE(script.find("mol = loadPdb /tmp/in.pdb\n"),
              std::string::npos);
    EXPECT_NE(script.find("saveamberparm mol /tmp/out.prmtop /tmp/out.inpcrd\n"),
              std::string::npos);
    EXPECT_NE(script.find("quit\n"), std::string::npos);
    EXPECT_EQ(CountSubstring(script, "\nbond mol."), 0u);
}

TEST_F(AmberLeapInputTest, GeneratedLeapScriptWithDisulfideEmitsBondLine) {
    auto protein = BuildTwoCysWithSgDistance(2.05);
    AmberPreparedChargeSource src(
        *protein, MakeCfg().preparation_policy,
        MakeUnsupportedVerdict(), MakeCfg());
    (void)src.GeneratedPdb(protein->Conformation());
    const std::string script = src.GeneratedLeapScript(
        "/tmp/in.pdb", "/tmp/out.prmtop", "/tmp/out.inpcrd");

    // Expect "bond mol.1.SG mol.2.SG" — 1-based PRMTOP residue indices.
    EXPECT_NE(script.find("bond mol.1.SG mol.2.SG\n"),
              std::string::npos) << script;
}


// ============================================================================
// Hard-precondition tests (re-tleap protection)
// ============================================================================

class AmberPreparedChargePreconditionDeathTest : public ::testing::Test {};

TEST_F(AmberPreparedChargePreconditionDeathTest,
       AbortsWhenChargesAlreadyPresent) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(r.Ok()) << r.error;
    // r.protein already has ForceFieldCharges (BuildFromProtonatedPdb
    // ran PrepareForceFieldCharges via the resolver).

    AmberSourceConfig cfg;
    cfg.flat_table_path = nmr::test::TestEnvironment::Ff14sbParams();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;
    cfg.tleap_path = "/dev/null/placeholder-tleap";

    AmberFlatTableCoverageVerdict verdict;
    verdict.ok = false;
    AmberFlatTableCoverageFailure f;
    f.kind = AmberFlatTableCoverageKind::UnsupportedTerminalVariant;
    f.terminal_token = "NTERM";
    f.ff_residue_name = "ASH";
    verdict.failures.push_back(f);

    AmberPreparedChargeSource src(*r.protein, cfg.preparation_policy,
                                  verdict, cfg);
    std::string err;
    EXPECT_DEATH(
        src.LoadCharges(*r.protein, r.protein->Conformation(), err),
        "re-tleap protection");
}

TEST_F(AmberPreparedChargePreconditionDeathTest,
       AbortsWhenBuildContextHasPrmtopPath) {
    auto protein = std::make_unique<Protein>();
    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.terminal_state = ResidueTerminalState::Internal;
    protein->AddResidue(res);

    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->prmtop_path = "/tmp/some_upstream.prmtop";
    protein->SetBuildContext(std::move(ctx));

    std::vector<Vec3> positions;
    AppendResidueAtomsFromAaType(*protein, AminoAcid::ALA, 0, positions);
    protein->AddConformation(std::move(positions), "test");

    AmberSourceConfig cfg;
    cfg.flat_table_path = nmr::test::TestEnvironment::Ff14sbParams();
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;
    cfg.tleap_path = "/dev/null/placeholder-tleap";
    AmberFlatTableCoverageVerdict verdict;
    verdict.ok = false;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  verdict, cfg);
    std::string err;
    EXPECT_DEATH(
        src.LoadCharges(*protein, protein->Conformation(), err),
        "re-tleap protection");
}

TEST(AmberPreparedChargeStep5NegativeTest,
     LoadChargesFailsWhenTleapPathDoesNotExist) {
    // When config.tleap_path is set to a non-existent path AND
    // RuntimeEnvironment::Tleap() is also unset/non-existent, LoadCharges
    // should fail with a named error mentioning "no usable tleap binary".
    auto protein = std::make_unique<Protein>();
    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.terminal_state = ResidueTerminalState::Internal;
    protein->AddResidue(res);
    std::vector<Vec3> positions;
    AppendResidueAtomsFromAaType(*protein, AminoAcid::ALA, 0, positions);
    protein->AddConformation(std::move(positions), "test");

    AmberSourceConfig cfg;
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;
    cfg.tleap_path = "/nonexistent/path/to/tleap_binary_for_test";
    AmberFlatTableCoverageVerdict verdict;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  verdict, cfg);
    std::string err;
    auto rows = src.LoadCharges(*protein, protein->Conformation(), err);

    // Outcome depends on whether RuntimeEnvironment::Tleap() is set on
    // this machine. Either result is acceptable for this test, but in
    // both cases LoadCharges should NOT produce real charges from a
    // non-existent config path; the rows must be empty.
    if (RuntimeEnvironment::Tleap().empty()) {
        EXPECT_TRUE(rows.empty());
        EXPECT_NE(err.find("no usable tleap binary"), std::string::npos) << err;
    } else {
        // RuntimeEnvironment fallback found a real tleap; the run will
        // attempt the full pipeline. We still require empty result OR a
        // legitimate parse outcome — but we don't assert success here
        // because the synthetic ALA protein may or may not survive a
        // real tleap roundtrip. The point of this test is the
        // "no usable tleap binary" path; if RuntimeEnvironment has one,
        // we skip that assertion.
        GTEST_SKIP() << "RuntimeEnvironment::Tleap() is set on this "
                        "machine; this test exercises only the no-tleap "
                        "branch.";
    }
}
