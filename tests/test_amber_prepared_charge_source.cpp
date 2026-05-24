// Step-5 integration tests: full AmberPreparedChargeSource roundtrip
// through tleap. Gated by RuntimeEnvironment::Tleap() availability —
// SKIPped when no tleap binary is found.
//
// Step-5 scope: verify the pipeline (PDB gen → tleap → PRMTOP parse →
// atom mapping) wires correctly and fails loudly on bad inputs. The
// successful round-trip on a stock protein (1UBQ-like) requires a load
// path that produces a Protein WITHOUT pre-loaded charges; the existing
// BuildFromProtonatedPdb pre-charges, which trips the re-tleap death
// guard. That round-trip is exercised in step 6's capping integration
// tests, where realistic minimal proteins drive the prepared path
// naturally.
//

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

#include <cmath>
#include <filesystem>
#include <vector>

using namespace nmr;

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

AmberFlatTableCoverageVerdict FakeUnsupportedVerdict() {
    AmberFlatTableCoverageVerdict v;
    v.ok = false;
    AmberFlatTableCoverageFailure f;
    f.kind = AmberFlatTableCoverageKind::UnsupportedTerminalVariant;
    f.terminal_token = "TEST";
    f.ff_residue_name = "TEST";
    v.failures.push_back(f);
    return v;
}

}  // namespace

class AmberPreparedChargeIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (RuntimeEnvironment::Tleap().empty()) {
            GTEST_SKIP() << "tleap binary not available; AMBER_REQUIRES_TLEAP "
                            "tests SKIPped.";
        }
    }

    AmberSourceConfig MakeCfg() {
        AmberSourceConfig cfg;
        cfg.flat_table_path = nmr::test::TestEnvironment::Ff14sbParams();
        cfg.preparation_policy =
            AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;
        // Leave tleap_path empty; LoadCharges falls back to RuntimeEnvironment.
        return cfg;
    }
};


TEST_F(AmberPreparedChargeIntegrationTest, LoadChargesFailsLoudlyOnInvalidPdbInput) {
    // Build a synthetic ALA protein with a stray atom name that doesn't
    // exist in the AMBER ALA template. tleap will reject the input.
    // LoadCharges must surface that as a named error (referencing the
    // tleap log path), not silent zero-charge rows.
    auto protein = std::make_unique<Protein>();
    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.terminal_state = ResidueTerminalState::Internal;
    size_t ri = protein->AddResidue(res);

    std::vector<Vec3> positions;
    AppendResidueAtomsFromAaType(*protein, AminoAcid::ALA, ri, positions);

    // Inject a stray atom not in the ALA template.
    auto stray = Atom::Create(Element::C);
    stray->pdb_atom_name = "ZZZZ";
    stray->residue_index = ri;
    size_t stray_ai = protein->AddAtom(std::move(stray));
    protein->MutableResidueAt(ri).atom_indices.push_back(stray_ai);
    positions.push_back(Vec3(0.0, 0.0, 0.0));

    protein->AddConformation(std::move(positions), "stray-atom-test");

    AmberPreparedChargeSource src(*protein, MakeCfg().preparation_policy,
                                  FakeUnsupportedVerdict(), MakeCfg());
    std::string err;
    auto rows = src.LoadCharges(*protein, protein->Conformation(), err);

    // The pipeline must NOT silently produce zeros. Either tleap
    // rejects (rc != 0 → "tleap failed") or PRMTOP atom mapping fails
    // (→ "no counterpart"); either error is acceptable. The point is:
    // empty rows + named error.
    EXPECT_TRUE(rows.empty());
    EXPECT_FALSE(err.empty());
    const bool tleap_rejected =
        err.find("tleap failed") != std::string::npos;
    const bool mapping_rejected =
        err.find("no counterpart") != std::string::npos ||
        err.find("ZZZZ") != std::string::npos;
    EXPECT_TRUE(tleap_rejected || mapping_rejected) << err;
}

TEST_F(AmberPreparedChargeIntegrationTest, GeneratedPdbIsByteIdenticalAcrossRuns) {
    // Determinism: building the source twice on the same Protein
    // produces identical PDB bodies. tleap is the source of any
    // numerical noise; our generators are pure functions of typed
    // state.
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(r.Ok()) << r.error;

    AmberPreparedChargeSource src1(
        *r.protein, MakeCfg().preparation_policy,
        FakeUnsupportedVerdict(), MakeCfg());
    AmberPreparedChargeSource src2(
        *r.protein, MakeCfg().preparation_policy,
        FakeUnsupportedVerdict(), MakeCfg());

    const std::string pdb1 = src1.GeneratedPdb(r.protein->Conformation());
    const std::string pdb2 = src2.GeneratedPdb(r.protein->Conformation());
    EXPECT_EQ(pdb1, pdb2);
    EXPECT_GT(pdb1.size(), 0u);
    EXPECT_NE(pdb1.find("END\n"), std::string::npos);
}


// Negative path: tleap binary genuinely missing (no config, no
// RuntimeEnvironment fallback). This test is NOT gated by
// RuntimeEnvironment::Tleap() because it is the test for the
// no-tleap branch.
TEST(AmberPreparedChargeStep5NegativeTest,
     LoadChargesFailsWhenNoTleapBinaryAnywhere) {
    if (!RuntimeEnvironment::Tleap().empty()) {
        GTEST_SKIP() << "RuntimeEnvironment::Tleap() is set; cannot exercise "
                        "the no-tleap branch on this machine.";
    }

    auto protein = std::make_unique<Protein>();
    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.terminal_state = ResidueTerminalState::Internal;
    protein->AddResidue(res);
    std::vector<Vec3> positions;
    AppendResidueAtomsFromAaType(*protein, AminoAcid::ALA, 0, positions);
    protein->AddConformation(std::move(positions), "no-tleap-test");

    AmberSourceConfig cfg;
    cfg.preparation_policy =
        AmberPreparationPolicy::UseCappedFragmentsForUnsupportedTerminalVariants;
    cfg.tleap_path = "/nonexistent/tleap";
    AmberFlatTableCoverageVerdict verdict;
    verdict.ok = false;

    AmberPreparedChargeSource src(*protein, cfg.preparation_policy,
                                  verdict, cfg);
    std::string err;
    auto rows = src.LoadCharges(*protein, protein->Conformation(), err);
    EXPECT_TRUE(rows.empty());
    EXPECT_NE(err.find("no usable tleap binary"), std::string::npos) << err;
}
