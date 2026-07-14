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
#include <cstdlib>
#include <filesystem>
#include <vector>

using namespace nmr;

namespace {

size_t AppendFinalizableAlaAtoms(Protein& protein,
                                size_t residue_index,
                                std::vector<Vec3>& positions_out) {
    // Coordinates are residue 28 of the protonated 1UBQ fixture, in the
    // canonical AminoAcidType ALA atom order. Keeping real local geometry
    // lets this full-pipeline failure fixture cross the topology boundary
    // before its outbound name is deliberately corrupted.
    const std::vector<Vec3> positions = {
        Vec3(37.499, 25.743, 14.571),  // N
        Vec3(38.794, 25.761, 13.880),  // CA
        Vec3(38.728, 26.591, 12.611),  // C
        Vec3(39.704, 27.346, 12.277),  // O
        Vec3(37.122, 24.970, 14.592),  // H
        Vec3(39.434, 26.174, 14.481),  // HA
        Vec3(39.285, 24.336, 13.566),  // CB
        Vec3(40.140, 24.381, 13.111),  // HB1
        Vec3(39.385, 23.838, 14.393),  // HB2
        Vec3(38.639, 23.889, 12.996),  // HB3
    };

    const AminoAcidType& aa_type = GetAminoAcidType(AminoAcid::ALA);
    if (aa_type.atoms.size() != positions.size()) {
        std::abort();
    }

    size_t cb_atom_index = 0;
    for (size_t i = 0; i < aa_type.atoms.size(); ++i) {
        const auto& templ = aa_type.atoms[i];
        auto atom = Atom::Create(templ.element);
        atom->pdb_atom_name = templ.name;
        atom->residue_index = residue_index;
        const size_t ai = protein.AddAtom(std::move(atom));
        protein.MutableResidueAt(residue_index).atom_indices.push_back(ai);
        positions_out.push_back(positions[i]);
        if (std::string(templ.name) == "CB") {
            cb_atom_index = ai;
        }
    }
    return cb_atom_index;
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
    const size_t cb_ai = AppendFinalizableAlaAtoms(*protein, ri, positions);
    protein->FinalizeConstruction(positions);
    protein->AddConformation(std::move(positions), "stray-atom-test");

    // Corrupt one outbound atom name only after the Protein has crossed
    // the topology boundary. LoadCharges must still surface tleap's named
    // failure rather than silently producing zero-charge rows.
    protein->MutableAtomAt(cb_ai).pdb_atom_name = "ZZZZ";

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
