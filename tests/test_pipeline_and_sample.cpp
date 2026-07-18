#include "TestEnvironment.h"
//
// test_pipeline_and_sample.cpp
//
// Validates the viewer's library surface:
//   1. Pipeline: OperationRunner::Run attaches all results
//   2. Traversal: per-atom fields are populated and consistent
//   3. SampleAt: grid evaluation matches atom-position values
//   4. Legal tensor calculators: SampleKernelAt returns sensible values
//
// This test is the contract between the library and the viewer.
// If it passes, the viewer can use the library without touching it.
//

#include <gtest/gtest.h>
#include <filesystem>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "OperationRunner.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "Of3Loader.h"
#include "OrcaShieldingResult.h"
#include "ForceFieldChargeTable.h"
#include "ProteinBuildContext.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "BiotSavartResult.h"
#include "HaighMallionResult.h"
#include "McConnellResult.h"
#include "RingSusceptibilityResult.h"
#include "PiQuadrupoleResult.h"
#include "DispersionResult.h"
#include "CoulombResult.h"
#include "HBondResult.h"
#include "ChargeSource.h"

namespace fs = std::filesystem;


// Helper: load a protonated protein with charges from prmtop
static std::unique_ptr<nmr::Protein> LoadTestProtein(const std::string& protein_id) {
    std::string dir = std::string(nmr::test::TestEnvironment::Consolidated()) + protein_id + "/";

    nmr::OrcaRunFiles files;
    files.pdb_path = dir + protein_id + "_WT.pdb";
    files.xyz_path = dir + protein_id + "_WT.xyz";
    files.prmtop_path = dir + protein_id + "_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        return nullptr;

    auto load = nmr::BuildFromOrca(files);
    if (!load.Ok()) return nullptr;
    return std::move(load.protein);
}


// BuildFromOrca records the preparation provenance on the Protein's
// build context. Relocated here from the retired test_protonation_pipeline
// (the only coverage of build-context provenance — codex review 2026-05-25).
TEST(OrcaLoaderProvenanceTest, BuildContextRecordsTleapFf14SB) {
    nmr::OrcaRunFiles files;
    const std::string dir = std::string(nmr::test::TestEnvironment::OrcaDir());
    files.pdb_path          = dir + "A0A7C5FAR6_WT.pdb";
    files.xyz_path          = dir + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path       = dir + "A0A7C5FAR6_WT.prmtop";
    files.tleap_script_path = dir + "A0A7C5FAR6_WT_tleap.in";
    if (!fs::exists(files.prmtop_path)) GTEST_SKIP() << "ORCA test data not found";

    auto result = nmr::BuildFromOrca(files);
    ASSERT_TRUE(result.Ok()) << result.error;

    const auto& ctx = result.protein->BuildContext();
    EXPECT_FALSE(ctx.pdb_source.empty());
    EXPECT_EQ(ctx.force_field, "ff14SB");
    EXPECT_EQ(ctx.protonation_tool, "tleap");
    EXPECT_FALSE(ctx.prmtop_path.empty());

    // Regression: parameterizing the shared heavy parse (LoadWithPrmtop) for
    // --of3 must NOT change ORCA's provenance. The ORCA conformation still
    // records the "AlphaFold+tleap" convention, unchanged.
    EXPECT_EQ(result.protein->PredictionAt(0).PredictionMethod(),
              "AlphaFold+tleap");
}


// BuildFromOf3 is the --of3 structural loader: it shares the ORCA heavy parse
// (prmtop topology + Branch-1 prmtop charges) but reads geometry from the
// tleap-emitted AMBER restart (rst7/inpcrd), records honest OpenFold
// provenance, and attaches NO DFT result. This round trip uses a committed
// REAL tleap prep (bmr10013.{prmtop,inpcrd}, tracked in git under tests/data/
// orca/) — genuine instrument output, not a hand-generated fixture.
TEST(Of3LoaderRoundTrip, StructureChargesProvenanceNoDft) {
    const std::string dir = std::string(nmr::test::TestEnvironment::OrcaDir());
    nmr::Of3Input input;
    input.prmtop_path       = dir + "bmr10013.prmtop";
    input.inpcrd_path       = dir + "bmr10013.inpcrd";
    input.prediction_method = "OpenFold+tleap";
    if (!fs::exists(input.prmtop_path) || !fs::exists(input.inpcrd_path))
        GTEST_SKIP() << "OF3 inpcrd fixture not found at " << dir;

    auto build = nmr::BuildFromOf3(input);
    ASSERT_TRUE(build.Ok()) << build.error;

    // Structure: prmtop topology + rst7 geometry reached the conformation in
    // prmtop atom order (proves ReadInpcrd parsed the 6F12.7 coordinates).
    EXPECT_EQ(build.protein->AtomCount(), 1284u);
    EXPECT_EQ(build.protein->ResidueCount(), 89u);
    EXPECT_EQ(build.protein->ConformationCount(), 1u);

    const auto& conf = build.protein->Conformation();
    EXPECT_NEAR(conf.AtomAt(0).Position().x(),  -0.675000, 1e-4);
    EXPECT_NEAR(conf.AtomAt(0).Position().y(),  -1.661000, 1e-4);
    EXPECT_NEAR(conf.AtomAt(0).Position().z(), -20.494000, 1e-4);
    const size_t last = build.protein->AtomCount() - 1;
    EXPECT_NEAR(conf.AtomAt(last).Position().x(),  -0.041798, 1e-4);
    EXPECT_NEAR(conf.AtomAt(last).Position().y(), -18.703585, 1e-4);
    EXPECT_NEAR(conf.AtomAt(last).Position().z(),  27.816413, 1e-4);

    // Provenance: ff14SB/tleap build context, exact prmtop path, honest
    // OpenFold prediction method (NOT AlphaFold).
    const auto& ctx = build.protein->BuildContext();
    EXPECT_EQ(ctx.prmtop_path, input.prmtop_path);
    EXPECT_EQ(ctx.force_field, "ff14SB");
    EXPECT_EQ(ctx.protonation_tool, "tleap");
    EXPECT_EQ(build.protein->PredictionAt(0).PredictionMethod(),
              "OpenFold+tleap");

    // Charges: Branch-1 prmtop authority, ff14SB, materialized onto the
    // protein's owned table.
    ASSERT_NE(build.charges, nullptr);
    EXPECT_EQ(build.charges->Kind(), nmr::ChargeModelKind::AmberPrmtop);
    EXPECT_EQ(build.charges->SourceForceField(), nmr::ForceField::Amber_ff14SB);
    ASSERT_TRUE(build.protein->HasForceFieldCharges());
    EXPECT_EQ(build.protein->ForceFieldCharges().Kind(),
              nmr::ChargeModelKind::AmberPrmtop);
    EXPECT_EQ(build.protein->ForceFieldCharges().AtomCount(), 1284u);
    EXPECT_EQ(build.net_charge, -9);

    // No DFT attached: OF3 cannot carry an nmr path (see the HasNmrOutPath
    // static_asserts in test_cli_parse.cpp) and the loader builds no
    // OrcaShieldingResult.
    EXPECT_FALSE(conf.HasResult<nmr::OrcaShieldingResult>());
}

// Honest-provenance guard: BuildFromOf3 refuses an empty prediction_method
// rather than silently recording blank provenance.
TEST(Of3LoaderRoundTrip, EmptyPredictionMethodRejected) {
    const std::string dir = std::string(nmr::test::TestEnvironment::OrcaDir());
    nmr::Of3Input input;
    input.prmtop_path       = dir + "bmr10013.prmtop";
    input.inpcrd_path       = dir + "bmr10013.inpcrd";
    input.prediction_method = "";  // empty on purpose
    if (!fs::exists(input.prmtop_path)) GTEST_SKIP() << "OF3 fixture not found";

    auto build = nmr::BuildFromOf3(input);
    EXPECT_FALSE(build.Ok());
    EXPECT_NE(build.error.find("prediction_method"), std::string::npos)
        << build.error;
}


// ============================================================================
// ReadInpcrd fail-loud hardening (adversarial Finding #2): the rst7 reader must
// REJECT wrong geometry rather than report success on it. These craft malformed
// rst7 files and drive them through BuildFromOf3 (ReadInpcrd is a file-local
// static; the public entry is the loader). Failure cases are rejected inside
// ReadInpcrd, before the prmtop is even consulted, so the diagnostic names the
// inpcrd reason. bmr10013.prmtop is supplied only so a REGRESSED reader (that
// wrongly succeeded) would still fail loudly on the count mismatch.
// ============================================================================
namespace {

namespace of3fs = std::filesystem;

std::string Of3OrcaDir() { return std::string(nmr::test::TestEnvironment::OrcaDir()); }

// One F12.7 coordinate field (right-justified, exactly 12 chars for |v|<1000).
std::string F12(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%12.7f", v);
    return std::string(buf);
}

// Right-justify a raw token into a 12-char field (for injecting non-numeric /
// non-finite fields). Returns tok unchanged if already >= 12 (over-width).
std::string Field12(const std::string& tok) {
    if (tok.size() >= 12) return tok;
    return std::string(12 - tok.size(), ' ') + tok;
}

// Write a crafted rst7: title, atom count, then the given lines verbatim.
of3fs::path WriteRst7(const std::string& tag, int natom,
                      const std::vector<std::string>& lines) {
    const of3fs::path dir = of3fs::path(testing::TempDir()) /
        ("of3_rst7_" + tag + "_" + std::to_string(::getpid()));
    of3fs::create_directories(dir);
    const of3fs::path p = dir / "model.inpcrd";
    std::ofstream out(p);
    out << "title\n" << natom << "\n";
    for (const auto& l : lines) out << l << "\n";
    return p;
}

nmr::Of3Input CraftedOf3(const of3fs::path& inpcrd) {
    nmr::Of3Input in;
    in.inpcrd_path       = inpcrd.string();
    in.prmtop_path       = Of3OrcaDir() + "bmr10013.prmtop";  // only a backstop
    in.prediction_method = "OpenFold+tleap";
    return in;
}

}  // namespace

TEST(Of3InpcrdHardening, OverWidthColumnRejected) {
    // 2 atoms = 6 coords on one line; field 0 is 13 chars (F12.7 overflow, e.g.
    // a negative value <= -1000), which shifts every fixed 12-char column.
    const std::string over = "-1234.5678901";  // 13 chars
    const std::string line = over + F12(0.1) + F12(0.2) + F12(0.3) + F12(0.4) + F12(0.5);
    const auto b = nmr::BuildFromOf3(CraftedOf3(WriteRst7("overwidth", 2, {line})));
    EXPECT_FALSE(b.Ok());
    EXPECT_NE(b.error.find("inpcrd"), std::string::npos) << b.error;
    EXPECT_NE(b.error.find("width"),  std::string::npos) << b.error;
}

TEST(Of3InpcrdHardening, NonFiniteRejected) {
    // "nan" parses numerically but is not finite.
    const std::string line = Field12("nan") + F12(0.1) + F12(0.2) +
                             F12(0.3) + F12(0.4) + F12(0.5);
    const auto b = nmr::BuildFromOf3(CraftedOf3(WriteRst7("nan", 2, {line})));
    EXPECT_FALSE(b.Ok());
    EXPECT_NE(b.error.find("inpcrd"), std::string::npos) << b.error;
    EXPECT_NE(b.error.find("finite"), std::string::npos) << b.error;
}

TEST(Of3InpcrdHardening, PartialNumericParseRejected) {
    // Fortran "D"-exponent: std::stod reads "1.234" and stops at 'D', leaving a
    // partial parse the naive reader would silently accept.
    const std::string line = Field12("1.234D+02") + F12(0.1) + F12(0.2) +
                             F12(0.3) + F12(0.4) + F12(0.5);
    const auto b = nmr::BuildFromOf3(CraftedOf3(WriteRst7("dexp", 2, {line})));
    EXPECT_FALSE(b.Ok());
    EXPECT_NE(b.error.find("inpcrd"),  std::string::npos) << b.error;
    EXPECT_NE(b.error.find("numeric"), std::string::npos) << b.error;
}

TEST(Of3InpcrdHardening, TrailingGarbageParseRejected) {
    // A numeric prefix with trailing non-numeric junk ("1.0000000gx").
    const std::string line = Field12("1.0000000gx") + F12(0.1) + F12(0.2) +
                             F12(0.3) + F12(0.4) + F12(0.5);
    const auto b = nmr::BuildFromOf3(CraftedOf3(WriteRst7("garbage", 2, {line})));
    EXPECT_FALSE(b.Ok());
    EXPECT_NE(b.error.find("inpcrd"),  std::string::npos) << b.error;
    EXPECT_NE(b.error.find("numeric"), std::string::npos) << b.error;
}

TEST(Of3InpcrdHardening, ShortRecordRejected) {
    // 2 atoms need 6 fields; provide only 5 (60 chars) -> width mismatch.
    const std::string line = F12(0.1) + F12(0.2) + F12(0.3) + F12(0.4) + F12(0.5);
    const auto b = nmr::BuildFromOf3(CraftedOf3(WriteRst7("short", 2, {line})));
    EXPECT_FALSE(b.Ok());
    EXPECT_NE(b.error.find("inpcrd"), std::string::npos) << b.error;
    EXPECT_NE(b.error.find("width"),  std::string::npos) << b.error;
}

TEST(Of3InpcrdHardening, CrossBlockFillRejected) {
    // 4 atoms = 12 coords over ceil(4/2)=2 lines. Line 2 is SHORT (atom 3 only,
    // 3 fields = 36 chars); a following velocity/box line would, in a reader
    // that just accumulates until 3*N, silently supply atom 4's coords from the
    // wrong block. The hardened reader rejects the short record and never reads
    // the later line.
    const std::string full  = F12(1)+F12(2)+F12(3)+F12(4)+F12(5)+F12(6);
    const std::string shrt  = F12(7)+F12(8)+F12(9);                  // atom 3 only
    const std::string velo  = F12(70)+F12(80)+F12(90)+F12(100)+F12(110)+F12(120);
    const auto b = nmr::BuildFromOf3(
        CraftedOf3(WriteRst7("crossblock", 4, {full, shrt, velo})));
    EXPECT_FALSE(b.Ok());
    EXPECT_NE(b.error.find("inpcrd"), std::string::npos) << b.error;
    EXPECT_NE(b.error.find("width"),  std::string::npos) << b.error;
    // The velocity-block magnitudes (70..120) must never have been read as coords.
    EXPECT_EQ(b.error.find("70.0000000"), std::string::npos) << b.error;
}

TEST(Of3InpcrdHardening, TrailingBoxBlockIgnoredValidSucceeds) {
    // A COMPLETE, valid coordinate block followed by a box line must load
    // correctly — the trailing block is ignored, not read as coordinates.
    const std::string dir = Of3OrcaDir();
    const std::string src = dir + "bmr10013.inpcrd";
    if (!of3fs::exists(src)) GTEST_SKIP() << "bmr10013 fixture not found";

    std::ifstream in(src);
    std::stringstream ss; ss << in.rdbuf();
    std::string content = ss.str();
    if (!content.empty() && content.back() != '\n') content.push_back('\n');
    content += "  10.0000000  10.0000000  10.0000000"
               "  90.0000000  90.0000000  90.0000000\n";  // box line
    const of3fs::path p = of3fs::path(testing::TempDir()) /
        ("of3_box_" + std::to_string(::getpid()) + ".inpcrd");
    { std::ofstream out(p); out << content; }

    nmr::Of3Input input;
    input.inpcrd_path       = p.string();
    input.prmtop_path       = dir + "bmr10013.prmtop";
    input.prediction_method = "OpenFold+tleap";
    auto b = nmr::BuildFromOf3(input);
    ASSERT_TRUE(b.Ok()) << b.error;
    EXPECT_EQ(b.protein->AtomCount(), 1284u);
    const auto& conf = b.protein->Conformation();
    EXPECT_NEAR(conf.AtomAt(0).Position().x(), -0.675000, 1e-4);
    const size_t last = b.protein->AtomCount() - 1;
    EXPECT_NEAR(conf.AtomAt(last).Position().z(), 27.816413, 1e-4);
    EXPECT_EQ(b.net_charge, -9);
}


// ============================================================================
// Test 1: Pipeline attaches all results in correct order
// ============================================================================

TEST(OperationRunnerTest, AttachesAllResults) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();

    nmr::PrmtopChargeSource prmtop(
        std::string(nmr::test::TestEnvironment::Consolidated()) + "P84477/P84477_WT.prmtop",
        nmr::ForceField::Amber_ff14SB);
    nmr::RunOptions opts;
    opts.charge_source = &prmtop;

    auto run_result = nmr::OperationRunner::Run(conf, opts);

    // Should have: Geometry, SpatialIndex, Enrichment, Dssp,
    // ChargeAssignment, BS, HM, MC, RingSusc scalar, PQ scalar,
    // Disp scalar, Coulomb, HBond scalars
    EXPECT_GE(run_result.attached.size(), 13u)
        << "Expected 13+ results, got " << run_result.attached.size();

    // Verify key results are attached
    EXPECT_TRUE(conf.HasResult<nmr::BiotSavartResult>());
    EXPECT_TRUE(conf.HasResult<nmr::HaighMallionResult>());
    EXPECT_TRUE(conf.HasResult<nmr::McConnellResult>());
    EXPECT_TRUE(conf.HasResult<nmr::CoulombResult>());
    EXPECT_TRUE(conf.HasResult<nmr::RingSusceptibilityResult>());
    EXPECT_TRUE(conf.HasResult<nmr::PiQuadrupoleResult>());
    EXPECT_TRUE(conf.HasResult<nmr::DispersionResult>());
    EXPECT_TRUE(conf.HasResult<nmr::HBondResult>());
}


// ============================================================================
// Test 2: Per-atom traversal produces populated fields
// ============================================================================

TEST(OperationRunnerTest, AtomFieldsPopulated) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();

    nmr::PrmtopChargeSource prmtop(
        std::string(nmr::test::TestEnvironment::Consolidated()) + "P84477/P84477_WT.prmtop",
        nmr::ForceField::Amber_ff14SB);
    nmr::RunOptions opts;
    opts.charge_source = &prmtop;

    nmr::OperationRunner::Run(conf, opts);

    // Traverse atoms — the viewer pattern
    int atoms_with_bs = 0;
    int atoms_with_mc = 0;
    int atoms_with_coulomb = 0;
    int atoms_with_rings = 0;
    double max_bs_t0 = 0;
    double max_mc_t0 = 0;

    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& atom = conf.AtomAt(i);

        // Position is non-zero (loaded correctly)
        EXPECT_GT(atom.Position().norm(), 0.0);

        // BS shielding
        if (std::abs(atom.bs_shielding_contribution.T0) > 1e-10) {
            atoms_with_bs++;
            max_bs_t0 = std::max(max_bs_t0,
                std::abs(atom.bs_shielding_contribution.T0));
        }

        // McConnell shielding
        if (std::abs(atom.mc_shielding_contribution.T0) > 1e-10) {
            atoms_with_mc++;
            max_mc_t0 = std::max(max_mc_t0,
                std::abs(atom.mc_shielding_contribution.T0));
        }

        // Coulomb E-field
        if (atom.coulomb_E_total.norm() > 1e-10) {
            atoms_with_coulomb++;
        }

        // Ring neighbourhood
        if (!atom.ring_neighbours.empty()) {
            atoms_with_rings++;
            // Check that ring neighbourhood has BS and HM tensors
            for (const auto& rn : atom.ring_neighbours) {
                EXPECT_GT(rn.distance_to_center, 0.0);
                // G_tensor should be nonzero for nearby rings
                if (rn.distance_to_center < 10.0) {
                    EXPECT_GT(rn.G_tensor.norm(), 0.0)
                        << "BS G_tensor zero at dist=" << rn.distance_to_center;
                }
            }
        }
    }

    // P84477 is a small protein with aromatic rings
    EXPECT_GT(atoms_with_bs, 0) << "No atoms with BS shielding";
    EXPECT_GT(atoms_with_mc, 0) << "No atoms with MC shielding";
    EXPECT_GT(atoms_with_coulomb, 0) << "No atoms with Coulomb E-field";
    EXPECT_GT(atoms_with_rings, 0) << "No atoms with ring neighbours";

    std::cout << "  P84477 traversal summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << protein->RingCount() << "\n"
              << "    with BS: " << atoms_with_bs
              << " (max T0=" << max_bs_t0 << ")\n"
              << "    with MC: " << atoms_with_mc
              << " (max T0=" << max_mc_t0 << ")\n"
              << "    with Coulomb: " << atoms_with_coulomb << "\n"
              << "    with ring neighbours: " << atoms_with_rings << "\n";
}


// ============================================================================
// Test 3: SampleAt agrees with atom-position results
// ============================================================================

TEST(SampleAtTest, BSMatchesAtomValues) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::OperationRunner::Run(conf, {});

    const auto& bs = conf.Result<nmr::BiotSavartResult>();

    // For atoms near rings, SampleKernelAt should approximately match
    // the stored bs_shielding_contribution. Not exact — the atom-position
    // result uses filters (RingBondedExclusion) that SampleAt doesn't.
    // But for atoms NOT bonded to ring vertices, they should be close.
    int tested = 0;
    int matched = 0;

    for (size_t i = 0; i < conf.AtomCount() && tested < 50; ++i) {
        const auto& atom = conf.AtomAt(i);
        double atom_t0 = atom.bs_shielding_contribution.T0;
        if (std::abs(atom_t0) < 0.001) continue;  // skip negligible

        // Skip atoms that are ring vertices or bonded to ring vertices
        // (filters differ between Compute and SampleAt)
        bool is_ring_atom = false;
        for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
            const auto& ring = protein->RingAt(ri);
            for (size_t vi : ring.atom_indices) {
                if (vi == i) { is_ring_atom = true; break; }
                // Check if bonded to vertex
                for (const auto& bond : protein->Bonds()) {
                    if ((bond.atom_index_a == i && bond.atom_index_b == vi) ||
                        (bond.atom_index_a == vi && bond.atom_index_b == i)) {
                        is_ring_atom = true;
                        break;
                    }
                }
                if (is_ring_atom) break;
            }
            if (is_ring_atom) break;
        }
        if (is_ring_atom) continue;

        nmr::SphericalTensor sampled = bs.SampleKernelAt(atom.Position());
        tested++;

        // Relative tolerance: these should match well for non-excluded atoms
        double rel_diff = std::abs(sampled.T0 - atom_t0)
                        / std::max(std::abs(atom_t0), 1e-6);
        if (rel_diff < 0.05) matched++;  // 5% tolerance

        EXPECT_NEAR(sampled.T0, atom_t0, std::abs(atom_t0) * 0.05 + 1e-6)
            << "BS SampleAt mismatch at atom " << i
            << " sampled=" << sampled.T0 << " stored=" << atom_t0;
    }

    EXPECT_GT(tested, 5) << "Too few testable atoms";
    std::cout << "  BS SampleAt: tested=" << tested
              << " matched=" << matched << "\n";
}


// ============================================================================
// Test 4: SampleAt on a grid produces physically sensible values
// ============================================================================

TEST(SampleAtTest, BSGridAboveRing) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::OperationRunner::Run(conf, {});

    if (protein->RingCount() == 0) GTEST_SKIP() << "No rings";

    const auto& bs = conf.Result<nmr::BiotSavartResult>();
    const auto& geom = conf.ring_geometries[0];

    // Sample on the ring normal at 3A above center — should be shielded (T0 > 0)
    nmr::Vec3 above = geom.center + 3.0 * geom.normal;
    auto st_above = bs.SampleKernelAt(above);

    // And in the ring plane at 5A — should be deshielded (T0 < 0)
    // Find a direction perpendicular to normal
    nmr::Vec3 perp = geom.normal.cross(
        std::abs(geom.normal.x()) < 0.9 ? nmr::Vec3(1,0,0) : nmr::Vec3(0,1,0)
    ).normalized();
    nmr::Vec3 inplane = geom.center + 5.0 * perp;
    auto st_inplane = bs.SampleKernelAt(inplane);

    std::cout << "  Ring 0: T0 at +3A normal = " << st_above.T0
              << ", T0 at 5A in-plane = " << st_inplane.T0 << "\n";

    // Classic ring current pattern: shielded above, deshielded in-plane
    // (with negative intensity, the signs flip, but the pattern holds)
    EXPECT_NE(st_above.T0, 0.0) << "Zero shielding above ring";
    // The above and in-plane should have opposite signs
    if (std::abs(st_above.T0) > 1e-4 && std::abs(st_inplane.T0) > 1e-4) {
        EXPECT_LT(st_above.T0 * st_inplane.T0, 0.0)
            << "Above and in-plane should have opposite T0 signs";
    }
}


// ============================================================================
// Test 5: B-field sampling for butterfly visualization
// ============================================================================

TEST(SampleAtTest, BSButterflyField) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::OperationRunner::Run(conf, {});

    if (protein->RingCount() == 0) GTEST_SKIP() << "No rings";

    const auto& bs = conf.Result<nmr::BiotSavartResult>();
    const auto& geom = conf.ring_geometries[0];

    // B-field on the ring axis should be parallel to the normal
    nmr::Vec3 above = geom.center + 3.0 * geom.normal;
    nmr::Vec3 B = bs.SampleBFieldAt(above);

    EXPECT_GT(B.norm(), 0.0) << "Zero B-field above ring";

    // B should be roughly along the normal direction
    double cos_angle = B.normalized().dot(geom.normal);
    EXPECT_GT(std::abs(cos_angle), 0.8)
        << "B-field not aligned with ring normal at 3A above";

    std::cout << "  B-field at +3A: |B|=" << B.norm()
              << " cos(B,n)=" << cos_angle << "\n";
}


// ============================================================================
// Test 6: Legal tensor SampleAt methods return non-zero for appropriate queries
// ============================================================================

TEST(SampleAtTest, AllCalculatorsSample) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::PrmtopChargeSource prmtop(
        std::string(nmr::test::TestEnvironment::Consolidated()) + "P84477/P84477_WT.prmtop",
        nmr::ForceField::Amber_ff14SB);
    nmr::RunOptions opts;
    opts.charge_source = &prmtop;
    nmr::OperationRunner::Run(conf, opts);

    if (protein->RingCount() == 0) GTEST_SKIP() << "No rings";

    const auto& geom = conf.ring_geometries[0];
    nmr::Vec3 test_point = geom.center + 3.0 * geom.normal;

    // HM
    auto hm_st = conf.Result<nmr::HaighMallionResult>().SampleKernelAt(test_point);
    EXPECT_NE(hm_st.T0, 0.0) << "HM SampleAt returned zero";

    // McConnell — sample near a bond midpoint
    nmr::Vec3 bond_test = conf.bond_midpoints[0] + nmr::Vec3(2.0, 0.0, 0.0);
    auto mc_st = conf.Result<nmr::McConnellResult>().SampleKernelAt(bond_test);
    // McConnell is pure T2 (T0 ≈ 0), so check T2 magnitude
    EXPECT_GT(mc_st.T2Magnitude(), 0.0) << "MC SampleAt returned zero T2";

    // Coulomb E-field
    auto E = conf.Result<nmr::CoulombResult>().SampleEFieldAt(test_point);
    EXPECT_GT(E.norm(), 0.0) << "Coulomb SampleEFieldAt returned zero";

    // Scalar-only calculators still attach and populate retained fields.
    const auto& atom0 = conf.AtomAt(0);
    (void)conf.Result<nmr::RingSusceptibilityResult>();
    (void)conf.Result<nmr::PiQuadrupoleResult>();
    (void)conf.Result<nmr::DispersionResult>();
    (void)conf.Result<nmr::HBondResult>();
    if (!atom0.ring_neighbours.empty()) {
        (void)atom0.ring_neighbours.front().chi_scalar;
        (void)atom0.ring_neighbours.front().quad_scalar;
        (void)atom0.ring_neighbours.front().disp_scalar;
    }
    (void)atom0.hbond_mcconnell_scalar;

    std::cout << "  All SampleAt methods exercised successfully\n"
              << "    BS T0=" << conf.Result<nmr::BiotSavartResult>().SampleKernelAt(test_point).T0 << "\n"
              << "    HM T0=" << hm_st.T0 << "\n"
              << "    MC T2_mag=" << mc_st.T2Magnitude() << "\n"
              << "    |E|=" << E.norm() << " V/A\n";
}


// ============================================================================
// Test 7: Pipeline without charges gracefully skips Coulomb
// ============================================================================

TEST(OperationRunnerTest, SkipsCoulombWithoutCharges) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();

    // No charge source — Coulomb should be skipped
    auto run_result = nmr::OperationRunner::Run(conf, {});

    EXPECT_FALSE(conf.HasResult<nmr::CoulombResult>())
        << "Coulomb should not be attached without charges";
    EXPECT_TRUE(conf.HasResult<nmr::BiotSavartResult>())
        << "BS should still be attached";
}


// ============================================================================
// Test 8: Ring and bond geometry accessible for viewer rendering
// ============================================================================

TEST(OperationRunnerTest, GeometryAccessible) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::OperationRunner::Run(conf, {});

    // Ring geometries
    EXPECT_EQ(conf.ring_geometries.size(), protein->RingCount());
    for (size_t i = 0; i < protein->RingCount(); ++i) {
        const auto& geo = conf.ring_geometries[i];
        EXPECT_GT(geo.radius, 0.0);
        EXPECT_NEAR(geo.normal.norm(), 1.0, 1e-6);
        EXPECT_GT(geo.vertices.size(), 2u);
    }

    // Bond geometries
    EXPECT_EQ(conf.bond_midpoints.size(), protein->BondCount());
    EXPECT_EQ(conf.bond_directions.size(), protein->BondCount());
    EXPECT_EQ(conf.bond_lengths.size(), protein->BondCount());

    // Global geometry
    EXPECT_GT(conf.radius_of_gyration, 0.0);

    std::cout << "  Geometry: "
              << protein->RingCount() << " rings, "
              << protein->BondCount() << " bonds, "
              << "Rg=" << conf.radius_of_gyration << " A\n";
}
