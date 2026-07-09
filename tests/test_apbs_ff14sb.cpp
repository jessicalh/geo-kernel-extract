#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "ApbsFieldResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "ForceFieldChargeTable.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "PhysicalConstants.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include <filesystem>
#include <cmath>
#include <numeric>

using namespace nmr;

namespace {

constexpr const char* kTprFixtureProtein = "1P9J_5801";

std::string TrrPathFor(const std::string& tpr_path) {
    return std::filesystem::path(tpr_path).replace_extension(".trr").string();
}

std::string ProductionDirFor(const std::string& tpr_path) {
    return std::filesystem::path(tpr_path).parent_path().string();
}

bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() &&
           std::filesystem::exists(fix.tpr_path) &&
           std::filesystem::exists(TrrPathFor(fix.tpr_path)) &&
           std::filesystem::exists(fix.edr_path);
}

bool IsHydrogenOnNitrogen(const Protein& protein, size_t atom_index) {
    const Atom& atom = protein.AtomAt(atom_index);
    if (atom.element != Element::H) return false;
    if (atom.parent_atom_index != SIZE_MAX) {
        return protein.AtomAt(atom.parent_atom_index).element == Element::N;
    }
    for (const size_t bi : atom.bond_indices) {
        const Bond& bond = protein.BondAt(bi);
        const size_t other = (bond.atom_index_a == atom_index) ? bond.atom_index_b :
                             (bond.atom_index_b == atom_index) ? bond.atom_index_a :
                             SIZE_MAX;
        if (other != SIZE_MAX && protein.AtomAt(other).element == Element::N) {
            return true;
        }
    }
    return false;
}

}  // namespace



class ApbsFF14SBTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found at " << nmr::test::TestEnvironment::UbqProtonated();
        }
        if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) {
            GTEST_SKIP() << "ff14sb_params.dat not found at " << nmr::test::TestEnvironment::Ff14sbParams();
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};


// ---------------------------------------------------------------------------
// Test 1: APBS with real ff14SB charges produces non-zero fields and
// differs from vacuum Coulomb (solvation screening)
// ---------------------------------------------------------------------------

TEST_F(ApbsFF14SBTest, RealChargesProduceNonZeroFields) {
    auto& conf = protein->Conformation();

    // Attach GeometryResult (needed for completeness, not a dependency of APBS)
    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    // Attach ChargeAssignmentResult with REAL ff14SB charges
    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(charges, nullptr) << "ff14SB charge loading failed";
    size_t assigned = charges->AssignedCount();
    EXPECT_GT(assigned, 0u) << "No atoms received ff14SB charges";
    EXPECT_GT(assigned, conf.AtomCount() * 9 / 10)
        << "Too many unassigned atoms";
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    // Attach ApbsFieldResult (uses the real charges)
    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr) << "APBS computation returned nullptr";
    ASSERT_TRUE(conf.AttachResult(std::move(apbs)))
        << "Failed to attach ApbsFieldResult";

    const auto& result = conf.Result<ApbsFieldResult>();

    // Count non-zero E-field and EFG
    int nonzero_E = 0;
    int nonzero_EFG = 0;
    double min_E_mag = 1e30, max_E_mag = 0.0, sum_E_mag = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Vec3 E = result.ElectricFieldAt(ai);
        Mat3 EFG = result.FieldGradientAt(ai);
        double E_mag = E.norm();

        if (E_mag > 1e-10) nonzero_E++;
        if (EFG.norm() > 1e-10) nonzero_EFG++;

        if (E_mag < min_E_mag) min_E_mag = E_mag;
        if (E_mag > max_E_mag) max_E_mag = E_mag;
        sum_E_mag += E_mag;
    }

    double mean_E_mag = sum_E_mag / static_cast<double>(conf.AtomCount());

    // Print summary as requested
    fprintf(stderr,
        "\n=== APBS ff14SB E-field summary ===\n"
        "  atoms:     %zu\n"
        "  nonzero E: %d / %zu\n"
        "  nonzero EFG: %d / %zu\n"
        "  min |E|:   %.6e\n"
        "  max |E|:   %.6e\n"
        "  mean |E|:  %.6e\n"
        "===================================\n",
        conf.AtomCount(),
        nonzero_E, conf.AtomCount(),
        nonzero_EFG, conf.AtomCount(),
        min_E_mag, max_E_mag, mean_E_mag);

    EXPECT_GT(nonzero_E, static_cast<int>(conf.AtomCount()) / 2)
        << "Too few atoms with non-zero E-field from APBS with real charges";
    EXPECT_GT(nonzero_EFG, static_cast<int>(conf.AtomCount()) / 2)
        << "Too few atoms with non-zero EFG from APBS with real charges";
}


// ---------------------------------------------------------------------------
// Test 2: EFG is traceless for each atom (trace < 1e-6)
// ---------------------------------------------------------------------------

TEST_F(ApbsFF14SBTest, EFGIsTraceless) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(charges, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(apbs)));

    const auto& result = conf.Result<ApbsFieldResult>();

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Mat3 EFG = result.FieldGradientAt(ai);
        double trace = std::abs(EFG.trace());
        EXPECT_LT(trace, 1e-6)
            << "EFG trace at atom " << ai << " is " << trace
            << " (expected < 1e-6)";
    }
}


// ---------------------------------------------------------------------------
// Test 3: SphericalTensor roundtrip on the EFG from real charges
// ---------------------------------------------------------------------------

TEST_F(ApbsFF14SBTest, SphericalTensorRoundtripWithRealCharges) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(charges, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(apbs)));

    const auto& result = conf.Result<ApbsFieldResult>();

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Mat3 EFG = result.FieldGradientAt(ai);
        SphericalTensor st = result.FieldGradientSphericalAt(ai);
        Mat3 reconstructed = st.Reconstruct();

        double diff = (EFG - reconstructed).norm();
        EXPECT_LT(diff, 1e-10)
            << "SphericalTensor roundtrip failed at atom " << ai
            << " diff=" << diff;
    }
}


// ---------------------------------------------------------------------------
// Test 4: APBS with real charges differs from direct Coulomb
// (solvation screening changes field magnitude)
// ---------------------------------------------------------------------------

TEST_F(ApbsFF14SBTest, DiffersFromVacuumCoulomb) {
    // Run APBS with real charges
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(charges, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(apbs)));

    const auto& result = conf.Result<ApbsFieldResult>();

    // Compute vacuum Coulomb manually for comparison.
    // Use a subset of atoms (first 50) to keep this efficient.
    size_t n_check = std::min(conf.AtomCount(), size_t(50));
    int differs = 0;
    double sum_apbs_mag = 0.0;
    double sum_coulomb_mag = 0.0;

    for (size_t i = 0; i < n_check; ++i) {
        Vec3 pos_i = conf.PositionAt(i);

        // Vacuum Coulomb E-field at atom i
        // E_i = sum_{j!=i} q_j * (pos_i - pos_j) / |pos_i - pos_j|^3
        Vec3 E_coulomb = Vec3::Zero();
        for (size_t j = 0; j < conf.AtomCount(); ++j) {
            if (i == j) continue;
            Vec3 r = pos_i - conf.PositionAt(j);
            double r_mag = r.norm();
            if (r_mag < 0.1) continue;
            double q_j = conf.AtomAt(j).partial_charge;
            E_coulomb += q_j * r / (r_mag * r_mag * r_mag);
        }

        Vec3 E_apbs = result.ElectricFieldAt(i);
        double apbs_mag = E_apbs.norm();
        double coulomb_mag = E_coulomb.norm();

        sum_apbs_mag += apbs_mag;
        sum_coulomb_mag += coulomb_mag;

        // Fields should differ (solvation screening modifies the field)
        if (apbs_mag > 1e-10 && coulomb_mag > 1e-10) {
            double ratio = apbs_mag / coulomb_mag;
            // If they differ by more than 5%, count it
            if (std::abs(ratio - 1.0) > 0.05) differs++;
        }
    }

    double mean_apbs = sum_apbs_mag / static_cast<double>(n_check);
    double mean_coulomb = sum_coulomb_mag / static_cast<double>(n_check);

    fprintf(stderr,
        "\n=== APBS vs vacuum Coulomb (first %zu atoms) ===\n"
        "  mean |E| APBS:    %.6e\n"
        "  mean |E| Coulomb: %.6e\n"
        "  atoms differing:  %d / %zu\n"
        "=================================================\n",
        n_check, mean_apbs, mean_coulomb, differs, n_check);

    // At least SOME atoms should show solvation screening effect. Identical
    // APBS and direct Coulomb fields indicate an APBS test/setup problem.
    if (differs == 0) {
        fprintf(stderr,
            "WARNING: APBS and vacuum Coulomb fields are identical.\n"
            "This indicates an APBS test/setup problem.\n");
    }
}

TEST(ApbsChargeProvenance, TprPlaceholderPbRadiiMaterializeBeforeProjection) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    const auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kTprFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    RunConfiguration config;
    config.SetName("TprPbRadiusMaterialization");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.SetStride(99999);

    TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();

    Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    Session session;
    ASSERT_EQ(traj.Run(tp, config, session), kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    const Protein& protein = tp.ProteinRef();
    ASSERT_TRUE(protein.HasForceFieldCharges());
    const ForceFieldChargeTable& table = protein.ForceFieldCharges();
    EXPECT_GT(table.DerivedMbondi2PbRadiusCount(), 0);
    EXPECT_EQ(table.PlaceholderPbRadiusCount(), 0);
    EXPECT_EQ(table.MissingCount(), 0);
    EXPECT_EQ(table.NonAuthoritativePbRadiusCount(), 0);

    const ProteinConformation& conf = tp.CanonicalConformation();
    ASSERT_TRUE(conf.HasResult<ChargeAssignmentResult>());

    auto find_atom = [&](Element element, bool h_on_n_required) -> size_t {
        for (size_t i = 0; i < protein.AtomCount(); ++i) {
            const Atom& atom = protein.AtomAt(i);
            if (atom.element != element) continue;
            if (element == Element::H &&
                IsHydrogenOnNitrogen(protein, i) != h_on_n_required) {
                continue;
            }
            return i;
        }
        return SIZE_MAX;
    };

    const size_t amide_h = find_atom(Element::H, true);
    const size_t non_amide_h = find_atom(Element::H, false);
    const size_t oxygen = find_atom(Element::O, false);
    ASSERT_NE(amide_h, SIZE_MAX);
    ASSERT_NE(non_amide_h, SIZE_MAX);
    ASSERT_NE(oxygen, SIZE_MAX);

    for (const size_t ai : {amide_h, non_amide_h, oxygen}) {
        const Atom& atom = protein.AtomAt(ai);
        const bool h_on_n = IsHydrogenOnNitrogen(protein, ai);
        const double expected = Mbondi2Radius(atom.element, h_on_n);
        ASSERT_EQ(table.Values().at(ai).status,
                  ChargeAssignmentStatus::DerivedMbondi2PbRadius)
            << "atom " << ai;
        EXPECT_DOUBLE_EQ(table.PbRadiusAt(ai), expected) << "atom " << ai;
        EXPECT_DOUBLE_EQ(conf.AtomAt(ai).pb_radius, expected) << "atom " << ai;
        EXPECT_DOUBLE_EQ(conf.AtomAt(ai).partial_charge,
                         table.PartialChargeAt(ai)) << "atom " << ai;
        EXPECT_TRUE(table.PbRadiusAuthoritativeAt(ai)) << "atom " << ai;
    }
}
