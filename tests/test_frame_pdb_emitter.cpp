//
// test_frame_pdb_emitter -- discipline tests for the opt-in trajectory
// frame-PDB writer. Mirrors the per-TR test pattern in
// test_amber_streaming.cpp (fleet_amber/1P9J_5801 fixture) but tests
// the emitter directly: not-configured-no-op, stride/window gating,
// decorator filename, content sanity.
//
// Singleton state is shared across test cases; the fixture calls
// FramePdbEmitter::Reset() in SetUp + TearDown to keep cases isolated.
//

#include "AIMNet2Result.h"

// Nanoflann-using headers (SpatialIndexResult via BiotSavartResult)
// must precede gromacs headers — gromacs vectypes.h #defines DIM 3
// which collides with nanoflann's DIM template parameter.
#include "GeometryResult.h"
#include "SpatialIndexResult.h"

#include "FramePdbEmitter.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "Protein.h"
#include "Bond.h"
#include "LegacyAmberTopology.h"
#include "CalculatorConfig.h"
#include "OperationLog.h"
#include "TestEnvironment.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

const std::string kFixtureProtein = "1P9J_5801";

std::string TrrPathFor(const std::string& tpr_path) {
    return fs::path(tpr_path).replace_extension(".trr").string();
}

std::string ProductionDirFor(const std::string& tpr_path) {
    return fs::path(tpr_path).parent_path().string();
}

bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() &&
           fs::exists(fix.tpr_path) &&
           fs::exists(TrrPathFor(fix.tpr_path)) &&
           fs::exists(fix.edr_path);
}

// Minimal config: GeometryResult + SpatialIndexResult only, calculators
// off. Same pattern as test_amber_streaming TrajectoryBuildAndRun.
void ConfigureMinimalRun(nmr::RunConfiguration& config,
                         const std::string& test_name) {
    config.SetName(test_name);
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
}

void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}

// Make a deterministic per-test temp directory under /tmp; clear it
// first so prior runs don't pollute counts.
fs::path FreshTempDir(const std::string& test_name) {
    fs::path dir = fs::temp_directory_path() /
                   ("frame_pdb_emitter_test_" + test_name + "_" +
                    std::to_string(::getpid()));
    fs::remove_all(dir);
    return dir;
}

size_t CountPdbsIn(const fs::path& dir) {
    if (!fs::exists(dir)) return 0;
    size_t n = 0;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (e.path().extension() == ".pdb") ++n;
    }
    return n;
}

}  // namespace


// ============================================================================
// Test fixture: enforces singleton-Reset isolation between cases.
// ============================================================================

class FramePdbEmitterTest : public ::testing::Test {
protected:
    void SetUp() override {
        nmr::FramePdbEmitter::Reset();
        LoadCalculatorConfig();
    }
    void TearDown() override {
        nmr::FramePdbEmitter::Reset();
    }
};


// ============================================================================
// Without Configure(), OnFrame is a no-op even when Trajectory::Run
// drives the per-frame loop.
// ============================================================================

TEST_F(FramePdbEmitterTest, NotConfiguredEmitsNothing) {
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    fs::path outdir = FreshTempDir("not_configured");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();

    // Deliberately do NOT call FramePdbEmitter::Configure.
    EXPECT_FALSE(nmr::FramePdbEmitter::IsActive());

    nmr::RunConfiguration config;
    ConfigureMinimalRun(config, "FramePdbEmitterNotConfigured");
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    EXPECT_EQ(CountPdbsIn(outdir), 0u)
        << "no Configure -> no PDBs anywhere";
}


// ============================================================================
// Configured with stride=1: emits one PDB per dispatched frame.
// ============================================================================

TEST_F(FramePdbEmitterTest, ConfiguredEmitsOnePerFrame) {
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    fs::path outdir = FreshTempDir("emit_per_frame");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();

    nmr::FramePdbEmitter::Config cfg;
    cfg.output_dir = outdir;
    cfg.stem       = "p1p9j";
    cfg.stride     = 1;
    nmr::FramePdbEmitter::Configure(tp.ProteinRef(), cfg);
    EXPECT_TRUE(nmr::FramePdbEmitter::IsActive());

    nmr::RunConfiguration config;
    ConfigureMinimalRun(config, "FramePdbEmitterConfiguredEmits");
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    EXPECT_GT(traj.FrameCount(), 1u);
    EXPECT_EQ(CountPdbsIn(outdir), traj.FrameCount())
        << "stride=1 -> one PDB per frame dispatched";
}


// ============================================================================
// Stride larger than fixture length leaves only frame 0 emitted.
// ============================================================================

TEST_F(FramePdbEmitterTest, StrideGatesToFrameZero) {
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    fs::path outdir = FreshTempDir("stride_gates");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();

    nmr::FramePdbEmitter::Config cfg;
    cfg.output_dir = outdir;
    cfg.stem       = "p1p9j";
    cfg.stride     = 99999;   // > fixture length: only frame 0 hits
    nmr::FramePdbEmitter::Configure(tp.ProteinRef(), cfg);

    nmr::RunConfiguration config;
    ConfigureMinimalRun(config, "FramePdbEmitterStrideGate");
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    EXPECT_EQ(CountPdbsIn(outdir), 1u)
        << "stride beyond fixture length must leave only frame 0";
}


// ============================================================================
// Time window: from_ps so high that no frame's time clears it.
// ============================================================================

TEST_F(FramePdbEmitterTest, WindowGatesNothing) {
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    fs::path outdir = FreshTempDir("window_gates");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();

    nmr::FramePdbEmitter::Config cfg;
    cfg.output_dir = outdir;
    cfg.stem       = "p1p9j";
    cfg.stride     = 1;
    cfg.from_ps    = 1.0e12;  // far in the future
    nmr::FramePdbEmitter::Configure(tp.ProteinRef(), cfg);

    nmr::RunConfiguration config;
    ConfigureMinimalRun(config, "FramePdbEmitterWindowGate");
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    EXPECT_EQ(CountPdbsIn(outdir), 0u)
        << "from_ps beyond every frame's time must gate them all out";
}


// ============================================================================
// Decorator appears in the emitted filenames.
// ============================================================================

TEST_F(FramePdbEmitterTest, DecoratorInFilename) {
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    fs::path outdir = FreshTempDir("decorator");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();

    nmr::FramePdbEmitter::Config cfg;
    cfg.output_dir = outdir;
    cfg.stem       = "p1p9j";
    cfg.decorator  = "chain-test";
    cfg.stride     = 99999;   // emit only frame 0 to keep test cheap
    nmr::FramePdbEmitter::Configure(tp.ProteinRef(), cfg);

    nmr::RunConfiguration config;
    ConfigureMinimalRun(config, "FramePdbEmitterDecorator");
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    ASSERT_EQ(CountPdbsIn(outdir), 1u);
    bool saw_decorator = false;
    for (const auto& e : fs::directory_iterator(outdir)) {
        if (e.path().filename().string().find("chain-test") != std::string::npos) {
            saw_decorator = true;
        }
    }
    EXPECT_TRUE(saw_decorator)
        << "expected decorator 'chain-test' in emitted PDB filename";
}


// ============================================================================
// Content sanity on a single emitted frame: HEADER, ATOM count matches
// protein.AtomCount, TER present, CONECT count matches disulfide count,
// END present, CRYST1 present (TRR carries box matrix on this fixture).
// ============================================================================

TEST_F(FramePdbEmitterTest, ContentSanity) {
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    fs::path outdir = FreshTempDir("content_sanity");

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();

    nmr::FramePdbEmitter::Config cfg;
    cfg.output_dir = outdir;
    cfg.stem       = "p1p9j";
    cfg.stride     = 99999;   // only frame 0
    nmr::FramePdbEmitter::Configure(tp.ProteinRef(), cfg);

    nmr::RunConfiguration config;
    ConfigureMinimalRun(config, "FramePdbEmitterContentSanity");
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    ASSERT_EQ(CountPdbsIn(outdir), 1u);
    fs::path pdb_path;
    for (const auto& e : fs::directory_iterator(outdir)) {
        pdb_path = e.path();
        break;
    }

    std::ifstream in(pdb_path);
    ASSERT_TRUE(in.is_open()) << "could not read " << pdb_path;
    std::stringstream ss; ss << in.rdbuf();
    const std::string body = ss.str();

    auto count_records = [&](const std::string& tag) {
        size_t n = 0; size_t pos = 0;
        while ((pos = body.find(tag, pos)) != std::string::npos) {
            // Match only at start of line: prev char is '\n' or pos == 0.
            if (pos == 0 || body[pos - 1] == '\n') ++n;
            pos += tag.size();
        }
        return n;
    };

    EXPECT_GE(count_records("HEADER"), 1u);
    EXPECT_GT(count_records("REMARK"), 0u);
    EXPECT_EQ(count_records("ATOM  "), tp.ProteinRef().AtomCount());
    EXPECT_GE(count_records("TER   "), 1u);
    EXPECT_GE(count_records("END"), 1u);

    // CRYST1: TRR carries box matrix on this fixture, so we expect one.
    EXPECT_GE(count_records("CRYST1"), 1u)
        << "TRR fixture should yield a CRYST1 record from box matrix";

    // CONECT: one per disulfide bond.
    size_t expected_conect = 0;
    for (const auto& b : tp.ProteinRef().LegacyAmber().BondList()) {
        if (b.category == nmr::BondCategory::Disulfide) ++expected_conect;
    }
    EXPECT_EQ(count_records("CONECT"), expected_conect)
        << "one CONECT per disulfide bond";
}
