/// @file
/// Behaviour tests for nmr::cli::Parse.

#include <gtest/gtest.h>

#include "Cli/Parse.h"

#include <string>
#include <vector>

namespace {

/// Build an argv from a brace-list of literals.
class Argv {
public:
    Argv(std::initializer_list<const char*> args) {
        storage_.reserve(args.size());
        for (const char* a : args) storage_.emplace_back(a);
        ptrs_.reserve(storage_.size() + 1);
        for (auto& s : storage_) ptrs_.push_back(s.data());
        ptrs_.push_back(nullptr);
    }

    int    argc()       const { return static_cast<int>(storage_.size()); }
    char** argv()             { return ptrs_.data(); }

private:
    std::vector<std::string> storage_;
    std::vector<char*>       ptrs_;
};

}  // namespace


// ---- Successful parses ----

TEST(CliParse, PdbBare) {
    Argv a{"nmr_extract", "--pdb", "x.pdb"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    ASSERT_TRUE(std::holds_alternative<nmr::cli::PdbMode>(*r.spec));
    EXPECT_EQ(std::get<nmr::cli::PdbMode>(*r.spec).pdb.string(), "x.pdb");
}

TEST(CliParse, ProtonatedPdbBare) {
    Argv a{"nmr_extract", "--protonated-pdb", "p.pdb"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_TRUE(std::holds_alternative<nmr::cli::ProtonatedPdbMode>(*r.spec));
}

TEST(CliParse, OrcaRoot) {
    Argv a{"nmr_extract", "--orca", "--root", "/r/WT"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    ASSERT_TRUE(std::holds_alternative<nmr::cli::OrcaMode>(*r.spec));
    const auto& m = std::get<nmr::cli::OrcaMode>(*r.spec);
    EXPECT_EQ(m.files.xyz_path,    "/r/WT.xyz");
    EXPECT_EQ(m.files.prmtop_path, "/r/WT.prmtop");
    EXPECT_EQ(m.files.nmr_out_path, "/r/WT_nmr.out");
}

TEST(CliParse, MutantPair) {
    Argv a{"nmr_extract", "--mutant", "--wt", "/r/WT", "--ala", "/r/ALA"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    ASSERT_TRUE(std::holds_alternative<nmr::cli::MutantMode>(*r.spec));
    const auto& m = std::get<nmr::cli::MutantMode>(*r.spec);
    EXPECT_EQ(m.wt.xyz_path,  "/r/WT.xyz");
    EXPECT_EQ(m.ala.xyz_path, "/r/ALA.xyz");
}

TEST(CliParse, TrajectoryDir) {
    Argv a{"nmr_extract", "--trajectory", "/data/run"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    ASSERT_TRUE(std::holds_alternative<nmr::cli::TrajectoryMode>(*r.spec));
    EXPECT_EQ(std::get<nmr::cli::TrajectoryMode>(*r.spec).dir.string(), "/data/run");
}

TEST(CliParse, TrajectoryPositionalAfterFlag) {
    // --trajectory DIR must not be confused by an intervening flag's value.
    Argv a{"nmr_extract",
           "--aimnet2", "/m/aimnet.jpt",
           "--trajectory", "/data/run",
           "--output", "/tmp/out"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_EQ(std::get<nmr::cli::TrajectoryMode>(*r.spec).dir.string(), "/data/run");
}


// ---- Per-mode flags ----

TEST(CliParse, PdbWithpH) {
    Argv a{"nmr_extract", "--pdb", "x.pdb", "--pH", "6.5"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_DOUBLE_EQ(std::get<nmr::cli::PdbMode>(*r.spec).pH, 6.5);
}

TEST(CliParse, PdbNoMopacFlipsDefault) {
    Argv a{"nmr_extract", "--pdb", "x.pdb", "--no-mopac"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_FALSE(std::get<nmr::cli::PdbMode>(*r.spec).mopac);
}

TEST(CliParse, SingleConfMopacDefaultsOn) {
    // MOPAC is the one single-conf toggle; APBS is always on and
    // home-rolled Coulomb is retired, so neither is parseable.
    Argv a{"nmr_extract", "--pdb", "x.pdb"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_TRUE(std::get<nmr::cli::PdbMode>(*r.spec).mopac);
}

TEST(CliParse, TrajectoryMopacDefaultsOff) {
    Argv a{"nmr_extract", "--trajectory", "/data/run"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_FALSE(std::get<nmr::cli::TrajectoryMode>(*r.spec).mopac);
}

TEST(CliParse, TrajectoryMopacFlagOptsIn) {
    Argv a{"nmr_extract", "--trajectory", "/data/run", "--mopac"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_TRUE(std::get<nmr::cli::TrajectoryMode>(*r.spec).mopac);
}


// ---- Trajectory stride (the single cadence knob) ----

TEST(CliParse, TrajectoryStrideDefaultsToOne) {
    Argv a{"nmr_extract", "--trajectory", "/data/run"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_EQ(std::get<nmr::cli::TrajectoryMode>(*r.spec).stride, 1u);
}

TEST(CliParse, TrajectoryStrideFlag) {
    Argv a{"nmr_extract", "--trajectory", "/data/run", "--stride", "5"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_EQ(std::get<nmr::cli::TrajectoryMode>(*r.spec).stride, 5u);
}

TEST(CliParse, TrajectoryStrideZeroNormalisedToOne) {
    Argv a{"nmr_extract", "--trajectory", "/data/run", "--stride", "0"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_EQ(std::get<nmr::cli::TrajectoryMode>(*r.spec).stride, 1u);
}


// ---- Common options ----

TEST(CliParse, CommonOptionsCarryThrough) {
    Argv a{"nmr_extract", "--pdb", "x.pdb",
           "--output", "/tmp/out",
           "--config", "params.toml",
           "--aimnet2", "/m/aimnet.jpt"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.common.has_value());
    EXPECT_EQ(r.common->output_dir.string(),         "/tmp/out");
    EXPECT_EQ(r.common->config_path.string(),        "params.toml");
    EXPECT_EQ(r.common->aimnet2_model_path.string(), "/m/aimnet.jpt");
}


// ---- Error / help surface ----

TEST(CliParse, RejectsMultipleModeFlags) {
    Argv a{"nmr_extract", "--pdb", "x.pdb", "--protonated-pdb", "y.pdb"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("multiple mode flags"), std::string::npos);
}

TEST(CliParse, HelpRequestedWithDashH) {
    Argv a{"nmr_extract", "-h"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_TRUE(r.help_requested);
    EXPECT_FALSE(r.spec.has_value());
}

TEST(CliParse, HelpRequestedWithNoArgs) {
    Argv a{"nmr_extract"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_TRUE(r.help_requested);
}

TEST(CliParse, NoModeFlagIsError) {
    Argv a{"nmr_extract", "--output", "/tmp/out"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("no mode flag"), std::string::npos);
}

TEST(CliParse, PdbWithoutPathIsError) {
    Argv a{"nmr_extract", "--pdb"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("--pdb requires"), std::string::npos);
}

TEST(CliParse, OrcaWithoutRootIsError) {
    Argv a{"nmr_extract", "--orca"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("--orca requires --root"), std::string::npos);
}

TEST(CliParse, MutantNeedsBothRoots) {
    Argv a{"nmr_extract", "--mutant", "--wt", "/r/WT"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("--mutant requires"), std::string::npos);
}

TEST(CliParse, TrajectoryNeedsDir) {
    Argv a{"nmr_extract", "--trajectory"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("--trajectory requires"), std::string::npos);
}
