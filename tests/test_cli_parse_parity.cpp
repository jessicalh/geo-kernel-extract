/// @file
/// Behaviour-equivalence tests: nmr::cli::Parse vs the legacy ParseJobSpec.
///
/// For every documented CLI surface, the two parsers must agree on:
///   - the mode identified
///   - the input files / dir
///   - the boolean toggles (mopac, apbs, coulomb)
///   - the common options (output, config, aimnet2)
///   - the trajectory emission settings (PDB stride, NPY stride, etc.)
///
/// Differences that are intentional (new parser is stricter on
/// multiple-mode-flag input; new trajectory mode does not parse
/// --no-apbs / --no-coulomb because the trajectory pipeline ignores
/// them) are tested separately under the @c CliParseStrictness suite.

#include <gtest/gtest.h>

#include "Cli/Parse.h"
#include "JobSpec.h"

#include <cstring>
#include <string>
#include <variant>
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

void ExpectParity(Argv& a) {
    const auto old_spec = nmr::ParseJobSpec(a.argc(), a.argv());
    const auto new_res  = nmr::cli::Parse(a.argc(), a.argv());

    ASSERT_TRUE(new_res.spec.has_value())
        << "new parser rejected: " << new_res.error;
    ASSERT_TRUE(new_res.common.has_value());
    ASSERT_NE(old_spec.mode, nmr::JobMode::None) << "old parser rejected";

    // Common options.
    EXPECT_EQ(new_res.common->output_dir.string(),         old_spec.output_dir);
    EXPECT_EQ(new_res.common->config_path.string(),        old_spec.config_path);
    EXPECT_EQ(new_res.common->aimnet2_model_path.string(), old_spec.aimnet2_model_path);

    std::visit(
        [&](const auto& m) {
            using T = std::decay_t<decltype(m)>;
            if constexpr (std::is_same_v<T, nmr::cli::PdbMode>) {
                EXPECT_EQ(old_spec.mode, nmr::JobMode::Pdb);
                EXPECT_EQ(m.pdb.string(), old_spec.pdb_path);
                EXPECT_DOUBLE_EQ(m.pH, old_spec.pH);
                EXPECT_EQ(m.mopac,   !old_spec.skip_mopac);
                EXPECT_EQ(m.apbs,    !old_spec.skip_apbs);
                EXPECT_EQ(m.coulomb, !old_spec.skip_coulomb);
            } else if constexpr (std::is_same_v<T, nmr::cli::ProtonatedPdbMode>) {
                EXPECT_EQ(old_spec.mode, nmr::JobMode::ProtonatedPdb);
                EXPECT_EQ(m.pdb.string(), old_spec.pdb_path);
                EXPECT_EQ(m.mopac,   !old_spec.skip_mopac);
                EXPECT_EQ(m.apbs,    !old_spec.skip_apbs);
                EXPECT_EQ(m.coulomb, !old_spec.skip_coulomb);
            } else if constexpr (std::is_same_v<T, nmr::cli::OrcaMode>) {
                EXPECT_EQ(old_spec.mode, nmr::JobMode::Orca);
                EXPECT_EQ(m.files.xyz_path,     old_spec.orca_files.xyz_path);
                EXPECT_EQ(m.files.prmtop_path,  old_spec.orca_files.prmtop_path);
                EXPECT_EQ(m.files.nmr_out_path, old_spec.orca_files.nmr_out_path);
                EXPECT_EQ(m.mopac,   !old_spec.skip_mopac);
            } else if constexpr (std::is_same_v<T, nmr::cli::MutantMode>) {
                EXPECT_EQ(old_spec.mode, nmr::JobMode::Mutant);
                EXPECT_EQ(m.wt.xyz_path,     old_spec.wt_files.xyz_path);
                EXPECT_EQ(m.wt.prmtop_path,  old_spec.wt_files.prmtop_path);
                EXPECT_EQ(m.ala.xyz_path,    old_spec.ala_files.xyz_path);
                EXPECT_EQ(m.ala.prmtop_path, old_spec.ala_files.prmtop_path);
                EXPECT_EQ(m.mopac, !old_spec.skip_mopac);
            } else if constexpr (std::is_same_v<T, nmr::cli::TrajectoryMode>) {
                EXPECT_EQ(old_spec.mode, nmr::JobMode::Trajectory);
                EXPECT_EQ(m.dir.string(), old_spec.traj_dir);
                EXPECT_EQ(m.mopac, old_spec.enable_mopac);
                if (m.emit_pdbs.has_value()) {
                    EXPECT_EQ(m.emit_pdbs->output_dir.string(), old_spec.emit_frame_pdbs_dir);
                    EXPECT_EQ(m.emit_pdbs->stride,              old_spec.pdb_stride);
                    EXPECT_EQ(m.emit_pdbs->from_ps,             old_spec.pdb_from_ps);
                    EXPECT_EQ(m.emit_pdbs->to_ps,               old_spec.pdb_to_ps);
                    EXPECT_EQ(m.emit_pdbs->decorator,           old_spec.pdb_decorator);
                } else {
                    EXPECT_TRUE(old_spec.emit_frame_pdbs_dir.empty());
                }
                if (m.emit_npys.has_value()) {
                    EXPECT_EQ(m.emit_npys->output_dir.string(), old_spec.emit_frame_npys_dir);
                    EXPECT_EQ(m.emit_npys->stride,              old_spec.npy_stride);
                    EXPECT_EQ(m.emit_npys->from_ps,             old_spec.npy_from_ps);
                    EXPECT_EQ(m.emit_npys->to_ps,               old_spec.npy_to_ps);
                } else {
                    EXPECT_TRUE(old_spec.emit_frame_npys_dir.empty());
                }
            }
        },
        *new_res.spec);
}

}  // namespace


// ---- Pdb mode ----

TEST(CliParseParity, PdbBare)            { Argv a{"nmr_extract", "--pdb", "x.pdb"}; ExpectParity(a); }
TEST(CliParseParity, PdbWithpH)          { Argv a{"nmr_extract", "--pdb", "x.pdb", "--pH", "6.5"}; ExpectParity(a); }
TEST(CliParseParity, PdbNoMopac)         { Argv a{"nmr_extract", "--pdb", "x.pdb", "--no-mopac"}; ExpectParity(a); }
TEST(CliParseParity, PdbNoApbs)          { Argv a{"nmr_extract", "--pdb", "x.pdb", "--no-apbs"}; ExpectParity(a); }
TEST(CliParseParity, PdbNoCoulomb)       { Argv a{"nmr_extract", "--pdb", "x.pdb", "--no-coulomb"}; ExpectParity(a); }
TEST(CliParseParity, PdbAllCommonFlags)  {
    Argv a{"nmr_extract", "--pdb", "x.pdb",
           "--output", "/tmp/out",
           "--config", "params.toml",
           "--aimnet2", "/m/aimnet.jpt"};
    ExpectParity(a);
}

// ---- ProtonatedPdb mode ----

TEST(CliParseParity, ProtonatedPdbBare)    { Argv a{"nmr_extract", "--protonated-pdb", "p.pdb"}; ExpectParity(a); }
TEST(CliParseParity, ProtonatedPdbNoMopac) { Argv a{"nmr_extract", "--protonated-pdb", "p.pdb", "--no-mopac"}; ExpectParity(a); }

// ---- Orca mode ----

TEST(CliParseParity, OrcaBare)    { Argv a{"nmr_extract", "--orca", "--root", "/r/WT"}; ExpectParity(a); }
TEST(CliParseParity, OrcaNoMopac) { Argv a{"nmr_extract", "--orca", "--root", "/r/WT", "--no-mopac"}; ExpectParity(a); }
TEST(CliParseParity, OrcaAllFlags) {
    Argv a{"nmr_extract", "--orca", "--root", "/r/WT",
           "--no-apbs", "--no-coulomb",
           "--output", "/tmp/out",
           "--aimnet2", "/m/aimnet.jpt"};
    ExpectParity(a);
}

// ---- Mutant mode ----

TEST(CliParseParity, MutantBare) {
    Argv a{"nmr_extract", "--mutant", "--wt", "/r/WT", "--ala", "/r/ALA"};
    ExpectParity(a);
}
TEST(CliParseParity, MutantNoMopac) {
    Argv a{"nmr_extract", "--mutant", "--wt", "/r/WT", "--ala", "/r/ALA", "--no-mopac"};
    ExpectParity(a);
}

// ---- Trajectory mode ----

TEST(CliParseParity, TrajectoryBare) {
    Argv a{"nmr_extract", "--trajectory", "/data/run"};
    ExpectParity(a);
}
TEST(CliParseParity, TrajectoryMopac) {
    Argv a{"nmr_extract", "--trajectory", "/data/run", "--mopac"};
    ExpectParity(a);
}
TEST(CliParseParity, TrajectoryPdbEmission) {
    Argv a{"nmr_extract", "--trajectory", "/data/run",
           "--emit-frame-pdbs", "/tmp/pdbs",
           "--pdb-stride", "100",
           "--pdb-from-ps", "1000",
           "--pdb-to-ps", "5000",
           "--pdb-decorator", "run1"};
    ExpectParity(a);
}
TEST(CliParseParity, TrajectoryNpyEmission) {
    Argv a{"nmr_extract", "--trajectory", "/data/run",
           "--emit-frame-npys", "/tmp/npys",
           "--npy-stride", "500",
           "--npy-from-ps", "0",
           "--npy-to-ps", "10000"};
    ExpectParity(a);
}
TEST(CliParseParity, TrajectoryBothEmissions) {
    Argv a{"nmr_extract", "--trajectory", "/data/run",
           "--emit-frame-pdbs", "/tmp/pdbs", "--pdb-stride", "100",
           "--emit-frame-npys", "/tmp/npys", "--npy-stride", "500",
           "--aimnet2", "/m/aimnet.jpt",
           "--output", "/tmp/out"};
    ExpectParity(a);
}
TEST(CliParseParity, TrajectoryAfterCommonArgs) {
    // --trajectory positional must not be confused by an intervening flag's value.
    Argv a{"nmr_extract",
           "--aimnet2", "/m/aimnet.jpt",
           "--trajectory", "/data/run",
           "--output", "/tmp/out"};
    ExpectParity(a);
}


// ---- New parser stricter (no parity expected; documents the divergence) ----

TEST(CliParseStrictness, RejectsMultipleModeFlags) {
    Argv a{"nmr_extract", "--pdb", "x.pdb", "--protonated-pdb", "y.pdb"};
    const auto res = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(res.spec.has_value());
    EXPECT_NE(res.error.find("multiple mode flags"), std::string::npos);
}

TEST(CliParseStrictness, HelpRequestedWithDashH) {
    Argv a{"nmr_extract", "-h"};
    const auto res = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_TRUE(res.help_requested);
    EXPECT_FALSE(res.spec.has_value());
}

TEST(CliParseStrictness, HelpRequestedWithNoArgs) {
    Argv a{"nmr_extract"};
    const auto res = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_TRUE(res.help_requested);
}

TEST(CliParseStrictness, NoModeFlagIsError) {
    Argv a{"nmr_extract", "--output", "/tmp/out"};
    const auto res = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(res.spec.has_value());
    EXPECT_NE(res.error.find("no mode flag"), std::string::npos);
}

TEST(CliParseStrictness, PdbWithoutPathIsError) {
    Argv a{"nmr_extract", "--pdb"};
    const auto res = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(res.spec.has_value());
    EXPECT_NE(res.error.find("--pdb requires"), std::string::npos);
}

TEST(CliParseStrictness, OrcaWithoutRootIsError) {
    Argv a{"nmr_extract", "--orca"};
    const auto res = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(res.spec.has_value());
    EXPECT_NE(res.error.find("--orca requires --root"), std::string::npos);
}

TEST(CliParseStrictness, MutantNeedsBothRoots) {
    Argv a{"nmr_extract", "--mutant", "--wt", "/r/WT"};
    const auto res = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(res.spec.has_value());
    EXPECT_NE(res.error.find("--mutant requires"), std::string::npos);
}

TEST(CliParseStrictness, TrajectoryNeedsDir) {
    Argv a{"nmr_extract", "--trajectory"};
    const auto res = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(res.spec.has_value());
    EXPECT_NE(res.error.find("--trajectory requires"), std::string::npos);
}
