/// @file
/// Behaviour tests for nmr::cli::Parse.

#include <gtest/gtest.h>

#include "Cli/Parse.h"
#include "Cli/Validate.h"
#include "Of3Loader.h"
#include "OrcaRunLoader.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

// Compile-time no-DFT contract for --of3. `Of3Input` (the neutral structural
// carrier) and `Of3Mode` must not carry an nmr_out_path member; `OrcaRunFiles`
// (the ORCA carrier) must. These static_asserts make the mode boundary a
// compile-time fact, not a runtime hope, without raising the project's C++17
// language level.
template <class T, class = void>
struct HasNmrOutPath : std::false_type {};

template <class T>
struct HasNmrOutPath<
    T, std::void_t<decltype(std::declval<T&>().nmr_out_path)>>
    : std::true_type {};

static_assert(!HasNmrOutPath<nmr::Of3Input>::value,
              "Of3Input must not carry a DFT/NMR path");
static_assert(!HasNmrOutPath<nmr::cli::Of3Mode>::value,
              "Of3Mode must not carry a DFT/NMR path");
static_assert(std::is_same_v<
                  decltype((std::declval<nmr::cli::Of3Mode&>().input)),
                  nmr::Of3Input&>,
              "Of3Mode.input must be the neutral Of3Input carrier");
// Contrast: the ORCA carrier DOES own the DFT path (documents the boundary).
static_assert(HasNmrOutPath<nmr::OrcaRunFiles>::value,
              "OrcaRunFiles owns the DFT path; --of3 must not");

// Recording overload set for the typed-dispatch wiring check (spec §6.3):
// visiting a parsed ModeSpec must select the alternative matching the flag.
struct RecordingVisitor {
    std::string operator()(const nmr::cli::PdbMode&)           const { return "pdb"; }
    std::string operator()(const nmr::cli::ProtonatedPdbMode&) const { return "protonated_pdb"; }
    std::string operator()(const nmr::cli::Of3Mode&)           const { return "of3"; }
    std::string operator()(const nmr::cli::OrcaMode&)          const { return "orca"; }
    std::string operator()(const nmr::cli::MutantMode&)        const { return "mutant"; }
    std::string operator()(const nmr::cli::TrajectoryMode&)    const { return "trajectory"; }
};

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


// ---- Trajectory window (--window-start N --window-len M; loud, never silent) ----

TEST(CliParse, TrajectoryNoWindowByDefault) {
    Argv a{"nmr_extract", "--trajectory", "/data/run"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    const auto& m = std::get<nmr::cli::TrajectoryMode>(*r.spec);
    EXPECT_EQ(m.window_start, 0u);
    EXPECT_EQ(m.window_len, 0u);  // 0 = no window = full trajectory
}

TEST(CliParse, TrajectoryWindowBothFlags) {
    Argv a{"nmr_extract", "--trajectory", "/data/run",
           "--window-start", "100", "--window-len", "20"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    const auto& m = std::get<nmr::cli::TrajectoryMode>(*r.spec);
    EXPECT_EQ(m.window_start, 100u);
    EXPECT_EQ(m.window_len, 20u);
}

TEST(CliParse, TrajectoryWindowStartAloneIsError) {
    // Both-or-neither: a lone --window-start must not silently run full.
    Argv a{"nmr_extract", "--trajectory", "/data/run", "--window-start", "100"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("both or neither"), std::string::npos);
}

TEST(CliParse, TrajectoryWindowLenMissingValueIsError) {
    // --window-len as the final token has no value; GetArg reports empty and
    // strtoull("")==0 would silently mean "no window". Must be a loud error.
    Argv a{"nmr_extract", "--trajectory", "/data/run",
           "--window-start", "100", "--window-len"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("require a non-negative integer"), std::string::npos);
}

TEST(CliParse, TrajectoryWindowLenZeroIsError) {
    // An explicit zero-length window is empty, not "full run" — reject loudly.
    Argv a{"nmr_extract", "--trajectory", "/data/run",
           "--window-start", "100", "--window-len", "0"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find(">= 1"), std::string::npos);
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


// ============================================================================
// --of3 mode: parse shape, flags, errors, typed dispatch, validation.
// ============================================================================

TEST(CliParse, Of3Root) {
    Argv a{"nmr_extract", "--of3", "--root", "/r/model"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    ASSERT_TRUE(std::holds_alternative<nmr::cli::Of3Mode>(*r.spec));
    const auto& m = std::get<nmr::cli::Of3Mode>(*r.spec);
    // Root expands to EXACTLY .prmtop + .inpcrd — never {root}.xyz, never
    // {root}_nmr.out. Geometry comes from the tleap rst7 restart.
    EXPECT_EQ(m.input.prmtop_path, "/r/model.prmtop");
    EXPECT_EQ(m.input.inpcrd_path, "/r/model.inpcrd");
    // Honest OpenFold provenance, supplied at the mode boundary.
    EXPECT_EQ(m.input.prediction_method, "OpenFold+tleap");
    // MOPAC defaults on for single-pose modes.
    EXPECT_TRUE(m.mopac);
}

TEST(CliParse, Of3NoMopac) {
    Argv a{"nmr_extract", "--of3", "--root", "/r/model", "--no-mopac"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    ASSERT_TRUE(std::holds_alternative<nmr::cli::Of3Mode>(*r.spec));
    // The shared single-conformation flag turns only `mopac` off.
    EXPECT_FALSE(std::get<nmr::cli::Of3Mode>(*r.spec).mopac);
}

TEST(CliParse, Of3WithoutRootIsError) {
    Argv a{"nmr_extract", "--of3"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("--of3 requires --root"), std::string::npos);
}

TEST(CliParse, RejectsOf3CombinedWithOrca) {
    // Two mode flags => the existing exclusive-mode error, unchanged.
    Argv a{"nmr_extract", "--of3", "--root", "/r/a", "--orca", "--root", "/r/b"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    EXPECT_FALSE(r.spec.has_value());
    EXPECT_NE(r.error.find("multiple mode flags"), std::string::npos);
}

// Regression: the ORCA parse shape is untouched by the OF3 additions — its
// carrier still owns all three paths including {root}_nmr.out.
TEST(CliParse, OrcaRootUnaffectedByOf3) {
    Argv a{"nmr_extract", "--orca", "--root", "/r/WT"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    ASSERT_TRUE(std::holds_alternative<nmr::cli::OrcaMode>(*r.spec));
    const auto& m = std::get<nmr::cli::OrcaMode>(*r.spec);
    EXPECT_EQ(m.files.xyz_path,     "/r/WT.xyz");
    EXPECT_EQ(m.files.prmtop_path,  "/r/WT.prmtop");
    EXPECT_EQ(m.files.nmr_out_path, "/r/WT_nmr.out");
}

// Typed-dispatch wiring (spec §6.3 check 1): a parsed --of3 spec selects the
// Of3Mode overload under std::visit.
TEST(CliDispatch, Of3ModeSelectedByVisitor) {
    Argv a{"nmr_extract", "--of3", "--root", "/r/model"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_EQ(std::visit(RecordingVisitor{}, *r.spec), "of3");
}

TEST(CliDispatch, OrcaModeStillSelectedByVisitor) {
    Argv a{"nmr_extract", "--orca", "--root", "/r/WT"};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    EXPECT_EQ(std::visit(RecordingVisitor{}, *r.spec), "orca");
}


// ---- Validation wiring (real temp files; never derives an _nmr.out) ----

namespace {

namespace fs = std::filesystem;

// Unique temp dir for a validation case; touch requested files under it.
fs::path MakeOf3Case(const std::string& tag, bool inpcrd, bool prmtop) {
    const fs::path dir = fs::path(testing::TempDir()) /
        ("of3_validate_" + tag + "_" + std::to_string(::getpid()));
    fs::create_directories(dir);
    const fs::path root = dir / "model";
    if (prmtop) std::ofstream(root.string() + ".prmtop") << "%VERSION\n";
    if (inpcrd) std::ofstream(root.string() + ".inpcrd") << "title\n0\n";
    return root;
}

}  // namespace

TEST(CliValidate, Of3PresentPasses) {
    const fs::path root = MakeOf3Case("ok", /*inpcrd=*/true, /*prmtop=*/true);
    Argv a{"nmr_extract", "--of3", "--root", root.string().c_str()};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    // No _nmr.out was created; OF3 validation must still pass (it does not
    // derive or inspect an NMR sibling). This catches accidental reuse of the
    // ORCA validator.
    EXPECT_FALSE(fs::exists(root.string() + "_nmr.out"));
    EXPECT_EQ(nmr::cli::Validate(*r.spec, nmr::cli::CommonOptions{}), "");
    // No fs::remove_all cleanup: torch's bundled libstdc++ ships a broken
    // std::filesystem::remove_all (segfaults). The temp dir is unique per
    // pid+tag under the ephemeral TempDir; leaving it is harmless.
}

TEST(CliValidate, Of3MissingInpcrdFails) {
    const fs::path root = MakeOf3Case("noinpcrd", /*inpcrd=*/false, /*prmtop=*/true);
    Argv a{"nmr_extract", "--of3", "--root", root.string().c_str()};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    const std::string e = nmr::cli::Validate(*r.spec, nmr::cli::CommonOptions{});
    EXPECT_NE(e.find("--of3 .inpcrd"), std::string::npos) << e;
    // No fs::remove_all cleanup (torch broken remove_all); see note above.
}

TEST(CliValidate, Of3MissingPrmtopFails) {
    const fs::path root = MakeOf3Case("noprm", /*inpcrd=*/true, /*prmtop=*/false);
    Argv a{"nmr_extract", "--of3", "--root", root.string().c_str()};
    const auto r = nmr::cli::Parse(a.argc(), a.argv());
    ASSERT_TRUE(r.spec.has_value());
    const std::string e = nmr::cli::Validate(*r.spec, nmr::cli::CommonOptions{});
    EXPECT_NE(e.find("--of3 .prmtop"), std::string::npos) << e;
    // No fs::remove_all cleanup (torch broken remove_all); see note above.
}
