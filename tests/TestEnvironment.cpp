#include "TestEnvironment.h"
#include "OperationLog.h"
#include "CalculatorConfig.h"

#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <cstdio>

namespace fs = std::filesystem;

namespace nmr {
namespace test {

std::string TestEnvironment::ubq_protonated_;
std::string TestEnvironment::ubq_crystal_;
std::string TestEnvironment::gmx_protonated_;
std::string TestEnvironment::orca_dir_;
std::string TestEnvironment::consolidated_;
std::string TestEnvironment::ff14sb_params_;
std::string TestEnvironment::aimnet2_model_;
std::string TestEnvironment::baseline_features_;
std::string TestEnvironment::fleet_amber_;
std::string TestEnvironment::fleet_amber_1p9j_5801_subpath_;
std::string TestEnvironment::fleet_amber_1z9b_6577_subpath_;
std::string TestEnvironment::larsen_1ubq_pm6_pdb_;
bool TestEnvironment::loaded_ = false;


bool TestEnvironment::RequireLoaded() {
    if (loaded_) return true;
    fprintf(stderr,
        "FATAL: TestEnvironment::Load() was not called. "
        "Call it in test_main.cpp before RUN_ALL_TESTS.\n");
    std::abort();
    return false;
}


void TestEnvironment::Load() {
    // Find testpaths.toml: NMR_TEST_DATA_DIR/../testpaths.toml
    // i.e., in the tests/ directory alongside the test sources.
    std::string toml_path;
#ifdef NMR_TEST_DATA_DIR
    toml_path = std::string(NMR_TEST_DATA_DIR) + "/../testpaths.toml";
#endif

    // Env var override
    const char* env = std::getenv("NMR_TESTPATHS_TOML");
    if (env) toml_path = env;

    // Repo root = parent of the tests/ dir that holds NMR_TEST_DATA_DIR.
    // Used to resolve repo-root-relative fixture paths (C10).
    std::string repo_root;
#ifdef NMR_TEST_DATA_DIR
    repo_root = fs::path(NMR_TEST_DATA_DIR).parent_path().parent_path().string();
#endif
    // A value that does not start with '/' is repo-root-relative; absolute
    // values (external data) are used verbatim. The *_subpath fragments are
    // deliberately relative and joined to fleet_amber at the call site, so
    // they are excluded here.
    auto resolveRelative = [&](const std::string& v) -> std::string {
        if (v.empty() || v.front() == '/' || repo_root.empty()) return v;
        return repo_root + "/" + v;
    };
    auto applyEnvPath = [&](std::string& slot, const char* env_name) {
        const char* env = std::getenv(env_name);
        if (env && *env) slot = resolveRelative(env);
    };
    auto ensureTrailingSlash = [](std::string& v) {
        if (!v.empty() && v.back() != '/') v.push_back('/');
    };

    if (!toml_path.empty() && fs::exists(toml_path)) {
        std::ifstream in(toml_path);
        std::string line;
        while (std::getline(in, line)) {
            auto pos = line.find('#');
            if (pos != std::string::npos) line = line.substr(0, pos);
            if (line.find('=') == std::string::npos) continue;

            auto eq = line.find('=');
            std::string key = line.substr(0, eq);
            std::string val = line.substr(eq + 1);

            auto trim = [](std::string& s) {
                while (!s.empty() && (s.front() == ' ' || s.front() == '\t' || s.front() == '"'))
                    s.erase(s.begin());
                while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '"'))
                    s.pop_back();
            };
            trim(key);
            trim(val);

            if (key.find("subpath") == std::string::npos)
                val = resolveRelative(val);

            if      (key == "ubq_protonated")  ubq_protonated_ = val;
            else if (key == "ubq_crystal")     ubq_crystal_ = val;
            else if (key == "gmx_protonated")  gmx_protonated_ = val;
            else if (key == "orca_dir")        orca_dir_ = val;
            else if (key == "consolidated")    consolidated_ = val;
            else if (key == "ff14sb_params")   ff14sb_params_ = val;
            else if (key == "aimnet2_model")   aimnet2_model_ = val;
            else if (key == "baseline_features") baseline_features_ = val;
            else if (key == "fleet_amber")     fleet_amber_ = val;
            else if (key == "fleet_amber_1p9j_5801_subpath")
                fleet_amber_1p9j_5801_subpath_ = val;
            else if (key == "fleet_amber_1z9b_6577_subpath")
                fleet_amber_1z9b_6577_subpath_ = val;
            else if (key == "larsen_1ubq_pm6_pdb")
                larsen_1ubq_pm6_pdb_ = val;
        }
        OperationLog::Info("TestEnvironment::Load", "read " + toml_path);
    } else {
        OperationLog::Error("TestEnvironment::Load",
            "testpaths.toml not found at " + toml_path +
            " — test data paths will be empty. Fix testpaths.toml.");
    }

    applyEnvPath(ubq_protonated_, "NMR_UBQ_PROTONATED");
    applyEnvPath(ubq_crystal_, "NMR_UBQ_CRYSTAL");
    applyEnvPath(gmx_protonated_, "NMR_GMX_PROTONATED");
    applyEnvPath(orca_dir_, "NMR_ORCA_TEST_DIR");
    applyEnvPath(consolidated_, "NMR_CONSOLIDATED_DIR");
    applyEnvPath(ff14sb_params_, "NMR_FF14SB_PARAMS");
    applyEnvPath(aimnet2_model_, "NMR_AIMNET2_MODEL");
    applyEnvPath(baseline_features_, "NMR_BASELINE_FEATURES");
    applyEnvPath(fleet_amber_, "NMR_FLEET_AMBER_DIR");
    applyEnvPath(larsen_1ubq_pm6_pdb_, "NMR_LARSEN_1UBQ_PM6_PDB");
    ensureTrailingSlash(orca_dir_);
    ensureTrailingSlash(consolidated_);

    loaded_ = true;

    auto status = [](const std::string& v) -> std::string {
        if (v.empty()) return "<not set>";
        return v;
    };

    OperationLog::Info("TestEnvironment::Load",
        "ubq=" + status(ubq_protonated_) +
        " ubq_crystal=" + status(ubq_crystal_) +
        " gmx=" + status(gmx_protonated_) +
        " orca=" + status(orca_dir_) +
        " consolidated=" + status(consolidated_) +
        " ff14sb=" + status(ff14sb_params_) +
        " aimnet2=" + status(aimnet2_model_) +
        " baseline=" + status(baseline_features_) +
        " fleet_amber=" + status(fleet_amber_) +
        " fleet_amber_1p9j_subpath=" + status(fleet_amber_1p9j_5801_subpath_) +
        " fleet_amber_1z9b_subpath=" + status(fleet_amber_1z9b_6577_subpath_) +
        " larsen_1ubq_pm6=" + status(larsen_1ubq_pm6_pdb_));
}


const std::string& TestEnvironment::UbqProtonated()   { RequireLoaded(); return ubq_protonated_; }
const std::string& TestEnvironment::UbqCrystal()     { RequireLoaded(); return ubq_crystal_; }
const std::string& TestEnvironment::GmxProtonated()   { RequireLoaded(); return gmx_protonated_; }
const std::string& TestEnvironment::OrcaDir()         { RequireLoaded(); return orca_dir_; }
const std::string& TestEnvironment::Consolidated()    { RequireLoaded(); return consolidated_; }
const std::string& TestEnvironment::Ff14sbParams()    { RequireLoaded(); return ff14sb_params_; }
const std::string& TestEnvironment::Aimnet2Model()    { RequireLoaded(); return aimnet2_model_; }
const std::string& TestEnvironment::BaselineFeatures() { RequireLoaded(); return baseline_features_; }
const std::string& TestEnvironment::FleetAmberData()  { RequireLoaded(); return fleet_amber_; }
const std::string& TestEnvironment::Larsen1UbqPm6Pdb() { RequireLoaded(); return larsen_1ubq_pm6_pdb_; }


std::string TestEnvironment::TempPath(const std::string& stem) {
    return (fs::temp_directory_path() / stem).string();
}


void TestEnvironment::LoadCalculatorConfig() {
    RequireLoaded();
    OperationLog::SetChannelMask(0xFFFFFFFF);
#ifdef NMR_TEST_DATA_DIR
    CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");
#else
    CalculatorConfig::Load();
#endif
}


AmberTrajectoryFixture TestEnvironment::FleetAmberTrajectory(
        const std::string& protein_id) {
    RequireLoaded();
    AmberTrajectoryFixture fix;
    fix.protein_id = protein_id;
    if (fleet_amber_.empty()) return fix;

    const std::string* subpath = nullptr;
    if (protein_id == "1P9J_5801")      subpath = &fleet_amber_1p9j_5801_subpath_;
    else if (protein_id == "1Z9B_6577") subpath = &fleet_amber_1z9b_6577_subpath_;
    if (!subpath || subpath->empty()) return fix;

    const std::string base = fleet_amber_ + "/" + *subpath;
    fix.tpr_path = base + ".tpr";
    fix.xtc_path = base + ".xtc";
    fix.edr_path = base + ".edr";
    return fix;
}

}  // namespace test
}  // namespace nmr
