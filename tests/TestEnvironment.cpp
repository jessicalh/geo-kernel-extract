#include "TestEnvironment.h"
#include "OperationLog.h"

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
std::string TestEnvironment::fleet_data_;
std::string TestEnvironment::ff14sb_params_;
std::string TestEnvironment::baseline_features_;
std::string TestEnvironment::fleet_amber_;
std::string TestEnvironment::fleet_amber_1p9j_5801_subpath_;
std::string TestEnvironment::fleet_amber_1z9b_6577_subpath_;
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

            if      (key == "ubq_protonated")  ubq_protonated_ = val;
            else if (key == "ubq_crystal")     ubq_crystal_ = val;
            else if (key == "gmx_protonated")  gmx_protonated_ = val;
            else if (key == "orca_dir")        orca_dir_ = val;
            else if (key == "consolidated")    consolidated_ = val;
            else if (key == "fleet_data")      fleet_data_ = val;
            else if (key == "ff14sb_params")   ff14sb_params_ = val;
            else if (key == "baseline_features") baseline_features_ = val;
            else if (key == "fleet_amber")     fleet_amber_ = val;
            else if (key == "fleet_amber_1p9j_5801_subpath")
                fleet_amber_1p9j_5801_subpath_ = val;
            else if (key == "fleet_amber_1z9b_6577_subpath")
                fleet_amber_1z9b_6577_subpath_ = val;
        }
        OperationLog::Info("TestEnvironment::Load", "read " + toml_path);
    } else {
        OperationLog::Error("TestEnvironment::Load",
            "testpaths.toml not found at " + toml_path +
            " — test data paths will be empty. Fix testpaths.toml.");
    }

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
        " fleet=" + status(fleet_data_) +
        " ff14sb=" + status(ff14sb_params_) +
        " baseline=" + status(baseline_features_) +
        " fleet_amber=" + status(fleet_amber_) +
        " fleet_amber_1p9j_subpath=" + status(fleet_amber_1p9j_5801_subpath_) +
        " fleet_amber_1z9b_subpath=" + status(fleet_amber_1z9b_6577_subpath_));
}


const std::string& TestEnvironment::UbqProtonated()   { RequireLoaded(); return ubq_protonated_; }
const std::string& TestEnvironment::UbqCrystal()     { RequireLoaded(); return ubq_crystal_; }
const std::string& TestEnvironment::GmxProtonated()   { RequireLoaded(); return gmx_protonated_; }
const std::string& TestEnvironment::OrcaDir()         { RequireLoaded(); return orca_dir_; }
const std::string& TestEnvironment::Consolidated()    { RequireLoaded(); return consolidated_; }
const std::string& TestEnvironment::FleetData()       { RequireLoaded(); return fleet_data_; }
const std::string& TestEnvironment::Ff14sbParams()    { RequireLoaded(); return ff14sb_params_; }
const std::string& TestEnvironment::BaselineFeatures() { RequireLoaded(); return baseline_features_; }
const std::string& TestEnvironment::FleetAmberData()  { RequireLoaded(); return fleet_amber_; }


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
