#pragma once
//
// TestEnvironment: where test data lives on this machine.
//
// Same pattern as RuntimeEnvironment: explicit Load(), RequireLoaded()
// on every accessor, complete state logged. Reads testpaths.toml from
// the project tests/ directory (path set by NMR_TEST_DATA_DIR CMake macro).
//
// testpaths.toml keys:
//   ubq_protonated  = "tests/data/1ubq_protonated.pdb"
//   orca_dir        = "path/to/orca/test/pair"
//   consolidated    = "path/to/734/mutant/pairs"
//   fleet_data      = "tests/data/fleet"
//   ff14sb_params   = "path/to/ff14sb_params.dat"  (test override)
//

#include <string>

namespace nmr {
namespace test {

class TestEnvironment {
public:
    static void Load();

    static const std::string& UbqProtonated();
    static const std::string& UbqCrystal();
    static const std::string& GmxProtonated();
    static const std::string& OrcaDir();
    static const std::string& Consolidated();
    static const std::string& FleetData();
    static const std::string& Ff14sbParams();
    static const std::string& BaselineFeatures();

    static bool RequireLoaded();

private:
    static std::string ubq_protonated_;
    static std::string ubq_crystal_;
    static std::string gmx_protonated_;
    static std::string orca_dir_;
    static std::string consolidated_;
    static std::string fleet_data_;
    static std::string ff14sb_params_;
    static std::string baseline_features_;
    static bool loaded_;
};

}  // namespace test
}  // namespace nmr
