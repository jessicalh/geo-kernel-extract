#pragma once
//
// TestEnvironment: where test data lives on this machine.
//
// Same pattern as RuntimeEnvironment: explicit Load(), RequireLoaded()
// on every accessor, complete state logged. Reads testpaths.toml from
// the project tests/ directory (path set by NMR_TEST_DATA_DIR CMake macro).
//
// testpaths.toml keys:
//   ubq_protonated                 = "tests/data/1ubq_protonated.pdb"
//   orca_dir                       = "path/to/orca/test/pair"
//   consolidated                   = "path/to/734/mutant/pairs"
//   fleet_data                     = "tests/data/fleet"
//   ff14sb_params                  = "path/to/ff14sb_params.dat"
//   fleet_amber                    = "tests/data/fleet_amber"
//   fleet_amber_<id>_subpath       = pinned per-protein production sub-path
//

#include <string>

namespace nmr {
namespace test {

// Trajectory file paths pinned for test fixtures. No discovery; the
// testpaths.toml entries name the production sub-path under
// <fleet_amber>/<protein_id>/ explicitly.
struct AmberTrajectoryFixture {
    std::string protein_id;
    std::string tpr_path;
    std::string xtc_path;
    std::string edr_path;
};

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

    // AMBER-ff GROMACS trajectory fixtures. FleetAmberData() is the
    // root; FleetAmberTrajectory(protein_id) builds a typed fixture
    // with explicit paths from the pinned subpath in testpaths.toml.
    // Returns a fixture with empty paths if the protein_id is unknown
    // (caller should GTEST_SKIP).
    static const std::string& FleetAmberData();
    static AmberTrajectoryFixture FleetAmberTrajectory(
        const std::string& protein_id);

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
    static std::string fleet_amber_;
    static std::string fleet_amber_1p9j_5801_subpath_;
    static std::string fleet_amber_1z9b_6577_subpath_;
    static bool loaded_;
};

}  // namespace test
}  // namespace nmr
