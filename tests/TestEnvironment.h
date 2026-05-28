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
//   ff14sb_params                  = "path/to/ff14sb_params.dat"
//   aimnet2_model                  = "path/to/aimnet2_wb97m_0.jpt"
//   fleet_amber                    = "tests/data/fleet_amber"
//   fleet_amber_<id>_subpath       = pinned per-protein production sub-path
//   larsen_1ubq_pm6_pdb            = external Larsen-archive PM6-D3H+ 1UBQ
//
// Path resolution (C10): a value that does not begin with '/' is treated
// as repo-root-relative and resolved against the repo root (the parent of
// NMR_TEST_DATA_DIR's tests/ directory). External/local fixture roots are
// supplied per machine via testpaths.toml or env vars such as
// NMR_CONSOLIDATED_DIR, NMR_AIMNET2_MODEL, and NMR_LARSEN_1UBQ_PM6_PDB.
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
    static const std::string& Ff14sbParams();
    static const std::string& Aimnet2Model();
    static const std::string& BaselineFeatures();

    // AMBER-ff GROMACS trajectory fixtures. FleetAmberData() is the
    // root; FleetAmberTrajectory(protein_id) builds a typed fixture
    // with explicit paths from the pinned subpath in testpaths.toml.
    // Returns a fixture with empty paths if the protein_id is unknown
    // (caller should GTEST_SKIP).
    static const std::string& FleetAmberData();
    static AmberTrajectoryFixture FleetAmberTrajectory(
        const std::string& protein_id);

    // External Larsen-archive PM6-D3H+ optimised 1UBQ geometry. Empty if
    // the key is unset; callers GTEST_SKIP when the file is absent.
    static const std::string& Larsen1UbqPm6Pdb();

    // A path under the system temp directory, for test debug artifacts.
    // No discovery, no fixed /tmp literals — std::filesystem chooses TMPDIR.
    static std::string TempPath(const std::string& stem);

    // Per-test setup every trajectory test previously inlined: enable the
    // full log channel mask and load the canonical data/calculator_params.toml
    // (the same file test_main loads globally). Replaces the byte-identical
    // local LoadCalculatorConfig() helpers (G3); the prior local path was
    // tests/data/../data/... which did not exist, so those silently fell
    // back to compiled defaults.
    static void LoadCalculatorConfig();

    static bool RequireLoaded();

private:
    static std::string ubq_protonated_;
    static std::string ubq_crystal_;
    static std::string gmx_protonated_;
    static std::string orca_dir_;
    static std::string consolidated_;
    static std::string ff14sb_params_;
    static std::string aimnet2_model_;
    static std::string baseline_features_;
    static std::string fleet_amber_;
    static std::string fleet_amber_1p9j_5801_subpath_;
    static std::string fleet_amber_1z9b_6577_subpath_;
    static std::string larsen_1ubq_pm6_pdb_;
    static bool loaded_;
};

}  // namespace test
}  // namespace nmr
