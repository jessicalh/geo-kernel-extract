#pragma once
//
// RuntimeEnvironment: the running program's configuration.
//
// PRECONDITION: Call RuntimeEnvironment::Load() once at program startup.
// Every accessor checks this and aborts if you forgot. There is
// no silent fallback. If the TOML doesn't exist, Load() still
// runs — it reads env vars and PATH, logs what it found, and
// records missing values as <not set>. But it must be called.
//
// TOML location: ~/.nmr_tools.toml (or pass explicit path).
// Every resolved value is logged at Load() time so you can see
// exactly what the program thinks its environment is.
//
// Live surface:
//   mopac         — MOPAC binary path for MopacResult subprocess calls
//   ff14sb_params — BuildFromPdb charge assignment
//   tleap          — AmberPreparedChargeSource runtime topology generation
//                    (resolved from TOML / AMBERHOME / PATH)
//   tmpdir         — temp files for MOPAC and AmberPreparedChargeSource
//                    work directories
//   bmrb_atom_nom  — data/bmrb_atom_nom.tbl path; consumed by
//                    CategoryInfoProjection at session start. Empty =
//                    no atom_nom.tbl-driven names; projection emits AMBER
//                    names as lookup fallbacks.
//   tensorcs15_dsn — libpq connection string for the local tensorcs15
//                    Postgres replica (ProCS15 tripeptide DFT data).
//                    Read from [databases].tensorcs15 in the TOML.
//                    Consumed by Session::LoadTripeptideDftTable.
//                    Empty = Session leaves the table unloaded, so
//                    OperationRunner skips tripeptide DFT calculators.
//   larsen_hbond_grids — directory holding the 6 dense.h5 grids
//                    produced by scripts/larsen_hbond_grid_parse/.
//                    Read from top-level `larsen_hbond_grids` key in
//                    the TOML. Consumed by Session::LoadLarsenHBondGrid.
//                    Empty = Session leaves the grid unloaded, so
//                    OperationRunner skips LarsenHBondShieldingResult.
//
// MOPAC and tleap are binary paths used for subprocess invocation.
//

#include <string>
#include <vector>

namespace nmr {

class RuntimeEnvironment {
public:
    // Load environment configuration. MUST be called before anything
    // else in this library. Reads TOML, falls back to env vars and
    // PATH. Logs the complete resolved state.
    // Empty path = use ~/.nmr_tools.toml
    static void Load(const std::string& tomlPath = "");

    // --- Live accessors ---

    static const std::string& Mopac();
    static const std::string& Tleap();
    static const std::string& Ff14sbParams();
    static const std::string& TmpDir();
    static const std::string& BmrbAtomNom();
    static const std::string& TensorCs15Dsn();
    static const std::string& LarsenHBondGridDir();

    // Verify required runtime files exist. Returns list of missing
    // (name + path) entries for MOPAC and ff14sb_params.
    static std::vector<std::string> Verify();

    // Generate a temp file path: tmpdir/guid_proteinName_suffix
    static std::string TempFilePath(const std::string& proteinName,
                                     const std::string& suffix);

    // Check that Load() has been called. Logs and aborts if not.
    // Goes in every accessor — same slot InitDefaults() occupied.
    // Future implementations copy this pattern and get the check.
    // PUBLIC so calculator implementations can use it too.
    static bool RequireLoaded();

private:
    static std::string mopac_;
    static std::string tleap_;
    static std::string ff14sb_params_;
    static std::string tmpDir_;
    static std::string bmrb_atom_nom_;
    static std::string tensorcs15_dsn_;
    static std::string larsen_hbond_grid_dir_;
    static std::string processGuid_;
    static bool loaded_;
};

}  // namespace nmr
