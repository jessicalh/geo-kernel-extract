#pragma once
//
// RuntimeEnvironment: the running program's configuration.
//
// PRECONDITION: Call RuntimeEnvironment::Load() once at program startup.
// Every accessor checks this and aborts if you forgot. There is
// no silent fallback. If the TOML doesn't exist, Load() still
// runs — it reads env vars and PATH, logs what it found, and
// warns about what's missing. But it must be called.
//
// TOML location: ~/.nmr_tools.toml (or pass explicit path).
// Every resolved value is logged at Load() time so you can see
// exactly what the program thinks its environment is.
//
// Live surface (2026-04-05):
//   ff14sb_params — BuildFromPdb charge assignment
//   tmpdir        — temp files for MOPAC work directories
//
// MOPAC is linked (libmopac.so), not a binary path. No accessor needed.
//
// Removed:
//   xtb           — replaced by linked libmopac (MopacResult)
//   propka, pdb2pqr, tleap, kaml, gmx, amberhome, apbs (binary)
//   These tools exist but are not called from any live code path.
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
    static const std::string& Ff14sbParams();
    static const std::string& TmpDir();

    // Verify live tools exist. Returns list of missing (name + path).
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
    static std::string ff14sb_params_;
    static std::string tmpDir_;
    static std::string processGuid_;
    static bool loaded_;
};

}  // namespace nmr
