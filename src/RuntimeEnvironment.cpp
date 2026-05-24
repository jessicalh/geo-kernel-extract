#include "RuntimeEnvironment.h"
#include "OperationLog.h"

#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <sstream>
#include <iomanip>

namespace fs = std::filesystem;

namespace nmr {

// Static members
std::string RuntimeEnvironment::mopac_;
std::string RuntimeEnvironment::tleap_;
std::string RuntimeEnvironment::ff14sb_params_;
std::string RuntimeEnvironment::tmpDir_;
std::string RuntimeEnvironment::bmrb_atom_nom_;
std::string RuntimeEnvironment::tensorcs15_dsn_;
std::string RuntimeEnvironment::larsen_hbond_grid_dir_;
std::string RuntimeEnvironment::processGuid_;
bool RuntimeEnvironment::loaded_ = false;


// ============================================================================
// RequireLoaded: the precondition check.
//
// Goes in every accessor, same slot the old InitDefaults() occupied.
// Returns true if loaded. Logs and aborts if not. Future code that
// copies the accessor pattern gets this check for free.
// ============================================================================

bool RuntimeEnvironment::RequireLoaded() {
    if (loaded_) return true;
    OperationLog::Error("RuntimeEnvironment",
        "FATAL: RuntimeEnvironment::Load() was not called before use. "
        "Call Load() at program startup before any library operations.");
    fprintf(stderr,
        "FATAL: RuntimeEnvironment::Load() was not called. "
        "Call it at program startup.\n");
    std::abort();
    return false;  // unreachable, but satisfies compilers
}


// Generate a short hex GUID for this process.
static std::string MakeGuid() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(0, 0xFFFFFFFF);
    std::ostringstream oss;
    oss << std::hex << std::setfill('0') << std::setw(8) << dist(gen);
    return oss.str();
}


// Resolve a binary: check TOML value first, then PATH.
// Returns empty string if not found anywhere — caller decides severity.
static std::string ResolveBinary(const std::string& toml_value,
                                  const std::string& bare_name) {
    if (!toml_value.empty() && fs::exists(toml_value))
        return toml_value;

    std::string which_cmd = "which " + bare_name + " 2>/dev/null";
    FILE* pipe = popen(which_cmd.c_str(), "r");
    if (pipe) {
        char buf[512];
        std::string result;
        while (fgets(buf, sizeof(buf), pipe))
            result += buf;
        pclose(pipe);
        while (!result.empty() && (result.back() == '\n' || result.back() == '\r'))
            result.pop_back();
        if (!result.empty() && fs::exists(result))
            return result;
    }

    return {};
}


// ============================================================================
// Load: read TOML, resolve everything, log the complete state.
// ============================================================================

void RuntimeEnvironment::Load(const std::string& tomlPath) {
    processGuid_ = MakeGuid();

    // --- Read TOML ---

    std::string path = tomlPath;
    if (path.empty()) {
        const char* home = std::getenv("HOME");
        if (home) path = std::string(home) + "/.nmr_tools.toml";
    }

    std::string toml_mopac, toml_tleap, toml_ff14sb, toml_tmpdir,
                toml_bmrb_atom_nom, toml_tensorcs15_dsn,
                toml_larsen_hbond_grid_dir;

    if (!path.empty() && fs::exists(path)) {
        std::ifstream in(path);
        std::string line;
        std::string current_section;   // empty = top-level
        auto trim = [](std::string& s) {
            while (!s.empty() && (s.front() == ' ' || s.front() == '\t' ||
                                  s.front() == '"'))
                s.erase(s.begin());
            while (!s.empty() && (s.back() == ' ' || s.back() == '\t' ||
                                  s.back() == '"'))
                s.pop_back();
        };
        while (std::getline(in, line)) {
            auto pos = line.find('#');
            if (pos != std::string::npos) line = line.substr(0, pos);

            // Section header [name]: track and continue
            std::string stripped = line;
            trim(stripped);
            if (stripped.size() >= 2 &&
                stripped.front() == '[' && stripped.back() == ']') {
                current_section = stripped.substr(1, stripped.size() - 2);
                trim(current_section);
                continue;
            }

            if (line.find('=') == std::string::npos) continue;

            auto eq = line.find('=');
            std::string key = line.substr(0, eq);
            std::string val = line.substr(eq + 1);
            trim(key);
            trim(val);

            if (current_section.empty()) {
                if      (key == "mopac")         toml_mopac = val;
                else if (key == "tleap")         toml_tleap = val;
                else if (key == "ff14sb_params") toml_ff14sb = val;
                else if (key == "tmpdir")        toml_tmpdir = val;
                else if (key == "bmrb_atom_nom") toml_bmrb_atom_nom = val;
                else if (key == "larsen_hbond_grids") toml_larsen_hbond_grid_dir = val;
            } else if (current_section == "databases") {
                if (key == "tensorcs15") toml_tensorcs15_dsn = val;
            }
        }
        OperationLog::Info("RuntimeEnvironment::Load", "read " + path);
    } else {
        OperationLog::Warn("RuntimeEnvironment::Load",
            "no TOML config at " + path + " — using env vars and PATH only");
    }

    // --- Resolve mopac: TOML → PATH → conda default ---

    mopac_ = ResolveBinary(toml_mopac, "mopac");
    if (mopac_.empty()) {
        // Try conda default location
        std::string conda_mopac = "/home/jessica/micromamba/envs/mm/bin/mopac";
        if (fs::exists(conda_mopac)) mopac_ = conda_mopac;
    }

    // --- Resolve tleap: TOML → AMBERHOME/bin/tleap → PATH → conda ---
    tleap_ = toml_tleap;
    if (tleap_.empty() || !fs::exists(tleap_)) {
        const char* amberhome = std::getenv("AMBERHOME");
        if (amberhome) {
            std::string ah_tleap = std::string(amberhome) + "/bin/tleap";
            if (fs::exists(ah_tleap)) tleap_ = ah_tleap;
        }
    }
    if (tleap_.empty()) {
        tleap_ = ResolveBinary("", "tleap");
    }
    if (tleap_.empty()) {
        std::string conda_tleap = "/home/jessica/micromamba/envs/mm/bin/tleap";
        if (fs::exists(conda_tleap)) tleap_ = conda_tleap;
    }

    // --- Resolve data files ---

    // ff14SB params: TOML → env var → NMR_DATA_DIR/ff14sb_params.dat
    if (!toml_ff14sb.empty() && fs::exists(toml_ff14sb)) {
        ff14sb_params_ = toml_ff14sb;
    } else {
        const char* ff_env = std::getenv("NMR_FF14SB_PARAMS");
        if (ff_env && fs::exists(ff_env)) {
            ff14sb_params_ = ff_env;
        } else {
#ifdef NMR_DATA_DIR
            std::string data_path = std::string(NMR_DATA_DIR) + "/ff14sb_params.dat";
            if (fs::exists(data_path))
                ff14sb_params_ = data_path;
#endif
        }
    }

    // --- Resolve directories ---

    if (!toml_tmpdir.empty()) {
        tmpDir_ = toml_tmpdir;
    } else {
        const char* tmp_env = std::getenv("NMR_TMPDIR");
        tmpDir_ = tmp_env ? tmp_env : "/tmp/nmr_shielding";
    }
    fs::create_directories(tmpDir_);

    // bmrb_atom_nom: TOML → env var → empty.
    // Empty is OK; CategoryInfoProjection runs inert in that case.
    if (!toml_bmrb_atom_nom.empty() && fs::exists(toml_bmrb_atom_nom)) {
        bmrb_atom_nom_ = toml_bmrb_atom_nom;
    } else {
        const char* env = std::getenv("NMR_BMRB_ATOM_NOM");
        if (env && fs::exists(env)) bmrb_atom_nom_ = env;
    }

    // tensorcs15 connection string: TOML → env var → empty. Empty is
    // OK; Session::LoadTripeptideDftTable will skip and the calculator
    // returns nullptr at Compute. The DSN is a libpq kv-pair string,
    // not a path, so we don't fs::exists-check it.
    if (!toml_tensorcs15_dsn.empty()) {
        tensorcs15_dsn_ = toml_tensorcs15_dsn;
    } else {
        const char* env = std::getenv("NMR_TENSORCS15_DSN");
        if (env) tensorcs15_dsn_ = env;
    }

    // larsen_hbond_grids directory: TOML → env var → empty. Empty is
    // OK; Session::LoadLarsenHBondGrid will skip. Fs-check the path
    // and only adopt it when it exists (mirrors bmrb_atom_nom).
    if (!toml_larsen_hbond_grid_dir.empty() &&
        fs::exists(toml_larsen_hbond_grid_dir)) {
        larsen_hbond_grid_dir_ = toml_larsen_hbond_grid_dir;
    } else {
        const char* env = std::getenv("NMR_LARSEN_HBOND_GRIDS");
        if (env && fs::exists(env)) larsen_hbond_grid_dir_ = env;
    }

    // --- Mark loaded ---

    loaded_ = true;

    // --- Log complete resolved state ---

    auto status = [](const std::string& v) -> std::string {
        if (v.empty()) return "<not set>";
        return v;
    };

    OperationLog::Info("RuntimeEnvironment::Load",
        "mopac=" + status(mopac_) +
        " tleap=" + status(tleap_) +
        " ff14sb_params=" + status(ff14sb_params_) +
        " tmpdir=" + status(tmpDir_) +
        " bmrb_atom_nom=" + status(bmrb_atom_nom_) +
        " tensorcs15_dsn=" + status(tensorcs15_dsn_) +
        " larsen_hbond_grids=" + status(larsen_hbond_grid_dir_) +
        " guid=" + processGuid_);
}


std::vector<std::string> RuntimeEnvironment::Verify() {
    RequireLoaded();
    std::vector<std::string> missing;
    auto check = [&](const std::string& name, const std::string& val) {
        if (val.empty() || !fs::exists(val))
            missing.push_back(name + " (" + (val.empty() ? "<not set>" : val) + ")");
    };
    check("mopac", mopac_);
    check("ff14sb_params", ff14sb_params_);
    return missing;
}


std::string RuntimeEnvironment::TempFilePath(const std::string& proteinName,
                                              const std::string& suffix) {
    RequireLoaded();
    return tmpDir_ + "/" + processGuid_ + "_" + proteinName + "_" + suffix;
}


const std::string& RuntimeEnvironment::Mopac()          { RequireLoaded(); return mopac_; }
const std::string& RuntimeEnvironment::Tleap()          { RequireLoaded(); return tleap_; }
const std::string& RuntimeEnvironment::Ff14sbParams()  { RequireLoaded(); return ff14sb_params_; }
const std::string& RuntimeEnvironment::TmpDir()        { RequireLoaded(); return tmpDir_; }
const std::string& RuntimeEnvironment::BmrbAtomNom()   { RequireLoaded(); return bmrb_atom_nom_; }
const std::string& RuntimeEnvironment::TensorCs15Dsn() { RequireLoaded(); return tensorcs15_dsn_; }
const std::string& RuntimeEnvironment::LarsenHBondGridDir() { RequireLoaded(); return larsen_hbond_grid_dir_; }

}  // namespace nmr
