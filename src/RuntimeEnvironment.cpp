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
std::string RuntimeEnvironment::ff14sb_params_;
std::string RuntimeEnvironment::tmpDir_;
std::string RuntimeEnvironment::ccd_path_;
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

    std::string toml_mopac, toml_ff14sb, toml_tmpdir, toml_ccd;

    if (!path.empty() && fs::exists(path)) {
        std::ifstream in(path);
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

            if      (key == "mopac")          toml_mopac = val;
            else if (key == "ff14sb_params") toml_ff14sb = val;
            else if (key == "tmpdir")        toml_tmpdir = val;
            else if (key == "ccd_path")      toml_ccd   = val;
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

    // CCD: TOML > env var > NMR_DATA_DIR/ccd/components.cif.gz.
    // No file discovery — if none of these resolve to an existing file,
    // ccd_path_ stays empty and CcdValidator reports ccd_loaded=false.
    if (!toml_ccd.empty() && fs::exists(toml_ccd)) {
        ccd_path_ = toml_ccd;
    } else {
        const char* ccd_env = std::getenv("NMR_CCD_PATH");
        if (ccd_env && fs::exists(ccd_env)) {
            ccd_path_ = ccd_env;
        } else {
#ifdef NMR_DATA_DIR
            std::string default_path = std::string(NMR_DATA_DIR) + "/ccd/components.cif";
            if (fs::exists(default_path))
                ccd_path_ = default_path;
#endif
        }
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
        " ff14sb_params=" + status(ff14sb_params_) +
        " tmpdir=" + status(tmpDir_) +
        " ccd_path=" + status(ccd_path_) +
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
const std::string& RuntimeEnvironment::Ff14sbParams()  { RequireLoaded(); return ff14sb_params_; }
const std::string& RuntimeEnvironment::TmpDir()        { RequireLoaded(); return tmpDir_; }
const std::string& RuntimeEnvironment::CcdPath()       { RequireLoaded(); return ccd_path_; }

}  // namespace nmr
