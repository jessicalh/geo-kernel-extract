/// @file
/// nmr::cli::Validate implementation.

#include "Cli/Validate.h"

#include <filesystem>
#include <variant>

namespace nmr::cli {

namespace {

namespace fs = std::filesystem;

std::string FileError(const fs::path& p, const char* role) {
    return std::string(role) + ": file not found: " + p.string();
}

std::string DirError(const fs::path& p, const char* role) {
    return std::string(role) + ": directory not found: " + p.string();
}

std::string CheckFile(const fs::path& p, const char* role) {
    if (p.empty())        return std::string(role) + ": no path specified";
    if (!fs::exists(p))   return FileError(p, role);
    if (!fs::is_regular_file(p)) return std::string(role) + ": not a regular file: " + p.string();
    return "";
}

std::string CheckDir(const fs::path& p, const char* role) {
    if (p.empty())                return std::string(role) + ": no path specified";
    if (!fs::exists(p))           return DirError(p, role);
    if (!fs::is_directory(p))     return std::string(role) + ": not a directory: " + p.string();
    return "";
}

std::string CheckOrcaFiles(const OrcaRunFiles& f, const char* side) {
    std::string e;
    e = CheckFile(f.xyz_path,    (std::string(side) + " .xyz").c_str());      if (!e.empty()) return e;
    e = CheckFile(f.prmtop_path, (std::string(side) + " .prmtop").c_str());   if (!e.empty()) return e;
    // nmr.out is optional; absence is not an error here.
    return "";
}

std::string CheckCommon(const CommonOptions& c) {
    // output_dir may be created at runtime; only check if non-empty
    // that it is not a path to a non-directory.
    if (!c.output_dir.empty() && fs::exists(c.output_dir) && !fs::is_directory(c.output_dir)) {
        return "--output: exists but is not a directory: " + c.output_dir.string();
    }
    if (!c.config_path.empty()) {
        const std::string e = CheckFile(c.config_path, "--config");
        if (!e.empty()) return e;
    }
    if (!c.aimnet2_model_path.empty()) {
        const std::string e = CheckFile(c.aimnet2_model_path, "--aimnet2 model");
        if (!e.empty()) return e;
    }
    return "";
}

std::string CheckTrajectory(const TrajectoryMode& m) {
    std::string e = CheckDir(m.dir, "--trajectory DIR");
    if (!e.empty()) return e;
    const auto files = TrajectoryInputFiles::FromProductionDir(m.dir);
    e = CheckFile(files.tpr, "trajectory production.tpr"); if (!e.empty()) return e;
    e = CheckFile(files.trr, "trajectory production.trr"); if (!e.empty()) return e;
    e = CheckFile(files.edr, "trajectory production.edr"); if (!e.empty()) return e;
    return "";
}

}  // namespace


std::string Validate(const ModeSpec& spec, const CommonOptions& common) {
    if (const std::string e = CheckCommon(common); !e.empty()) return e;

    return std::visit(
        [&](const auto& m) -> std::string {
            using T = std::decay_t<decltype(m)>;
            if constexpr (std::is_same_v<T, PdbMode>) {
                return CheckFile(m.pdb, "--pdb");
            } else if constexpr (std::is_same_v<T, ProtonatedPdbMode>) {
                return CheckFile(m.pdb, "--protonated-pdb");
            } else if constexpr (std::is_same_v<T, OrcaMode>) {
                return CheckOrcaFiles(m.files, "--orca");
            } else if constexpr (std::is_same_v<T, MutantMode>) {
                std::string e = CheckOrcaFiles(m.wt,  "--mutant --wt");
                if (!e.empty()) return e;
                return CheckOrcaFiles(m.ala, "--mutant --ala");
            } else if constexpr (std::is_same_v<T, TrajectoryMode>) {
                return CheckTrajectory(m);
            } else {
                return "internal: unhandled mode variant";
            }
        },
        spec);
}

}  // namespace nmr::cli
