#include "CalcsetManifestEmitter.h"
#include "AppVersion.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace nmr {
namespace {

namespace fs = std::filesystem;

constexpr int kSchemaVersion = 1;
constexpr const char* kLgsWriter = "nmr_extract lgs draft 0.1";
constexpr const char* kDefaultDftMethod =
    "r2SCAN def2-SVP def2/J NMR  CPCM(Water)";

struct Timestamp {
    std::string iso_utc;
    std::string filename_utc;
};

struct DftFrame {
    int frame_index = 0;
    fs::path meta_json;
};

Timestamp MakeTimestamp() {
    using namespace std::chrono;
    const auto now = system_clock::now();
    const auto t = system_clock::to_time_t(now);
    std::tm tm{};
    gmtime_r(&t, &tm);

    char iso[32];
    std::strftime(iso, sizeof(iso), "%Y-%m-%dT%H:%M:%SZ", &tm);

    char file[32];
    std::strftime(file, sizeof(file), "%Y%m%dT%H%M%SZ", &tm);

    return {iso, file};
}

std::string FilenameStemFromProteinId(const std::string& protein_id) {
    std::string out = protein_id.empty() ? std::string("calcset") : protein_id;
    for (char& ch : out) {
        if (ch == '/' || ch == '\\') ch = '_';
    }
    return out;
}

std::string GenericString(const fs::path& p) {
    const std::string s = p.generic_string();
    return s.empty() ? std::string(".") : s;
}

fs::path RootRelativeInput(const fs::path& root, const fs::path& target) {
    if (target.empty()) return target;
    if (target.is_absolute()) return target;
    return root / target;
}

std::string RelPath(const fs::path& root, const fs::path& target) {
    const fs::path abs_root = fs::absolute(root);
    const fs::path abs_target = fs::absolute(RootRelativeInput(root, target));

    std::error_code ec;
    fs::path rel = fs::relative(abs_target, abs_root, ec);
    if (ec || rel.empty()) {
        rel = abs_target.lexically_relative(abs_root);
    }
    if (rel.empty()) rel = ".";
    return GenericString(rel);
}

bool RequireDirectory(const fs::path& p, const char* label, std::string& error) {
    if (p.empty()) {
        error = std::string(label) + " path is empty";
        return false;
    }
    std::error_code ec;
    if (!fs::is_directory(p, ec)) {
        error = std::string(label) + " is not a directory: " + p.string();
        return false;
    }
    return true;
}

bool RequireFile(const fs::path& p, const char* label, std::string& error) {
    if (p.empty()) {
        error = std::string(label) + " path is empty";
        return false;
    }
    std::error_code ec;
    if (!fs::is_regular_file(p, ec)) {
        error = std::string(label) + " is not a regular file: " + p.string();
        return false;
    }
    return true;
}

nlohmann::ordered_json BaseManifest(
        const CalcsetManifestEmitter::Identity& identity,
        const std::string& kind) {
    nlohmann::ordered_json j;
    j["schema_version"] = kSchemaVersion;
    j["kind"] = kind;
    j["dataset_id"] = identity.dataset_id;
    j["protein_id"] = identity.protein_id;
    j["human_name"] = identity.human_name;
    return j;
}

void AddMetadata(nlohmann::ordered_json& j,
                 const Timestamp& timestamp) {
    j["metadata"] = nlohmann::ordered_json{
        {"generated_at_utc", timestamp.iso_utc},
        {"lgs_writer", kLgsWriter},
        {"app_version", AppVersion()},
        {"app_version_scope", AppVersionScope()},
    };
}

CalcsetManifestEmitter::Result WriteJson(
        const fs::path& root,
        const CalcsetManifestEmitter::Identity& identity,
        const Timestamp& timestamp,
        const nlohmann::ordered_json& j) {
    std::error_code ec;
    fs::create_directories(root, ec);
    if (ec) {
        return {false, {}, "could not create calcset root " + root.string() +
                           ": " + ec.message()};
    }

    const fs::path path = root /
        (FilenameStemFromProteinId(identity.protein_id) + "_" +
         timestamp.filename_utc + ".lgs");
    std::ofstream out(path);
    if (!out) {
        return {false, path, "could not open " + path.string() + " for write"};
    }

    out << j.dump(2, ' ', false,
                  nlohmann::ordered_json::error_handler_t::replace) << "\n";
    if (!out.good()) {
        return {false, path, "write failed for " + path.string()};
    }
    return {true, path, ""};
}

std::vector<DftFrame> CollectDftFrames(const fs::path& jobs_dir) {
    std::vector<DftFrame> frames;
    std::error_code ec;
    if (!fs::is_directory(jobs_dir, ec)) return frames;

    std::vector<fs::path> job_dirs;
    for (const auto& entry : fs::directory_iterator(jobs_dir, ec)) {
        if (ec) break;
        if (entry.is_directory()) job_dirs.push_back(entry.path());
    }
    std::sort(job_dirs.begin(), job_dirs.end());

    for (const fs::path& job_dir : job_dirs) {
        const std::string job_id = job_dir.filename().string();
        const fs::path meta_path = job_dir / (job_id + "_meta.json");
        if (!fs::is_regular_file(meta_path, ec)) continue;

        std::ifstream in(meta_path);
        if (!in) continue;

        try {
            const auto meta = nlohmann::ordered_json::parse(in);
            if (!meta.contains("frame_index") ||
                !meta["frame_index"].is_number_integer()) {
                continue;
            }
            if (!meta.contains("orca_exit_code") ||
                !meta["orca_exit_code"].is_number_integer() ||
                meta["orca_exit_code"].get<int>() != 0) {
                continue;
            }
            frames.push_back({meta["frame_index"].get<int>(), meta_path});
        } catch (...) {
            continue;
        }
    }

    std::sort(frames.begin(), frames.end(),
        [](const DftFrame& a, const DftFrame& b) {
            if (a.frame_index != b.frame_index) {
                return a.frame_index < b.frame_index;
            }
            return a.meta_json < b.meta_json;
        });
    return frames;
}

nlohmann::ordered_json ComputeFrameStride(const std::vector<int>& frame_indices) {
    if (frame_indices.empty()) {
        return nlohmann::ordered_json{
            {"first", 0}, {"last", 0}, {"step", 1},
        };
    }
    if (frame_indices.size() == 1) {
        return nlohmann::ordered_json{
            {"first", frame_indices.front()},
            {"last", frame_indices.front()},
            {"step", 1},
        };
    }

    std::vector<int> diffs;
    for (std::size_t i = 1; i < frame_indices.size(); ++i) {
        const int diff = frame_indices[i] - frame_indices[i - 1];
        if (diff > 0) diffs.push_back(diff);
    }
    int step = 1;
    if (!diffs.empty()) {
        std::sort(diffs.begin(), diffs.end());
        step = diffs[diffs.size() / 2];
    }

    return nlohmann::ordered_json{
        {"first", frame_indices.front()},
        {"last", frame_indices.back()},
        {"step", step},
    };
}

int ReadCampaignTarget(const fs::path& root, int fallback) {
    const fs::path snapshot = root / "dft" / "_consolidation_snapshot.json";
    std::ifstream in(snapshot);
    if (!in) return fallback;

    try {
        const auto j = nlohmann::ordered_json::parse(in);
        if (j.contains("campaign_target") &&
            j["campaign_target"].is_number_integer()) {
            const int target = j["campaign_target"].get<int>();
            if (target > 0) return target;
        }
    } catch (...) {
        return fallback;
    }
    return fallback;
}

std::optional<nlohmann::ordered_json> BuildDftBlock(const fs::path& root) {
    const fs::path jobs_dir = root / "dft" / "jobs";
    const std::vector<DftFrame> frames = CollectDftFrames(jobs_dir);
    if (frames.empty()) return std::nullopt;

    std::vector<int> frame_indices;
    frame_indices.reserve(frames.size());
    nlohmann::ordered_json frame_array = nlohmann::ordered_json::array();
    for (const DftFrame& frame : frames) {
        frame_indices.push_back(frame.frame_index);
        frame_array.push_back(nlohmann::ordered_json{
            {"frame_index", frame.frame_index},
            {"meta_json", RelPath(root, frame.meta_json)},
        });
    }

    return nlohmann::ordered_json{
        {"method", kDefaultDftMethod},
        {"campaign_target_frames",
            ReadCampaignTarget(root, static_cast<int>(frames.size()))},
        {"frame_stride", ComputeFrameStride(frame_indices)},
        {"frames", frame_array},
    };
}

}  // namespace

CalcsetManifestEmitter::Identity
CalcsetManifestEmitter::IdentityFromProteinId(const std::string& protein_id) {
    const std::string id = protein_id.empty() ? std::string("calcset") : protein_id;
    return {id, id, id + " calcset"};
}

std::optional<double> CalcsetManifestEmitter::FrameDtPsFromTimes(
        const std::vector<double>& frame_times) {
    if (frame_times.size() < 2) return std::nullopt;

    std::vector<double> diffs;
    diffs.reserve(frame_times.size() - 1);
    for (std::size_t i = 1; i < frame_times.size(); ++i) {
        const double diff = frame_times[i] - frame_times[i - 1];
        if (diff > 0.0) diffs.push_back(diff);
    }
    if (diffs.empty()) return std::nullopt;

    std::sort(diffs.begin(), diffs.end());
    return diffs[diffs.size() / 2];
}

CalcsetManifestEmitter::Result CalcsetManifestEmitter::WriteSinglePose(
        const fs::path& root,
        const Identity& identity,
        const SinglePoseArtifacts& artifacts) {
    std::string error;
    const fs::path pose_dir = RootRelativeInput(root, artifacts.pose_dir);
    const fs::path extraction_manifest =
        RootRelativeInput(root, artifacts.extraction_manifest);

    if (!RequireDirectory(root, "calcset root", error)) return {false, {}, error};
    if (!RequireDirectory(pose_dir, "single_pose.pose_dir", error)) {
        return {false, {}, error};
    }
    if (!RequireFile(extraction_manifest,
                     "single_pose.extraction_manifest", error)) {
        return {false, {}, error};
    }

    const Timestamp timestamp = MakeTimestamp();
    nlohmann::ordered_json j = BaseManifest(identity, "single_pose");
    j["single_pose"] = nlohmann::ordered_json{
        {"pose_kind", artifacts.pose_kind},
        {"pose_dir", RelPath(root, pose_dir)},
        {"extraction_manifest", RelPath(root, extraction_manifest)},
    };
    AddMetadata(j, timestamp);

    return WriteJson(root, identity, timestamp, j);
}

CalcsetManifestEmitter::Result CalcsetManifestEmitter::WriteTrajectory(
        const fs::path& root,
        const Identity& identity,
        const TrajectoryArtifacts& artifacts) {
    std::string error;
    const fs::path md_dir = RootRelativeInput(root, artifacts.md_dir);
    const fs::path topology_top = RootRelativeInput(root, artifacts.topology_top);
    const fs::path extraction_dir =
        RootRelativeInput(root, artifacts.extraction_dir);
    const fs::path trajectory_h5 =
        RootRelativeInput(root, artifacts.trajectory_h5);
    const fs::path extraction_manifest =
        RootRelativeInput(root, artifacts.extraction_manifest);

    if (!RequireDirectory(root, "calcset root", error)) return {false, {}, error};
    if (!RequireDirectory(md_dir, "trajectory.md_dir", error)) {
        return {false, {}, error};
    }
    if (!RequireFile(topology_top, "trajectory.topology_top", error)) {
        return {false, {}, error};
    }
    if (!RequireDirectory(extraction_dir, "trajectory.extraction_dir", error)) {
        return {false, {}, error};
    }
    if (!RequireFile(trajectory_h5, "trajectory.trajectory_h5", error)) {
        return {false, {}, error};
    }
    if (!RequireFile(extraction_manifest,
                     "trajectory.extraction_manifest", error)) {
        return {false, {}, error};
    }
    if (!artifacts.frame_dt_ps || *artifacts.frame_dt_ps <= 0.0) {
        return {false, {},
                "trajectory.frame_dt_ps could not be derived from recorded times"};
    }

    const Timestamp timestamp = MakeTimestamp();
    nlohmann::ordered_json j = BaseManifest(identity, "trajectory");
    j["trajectory"] = nlohmann::ordered_json{
        {"md_dir", RelPath(root, md_dir)},
        {"topology_top", RelPath(root, topology_top)},
        {"extraction_dir", RelPath(root, extraction_dir)},
        {"trajectory_h5", RelPath(root, trajectory_h5)},
        {"extraction_manifest", RelPath(root, extraction_manifest)},
        {"frame_dt_ps", *artifacts.frame_dt_ps},
        {"frame_index_basis", artifacts.frame_index_basis},
    };

    const fs::path reference_pdb = RootRelativeInput(root, artifacts.reference_pdb);
    std::error_code ec;
    if (!artifacts.reference_pdb.empty() && fs::is_regular_file(reference_pdb, ec)) {
        j["trajectory"]["reference_pdb"] = RelPath(root, reference_pdb);
    }

    if (auto dft = BuildDftBlock(root)) {
        j["dft"] = *dft;
    }
    AddMetadata(j, timestamp);

    return WriteJson(root, identity, timestamp, j);
}

CalcsetManifestEmitter::Result CalcsetManifestEmitter::WriteMutantPair(
        const fs::path& root,
        const Identity& identity,
        const MutantPairArtifacts& artifacts) {
    std::string error;
    const fs::path wt_lgs = RootRelativeInput(root, artifacts.wt_lgs);
    const fs::path ala_lgs = RootRelativeInput(root, artifacts.ala_lgs);

    if (!RequireDirectory(root, "calcset root", error)) return {false, {}, error};
    if (!RequireFile(wt_lgs, "mutant_pair.wt_lgs", error)) {
        return {false, {}, error};
    }
    if (!RequireFile(ala_lgs, "mutant_pair.ala_lgs", error)) {
        return {false, {}, error};
    }

    const Timestamp timestamp = MakeTimestamp();
    nlohmann::ordered_json j = BaseManifest(identity, "mutant_pair");
    j["mutant_pair"] = nlohmann::ordered_json{
        {"wt_lgs", RelPath(root, wt_lgs)},
        {"ala_lgs", RelPath(root, ala_lgs)},
    };
    AddMetadata(j, timestamp);

    return WriteJson(root, identity, timestamp, j);
}

}  // namespace nmr
