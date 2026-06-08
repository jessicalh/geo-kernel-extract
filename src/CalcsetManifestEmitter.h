#pragma once
//
// Producer-side .lgs writer for calcset manifests.
//
// The manifest is a pointer-only wrapper around artifacts the extractor
// already wrote. It does not duplicate sidecar contents.

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace nmr {

class CalcsetManifestEmitter {
public:
    struct Identity {
        std::string protein_id;
        std::string dataset_id;
        std::string human_name;
    };

    struct Result {
        bool ok = false;
        std::filesystem::path path;
        std::string error;
    };

    struct SinglePoseArtifacts {
        std::string pose_kind;
        std::filesystem::path pose_dir = ".";
        std::filesystem::path extraction_manifest;
    };

    struct TrajectoryArtifacts {
        std::filesystem::path md_dir;
        std::filesystem::path topology_top;
        std::filesystem::path extraction_dir = ".";
        std::filesystem::path trajectory_h5;
        std::filesystem::path extraction_manifest;
        std::optional<double> frame_dt_ps;
        std::string frame_index_basis = "trr_frame_index";
        std::filesystem::path reference_pdb;
    };

    struct MutantPairArtifacts {
        std::filesystem::path wt_lgs;
        std::filesystem::path ala_lgs;
    };

    static Identity IdentityFromProteinId(const std::string& protein_id);

    static std::optional<double> FrameDtPsFromTimes(
        const std::vector<double>& frame_times);

    static Result WriteSinglePose(const std::filesystem::path& root,
                                  const Identity& identity,
                                  const SinglePoseArtifacts& artifacts);

    static Result WriteTrajectory(const std::filesystem::path& root,
                                  const Identity& identity,
                                  const TrajectoryArtifacts& artifacts);

    static Result WriteMutantPair(const std::filesystem::path& root,
                                  const Identity& identity,
                                  const MutantPairArtifacts& artifacts);
};

}  // namespace nmr
