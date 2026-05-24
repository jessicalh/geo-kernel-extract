#pragma once
/// @file
/// Per-frame sidecar emission options for trajectory mode.

#include <cstddef>
#include <filesystem>
#include <limits>
#include <string>

namespace nmr::cli {

/// @brief Per-frame PDB sidecar emission for trajectory runs.
///
/// Consumed by @c FramePdbEmitter::Configure. Stem is derived from the
/// trajectory directory's basename, not supplied here.
struct FramePdbEmission {
    std::filesystem::path output_dir;          ///< Sidecar output directory.
    std::size_t           stride = 1;          ///< Emit every N-th frame as read.
    double                from_ps =
        -std::numeric_limits<double>::infinity();  ///< Lower time bound (inclusive).
    double                to_ps =
         std::numeric_limits<double>::infinity();  ///< Upper time bound (exclusive).
    std::string           decorator;           ///< Optional tag in the filename stem.
};

/// @brief Per-frame NPY sidecar emission for trajectory runs.
///
/// Each accepted frame becomes @c output_dir/frame_NNNNNN/ with the full
/// per-conformation NPY set. Per-protein topology sidecars land once at
/// @c output_dir. Consumed by @c FrameNpyEmitter::Configure.
struct FrameNpyEmission {
    std::filesystem::path output_dir;          ///< Sidecar output directory.
    std::size_t           stride = 1;          ///< Emit every N-th frame as read.
    double                from_ps =
        -std::numeric_limits<double>::infinity();
    double                to_ps =
         std::numeric_limits<double>::infinity();
};

}  // namespace nmr::cli
