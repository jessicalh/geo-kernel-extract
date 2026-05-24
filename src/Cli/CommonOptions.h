#pragma once
/// @file
/// Mode-agnostic command-line options shared by every nmr_extract mode.

#include <filesystem>

namespace nmr::cli {

/// @brief Options every mode accepts.
///
/// Distinct from the per-mode @ref ModeSpec variant: these settings
/// describe how the program runs rather than what work it performs.
/// @see ModeSpec
struct CommonOptions {
    /// Output directory for NPY feature arrays and trajectory H5.
    /// Empty in viewer-only contexts.
    std::filesystem::path output_dir;

    /// TOML file with calculator parameter overrides.
    /// Empty means use built-in literature defaults.
    std::filesystem::path config_path;

    /// AIMNet2 TorchScript (.jpt) model path.
    /// Empty disables AIMNet2-derived calculators.
    std::filesystem::path aimnet2_model_path;
};

}  // namespace nmr::cli
