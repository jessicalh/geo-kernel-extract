#pragma once
/// @file
/// File-existence validation for a parsed @ref nmr::cli::ModeSpec.

#include "CommonOptions.h"
#include "ModeSpec.h"

#include <string>

namespace nmr::cli {

/// @brief Verify every required file referenced by @c spec exists.
///
/// Returns an empty string on success. On failure, returns a single
/// diagnostic naming the first missing or invalid file with its full
/// path. Call after @ref Parse, before any pipeline construction.
///
/// Optional files (e.g. ORCA NMR output) that are absent do not fail
/// here; the calling pipeline decides how to handle their absence.
std::string Validate(const ModeSpec& spec, const CommonOptions& common);

}  // namespace nmr::cli
