#pragma once
/// @file
/// Command-line argument parsing for nmr_extract.

#include "CommonOptions.h"
#include "ModeSpec.h"

#include <optional>
#include <string>
#include <vector>

namespace nmr::cli {

/// @brief Outcome of parsing argv.
///
/// On success @c spec and @c common are populated and @c error is empty.
/// On parse failure @c error names the problem and @c spec is unset.
/// On @c --help or no arguments @c help_requested is true and the
/// caller should print usage and exit cleanly.
struct ParseResult {
    std::optional<ModeSpec>      spec;
    std::optional<CommonOptions> common;
    std::string                  error;
    std::vector<std::string>     warnings;
    bool                         help_requested = false;
};

/// @brief Parse argv into a @ref ParseResult.
///
/// Does not check file existence; that is @ref Validate's job.
/// @param argc Argument count, unchanged.
/// @param argv Argument vector, unchanged.
/// @return ParseResult describing what was found.
ParseResult Parse(int argc, char* argv[]);

}  // namespace nmr::cli
