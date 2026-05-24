#pragma once
/// @file
/// Help-text emission for nmr_extract.

namespace nmr::cli {

/// @brief Print the program usage to stderr.
///
/// Organised per-mode: each supported mode appears with its required
/// arguments, optional flags, and a one-line description. Common
/// options appear once at the end.
/// @param program_name The argv[0] value to embed in the synopsis.
void PrintUsage(const char* program_name);

}  // namespace nmr::cli
