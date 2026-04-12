#pragma once
//
// CalculatorConfig: tuneable calculator parameters from TOML.
//
// Load() once at startup. Get() returns a double — the TOML value
// if present, else the compiled literature default. That is all.
//
// The TOML file is flat key = value. No sections, no nesting.
// Missing file → all defaults. Missing key → default for that key.
//
// Validate() checks all expected keys, logs every value and its
// source (toml or default), returns a list of any unknown keys
// found in the file. Call once before calculators run.
//

#include <string>
#include <vector>
#include <unordered_map>

namespace nmr {

class CalculatorConfig {
public:
    // Load parameter overrides from a TOML file.
    // Empty path → defaults only (no file read).
    static void Load(const std::string& path = "");

    // Get a numeric parameter value. Returns override if loaded, else default.
    // Aborts on unknown key (programming error).
    static double Get(const std::string& key);

    // Get a string value from the config file (e.g. model paths).
    // Returns defaultVal if the key was not set.
    static std::string GetString(const std::string& key,
                                 const std::string& defaultVal = "");

    // Validate the loaded configuration.
    // Logs every parameter: key, value, source (toml/default), unit.
    // Returns list of unknown keys found in the TOML file (typos).
    // Empty return = all clean.
    static std::vector<std::string> Validate();

    static bool IsLoaded() { return loaded_; }

private:
    static void InitDefaults();

    struct ParamEntry {
        double value;
        const char* unit;
        const char* description;
    };

    static std::unordered_map<std::string, ParamEntry> defaults_;
    static std::unordered_map<std::string, double> overrides_;
    static std::unordered_map<std::string, std::string> string_overrides_;
    static bool loaded_;
    static bool defaults_initialised_;
};

}  // namespace nmr
