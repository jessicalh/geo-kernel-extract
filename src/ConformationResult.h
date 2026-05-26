#pragma once
//
// ConformationResult: base class for all computed results that attach to
// a ProteinConformation. Each result type is a named singleton.
//
// Access by template: conf.Result<T>(), conf.HasResult<T>()
// Dependency checking at attach time.
//
// WriteFeatures: each result knows what it computed and writes its own
// contribution to disk. The static WriteAllFeatures traverses the
// conformation's accumulated results and calls each one.
//

#include <string>
#include <vector>
#include <typeindex>
#include <memory>

namespace nmr {

class ProteinConformation;

class ConformationResult {
public:
    virtual ~ConformationResult() = default;

    // Human-readable name for logging and diagnostics
    virtual std::string Name() const = 0;

    // Dependencies: type_index values of results that must be attached first.
    virtual std::vector<std::type_index> Dependencies() const = 0;

    // Write this result's features to .npy files in output_dir.
    // Returns the number of arrays written.
    virtual int WriteFeatures(const ProteinConformation& conf,
                              const std::string& output_dir) const { return 0; }

    // Calls each attached result's WriteFeatures (iterated in
    // unordered_map order; dependency order is enforced at attach time,
    // not here). Also writes the identity arrays and ring tables that
    // belong to no single result.
    static int WriteAllFeatures(const ProteinConformation& conf,
                                const std::string& output_dir);
};

}  // namespace nmr
