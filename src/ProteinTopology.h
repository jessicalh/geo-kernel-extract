#pragma once
//
// ProteinTopology: intentionally minimal abstract base — calculators
// bind to a concrete topology (e.g. LegacyAmberTopology), not to a
// broad virtual interface here.
//

#include <cstddef>
#include <string_view>

namespace nmr {

class IupacAtomMap;

enum class ProteinTopologyKind {
    LegacyAmber
};

class ProteinTopology {
public:
    virtual ~ProteinTopology() = default;

    virtual ProteinTopologyKind Kind() const = 0;
    virtual std::string_view Name() const = 0;
    virtual size_t AtomCount() const = 0;
    virtual size_t ResidueCount() const = 0;

    virtual const IupacAtomMap* IupacOrNull() const { return nullptr; }
};

}  // namespace nmr
