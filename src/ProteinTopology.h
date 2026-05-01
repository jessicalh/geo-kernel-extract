#pragma once
//
// ProteinTopology: abstract family for named protein topology contracts.
//
// This base class is intentionally small. Concrete topology classes expose
// their own calculator-language accessors; calculators bind to those concrete
// contracts rather than to a broad virtual interface.
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
