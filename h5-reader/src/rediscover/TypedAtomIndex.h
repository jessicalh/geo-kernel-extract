// TypedAtomIndex — scoped, set-returning atom lookup over QtAtom identity.
//
// This is deliberately not a name or position lookup. Selectors are partial
// typed predicates, scopes are explicit atom-index sets, and asserted-unique
// lookup fails loud when chemistry yields 0 or N matches.

#pragma once

#include "../model/QtSemanticEnums.h"
#include "../model/Types.h"

#include <QString>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace h5reader::model {
class QtProtein;
}

namespace h5reader::rediscover {

struct TypedAtomSelector {
    std::optional<model::Element> element;
    std::optional<model::Locant> locant;
    std::optional<model::BranchAddress> branch;
    std::optional<model::DiastereotopicIndex> diIndex;
    std::optional<model::BackboneRole> backboneRole;
};

class TypedAtomIndex {
public:
    TypedAtomIndex() = default;
    explicit TypedAtomIndex(const model::QtProtein& protein);

    std::vector<int32_t> residueScope(std::size_t residueIndex) const;
    std::vector<int32_t> select(const std::vector<int32_t>& scope,
                                const TypedAtomSelector& selector) const;
    std::optional<int32_t> selectUnique(const std::vector<int32_t>& scope,
                                        const TypedAtomSelector& selector,
                                        QString* err_out = nullptr) const;

private:
    const model::QtProtein* protein_ = nullptr;
};

}  // namespace h5reader::rediscover
