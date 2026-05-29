// QtTopology implementation — reverse-index caches built at
// construction time so per-atom / per-ring lookups are O(1) at runtime.

#include "QtTopology.h"

#include <algorithm>

namespace h5reader::model {

QtTopology::QtTopology(std::size_t atomCount,
                       std::vector<QtBond> bonds,
                       std::vector<std::unique_ptr<QtRing>> rings,
                       std::vector<QtRingMembership> ringMemberships,
                       std::size_t aromaticRingCount,
                       std::size_t saturatedRingCount)
    : bonds_(std::move(bonds))
    , rings_(std::move(rings))
    , ringMemberships_(std::move(ringMemberships))
    , aromaticRingCount_(aromaticRingCount)
    , saturatedRingCount_(saturatedRingCount) {
    // bondsByAtom_: walk bonds, push bond index into each endpoint's bucket.
    bondsByAtom_.resize(atomCount);
    for (std::size_t bi = 0; bi < bonds_.size(); ++bi) {
        const auto& b = bonds_[bi];
        if (b.atomIndexA >= 0 && static_cast<std::size_t>(b.atomIndexA) < atomCount)
            bondsByAtom_[static_cast<std::size_t>(b.atomIndexA)].push_back(static_cast<int32_t>(bi));
        if (b.atomIndexB >= 0 && static_cast<std::size_t>(b.atomIndexB) < atomCount)
            bondsByAtom_[static_cast<std::size_t>(b.atomIndexB)].push_back(static_cast<int32_t>(bi));
    }

    // ringMembershipsByAtom_: walk memberships, push membership index into
    // the atom's bucket.
    ringMembershipsByAtom_.resize(atomCount);
    ringMembershipsByRing_.resize(rings_.size());
    for (std::size_t mi = 0; mi < ringMemberships_.size(); ++mi) {
        const auto& m = ringMemberships_[mi];
        if (m.atomIndex >= 0 && static_cast<std::size_t>(m.atomIndex) < atomCount)
            ringMembershipsByAtom_[static_cast<std::size_t>(m.atomIndex)].push_back(static_cast<int32_t>(mi));
        if (m.ringId >= 0 && static_cast<std::size_t>(m.ringId) < rings_.size())
            ringMembershipsByRing_[static_cast<std::size_t>(m.ringId)].push_back(static_cast<int32_t>(mi));
    }

    // Sort ring-membership-by-ring lists by ring_atom_order to preserve
    // canonical walk order (matters for ring polygon overlay vertex order
    // + Johnson-Bovey loop integral orientation).
    for (auto& bucket : ringMembershipsByRing_) {
        std::sort(bucket.begin(), bucket.end(), [this](int32_t a, int32_t b) {
            return ringMemberships_[static_cast<std::size_t>(a)].ringAtomOrder
                   < ringMemberships_[static_cast<std::size_t>(b)].ringAtomOrder;
        });
    }
}

const std::vector<int32_t>& QtTopology::EmptyIndexVector() {
    static const std::vector<int32_t> kEmpty;
    return kEmpty;
}

const std::vector<int32_t>& QtTopology::bondIndicesForAtom(std::size_t atomIdx) const {
    return atomIdx < bondsByAtom_.size() ? bondsByAtom_[atomIdx] : EmptyIndexVector();
}

std::optional<std::size_t> QtTopology::absoluteRingIndex(QtRingAxis axis, std::size_t index) const {
    switch (axis) {
    case QtRingAxis::Ring:
        return index < rings_.size() ? std::optional<std::size_t>(index) : std::nullopt;
    case QtRingAxis::AromaticRing:
        return index < aromaticRingCount_ ? std::optional<std::size_t>(index) : std::nullopt;
    case QtRingAxis::SaturatedRing: {
        if (index >= saturatedRingCount_)
            return std::nullopt;
        const std::size_t absolute = aromaticRingCount_ + index;
        return absolute < rings_.size() ? std::optional<std::size_t>(absolute) : std::nullopt;
    }
    }
    return std::nullopt;
}

const std::vector<int32_t>& QtTopology::ringMembershipsForAtom(std::size_t atomIdx) const {
    return atomIdx < ringMembershipsByAtom_.size() ? ringMembershipsByAtom_[atomIdx] : EmptyIndexVector();
}

const std::vector<int32_t>& QtTopology::ringMembershipsForRing(std::size_t ringId) const {
    return ringId < ringMembershipsByRing_.size() ? ringMembershipsByRing_[ringId] : EmptyIndexVector();
}

}  // namespace h5reader::model
