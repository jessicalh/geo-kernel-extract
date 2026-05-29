// QtTopology — wraps the bonds + rings + ring_membership topology
// surfaces, mirroring `nmr::LegacyAmberTopology` shape (without the
// FF-specific "LegacyAmber" branding — the reader doesn't carry the
// per-FF substrate generation history, so the name is reduced).
//
// Library mirror: src/LegacyAmberTopology.h. Owns:
//   - QtBond vector (bonds.npy)
//   - QtRing polymorphic vector (rings.npy + ring_membership.npy join)
//   - QtRingMembership vector (ring_membership.npy)
// And builds reverse-index caches at construction:
//   - per-atom bond indices (CSR-equivalent of Atom.bond_indices)
//   - per-atom ring-membership indices
//   - per-ring membership indices (canonical-walk order preserved)
//
// AtomSemantic substrate is NOT held here — it lives on QtAtom
// directly (the projection collapses the library's separate
// AtomSemanticTable array into per-atom fields).
//
// Move-only. Owned by QtProtein; held as `std::unique_ptr<QtTopology>`.

#pragma once

#include "QtBond.h"
#include "QtRing.h"
#include "QtRingMembership.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace h5reader::model {

enum class QtRingAxis : std::uint8_t {
    Ring = 0,
    AromaticRing,
    SaturatedRing,
};

class QtTopology {
public:
    // Construct from typed vectors. Atom count is needed to size the
    // reverse-index caches. Rings vector ownership transfers.
    QtTopology(std::size_t atomCount,
               std::vector<QtBond> bonds,
               std::vector<std::unique_ptr<QtRing>> rings,
               std::vector<QtRingMembership> ringMemberships,
               std::size_t aromaticRingCount,
               std::size_t saturatedRingCount);

    ~QtTopology() = default;
    QtTopology(const QtTopology&) = delete;
    QtTopology& operator=(const QtTopology&) = delete;
    QtTopology(QtTopology&&) = default;
    QtTopology& operator=(QtTopology&&) = default;

    // ----- Bond axis -----
    std::size_t bondCount() const { return bonds_.size(); }
    const QtBond& bondAt(std::size_t i) const { return bonds_[i]; }
    const std::vector<QtBond>& bonds() const { return bonds_; }
    // Bond indices for one atom — CSR-equivalent of Atom.bond_indices.
    // Returns empty span if atomIdx out of range.
    const std::vector<int32_t>& bondIndicesForAtom(std::size_t atomIdx) const;

    // ----- Ring axis (aromatic first, then saturated; rings.npy row order) -----
    std::size_t ringCount() const { return rings_.size(); }
    std::size_t aromaticRingCount() const { return aromaticRingCount_; }
    std::size_t saturatedRingCount() const { return saturatedRingCount_; }
    std::optional<std::size_t> absoluteRingIndex(QtRingAxis axis, std::size_t index) const;
    const QtRing& ringAt(std::size_t i) const { return *rings_[i]; }
    const QtRing* ringAtPtr(std::size_t i) const { return rings_[i].get(); }

    // ----- Ring membership axis -----
    std::size_t ringMembershipCount() const { return ringMemberships_.size(); }
    const QtRingMembership& ringMembershipAt(std::size_t i) const { return ringMemberships_[i]; }
    // Membership indices for one atom (atom may belong to up to 2-3 rings
    // for fused TRP bridgeheads).
    const std::vector<int32_t>& ringMembershipsForAtom(std::size_t atomIdx) const;
    // Membership indices for one ring, ordered by ring_atom_order.
    const std::vector<int32_t>& ringMembershipsForRing(std::size_t ringId) const;

private:
    std::vector<QtBond> bonds_;
    std::vector<std::unique_ptr<QtRing>> rings_;
    std::vector<QtRingMembership> ringMemberships_;
    std::size_t aromaticRingCount_ = 0;
    std::size_t saturatedRingCount_ = 0;

    // Reverse-index caches built at construction
    std::vector<std::vector<int32_t>> bondsByAtom_;
    std::vector<std::vector<int32_t>> ringMembershipsByAtom_;
    std::vector<std::vector<int32_t>> ringMembershipsByRing_;

    static const std::vector<int32_t>& EmptyIndexVector();
};

}  // namespace h5reader::model
