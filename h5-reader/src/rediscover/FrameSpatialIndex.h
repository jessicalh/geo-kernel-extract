// FrameSpatialIndex — a per-frame nanoflann KD-tree over anisotropic-bond
// midpoints, for McConnell source discovery. Built lazily per frame (cheap:
// a few hundred bonds). Discovery is by documented cutoff (no glob, no
// try-and-fail); the cutoff is a required parameter at the call site per the
// substrate conventions (no default cutoffs).
//
// The cloud is the subset of bonds that carry magnetic anisotropy in the
// McConnell sense — the four bond categories the design names (PeptideCO,
// PeptideCN, SidechainCO, Aromatic). Each kept bond contributes its midpoint
// + its endpoint positions so the caller can form the bond axis on the fly.

#pragma once

#include <nanoflann/nanoflann.hpp>

#include "../model/Types.h"  // Vec3

#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::model {
class QtProtein;
class Conformation;
}

namespace h5reader::rediscover {

using model::Vec3;

// One anisotropic bond's per-frame geometry, kept alongside its identity.
struct AnisoBond {
    int32_t bondIndex = -1;
    int32_t atomIndexA = -1;  // endpoint atom indices (provenance + equivariant fit)
    int32_t atomIndexB = -1;
    int     category  = -1;  // BondCategory ordinal
    int     order     = -1;  // BondOrder ordinal
    int     elemA     = -1;  // Element ordinal of endpoint A
    int     elemB     = -1;
    Vec3    posA = Vec3::Zero();
    Vec3    posB = Vec3::Zero();
    Vec3    midpoint = Vec3::Zero();
};

class FrameSpatialIndex {
public:
    // Build the index for one frame: collects the anisotropic bonds (the four
    // McConnell categories) and their midpoints at this frame's positions,
    // then builds the KD-tree.
    FrameSpatialIndex(const model::QtProtein& protein, const model::Conformation& conf,
                      std::size_t frame);

    // Bond cloud-indices whose midpoint is within `cutoff` Å of `query`.
    // cutoff is REQUIRED (no default) per the substrate conventions.
    std::vector<std::size_t> Within(const Vec3& query, double cutoff) const;

    const AnisoBond& bond(std::size_t cloudIndex) const { return bonds_[cloudIndex]; }
    std::size_t bondCount() const { return bonds_.size(); }

    // ── nanoflann point-cloud adaptor surface ──
    std::size_t kdtree_get_point_count() const { return bonds_.size(); }
    double kdtree_get_pt(std::size_t idx, std::size_t dim) const {
        return bonds_[idx].midpoint[static_cast<int>(dim)];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }

private:
    using KdTree = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, FrameSpatialIndex>, FrameSpatialIndex, 3, std::size_t>;

    std::vector<AnisoBond> bonds_;
    std::unique_ptr<KdTree> tree_;
};

}  // namespace h5reader::rediscover
