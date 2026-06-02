// SpatialIndexSet — day-one KD trees for all source-cloud kinds.

#pragma once

#include "../model/Types.h"

#include <nanoflann/nanoflann.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace h5reader::rediscover {

class RunData;

using model::Vec3;

enum class CloudKind : int {
    Atoms = 0,
    BondMidpoints = 1,     // anisotropic McConnell subset, legacy callers
    RingCenters = 2,
    ChargeSites = 3,
    AllBondMidpoints = 4,  // full producer BondCategory enumeration
};

struct SourceRef {
    CloudKind kind = CloudKind::Atoms;
    int32_t cloud_index = -1;
    int32_t entity_index = -1;  // atom index, bond index, ring index, or charge-site atom index
};

struct CloudPoint {
    Vec3 point = Vec3::Zero();
    SourceRef ref;
};

class CloudTree {
public:
    CloudTree() = default;
    CloudTree(CloudKind kind, std::vector<CloudPoint> points);

    std::vector<SourceRef> near(const Vec3& query, double cutoff) const;
    std::vector<SourceRef> range(const Vec3& query, double innerCutoff, double outerCutoff) const;
    std::vector<SourceRef> nearAtom(std::size_t atom, double cutoff) const;
    std::vector<SourceRef> rangeAtom(std::size_t atom, double innerCutoff, double outerCutoff) const;

    const Vec3& pointAt(std::size_t cloudIndex) const { return points_[cloudIndex].point; }
    std::size_t size() const { return points_.size(); }

    std::size_t kdtree_get_point_count() const { return points_.size(); }
    double kdtree_get_pt(std::size_t idx, std::size_t dim) const {
        return points_[idx].point[static_cast<int>(dim)];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }

private:
    using KdTree = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, CloudTree>, CloudTree, 3, std::size_t>;

    CloudKind kind_ = CloudKind::Atoms;
    std::vector<CloudPoint> points_;
    std::unique_ptr<KdTree> tree_;
};

class SpatialIndexSet {
public:
    SpatialIndexSet() = default;
    explicit SpatialIndexSet(const RunData& run);

    const CloudTree& tree(CloudKind kind, std::size_t frame) const;
    std::vector<SourceRef> near(CloudKind kind, std::size_t frame, const Vec3& query,
                                double cutoff) const;
    std::vector<SourceRef> range(CloudKind kind, std::size_t frame, const Vec3& query,
                                 double innerCutoff, double outerCutoff) const;

    std::size_t frameCount() const { return frameCount_; }

private:
    std::size_t frameCount_ = 0;
    std::array<std::vector<CloudTree>, 5> trees_;
};

}  // namespace h5reader::rediscover
