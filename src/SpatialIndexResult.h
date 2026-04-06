#pragma once
//
// SpatialIndexResult: nanoflann KD-tree spatial index over atom positions.
//
// Dependencies: GeometryResult.
//
// Builds KD-trees for: atom positions, ring centers, bond midpoints.
// Provides: radius search, K-nearest search.
// Builds 15A neighbour lists: for each atom, all atoms within 15A with
// stored distance (double) and direction (Vec3, normalised).
//
// Neighbours stored on ConformationAtom::spatial_neighbours.
//
// Constitution: 15A cutoff for ring current calculations.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "GeometryResult.h"
#include "nanoflann.hpp"
#include <vector>
#include <typeindex>

namespace nmr {

// Constitution: 15A cutoff where ring current effects become negligible
constexpr double SPATIAL_NEIGHBOUR_CUTOFF_A = 15.0;

// ============================================================================
// Point cloud adaptors for nanoflann
// ============================================================================

struct PointCloud {
    std::vector<Vec3> points;

    size_t kdtree_get_point_count() const { return points.size(); }

    double kdtree_get_pt(size_t idx, size_t dim) const {
        return points[idx](static_cast<int>(dim));
    }

    template<class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud>,
    PointCloud, 3, size_t>;


class SpatialIndexResult : public ConformationResult {
public:
    std::string Name() const override { return "SpatialIndexResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { typeid(GeometryResult) };
    }

    // Factory: build KD-trees and neighbour lists
    static std::unique_ptr<SpatialIndexResult> Compute(ProteinConformation& conf);

    // Radius search on atom positions: returns indices within radius
    std::vector<size_t> AtomsWithinRadius(Vec3 point, double radius) const;

    // K-nearest atoms to a point
    std::vector<size_t> KNearestAtoms(Vec3 point, size_t k) const;

    // Radius search on ring centers
    std::vector<size_t> RingsWithinRadius(Vec3 point, double radius) const;

    // Radius search on bond midpoints
    std::vector<size_t> BondsWithinRadius(Vec3 point, double radius) const;

private:
    ProteinConformation* conf_ = nullptr;

    PointCloud atom_cloud_;
    PointCloud ring_cloud_;
    PointCloud bond_cloud_;

    std::unique_ptr<KDTree> atom_tree_;
    std::unique_ptr<KDTree> ring_tree_;
    std::unique_ptr<KDTree> bond_tree_;
};

}  // namespace nmr
