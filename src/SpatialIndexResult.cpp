#include "SpatialIndexResult.h"
#include "Protein.h"
#include "OperationLog.h"
#include <cmath>

namespace nmr {

std::unique_ptr<SpatialIndexResult> SpatialIndexResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("SpatialIndexResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    auto result = std::make_unique<SpatialIndexResult>();
    result->conf_ = &conf;

    const Protein& protein = conf.ProteinRef();
    const auto& positions = conf.Positions();

    // ---------------------------------------------------------------
    // Build atom position KD-tree
    // ---------------------------------------------------------------
    result->atom_cloud_.points = positions;
    result->atom_tree_ = std::make_unique<KDTree>(
        3, result->atom_cloud_,
        nanoflann::KDTreeSingleIndexAdaptorParams(10 /* leaf_max_size */));
    result->atom_tree_->buildIndex();

    // ---------------------------------------------------------------
    // Build ring center KD-tree
    // ---------------------------------------------------------------
    result->ring_cloud_.points.resize(protein.RingCount());
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        result->ring_cloud_.points[ri] = conf.ring_geometries[ri].center;
    }
    if (!result->ring_cloud_.points.empty()) {
        result->ring_tree_ = std::make_unique<KDTree>(
            3, result->ring_cloud_,
            nanoflann::KDTreeSingleIndexAdaptorParams(10));
        result->ring_tree_->buildIndex();
    }

    // ---------------------------------------------------------------
    // Build bond midpoint KD-tree
    // ---------------------------------------------------------------
    result->bond_cloud_.points = conf.bond_midpoints;
    if (!result->bond_cloud_.points.empty()) {
        result->bond_tree_ = std::make_unique<KDTree>(
            3, result->bond_cloud_,
            nanoflann::KDTreeSingleIndexAdaptorParams(10));
        result->bond_tree_->buildIndex();
    }

    // ---------------------------------------------------------------
    // Build 15A neighbour lists for every atom.
    // For each atom, find all other atoms within SPATIAL_NEIGHBOUR_CUTOFF_A.
    // Store distance and normalised direction.
    // Self is NOT included.
    // ---------------------------------------------------------------
    const double cutoff_sq = SPATIAL_NEIGHBOUR_CUTOFF_A * SPATIAL_NEIGHBOUR_CUTOFF_A;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const Vec3& pos = positions[ai];
        double query_pt[3] = { pos.x(), pos.y(), pos.z() };

        // nanoflann radius search: returns pairs of (index, squared_distance)
        std::vector<nanoflann::ResultItem<size_t, double>> matches;
        result->atom_tree_->radiusSearch(query_pt, cutoff_sq, matches);

        auto& neighbours = conf.MutableAtomAt(ai).spatial_neighbours;
        neighbours.clear();
        neighbours.reserve(matches.size());

        for (const auto& match : matches) {
            if (match.first == ai) continue;  // skip self

            double dist = std::sqrt(match.second);
            if (dist < 1e-12) continue;  // degenerate

            Vec3 direction = (positions[match.first] - pos) / dist;

            AtomNeighbour nb;
            nb.atom_index = match.first;
            nb.distance = dist;
            nb.direction = direction;
            neighbours.push_back(nb);
        }
    }

    // Diagnostics
    size_t total_neighbours = 0;
    size_t min_neighbours = SIZE_MAX;
    size_t max_neighbours = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        size_t n = conf.AtomAt(ai).spatial_neighbours.size();
        total_neighbours += n;
        if (n < min_neighbours) min_neighbours = n;
        if (n > max_neighbours) max_neighbours = n;
    }

    OperationLog::Info(LogResultAttach, "SpatialIndexResult::Compute",
        "total_neighbour_pairs=" + std::to_string(total_neighbours) +
        " min=" + std::to_string(min_neighbours) +
        " max=" + std::to_string(max_neighbours) +
        " cutoff=" + std::to_string(SPATIAL_NEIGHBOUR_CUTOFF_A) + "A");

    return result;
}


std::vector<size_t> SpatialIndexResult::AtomsWithinRadius(
        Vec3 point, double radius) const {
    double query_pt[3] = { point.x(), point.y(), point.z() };
    std::vector<nanoflann::ResultItem<size_t, double>> matches;
    atom_tree_->radiusSearch(query_pt, radius * radius, matches);

    std::vector<size_t> indices;
    indices.reserve(matches.size());
    for (const auto& m : matches) indices.push_back(m.first);
    return indices;
}


std::vector<size_t> SpatialIndexResult::KNearestAtoms(
        Vec3 point, size_t k) const {
    double query_pt[3] = { point.x(), point.y(), point.z() };
    std::vector<size_t> indices(k);
    std::vector<double> dists(k);
    size_t found = atom_tree_->knnSearch(query_pt, k, indices.data(), dists.data());
    indices.resize(found);
    return indices;
}


std::vector<size_t> SpatialIndexResult::RingsWithinRadius(
        Vec3 point, double radius) const {
    if (!ring_tree_) return {};
    double query_pt[3] = { point.x(), point.y(), point.z() };
    std::vector<nanoflann::ResultItem<size_t, double>> matches;
    ring_tree_->radiusSearch(query_pt, radius * radius, matches);

    std::vector<size_t> indices;
    indices.reserve(matches.size());
    for (const auto& m : matches) indices.push_back(m.first);
    return indices;
}


std::vector<size_t> SpatialIndexResult::BondsWithinRadius(
        Vec3 point, double radius) const {
    if (!bond_tree_) return {};
    double query_pt[3] = { point.x(), point.y(), point.z() };
    std::vector<nanoflann::ResultItem<size_t, double>> matches;
    bond_tree_->radiusSearch(query_pt, radius * radius, matches);

    std::vector<size_t> indices;
    indices.reserve(matches.size());
    for (const auto& m : matches) indices.push_back(m.first);
    return indices;
}

}  // namespace nmr
