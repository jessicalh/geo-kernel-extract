#include "GeometryResult.h"
#include <limits>
#include <cmath>

namespace nmr {

std::unique_ptr<GeometryResult> GeometryResult::Compute(ProteinConformation& conf) {
    auto result = std::make_unique<GeometryResult>();
    result->conf_ = &conf;

    const Protein& protein = conf.ProteinRef();
    const auto& positions = conf.Positions();

    // ---------------------------------------------------------------
    // Ring geometry: center (centroid), normal (SVD), radius, vertices
    // ---------------------------------------------------------------
    conf.ring_geometries.resize(protein.RingCount());
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        conf.ring_geometries[ri] = protein.RingAt(ri).ComputeGeometry(positions);
    }

    // ---------------------------------------------------------------
    // Bond geometry: length, direction, midpoint
    // ---------------------------------------------------------------
    conf.bond_lengths.resize(protein.BondCount());
    conf.bond_directions.resize(protein.BondCount());
    conf.bond_midpoints.resize(protein.BondCount());
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        conf.bond_lengths[bi] = bond.Length(positions);
        conf.bond_directions[bi] = bond.Direction(positions);
        conf.bond_midpoints[bi] = bond.Midpoint(positions);
    }

    // ---------------------------------------------------------------
    // Global geometry: bounding box, center of geometry, radius of gyration
    // ---------------------------------------------------------------
    if (!positions.empty()) {
        Vec3 bmin = positions[0];
        Vec3 bmax = positions[0];
        Vec3 center = Vec3::Zero();

        for (const auto& p : positions) {
            bmin = bmin.cwiseMin(p);
            bmax = bmax.cwiseMax(p);
            center += p;
        }
        center /= static_cast<double>(positions.size());

        conf.bounding_min = bmin;
        conf.bounding_max = bmax;
        conf.center_of_geometry = center;

        double rg_sum = 0.0;
        for (const auto& p : positions)
            rg_sum += (p - center).squaredNorm();
        conf.radius_of_gyration = std::sqrt(rg_sum / positions.size());
    }

    // ---------------------------------------------------------------
    // Pre-built collections: rings by type, bonds by category, residues by type
    // ---------------------------------------------------------------
    for (size_t ri = 0; ri < protein.RingCount(); ++ri)
        conf.rings_by_type[protein.RingAt(ri).type_index].push_back(ri);

    for (size_t bi = 0; bi < protein.BondCount(); ++bi)
        conf.bonds_by_category[protein.BondAt(bi).category].push_back(bi);

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri)
        conf.residues_by_type[protein.ResidueAt(ri).type].push_back(ri);

    // ---------------------------------------------------------------
    // Ring pair properties
    // ---------------------------------------------------------------
    for (size_t i = 0; i < protein.RingCount(); ++i) {
        for (size_t j = i + 1; j < protein.RingCount(); ++j) {
            const auto& gi = conf.ring_geometries[i];
            const auto& gj = conf.ring_geometries[j];

            ProteinConformation::RingPair pair;
            pair.ring_a = i;
            pair.ring_b = j;
            pair.center_distance = (gi.center - gj.center).norm();
            pair.normal_dot = gi.normal.dot(gj.normal);
            pair.normal_cross_mag = gi.normal.cross(gj.normal).norm();
            pair.is_fused = (protein.RingAt(i).fused_partner_index == j);

            conf.ring_pairs.push_back(pair);
        }
    }

    return result;
}


const RingGeometry& GeometryResult::RingGeometryAt(size_t ring_index) const {
    return conf_->ring_geometries[ring_index];
}

double GeometryResult::BondLengthAt(size_t bond_index) const {
    return conf_->bond_lengths[bond_index];
}

Vec3 GeometryResult::BondMidpointAt(size_t bond_index) const {
    return conf_->bond_midpoints[bond_index];
}

Vec3 GeometryResult::BondDirectionAt(size_t bond_index) const {
    return conf_->bond_directions[bond_index];
}

}  // namespace nmr
