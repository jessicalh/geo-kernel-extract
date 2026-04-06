#include "Ring.h"
#include <Eigen/SVD>

namespace nmr {

RingGeometry Ring::ComputeGeometry(const std::vector<Vec3>& positions) const {
    RingGeometry geo;
    if (atom_indices.empty()) return geo;

    // Collect vertex positions
    geo.vertices.reserve(atom_indices.size());
    for (size_t idx : atom_indices)
        geo.vertices.push_back(positions[idx]);

    // Center: centroid of vertex positions
    geo.center = Vec3::Zero();
    for (const auto& v : geo.vertices)
        geo.center += v;
    geo.center /= static_cast<double>(geo.vertices.size());

    // Normal: best-fit plane via SVD of centered vertex coordinates.
    // The smallest singular value's corresponding right singular vector
    // is the normal to the plane of best fit.
    if (geo.vertices.size() >= 3) {
        Eigen::MatrixXd coords(geo.vertices.size(), 3);
        for (size_t i = 0; i < geo.vertices.size(); ++i)
            coords.row(i) = (geo.vertices[i] - geo.center).transpose();

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(coords, Eigen::ComputeFullV);
        geo.normal = svd.matrixV().col(2);

        // Consistent orientation: normal in same direction as cross product
        // of first two edges (right-hand rule).
        Vec3 edge01 = geo.vertices[1] - geo.vertices[0];
        Vec3 edge02 = geo.vertices[2] - geo.vertices[0];
        if (geo.normal.dot(edge01.cross(edge02)) < 0)
            geo.normal = -geo.normal;
    }

    // Radius: mean distance from center to vertices
    geo.radius = 0.0;
    for (const auto& v : geo.vertices)
        geo.radius += (v - geo.center).norm();
    geo.radius /= static_cast<double>(geo.vertices.size());

    return geo;
}


std::unique_ptr<Ring> CreateRing(RingTypeIndex type) {
    switch (type) {
        case RingTypeIndex::PheBenzene:   return std::make_unique<PheBenzeneRing>();
        case RingTypeIndex::TyrPhenol:    return std::make_unique<TyrPhenolRing>();
        case RingTypeIndex::TrpBenzene:   return std::make_unique<TrpBenzeneRing>();
        case RingTypeIndex::TrpPyrrole:   return std::make_unique<TrpPyrroleRing>();
        case RingTypeIndex::TrpPerimeter: return std::make_unique<IndolePerimeterRing>();
        case RingTypeIndex::HisImidazole: return std::make_unique<HisImidazoleRing>();
        case RingTypeIndex::HidImidazole: return std::make_unique<HidImidazoleRing>();
        case RingTypeIndex::HieImidazole: return std::make_unique<HieImidazoleRing>();
        default: return std::make_unique<PheBenzeneRing>();
    }
}

}  // namespace nmr
