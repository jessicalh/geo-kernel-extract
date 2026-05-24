#include "DemoResult.h"
#include "Protein.h"
#include "PhysicalConstants.h"
#include <limits>
#include <cmath>

namespace nmr {

std::unique_ptr<DemoResult> DemoResult::Compute(ProteinConformation& conf) {
    auto result = std::make_unique<DemoResult>();
    result->conf_ = &conf;

    // DemoResult depends on GeometryResult for ring geometry.
    // Compute produces the result; AttachResult checks dependencies.
    // But Compute must not crash if called before dependencies are met.
    if (!conf.HasResult<GeometryResult>()) {
        // Return a valid but empty result — AttachResult will reject it
        return result;
    }

    const Protein& protein = conf.ProteinRef();

    // ---------------------------------------------------------------
    // Per-atom: distance and direction to nearest ring center
    // ---------------------------------------------------------------
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Vec3 pos = conf.PositionAt(ai);
        double min_dist = std::numeric_limits<double>::max();
        Vec3 min_dir = Vec3::Zero();

        for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
            const auto& geo = conf.ring_geometries[ri];
            Vec3 delta = geo.center - pos;
            double dist = delta.norm();
            if (dist < min_dist) {
                min_dist = dist;
                min_dir = (dist > 1e-10) ? Vec3(delta / dist) : Vec3::Zero();
            }
        }

        // Store on ConformationAtom
        auto& ca = conf.MutableAtomAt(ai);
        ca.demo_nearest_ring_distance = (min_dist < 1e30) ? min_dist : NO_DATA_SENTINEL;
        ca.demo_nearest_ring_direction = min_dir;
    }

    // ---------------------------------------------------------------
    // SphericalTensor demonstration: decompose a known Mat3
    // ---------------------------------------------------------------
    // Test matrix: asymmetric with known properties
    Mat3 test;
    test << 1.0, 0.5, 0.3,
            0.2, 2.0, 0.7,
            0.1, 0.4, 3.0;

    result->test_decomposition_ = SphericalTensor::Decompose(test);
    result->test_reconstructed_ = result->test_decomposition_.Reconstruct();

    return result;
}


double DemoResult::NearestRingDistance(size_t atom_index) const {
    return conf_->AtomAt(atom_index).demo_nearest_ring_distance;
}

Vec3 DemoResult::NearestRingDirection(size_t atom_index) const {
    return conf_->AtomAt(atom_index).demo_nearest_ring_direction;
}

}  // namespace nmr
