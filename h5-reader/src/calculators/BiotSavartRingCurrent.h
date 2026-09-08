#pragma once

#include "../model/QtRing.h"
#include "../model/Types.h"

#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

namespace h5reader::calculators {

enum class BiotSavartSampleStatus {
    Evaluated,
    OutsideSpatialCutoff,
    WireSingularity,
    InvalidInput,
};

const char* BiotSavartSampleStatusName(BiotSavartSampleStatus status);

struct BiotSavartSegmentContribution {
    std::size_t edge = 0;
    int loopSide = 0;
    model::Vec3 startA = model::Vec3::Zero();
    model::Vec3 endA = model::Vec3::Zero();
    model::Vec3 inducedFieldTeslaPerTesla = model::Vec3::Zero();
};

struct BiotSavartRingSample {
    BiotSavartSampleStatus status = BiotSavartSampleStatus::InvalidInput;
    model::Vec3 pointA = model::Vec3::Zero();
    double ringCurrentIntensityNanoamperesPerTesla = 0.0;
    double distanceToCenterA = std::numeric_limits<double>::quiet_NaN();
    double minimumWireDistanceA = std::numeric_limits<double>::quiet_NaN();
    model::Vec3 inducedFieldTeslaPerTesla = model::Vec3::Zero();
    model::Mat3 shieldingTensorPpm = model::Mat3::Zero();
    model::SphericalTensor spherical;
    std::vector<BiotSavartSegmentContribution> segments;

    bool evaluated() const { return status == BiotSavartSampleStatus::Evaluated; }
};

// One finite-polygon Johnson-Bovey source, built and sampled by the same
// numerical construction as nmr_extract's BiotSavartResult.
class BiotSavartRingCurrent {
public:
    static std::optional<BiotSavartRingCurrent>
    Build(const std::vector<model::Vec3>& orderedVerticesA, double lobeOffsetA);

    // The intensity is induced ring current per applied field. Passing 1 nA/T
    // reproduces nmr_extract's unit-current kernel numerically; a literature
    // intensity produces the corresponding shielding response in ppm.
    BiotSavartRingSample evaluate(
        const model::Vec3& pointA,
        double ringCurrentIntensityNanoamperesPerTesla,
        bool includeSegments = false) const;

    const std::vector<model::Vec3>& verticesA() const { return verticesA_; }
    const model::RingGeometry& geometry() const { return geometry_; }
    double lobeOffsetA() const { return lobeOffsetA_; }

    std::vector<model::Vec3> upperLoopVerticesA() const;
    std::vector<model::Vec3> lowerLoopVerticesA() const;

private:
    BiotSavartRingCurrent(std::vector<model::Vec3> verticesA,
                          model::RingGeometry geometry,
                          double lobeOffsetA);

    std::vector<model::Vec3> verticesA_;
    model::RingGeometry geometry_;
    double lobeOffsetA_ = 0.0;
};

}  // namespace h5reader::calculators
