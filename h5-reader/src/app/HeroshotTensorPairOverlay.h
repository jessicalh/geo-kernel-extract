#pragma once

#include "../model/Types.h"

#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <cstddef>
#include <vector>

namespace h5reader::app {

// Resthero-only comparison of two symmetric-traceless tensors at one physical
// site. The candidate is a translucent signed surface; the reference is the
// same surface in wireframe. Both use one fixed Angstrom-per-ppm scale.
class HeroshotTensorPairOverlay {
public:
    struct Sample {
        model::Vec3 center = model::Vec3::Zero();
        model::Mat3 candidate = model::Mat3::Zero();
        model::Mat3 reference = model::Mat3::Zero();
        model::Vec3 referenceNormal = model::Vec3::Zero();
    };

    struct Style {
        double scaleAperPpm = 0.42;
        double minimumRadiusA = 0.015;
        double candidateOpacity = 0.44;
        double referenceOpacity = 0.95;
        double referenceLineWidth = 2.2;
        int thetaResolution = 96;
        int phiResolution = 96;
        bool showCandidate = true;
        bool showReference = true;
        bool showSurfaces = true;
        bool showDirectors = false;
        bool showReferenceNormal = false;
        double candidateDirectorHalfLengthA = 1.2;
        double referenceDirectorHalfLengthA = 1.344;
        double referenceNormalHalfLengthA = 1.488;
        double directorRadiusA = 0.025;
    };

    explicit HeroshotTensorPairOverlay(vtkSmartPointer<vtkRenderer> renderer);
    ~HeroshotTensorPairOverlay();

    HeroshotTensorPairOverlay(const HeroshotTensorPairOverlay&) = delete;
    HeroshotTensorPairOverlay& operator=(const HeroshotTensorPairOverlay&) = delete;

    bool show(const std::vector<Sample>& samples, const Style& style);
    void clear();
    std::size_t size() const { return sampleCount_; }

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkActor>> actors_;
    std::size_t sampleCount_ = 0;
};

}  // namespace h5reader::app
