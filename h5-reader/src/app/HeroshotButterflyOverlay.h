// HeroshotButterflyOverlay -- high-resolution transient ring-current surfaces.
//
// This is deliberately separate from QtFieldGridOverlay. The normal viewer
// butterfly remains the playback/UI overlay. This class is a resthero/export
// layer: sample a declared ring or fused-ring field more densely, contour its
// signed T0 isovalues, and clear it afterward.

#pragma once

#include "../model/Types.h"

#include <vtkActor.h>
#include <vtkContourFilter.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#include <cstddef>
#include <vector>

namespace h5reader::model {
class Conformation;
class QtProtein;
}  // namespace h5reader::model

namespace h5reader::app {

class HeroshotButterflyOverlay {
public:
    enum class Mode {
        BiotSavart,
        HaighMallion,
        Sum,
        CircularCandidateA,
    };

    struct Style {
        int gridDim = 40;
        double thresholdPpm = 0.75;
        double gaussianExtentA = 7.0;
        double gaussianPeak = 1.0;
        double opacity = 0.24;
        bool showShielded = true;
        bool showDeshielded = true;
        bool showSourceLoops = false;
        double sourceLoopTubeRadiusA = 0.018;
        double sourceLoopOpacity = 0.88;
        int sourceLoopResolution = 128;
        Mode mode = Mode::BiotSavart;
    };

    struct CircularSource {
        std::size_t ring = 0;
        model::Vec3 center = model::Vec3::Zero();
        model::Vec3 normal = model::Vec3::Zero();
        double planeRmsA = 0.0;
        double radiusA = 0.0;
        double lobeOffsetA = 0.0;
        double currentNanoamperePerTesla = 0.0;
    };

    struct Stats {
        double minT0 = 0.0;
        double maxT0 = 0.0;
        std::size_t shieldedPoints = 0;
        std::size_t deshieldedPoints = 0;
        std::size_t shieldedCells = 0;
        std::size_t deshieldedCells = 0;
    };

    explicit HeroshotButterflyOverlay(vtkSmartPointer<vtkRenderer> renderer);
    ~HeroshotButterflyOverlay();

    HeroshotButterflyOverlay(const HeroshotButterflyOverlay&) = delete;
    HeroshotButterflyOverlay& operator=(const HeroshotButterflyOverlay&) = delete;

    bool show(const model::QtProtein& protein, const model::Conformation& conformation, std::size_t ring, std::size_t frame);
    bool show(const model::QtProtein& protein,
              const model::Conformation& conformation,
              std::size_t ring,
              std::size_t frame,
              const Style& style);
    bool show(const model::QtProtein& protein,
              const model::Conformation& conformation,
              const std::vector<std::size_t>& rings,
              std::size_t frame,
              const Style& style);
    void clear();
    std::size_t size() const { return pipelines_.size(); }
    const Stats& stats() const { return stats_; }
    const std::vector<CircularSource>& circularSources() const { return circularSources_; }
    std::size_t sourceLoopActorCount() const;

private:
    struct Pipeline {
        vtkSmartPointer<vtkImageData> imageData;
        vtkSmartPointer<vtkFloatArray> scalars;
        vtkSmartPointer<vtkTrivialProducer> producer;
        vtkSmartPointer<vtkContourFilter> contourShielded;
        vtkSmartPointer<vtkContourFilter> contourDeshielded;
        vtkSmartPointer<vtkPolyDataMapper> mapperShielded;
        vtkSmartPointer<vtkPolyDataMapper> mapperDeshielded;
        vtkSmartPointer<vtkActor> actorShielded;
        vtkSmartPointer<vtkActor> actorDeshielded;
        std::vector<vtkSmartPointer<vtkActor>> sourceLoopActors;
    };

    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<Pipeline> pipelines_;
    std::vector<CircularSource> circularSources_;
    Stats stats_;
};

}  // namespace h5reader::app
