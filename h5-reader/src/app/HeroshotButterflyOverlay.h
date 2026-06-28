// HeroshotButterflyOverlay -- high-resolution transient ring-current surfaces.
//
// This is deliberately separate from QtFieldGridOverlay. The normal viewer
// butterfly remains the playback/UI overlay. This class is a resthero/export
// layer: sample the same closed-form field equations more densely for one
// selected ring, contour the same signed T0 isovalues, and clear it afterward.

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
}

namespace h5reader::app {

class HeroshotButterflyOverlay {
public:
    enum class Mode {
        BiotSavart,
        HaighMallion,
        Sum,
    };

    struct Style {
        int gridDim = 40;
        double thresholdPpm = 0.75;
        double gaussianExtentA = 7.0;
        double gaussianPeak = 1.0;
        double opacity = 0.24;
        bool showShielded = true;
        bool showDeshielded = true;
        Mode mode = Mode::BiotSavart;
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

    bool show(const model::QtProtein& protein,
              model::Conformation& conformation,
              std::size_t ring,
              std::size_t frame,
              const Style& style = Style{});
    void clear();
    std::size_t size() const { return pipelines_.size(); }
    const Stats& stats() const { return stats_; }

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
    };

    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<Pipeline> pipelines_;
    Stats stats_;
};

}  // namespace h5reader::app
