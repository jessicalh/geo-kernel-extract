// AtomTrackOverlay -- transient heroshot signal geometry for one atom path.
//
// Draws sampled positions as signal marks, optionally connected by hairline
// sample edges. The default is screen-space points; heroshots may request tiny
// 3-D dot glyphs. It deliberately avoids atom/bond visual grammar: no tubes, no
// chemistry-sized spheres. The caller supplies the samples and their scalar
// intensity; this class only maps them to VTK geometry. No text, no labels, no
// interpolation claims.

#pragma once

#include "../model/Types.h"

#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <array>
#include <vector>

namespace h5reader::app {

class AtomTrackOverlay {
public:
    enum class ColorMode {
        Signed,
        Absolute,
    };
    enum class PointShape {
        ScreenPoint,
        Sphere,
    };

    struct Style {
        double pointSizePixels = 4.5;
        double sphereRadiusA = 0.03;
        double currentPointScale = 1.8;
        double pointOpacity = 0.82;
        double haloScale = 3.2;
        double haloOpacity = 0.10;
        double lineWidthPixels = 1.0;
        double lineOpacity = 0.24;
        double colorScale = 0.0;
        double colorGamma = 0.65;
        double minColorFraction = 0.10;
        bool showPoints = true;
        bool showHalos = true;
        bool showLines = false;
        bool highlightCurrent = false;
        ColorMode colorMode = ColorMode::Signed;
        PointShape pointShape = PointShape::ScreenPoint;
    };

    struct Sample {
        model::Vec3 position = model::Vec3::Zero();
        double intensity = 0.0;
        bool current = false;
    };

    explicit AtomTrackOverlay(vtkSmartPointer<vtkRenderer> sceneRenderer);
    ~AtomTrackOverlay();

    AtomTrackOverlay(const AtomTrackOverlay&) = delete;
    AtomTrackOverlay& operator=(const AtomTrackOverlay&) = delete;

    void show(const std::vector<Sample>& samples, const Style& style = Style{});
    void clear();
    std::size_t size() const { return actors_.size(); }

private:
    std::array<unsigned char, 3> colorFor(double value,
                                          double scale,
                                          ColorMode mode,
                                          double gamma,
                                          double minFraction) const;
    void addActor(vtkSmartPointer<vtkActor> actor);

    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkActor>> actors_;
};

}  // namespace h5reader::app
