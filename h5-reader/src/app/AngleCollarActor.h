// AngleCollarActor -- transient resthero geometry for a true dihedral angle.
//
// Draws an Elizabethan-collar style cone around a hinge axis, plus the true
// signed angle arcs and seams on its rim. The caller supplies the real
// axis/reference geometry and the real angle values; this class only turns
// them into transient VTK geometry. Resthero-only, no UI ownership.

#pragma once

#include "../model/Types.h"

#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <array>
#include <vector>

namespace h5reader::app {

class AngleCollarActor {
public:
    struct Style {
        double radius = 1.25;
        double tubeRadius = 0.035;
        double axisPadding = 0.35;
        double coneLength = 0.0;
        double neckRadius = 0.0;
        double rimRadius = 0.0;
        double coneOpacity = 0.26;
        double coneDirection = -1.0;
        int arcSegments = 96;
    };

    struct Arc {
        double angleRadians = 0.0;
        double opacity = 1.0;
        double radiusScale = 1.0;
        std::array<double, 3> color{{1.0, 0.72, 0.18}};
    };

    explicit AngleCollarActor(vtkSmartPointer<vtkRenderer> sceneRenderer);
    ~AngleCollarActor();

    AngleCollarActor(const AngleCollarActor&) = delete;
    AngleCollarActor& operator=(const AngleCollarActor&) = delete;

    bool show(const model::Vec3& axisStart,
              const model::Vec3& axisEnd,
              const model::Vec3& center,
              const model::Vec3& referenceDirection,
              const std::vector<Arc>& arcs);
    bool show(const model::Vec3& axisStart,
              const model::Vec3& axisEnd,
              const model::Vec3& center,
              const model::Vec3& referenceDirection,
              const std::vector<Arc>& arcs,
              const Style& style);
    void clear();

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkActor>> actors_;
};

}  // namespace h5reader::app
