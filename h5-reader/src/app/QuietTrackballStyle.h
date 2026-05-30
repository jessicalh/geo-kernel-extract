// QuietTrackballStyle — minimal subclass of vtkInteractorStyleTrackballCamera
// that no-ops every camera manipulator.
//
// Per spec/viewport_pipeline_2026-05-30.md §2.7 / §4.2: the adapter's paint
// chain (QVTKRenderWindowAdapter::paint at .cxx:241-266) calls iren->Render()
// and has no fallback when the interactor is present, so we must keep
// iren->EnableRender ON or paintGL never actually draws. We can't switch
// off the trackball that way; instead this subclass overrides every
// manipulator (Rotate/Spin/Pan/Dolly/EnvironmentRotate/OnMouseWheel) to a
// no-op so VTK's stock trackball never mutates the camera or schedules
// renders from event paths the Qt eventFilter doesn't intercept (3D-mouse,
// touch, future extensions).
//
// All actual camera input flows through CameraInputFilter installed on
// QVTKOpenGLNativeWidget (Stage X in the spec). The trackball still gets
// called by VTK so the paint chain works, but does nothing on its own
// initiative.

#pragma once

#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>

namespace h5reader::app {

class QuietTrackballStyle : public vtkInteractorStyleTrackballCamera {
public:
    static QuietTrackballStyle* New();
    vtkTypeMacro(QuietTrackballStyle, vtkInteractorStyleTrackballCamera)

    void Rotate() override {}
    void Spin() override {}
    void Pan() override {}
    void Dolly() override {}
    void EnvironmentRotate() override {}
    void OnMouseWheelForward() override {}
    void OnMouseWheelBackward() override {}

protected:
    QuietTrackballStyle() = default;
    ~QuietTrackballStyle() override = default;
};

}  // namespace h5reader::app
