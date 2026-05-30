// CameraInputFilter — Qt eventFilter that intercepts mouse + wheel
// events on the QVTKOpenGLNativeWidget and routes them through
// CameraComposer::applyGesture instead of letting VTK's trackball
// process them.
//
// Pattern from QtAtomPicker (QtAtomPicker.cpp:48 installs an event
// filter the same way). Per spec/viewport_pipeline_2026-05-30.md §4.1:
//   * MouseButtonPress (Left/Middle/Right) starts a gesture
//   * MouseMove drives the active gesture, calls
//     CameraComposer::applyGesture, asks MoleculeScene to render
//   * MouseButtonRelease ends the gesture
//   * Wheel events dolly
//   * MouseButtonDblClick is NOT intercepted (the picker owns it)
//
// Install AFTER the picker so Qt's filter chain calls THIS filter
// first; double-click events fall through to the picker.
//
// Gesture mapping (matches stock trackball conventions in this codebase):
//   Left drag       — rotate (azimuth + elevation)
//   Middle drag     — pan
//   Right drag      — dolly (vertical = zoom)
//   Shift+Left drag — pan
//   Wheel           — dolly
//
// No touch / 3D-mouse arms per the implementation prompt §1 (decision:
// no touch, no SpaceMouse). Future input devices add their arms here
// without changing the rest of the pipeline.

#pragma once

#include "CameraComposer.h"

#include "../diagnostics/ThreadGuard.h"

#include <QObject>
#include <QPoint>
#include <QPointF>
#include <QPointer>

class QMouseEvent;
class QVTKOpenGLNativeWidget;
class QWheelEvent;

namespace h5reader::app {

class MoleculeScene;

class CameraInputFilter final : public QObject {
    Q_OBJECT
public:
    CameraInputFilter(QVTKOpenGLNativeWidget* widget,
                       MoleculeScene*          scene,
                       CameraComposer*         composer,
                       QObject*                parent = nullptr);
    ~CameraInputFilter() override;

protected:
    bool eventFilter(QObject* obj, QEvent* event) override;

private:
    enum class Gesture { None, Rotate, Pan, Dolly };

    void handleMouseDown(QMouseEvent* me);
    void handleMouseMove(QMouseEvent* me);
    void handleMouseUp(QMouseEvent* me);
    void handleWheel(QWheelEvent* we);

    QPointer<QVTKOpenGLNativeWidget> widget_;
    QPointer<MoleculeScene>          scene_;
    QPointer<CameraComposer>         composer_;

    Gesture activeGesture_ = Gesture::None;
    QPointF lastPos_;
};

}  // namespace h5reader::app
