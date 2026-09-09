// CameraInputFilter — Qt eventFilter that intercepts mouse + wheel
// events on the QVTKOpenGLNativeWidget and routes them through
// CameraComposer::applyGesture instead of letting VTK's trackball
// process them.
//
// QtAtomPicker installs its event filter in the same way:
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
// Touch and 3D-mouse input are not handled here.

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

signals:
    // A left-button press+release with no meaningful drag — a plain click, not
    // a rotate. pos is in widget coordinates. ReaderMainWindow uses it to toggle
    // playback when the click misses every atom.
    void viewportClicked(QPointF posInWidget);

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

    // Click-vs-drag discrimination for viewportClicked.
    QPointF         pressPos_;
    bool            moved_       = false;
    Qt::MouseButton pressButton_ = Qt::NoButton;
};

}  // namespace h5reader::app
