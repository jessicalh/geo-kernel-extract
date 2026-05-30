#include "CameraInputFilter.h"

#include "MoleculeScene.h"

#include "../diagnostics/ObjectCensus.h"

#include <QEvent>
#include <QLoggingCategory>
#include <QMouseEvent>
#include <QVTKOpenGLNativeWidget.h>
#include <QWheelEvent>

#include <cmath>

namespace h5reader::app {

namespace {
// Gesture sensitivity. Chosen to match the stock vtkInteractorStyleTrackballCamera
// feel: 200 pixels of mouse drag at the default widget size sweeps roughly
// 90 degrees of rotation. Empirical — tuned alongside the harness.
constexpr double kRotateRadiansPerPixel = 0.005;
constexpr double kDollyPerWheelTick     = 1.0 / 1200.0;  // 1 wheel notch = 120 angleDelta
constexpr double kDollyPerDragPixel     = 1.0 / 200.0;
}  // namespace

CameraInputFilter::CameraInputFilter(QVTKOpenGLNativeWidget* widget,
                                      MoleculeScene*          scene,
                                      CameraComposer*         composer,
                                      QObject*                parent)
    : QObject(parent),
      widget_(widget),
      scene_(scene),
      composer_(composer) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("CameraInputFilter"));
    if (widget_) widget_->installEventFilter(this);
}

CameraInputFilter::~CameraInputFilter() {
    if (widget_) widget_->removeEventFilter(this);
}

bool CameraInputFilter::eventFilter(QObject* obj, QEvent* event) {
    if (obj != widget_.data())
        return QObject::eventFilter(obj, event);
    if (!composer_ || !scene_)
        return QObject::eventFilter(obj, event);

    switch (event->type()) {
        case QEvent::MouseButtonPress: {
            handleMouseDown(static_cast<QMouseEvent*>(event));
            // Consume; the trackball is a QuietTrackballStyle so it's a
            // no-op anyway, but this also stops the picker (filter
            // installed before us) from misinterpreting single clicks.
            return activeGesture_ != Gesture::None;
        }
        case QEvent::MouseMove: {
            if (activeGesture_ == Gesture::None)
                return false;
            handleMouseMove(static_cast<QMouseEvent*>(event));
            return true;
        }
        case QEvent::MouseButtonRelease: {
            if (activeGesture_ == Gesture::None)
                return false;
            handleMouseUp(static_cast<QMouseEvent*>(event));
            return true;
        }
        case QEvent::Wheel: {
            handleWheel(static_cast<QWheelEvent*>(event));
            return true;
        }
        default:
            return QObject::eventFilter(obj, event);
    }
}

void CameraInputFilter::handleMouseDown(QMouseEvent* me) {
    ASSERT_THREAD(this);
    if (!me) return;
    lastPos_ = me->position();
    const auto mods = me->modifiers();
    const auto buttons = me->button();
    if (buttons == Qt::LeftButton) {
        activeGesture_ = (mods & Qt::ShiftModifier) ? Gesture::Pan : Gesture::Rotate;
    } else if (buttons == Qt::MiddleButton) {
        activeGesture_ = Gesture::Pan;
    } else if (buttons == Qt::RightButton) {
        activeGesture_ = Gesture::Dolly;
    } else {
        activeGesture_ = Gesture::None;
    }
}

void CameraInputFilter::handleMouseMove(QMouseEvent* me) {
    ASSERT_THREAD(this);
    if (!me || !composer_) return;
    const QPointF p = me->position();
    const double dx = p.x() - lastPos_.x();
    const double dy = p.y() - lastPos_.y();
    lastPos_ = p;

    switch (activeGesture_) {
        case Gesture::Rotate: {
            CameraGesture g;
            g.kind = CameraGesture::Kind::Azimuth;
            g.deltaRadians = -dx * kRotateRadiansPerPixel;
            composer_->applyGesture(g);
            g.kind = CameraGesture::Kind::Elevation;
            g.deltaRadians = -dy * kRotateRadiansPerPixel;
            composer_->applyGesture(g);
            break;
        }
        case Gesture::Pan: {
            CameraGesture g;
            g.kind = CameraGesture::Kind::Pan;
            g.dxScreenPx = dx;
            g.dyScreenPx = dy;
            composer_->applyGesture(g);
            break;
        }
        case Gesture::Dolly: {
            CameraGesture g;
            g.kind = CameraGesture::Kind::Dolly;
            const double scale = 1.0 + dy * kDollyPerDragPixel;
            g.dollyFactor = scale > 1e-3 ? scale : 1e-3;
            composer_->applyGesture(g);
            break;
        }
        case Gesture::None:
            break;
    }
    if (scene_) scene_->requestRender(MoleculeScene::RenderSource::CameraInput);
}

void CameraInputFilter::handleMouseUp(QMouseEvent* /*me*/) {
    ASSERT_THREAD(this);
    activeGesture_ = Gesture::None;
}

void CameraInputFilter::handleWheel(QWheelEvent* we) {
    ASSERT_THREAD(this);
    if (!we || !composer_) return;
    const int angle = we->angleDelta().y();
    if (angle == 0) return;
    CameraGesture g;
    g.kind = CameraGesture::Kind::Dolly;
    g.dollyFactor = std::exp(angle * kDollyPerWheelTick);
    composer_->applyGesture(g);
    if (scene_) scene_->requestRender(MoleculeScene::RenderSource::CameraInput);
}

}  // namespace h5reader::app
