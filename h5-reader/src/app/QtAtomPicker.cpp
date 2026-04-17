#include "QtAtomPicker.h"

#include "QtPlaybackController.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/QtConformation.h"
#include "../model/QtFrame.h"
#include "../model/QtProtein.h"

#include <QEvent>
#include <QLoggingCategory>
#include <QMouseEvent>

#include <QVTKOpenGLNativeWidget.h>

#include <vtkCamera.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cPicker, "h5reader.picker")

// Maximum ray-to-atom distance for a pick to register, in Angstroms.
// Matches the library viewer's tolerance.
constexpr double kMaxPickDistanceA = 2.0;
}

QtAtomPicker::QtAtomPicker(QVTKOpenGLNativeWidget*                vtkWidget,
                            vtkSmartPointer<vtkRenderer>           renderer,
                            const model::QtProtein*                 protein,
                            const model::QtConformation*            conformation,
                            const QtPlaybackController*             playback,
                            QObject*                                parent)
    : QObject(parent),
      vtkWidget_(vtkWidget),
      renderer_(std::move(renderer)),
      protein_(protein),
      conformation_(conformation),
      playback_(playback)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtAtomPicker"));
    if (vtkWidget_) vtkWidget_->installEventFilter(this);
}

QtAtomPicker::~QtAtomPicker() {
    if (vtkWidget_) vtkWidget_->removeEventFilter(this);
}

bool QtAtomPicker::eventFilter(QObject* obj, QEvent* event) {
    if (obj == vtkWidget_.data()
        && event->type() == QEvent::MouseButtonDblClick) {
        auto* me = static_cast<QMouseEvent*>(event);
        doPick(me->position().x(), me->position().y());
        return true;
    }
    return QObject::eventFilter(obj, event);
}

void QtAtomPicker::doPick(int displayX, int displayY) {
    ASSERT_THREAD(this);
    if (!protein_ || !conformation_ || !renderer_ || !vtkWidget_) return;

    // Convert Qt widget coords → VTK widget coords. Qt origin is top-
    // left; VTK is bottom-left. Device pixel ratio for Hi-DPI displays.
    const double dpr = vtkWidget_->devicePixelRatioF();
    const int    vtkX = static_cast<int>(displayX * dpr);
    const int    vtkY = static_cast<int>((vtkWidget_->height() - displayY) * dpr);

    qCDebug(cPicker).noquote()
        << "qt=(" << displayX << "," << displayY << ")"
        << "vtk=(" << vtkX << "," << vtkY << ") dpr=" << dpr;

    auto* camera = renderer_->GetActiveCamera();
    double camPos[3]; camera->GetPosition(camPos);
    const model::Vec3 rayOrigin(camPos[0], camPos[1], camPos[2]);

    renderer_->SetDisplayPoint(vtkX, vtkY, 0.0);
    renderer_->DisplayToWorld();
    double worldPt[4]; renderer_->GetWorldPoint(worldPt);
    const double w = worldPt[3];
    if (std::abs(w) < 1e-12) return;
    const model::Vec3 clickWorld(worldPt[0] / w, worldPt[1] / w, worldPt[2] / w);
    const model::Vec3 rayDir = (clickWorld - rayOrigin).normalized();

    // Walk all atoms at the current frame; take the one whose closest
    // approach to the ray is smallest, provided it's in front of the
    // camera (projLen >= 0) and within the pick tolerance.
    const int t = playback_ ? playback_->currentFrame() : 0;
    const auto& frame = conformation_->frame(
        static_cast<size_t>(std::max(0, t)));

    double bestDist = std::numeric_limits<double>::infinity();
    size_t bestAtom = 0;
    bool   found    = false;

    const size_t N = protein_->atomCount();
    for (size_t i = 0; i < N; ++i) {
        const model::Vec3 pos = frame.position(i);
        const model::Vec3 toAtom = pos - rayOrigin;
        const double projLen = toAtom.dot(rayDir);
        if (projLen < 0.0) continue;                  // behind camera
        const model::Vec3 closest = rayOrigin + projLen * rayDir;
        const double d = (pos - closest).norm();
        if (d < bestDist) {
            bestDist = d;
            bestAtom = i;
            found    = true;
        }
    }

    if (!found || bestDist > kMaxPickDistanceA) {
        qCDebug(cPicker).noquote()
            << "no pick | best dist=" << bestDist << "Å";
        return;
    }

    qCInfo(cPicker).noquote()
        << "atom" << bestAtom << "| dist=" << bestDist << "Å"
        << "| frame=" << t;
    emit atomPicked(bestAtom);
}

}  // namespace h5reader::app
