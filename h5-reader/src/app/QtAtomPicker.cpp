#include "QtAtomPicker.h"

#include "QtPlaybackController.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/Conformation.h"
#include "../model/QtProtein.h"

#include <QEvent>
#include <QLoggingCategory>
#include <QMouseEvent>

#include <QVTKOpenGLNativeWidget.h>

#include <vtkCamera.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

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
                            model::Conformation*                    conformation,
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
        doPick(me->position().x(), me->position().y(), me->modifiers());
        return true;
    }
    return QObject::eventFilter(obj, event);
}

std::optional<std::size_t> QtAtomPicker::atomAt(int displayX, int displayY) const {
    if (!protein_ || !conformation_ || !renderer_ || !vtkWidget_)
        return std::nullopt;

    // Convert Qt widget coords → VTK widget coords. Qt origin is top-
    // left; VTK is bottom-left. Device pixel ratio for Hi-DPI displays.
    const double dpr = vtkWidget_->devicePixelRatioF();
    const int    vtkX = static_cast<int>(displayX * dpr);
    const int    vtkY = static_cast<int>((vtkWidget_->height() - displayY) * dpr);

    auto* camera = renderer_->GetActiveCamera();
    double camPos[3]; camera->GetPosition(camPos);
    const model::Vec3 rayOrigin(camPos[0], camPos[1], camPos[2]);

    renderer_->SetDisplayPoint(vtkX, vtkY, 0.0);
    renderer_->DisplayToWorld();
    double worldPt[4]; renderer_->GetWorldPoint(worldPt);
    const double w = worldPt[3];
    if (std::abs(w) < 1e-12) return std::nullopt;
    const model::Vec3 clickWorld(worldPt[0] / w, worldPt[1] / w, worldPt[2] / w);
    const model::Vec3 rayDir = (clickWorld - rayOrigin).normalized();

    // Walk all atoms at the current frame; take the one whose closest
    // approach to the ray is smallest, provided it's in front of the
    // camera (projLen >= 0) and within the pick tolerance.
    const int t = playback_ ? playback_->currentFrame() : 0;
    const size_t st = static_cast<size_t>(std::max(0, t));

    double bestDist = std::numeric_limits<double>::infinity();
    size_t bestAtom = 0;
    bool   found    = false;

    const size_t N = protein_->atomCount();
    for (size_t i = 0; i < N; ++i) {
        const model::Vec3 pos = conformation_->atomPosition(st, i);
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

    if (!found || bestDist > kMaxPickDistanceA)
        return std::nullopt;
    return bestAtom;
}

void QtAtomPicker::doPick(int displayX, int displayY, Qt::KeyboardModifiers mods) {
    ASSERT_THREAD(this);
    const auto hit = atomAt(displayX, displayY);
    if (!hit) {
        qCDebug(cPicker).noquote() << "no pick";
        return;
    }
    qCInfo(cPicker).noquote()
        << "atom" << *hit
        << "| frame=" << (playback_ ? playback_->currentFrame() : 0);
    emit atomPicked(*hit, mods);
}

}  // namespace h5reader::app
