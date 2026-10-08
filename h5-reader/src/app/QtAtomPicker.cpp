#include "QtAtomPicker.h"

#include "MoleculeScene.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QEvent>
#include <QLoggingCategory>
#include <QMouseEvent>

#include <QVTKOpenGLNativeWidget.h>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cPicker, "h5reader.picker")
}

QtAtomPicker::QtAtomPicker(QVTKOpenGLNativeWidget* vtkWidget,
                         MoleculeScene* scene, QObject* parent)
    : QObject(parent), vtkWidget_(vtkWidget), scene_(scene)
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
        ASSERT_THREAD(this);
        auto* mouse = static_cast<QMouseEvent*>(event);
        if (mouse->button() != Qt::LeftButton || !scene_)
            return QObject::eventFilter(obj, event);
        const auto hit = scene_->pickAt(mouse->position());
        if (hit && hit->atom) {
            qCInfo(cPicker) << "atom" << *hit->atom;
            emit atomPicked(*hit->atom, mouse->modifiers());
        }
        return true;
    }
    return QObject::eventFilter(obj, event);
}

}  // namespace h5reader::app
