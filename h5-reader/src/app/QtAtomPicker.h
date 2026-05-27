// QtAtomPicker — ray-cast double-click picker for the VTK viewport.
//
// Installs an event filter on the QVTKOpenGLNativeWidget, catches
// double-clicks, converts the click point to a world-space ray via
// the renderer's camera, loops protein atoms, picks the one nearest
// to the ray (with a tolerance threshold), and emits atomPicked(idx).
//
// Pattern ported from ui/src/MainWindow.cpp::pickAtom (library viewer).
// Uses screen-space → world projection + nearest-to-ray; vtkCellPicker
// and friends don't play well with vtkOpenGLMoleculeMapper's GPU
// imposter rendering, so we pick against atom positions directly.
//
// Hi-DPI correct via QWidget::devicePixelRatioF(). Y-axis flipped
// between Qt (top-left origin) and VTK (bottom-left).

#pragma once

#include <QObject>
#include <QPointer>

#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <cstddef>

class QVTKOpenGLNativeWidget;

namespace h5reader::model {
class Conformation;
class QtProtein;
}

namespace h5reader::app {

class QtPlaybackController;

class QtAtomPicker final : public QObject {
    Q_OBJECT
public:
    QtAtomPicker(QVTKOpenGLNativeWidget*                vtkWidget,
                 vtkSmartPointer<vtkRenderer>           renderer,
                 const model::QtProtein*                 protein,
                 model::Conformation*                    conformation,
                 const QtPlaybackController*             playback,
                 QObject*                                parent = nullptr);
    ~QtAtomPicker() override;

signals:
    // Emitted when the user double-clicks over an atom (within the
    // 2 Å ray-distance tolerance). atomIdx is an index into QtProtein;
    // modifiers carries the keyboard state at click time (Shift =
    // add-to-selection) for AtomSelection to interpret — the picker stays
    // dumb and does not act on it. If no atom is close enough, none fires.
    void atomPicked(std::size_t atomIdx, Qt::KeyboardModifiers modifiers);

protected:
    bool eventFilter(QObject* obj, QEvent* event) override;

private:
    void doPick(int displayX, int displayY, Qt::KeyboardModifiers mods);

    QPointer<QVTKOpenGLNativeWidget> vtkWidget_;
    vtkSmartPointer<vtkRenderer>     renderer_;
    const model::QtProtein*          protein_;
    QPointer<model::Conformation>    conformation_;
    const QtPlaybackController*      playback_;
};

}  // namespace h5reader::app
