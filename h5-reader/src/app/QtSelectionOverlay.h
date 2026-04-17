// QtSelectionOverlay — yellow translucent sphere at the picked atom.
//
// Listens for atomPicked(idx) from QtAtomPicker and for frameChanged(t)
// from QtPlaybackController. On either signal, moves the sphere to
// the picked atom's current-frame position.
//
// Sphere: 1.0 Å radius, yellow, 0.3 opacity — matches the library
// viewer's selection highlight.

#pragma once

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>

#include <cstddef>

namespace h5reader::app {

class QtSelectionOverlay final : public QObject {
    Q_OBJECT
public:
    explicit QtSelectionOverlay(
        vtkSmartPointer<vtkRenderer>                  renderer,
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
        QObject* parent = nullptr);
    ~QtSelectionOverlay() override;

    void Build(const model::QtProtein&      protein,
               const model::QtConformation& conformation);

public slots:
    // Move the sphere to this atom and show it. Called from
    // QtAtomPicker::atomPicked.
    void setPickedAtom(std::size_t atomIdx);

    // Clear the selection (hide the sphere).
    void clearSelection();

    // Called by MoleculeScene when the frame advances. Moves the sphere
    // to the picked atom's new position, or no-ops if no atom is picked.
    void setFrame(int t);

private:
    void applyCurrentPosition(int t);

    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    vtkSmartPointer<vtkSphereSource>              sphere_;
    vtkSmartPointer<vtkActor>                     actor_;

    const model::QtProtein*      protein_      = nullptr;
    const model::QtConformation* conformation_ = nullptr;

    bool        hasSelection_ = false;
    std::size_t pickedAtom_   = 0;
};

}  // namespace h5reader::app
