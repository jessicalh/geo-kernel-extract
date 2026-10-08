// Double-click event handling. MoleculeScene picks the rendered atom;
// AtomSelection interprets the keyboard modifiers.

#pragma once

#include <QObject>
#include <QPointer>

#include <cstddef>

class QVTKOpenGLNativeWidget;

namespace h5reader::app {

class MoleculeScene;

class QtAtomPicker final : public QObject {
    Q_OBJECT
public:
    QtAtomPicker(QVTKOpenGLNativeWidget* vtkWidget,
                 MoleculeScene* scene, QObject* parent = nullptr);
    ~QtAtomPicker() override;

signals:
    // atomIdx is the protein atom index, including when the scene is filtered.
    void atomPicked(std::size_t atomIdx, Qt::KeyboardModifiers modifiers);

protected:
    bool eventFilter(QObject* obj, QEvent* event) override;

private:
    QPointer<QVTKOpenGLNativeWidget> vtkWidget_;
    QPointer<MoleculeScene>          scene_;
};

}  // namespace h5reader::app
