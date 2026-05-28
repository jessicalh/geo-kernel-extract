// SelectionDock — the QListView panel bound to the AtomSelection model.
//
// The selection's first VIEW: it makes the QAbstractListModel real (a model
// with no view never runs its data()/roleNames()) and gives the ≤4-atom set a
// 2-D readout — slot colour swatch (DecorationRole) + "TRP42:CA" label
// (DisplayRole) — alongside the 3-D spheres. A header line reports the count
// and the geometry the set defines (Distance / Angle / Dihedral).
//
// Increment 1: display-only (no row selection). Later increments make a row
// the focus on click and host the comparison/value columns.

#pragma once

#include <QDockWidget>
#include <QPointer>

class QLabel;
class QListView;

namespace h5reader::model {
class AtomSelection;
}

namespace h5reader::app {

class SelectionDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit SelectionDock(QWidget* parent = nullptr);
    ~SelectionDock() override = default;

    // Bind the model. The view shows its rows; the header tracks its
    // changed()/cleared() signals.
    void setModel(model::AtomSelection* selection);

private slots:
    // Update the count + geometry-kind caption.
    void refreshHeader();

private:
    QPointer<QLabel>               header_;
    QPointer<QLabel>               detail_;
    QPointer<QListView>            list_;
    QPointer<model::AtomSelection> selection_;
};

}  // namespace h5reader::app
