// QtAtomInspectorDock — tabified dock displaying the full typed state
// of the picked atom at the current frame.
//
// Tree structure:
//   Identity                (element, role, residue, chain, flags — static)
//   Position                (per frame)
//   Ring current            (bs/hm/rs shielding, proximity counts, B-field)
//   Bond anisotropy         (mc shielding, nearest C=O / C-N stats)
//   Quadrupole / Dispersion (pq, disp shielding)
//   Electrostatics          (coulomb, apbs, aimnet2 shielding + E-field)
//   H-bond                  (nearest partner, counts, donor/acceptor flags)
//   SASA                    (Å² + outward normal)
//   Water environment       (E-field + shell counts + dipole cos)
//   Charges                 (AIMNet2, EEQ, CN)
//
// Updates on atomPicked(idx) from QtAtomPicker AND on frameChanged(t)
// from QtPlaybackController.

#pragma once

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QDockWidget>
#include <QPointer>

#include <cstddef>

class QTreeWidget;
class QTreeWidgetItem;

namespace h5reader::app {

class QtAtomInspectorDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit QtAtomInspectorDock(QWidget* parent = nullptr);
    ~QtAtomInspectorDock() override = default;

    // Bind the typed model. Call once after H5 load, before any
    // setPickedAtom / setFrame.
    void setContext(const model::QtProtein*      protein,
                    const model::QtConformation* conformation);

public slots:
    // The dock's two inputs: which atom and which frame. Both cause a
    // full tree rebuild because per-frame values affect many fields.
    void setPickedAtom(std::size_t atomIdx);
    void setFrame(int t);

    // Clear the tree (e.g. load unmounted or picker cleared).
    void clearSelection();

private:
    void rebuild();
    void populateIdentity(QTreeWidgetItem* parent);
    void populatePerFrame(QTreeWidgetItem* root);

    QPointer<QTreeWidget>        tree_;
    const model::QtProtein*      protein_      = nullptr;
    const model::QtConformation* conformation_ = nullptr;
    bool                         hasSelection_ = false;
    std::size_t                  atomIdx_      = 0;
    int                          frame_        = 0;
};

}  // namespace h5reader::app
