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

#include "../model/Conformation.h"
#include "../model/QtProtein.h"
#include "../model/TrajectoryFieldAvailability.h"

#include <QDockWidget>
#include <QPointer>
#include <QString>

#include <cstddef>
#include <memory>

class QTreeWidget;
class QTreeWidgetItem;

namespace h5reader::app {

// Light, Qt-only carrier for the focused atom's DFT CSA tensor shape so this
// header stays free of the heavy CsaProbe / DftShieldingStore graph.
// ReaderMainWindow fills it from the AtomCsaResult it already computes for the
// glyph; the panel shows it as the "CSA shielding tensor (DFT)" section. The
// per-axis colours match CsaTensorOverlay's arrows (amber/teal/violet) so the
// in-scene colour-coded arrows stay decodable without in-scene labels.
struct CsaTensorInfo {
    bool    framed = false;
    QString frameKind;        // human label; "unframed (raw PAS)" when not framed
    double  sigmaIso = 0.0;   // absolute shielding (ppm), NOT chemical shift
    double  span = 0.0;       // sigma33 - sigma11 (ppm)
    double  skew = 0.0;       // 3 (sigma22 - iso) / span
    double  eta = 0.0;        // asymmetry [0,1]
    double  sigma11 = 0.0;    // ascending principal values (ppm)
    double  sigma22 = 0.0;
    double  sigma33 = 0.0;
};

class QtAtomInspectorDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit QtAtomInspectorDock(QWidget* parent = nullptr);
    ~QtAtomInspectorDock() override = default;

    // Bind the typed model. Call once after H5 load, before any
    // setPickedAtom / setFrame.
    void setContext(const model::QtProtein* protein,
                    model::Conformation*    conformation);
    void setFieldAvailability(std::shared_ptr<const model::TrajectoryFieldAvailability> availability);

    // The focused atom's DFT CSA tensor for the current frame, mirrored from the
    // glyph driver (ReaderMainWindow::updateCsaGlyph) so picture and numbers
    // agree. setCsaTensor shows the section (iff this stays the focused atom);
    // clearCsaTensor hides it.
    void setCsaTensor(std::size_t atom, const CsaTensorInfo& info);
    void clearCsaTensor();

public slots:
    // The dock's two inputs: which atom and which frame. Both cause a
    // full tree rebuild because per-frame values affect many fields.
    void setPickedAtom(std::size_t atomIdx);
    void setFrame(int t);

    // Clear the tree (e.g. load unmounted or picker cleared).
    void clearSelection();

private slots:
    // The conformation finished loading `frame`'s snapshot; if it is the
    // parked frame, rebuild to show the full per-frame detail. Async-shaped:
    // v1 loads synchronously so this fires inside requestSnapshot, but the
    // committed prefetch increment will fire it from a worker handoff.
    void onSnapshotReady(std::size_t frame);

private:
    void rebuild();
    void populateIdentity(QTreeWidgetItem* parent);
    void populatePerFrame(QTreeWidgetItem* root);
    void populateCsa(QTreeWidgetItem* root);

    QPointer<QTreeWidget>         tree_;
    const model::QtProtein*       protein_      = nullptr;
    QPointer<model::Conformation> conformation_;
    std::shared_ptr<const model::TrajectoryFieldAvailability> availability_;
    bool                         hasSelection_ = false;
    std::size_t                  atomIdx_      = 0;
    int                          frame_        = 0;
    // CSA tensor mirror (fed by ReaderMainWindow); shown iff csaAtom_ == atomIdx_.
    bool                         hasCsa_       = false;
    std::size_t                  csaAtom_      = 0;
    CsaTensorInfo                csa_;
};

}  // namespace h5reader::app
