// NewmanPanel -- a dashboard display form for a residue's backbone dihedral
// (the "Residue dihedral time series" metric shown as Newman dials instead of a
// line strip). Draws the residue's phi + psi as two d'Arsonval dials: the live
// current-frame torsion is the needle, behind it a translucent pastel haze of
// the whole-trajectory angle occupancy. Paints into the dark StripStack canvas;
// the projection is recomputed each paint from context.currentFrame, so it
// animates for free as playback advances.

#pragma once

#include "AbstractStripPanel.h"
#include "NewmanPaint.h"
#include "NewmanProjection.h"

#include "../model/SignalDictionary.h"

#include <QString>

#include <cstddef>

namespace h5reader::model {
class QtProtein;
class Conformation;
}

namespace h5reader::app {

class NewmanPanel final : public AbstractStripPanel {
public:
    NewmanPanel(QString title,
                const model::QtProtein* protein,
                const model::Conformation* conf,
                std::size_t residue,
                NewmanOccupancy occPhi,
                NewmanOccupancy occPsi,
                model::SignalBinding revealBinding);

    void paint(QPainter& p, const PanelGeometry& geometry, const PaintContext& context) const override;

    bool hasRevealBinding() const override { return hasReveal_; }
    model::SignalBinding revealBinding() const override { return reveal_; }

    // The two dials want a wide-ish, roughly 2:1 area (phi | psi).
    std::optional<double> preferredAspect() const override { return 2.0; }

    PanelDisplayData displayData() const override;

private:
    QString                    title_;
    const model::QtProtein*    protein_ = nullptr;
    const model::Conformation* conf_    = nullptr;
    std::size_t                residue_ = 0;
    NewmanOccupancy            occPhi_;
    NewmanOccupancy            occPsi_;
    model::SignalBinding       reveal_;
    bool                       hasReveal_ = false;
};

}  // namespace h5reader::app
