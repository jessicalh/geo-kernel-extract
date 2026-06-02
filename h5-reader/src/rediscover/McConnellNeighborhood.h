// McConnellNeighborhood — the McConnell bond-anisotropy extraction.
//
// Stratum: backbone amide HN (QtAtom::IsBackboneAmideHydrogen()). Sources:
// nearby anisotropic bonds discovered PER FRAME via a nanoflann KD-tree over
// bond midpoints (FrameSpatialIndex). Per source: the BOND AXIS (unit B−A, NOT
// "C=O axis") and the form (3cos²θ − 1)/r³ where θ is the angle between the
// bond axis and the (target ← bond-midpoint) vector, r the midpoint distance.
// The Δχ susceptibility lives in the downstream parameter, NOT in the kernel —
// the kernel is /r³. Categories: PeptideCO, PeptideCN, SidechainCO, Aromatic.
//
// The atom's own local frame is the HN amide-plane frame (z = N→H; x in the
// amide plane via C_prev; y = z × x).

#pragma once

#include "RediscoveryExtraction.h"

namespace h5reader::rediscover {

class McConnellNeighborhood final : public RediscoveryExtraction {
public:
    QString name() const override { return QStringLiteral("mcconnell"); }
    FeatureSchema schema() const override;
    std::size_t extract(const Body& body, RecordSink& sink) const override;

    // Source-discovery cutoff (Å). Required at the CLI and recorded in the
    // emit; the default matches the producer McConnell bond-anisotropy cutoff.
    double cutoff_A = 10.0;
};

}  // namespace h5reader::rediscover
