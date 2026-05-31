// RingCurrentNeighborhood — the Pople ring-current extraction.
//
// Stratum: aromatic ring-facing hydrogens (QtAtom::IsAromaticRingHydrogen()).
// Sources: nearby aromatic rings, their per-frame geometry READ from the H5
// ring-neighbourhood time-series (distance, ρ, z, in-plane angle) so little is
// recomputed. Ring identity comes from ring_membership_per_atom — the FRAME-0
// snapshot the H5 freezes (NOT "every ring that ever enters cutoff"). Ring
// physics (intensity, JB offset, nitrogen count, aromaticity, size, fused)
// comes from the typed QtRing virtuals. cosθ = z / distance about the ring
// normal; dipolar term (3cos²θ − 1)/r³.
//
// The atom's own local frame is the aromatic-H ring-normal frame built from
// the ring the H belongs to (its parent heavy atom's ring membership).

#pragma once

#include "RediscoveryExtraction.h"

namespace h5reader::rediscover {

class RingCurrentNeighborhood final : public RediscoveryExtraction {
public:
    QString name() const override { return QStringLiteral("ring_current"); }
    FeatureSchema schema() const override;
    std::size_t extract(const Body& body, RecordSink& sink) const override;
};

}  // namespace h5reader::rediscover
