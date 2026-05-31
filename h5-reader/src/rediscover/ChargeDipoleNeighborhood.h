// ChargeDipoleNeighborhood — FF14SB charge-site dipole extraction.
//
// Stratum: backbone amide HN (same predicate as McConnell). Sources: charged
// protein atoms in the charge-sites KD cloud within cutoff_A, excluding the
// target atom's own residue for a through-space sum. Per source: FF14SB charge
// q_i and displacement (r_i - r_atom) in the target HN amide-plane frame.
// Reducer: mu = Σ q_i (r_i - r_atom), emitted as a local-frame Vec3.

#pragma once

#include "RediscoveryExtraction.h"

namespace h5reader::rediscover {

class ChargeDipoleNeighborhood final : public RediscoveryExtraction {
public:
    QString name() const override { return QStringLiteral("charge_dipole"); }
    FeatureSchema schema() const override;
    std::size_t extract(const Body& body, RecordSink& sink) const override;

    // Required and recorded source-discovery cutoff (Å), per conventions.
    double cutoff_A = 6.0;
    QString charge_source = QStringLiteral("ff14sb");
    bool exclude_residue = true;
};

}  // namespace h5reader::rediscover
