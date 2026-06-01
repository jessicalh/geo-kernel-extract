// EfgFeature — focused per_atom_feature cell for APBS EFG -> DFT T2.
//
// Shape: one row per DFT-present (atom, frame), APBS EFG T2 sidecar, DFT
// target T2 sidecar. No source set, no reducer, no registry. This is #29's
// second non-source_sum data point after broad_backbone; unification is due,
// but this spike must not disturb the frozen ring/mc oracle path.

#pragma once

#include "AnalysisBody.h"
#include "EfgFeatureSink.h"

namespace h5reader::rediscover {

EfgFeatureStats RunEfgPerAtomFeature(const Body& body, EfgFeatureSink& sink);

}  // namespace h5reader::rediscover
