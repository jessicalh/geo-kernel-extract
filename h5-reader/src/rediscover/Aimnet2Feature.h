// Aimnet2Feature -- focused per_atom_feature cell for AIMNet2 charge,
// charge-response-gradient (CRG), and 256-d embedding.

#pragma once

#include "Aimnet2FeatureSink.h"
#include "AnalysisBody.h"

namespace h5reader::rediscover {

Aimnet2FeatureStats RunAimnet2PerAtomFeature(const Body& body, Aimnet2FeatureSink& sink);

}  // namespace h5reader::rediscover
