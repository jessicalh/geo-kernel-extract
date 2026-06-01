// BuckinghamEfield — focused per_atom_feature cell for APBS E-field -> DFT T0.
//
// The spine builds the local backbone frame and projects apbs_efield into it.
// Python fits only the emitted scalar invariants: E_proj and |E|^2.

#pragma once

#include "AnalysisBody.h"
#include "BuckinghamEfieldSink.h"

namespace h5reader::rediscover {

BuckinghamEfieldStats RunBuckinghamEfieldPerAtomFeature(const Body& body,
                                                        BuckinghamEfieldSink& sink);

}  // namespace h5reader::rediscover
