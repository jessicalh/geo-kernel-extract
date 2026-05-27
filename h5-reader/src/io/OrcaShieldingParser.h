// OrcaShieldingParser — parse the CHEMICAL SHIELDINGS section of an ORCA
// *_nmr.out into per-atom DFT shielding (total / diamagnetic / paramagnetic,
// each a SphericalTensor).
//
// Line-based state machine, NOT regex (structural-data discipline; ORCA's block
// layout is rigid). Qt-free + std::istream-based so the parse is directly
// testable; the app wraps it with QFile. The ORCA identity total == dia + para
// is the natural cross-check (done by the caller / test).

#pragma once

#include "../model/DftShielding.h"

#include <istream>

namespace h5reader::io {

// Parse a stream of ORCA *_nmr.out text. Returns one DftAtomShielding per
// nucleus, indexed by ORCA nucleus index (== emitted-PDB / topology order).
model::DftShieldingFrame ParseOrcaNmrShielding(std::istream& in);

}  // namespace h5reader::io
