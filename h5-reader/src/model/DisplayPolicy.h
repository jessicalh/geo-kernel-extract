// DisplayPolicy -- per-field PRESENTATION limits as data, not mechanical-by-shape.
//
// The display catalog historically assigned modes by coarse value_shape in hand
// helpers (TrajectorySignalCatalog.cpp tensorStripModes()/categoryStripModes()/…),
// divorced from each field's physics -- so static topology TABLES got time-strip
// modes nothing renders, and a 256-d embedding got a meaningless strip
// (GET /catalog/display-audit, 2026-06-22). This is the seam where presentation
// limits become per-field DATA. Pure functions over SignalDescriptor; shared by
// the catalog, REST, and the headless test.
//
// A2-step-1 (here): displayability -- is a field a dashboard signal at all.
// A2-step-2 (next): derive the OFFERED MODES from the dictionary (FieldSpec irreps
// 0e/1e/2e -> tensor components; shape+axis -> surfaces), replacing the mechanical
// helpers. See notes/DISPLAY_COHERENCE_RELAY_2026-06-22.md.

#pragma once

#include "DashboardSignal.h"

#include <QLatin1String>

namespace h5reader::model {

// Is this descriptor a dashboard DISPLAYABLE signal at all? Two field kinds are
// intentionally NOT plottable dashboard signals (the mechanical-by-shape pass
// offered them dead/nonsense modes):
//   * the structural topology TABLES (atoms / residues / bonds / rings /
//     ring_membership) -- they ARE the molecule, shown in the 3-D scene, not a
//     per-frame series. (topology:bond_length is a real per-bond Scalar and STAYS
//     displayable -- hence the Category gate, not a blanket family check.)
//   * the 256-d AIMNet2 Embedding -- an ML feature vector, not a plottable curve.
inline bool IsDashboardDisplayable(const SignalDescriptor& descriptor) {
    if (descriptor.valueShape == SignalValueShape::Embedding)
        return false;
    if (descriptor.family == QLatin1String("topology")
        && descriptor.valueShape == SignalValueShape::Category)
        return false;
    return true;
}

}  // namespace h5reader::model
