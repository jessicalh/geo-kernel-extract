// ExtractionSupport — shared helpers both extractions use to build the typed
// record: the DFT target (raw + library-basis decomposition + local-frame
// rotation), the identity columns, and the shared schema column groups. Keeps
// RingCurrentNeighborhood and McConnellNeighborhood thin.

#pragma once

#include "LocalFrameBasis.h"
#include "RecordSink.h"
#include "RediscoverTypes.h"
#include "RunData.h"

#include <vector>

namespace h5reader::rediscover {

// Fill the identity / frame / local-frame fields of `rec` from the run.
// h5Row is the H5 frame row; the frame map gives the original index.
void FillIdentity(NeighborhoodRecord& rec, const RunData& run, std::size_t atomIdx,
                  std::size_t h5Row, const QString& stratum, const LocalFrame& frame);

// Build the DFT target for (atom, original index), decomposing the raw total
// in the LIBRARY T2 order and rotating it into `frame`. present=false when the
// frame has no DFT job for this atom.
DftTarget BuildTarget(const RunData& run, std::size_t atomIdx, std::size_t originalIndex,
                      const LocalFrame& frame);

// Column groups shared by both extractions (identity + bare-kernel + target).
// Appended to each extraction's per-source / aggregated schemas so the column
// order matches RecordSink's writeIdentity / writeBareKernel / writeTarget.
std::vector<FeatureColumn> IdentityColumns();
std::vector<FeatureColumn> BareKernelColumns();
std::vector<FeatureColumn> TargetColumns();

// DFT frame-alignment diagnostic (resolves the T2 Cartesian-frame caveat).
// Across the DFT frames, the optimal rotation (Kabsch) between the ORCA-input
// geometry — the orientation the DFT tensors live in — and the H5 positions the
// extractor already holds. A rotation ~0° means the raw-tensor T2 components are
// in the H5 frame and are comparable to the H5 kernel-T2; a large rotation means
// the T2 components must be rotated first (T0 / |T2| are rotation-invariant and
// safe regardless). All reading stays inside the reader — no parallel H5 path.
struct DftFrameAlignment {
    int    n_frames       = 0;
    double mean_angle_deg = 0.0;
    double max_angle_deg  = 0.0;
    double mean_rmsd_A    = 0.0;  // RMSD after the optimal rotation
    double max_rmsd_A     = 0.0;
    int    n_atoms_used   = 0;    // matched atoms in the last frame (sanity)
};
DftFrameAlignment CheckDftFrameAlignment(const RunData& run);

}  // namespace h5reader::rediscover
