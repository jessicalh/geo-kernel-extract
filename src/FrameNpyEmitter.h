#pragma once
//
// FrameNpyEmitter -- opt-in per-frame NPY writer for trajectory runs.
//
// Symmetric to FramePdbEmitter. Singleton, fixed-shape, all state
// file-local in the .cpp; the public surface is three static methods.
//
//   Configure(protein, config) -- stash a const-Protein handle + the
//                                  Config (output dir). Called once at
//                                  startup from nmr_extract in
//                                  trajectory mode.
//                                  Does NOT touch the protein: at this
//                                  point in trajectory mode the protein
//                                  is not yet finalized (no rings,
//                                  bonds, LegacyAmberTopology).
//   OnFrame(conf, frame_idx, time_ps)
//                                  -- called per-frame from
//                                  Trajectory::Run. If not configured,
//                                  early-returns. On first call (the
//                                  protein is now finalized by
//                                  tp.Seed) emits the per-protein
//                                  CategoryInfoProjection and
//                                  TopologySidecar into output_dir.
//                                  Per dispatched frame: creates
//                                  output_dir/frame_NNNNNN/ and calls
//                                  ConformationResult::WriteAllFeatures,
//                                  mirroring the single-conformation
//                                  output of modes 1-4. No local gate;
//                                  the single --stride already selected
//                                  the frame.
//   Reset()                       -- clear configuration. For tests.
//
// Cost reality. WriteAllFeatures emits many NPY files per accepted
// frame, and wide per-atom arrays such as AIMNet2 embeddings dominate
// disk use. The H5 already carries the trajectory-scope time series;
// per-frame NPYs serve a different consumer (calibration scripts that
// read NPY without h5py, cross-tool numpy ingestion).
//
// Deliberately not a TrajectoryResult / ConformationResult. Holds no
// Welford / DenseBuffer / Selection state; participates in no
// dependency graph; emits no H5. Projection-only output, the
// per-frame NPY analog of FramePdbEmitter.
//

#include <cstddef>
#include <filesystem>

namespace nmr {

class Protein;
class ProteinConformation;

class FrameNpyEmitter {
public:
    struct Config {
        std::filesystem::path output_dir;  // empty = inert
    };

    static void Configure(const Protein& protein, Config config);
    static void OnFrame(const ProteinConformation& conf, std::size_t frame_idx, double time_ps);
    static void Reset();

    static bool IsActive();

    FrameNpyEmitter() = delete;
    FrameNpyEmitter(const FrameNpyEmitter&) = delete;
    FrameNpyEmitter& operator=(const FrameNpyEmitter&) = delete;
};

}  // namespace nmr
