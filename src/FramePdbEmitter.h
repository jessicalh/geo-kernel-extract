#pragma once
//
// FramePdbEmitter -- opt-in per-frame PDB writer for trajectory runs.
//
// Singleton, fixed-shape (no inheritance, no virtuals). All state is
// file-local in the .cpp; the public surface is three static methods.
//
//   Configure(protein, config) -- assemble: stash a const-Protein handle
//                                  + the Config (output dir, stem).
//                                  Called once at startup from nmr_extract
//                                  in trajectory mode.
//   OnFrame(conf, frame_idx, time_ps, box_matrix?)
//                                  -- called per-frame from
//                                  Trajectory::Run for each dispatched
//                                  frame. If not configured, early-returns;
//                                  else emits a single PDB file. Cadence is
//                                  the single --stride; no local gate.
//   Reset()                       -- clear configuration; OnFrame
//                                  becomes inert again. For tests.
//
// Reads ONLY (no model mutation):
//   - protein.AtomAt(i) / ResidueAt(j) for identity, chain layout,
//     atom names, sequence numbers
//   - protein.LegacyAmber().BondList() for disulfides (CONECT records)
//   - per-frame ProteinConformation positions
//   - optional box matrix for CRYST1 (TRR path supplies it; PDB path
//     does not -- CRYST1 omitted when null or zero)
//
// Writes:
//   - one PDB file per dispatched frame, named
//     {stem}_f{NNNNNN}_t{ps:.1f}.pdb in output_dir
//   - HEADER + REMARK provenance, optional CRYST1, ATOM with hydrogens,
//     TER between biological chains, CONECT for disulfides only, END
//
// Deliberately not a TrajectoryResult / ConformationResult. Holds no
// Welford / DenseBuffer / Selection state; participates in no
// dependency graph; emits no H5. It is a projection-only output.
//

#include <Eigen/Dense>
#include <cstddef>
#include <filesystem>
#include <string>

namespace nmr {

class Protein;
class ProteinConformation;

class FramePdbEmitter {
public:
    struct Config {
        std::filesystem::path output_dir;  // empty = inert
        std::string           stem;        // from trajectory dir basename
    };

    static void Configure(const Protein& protein, Config config);
    static void OnFrame(const ProteinConformation& conf,
                        std::size_t frame_idx,
                        double time_ps,
                        const Eigen::Matrix3d* box_matrix = nullptr);
    static void Reset();

    // Test introspection.
    static bool IsActive();

    FramePdbEmitter() = delete;
    FramePdbEmitter(const FramePdbEmitter&) = delete;
    FramePdbEmitter& operator=(const FramePdbEmitter&) = delete;
};

}  // namespace nmr
