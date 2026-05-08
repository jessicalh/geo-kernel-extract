#pragma once
//
// FramePdbEmitter -- opt-in per-frame PDB writer for trajectory runs.
//
// Singleton, fixed-shape (no inheritance, no virtuals). All state is
// file-local in the .cpp; the public surface is three static methods.
//
//   Configure(protein, config) -- assemble: stash a const-Protein handle
//                                  + the Config (output dir, stem,
//                                  decorator, stride, time window).
//                                  Called once at startup from JobSpec
//                                  when --emit-frame-pdbs is present.
//   OnFrame(conf, frame_idx, time_ps, box_matrix?)
//                                  -- called per-frame from
//                                  Trajectory::Run. If not configured,
//                                  early-returns; else gates on window
//                                  + stride and emits a single PDB
//                                  file when all gates pass.
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
//   - one PDB file per accepted frame, named
//     {stem}{_decorator?}_f{NNNNNN}_t{ps:.1f}.pdb in output_dir
//   - HEADER + REMARK provenance, optional CRYST1, ATOM with hydrogens,
//     TER between biological chains, CONECT for disulfides only, END
//
// Deliberately not a TrajectoryResult / ConformationResult. Holds no
// Welford / DenseBuffer / Selection state; participates in no
// dependency graph; emits no H5. It is a projection-only output, the
// PDB analog of the planned IUPAC / BMRB string projections on
// LegacyAmberTopology.
//

#include <Eigen/Dense>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <string>

namespace nmr {

class Protein;
class ProteinConformation;

class FramePdbEmitter {
public:
    struct Config {
        std::filesystem::path output_dir;        // empty = inert
        std::string           stem;              // from trajectory dir basename
        std::string           decorator;         // optional run tag; empty = none
        std::size_t           stride = 1;        // frames-as-read
        double                from_ps =
            -std::numeric_limits<double>::infinity();
        double                to_ps =
             std::numeric_limits<double>::infinity();
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
