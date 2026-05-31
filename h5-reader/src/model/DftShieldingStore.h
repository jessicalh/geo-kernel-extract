// DftShieldingStore — per-frame DFT shielding provider for the strip chart.
//
// Maps an ORIGINAL trajectory frame index (the key shared by the H5
// frame_indices, the per-frame npys/frame_NNNNNN dirs, and the ORCA job
// dirs) to that frame's ORCA shielding, parsed lazily from the
// successful .out (meta.json -> files.out_primary, NOT a glob — a frame
// may have retry .out files). The DFT campaign is partial (≈500 of 751
// frames computed at present), so a missing job is an honest GAP, never
// a faked value.
//
// Post-2026-05-31 SIMPLIFY: the `originalIndex -> meta.json` map is
// taken directly from the `.LGS`'s `dft.frames[]` array; this store no
// longer parses `_fNNNNNN_t<ps>` from job-dir names. The store accepts
// a typed `DftFrame` list at construction; the rest of the lazy/single-
// resident behaviour is unchanged.
//
// A full Qt citizen (QObject), deliberately — it owns lazy, eventually-
// async file I/O and emits a readiness signal. Architecturally this is
// just one per-frame source provider, orthogonal to FrameNpyLoader and
// any other reader: load one frame's source data, let observers sample
// it, then release it. It mirrors the frame-source contract:
//   * sample()/frame() are CHEAP — current-frame-or-null, they never parse or block;
//   * requestFrame() parses/validates one frame, makes it resident for observers,
//     then emits frameReady().
// The dashboard strips own persistent display history. This store does
// not accumulate parsed DFT frames; once another frame is requested,
// the prior parsed frame is released and the source data is effectively
// back on disk.
//
// Validation before a frame is exposed (the loader is strict even
// though the parser is permissive — "be smart, not fluff it", user
// 2026-05-27):
//   * atom count == topology atom count;
//   * no parser holes (every atom has a real element, not the default Unknown);
//   * the ORCA identity total == dia + para holds (T0 suffices:
//     decomposition is linear). A frame that fails is logged at the
//     seam and treated as absent.

#pragma once

#include "DftShielding.h"

#include "../io/CalcsetManifest.h"

#include <QObject>
#include <QString>

#include <cstddef>
#include <memory>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace h5reader::model {

class QtProtein;

// Which shielding part a chart channel reads. total = dia + para.
enum class DftPart { Total, Dia, Para };

// Which scalar of the (rank-2) shielding tensor a channel plots. T0 is the
// isotropic shielding (the headline NMR number, ppm); |T2| is the anisotropy.
// The full tensor is kept end-to-end; a channel selects ONE scalar to draw.
enum class DftScalar { IsotropicT0, AnisotropyT2 };

class DftShieldingStore final : public QObject {
    Q_OBJECT
public:
    // frames: the `.LGS`'s typed `dft.frames[]` list (frame_index +
    //         resolved meta.json abspath per entry). protein: topology
    //         spine (atom-count validation; outlives the store — both
    //         window-lived). The ctor copies the frames vector into a
    //         hash for O(1) lookup; .out files are parsed lazily on
    //         requestFrame().
    DftShieldingStore(const QtProtein* protein,
                      const std::vector<h5reader::io::DftFrame>& frames,
                      QObject* parent = nullptr);
    ~DftShieldingStore() override = default;

    std::size_t jobCount() const { return metaByOriginal_.size(); }

    // Does a DFT job exist on disk for this original frame index? (cheap map
    // lookup — distinguishes a "not computed" gap from "not yet parsed".)
    bool hasJob(std::size_t originalIndex) const;

    // Current-resident-or-null; NEVER parses or blocks. null == not resident:
    // call requestFrame() and react to frameReady().
    const DftShieldingFrame* frame(std::size_t originalIndex) const;

    // Parse + validate `originalIndex`, make that single parsed frame resident,
    // then emit frameReady(). v1: SYNCHRONOUS (parses on the calling thread).
    // Idempotent for the resident frame or a known-absent frame. A job that does
    // not exist, or fails validation, is remembered as absent so it is not
    // re-attempted every frame.
    void requestFrame(std::size_t originalIndex);

    // Cheap chart sample: resident value for (atom, part, scalar), or nullopt when
    // the frame is not resident / absent / the atom is out of range. Never
    // parses — the caller drives loading with requestFrame().
    std::optional<double> sample(std::size_t originalIndex, std::size_t atom,
                                 DftPart part, DftScalar scalar) const;

signals:
    // Emitted when frame(originalIndex) has become resolved (a valid frame is
    // now resident, OR it was determined absent/invalid — check frame()/hasJob).
    void frameReady(std::size_t originalIndex);

private:
    // Read meta.json -> files.out_primary, parse that .out, validate against the
    // topology. Returns null (and logs at the seam) on any failure.
    std::shared_ptr<const DftShieldingFrame> loadAndValidate(std::size_t originalIndex) const;

    const QtProtein* protein_ = nullptr;

    // originalIndex -> absolute meta.json path (built once at construction).
    std::unordered_map<std::size_t, QString> metaByOriginal_;

    // The single parsed frame currently exposed to observers. Persistent chart
    // history lives in strip ChannelBuffers, not here.
    std::optional<std::size_t> residentOriginal_;
    std::shared_ptr<const DftShieldingFrame> residentFrame_;

    // Negative cache only: no parsed data is retained for these frames.
    std::unordered_set<std::size_t> resolvedAbsent_;
};

}  // namespace h5reader::model
