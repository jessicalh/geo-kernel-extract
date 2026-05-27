// DftShieldingStore — per-frame DFT shielding provider for the strip chart.
//
// Maps an ORIGINAL trajectory frame index (the key shared by the H5
// frame_indices, the per-frame npys/frame_NNNNNN dirs, and the ORCA job dirs
// ..._fNNNNNN_t<ps>) to that frame's ORCA shielding, parsed lazily from the
// successful .out (meta.json -> files.out_primary, NOT a glob — a frame may have
// retry .out files). The DFT campaign is partial (≈500 of 751 frames computed
// at present), so a missing job is an honest GAP, never a faked value.
//
// A full Qt citizen (QObject), deliberately — it owns lazy, eventually-async
// file I/O and emits a readiness signal, exactly the manager role the
// Conformation snapshot facade plays (memory project_h5reader_buffering_decision
// _20260526). It mirrors that facade's contract:
//   * sample()/frame() are CHEAP — cached-or-null, they never parse or block;
//   * requestFrame() does the parse+validate+cache, then emits frameReady().
// The committed-later prefetch worker fills frames behind this same signature.
// The parsed DftShieldingFrame it caches is plain immutable data (DftShielding.h)
// behind shared_ptr<const> — no QObject, lock-free to hand from a future worker.
//
// Validation before a frame is exposed (the loader is strict even though the
// parser is permissive — "be smart, not fluff it", user 2026-05-27):
//   * atom count == topology atom count;
//   * no parser holes (every atom has a real element, not the default Unknown);
//   * the ORCA identity total == dia + para holds (T0 suffices: decomposition
//     is linear). A frame that fails is logged at the seam and treated as absent.

#pragma once

#include "DftShielding.h"

#include <QObject>
#include <QString>

#include <cstddef>
#include <memory>
#include <optional>
#include <unordered_map>

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
    // jobsDir: the dft/jobs directory. protein: topology spine (atom-count
    // validation; outlives the store — both window-lived). The ctor scans the
    // job-dir NAMES once (cheap, documented-convention parse of the _fNNNNNN_
    // token) to build the original-index -> dir map; .out files are parsed
    // lazily on requestFrame().
    DftShieldingStore(const QtProtein* protein, const QString& jobsDir, QObject* parent = nullptr);
    ~DftShieldingStore() override = default;

    std::size_t jobCount() const { return dirByOriginal_.size(); }

    // Does a DFT job exist on disk for this original frame index? (cheap map
    // lookup — distinguishes a "not computed" gap from "not yet parsed".)
    bool hasJob(std::size_t originalIndex) const;

    // Cached-or-null; NEVER parses or blocks. null == not resident: call
    // requestFrame() and react to frameReady().
    const DftShieldingFrame* frame(std::size_t originalIndex) const;

    // Parse + validate + cache `originalIndex`, then emit frameReady(). v1:
    // SYNCHRONOUS (parses on the calling thread). Idempotent for a resident or
    // known-absent frame. A job that does not exist, or fails validation, is
    // cached as absent so it is not re-attempted every frame.
    void requestFrame(std::size_t originalIndex);

    // Cheap chart sample: cached value for (atom, part, scalar), or nullopt when
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
    QString          jobsDir_;

    // originalIndex -> absolute job directory path (built once at construction).
    std::unordered_map<std::size_t, QString> dirByOriginal_;

    // originalIndex -> parsed frame; a null shared_ptr value means "resolved as
    // absent/invalid" (cached so we do not re-parse a known-bad frame).
    mutable std::unordered_map<std::size_t, std::shared_ptr<const DftShieldingFrame>> cache_;
};

}  // namespace h5reader::model
