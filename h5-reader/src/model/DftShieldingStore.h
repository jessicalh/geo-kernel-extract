// Lazy per-frame ORCA shielding provider. Original trajectory frame indices
// map to the exact meta.json paths declared by the LGS manifest; no directory
// discovery is performed. One validated frame is resident at a time, while
// dashboard channels retain any longer display history they need.
//
// A frame is exposed only when its atom count matches the topology, every atom
// was parsed, and total = diamagnetic + paramagnetic isotropic shielding. A
// missing or invalid frame remains an honest gap and is logged at the loader.

#pragma once

#include "DftShielding.h"

#include "../io/CalcsetManifest.h"

#include <QObject>
#include <QString>

#include <cstddef>
#include <memory>
#include <optional>
#include <thread>
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
    // The protein outlives this store and is used to validate parsed frames.
    DftShieldingStore(const QtProtein* protein,
                      const std::vector<h5reader::io::DftFrame>& frames,
                      QObject* parent = nullptr);
    ~DftShieldingStore() override;

    std::size_t jobCount() const { return metaByOriginal_.size(); }

    // Does a DFT job exist on disk for this original frame index? (cheap map
    // lookup — distinguishes a "not computed" gap from "not yet parsed".)
    bool hasJob(std::size_t originalIndex) const;

    // Current-resident-or-null; NEVER parses or blocks. null == not resident:
    // call requestFrame() and react to frameReady().
    const DftShieldingFrame* frame(std::size_t originalIndex) const;

    // Parse and validate `originalIndex`, make that frame resident, then emit
    // frameReady(). This method blocks on the calling thread.
    // Idempotent for the resident frame or a known-absent frame. A job that does
    // not exist, or fails validation, is remembered as absent so it is not
    // re-attempted every frame.
    void requestFrame(std::size_t originalIndex);

    // Same residency contract as requestFrame(), but returns immediately and
    // emits frameReady() from the GUI thread after the parser finishes. This is
    // for live display refreshes that must not make frame advance wait on ORCA
    // output parsing. Dashboard history sampling and REST probes still use the
    // blocking requestFrame() path when they need a deterministic sample now.
    void requestFrameAsync(std::size_t originalIndex);

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
    void finishAsyncFrame(std::size_t originalIndex,
                          std::shared_ptr<const DftShieldingFrame> frame);

    const QtProtein* protein_ = nullptr;

    // originalIndex -> absolute meta.json path (built once at construction).
    std::unordered_map<std::size_t, QString> metaByOriginal_;

    // The single parsed frame currently exposed to observers. Persistent chart
    // history lives in strip ChannelBuffers, not here.
    std::optional<std::size_t> residentOriginal_;
    std::shared_ptr<const DftShieldingFrame> residentFrame_;

    // Negative cache only: no parsed data is retained for these frames.
    std::unordered_set<std::size_t> resolvedAbsent_;

    // One worker serves live-glyph requests. A later request replaces the
    // pending frame, and the completed worker is joined before reuse.
    std::unordered_set<std::size_t> asyncInFlight_;
    std::optional<std::size_t> asyncPendingOriginal_;
    std::thread asyncThread_;
};

}  // namespace h5reader::model
