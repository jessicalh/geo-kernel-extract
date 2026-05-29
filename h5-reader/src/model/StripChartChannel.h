// StripChartChannel — the strip chart's value-buffer seam.
//
// The user's design (2026-05-27): "a graph would load storage objects for its
// displayed data and we'd feed those ... a vector of the thing we are loading,
// not a fancy rewrite of it — the graph would know how to render that loaded
// object we already have, and the vector would be our buffer of data. You still
// cannot see the future and need a real strip chart, but the past would remain."
//
// Two plain data/logic types, NO QObject and NO `Qt` prefix (born honest — they
// are DATA, the way QtConformationSnapshot is data: no signals, no thread
// affinity, no lifetime in the object graph, so nothing to gain from QObject):
//
//   ChannelBuffer — the authoritative history of ONE plotted series. A
//     std::vector<double> we own (≈ 8 MB at 1e6 frames — held whole, NO
//     decimation), grown one append per playback frame. The PAST stays; the
//     FUTURE is never written (a real strip chart, not a cursor sliding over a
//     pre-built curve). An invalid frame appends NaN with valid=0 (NOT a
//     synthetic zero — fixes the old dock's leading-zero bug) so the renderer
//     can break the line at a gap (e.g. an uncomputed DFT frame).
//
//   ChannelSource — how a buffer fills: a small functor returning ONE value for
//     ONE frame (nullopt = no value at this frame → a gap). The geometry source
//     wraps model::Measure over the selection; the DFT source queries the
//     DftShieldingStore for the focus atom. DashboardDisplayController builds
//     the functor because it coordinates the Conformation / selection / DFT store.
//
// The presentation (Qt Charts now, a richer / VTK widget later) is FED these
// buffers and renders the visible window — it does not own or transform them.
// That is the swap-the-view-keep-the-data seam; these two types are the durable
// half, the renderer the throwaway half.

#pragma once

#include <QString>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <vector>

namespace h5reader::model {

// One plotted series' authoritative history. Frame-indexed and contiguous from
// frame 0: values[f] is the value at frame f, so values.size() == #frames
// appended and lastFrame() == size()-1. The owner appends contiguously (it
// backfills lastFrame()+1 .. t on a forward jump) so this invariant holds; this
// type just stores and tracks the data range.
struct ChannelBuffer {
    QString id;     // stable key, e.g. "geometry" / "dft.total.T0"
    QString label;  // axis / readout label, e.g. "dihedral" / "DFT σ (total)"
    QString unit;   // "Å" / "°" / "ppm"

    std::vector<double>  values;  // one per appended frame; NaN where !valid
    std::vector<uint8_t> valid;   // 1 == values[i] is a real measurement

    // Data range over the VALID samples only (drives the y-axis). hasRange is
    // false until the first valid sample, so a leading run of gaps does not pin
    // the range to a fake 0 (the old dock's bug).
    double yMin     = 0.0;
    double yMax     = 0.0;
    bool   hasRange = false;

    bool        empty() const { return values.empty(); }
    std::size_t size() const { return values.size(); }
    // Highest frame index present; -1 when empty. long long so the empty case
    // is a clean sentinel rather than an unsigned wrap.
    long long lastFrame() const { return static_cast<long long>(values.size()) - 1; }

    void clear() {
        values.clear();
        valid.clear();
        yMin = yMax = 0.0;
        hasRange = false;
    }

    // Append the next frame's value. nullopt or non-finite → a gap (NaN, valid
    // 0) that leaves the data range untouched. The frame index is implicit
    // (== the prior size), keeping the contiguous-from-0 invariant.
    void append(std::optional<double> v) {
        if (v && std::isfinite(*v)) {
            const double x = *v;
            values.push_back(x);
            valid.push_back(1);
            if (!hasRange) {
                yMin = yMax = x;
                hasRange = true;
            } else {
                yMin = std::min(yMin, x);
                yMax = std::max(yMax, x);
            }
        } else {
            values.push_back(std::numeric_limits<double>::quiet_NaN());
            valid.push_back(0);
        }
    }
};

// How a ChannelBuffer fills, one frame at a time. nullopt == no value at this
// frame (a gap). Kept deliberately thin: the owner supplies the closure, which
// captures whatever it needs (Conformation + selection for geometry; the DFT
// store + focus atom for shielding).
struct ChannelSource {
    QString                                              id;
    QString                                              label;
    QString                                              unit;
    std::function<std::optional<double>(std::size_t frame)> sample;
};

}  // namespace h5reader::model
