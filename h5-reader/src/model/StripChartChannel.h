// Value storage and sampling for strip charts. These are plain data types with
// no QObject lifetime or thread affinity.
//
//   ChannelBuffer — the retained history of one plotted series. It grows by one
//     sample per playback frame without decimation. An invalid frame appends
//     NaN with valid=0 so the renderer can break the line at a gap.
//
//   ChannelSource — a functor returning one value for one frame (nullopt means a
//     gap). The geometry source wraps model::Measure over the selection; the DFT
//     source queries DftShieldingStore for the focus atom.
//     DashboardDisplayController builds the functor because it coordinates the
//     conformation, selection, and DFT store.
//
// Renderers consume these buffers without owning or transforming them.

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

    void replace(std::size_t index, std::optional<double> value) {
        if (index >= values.size() || index >= valid.size())
            return;
        if (value && std::isfinite(*value)) {
            values[index] = *value;
            valid[index] = 1;
        } else {
            values[index] = std::numeric_limits<double>::quiet_NaN();
            valid[index] = 0;
        }
        hasRange = false;
        yMin = yMax = 0.0;
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i >= valid.size() || valid[i] == 0 || !std::isfinite(values[i]))
                continue;
            if (!hasRange) {
                yMin = yMax = values[i];
                hasRange = true;
            } else {
                yMin = std::min(yMin, values[i]);
                yMax = std::max(yMax, values[i]);
            }
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
