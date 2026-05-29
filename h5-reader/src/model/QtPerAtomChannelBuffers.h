// QtPerAtomChannelBuffers.h — per-atom × per-channel typed buffers.
//
// Used by KernelDynamicsTrajectoryResult's 5 outputs (acf, power_spectrum,
// decay_time, peak_freq, spectral_centroid) and KernelCoherence
// (Phase G — Matrix shape uses a sibling QtPerAtomMatrix struct).
//
// Two row shapes:
//   - PerAtomChannelCurve   = (N atoms, C channels, S samples)  rho/lag or PSD/freq
//   - PerAtomChannelScalar  = (N atoms, C channels)             per-channel reduction
//
// Both carry channel metadata (names + units) so the catalog descriptor
// channels[] list and the H5's channel_names string array stay aligned.
// External axis values (lag_times_ps or frequencies_per_ps) ride on the
// curve struct so consumers don't have to re-derive them.

#pragma once

#include <QString>
#include <QStringList>

#include <cstddef>
#include <vector>

namespace h5reader::model {

struct QtPerAtomChannelCurve {
    std::size_t n_atoms = 0;
    std::size_t n_channels = 0;
    std::size_t n_samples = 0;          // L for ACF curves, F for PSD curves
    std::vector<double> data;           // (N * C * n_samples,) row-major
    std::vector<double> axis_values;    // (n_samples,) — lag_ps or freq_per_ps
    QStringList channel_names;          // (n_channels,)
    QStringList channel_units;          // (n_channels,)
    QString axis_unit;                  // "ps" (lag) or "1/ps" (freq)
    QString axis_label;                 // "lag" or "frequency"
    QString units;                      // value units (dimensionless ACF, units² ps PSD)
    QString result_name;

    double at(std::size_t atom, std::size_t channel, std::size_t sample) const {
        if (atom >= n_atoms || channel >= n_channels || sample >= n_samples
            || data.empty()) {
            return 0.0;
        }
        return data[(atom * n_channels + channel) * n_samples + sample];
    }
};

struct QtPerAtomChannelScalar {
    std::size_t n_atoms = 0;
    std::size_t n_channels = 0;
    std::vector<double> data;           // (N * C,) row-major
    QStringList channel_names;          // (n_channels,)
    QStringList channel_units;          // (n_channels,)
    QString units;
    QString result_name;

    double at(std::size_t atom, std::size_t channel) const {
        if (atom >= n_atoms || channel >= n_channels || data.empty())
            return 0.0;
        return data[atom * n_channels + channel];
    }
};

// Per-atom NxN matrix (KernelCoherence's 13×13 Pearson per atom).
// Symmetric, diagonal=1; channel-constant rows/cols are NaN.
struct QtPerAtomMatrix {
    std::size_t n_atoms = 0;
    std::size_t n_channels = 0;
    std::vector<double> data;           // (N * C * C,) row-major
    QStringList channel_names;          // (n_channels,)
    QStringList channel_units;          // (n_channels,)
    QString units;
    QString result_name;

    double at(std::size_t atom, std::size_t a, std::size_t b) const {
        if (atom >= n_atoms || a >= n_channels || b >= n_channels || data.empty())
            return 0.0;
        return data[(atom * n_channels + a) * n_channels + b];
    }
};

// KernelCoherence composite — one matrix per atom + the shared channel
// metadata. Static (no time axis).
struct QtKernelCoherence {
    QtPerAtomMatrix matrix;
};

// Composite TR-level buffer for /trajectory/kernel_dynamics. One reader
// fills all 5 sub-buffers in one pass and stamps them with the shared
// channel/lag/freq metadata.
struct QtKernelDynamics {
    std::size_t n_atoms = 0;
    std::size_t n_channels = 0;
    double sample_interval_ps = 0.0;
    QStringList channel_names;
    QStringList channel_units;

    QtPerAtomChannelCurve   acf;                       // lag domain
    QtPerAtomChannelCurve   power_spectrum;            // frequency domain
    QtPerAtomChannelScalar  decay_time_ps;
    QtPerAtomChannelScalar  peak_freq_per_ps;
    QtPerAtomChannelScalar  spectral_centroid_per_ps;
};

}  // namespace h5reader::model
