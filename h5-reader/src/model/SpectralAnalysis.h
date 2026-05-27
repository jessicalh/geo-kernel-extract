// SpectralAnalysis — power spectrum of a uniformly-sampled scalar series.
//
// The killer app charts a geometric observable (distance / angle / dihedral) of
// a selected atom tuple over a trajectory. Its FORWARD FOURIER TRANSFORM turns
// that time record into a frequency spectrum: a periodic motion (an aromatic
// ring flip, a backbone libration, a sidechain rotamer hop) shows up as a peak
// at its characteristic frequency. This is how the reader SEES periodicity in a
// dihedral — directly observable and assumption-free, which is the point.
//
// (An orientational-ACF / Lipari-Szabo S²/τc path was prototyped and pulled:
// in the lab frame, without removing global tumbling, those relaxation numbers
// may not validly apply to our data — dihedral periodicity is the honest, more
// compelling demonstration. See notes/PLANNED_ANALYSIS_METHODS.md.)
//
// Pure math over a plain series — no Qt, no VTK, no rendering. Eigen's FFT
// (unsupported/Eigen/FFT, kissfft backend) does the transform; the model layer
// already depends on Eigen, so this adds no new dependency.

#pragma once

#include <vector>

namespace h5reader::model {

// One-sided power spectrum of a real series. Index 0 is the DC bin (≈ 0 after
// mean subtraction); the rest are the physical frequencies. `valid` is false
// when the record is too short (< 4 samples) or the sample spacing is not
// positive.
struct PowerSpectrum {
    std::vector<double> frequencyPerNs;  // one-sided bins, ns^-1 (length N/2 + 1)
    std::vector<double> power;           // |X_k|^2 of the windowed series
    double dominantPeriodPs = 0.0;       // period (ps) of the strongest non-DC bin
    bool   valid = false;
};

// Compute the power spectrum of `series` sampled every `dtPicoseconds`.
//
// The series is MEAN-SUBTRACTED (the DC term carries no periodic information and
// would otherwise dominate bin 0) and HANN-WINDOWED (the trajectory is a finite,
// generally non-periodic record; tapering the ends stops a real period from
// smearing across neighbouring bins — spectral leakage) before the transform.
// Frequencies are the one-sided bins f_k = k / (N * dt), reported in ns^-1 for
// readable MD numbers; the dominant period is 1 / f_peak over the non-DC bins.
PowerSpectrum ComputePowerSpectrum(const std::vector<double>& series, double dtPicoseconds);

}  // namespace h5reader::model
