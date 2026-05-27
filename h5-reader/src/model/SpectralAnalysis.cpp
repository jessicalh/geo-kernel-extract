#include "SpectralAnalysis.h"

#include <unsupported/Eigen/FFT>

#include <cmath>
#include <complex>
#include <cstddef>

namespace h5reader::model {

namespace {
constexpr double kPi = 3.14159265358979323846;  // local: portable vs <cmath> M_PI
}  // namespace

PowerSpectrum ComputePowerSpectrum(const std::vector<double>& series, double dtPicoseconds) {
    PowerSpectrum s;
    const std::size_t N = series.size();
    if (N < 4 || dtPicoseconds <= 0.0)
        return s;  // too short to resolve any period, or no time base

    // 1. Mean: the DC component (the average distance/angle) is not periodic and
    //    would otherwise swamp bin 0.
    double mean = 0.0;
    for (double v : series)
        mean += v;
    mean /= static_cast<double>(N);

    // 2. Hann window over the mean-subtracted record: w_n = 0.5 (1 - cos(2 pi n /
    //    (N-1))). The finite, non-periodic trajectory would otherwise leak a real
    //    period across neighbouring bins; the taper concentrates it.
    std::vector<double> windowed(N);
    const double twoPiOverNm1 = 2.0 * kPi / static_cast<double>(N - 1);
    for (std::size_t n = 0; n < N; ++n) {
        const double w = 0.5 * (1.0 - std::cos(twoPiOverNm1 * static_cast<double>(n)));
        windowed[n] = (series[n] - mean) * w;
    }

    // 3. Forward FFT (Eigen, kissfft backend). Real input -> full complex
    //    spectrum of length N; the real signal makes it conjugate-symmetric, so
    //    only the one-sided half 0..N/2 is independent.
    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> spectrum;
    fft.fwd(spectrum, windowed);

    // 4. One-sided power + frequency bins. f_k = k / (N dt); reported in ns^-1
    //    (dt is in ps, so * 1000). The strongest non-DC bin gives the dominant
    //    period (N dt) / k_peak, in ps.
    const std::size_t half     = N / 2;
    const double      recordPs = static_cast<double>(N) * dtPicoseconds;  // N * dt
    s.frequencyPerNs.resize(half + 1);
    s.power.resize(half + 1);

    // Strictly-positive peak gate (init 0.0, not -1.0): a flat / constant series
    // is all-zero power after mean subtraction, so NO non-DC bin clears 0 and
    // peakK stays 0 — dominantPeriodPs keeps its 0.0 default, which the readout
    // shows as "no dominant period". With the old -1.0 init the first bin always
    // "won" and a featureless series reported a bogus record-length period.
    double      peakPower = 0.0;
    std::size_t peakK     = 0;
    for (std::size_t k = 0; k <= half; ++k) {
        s.frequencyPerNs[k] = (static_cast<double>(k) / recordPs) * 1000.0;  // ps^-1 -> ns^-1
        s.power[k]          = std::norm(spectrum[k]);                        // |X_k|^2
        if (k >= 1 && s.power[k] > peakPower) {
            peakPower = s.power[k];
            peakK     = k;
        }
    }
    if (peakK >= 1)
        s.dominantPeriodPs = recordPs / static_cast<double>(peakK);

    s.valid = true;
    return s;
}

}  // namespace h5reader::model
