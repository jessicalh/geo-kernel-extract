#pragma once
//
// TrajectorySpectral: inline primitives for time-correlation and
// spectral-density estimation of a per-frame signal sampled along a
// trajectory. Free functions plus a small accumulator struct, not a
// helper class -- the same posture TrajectoryMoments.h takes for the
// Welford discipline, which PATTERNS.md's "utility namespace" rule
// sanctions for shared numerical primitives. Defining the estimator once
// here keeps the bias convention, the window, and the bounded-memory
// bookkeeping from drifting across the TrajectoryResults that
// autocorrelate kernels, bond vectors, and dihedral angles.
//
// --- Biased autocorrelation, full-range mean ---
// After Berne & Harp 1970 (Adv. Chem. Phys. 17, 63) and identical to the
// estimator in BsT0AutocorrelationTrajectoryResult:
//
//   mu     = (1/T) sum_t x_t
//   C(k)   = (1/T) sum_{t=0}^{T-k-1} (x_t - mu)(x_{t+k} - mu)   (autocovariance)
//   rho(k) = C(k) / C(0)                                        (in [-1, 1])
//
// The 1/T normalisation (not 1/(T-k)) tapers the high-variance long-lag
// tail and keeps |rho(k)| <= 1 by Cauchy-Schwarz. The autocovariance
// C(k) is the fundamental output (power units); callers normalise to
// rho(k) for the interpretable memory curve and feed C(k) to the power
// spectrum below.
//
// --- Bounded memory, EXACT (not an approximation) ---
// The full-range mean needs every sample, yet C(k) is recoverable in
// O(L) memory by accumulating, per lag k, the raw cross-product
// P_k = sum_t x_t x_{t+k}, the total sum S and count T, and the first-L
// and last-L samples for the tail corrections (codex CRITICAL,
// 2026-05-28; the refactor the BsT0 header sketches):
//
//   left_k  = sum_{t=0}^{T-k-1} x_t = S - (sum of the last  k samples)
//   right_k = sum_{t=k}^{T-1}  x_t = S - (sum of the first k samples)
//   C(k)    = ( P_k - mu (left_k + right_k) + (T-k) mu^2 ) / T
//
// A running-mean variant (subtracting a mean that drifts as samples
// arrive) is biased and is NOT used. BruteForceBiasedAutocovariance() is
// the full-history reference the bounded accumulator is unit-tested
// against for bit-level agreement.
//

#include <cmath>
#include <cstddef>
#include <vector>

namespace nmr {

// Second Legendre polynomial. P2(1) = 1, P2(0) = -1/2, P2(-1) = 1. The
// rank-2 reorientational correlation function uses this; dipolar and CSA
// relaxation are rank-2 interactions (Lipari & Szabo 1982).
inline double P2(double u) { return 0.5 * (3.0 * u * u - 1.0); }


// Full-history biased autocovariance C(k), k = 0..n_lags-1, with the
// full-range mean. O(T * n_lags) time, O(T) memory. This is the
// reference implementation; the bounded BiasedAcfAccumulator must
// reproduce it bit-for-bit. A constant signal (C(0) below eps) is the
// caller's to detect via C(0); this returns the raw covariances.
inline std::vector<double> BruteForceBiasedAutocovariance(
        const std::vector<double>& x, std::size_t n_lags) {
    const std::size_t T = x.size();
    std::vector<double> cov(n_lags, 0.0);
    if (T == 0) return cov;
    double sum = 0.0;
    for (double v : x) sum += v;
    const double mu = sum / static_cast<double>(T);
    for (std::size_t k = 0; k < n_lags && k < T; ++k) {
        double ck = 0.0;
        for (std::size_t t = 0; t + k < T; ++t) {
            ck += (x[t] - mu) * (x[t + k] - mu);
        }
        cov[k] = ck / static_cast<double>(T);
    }
    return cov;
}


// Exact O(L) streaming biased autocovariance with the full-range mean.
// One instance per signal channel. Push() each frame's scalar in order;
// Finalize() returns C(k) for k = 0..n_lags-1. Reproduces
// BruteForceBiasedAutocovariance bit-for-bit (the streaming unit test).
//
// State (all O(L)): per-lag cross-products P_k, the first-L samples
// (prefix_, for right_k), a ring buffer of the last-L samples (window_,
// for left_k and for forming cross-products as samples arrive), the
// running sum, and the count.
class BiasedAcfAccumulator {
public:
    explicit BiasedAcfAccumulator(std::size_t n_lags)
        : n_lags_(n_lags),
          p_(n_lags, 0.0),
          window_(n_lags, 0.0) {
        prefix_.reserve(n_lags);
    }

    // Add the next frame's sample. Cross-products pair the incoming x_j
    // with earlier samples x_{j-k} held in the ring buffer, so they are
    // accumulated BEFORE x_j enters the buffer (k = 0 is x_j * x_j).
    void Push(double x) {
        if (n_lags_ == 0) { sum_ += x; ++count_; return; }
        p_[0] += x * x;                                  // lag 0
        const std::size_t kmax = (window_count_ < n_lags_ - 1)
                                     ? window_count_ : n_lags_ - 1;
        for (std::size_t k = 1; k <= kmax; ++k) {
            // sample k steps back: most-recent is at window_head_-1.
            const std::size_t idx =
                (window_head_ + n_lags_ - k) % n_lags_;
            p_[k] += window_[idx] * x;
        }
        // x_j enters the ring (capacity n_lags_, overwriting oldest).
        window_[window_head_] = x;
        window_head_ = (window_head_ + 1) % n_lags_;
        if (window_count_ < n_lags_) ++window_count_;
        if (prefix_.size() < n_lags_) prefix_.push_back(x);
        sum_ += x;
        ++count_;
    }

    // C(k), k = 0..n_lags-1. Lags >= T contribute no pairs -> 0.
    std::vector<double> Finalize() const {
        std::vector<double> cov(n_lags_, 0.0);
        const std::size_t T = count_;
        if (T == 0 || n_lags_ == 0) return cov;
        const double mu = sum_ / static_cast<double>(T);

        // Suffix (last-k) sums from the ring buffer, oldest-to-newest
        // order reconstructed; prefix (first-k) sums from prefix_.
        for (std::size_t k = 0; k < n_lags_ && k < T; ++k) {
            double first_k = 0.0;            // sum of the first k samples
            for (std::size_t i = 0; i < k; ++i) first_k += prefix_[i];
            double last_k = 0.0;             // sum of the last  k samples
            for (std::size_t i = 0; i < k; ++i) {
                const std::size_t idx =
                    (window_head_ + n_lags_ - 1 - i) % n_lags_;
                last_k += window_[idx];
            }
            const double left_k  = sum_ - last_k;
            const double right_k = sum_ - first_k;
            cov[k] = (p_[k] - mu * (left_k + right_k)
                      + static_cast<double>(T - k) * mu * mu)
                     / static_cast<double>(T);
        }
        return cov;
    }

    std::size_t Count() const { return count_; }

private:
    std::size_t n_lags_;
    std::vector<double> p_;          // per-lag cross-products P_k
    std::vector<double> window_;     // ring buffer of the last n_lags_ samples
    std::size_t window_head_  = 0;   // next write index into window_
    std::size_t window_count_ = 0;   // filled entries (<= n_lags_)
    std::vector<double> prefix_;     // first n_lags_ samples, in order
    double sum_       = 0.0;
    std::size_t count_ = 0;
};


// Parzen lag window, w(0) = 1, w(M) = 0, with M = max lag (Parzen 1961,
// Ann. Math. Statist. 32, 329; de la Vallee Poussin window). Its spectral
// window is non-negative, so the Blackman-Tukey estimate below is
// guaranteed non-negative -- unlike a Hann lag window, which can produce
// negative bins (codex MAJOR, 2026-05-28).
//
//   w(k) = 1 - 6(k/M)^2 + 6(k/M)^3      0   <= k <= M/2
//        = 2(1 - k/M)^3                 M/2 <  k <= M
//        = 0                            k   >  M
inline double ParzenLagWindow(std::size_t k, std::size_t M) {
    if (M == 0 || k > M) return (k == 0) ? 1.0 : 0.0;
    const double r = static_cast<double>(k) / static_cast<double>(M);
    if (k <= M / 2) return 1.0 - 6.0 * r * r + 6.0 * r * r * r;
    const double s = 1.0 - r;
    return 2.0 * s * s * s;
}


// Blackman-Tukey power spectral density from an autocovariance sequence
// C(k), via the Parzen lag window (guaranteed S(f) >= 0). Direct cosine
// transform -- no FFT dependency; O(n_freq * n_lags), negligible at the
// lag counts here. Frequencies f_m = m / (2 L dt) for m = 0..n_freq-1,
// spanning 0 .. Nyquist = 1/(2 dt). One-sided convention: interior bins
// doubled, the f=0 and Nyquist bins carry single weight. dt is the frame
// interval (ps); S has units (signal-unit)^2 * ps.
//
//   S(f_m) = dt [ C(0) + 2 sum_{k=1}^{L-1} w(k) C(k) cos(2 pi f_m k dt) ]
//
// Feeding the autocovariance C(k) gives a power density; feeding the
// normalised rho(k) gives a unit-normalised spectral shape (both >= 0).
inline std::vector<double> ParzenPowerSpectrum(
        const std::vector<double>& cov, double dt, std::size_t n_freq) {
    const std::size_t L = cov.size();
    std::vector<double> spectrum(n_freq, 0.0);
    if (L == 0 || n_freq == 0 || !(dt > 0.0)) return spectrum;
    const double two_pi = 2.0 * M_PI;
    for (std::size_t m = 0; m < n_freq; ++m) {
        const double f = static_cast<double>(m)
                         / (2.0 * static_cast<double>(L) * dt);
        double acc = cov[0];                       // k = 0 term, w(0) = 1
        for (std::size_t k = 1; k < L; ++k) {
            acc += 2.0 * ParzenLagWindow(k, L) * cov[k]
                   * std::cos(two_pi * f * static_cast<double>(k) * dt);
        }
        // One-sided PSD: double the interior bins so integrating over
        // [0, Nyquist] recovers the total power; DC (m==0) and Nyquist
        // (m==L) carry single weight.
        const double side = (m == 0 || m == L) ? 1.0 : 2.0;
        spectrum[m] = side * dt * acc;
    }
    return spectrum;
}


// Frequency grid matching ParzenPowerSpectrum: f_m = m / (2 L dt), in
// units of 1/ps (so multiply by 1000 for 1/ns). Emitted alongside the
// spectrum so consumers read the axis rather than reconstruct it.
inline std::vector<double> SpectrumFrequenciesPerPs(
        std::size_t n_lags, double dt, std::size_t n_freq) {
    std::vector<double> f(n_freq, 0.0);
    if (n_lags == 0 || !(dt > 0.0)) return f;
    for (std::size_t m = 0; m < n_freq; ++m) {
        f[m] = static_cast<double>(m)
               / (2.0 * static_cast<double>(n_lags) * dt);
    }
    return f;
}


// Second-rank reorientational time-correlation function of a unit vector,
//
//   C_I(k) = < P2( u(t) . u(t+k) ) >,    P2(x) = (3 x^2 - 1) / 2,
//
// the orientational autocorrelation that governs dipolar and CSA
// relaxation (Lipari & Szabo 1982, J. Am. Chem. Soc. 104, 4546). UNBIASED
// estimator -- averaged over the T-k pairs available at each lag -- so
// C_I(0) = 1 and the long-time value is the order parameter S^2 (the
// plateau). NOT mean-subtracted: P2 already carries the orientational
// structure (distinct from the mean-subtracted scalar ACF above). Bounded
// O(L) memory: a ring buffer of the last L unit vectors plus per-lag P2
// sums and pair counts. The unit vector is passed as components so this
// header carries no Eigen dependency; the caller normalises. Push the
// body-frame (tumbling-removed) vector for the internal C_I, the lab-frame
// vector for the total correlation function.
class LegendreTcfAccumulator {
public:
    explicit LegendreTcfAccumulator(std::size_t n_lags)
        : n_lags_(n_lags), p2_(n_lags, 0.0), npairs_(n_lags, 0),
          wx_(n_lags, 0.0), wy_(n_lags, 0.0), wz_(n_lags, 0.0) {}

    void Push(double ux, double uy, double uz) {
        ++count_;
        if (n_lags_ == 0) return;
        p2_[0] += 1.0; ++npairs_[0];          // P2(u . u) = P2(1) = 1
        const std::size_t kmax =
            (window_count_ < n_lags_ - 1) ? window_count_ : n_lags_ - 1;
        for (std::size_t k = 1; k <= kmax; ++k) {
            const std::size_t idx = (window_head_ + n_lags_ - k) % n_lags_;
            const double dot = ux * wx_[idx] + uy * wy_[idx] + uz * wz_[idx];
            p2_[k] += P2(dot);
            ++npairs_[k];
        }
        wx_[window_head_] = ux; wy_[window_head_] = uy; wz_[window_head_] = uz;
        window_head_ = (window_head_ + 1) % n_lags_;
        if (window_count_ < n_lags_) ++window_count_;
    }

    // C_I(k) = mean P2 over the T-k pairs; NaN where no pairs (k >= T).
    std::vector<double> Finalize() const {
        std::vector<double> c(n_lags_, std::nan(""));
        for (std::size_t k = 0; k < n_lags_; ++k) {
            if (npairs_[k] > 0)
                c[k] = p2_[k] / static_cast<double>(npairs_[k]);
        }
        return c;
    }
    std::size_t Count() const { return count_; }

private:
    std::size_t n_lags_;
    std::vector<double> p2_;
    std::vector<std::size_t> npairs_;
    std::vector<double> wx_, wy_, wz_;     // ring buffer of last-L unit vectors
    std::size_t window_head_  = 0;
    std::size_t window_count_ = 0;
    std::size_t count_        = 0;
};


// Circular autocorrelation of an angle,
//
//   C(k) = < cos( theta(t+k) - theta(t) ) >
//        = < cos(t) cos(t+k) + sin(t) sin(t+k) >,
//
// the proper periodic correlation -- the naive < theta(t) theta(0) > is
// wrong across the +-pi branch cut. UNBIASED (T-k pairs per lag); C(0) = 1,
// the long-time value is < cos theta >^2 + < sin theta >^2 (the circular
// order). Bounded O(L) ring buffers of cos and sin. For phi/psi/chi
// torsional decorrelation.
class CircularAcfAccumulator {
public:
    explicit CircularAcfAccumulator(std::size_t n_lags)
        : n_lags_(n_lags), c_(n_lags, 0.0), npairs_(n_lags, 0),
          wc_(n_lags, 0.0), ws_(n_lags, 0.0) {}

    void Push(double angle_rad) {
        ++count_;
        if (n_lags_ == 0) return;
        const double cs = std::cos(angle_rad);
        const double sn = std::sin(angle_rad);
        c_[0] += 1.0; ++npairs_[0];           // cos(0) = 1
        const std::size_t kmax =
            (window_count_ < n_lags_ - 1) ? window_count_ : n_lags_ - 1;
        for (std::size_t k = 1; k <= kmax; ++k) {
            const std::size_t idx = (window_head_ + n_lags_ - k) % n_lags_;
            c_[k] += cs * wc_[idx] + sn * ws_[idx];
            ++npairs_[k];
        }
        wc_[window_head_] = cs; ws_[window_head_] = sn;
        window_head_ = (window_head_ + 1) % n_lags_;
        if (window_count_ < n_lags_) ++window_count_;
    }

    std::vector<double> Finalize() const {
        std::vector<double> c(n_lags_, std::nan(""));
        for (std::size_t k = 0; k < n_lags_; ++k) {
            if (npairs_[k] > 0)
                c[k] = c_[k] / static_cast<double>(npairs_[k]);
        }
        return c;
    }
    std::size_t Count() const { return count_; }

private:
    std::size_t n_lags_;
    std::vector<double> c_;
    std::vector<std::size_t> npairs_;
    std::vector<double> wc_, ws_;          // ring buffer of last-L cos, sin
    std::size_t window_head_  = 0;
    std::size_t window_count_ = 0;
    std::size_t count_        = 0;
};

}  // namespace nmr
