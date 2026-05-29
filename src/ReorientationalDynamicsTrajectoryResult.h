#pragma once
//
// ReorientationalDynamicsTrajectoryResult: per bond-vector Lipari-Szabo
// model-free order parameters from MD, across the whole protein backbone.
// v1 covers the backbone vectors N-H, Calpha-Halpha, and C=O (per
// residue). On top of the assumption-light model-free core it adds the 15N
// relaxation observables -- the Lipari-Szabo spectral density J(omega) and
// the 15N R1/R2/{1H}-NOE -- for the N-H vectors only, at a configurable
// reporting field. These ride on the SAME global tau_m and per-vector
// S^2/tau_e, so they inherit the tau_m_converged flag: on a run too short
// for tau_m (1P9J at 15 ns) the rates are still computed and emitted, but
// flagged unreliable, never silently trusted. Gyromagnetic ratios, hbar,
// r_NH, and the 15N CSA are cited in PhysicalConstants.h; the reporting
// field (which genuinely varies per experiment) is the CalculatorConfig
// key relaxation_field_tesla.
//
// Global tumbling removed by per-frame Kabsch superposition of the
// backbone heavy atoms (N/CA/C/O) onto frame 0 (the KabschRotation block
// is cloned from RmsdTrackingTrajectoryResult, PATTERNS.md 17). The
// body-frame second-rank TCF gives the internal order parameter and the
// effective internal correlation time; the lab-frame TCF carries overall
// tumbling, from which a single global tau_m is estimated with a
// convergence diagnostic.
//
// S^2 -- Henry & Szabo 1985 (J. Chem. Phys. 82, 4753) order-tensor
// estimator on the body-frame unit-vector components averaged over the run:
//   S^2 = (3/2)[ <x^2>^2 + <y^2>^2 + <z^2>^2
//              + 2<xy>^2 + 2<xz>^2 + 2<yz>^2 ] - 1/2
// The 3/2 multiplies the WHOLE six-term bracket; off-diagonals weight 2.
// Computed in the body frame (the lab-frame average gives ~0).
//
// tau_e -- Lipari & Szabo 1982 area method on the converged part:
//   tau_e = integral_0^inf [C_I(t) - S^2] / (1 - S^2) dt, truncated at the
//   first lag where C_I <= S^2 (avoids integrating tail noise).
//
// tau_m -- the lab-frame TCF of a rigid vector (S^2 ~ 1) is ~ C_O(t), the
// overall-tumbling correlation function; tau_m is the area of C_O averaged
// over the high-S^2 N-H set. Honest on long runs (Trp-cage at 1 us is
// hundreds of tau_m), unreliable on short ones (1P9J at 15 ns is a few):
// emitted with trajectory_length_over_tau_m and a converged flag, never
// silently trusted.
//
// Deferred (documented): sidechain X-H, the methyl symmetry axis (with the
// 1/9 fast-rotation factor), and aromatic C-H; and the dipole-CSA
// cross-correlated rates (eta_xy/eta_z), which need the relative dipolar/CSA
// tensor orientation this v1 does not carry.
//
// Lifecycle: FO. Per vector, two Legendre TCF accumulators (body, lab) and
// six order-tensor running sums; Finalize derives S^2, tau_e, the
// orientation tensor, and the global tau_m. Result-owned arrays written
// directly (no DenseBuffer; codex).
//
// Emission /trajectory/reorientational_dynamics/ (per vector, V rows):
//   bond_vector_autocorrelation      (V, L) float64  internal C_I(k), C_I(0)=1
//   bond_vector_autocorrelation_lab  (V, L) float64  lab-frame TCF
//   order_parameter_S2               (V,)   float64  Henry-Szabo, body frame
//   lipari_szabo_tau_e               (V,)   float64  ps, area method
//   bond_orientation_tensor          (V, 3, 3) float64  body-frame <u (x) u>
//   vector_kind                      (V,)   uint8   1=NH, 2=CaHa, 3=CO
//   owning_atom, tail_atom, head_atom (V,)  int32
//   residue_index                    (V,)   int32
//   lag_frames (L,) uint64 ; lag_times_ps (L,) float64
//   -- 15N relaxation (N-H rows only; NaN for CaHa/CO) --
//   spectral_density_j               (V, 5) float64  J at the KTB combination
//       frequencies [0, wN, wH-wN, wH, wH+wN], seconds (SIGNED; 15N gamma_N<0
//       so |wH-wN| is the high combination and |wH+wN| the low)
//   relaxation_R1, relaxation_R2     (V,)   float64  s^-1
//   relaxation_NOE                   (V,)   float64  dimensionless
//   relaxation_larmor_freqs_rad_per_s (5,)  float64  the J sampling frequencies
//   attrs: tau_m_ps, tau_m_converged, trajectory_length_over_tau_m,
//          tau_m_provenance, superposition, reference, estimator,
//          relaxation_field_tesla, relaxation_proton_larmor_MHz,
//          relaxation_nh_bond_length_A, relaxation_n15_csa_ppm,
//          relaxation_spectral_density_model, relaxation_equations,
//          relaxation_reliability, ...
//

#include "TrajectoryResult.h"
#include "TrajectorySpectral.h"
#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ReorientationalDynamicsTrajectoryResult : public TrajectoryResult {
public:
    // Lag count from the CalculatorConfig parameter `dynamics_n_lags`
    // (default 120), read at Create into n_lags_.

    // Vector kinds, also the emitted vector_kind codes.
    enum class Kind : std::uint8_t { None = 0, NH = 1, CaHa = 2, CO = 3 };

    std::string Name() const override {
        return "ReorientationalDynamicsTrajectoryResult";
    }

    // No declared dependency: positions and the Residue backbone cache are
    // present after tp.Seed (PATTERNS.md 15).
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<ReorientationalDynamicsTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    std::size_t NumVectors() const { return n_vectors_; }

private:
    // Per-vector identity, captured at Create.
    std::vector<std::uint8_t>  kind_;          // Kind
    std::vector<std::size_t>   tail_atom_;     // vector = head - tail
    std::vector<std::size_t>   head_atom_;
    std::vector<std::size_t>   owning_atom_;   // H (NH), HA (CaHa), O (CO)
    std::vector<std::int32_t>  residue_index_;

    // Backbone-heavy alignment set (N/CA/C/O) + frame-0 reference.
    std::vector<std::size_t> align_atoms_;
    std::vector<Vec3>        reference_positions_;
    bool                     reference_captured_ = false;
    std::size_t              n_lags_ = 120;  // CalculatorConfig dynamics_n_lags

    // Per-vector accumulators: body-frame and lab-frame Legendre TCF, and
    // six body-frame order-tensor sums (xx, yy, zz, xy, xz, yz), flat.
    std::vector<LegendreTcfAccumulator> body_acc_;
    std::vector<LegendreTcfAccumulator> lab_acc_;
    std::vector<double> order_sum_;            // V * 6

    std::vector<double> frame_times_;

    // Finalized, result-owned.
    std::vector<double> acf_internal_;   // V * L
    std::vector<double> acf_lab_;        // V * L
    std::vector<double> s2_;             // V
    std::vector<double> tau_e_;          // V
    std::vector<double> orient_tensor_;  // V * 9
    double tau_m_ps_ = 0.0;
    double traj_len_over_tau_m_ = 0.0;
    bool   tau_m_converged_ = false;

    // Relaxation layer (15N-1H, NH vectors only; NaN elsewhere). Derived at
    // Finalize from the global tau_m + per-vector S^2/tau_e at the field
    // CalculatorConfig("relaxation_field_tesla"); reliability rides the
    // tau_m_converged flag above.
    std::vector<double> spectral_density_j_;  // V * 5 : J at [0, wN, wH-wN, wH, wH+wN], seconds
    std::vector<double> r1_;                  // V : 15N R1 (s^-1)
    std::vector<double> r2_;                  // V : 15N R2 (s^-1)
    std::vector<double> noe_;                 // V : 15N{1H} steady-state NOE
    std::vector<double> relax_freqs_;         // 5 : [0, wN, wH-wN, wH, wH+wN] rad/s
    double relaxation_field_tesla_ = 0.0;
    double proton_larmor_mhz_ = 0.0;

    std::size_t n_vectors_ = 0;
    std::size_t n_frames_  = 0;
    double sample_interval_ps_ = 0.0;
    bool finalized_ = false;
};

}  // namespace nmr
