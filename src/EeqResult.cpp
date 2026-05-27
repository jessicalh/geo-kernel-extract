#include "EeqResult.h"
#include "Protein.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <Eigen/Dense>
#include <cmath>

namespace nmr {

// D4 EEQ parameters (atomic units): values in PhysicalConstants.h
// (D4EeqParams, D4EeqParamsFor).
// Caldeweyher et al. 2019, DOI: 10.1063/1.5090222.


std::unique_ptr<EeqResult> EeqResult::Compute(ProteinConformation& conf) {

    const size_t N = conf.AtomCount();
    const Protein& protein = conf.ProteinRef();

    OperationLog::Scope scope("EeqResult::Compute",
        "atoms=" + std::to_string(N));

    // Empty conformation: nothing to solve. Guards the atom-0 stats seed
    // and the ÷N neutrality shift below against a zero-atom conformation.
    if (N == 0)
        return nullptr;

    // ── Configuration (TOML) ────────────────────────────────────────

    const double Q_total      = CalculatorConfig::Get("eeq_total_charge");
    const double cn_steepness = CalculatorConfig::Get("eeq_cn_steepness");
    const double cn_cutoff  = CalculatorConfig::Get("eeq_cn_cutoff");
    const double charge_clamp = CalculatorConfig::Get("eeq_charge_clamp");

    const double cn_cutoff_bohr = cn_cutoff * BOHR_PER_ANGSTROM;
    const double cn_cutoff_bohr_sq = cn_cutoff_bohr * cn_cutoff_bohr;

    // ── Pre-cache per-atom parameters and positions in Bohr ─────────

    std::vector<D4EeqParams> params(N);
    for (size_t i = 0; i < N; ++i) {
        params[i] = D4EeqParamsFor(protein.AtomAt(i).element);
    }

    Eigen::MatrixXd pos(N, 3);
    for (size_t i = 0; i < N; ++i) {
        Vec3 p = conf.PositionAt(i);
        pos(i, 0) = p.x() * BOHR_PER_ANGSTROM;
        pos(i, 1) = p.y() * BOHR_PER_ANGSTROM;
        pos(i, 2) = p.z() * BOHR_PER_ANGSTROM;
    }

    GeometryChoiceBuilder choices(conf);

    // GeometryChoice: parameter summary
    choices.Record(CalculatorId::EEQ, 0, "eeq_parameters",
        [Q_total, cn_steepness, cn_cutoff, charge_clamp, N]
        (GeometryChoice& gc) {
            AddNumber(gc, "total_charge", Q_total, "e");
            AddNumber(gc, "cn_steepness", cn_steepness, "");
            AddNumber(gc, "cn_cutoff", cn_cutoff, "A");
            AddNumber(gc, "charge_clamp", charge_clamp, "e");
            AddNumber(gc, "n_atoms", static_cast<double>(N), "count");
            AddNumber(gc, "method", 0.0, "D4_EEQ_Caldeweyher_2019");
        });

    // ── Step 1: coordination numbers (error function counting) ──────
    //
    // CN_i = Σ_{j≠i} ½ erfc(k · (R_ij/(r_cov_i + r_cov_j) - 1))
    //
    // Pairs beyond cn_cutoff are skipped — erfc is negligible there.

    Eigen::VectorXd cn = Eigen::VectorXd::Zero(N);
    long cn_pairs_counted = 0;
    long cn_pairs_skipped = 0;

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            double dx = pos(i, 0) - pos(j, 0);
            double dy = pos(i, 1) - pos(j, 1);
            double dz = pos(i, 2) - pos(j, 2);
            double dist_bohr_sq = dx * dx + dy * dy + dz * dz;
            if (dist_bohr_sq > cn_cutoff_bohr_sq) {
                ++cn_pairs_skipped;
                continue;
            }
            double dist_bohr = std::sqrt(dist_bohr_sq);
            double rcov_sum_bohr = params[i].rcov + params[j].rcov;
            double count =
                0.5 * std::erfc(cn_steepness * (dist_bohr / rcov_sum_bohr - 1.0));
            cn(i) += count;
            cn(j) += count;
            ++cn_pairs_counted;
        }
    }

    // ── Step 2: effective electronegativities ───────────────────────
    //
    // χ_eff = χ + κ · √(CN + ε)

    Eigen::VectorXd chi_eff(N);
    for (size_t i = 0; i < N; ++i)
        chi_eff(i) = params[i].chi
                   + params[i].kappa * std::sqrt(cn(i) + 1e-14);
                                                // 1e-14: dftd4 guard against sqrt(0)

    // ── Step 3: Coulomb matrix A (N×N, symmetric EEQ matrix) ────────
    //
    // Diagonal: A_ii = η_i + √(2/π)/r_i (hardness plus Gaussian self term)
    // Off-diagonal: A_ij = γ(R_ij) = 1/√(R² + 1/(η_i·η_j))
    //
    // Ohno-Klopman kernel (Ohno 1964, Klopman 1964):
    //   R→0:  γ → √(η_i·η_j)  (finite, no singularity)
    //   R→∞:  γ → 1/R          (bare Coulomb)
    //
    // Note: the dftd4 reference uses a Gaussian-charge form erf(αR)/R.
    // The Ohno-Klopman form is a deliberate divergence that preserves the
    // correct asymptotics; charges differ slightly from exact dftd4 — we
    // want geometry-responsive charges, not exact dftd4 reproduction.

    // Field aliases: gam = chemical hardness η_i; rad = Gaussian charge
    // radius r_i (rcov, used in Step 1, is the separate covalent radius).
    Eigen::MatrixXd A(N, N);
    for (size_t i = 0; i < N; ++i) {
        // Diagonal: η_i + √(2/π) / r_i (Caldeweyher 2019 Eq. 8). The second
        // term is the self-Coulomb energy of the Gaussian charge distribution.
        A(i, i) = params[i].gam + SQRT_2_OVER_PI / params[i].rad;
        for (size_t j = i + 1; j < N; ++j) {
            double dx = pos(i, 0) - pos(j, 0);
            double dy = pos(i, 1) - pos(j, 1);
            double dz = pos(i, 2) - pos(j, 2);
            double dist_bohr_sq = dx * dx + dy * dy + dz * dz;
            double hardness_product_inv = 1.0 / (params[i].gam * params[j].gam);
            double gamma = 1.0 / std::sqrt(dist_bohr_sq + hardness_product_inv);
            A(i, j) = gamma;
            A(j, i) = gamma;
        }
    }

    // ── Step 4: solve with net-charge constraint ───────────────────
    //
    // Block elimination of the augmented system:
    //   [A  1] [q] = [-χ_eff]
    //   [1' 0] [λ]   [Q     ]
    //
    // Solve A·u = χ_eff and A·v = 1, then:
    //   λ = -(Q + 1'u) / (1'v)
    //   q = -(u + λ·v)
    // In code: u ≡ A_inv_chi, v ≡ A_inv_ones, q ≡ charges.

    // Cholesky (LLT): this formulation is intended to be SPD; failure is
    // checked below instead of assuming every parameter set is benign.
    Eigen::LLT<Eigen::MatrixXd> chol(A);
    if (chol.info() != Eigen::Success) {
        OperationLog::Error("EeqResult",
            "Cholesky decomposition failed (N=" + std::to_string(N) +
            ") — Coulomb matrix not positive definite");
        return nullptr;
    }

    Eigen::VectorXd A_inv_chi = chol.solve(chi_eff);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd A_inv_ones = chol.solve(ones);

    double denom = ones.dot(A_inv_ones);
    if (std::abs(denom) < 1e-30) {
        OperationLog::Error("EeqResult", "charge constraint denominator is zero");
        return nullptr;
    }

    double lambda = -(Q_total + ones.dot(A_inv_chi)) / denom;
    Eigen::VectorXd charges = -(A_inv_chi + lambda * A_inv_ones);

    // Remove the Cholesky solve residual so the (pre-clamp) sum equals
    // Q_total.  The block elimination is algebraically exact; the solve
    // introduces a residual ~cond(A).  A uniform shift preserves per-atom
    // charge differences.  NOTE: the Step-5 clamp can move the stored sum
    // off Q_total when it fires (see n_clamped).
    double charge_residual = charges.sum() - Q_total;
    charges.array() -= charge_residual / static_cast<double>(N);

    // ── Step 5: store charges with clamp guard ──────────────────────

    auto result = std::make_unique<EeqResult>();
    result->conf_ = &conf;

    int n_clamped = 0;
    for (size_t i = 0; i < N; ++i) {
        auto& atom = conf.MutableAtomAt(i);
        atom.eeq_cn = cn(i);

        double qi = charges(i);
        if (std::abs(qi) > charge_clamp) {
            double original = qi;
            qi = (qi > 0) ? charge_clamp : -charge_clamp;
            ++n_clamped;

            // GeometryChoice: charge-clamp record
            choices.Record(CalculatorId::EEQ, i, "eeq charge clamp",
                [&conf, i, original, qi, &cn](GeometryChoice& gc) {
                    AddAtom(gc, &conf.AtomAt(i), i,
                            EntityRole::Target, EntityOutcome::Triggered);
                    AddNumber(gc, "original_charge", original, "e");
                    AddNumber(gc, "clamped_charge", qi, "e");
                    AddNumber(gc, "coordination_number", cn(i), "");
                });
        }
        atom.eeq_charge = qi;
    }

    // ── Charge statistics ───────────────────────────────────────────

    // Seed from the stored (post-clamp) charge so the reported min/max are
    // consistent with the values the loop compares and every consumer reads.
    double q_sum = 0.0;
    double q_min = conf.AtomAt(0).eeq_charge, q_max = q_min;
    for (size_t i = 0; i < N; ++i) {
        double qi = conf.AtomAt(i).eeq_charge;
        q_sum += qi;
        if (qi < q_min) q_min = qi;
        if (qi > q_max) q_max = qi;
    }

    // GeometryChoice: charge statistics summary
    choices.Record(CalculatorId::EEQ, 0, "eeq_charge_statistics",
        [q_sum, q_min, q_max, n_clamped, cn_pairs_counted, cn_pairs_skipped]
        (GeometryChoice& gc) {
            AddNumber(gc, "charge_sum", q_sum, "e");
            AddNumber(gc, "charge_min", q_min, "e");
            AddNumber(gc, "charge_max", q_max, "e");
            AddNumber(gc, "n_clamped", static_cast<double>(n_clamped), "count");
            AddNumber(gc, "cn_pairs_counted",
                      static_cast<double>(cn_pairs_counted), "count");
            AddNumber(gc, "cn_pairs_skipped",
                      static_cast<double>(cn_pairs_skipped), "count");
        });

    OperationLog::Info(LogCalcOther, "EeqResult",
        "N=" + std::to_string(N) +
        " Σq=" + std::to_string(q_sum) +
        " range=[" + std::to_string(q_min) + "," +
        std::to_string(q_max) + "]" +
        " clamped=" + std::to_string(n_clamped) +
        " cn_pairs=" + std::to_string(cn_pairs_counted));

    return result;
}


int EeqResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const size_t N = conf.AtomCount();

    // eeq_charges: (N,) — partial charges in elementary charges
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).eeq_charge;
        NpyWriter::WriteFloat64(output_dir + "/eeq_charges.npy", data.data(), N);
    }

    // eeq_cn: (N,) — coordination number (intermediate, for traceability)
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).eeq_cn;
        NpyWriter::WriteFloat64(output_dir + "/eeq_cn.npy", data.data(), N);
    }

    return 2;
}

}  // namespace nmr
