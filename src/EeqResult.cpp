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

// ============================================================================
// D4 EEQ parameters
//
// Source: Caldeweyher, Ehlert, Hansen, Neugebauer, Spicher, Bannwarth &
// Grimme, J. Chem. Phys. 150, 154122 (2019).  DOI: 10.1063/1.5090222.
// Reference implementation: dftd4/dftd4 (Apache-2.0), src/dftd4/data/.
//
// All values in atomic units (Hartree, Bohr).  Fitted to reproduce
// Hirshfeld charges from DFT/def2-TZVP.  Same parameters used in
// TURBOMOLE, ORCA D4, xTB, and DFTB+.
// ============================================================================

// D4 EEQ parameters: see PhysicalConstants.h (D4EeqParams, D4D4EeqParamsFor).
// Caldeweyher et al. 2019, DOI: 10.1063/1.5090222.


std::unique_ptr<EeqResult> EeqResult::Compute(ProteinConformation& conf) {

    const size_t N = conf.AtomCount();
    const Protein& protein = conf.ProteinRef();

    OperationLog::Scope scope("EeqResult::Compute",
        "atoms=" + std::to_string(N));

    // ── Named constants from TOML (no buried literals) ──────────────

    const double Q_total    = CalculatorConfig::Get("eeq_total_charge");
    const double cn_k       = CalculatorConfig::Get("eeq_cn_steepness");
    const double cn_cutoff  = CalculatorConfig::Get("eeq_cn_cutoff");
    const double charge_clamp = CalculatorConfig::Get("eeq_charge_clamp");

    const double cn_cutoff_bohr = cn_cutoff * BOHR_PER_ANGSTROM;
    const double cn_cutoff_bohr_sq = cn_cutoff_bohr * cn_cutoff_bohr;

    // ── Pre-cache per-atom parameters and positions in Bohr ─────────

    std::vector<D4EeqParams> params(N);
    int n_fallback = 0;
    for (size_t i = 0; i < N; ++i) {
        Element el = protein.AtomAt(i).element;
        params[i] = D4EeqParamsFor(el);
        if (el == Element::Unknown || el == Element::S)
            ++n_fallback;  // S uses real params; Unknown uses fallback
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
        [Q_total, cn_k, cn_cutoff, charge_clamp, N](GeometryChoice& gc) {
            AddNumber(gc, "total_charge", Q_total, "e");
            AddNumber(gc, "cn_steepness", cn_k, "");
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
            double R2 = dx * dx + dy * dy + dz * dz;
            if (R2 > cn_cutoff_bohr_sq) {
                ++cn_pairs_skipped;
                continue;
            }
            double R = std::sqrt(R2);
            double rcov_sum = params[i].rcov + params[j].rcov;
            double count = 0.5 * std::erfc(cn_k * (R / rcov_sum - 1.0));
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

    // ── Step 3: Coulomb matrix A (N×N, symmetric positive definite) ─
    //
    // Diagonal: A_ii = η_i (chemical hardness)
    // Off-diagonal: A_ij = γ(R_ij) = 1/√(R² + 1/(η_i·η_j))
    //
    // Ohno-Klopman kernel (Ohno 1964, Klopman 1964):
    //   R→0:  γ → √(η_i·η_j)  (finite, no singularity)
    //   R→∞:  γ → 1/R          (bare Coulomb)
    //
    // Note: the dftd4 reference uses a Gaussian-charge form erf(αR)/R.
    // The Ohno-Klopman form is a deliberate simplification that avoids
    // the erf() evaluation cost while preserving the correct asymptotic
    // behaviour.  Charges will differ slightly from exact dftd4 output.
    // For our purpose (geometry-responsive charge patterns, not exact
    // reproduction), this is acceptable.

    Eigen::MatrixXd A(N, N);
    for (size_t i = 0; i < N; ++i) {
        // Diagonal: η_i + √(2/π) / r_i (Caldeweyher 2019 Eq. 8)
        // The second term is the self-Coulomb energy of the Gaussian charge
        // distribution.  Makes the matrix diagonally dominant → SPD.
        A(i, i) = params[i].gam + SQRT_2_OVER_PI / params[i].rad;
        for (size_t j = i + 1; j < N; ++j) {
            double dx = pos(i, 0) - pos(j, 0);
            double dy = pos(i, 1) - pos(j, 1);
            double dz = pos(i, 2) - pos(j, 2);
            double R2 = dx * dx + dy * dy + dz * dz;
            double gam_prod_inv = 1.0 / (params[i].gam * params[j].gam);
            double gamma = 1.0 / std::sqrt(R2 + gam_prod_inv);
            A(i, j) = gamma;
            A(j, i) = gamma;
        }
    }

    // ── Step 4: solve with charge neutrality constraint ─────────────
    //
    // Block elimination of the augmented system:
    //   [A  1] [q] = [-χ_eff]
    //   [1' 0] [λ]   [Q     ]
    //
    // Solve A·u = χ_eff and A·v = 1, then:
    //   λ = -(Q + 1'u) / (1'v)
    //   q = -(u + λ·v)

    // Cholesky (LLT): the matrix is SPD because the diagonal includes
    // both the hardness η and the self-Coulomb √(2/π)/r correction,
    // which together dominate the off-diagonal Ohno-Klopman γ values.
    Eigen::LLT<Eigen::MatrixXd> chol(A);
    if (chol.info() != Eigen::Success) {
        OperationLog::Error("EeqResult",
            "Cholesky decomposition failed (N=" + std::to_string(N) +
            ") — Coulomb matrix not positive definite");
        return nullptr;
    }

    Eigen::VectorXd u = chol.solve(chi_eff);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd v = chol.solve(ones);

    double denom = ones.dot(v);
    if (std::abs(denom) < 1e-30) {
        OperationLog::Error("EeqResult", "charge constraint denominator is zero");
        return nullptr;
    }

    double lambda = -(Q_total + ones.dot(u)) / denom;
    Eigen::VectorXd q = -(u + lambda * v);

    // Enforce exact charge neutrality.  The block elimination is
    // algebraically exact but the Cholesky solve introduces residual
    // proportional to cond(A).  A uniform shift preserves per-atom
    // charge differences.  Standard practice in EEQ implementations.
    double q_residual = q.sum() - Q_total;
    q.array() -= q_residual / static_cast<double>(N);

    // ── Step 5: store charges with clamp guard ──────────────────────

    auto result = std::make_unique<EeqResult>();
    result->conf_ = &conf;

    int n_clamped = 0;
    for (size_t i = 0; i < N; ++i) {
        auto& atom = conf.MutableAtomAt(i);
        atom.eeq_cn = cn(i);

        double qi = q(i);
        if (std::abs(qi) > charge_clamp) {
            double original = qi;
            qi = (qi > 0) ? charge_clamp : -charge_clamp;
            ++n_clamped;

            // GeometryChoice: traceable decision — charge clamped
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

    double q_sum = 0.0, q_min = q(0), q_max = q(0);
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
