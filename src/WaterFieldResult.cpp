#include "WaterFieldResult.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "Types.h"

#include <cmath>

namespace nmr {

static void PackST(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

// Cutoffs for water shells and E-field summation.
static constexpr double FIRST_SHELL_CUTOFF  = 3.5;   // Angstroms
static constexpr double SECOND_SHELL_CUTOFF = 5.5;   // Angstroms
static constexpr double EFIELD_CUTOFF       = 15.0;  // Angstroms

// Coulomb constant: 1 / (4π ε₀) in V·Å/e = 14.3996 V·Å/e
static constexpr double KE = 14.3996;


std::unique_ptr<WaterFieldResult> WaterFieldResult::Compute(
        ProteinConformation& conf,
        const SolventEnvironment& solvent) {

    OperationLog::Scope scope("WaterFieldResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " waters=" + std::to_string(solvent.WaterCount()));

    if (solvent.Empty()) {
        OperationLog::Error("WaterFieldResult",
            "no solvent data — need full-system trajectory");
        return nullptr;
    }

    auto result = std::make_unique<WaterFieldResult>();
    result->conf_ = &conf;

    const size_t N = conf.AtomCount();
    const size_t W = solvent.WaterCount();
    const double cutoff_sq = EFIELD_CUTOFF * EFIELD_CUTOFF;
    const double first_sq  = FIRST_SHELL_CUTOFF * FIRST_SHELL_CUTOFF;
    const double second_sq = SECOND_SHELL_CUTOFF * SECOND_SHELL_CUTOFF;

    // For each protein atom, sum E-field from nearby water charges.
    // Water molecule has 3 charge sites: O (q_O), H1 (q_H), H2 (q_H).
    //
    // E_i = Σ_j  q_j * r_ij / |r_ij|³   (V/A)
    // V_ij = q_j * (3 r_ij r_ij^T / |r_ij|⁵  -  I / |r_ij|³)  (V/A²)

    for (size_t ai = 0; ai < N; ++ai) {
        auto& atom = conf.MutableAtomAt(ai);
        Vec3 pos_i = atom.Position();

        Vec3 E_total = Vec3::Zero();
        Mat3 V_total = Mat3::Zero();
        Vec3 E_first = Vec3::Zero();
        Mat3 V_first = Mat3::Zero();
        int n_first = 0;
        int n_second = 0;

        for (size_t wi = 0; wi < W; ++wi) {
            const auto& water = solvent.waters[wi];

            // Quick distance check on oxygen
            Vec3 r_O = water.O_pos - pos_i;
            double d_O_sq = r_O.squaredNorm();
            if (d_O_sq > cutoff_sq) continue;

            double d_O = std::sqrt(d_O_sq);
            bool in_first = (d_O_sq < first_sq);

            // Shell counts (based on O distance)
            if (in_first)
                ++n_first;
            else if (d_O_sq < second_sq)
                ++n_second;

            // Sum contribution from all 3 charge sites
            auto add_charge = [&](const Vec3& q_pos, double q) {
                Vec3 r = q_pos - pos_i;
                double r2 = r.squaredNorm();
                if (r2 < 0.01) return;  // skip self-overlap
                double r_mag = std::sqrt(r2);
                double r3 = r2 * r_mag;
                double r5 = r3 * r2;

                // E-field: E += q * r / r³
                Vec3 E_contrib = (KE * q / r3) * r;
                E_total += E_contrib;

                // EFG: V += q * (3*r*r^T/r⁵ - I/r³)
                Mat3 V_contrib = KE * q * (3.0 * r * r.transpose() / r5
                                           - Mat3::Identity() / r3);
                V_total += V_contrib;

                if (in_first) {
                    E_first += E_contrib;
                    V_first += V_contrib;
                }
            };

            add_charge(water.O_pos,  water.O_charge);
            add_charge(water.H1_pos, water.H_charge);
            add_charge(water.H2_pos, water.H_charge);
        }

        // Make EFG traceless (remove self-potential artifact)
        double trace_total = V_total.trace() / 3.0;
        V_total -= trace_total * Mat3::Identity();
        double trace_first = V_first.trace() / 3.0;
        V_first -= trace_first * Mat3::Identity();

        // Store on atom
        atom.water_efield       = E_total;
        atom.water_efg          = V_total;
        atom.water_efg_spherical = SphericalTensor::Decompose(V_total);
        atom.water_efield_first  = E_first;
        atom.water_efg_first     = V_first;
        atom.water_efg_first_spherical = SphericalTensor::Decompose(V_first);
        atom.water_n_first  = n_first;
        atom.water_n_second = n_second;
    }

    OperationLog::Info(LogCalcOther, "WaterFieldResult",
        "computed E-field + EFG for " + std::to_string(N) +
        " atoms from " + std::to_string(W) + " water molecules");

    return result;
}


int WaterFieldResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const size_t N = conf.AtomCount();
    int n_files = 0;

    // water_efield: (N, 3)
    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const auto& a = conf.AtomAt(i);
            data[i * 3 + 0] = a.water_efield.x();
            data[i * 3 + 1] = a.water_efield.y();
            data[i * 3 + 2] = a.water_efield.z();
        }
        NpyWriter::WriteFloat64(output_dir + "/water_efield.npy", data.data(), N, 3);
        ++n_files;
    }

    // water_efg: (N, 9) SphericalTensor
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i) {
            const auto& st = conf.AtomAt(i).water_efg_spherical;
            PackST(st, &data[i * 9]);
        }
        NpyWriter::WriteFloat64(output_dir + "/water_efg.npy", data.data(), N, 9);
        ++n_files;
    }

    // water_efg_first: (N, 9) SphericalTensor (first shell only)
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i) {
            const auto& st = conf.AtomAt(i).water_efg_first_spherical;
            PackST(st, &data[i * 9]);
        }
        NpyWriter::WriteFloat64(output_dir + "/water_efg_first.npy", data.data(), N, 9);
        ++n_files;
    }

    // water_efield_first: (N, 3) first shell E-field
    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const auto& a = conf.AtomAt(i);
            data[i * 3 + 0] = a.water_efield_first.x();
            data[i * 3 + 1] = a.water_efield_first.y();
            data[i * 3 + 2] = a.water_efield_first.z();
        }
        NpyWriter::WriteFloat64(output_dir + "/water_efield_first.npy", data.data(), N, 3);
        ++n_files;
    }

    // water_shell_counts: (N, 2) [n_first, n_second]
    {
        std::vector<double> data(N * 2);
        for (size_t i = 0; i < N; ++i) {
            const auto& a = conf.AtomAt(i);
            data[i * 2 + 0] = static_cast<double>(a.water_n_first);
            data[i * 2 + 1] = static_cast<double>(a.water_n_second);
        }
        NpyWriter::WriteFloat64(output_dir + "/water_shell_counts.npy", data.data(), N, 2);
        ++n_files;
    }

    return n_files;
}

}  // namespace nmr
