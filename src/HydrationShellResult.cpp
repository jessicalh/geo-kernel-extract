#include "HydrationShellResult.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <limits>

namespace nmr {

static constexpr double FIRST_SHELL  = 3.5;   // A
static constexpr double ION_CUTOFF   = 20.0;  // A


std::unique_ptr<HydrationShellResult> HydrationShellResult::Compute(
        ProteinConformation& conf,
        const SolventEnvironment& solvent) {

    OperationLog::Scope scope("HydrationShellResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " waters=" + std::to_string(solvent.WaterCount()) +
        " ions=" + std::to_string(solvent.IonCount()));

    if (solvent.Empty()) {
        OperationLog::Error("HydrationShellResult",
            "no solvent data — need full-system trajectory");
        return nullptr;
    }

    auto result = std::make_unique<HydrationShellResult>();
    result->conf_ = &conf;

    const size_t N = conf.AtomCount();
    const size_t W = solvent.WaterCount();
    const double first_sq = FIRST_SHELL * FIRST_SHELL;

    // Compute protein center of mass (for half-shell direction)
    Vec3 protein_com = Vec3::Zero();
    for (size_t i = 0; i < N; ++i)
        protein_com += conf.AtomAt(i).Position();
    protein_com /= static_cast<double>(N);

    for (size_t ai = 0; ai < N; ++ai) {
        auto& atom = conf.MutableAtomAt(ai);
        Vec3 pos_i = atom.Position();

        // Direction from this atom toward the protein interior
        Vec3 to_interior = (protein_com - pos_i).normalized();

        // Half-shell: count first-shell waters on exposed vs buried side
        int n_exposed = 0;
        int n_buried = 0;
        double dipole_cos_sum = 0.0;
        int dipole_count = 0;

        for (size_t wi = 0; wi < W; ++wi) {
            Vec3 r = solvent.waters[wi].O_pos - pos_i;
            double d_sq = r.squaredNorm();
            if (d_sq > first_sq) continue;

            double d = std::sqrt(d_sq);
            Vec3 r_hat = r / d;

            // Is this water on the interior side or the exposed side?
            double cos_interior = r_hat.dot(to_interior);
            if (cos_interior > 0)
                ++n_buried;
            else
                ++n_exposed;

            // Water dipole orientation relative to atom→water vector
            Vec3 dipole = solvent.waters[wi].Dipole();
            double dip_mag = dipole.norm();
            if (dip_mag > 1e-10) {
                dipole_cos_sum += dipole.dot(r_hat) / dip_mag;
                ++dipole_count;
            }
        }

        int n_total_shell = n_exposed + n_buried;
        atom.half_shell_asymmetry =
            (n_total_shell > 0)
            ? static_cast<double>(n_exposed) / static_cast<double>(n_total_shell)
            : 0.0;

        atom.mean_water_dipole_cos =
            (dipole_count > 0)
            ? dipole_cos_sum / static_cast<double>(dipole_count)
            : 0.0;

        // Nearest ion
        double best_ion_dist = std::numeric_limits<double>::max();
        double best_ion_charge = 0.0;
        for (size_t ii = 0; ii < solvent.IonCount(); ++ii) {
            double d = (solvent.ions[ii].pos - pos_i).norm();
            if (d < best_ion_dist) {
                best_ion_dist = d;
                best_ion_charge = solvent.ions[ii].charge;
            }
        }
        atom.nearest_ion_distance = (best_ion_dist < ION_CUTOFF) ? best_ion_dist : 0.0;
        atom.nearest_ion_charge   = (best_ion_dist < ION_CUTOFF) ? best_ion_charge : 0.0;
    }

    OperationLog::Info(LogCalcOther, "HydrationShellResult",
        "computed hydration shell for " + std::to_string(N) + " atoms");

    return result;
}


int HydrationShellResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const size_t N = conf.AtomCount();

    // hydration_shell: (N, 4)
    // [half_shell_asymmetry, mean_water_dipole_cos, nearest_ion_distance, nearest_ion_charge]
    std::vector<double> data(N * 4);
    for (size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        data[i * 4 + 0] = a.half_shell_asymmetry;
        data[i * 4 + 1] = a.mean_water_dipole_cos;
        data[i * 4 + 2] = a.nearest_ion_distance;
        data[i * 4 + 3] = a.nearest_ion_charge;
    }
    NpyWriter::WriteFloat64(output_dir + "/hydration_shell.npy", data.data(), N, 4);
    return 1;
}

}  // namespace nmr
