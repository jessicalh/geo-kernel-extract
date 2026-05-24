#include "HydrationGeometryResult.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>

namespace nmr {


std::unique_ptr<HydrationGeometryResult> HydrationGeometryResult::Compute(
        ProteinConformation& conf,
        const SolventEnvironment& solvent) {

    OperationLog::Scope scope("HydrationGeometryResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " waters=" + std::to_string(solvent.WaterCount()));

    if (solvent.Empty()) {
        OperationLog::Error("HydrationGeometryResult",
            "no solvent data — need full-system trajectory");
        return nullptr;
    }

    auto result = std::make_unique<HydrationGeometryResult>();
    result->conf_ = &conf;

    const size_t N = conf.AtomCount();
    const size_t W = solvent.WaterCount();

    // Named constants from TOML (no buried literals).
    const double first_shell_cutoff = CalculatorConfig::Get("water_first_shell_cutoff");
    const double first_sq = first_shell_cutoff * first_shell_cutoff;
    const double near_zero = CalculatorConfig::Get("near_zero_vector_norm_threshold");

    GeometryChoiceBuilder choices(conf);

    // GeometryChoice: one summary record for the parameters used
    choices.Record(CalculatorId::HydrationGeometry, 0, "hydration_geometry_parameters",
        [first_shell_cutoff, N, W](GeometryChoice& gc) {
            AddNumber(gc, "first_shell_cutoff", first_shell_cutoff, "A");
            AddNumber(gc, "n_atoms", static_cast<double>(N), "count");
            AddNumber(gc, "n_waters", static_cast<double>(W), "count");
            AddNumber(gc, "reference_frame", 0.0, "SASA_normal");
        });

    for (size_t ai = 0; ai < N; ++ai) {
        auto& atom = conf.MutableAtomAt(ai);

        // Surface normal from SasaResult (already computed and stored on atom)
        Vec3 normal = atom.sasa_normal;
        atom.water_surface_normal = normal;

        // Accumulate first-shell water dipoles
        Vec3 dipole_sum = Vec3::Zero();
        int n_exposed = 0;
        int n_buried = 0;
        int n_shell = 0;

        Vec3 pos_i = atom.Position();

        for (size_t wi = 0; wi < W; ++wi) {
            Vec3 r = solvent.waters[wi].O_pos - pos_i;
            double d_sq = r.squaredNorm();
            if (d_sq > first_sq) continue;

            ++n_shell;

            // Accumulate water dipole
            dipole_sum += solvent.waters[wi].Dipole();

            // Half-shell: is this water on the exposed side (along surface normal)
            // or the buried side (against normal)?
            // For buried atoms (normal == 0), all waters count as exposed.
            double d = std::sqrt(d_sq);
            Vec3 r_hat = r / d;
            if (normal.norm() > near_zero) {
                double cos_normal = r_hat.dot(normal);
                if (cos_normal > 0)
                    ++n_exposed;  // water is on the outward (solvent) side
                else
                    ++n_buried;   // water is on the interior side
            } else {
                // Buried atom: no meaningful normal, count all as exposed
                ++n_exposed;
            }
        }

        atom.sasa_first_shell_count = n_shell;

        // Half-shell asymmetry: fraction on exposed side
        int n_total = n_exposed + n_buried;
        atom.sasa_half_shell_asymmetry =
            (n_total > 0)
            ? static_cast<double>(n_exposed) / static_cast<double>(n_total)
            : 0.0;

        // Net water dipole vector
        atom.water_dipole_vector = dipole_sum;

        // Dipole alignment: cos(net dipole, surface normal)
        double dip_mag = dipole_sum.norm();
        double norm_mag = normal.norm();
        if (dip_mag > near_zero && norm_mag > near_zero)
            atom.sasa_dipole_alignment = dipole_sum.dot(normal) / (dip_mag * norm_mag);
        else
            atom.sasa_dipole_alignment = 0.0;

        // Dipole coherence: |Σ dᵢ| / n
        // Measures how ordered vs random the first-shell water dipoles are.
        // Fully aligned → coherence ≈ magnitude of one water dipole.
        // Random → coherence → 0 as n → ∞.
        atom.sasa_dipole_coherence =
            (n_shell > 0)
            ? dip_mag / static_cast<double>(n_shell)
            : 0.0;
    }

    OperationLog::Info(LogCalcOther, "HydrationGeometryResult",
        "computed water polarisation for " + std::to_string(N) + " atoms"
        " first_shell=" + std::to_string(first_shell_cutoff) + "A"
        " ref=SASA_normal");

    return result;
}


int HydrationGeometryResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const size_t N = conf.AtomCount();

    // water_polarization: (N, 10)
    // [dipole_x, dipole_y, dipole_z, normal_x, normal_y, normal_z,
    //  asymmetry, alignment, coherence, shell_count]
    std::vector<double> data(N * 10);
    for (size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        data[i * 10 + 0] = a.water_dipole_vector.x();
        data[i * 10 + 1] = a.water_dipole_vector.y();
        data[i * 10 + 2] = a.water_dipole_vector.z();
        data[i * 10 + 3] = a.water_surface_normal.x();
        data[i * 10 + 4] = a.water_surface_normal.y();
        data[i * 10 + 5] = a.water_surface_normal.z();
        data[i * 10 + 6] = a.sasa_half_shell_asymmetry;
        data[i * 10 + 7] = a.sasa_dipole_alignment;
        data[i * 10 + 8] = a.sasa_dipole_coherence;
        data[i * 10 + 9] = static_cast<double>(a.sasa_first_shell_count);
    }
    NpyWriter::WriteFloat64(output_dir + "/water_polarization.npy", data.data(), N, 10);
    return 1;
}

}  // namespace nmr
