#include "SasaResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"

#include <cmath>

namespace nmr {

// Bondi vdW radii: see PhysicalConstants.h (BondiVdwRadius). Bondi 1964.
static double BondiRadius(Element el) { return BondiVdwRadius(el); }


// Fibonacci lattice on unit sphere (Gonzalez 2010).
static std::vector<Vec3> FibonacciSphere(int n) {
    std::vector<Vec3> points(n);
    double golden = (1.0 + std::sqrt(5.0)) / 2.0;
    for (int i = 0; i < n; ++i) {
        double theta = std::acos(1.0 - 2.0 * (i + 0.5) / n);
        double phi = 2.0 * PI * i / golden;
        points[i] = Vec3(std::sin(theta) * std::cos(phi),
                         std::sin(theta) * std::sin(phi),
                         std::cos(theta));
    }
    return points;
}


std::vector<std::type_index> SasaResult::Dependencies() const {
    return { typeid(SpatialIndexResult) };
}


std::unique_ptr<SasaResult> SasaResult::Compute(ProteinConformation& conf) {
    OperationLog::Scope scope("SasaResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();
    const auto& spatial = conf.Result<SpatialIndexResult>();

    const double probe_radius = CalculatorConfig::Get("sasa_probe_radius");
    const int n_points = static_cast<int>(CalculatorConfig::Get("sasa_n_points"));

    auto unit_sphere = FibonacciSphere(n_points);

    auto result = std::make_unique<SasaResult>();
    result->conf_ = &conf;

    // GeometryChoice: record the parameters used
    GeometryChoiceBuilder choices(conf);

    for (size_t i = 0; i < N; ++i) {
        Vec3 pos_i = conf.PositionAt(i);
        double r_i = BondiRadius(protein.AtomAt(i).element) + probe_radius;
        double sphere_area = 4.0 * PI * r_i * r_i;

        // Find all atoms whose expanded spheres could overlap
        double max_vdw = BondiRadius(Element::S); // largest Bondi radius in table
        double search_radius = r_i + max_vdw + probe_radius;
        auto neighbours = spatial.AtomsWithinRadius(pos_i, search_radius);

        int exposed = 0;
        Vec3 normal_sum = Vec3::Zero();
        for (int p = 0; p < n_points; ++p) {
            Vec3 test_point = pos_i + unit_sphere[p] * r_i;

            bool occluded = false;
            for (size_t j : neighbours) {
                if (j == i) continue;
                double r_j = BondiRadius(protein.AtomAt(j).element) + probe_radius;
                double dist_sq = (test_point - conf.PositionAt(j)).squaredNorm();
                if (dist_sq < r_j * r_j) {
                    occluded = true;
                    break;
                }
            }
            if (!occluded) {
                ++exposed;
                normal_sum += unit_sphere[p];
            }
        }

        auto& ca = conf.MutableAtomAt(i);
        ca.atom_sasa = sphere_area * static_cast<double>(exposed) / n_points;

        // Surface normal: average direction of non-occluded test points.
        // For fully buried atoms (exposed == 0), normal remains zero.
        double normal_mag = normal_sum.norm();
        if (normal_mag > CalculatorConfig::Get("near_zero_vector_norm_threshold"))
            ca.sasa_normal = normal_sum / normal_mag;
        else
            ca.sasa_normal = Vec3::Zero();
    }

    // Record a single GeometryChoice summarising the SASA parameters
    choices.Record(CalculatorId::SASA, 0, "sasa_parameters",
        [probe_radius, n_points, N](GeometryChoice& gc) {
            AddNumber(gc, "probe_radius", probe_radius, "A");
            AddNumber(gc, "n_sphere_points", static_cast<double>(n_points), "count");
            AddNumber(gc, "n_atoms", static_cast<double>(N), "count");
            AddNumber(gc, "radii_source", 0.0, "Bondi_1964");
        });

    OperationLog::Info(LogCalcOther, "SasaResult::Compute",
        std::to_string(N) + " atoms, probe=" +
        std::to_string(probe_radius) + "A, points=" +
        std::to_string(n_points));

    return result;
}


double SasaResult::AtomSASA(size_t atom_index) const {
    return conf_->AtomAt(atom_index).atom_sasa;
}

const std::vector<double>& SasaResult::AllSASA() const {
    // Rebuild from atoms on demand (rare path — WriteFeatures is the normal one)
    thread_local std::vector<double> buf;
    buf.resize(conf_->AtomCount());
    for (size_t i = 0; i < buf.size(); ++i)
        buf[i] = conf_->AtomAt(i).atom_sasa;
    return buf;
}

int SasaResult::WriteFeatures(const ProteinConformation& conf,
                               const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    // atom_sasa: (N,)
    std::vector<double> sasa_data(N);
    for (size_t i = 0; i < N; ++i)
        sasa_data[i] = conf.AtomAt(i).atom_sasa;
    NpyWriter::WriteFloat64(output_dir + "/atom_sasa.npy", sasa_data.data(), N);

    // sasa_normal: (N, 3) — outward surface normal from non-occluded test points
    std::vector<double> normal_data(N * 3);
    for (size_t i = 0; i < N; ++i) {
        const auto& n = conf.AtomAt(i).sasa_normal;
        normal_data[i * 3 + 0] = n.x();
        normal_data[i * 3 + 1] = n.y();
        normal_data[i * 3 + 2] = n.z();
    }
    NpyWriter::WriteFloat64(output_dir + "/sasa_normal.npy", normal_data.data(), N, 3);

    return 2;
}

}  // namespace nmr
