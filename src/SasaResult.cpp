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

// Bondi vdW radii (Angstroms). Bondi 1964.
static double BondiRadius(Element el) {
    switch (el) {
        case Element::H: return 1.20;
        case Element::C: return 1.70;
        case Element::N: return 1.55;
        case Element::O: return 1.52;
        case Element::S: return 1.80;
        default:         return 1.70;
    }
}


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
    result->sasa_.resize(N, 0.0);

    // GeometryChoice: record the parameters used
    GeometryChoiceBuilder choices(conf);

    for (size_t i = 0; i < N; ++i) {
        Vec3 pos_i = conf.PositionAt(i);
        double r_i = BondiRadius(protein.AtomAt(i).element) + probe_radius;
        double sphere_area = 4.0 * PI * r_i * r_i;

        // Find all atoms whose expanded spheres could overlap
        double max_vdw = 1.80; // S is largest Bondi radius
        double search_radius = r_i + max_vdw + probe_radius;
        auto neighbours = spatial.AtomsWithinRadius(pos_i, search_radius);

        int exposed = 0;
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
            if (!occluded) ++exposed;
        }

        result->sasa_[i] = sphere_area * static_cast<double>(exposed) / n_points;
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
    if (atom_index >= sasa_.size()) return 0.0;
    return sasa_[atom_index];
}


int SasaResult::WriteFeatures(const ProteinConformation& conf,
                               const std::string& output_dir) const {
    NpyWriter::WriteFloat64(output_dir + "/atom_sasa.npy",
                            sasa_.data(), sasa_.size());
    return 1;
}

}  // namespace nmr
