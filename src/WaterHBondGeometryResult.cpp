#include "WaterHBondGeometryResult.h"

#include "CalculatorConfig.h"
#include "LegacyAmberTopology.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "SemanticEnums.h"
#include "SpatialIndexResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace nmr {

namespace water_hbond_geometry_detail {

GeometryEvaluation EvaluateGeometry(const Vec3& donor_heavy,
                                    const Vec3& donor_hydrogen,
                                    const Vec3& acceptor) {
    GeometryEvaluation out;
    out.donor_heavy_acceptor_distance_A =
        (donor_heavy - acceptor).norm();
    out.h_acceptor_distance_A = (donor_hydrogen - acceptor).norm();
    const Vec3 h_to_donor = donor_heavy - donor_hydrogen;
    const Vec3 h_to_acceptor = acceptor - donor_hydrogen;
    const double denom = h_to_donor.norm() * h_to_acceptor.norm();
    if (denom > CalculatorConfig::Get("near_zero_vector_norm_threshold")) {
        const double cosine = std::clamp(
            h_to_donor.dot(h_to_acceptor) / denom, -1.0, 1.0);
        out.angle_deg = std::acos(cosine) * 180.0 / PI;
    } else {
        out.angle_deg = std::numeric_limits<double>::quiet_NaN();
    }
    out.passes_geometry = std::isfinite(out.angle_deg) &&
        out.angle_deg >= 120.0 && out.h_acceptor_distance_A <= 2.5;
    return out;
}

std::uint8_t ProteinAcceptorClass(const AtomSemanticTable& sem) {
    if (sem.formal_charge > 0) return 0;
    if (sem.IsBackboneCarbonylOxygen()) return 1;
    if (sem.IsSidechainAmideOxygen()) return 2;
    if (sem.IsSidechainCarboxylateOxygen()) return 3;
    if (sem.element == Element::O &&
        (sem.planar_group == PlanarGroupKind::AromaticHydroxyl ||
         sem.planar_group == PlanarGroupKind::AromaticOxide ||
         sem.planar_group == PlanarGroupKind::None)) return 4;
    if (sem.element == Element::N &&
        sem.ring_position.primary.position ==
            RingPositionLabel::Heteroatom_NoH) return 5;
    if ((sem.element == Element::N &&
         sem.planar_group == PlanarGroupKind::None) ||
        sem.element == Element::S) return 6;
    return 0;
}

}  // namespace water_hbond_geometry_detail

std::vector<std::type_index> WaterHBondGeometryResult::Dependencies() const {
    return {std::type_index(typeid(SpatialIndexResult))};
}

std::unique_ptr<WaterHBondGeometryResult>
WaterHBondGeometryResult::Compute(ProteinConformation& conf,
                                  const SolventEnvironment& solvent) {
    auto result = std::make_unique<WaterHBondGeometryResult>();
    const Protein& protein = conf.ProteinRef();
    const LegacyAmberTopology& topology = protein.LegacyAmber();
    const size_t N = conf.AtomCount();
    result->counts_.assign(N, std::array<std::int32_t, 6>{0,0,0,0,-1,0});
    const double nan = std::numeric_limits<double>::quiet_NaN();
    result->nearest_.assign(
        N, std::array<double, 8>{nan,nan,nan,0.0,-1.0,0.0,0.0,0.0});

    if (!topology.HasAtomSemantic()) {
        OperationLog::Warn("WaterHBondGeometryResult",
            "typed atom semantics absent; emitting empty candidate table");
        return result;
    }

    const double candidate_cutoff =
        CalculatorConfig::Get("water_first_shell_cutoff");

    auto register_candidate = [&](Candidate candidate) {
        const size_t ai = candidate.protein_atom;
        auto& count = result->counts_[ai];
        const bool water_donates = candidate.mode == 1;
        count[water_donates ? 0 : 2]++;
        if (candidate.geometry.passes_geometry)
            count[water_donates ? 1 : 3]++;

        auto& nearest = result->nearest_[ai];
        nearest[6] += 1.0;
        if (candidate.geometry.passes_geometry) nearest[7] += 1.0;
        if (!std::isfinite(nearest[0]) ||
            candidate.water_O_distance_A < nearest[0]) {
            nearest[0] = candidate.water_O_distance_A;
            nearest[1] = candidate.geometry.h_acceptor_distance_A;
            nearest[2] = candidate.geometry.angle_deg;
            nearest[3] = static_cast<double>(candidate.mode);
            nearest[4] = static_cast<double>(candidate.water_index);
            nearest[5] = candidate.geometry.passes_geometry ? 1.0 : 0.0;
            count[4] = static_cast<std::int32_t>(candidate.water_index);
            count[5] = candidate.geometry.passes_geometry
                ? candidate.mode : 0;
        }
        result->candidates_.push_back(std::move(candidate));
    };

    // Mode 1: water O-H donates to a typed protein acceptor.
    for (size_t ai = 0; ai < N; ++ai) {
        const AtomSemanticTable& sem = topology.SemanticAt(ai);
        const std::uint8_t acceptor_class =
            water_hbond_geometry_detail::ProteinAcceptorClass(sem);
        if (acceptor_class == 0) continue;
        const Vec3 acceptor = conf.PositionAt(ai);
        for (size_t wi = 0; wi < solvent.waters.size(); ++wi) {
            const WaterMolecule& water = solvent.waters[wi];
            if ((water.O_pos - acceptor).norm() > candidate_cutoff) continue;
            const Vec3 chosen_h =
                (water.H1_pos - acceptor).squaredNorm() <=
                (water.H2_pos - acceptor).squaredNorm()
                    ? water.H1_pos : water.H2_pos;
            Candidate row;
            row.protein_atom = ai;
            row.protein_residue = protein.AtomAt(ai).residue_index;
            row.water_index = wi;
            row.mode = 1;
            row.protein_role = acceptor_class;
            row.water_O_distance_A = (water.O_pos - acceptor).norm();
            row.geometry = water_hbond_geometry_detail::EvaluateGeometry(
                water.O_pos, chosen_h, acceptor);
            row.water_O = water.O_pos;
            row.water_H = chosen_h;
            register_candidate(std::move(row));
        }
    }

    // Mode 2: a typed protein polar H donates to water oxygen.
    for (size_t ai = 0; ai < N; ++ai) {
        const AtomSemanticTable& sem = topology.SemanticAt(ai);
        if (!sem.IsPolarH()) continue;
        const size_t parent = protein.AtomAt(ai).parent_atom_index;
        if (parent == SIZE_MAX || parent >= N) continue;
        const Vec3 donor_h = conf.PositionAt(ai);
        const Vec3 donor_heavy = conf.PositionAt(parent);
        for (size_t wi = 0; wi < solvent.waters.size(); ++wi) {
            const WaterMolecule& water = solvent.waters[wi];
            if ((water.O_pos - donor_heavy).norm() > candidate_cutoff)
                continue;
            Candidate row;
            row.protein_atom = ai;
            row.protein_residue = protein.AtomAt(ai).residue_index;
            row.water_index = wi;
            row.mode = 2;
            row.protein_role = static_cast<int>(sem.polar_h);
            row.water_O_distance_A = (water.O_pos - donor_h).norm();
            row.geometry = water_hbond_geometry_detail::EvaluateGeometry(
                donor_heavy, donor_h, water.O_pos);
            row.water_O = water.O_pos;
            row.water_H = Vec3::Constant(nan);  // water is acceptor here
            register_candidate(std::move(row));
        }
    }

    OperationLog::Info(LogCalcOther, "WaterHBondGeometryResult",
        "candidates=" + std::to_string(result->candidates_.size()) +
        " waters=" + std::to_string(solvent.waters.size()));
    return result;
}

int WaterHBondGeometryResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const size_t P = candidates_.size();
    std::vector<double> candidates(P * 16);
    for (size_t i = 0; i < P; ++i) {
        const Candidate& row = candidates_[i];
        double* out = &candidates[i * 16];
        out[0] = static_cast<double>(row.protein_atom);
        out[1] = static_cast<double>(row.protein_residue);
        out[2] = static_cast<double>(row.water_index);
        out[3] = static_cast<double>(row.mode);
        out[4] = static_cast<double>(row.protein_role);
        out[5] = row.water_O_distance_A;
        out[6] = row.geometry.h_acceptor_distance_A;
        out[7] = row.geometry.donor_heavy_acceptor_distance_A;
        out[8] = row.geometry.angle_deg;
        out[9] = row.water_O.x(); out[10] = row.water_O.y();
        out[11] = row.water_O.z(); out[12] = row.water_H.x();
        out[13] = row.water_H.y(); out[14] = row.water_H.z();
        out[15] = row.geometry.passes_geometry ? 1.0 : 0.0;
    }
    std::vector<std::int32_t> counts(conf.AtomCount() * 6);
    std::vector<double> nearest(conf.AtomCount() * 8);
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        std::copy(counts_[i].begin(), counts_[i].end(), counts.begin() + i*6);
        std::copy(nearest_[i].begin(), nearest_[i].end(), nearest.begin() + i*8);
    }
    int written = 0;
    if (NpyWriter::WriteFloat64(output_dir + "/water_hbond_candidates.npy",
                                candidates.data(), P, 16)) ++written;
    if (NpyWriter::WriteInt32(output_dir + "/water_hbond_counts.npy",
                              counts.data(), conf.AtomCount(), 6)) ++written;
    if (NpyWriter::WriteFloat64(output_dir + "/water_hbond_nearest.npy",
                                nearest.data(), conf.AtomCount(), 8)) ++written;
    return written;
}

}  // namespace nmr
