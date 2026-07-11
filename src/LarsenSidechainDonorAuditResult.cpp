#include "LarsenSidechainDonorAuditResult.h"

#include "LarsenHBondGeometryCommon.h"
#include "LegacyAmberTopology.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "SemanticEnums.h"
#include "SpatialIndexResult.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <vector>

namespace fs = std::filesystem;

namespace nmr {
namespace {

// Standard donor-heavy/acceptor and H/acceptor geometric audit gates.
// Candidate enumeration is deliberately the heavy-heavy 3.5 A shell;
// the pass flag additionally applies the explicit H...A and D-H...A gates.
constexpr double kCandidateHeavyDistance_A = 3.5;
constexpr double kPassHydrogenDistance_A = 2.5;
constexpr double kPassDonorAngleDeg = 120.0;

std::int32_t ToInt32Index(std::size_t index) {
    if (index == SIZE_MAX ||
        index > static_cast<std::size_t>(
            std::numeric_limits<std::int32_t>::max())) {
        return -1;
    }
    return static_cast<std::int32_t>(index);
}

}  // namespace

std::vector<std::type_index>
LarsenSidechainDonorAuditResult::Dependencies() const {
    return {std::type_index(typeid(SpatialIndexResult))};
}

std::unique_ptr<LarsenSidechainDonorAuditResult>
LarsenSidechainDonorAuditResult::Compute(ProteinConformation& conf) {
    OperationLog::Scope scope(
        "LarsenSidechainDonorAuditResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    auto result = std::make_unique<LarsenSidechainDonorAuditResult>();
    result->donor_atoms_.resize(conf.AtomCount());

    const Protein& protein = conf.ProteinRef();
    const LegacyAmberTopology& topology = protein.LegacyAmber();
    if (!topology.HasAtomSemantic()) {
        OperationLog::Warn("LarsenSidechainDonorAuditResult::Compute",
            "typed atom semantics absent; emitting zero donor audit");
        return result;
    }

    const SpatialIndexResult& spatial = conf.Result<SpatialIndexResult>();

    for (std::size_t donor_idx = 0; donor_idx < conf.AtomCount();
         ++donor_idx) {
        DonorAtomRecord& donor_row = result->donor_atoms_[donor_idx];
        donor_row.residue_index = ToInt32Index(
            protein.AtomAt(donor_idx).residue_index);

        const AtomSemanticTable& donor_sem = topology.SemanticAt(donor_idx);
        donor_row.polar_h_kind =
            static_cast<std::int32_t>(donor_sem.polar_h);

        if (donor_sem.element != Element::H || !donor_sem.IsPolarH()) {
            continue;
        }

        const std::size_t parent_idx =
            larsen_hbond_geometry::BondedHeavyAtomOfHydrogen(
                protein, donor_idx);
        donor_row.parent_atom = ToInt32Index(parent_idx);
        if (parent_idx == SIZE_MAX) continue;

        // Sidechain is decided from the typed parent row.  This excludes
        // backbone amide and N-terminal cap hydrogens even when their own
        // PolarHKind is exchangeable, while retaining every typed sidechain
        // class (amide, indole, ammonium/amine, guanidinium, imidazole,
        // hydroxyl, carboxyl OH, and thiol).
        if (topology.SemanticAt(parent_idx).IsBackbone()) continue;
        donor_row.is_sidechain_polar_h = 1;

        const Vec3 donor_h_pos = conf.PositionAt(donor_idx);
        const Vec3 parent_pos = conf.PositionAt(parent_idx);
        std::vector<std::size_t> nearby = spatial.AtomsWithinRadius(
            parent_pos, kCandidateHeavyDistance_A);
        std::sort(nearby.begin(), nearby.end());

        for (std::size_t acceptor_idx : nearby) {
            if (acceptor_idx == parent_idx || acceptor_idx == donor_idx)
                continue;
            if (protein.AtomAt(acceptor_idx).element != Element::O)
                continue;

            const auto acceptor =
                larsen_hbond_geometry::ClassifyAcceptor(
                    protein, acceptor_idx);
            if (!acceptor) continue;

            CandidateRecord row;
            row.donor_atom = donor_idx;
            row.donor_residue = protein.AtomAt(donor_idx).residue_index;
            row.polar_h_kind =
                static_cast<std::uint8_t>(donor_sem.polar_h);
            row.parent_atom = parent_idx;
            row.acceptor_atom = acceptor_idx;
            row.acceptor_residue =
                protein.AtomAt(acceptor_idx).residue_index;
            row.acceptor_class = acceptor->class_;

            const Vec3 acceptor_pos = conf.PositionAt(acceptor_idx);
            row.h_acceptor_distance_A =
                (donor_h_pos - acceptor_pos).norm();
            row.parent_acceptor_distance_A =
                (parent_pos - acceptor_pos).norm();
            row.angle_parent_h_acceptor_deg =
                larsen_hbond_geometry::AngleDegrees(
                    parent_pos, donor_h_pos, acceptor_pos);
            row.rho_deg = std::numeric_limits<double>::quiet_NaN();
            if (acceptor->C_idx != SIZE_MAX &&
                acceptor->third_idx != SIZE_MAX) {
                const LarsenHBondGeometry geometry =
                    ComputeLarsenHBondGeometry(
                        donor_h_pos,
                        acceptor_pos,
                        conf.PositionAt(acceptor->C_idx),
                        conf.PositionAt(acceptor->third_idx));
                row.rho_deg = geometry.rho_deg;
            }

            row.passes_geometry =
                std::isfinite(row.h_acceptor_distance_A) &&
                std::isfinite(row.parent_acceptor_distance_A) &&
                std::isfinite(row.angle_parent_h_acceptor_deg) &&
                row.parent_acceptor_distance_A <=
                    kCandidateHeavyDistance_A &&
                row.h_acceptor_distance_A <= kPassHydrogenDistance_A &&
                row.angle_parent_h_acceptor_deg >= kPassDonorAngleDeg;
            row.modeled_by_larsen_table2 = false;

            ++donor_row.n_acceptor_candidates_3p5A;
            if (row.passes_geometry) ++donor_row.n_geometry_pass;
            result->candidates_.push_back(row);
        }
    }

    OperationLog::Info(LogCalcOther,
        "LarsenSidechainDonorAuditResult::Compute",
        "candidates=" + std::to_string(result->candidates_.size()));
    return result;
}

int LarsenSidechainDonorAuditResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const std::size_t N = conf.AtomCount();
    std::vector<std::int32_t> atoms(N * 6, 0);
    for (std::size_t i = 0; i < N; ++i) {
        const DonorAtomRecord& row = donor_atoms_[i];
        atoms[i * 6 + 0] = row.is_sidechain_polar_h;
        atoms[i * 6 + 1] = row.polar_h_kind;
        atoms[i * 6 + 2] = row.parent_atom;
        atoms[i * 6 + 3] = row.residue_index;
        atoms[i * 6 + 4] = row.n_acceptor_candidates_3p5A;
        atoms[i * 6 + 5] = row.n_geometry_pass;
    }

    const std::size_t P = candidates_.size();
    std::vector<double> candidates(P * 13, 0.0);
    for (std::size_t i = 0; i < P; ++i) {
        const CandidateRecord& row = candidates_[i];
        double* out = &candidates[i * 13];
        out[0] = static_cast<double>(row.donor_atom);
        out[1] = static_cast<double>(row.donor_residue);
        out[2] = static_cast<double>(row.polar_h_kind);
        out[3] = static_cast<double>(row.parent_atom);
        out[4] = static_cast<double>(row.acceptor_atom);
        out[5] = static_cast<double>(row.acceptor_residue);
        out[6] = static_cast<double>(row.acceptor_class);
        out[7] = row.h_acceptor_distance_A;
        out[8] = row.parent_acceptor_distance_A;
        out[9] = row.angle_parent_h_acceptor_deg;
        out[10] = row.rho_deg;
        out[11] = row.passes_geometry ? 1.0 : 0.0;
        out[12] = 0.0;  // unsupported donor classes are audit-only
    }

    const fs::path dir(output_dir);
    NpyWriter::WriteInt32(
        (dir / "larsen_sidechain_donor_atoms.npy").string(),
        atoms.data(), N, 6);
    NpyWriter::WriteFloat64(
        (dir / "larsen_sidechain_donor_candidates.npy").string(),
        candidates.data(), P, 13);
    return 2;
}

}  // namespace nmr

