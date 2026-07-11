#include "SidechainCarbonylAnisotropyResult.h"

#include "CalculatorConfig.h"
#include "GeometryResult.h"
#include "McConnellResult.h"
#include "MopacResult.h"
#include "NpyWriter.h"
#include "Protein.h"
#include "ProteinConformation.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace nmr {

namespace sidechain_carbonyl_anisotropy_detail {

namespace {

std::size_t OtherEndpoint(const Bond& bond, std::size_t atom_index) {
    if (bond.atom_index_a == atom_index) return bond.atom_index_b;
    if (bond.atom_index_b == atom_index) return bond.atom_index_a;
    return SIZE_MAX;
}

}  // namespace


SourceBond ClassifySourceBond(const Protein& protein,
                              std::size_t bond_index) {
    SourceBond source;
    source.bond_index = bond_index;
    if (bond_index >= protein.BondCount()) return source;

    const Bond& bond = protein.BondAt(bond_index);
    source.bond_category = bond.category;
    if (bond.category != BondCategory::SidechainCO) return source;
    if (bond.atom_index_a >= protein.AtomCount() ||
        bond.atom_index_b >= protein.AtomCount()) {
        return source;
    }

    const Atom& atom_a = protein.AtomAt(bond.atom_index_a);
    const Atom& atom_b = protein.AtomAt(bond.atom_index_b);
    if (atom_a.element == Element::C && atom_b.element == Element::O) {
        source.carbon_atom = bond.atom_index_a;
        source.oxygen_atom = bond.atom_index_b;
    } else if (atom_a.element == Element::O &&
               atom_b.element == Element::C) {
        source.carbon_atom = bond.atom_index_b;
        source.oxygen_atom = bond.atom_index_a;
    } else {
        return source;
    }

    const Atom& carbon = protein.AtomAt(source.carbon_atom);
    const Atom& oxygen = protein.AtomAt(source.oxygen_atom);
    if (carbon.residue_index >= protein.ResidueCount() ||
        oxygen.residue_index != carbon.residue_index) {
        return source;
    }
    source.residue_index = carbon.residue_index;

    const LegacyAmberTopology& topology = protein.LegacyAmber();
    if (!topology.HasAtomSemantic()) return source;
    const AtomSemanticTable& oxygen_semantic =
        topology.SemanticAt(source.oxygen_atom);
    const AtomSemanticTable& carbon_semantic =
        topology.SemanticAt(source.carbon_atom);
    source.planar_group_kind = oxygen_semantic.planar_group;

    if (oxygen_semantic.IsSidechainAmideOxygen()) {
        source.oxygen_semantic_class =
            OxygenSemanticClass::SidechainAmide;
    } else if (oxygen_semantic.IsSidechainCarboxylateOxygen()) {
        source.oxygen_semantic_class =
            OxygenSemanticClass::SidechainCarboxylate;
    } else {
        return source;
    }

    if (carbon_semantic.planar_group != source.planar_group_kind)
        return source;

    // The canonical plane reference is the lowest-index atom bonded to the
    // carbon that belongs to the same typed planar group, excluding the
    // source oxygen.  This selects amide N for ASN/GLN and the partner O for
    // ASP/GLU carboxylates without parsing atom names.
    for (std::size_t adjacent_bond_index : carbon.bond_indices) {
        if (adjacent_bond_index >= protein.BondCount()) continue;
        const std::size_t other = OtherEndpoint(
            protein.BondAt(adjacent_bond_index), source.carbon_atom);
        if (other == SIZE_MAX || other == source.oxygen_atom ||
            other >= protein.AtomCount()) {
            continue;
        }
        const Atom& other_atom = protein.AtomAt(other);
        if (other_atom.residue_index != source.residue_index) continue;
        if (topology.SemanticAt(other).planar_group !=
            source.planar_group_kind) {
            continue;
        }
        source.plane_reference_atom =
            std::min(source.plane_reference_atom, other);
    }

    source.source_valid = true;
    return source;
}


SourceFrame BuildSourceFrame(const ProteinConformation& conf,
                             const SourceBond& source) {
    SourceFrame frame;
    const Protein& protein = conf.ProteinRef();
    if (source.bond_index < conf.bond_midpoints.size()) {
        frame.origin = conf.bond_midpoints[source.bond_index];
    }
    if (source.bond_index < conf.bond_lengths.size()) {
        frame.bond_length = conf.bond_lengths[source.bond_index];
    }

    if (!source.source_valid ||
        source.carbon_atom >= conf.AtomCount() ||
        source.oxygen_atom >= conf.AtomCount() ||
        source.plane_reference_atom >= conf.AtomCount() ||
        source.carbon_atom >= protein.AtomCount() ||
        source.oxygen_atom >= protein.AtomCount()) {
        return frame;
    }

    const Vec3 carbon_position = conf.PositionAt(source.carbon_atom);
    const Vec3 co = conf.PositionAt(source.oxygen_atom) - carbon_position;
    const Vec3 plane_reference =
        conf.PositionAt(source.plane_reference_atom) - carbon_position;
    const double co_norm = co.norm();
    const Vec3 raw_normal = co.cross(plane_reference);
    frame.normal_norm = raw_normal.norm();

    const double guard =
        CalculatorConfig::Get("near_zero_vector_norm_threshold");
    if (!std::isfinite(co_norm) || co_norm <= guard ||
        !std::isfinite(frame.normal_norm) || frame.normal_norm <= guard) {
        return frame;
    }

    frame.x_axis = co / co_norm;                 // carbon -> oxygen
    frame.z_axis = raw_normal / frame.normal_norm;  // plane normal
    const Vec3 raw_y = frame.z_axis.cross(frame.x_axis);
    const double y_norm = raw_y.norm();
    if (!std::isfinite(y_norm) || y_norm <= guard) return frame;
    frame.y_axis = raw_y / y_norm;

    frame.orthogonality_error = std::max({
        std::abs(frame.x_axis.dot(frame.y_axis)),
        std::abs(frame.x_axis.dot(frame.z_axis)),
        std::abs(frame.y_axis.dot(frame.z_axis)),
        std::abs(frame.x_axis.norm() - 1.0),
        std::abs(frame.y_axis.norm() - 1.0),
        std::abs(frame.z_axis.norm() - 1.0)});
    frame.frame_valid = frame.x_axis.allFinite() &&
                        frame.y_axis.allFinite() &&
                        frame.z_axis.allFinite() &&
                        std::isfinite(frame.orthogonality_error);
    return frame;
}

}  // namespace sidechain_carbonyl_anisotropy_detail


std::vector<std::type_index>
SidechainCarbonylAnisotropyResult::Dependencies() const {
    return {
        std::type_index(typeid(GeometryResult)),
        std::type_index(typeid(McConnellResult)),
    };
}


std::unique_ptr<SidechainCarbonylAnisotropyResult>
SidechainCarbonylAnisotropyResult::Compute(ProteinConformation& conf) {
    using namespace sidechain_carbonyl_anisotropy_detail;

    auto result = std::make_unique<SidechainCarbonylAnisotropyResult>();
    const Protein& protein = conf.ProteinRef();
    const std::size_t atom_count = conf.AtomCount();

    for (std::size_t bond_index = 0;
         bond_index < protein.BondCount(); ++bond_index) {
        if (protein.BondAt(bond_index).category !=
            BondCategory::SidechainCO) {
            continue;
        }
        SourceBond source = ClassifySourceBond(protein, bond_index);
        result->frames_.push_back(BuildSourceFrame(conf, source));
        result->sources_.push_back(std::move(source));
    }

    result->has_mopac_bond_orders_ = conf.HasResult<MopacResult>();
    result->fixed_tensors_.resize(atom_count);
    result->bond_order_tensors_.resize(atom_count);
    result->scalar_audit_.resize(atom_count);

    const std::size_t category = static_cast<std::size_t>(
        McConnellSourceCategory::SidechainCO);
    const std::size_t fixed_channel = static_cast<std::size_t>(
        McConnellChannel::Fixed);
    const std::size_t bond_order_channel = static_cast<std::size_t>(
        McConnellChannel::BondOrder);
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (std::size_t atom_index = 0;
         atom_index < atom_count; ++atom_index) {
        const ConformationAtom& atom = conf.AtomAt(atom_index);
        result->fixed_tensors_[atom_index] =
            atom.mcconnell_source_tensors[category][fixed_channel];

        if (result->has_mopac_bond_orders_) {
            result->bond_order_tensors_[atom_index] =
                atom.mcconnell_source_tensors[category][bond_order_channel];
        } else {
            SphericalTensor missing;
            missing.T0 = nan;
            missing.T1.fill(nan);
            missing.T2.fill(nan);
            result->bond_order_tensors_[atom_index] = missing;
        }

        double nearest_distance = std::numeric_limits<double>::infinity();
        std::size_t accepted_source_count = 0;
        for (const BondNeighbourhood& neighbour : atom.bond_neighbours) {
            if (neighbour.bond_category != BondCategory::SidechainCO)
                continue;
            ++accepted_source_count;
            nearest_distance = std::min(
                nearest_distance, neighbour.distance_to_midpoint);
        }

        result->scalar_audit_[atom_index] = {
            result->fixed_tensors_[atom_index].T2Magnitude(),
            result->has_mopac_bond_orders_
                ? result->bond_order_tensors_[atom_index].T2Magnitude()
                : nan,
            static_cast<double>(accepted_source_count),
            std::isfinite(nearest_distance) ? nearest_distance : nan,
        };
    }

    return result;
}


int SidechainCarbonylAnisotropyResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    (void)conf;
    const std::size_t source_count = sources_.size();
    const std::size_t atom_count = fixed_tensors_.size();

    std::vector<std::int32_t> source_rows(source_count * 8, -1);
    std::vector<double> frame_rows(
        source_count * 12,
        std::numeric_limits<double>::quiet_NaN());
    std::vector<double> quality_rows(
        source_count * 4,
        std::numeric_limits<double>::quiet_NaN());

    auto index_or_missing = [](std::size_t value) -> std::int32_t {
        return value == SIZE_MAX ? -1 : static_cast<std::int32_t>(value);
    };

    for (std::size_t row = 0; row < source_count; ++row) {
        const auto& source = sources_[row];
        const auto& frame = frames_[row];
        std::int32_t* index = &source_rows[row * 8];
        index[0] = index_or_missing(source.bond_index);
        index[1] = index_or_missing(source.carbon_atom);
        index[2] = index_or_missing(source.oxygen_atom);
        index[3] = index_or_missing(source.residue_index);
        index[4] = static_cast<std::int32_t>(source.bond_category);
        index[5] = static_cast<std::int32_t>(source.planar_group_kind);
        index[6] = static_cast<std::int32_t>(
            source.oxygen_semantic_class);
        index[7] = source.source_valid ? 1 : 0;

        double* packed_frame = &frame_rows[row * 12];
        packed_frame[0] = frame.origin.x();
        packed_frame[1] = frame.origin.y();
        packed_frame[2] = frame.origin.z();
        packed_frame[3] = frame.x_axis.x();
        packed_frame[4] = frame.x_axis.y();
        packed_frame[5] = frame.x_axis.z();
        packed_frame[6] = frame.y_axis.x();
        packed_frame[7] = frame.y_axis.y();
        packed_frame[8] = frame.y_axis.z();
        packed_frame[9] = frame.z_axis.x();
        packed_frame[10] = frame.z_axis.y();
        packed_frame[11] = frame.z_axis.z();

        double* quality = &quality_rows[row * 4];
        quality[0] = frame.bond_length;
        quality[1] = frame.orthogonality_error;
        quality[2] = frame.normal_norm;
        quality[3] = frame.frame_valid ? 1.0 : 0.0;
    }

    std::vector<double> fixed_rows(atom_count * 9);
    std::vector<double> bond_order_rows(atom_count * 9);
    std::vector<double> audit_rows(atom_count * 4);
    for (std::size_t atom_index = 0;
         atom_index < atom_count; ++atom_index) {
        fixed_tensors_[atom_index].PackFull9(
            &fixed_rows[atom_index * 9]);
        bond_order_tensors_[atom_index].PackFull9(
            &bond_order_rows[atom_index * 9]);
        std::copy(scalar_audit_[atom_index].begin(),
                  scalar_audit_[atom_index].end(),
                  &audit_rows[atom_index * 4]);
    }

    NpyWriter::WriteInt32(
        output_dir + "/sidechain_co_source_bonds.npy",
        source_rows.data(), source_count, 8);
    NpyWriter::WriteFloat64(
        output_dir + "/sidechain_co_frame.npy",
        frame_rows.data(), source_count, 12);
    NpyWriter::WriteFloat64(
        output_dir + "/sidechain_co_frame_quality.npy",
        quality_rows.data(), source_count, 4);
    NpyWriter::WriteFloat64(
        output_dir + "/sidechain_co_fixed_T2.npy",
        fixed_rows.data(), atom_count, 9);
    NpyWriter::WriteFloat64(
        output_dir + "/sidechain_co_bo_T2.npy",
        bond_order_rows.data(), atom_count, 9);
    NpyWriter::WriteFloat64(
        output_dir + "/sidechain_co_scalar_audit.npy",
        audit_rows.data(), atom_count, 4);
    return 6;
}

}  // namespace nmr
