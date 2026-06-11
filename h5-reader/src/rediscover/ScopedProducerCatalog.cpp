#include "ScopedProducerCatalog.h"

namespace h5reader::rediscover {

namespace {

QString qsv(std::string_view s) {
    return QString::fromUtf8(s.data(), static_cast<qsizetype>(s.size()));
}

}  // namespace

bool IsExcludedProducerGroup(io::FieldGroup group) {
    switch (group) {
    case io::FieldGroup::LarsenHBond:
    case io::FieldGroup::Tripeptide:
    case io::FieldGroup::Delta:
    case io::FieldGroup::Rediscover:
    case io::FieldGroup::McConnellLegacy:
    case io::FieldGroup::MOPACMcConnellLegacy:
    case io::FieldGroup::MOPACCoulombLegacy:
    case io::FieldGroup::CoulombLegacy:
        return true;
    default:
        return false;
    }
}

bool IsScopedProducerField(const io::FieldSpec& spec) {
    return !IsExcludedProducerGroup(spec.group);
}

bool IsStructuredProducerField(const io::FieldSpec& spec) {
    return spec.kind == io::FieldKind::AtomsCategoryInfo
           || spec.group == io::FieldGroup::Topology;
}

std::size_t NominalComponentCount(const io::FieldSpec& spec) {
    return spec.cols > 0 ? static_cast<std::size_t>(spec.cols) : 1;
}

QString FieldStem(const io::FieldSpec& spec) { return qsv(spec.stem); }

QString FieldGroupName(io::FieldGroup group) {
    switch (group) {
    case io::FieldGroup::Identity: return QStringLiteral("identity");
    case io::FieldGroup::Enrichment: return QStringLiteral("enrichment");
    case io::FieldGroup::ChargeAssignment: return QStringLiteral("charge_assignment");
    case io::FieldGroup::BiotSavart: return QStringLiteral("biot_savart");
    case io::FieldGroup::HaighMallion: return QStringLiteral("haigh_mallion");
    case io::FieldGroup::PiQuadrupole: return QStringLiteral("pi_quadrupole");
    case io::FieldGroup::Dispersion: return QStringLiteral("dispersion");
    case io::FieldGroup::McConnell: return QStringLiteral("mcconnell");
    case io::FieldGroup::Coulomb: return QStringLiteral("coulomb");
    case io::FieldGroup::HBond: return QStringLiteral("hbond");
    case io::FieldGroup::DSSP: return QStringLiteral("dssp");
    case io::FieldGroup::SASA: return QStringLiteral("sasa");
    case io::FieldGroup::WaterField: return QStringLiteral("water_field");
    case io::FieldGroup::Hydration: return QStringLiteral("hydration");
    case io::FieldGroup::WaterPolarization: return QStringLiteral("water_polarization");
    case io::FieldGroup::EEQ: return QStringLiteral("eeq");
    case io::FieldGroup::Gromacs: return QStringLiteral("gromacs");
    case io::FieldGroup::Bonded: return QStringLiteral("bonded");
    case io::FieldGroup::MOPACCore: return QStringLiteral("mopac_core");
    case io::FieldGroup::MOPACCoulomb: return QStringLiteral("mopac_coulomb");
    case io::FieldGroup::APBS: return QStringLiteral("apbs");
    case io::FieldGroup::Orca: return QStringLiteral("orca");
    case io::FieldGroup::AIMNet2: return QStringLiteral("aimnet2");
    case io::FieldGroup::PlanarGeometry: return QStringLiteral("planar_geometry");
    case io::FieldGroup::Topology: return QStringLiteral("topology");
    default: return QStringLiteral("excluded");
    }
}

QString NativeAxisName(io::NativeAxis axis) {
    switch (axis) {
    case io::NativeAxis::Atom: return QStringLiteral("atom");
    case io::NativeAxis::RingContributionPair: return QStringLiteral("ring_contribution_pair");
    case io::NativeAxis::AromaticRing: return QStringLiteral("aromatic_ring");
    case io::NativeAxis::Protein: return QStringLiteral("protein");
    case io::NativeAxis::Bond: return QStringLiteral("bond");
    case io::NativeAxis::MOPACBondNeighborPair: return QStringLiteral("mopac_bond_neighbor_pair");
    case io::NativeAxis::MOPACUniquePair: return QStringLiteral("mopac_unique_pair");
    case io::NativeAxis::Residue: return QStringLiteral("residue");
    case io::NativeAxis::SaturatedRing: return QStringLiteral("saturated_ring");
    case io::NativeAxis::Ring: return QStringLiteral("ring");
    case io::NativeAxis::RingMembership: return QStringLiteral("ring_membership");
    default: return QStringLiteral("excluded");
    }
}

const std::vector<const io::FieldSpec*>& ScopedProducerCatalog() {
    static const std::vector<const io::FieldSpec*> fields = [] {
        std::vector<const io::FieldSpec*> out;
        for (const io::FieldSpec& spec : io::kFieldCatalog) {
            if (IsScopedProducerField(spec)) out.push_back(&spec);
        }
        return out;
    }();
    return fields;
}

}  // namespace h5reader::rediscover
