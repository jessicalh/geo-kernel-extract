#include "TrajectorySignalCatalog.h"

#include <QSet>
#include <QStringList>

namespace h5reader::model {

namespace {

UnitSpec unit(UnitDimension dimension,
              const char* sourceSymbol = "",
              const char* displaySymbol = nullptr,
              double scale = 1.0,
              double offset = 0.0,
              bool convertible = true) {
    UnitSpec spec;
    spec.dimension = dimension;
    spec.sourceSymbol = QString::fromLatin1(sourceSymbol);
    spec.displaySymbol = QString::fromLatin1(displaySymbol ? displaySymbol : sourceSymbol);
    spec.scaleToDisplay = scale;
    spec.offsetToDisplay = offset;
    spec.convertible = convertible;
    return spec;
}

ChannelDescriptor channel(const char* id,
                          const char* label,
                          SignalValueShape shape,
                          const UnitSpec& sourceUnits,
                          const UnitSpec& displayUnits = UnitSpec{}) {
    ChannelDescriptor descriptor;
    descriptor.id = QString::fromLatin1(id);
    descriptor.label = QString::fromLatin1(label);
    descriptor.valueShape = shape;
    descriptor.sourceUnits = sourceUnits;
    descriptor.defaultDisplayUnits = displayUnits.sourceSymbol.isEmpty() && displayUnits.displaySymbol.isEmpty()
                                         ? sourceUnits
                                         : displayUnits;
    return descriptor;
}

QVector<ChannelDescriptor> scalarChannels(const UnitSpec& units) {
    return {channel("value", "Value", SignalValueShape::Scalar, units)};
}

QVector<ChannelDescriptor> countChannels() {
    const UnitSpec count = unit(UnitDimension::Count, "count", "count");
    return {channel("count", "Count", SignalValueShape::Count, count)};
}

QVector<ChannelDescriptor> vectorChannels(const UnitSpec& units) {
    return {
        channel("x", "X", SignalValueShape::Scalar, units),
        channel("y", "Y", SignalValueShape::Scalar, units),
        channel("z", "Z", SignalValueShape::Scalar, units),
        channel("magnitude", "Magnitude", SignalValueShape::Scalar, units),
    };
}

QVector<ChannelDescriptor> sphericalTensorChannels(const UnitSpec& units) {
    return {
        channel("T0", "T0", SignalValueShape::Scalar, units),
        channel("T1", "T1", SignalValueShape::TensorComponents, units),
        channel("T2", "T2", SignalValueShape::EfgT2, units),
        channel("component", "Component", SignalValueShape::Scalar, units),
    };
}

QVector<ChannelDescriptor> efgChannels(const UnitSpec& units) {
    return {
        channel("T2", "T2", SignalValueShape::EfgT2, units),
        channel("magnitude", "Magnitude", SignalValueShape::Scalar, units),
        channel("component", "Component", SignalValueShape::Scalar, units),
    };
}

QVector<ChannelDescriptor> angleChannels() {
    const UnitSpec degrees = unit(UnitDimension::Angle, "degrees", "deg");
    return {channel("angle", "Angle", SignalValueShape::Scalar, degrees)};
}

QStringList tensorStripModes() {
    return {
        QStringLiteral("strip.tensor.T0"),
        QStringLiteral("strip.tensor.T1"),
        QStringLiteral("strip.tensor.T2"),
        QStringLiteral("strip.tensor.component"),
        QStringLiteral("strip.spectrum"),
    };
}

QStringList tensorStaticModes() {
    return {
        QStringLiteral("static.tensor"),
        QStringLiteral("static.scalar"),
        QStringLiteral("static.table"),
        QStringLiteral("static.atomColor"),
    };
}

QStringList efgStripModes() {
    return {
        QStringLiteral("strip.tensor.T2"),
        QStringLiteral("strip.tensor.component"),
        QStringLiteral("strip.spectrum"),
    };
}

QStringList efgStaticModes() {
    return {
        QStringLiteral("static.efg"),
        QStringLiteral("static.tensor"),
        QStringLiteral("static.table"),
        QStringLiteral("static.atomColor"),
    };
}

QStringList scalarStripModes() {
    return {
        QStringLiteral("strip.scalar"),
        QStringLiteral("strip.spectrum"),
    };
}

QStringList scalarStaticModes() {
    return {
        QStringLiteral("static.scalar"),
        QStringLiteral("static.table"),
        QStringLiteral("static.atomColor"),
    };
}

QStringList vectorStripModes() {
    return {
        QStringLiteral("strip.vector.component"),
        QStringLiteral("strip.vector.magnitude"),
        QStringLiteral("strip.spectrum"),
    };
}

QStringList vectorStaticModes() {
    return {
        QStringLiteral("static.vector"),
        QStringLiteral("static.vectorGlyph"),
        QStringLiteral("static.table"),
    };
}

QStringList categoryStripModes() {
    return {
        QStringLiteral("strip.category"),
        QStringLiteral("strip.count"),
    };
}

QStringList categoryStaticModes() {
    return {
        QStringLiteral("static.category"),
        QStringLiteral("static.table"),
    };
}

QStringList perClassStripModes() {
    return {
        QStringLiteral("strip.per-class"),
        QStringLiteral("strip.scalar"),
        QStringLiteral("strip.spectrum"),
    };
}

QStringList perClassStaticModes() {
    return {
        QStringLiteral("static.per-class"),
        QStringLiteral("static.table"),
        QStringLiteral("static.atomColor"),
    };
}

QStringList rollupStripModes() {
    return {
        QStringLiteral("strip.rollup"),
    };
}

QStringList rollupStaticModes() {
    return {
        QStringLiteral("static.rollup"),
        QStringLiteral("static.table"),
    };
}

QStringList eventStripModes() {
    return {
        QStringLiteral("strip.event"),
        QStringLiteral("strip.category"),
    };
}

QStringList eventStaticModes() {
    return {
        QStringLiteral("static.event"),
        QStringLiteral("static.table"),
    };
}

bool hasImplementedTemporalSampler(SignalSourceKind sourceKind, const QString& storagePath)
{
    if (storagePath.isEmpty())
        return false;

    switch (sourceKind) {
    case SignalSourceKind::DenseH5Trajectory: {
        static const QSet<QString> kDensePaths = {
            QStringLiteral("/trajectory/positions"),
            QStringLiteral("/trajectory/bs_shielding_time_series"),
            QStringLiteral("/trajectory/hm_shielding_time_series"),
            QStringLiteral("/trajectory/mc_shielding_time_series"),
            QStringLiteral("/trajectory/piquad_shielding_time_series"),
            QStringLiteral("/trajectory/ringchi_shielding_time_series"),
            QStringLiteral("/trajectory/disp_shielding_time_series"),
            QStringLiteral("/trajectory/hbond_shielding_time_series"),
            QStringLiteral("/trajectory/mopac_coulomb_shielding_time_series"),
            QStringLiteral("/trajectory/mopac_mc_shielding_time_series"),
            QStringLiteral("/trajectory/tripeptide_bb_shielding_time_series"),
            QStringLiteral("/trajectory/tripeptide_neighbor_shielding_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_1pHB_shielding_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_1pHaB_shielding_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_2pHB_shielding_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_2pHaB_shielding_time_series"),
            QStringLiteral("/trajectory/sasa_time_series"),
            QStringLiteral("/trajectory/aimnet2_charge_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_count_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_water_term_time_series"),
            QStringLiteral("/trajectory/bonded_energy_time_series"),
            QStringLiteral("/trajectory/mopac_vs_ff14sb_reconciliation"),
            QStringLiteral("/trajectory/water_field_time_series"),
            QStringLiteral("/trajectory/hydration_shell_time_series"),
            QStringLiteral("/trajectory/hydration_geometry_time_series"),
            QStringLiteral("/trajectory/apbs_efield_time_series"),
            QStringLiteral("/trajectory/apbs_efg_time_series"),
            QStringLiteral("/trajectory/aimnet2_embedding_time_series"),
            QStringLiteral("/trajectory/aimnet2_charge_response_gradient_time_series"),
            QStringLiteral("/trajectory/tripeptide_bb_residual_vec_time_series"),
            QStringLiteral("/trajectory/tripeptide_neighbor_residual_vec_prev_time_series"),
            QStringLiteral("/trajectory/tripeptide_neighbor_residual_vec_next_time_series"),
            QStringLiteral("/trajectory/tripeptide_bb_method_tag_time_series"),
            QStringLiteral("/trajectory/dihedral_time_series"),
            QStringLiteral("/trajectory/dssp8_time_series"),
            QStringLiteral("/trajectory/j_coupling_time_series"),
            QStringLiteral("/trajectory/ring_pucker_time_series"),
            QStringLiteral("/trajectory/ring_neighbourhood_trajectory_stats"),
            QStringLiteral("/trajectory/gromacs_energy_time_series"),
            QStringLiteral("/trajectory/rmsd_tracking"),
            QStringLiteral("/trajectory/bond_length_stats"),
            QStringLiteral("/trajectory/bs_welford"),
            QStringLiteral("/trajectory/hm_welford"),
            QStringLiteral("/trajectory/mc_welford"),
            QStringLiteral("/trajectory/sasa_welford"),
            QStringLiteral("/trajectory/eeq_welford"),
            QStringLiteral("/trajectory/hbond_count_welford"),
            QStringLiteral("/trajectory/mopac_charge_welford"),
            QStringLiteral("/trajectory/mopac_bond_order_welford"),
            QStringLiteral("/trajectory/water_field_welford"),
            QStringLiteral("/trajectory/aimnet2_charge_response_gradient_welford"),
            QStringLiteral("/trajectory/hydration_shell_welford"),
            QStringLiteral("/trajectory/hydration_geometry_welford"),
            QStringLiteral("/trajectory/bs_t0_autocorrelation"),
            QStringLiteral("/trajectory/ired_order_parameters"),
            QStringLiteral("/trajectory/kernel_dynamics"),
            QStringLiteral("/trajectory/reorientational_dynamics"),
            QStringLiteral("/trajectory/dihedral_autocorrelation"),
            QStringLiteral("/trajectory/kernel_coherence"),
            QStringLiteral("/trajectory/dssp8_transition"),
            QStringLiteral("/trajectory/dihedral_bin_transition"),
        };
        return kDensePaths.contains(storagePath);
    }
    case SignalSourceKind::FrameNpySnapshot:
        return true;
    case SignalSourceKind::OrcaDftFrame:
        return storagePath.startsWith(QStringLiteral("orca_"));
    case SignalSourceKind::Topology:
        return storagePath == QStringLiteral("bond_length");
    case SignalSourceKind::DerivedGeometry:
        return storagePath == QStringLiteral("distance")
               || storagePath == QStringLiteral("angle")
               || storagePath == QStringLiteral("dihedral")
               || storagePath == QStringLiteral("atom_displacement");
    case SignalSourceKind::SelectionEvents:
        return storagePath == QStringLiteral("/trajectory/selections")
               || storagePath == QStringLiteral("selection_timeline")
               || storagePath == QStringLiteral("selection_counts");
    }
    return false;
}

SignalDescriptor makeDescriptor(const char* id,
                                const char* conceptKey,
                                SignalSourceKind sourceKind,
                                const char* importSet,
                                const char* family,
                                const char* label,
                                SourceResidency residency,
                                SignalAxis nativeAxis,
                                SignalAxis requiredAnchor,
                                SignalValueShape valueShape,
                                const UnitSpec& units,
                                const QStringList& temporalModes,
                                const QStringList& staticModes,
                                const QVector<ChannelDescriptor>& channels,
                                const char* storagePath = "",
                                bool sourceAttachedMask = false,
                                bool frameLocalPayload = false,
                                SampleStatus samplingStatus = SampleStatus::Gap,
                                GapReason samplingGapReason = GapReason::Pending) {
    SignalDescriptor descriptor;
    descriptor.id = QString::fromLatin1(id);
    descriptor.conceptKey = QString::fromLatin1(conceptKey);
    descriptor.sourceKind = sourceKind;
    descriptor.importSet = QString::fromLatin1(importSet);
    descriptor.family = QString::fromLatin1(family);
    descriptor.label = QString::fromLatin1(label);
    descriptor.storagePath = QString::fromLatin1(storagePath);
    descriptor.sourceUnits = units;
    descriptor.defaultDisplayUnits = units;
    descriptor.residency = residency;
    descriptor.nativeAxis = nativeAxis;
    descriptor.requiredAnchor = requiredAnchor;
    descriptor.valueShape = valueShape;
    descriptor.temporalModes = temporalModes;
    descriptor.staticModes = staticModes;
    descriptor.temporal = !temporalModes.isEmpty();
    descriptor.staticDisplay = !staticModes.isEmpty();
    descriptor.sourceAttachedMask = sourceAttachedMask;
    descriptor.frameLocalPayload = frameLocalPayload;
    descriptor.finiteScalarRequired = valueShape != SignalValueShape::Category
                                      && valueShape != SignalValueShape::EventRecord
                                      && valueShape != SignalValueShape::Embedding;
    descriptor.samplingStatus = samplingStatus;
    descriptor.samplingGapReason = samplingGapReason;
    if (descriptor.samplingStatus == SampleStatus::Gap
        && descriptor.samplingGapReason == GapReason::Pending) {
        if (!descriptor.temporal) {
            descriptor.samplingStatus = SampleStatus::NotAvailable;
            descriptor.samplingGapReason = GapReason::NotApplicable;
        } else if (hasImplementedTemporalSampler(descriptor.sourceKind, descriptor.storagePath)) {
            descriptor.samplingStatus = SampleStatus::Valid;
            descriptor.samplingGapReason = GapReason::None;
        }
    }
    descriptor.channels = channels;
    descriptor.tags = {
        descriptor.importSet,
        descriptor.family,
        ToString(sourceKind),
        ToString(nativeAxis),
        ToString(valueShape),
    };
    return descriptor;
}

void add(QVector<SignalDescriptor>& descriptors, const SignalDescriptor& descriptor) {
    descriptors.push_back(descriptor);
}

void addDenseH5(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec none = unit(UnitDimension::Dimensionless, "", "");
    const UnitSpec length = unit(UnitDimension::Length, "A", "A");
    const UnitSpec shielding = unit(UnitDimension::MagneticShielding, "ppm", "ppm");
    const UnitSpec efield = unit(UnitDimension::ElectricField, "V/A", "V/A");
    const UnitSpec efg = unit(UnitDimension::ElectricFieldGradient, "V/A^2", "V/A^2");
    const UnitSpec charge = unit(UnitDimension::Charge, "e", "e");
    const UnitSpec energy = unit(UnitDimension::Energy, "kJ/mol", "kJ/mol");
    const UnitSpec angle = unit(UnitDimension::Angle, "radians", "deg", 57.29577951308232);
    const UnitSpec tag = unit(UnitDimension::Tag, "tag", "tag", 1.0, 0.0, false);

    add(descriptors,
        makeDescriptor("h5:positions",
                       "positions",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "identity",
                       "Atom positions",
                       SourceResidency::StartupLoaded,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::Vector3,
                       length,
                       vectorStripModes(),
                       {QStringLiteral("static.geometry"), QStringLiteral("static.topology"), QStringLiteral("static.table")},
                       vectorChannels(length),
                       "/trajectory/positions"));

    const struct TensorGroup {
        const char* id;
        const char* conceptKey;
        const char* family;
        const char* label;
        const char* path;
        SignalValueShape shape;
    } tensorGroups[] = {
        {"h5:bs_shielding_time_series", "bs_shielding", "biot_savart", "Biot-Savart shielding time series", "/trajectory/bs_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:hm_shielding_time_series", "hm_shielding", "haigh_mallion", "Haigh-Mallion shielding time series", "/trajectory/hm_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:mc_shielding_time_series", "mc_shielding", "mcconnell", "McConnell shielding time series", "/trajectory/mc_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:piquad_shielding_time_series", "pq_shielding", "pi_quadrupole", "Pi quadrupole shielding time series", "/trajectory/piquad_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:ringchi_shielding_time_series", "ringchi_shielding", "ring_susceptibility", "Ring susceptibility shielding time series", "/trajectory/ringchi_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:disp_shielding_time_series", "disp_shielding", "dispersion", "Dispersion shielding time series", "/trajectory/disp_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:hbond_shielding_time_series", "hbond_shielding", "hbond", "H-bond shielding time series", "/trajectory/hbond_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:mopac_coulomb_shielding_time_series", "mopac_coulomb_shielding", "mopac_coulomb", "MOPAC Coulomb T2 shielding time series", "/trajectory/mopac_coulomb_shielding_time_series", SignalValueShape::EfgT2},
        {"h5:mopac_mc_shielding_time_series", "mopac_mc_shielding", "mopac_mcconnell", "MOPAC McConnell shielding time series", "/trajectory/mopac_mc_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:tripeptide_bb_shielding_time_series", "tripeptide_bb_shielding", "tripeptide", "Tripeptide backbone shielding time series", "/trajectory/tripeptide_bb_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:tripeptide_neighbor_shielding_time_series", "tripeptide_neighbor_shielding", "tripeptide", "Tripeptide neighbor shielding time series", "/trajectory/tripeptide_neighbor_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:larsen_hbond_1pHB_shielding_time_series", "larsen_hbond_1pHB_shielding", "larsen_hbond", "Larsen 1pHB shielding time series", "/trajectory/larsen_hbond_1pHB_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:larsen_hbond_1pHaB_shielding_time_series", "larsen_hbond_1pHaB_shielding", "larsen_hbond", "Larsen 1pHaB shielding time series", "/trajectory/larsen_hbond_1pHaB_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:larsen_hbond_2pHB_shielding_time_series", "larsen_hbond_2pHB_shielding", "larsen_hbond", "Larsen 2pHB shielding time series", "/trajectory/larsen_hbond_2pHB_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:larsen_hbond_2pHaB_shielding_time_series", "larsen_hbond_2pHaB_shielding", "larsen_hbond", "Larsen 2pHaB shielding time series", "/trajectory/larsen_hbond_2pHaB_shielding_time_series", SignalValueShape::SphericalTensor},
    };

    for (const TensorGroup& group : tensorGroups) {
        const bool efgOnly = group.shape == SignalValueShape::EfgT2;
        add(descriptors,
            makeDescriptor(group.id,
                           group.conceptKey,
                           SignalSourceKind::DenseH5Trajectory,
                           "TrajectoryH5",
                           group.family,
                           group.label,
                           SourceResidency::StartupLoaded,
                           SignalAxis::Atom,
                           SignalAxis::Atom,
                           group.shape,
                           shielding,
                           efgOnly ? efgStripModes() : tensorStripModes(),
                           efgOnly ? efgStaticModes() : tensorStaticModes(),
                           efgOnly ? efgChannels(shielding) : sphericalTensorChannels(shielding),
                           group.path,
                           true));
    }

    add(descriptors, makeDescriptor("h5:sasa_time_series", "atom_sasa", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "sasa", "SASA time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, unit(UnitDimension::Length, "A^2", "A^2"), scalarStripModes(), scalarStaticModes(), scalarChannels(unit(UnitDimension::Length, "A^2", "A^2")), "/trajectory/sasa_time_series", true));
    add(descriptors, makeDescriptor("h5:aimnet2_charge_time_series", "aimnet2_charges", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "aimnet2", "AIMNet2 charge time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, charge, scalarStripModes(), scalarStaticModes(), scalarChannels(charge), "/trajectory/aimnet2_charge_time_series", true));
    add(descriptors, makeDescriptor("h5:larsen_hbond_count_time_series", "larsen_hbond_count", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "larsen_hbond", "Larsen H-bond count time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.scalar")}, {QStringLiteral("static.scalar"), QStringLiteral("static.atomColor"), QStringLiteral("static.table")}, countChannels(), "/trajectory/larsen_hbond_count_time_series", true));
    add(descriptors, makeDescriptor("h5:larsen_hbond_water_term_time_series", "larsen_hbond_water_term", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "larsen_hbond", "Larsen H-bond water term time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, shielding, scalarStripModes(), scalarStaticModes(), scalarChannels(shielding), "/trajectory/larsen_hbond_water_term_time_series", true));
    add(descriptors, makeDescriptor("h5:bonded_energy_time_series", "bonded_energy", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "bonded", "Bonded energy time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, energy, scalarStripModes(), scalarStaticModes(), scalarChannels(energy), "/trajectory/bonded_energy_time_series", true));
    add(descriptors, makeDescriptor("h5:mopac_vs_ff14sb_reconciliation", "mopac_vs_ff14sb_reconciliation", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "mopac_core", "MOPAC vs FF14SB reconciliation", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, energy, scalarStripModes(), {QStringLiteral("static.scalar"), QStringLiteral("static.table")}, scalarChannels(energy), "/trajectory/mopac_vs_ff14sb_reconciliation", true));

    add(descriptors, makeDescriptor("h5:water_field_efield_time_series", "water_efield", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "water_field", "Water electric field time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield), "/trajectory/water_field_time_series", true));
    add(descriptors, makeDescriptor("h5:water_field_efg_time_series", "water_efg", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "water_field", "Water EFG time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg), "/trajectory/water_field_time_series", true));
    add(descriptors, makeDescriptor("h5:water_shell_count_time_series", "water_shell_counts", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "water_field", "Water shell counts", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.scalar")}, {QStringLiteral("static.scalar"), QStringLiteral("static.atomColor"), QStringLiteral("static.table")}, countChannels(), "/trajectory/water_field_time_series", true));
    add(descriptors, makeDescriptor("h5:hydration_shell_time_series", "hydration_shell", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "hydration", "Hydration shell time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {}, "/trajectory/hydration_shell_time_series", true));
    add(descriptors, makeDescriptor("h5:hydration_geometry_time_series", "hydration_geometry", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "hydration", "Hydration geometry time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {}, "/trajectory/hydration_geometry_time_series", true));
    add(descriptors, makeDescriptor("h5:apbs_efield_time_series", "apbs_E", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "apbs", "APBS electric field time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield), "/trajectory/apbs_efield_time_series", true));
    add(descriptors, makeDescriptor("h5:apbs_efg_time_series", "apbs_efg", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "apbs", "APBS EFG time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg), "/trajectory/apbs_efg_time_series", true));
    add(descriptors, makeDescriptor("h5:aimnet2_embedding_time_series", "aimnet2_embedding", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "aimnet2", "AIMNet2 embedding time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Embedding, none, {QStringLiteral("strip.embedding.component")}, {QStringLiteral("static.embedding"), QStringLiteral("static.table")}, {}, "/trajectory/aimnet2_embedding_time_series", true));
    add(descriptors, makeDescriptor("h5:aimnet2_charge_response_gradient_time_series", "aimnet2_charge_response_gradient", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "aimnet2", "AIMNet2 charge-response gradient", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, charge, vectorStripModes(), vectorStaticModes(), vectorChannels(charge), "/trajectory/aimnet2_charge_response_gradient_time_series", true));
    add(descriptors, makeDescriptor("h5:tripeptide_bb_residual_vec_time_series", "tripeptide_bb_residual_vec", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "tripeptide", "Tripeptide backbone residual vector", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, shielding, vectorStripModes(), vectorStaticModes(), vectorChannels(shielding), "/trajectory/tripeptide_bb_residual_vec_time_series", true));
    add(descriptors, makeDescriptor("h5:tripeptide_neighbor_residual_vec_prev_time_series", "tripeptide_neighbor_residual_vec_prev", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "tripeptide", "Tripeptide previous-neighbor residual vector", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, shielding, vectorStripModes(), vectorStaticModes(), vectorChannels(shielding), "/trajectory/tripeptide_neighbor_residual_vec_prev_time_series", true));
    add(descriptors, makeDescriptor("h5:tripeptide_neighbor_residual_vec_next_time_series", "tripeptide_neighbor_residual_vec_next", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "tripeptide", "Tripeptide next-neighbor residual vector", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, shielding, vectorStripModes(), vectorStaticModes(), vectorChannels(shielding), "/trajectory/tripeptide_neighbor_residual_vec_next_time_series", true));
    add(descriptors, makeDescriptor("h5:tripeptide_bb_method_tag_time_series", "tripeptide_bb_method_tag", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "tripeptide", "Tripeptide method tag", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Category, tag, categoryStripModes(), categoryStaticModes(), {channel("method", "Method", SignalValueShape::Category, tag)}, "/trajectory/tripeptide_bb_method_tag_time_series", true));

    add(descriptors, makeDescriptor("h5:dihedral_time_series", "residue_dihedral", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "planar_geometry", "Residue dihedral time series", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::PerClassBlock, angle, perClassStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {channel("phi", "Phi", SignalValueShape::Scalar, angle), channel("psi", "Psi", SignalValueShape::Scalar, angle), channel("omega", "Omega", SignalValueShape::Scalar, angle), channel("chi", "Chi", SignalValueShape::Scalar, angle)}, "/trajectory/dihedral_time_series", true));
    add(descriptors, makeDescriptor("h5:dssp8_time_series", "dssp_ss8", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "dssp", "DSSP8 residue state time series", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.category"), QStringLiteral("static.table")}, {channel("ss8", "SS8", SignalValueShape::Category, tag)}, "/trajectory/dssp8_time_series", true));
    add(descriptors, makeDescriptor("h5:j_coupling_time_series", "j_coupling", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "j_coupling", "J-coupling time series", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::PerClassBlock, unit(UnitDimension::Frequency, "Hz", "Hz"), perClassStripModes(), {QStringLiteral("static.per-class"), QStringLiteral("static.table")}, {}, "/trajectory/j_coupling_time_series", true));
    add(descriptors, makeDescriptor("h5:ring_pucker_time_series", "ring_pucker", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "planar_geometry", "Ring pucker time series", SourceResidency::StartupLoaded, SignalAxis::Ring, SignalAxis::Ring, SignalValueShape::PerClassBlock, angle, perClassStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {}, "/trajectory/ring_pucker_time_series", true));
    add(descriptors, makeDescriptor("h5:ring_neighbourhood_trajectory_stats", "ring_neighbourhood", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "ring_current", "Ring neighbourhood trajectory stats", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::PerClassBlock, length, perClassStripModes(), {QStringLiteral("static.per-class"), QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {channel("distance", "Distance", SignalValueShape::Scalar, length), channel("rho", "Rho", SignalValueShape::Scalar, length), channel("z", "Z", SignalValueShape::Scalar, length), channel("in_plane_angle", "In-plane angle", SignalValueShape::Scalar, angle)}, "/trajectory/ring_neighbourhood_trajectory_stats", true));
    add(descriptors, makeDescriptor("h5:gromacs_energy_time_series", "gromacs_energy", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "gromacs", "Gromacs energy/runtime time series", SourceResidency::StartupLoaded, SignalAxis::System, SignalAxis::System, SignalValueShape::PerClassBlock, none, {QStringLiteral("strip.system"), QStringLiteral("strip.per-class"), QStringLiteral("strip.tensor.component")}, {QStringLiteral("static.system"), QStringLiteral("static.table")}, {channel("energy", "Energy", SignalValueShape::Scalar, energy), channel("temperature", "Temperature", SignalValueShape::Scalar, unit(UnitDimension::Temperature, "K", "K")), channel("pressure", "Pressure", SignalValueShape::Scalar, unit(UnitDimension::Pressure, "bar", "bar")), channel("volume", "Volume", SignalValueShape::Scalar, unit(UnitDimension::Volume, "nm^3", "nm^3"))}, "/trajectory/gromacs_energy_time_series", true));
    add(descriptors, makeDescriptor("h5:rmsd_tracking", "rmsd_tracking", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "rmsd", "RMSD tracking", SourceResidency::StartupLoaded, SignalAxis::System, SignalAxis::System, SignalValueShape::Scalar, length, {QStringLiteral("strip.system"), QStringLiteral("strip.scalar"), QStringLiteral("strip.event")}, {QStringLiteral("static.system"), QStringLiteral("static.table")}, scalarChannels(length), "/trajectory/rmsd_tracking", true));

    const struct RollupGroup {
        const char* id;
        const char* conceptKey;
        const char* family;
        const char* label;
        const char* path;
        SignalAxis axis;
        UnitSpec units;
    } rollups[] = {
        {"h5:bond_length_stats", "bond_length.stats", "topology", "Bond length statistics", "/trajectory/bond_length_stats", SignalAxis::Bond, length},
        {"h5:bs_welford", "bs_shielding.stats", "biot_savart", "Biot-Savart Welford rollup", "/trajectory/bs_welford", SignalAxis::Atom, shielding},
        {"h5:hm_welford", "hm_shielding.stats", "haigh_mallion", "Haigh-Mallion Welford rollup", "/trajectory/hm_welford", SignalAxis::Atom, shielding},
        {"h5:mc_welford", "mc_shielding.stats", "mcconnell", "McConnell Welford rollup", "/trajectory/mc_welford", SignalAxis::Atom, shielding},
        {"h5:sasa_welford", "atom_sasa.stats", "sasa", "SASA Welford rollup", "/trajectory/sasa_welford", SignalAxis::Atom, unit(UnitDimension::Length, "A^2", "A^2")},
        {"h5:eeq_welford", "eeq_charges.stats", "eeq", "EEQ charge Welford rollup", "/trajectory/eeq_welford", SignalAxis::Atom, charge},
        {"h5:hbond_count_welford", "hbond_count.stats", "hbond", "H-bond count Welford rollup", "/trajectory/hbond_count_welford", SignalAxis::Atom, unit(UnitDimension::Count, "count", "count")},
        {"h5:mopac_charge_welford", "mopac_charges.stats", "mopac_core", "MOPAC charge Welford rollup", "/trajectory/mopac_charge_welford", SignalAxis::Atom, charge},
        {"h5:mopac_bond_order_welford", "mopac_bond_orders.stats", "mopac_core", "MOPAC bond order Welford rollup", "/trajectory/mopac_bond_order_welford", SignalAxis::Bond, none},
        {"h5:water_field_welford", "water_field.stats", "water_field", "Water field Welford rollup", "/trajectory/water_field_welford", SignalAxis::Atom, efield},
        {"h5:aimnet2_charge_response_gradient_welford", "aimnet2_charge_response_gradient.stats", "aimnet2", "AIMNet2 CRG Welford rollup", "/trajectory/aimnet2_charge_response_gradient_welford", SignalAxis::Atom, charge},
        {"h5:hydration_shell_welford", "hydration_shell.stats", "hydration", "Hydration shell Welford rollup", "/trajectory/hydration_shell_welford", SignalAxis::Atom, none},
        {"h5:hydration_geometry_welford", "hydration_geometry.stats", "hydration", "Hydration geometry Welford rollup", "/trajectory/hydration_geometry_welford", SignalAxis::Atom, none},
        {"h5:bs_t0_autocorrelation", "bs_shielding.autocorrelation", "biot_savart", "Biot-Savart T0 autocorrelation", "/trajectory/bs_t0_autocorrelation", SignalAxis::Atom, shielding},
    };

    for (const RollupGroup& group : rollups) {
        add(descriptors,
            makeDescriptor(group.id,
                           group.conceptKey,
                           SignalSourceKind::DenseH5Trajectory,
                           "TrajectoryH5",
                           group.family,
                           group.label,
                           SourceResidency::StartupLoaded,
                           group.axis,
                           group.axis,
                           SignalValueShape::RollupMoments,
                           group.units,
                           rollupStripModes(),
                           rollupStaticModes(),
                           {},
                           group.path,
                           true));
    }

    // ── Per-atom × per-channel (KernelDynamics) ─────────────────────
    // The 13 kernel channels are stable across the codebase (producer
    // emits them in canonical thesis-narrative order). Catalog declares
    // them as ChannelDescriptors so the strip-mode dispatch + the
    // panel renderers can index by channel id.
    static const struct { const char* id; const char* label; } kKernelChannels[] = {
        {"bs_T0", "BS T0"},     {"bs_absT2", "BS |T2|"},
        {"hm_T0", "HM T0"},     {"hm_absT2", "HM |T2|"},
        {"mc_T0", "MC T0"},     {"mc_absT2", "MC |T2|"},
        {"ringchi_T0", "RingChi T0"}, {"ringchi_absT2", "RingChi |T2|"},
        {"hbond_T0", "H-bond T0"},    {"hbond_absT2", "H-bond |T2|"},
        {"piquad_absT2", "PiQuad |T2|"},
        {"disp_absT2", "Disp |T2|"},
        {"apbs_absT2", "APBS |T2|"},
    };
    QVector<ChannelDescriptor> kernelChannels;
    for (const auto& ch : kKernelChannels) {
        kernelChannels.push_back(channel(ch.id, ch.label, SignalValueShape::Scalar, none));
    }

    // Two curve-shaped descriptors (ACF over lag, PSD over frequency) +
    // three per-channel scalar reductions all sharing one storagePath.
    // The denseH5Plan branch dispatches by descriptor.conceptKey.
    add(descriptors,
        makeDescriptor("h5:kernel_dynamics_acf",
                       "kernel_dynamics.acf",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "kernel_dynamics",
                       "Kernel autocorrelation (ACF)",
                       SourceResidency::StartupLoaded,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::CurveOverLag,
                       none,
                       {},                       // not a temporal strip
                       {QStringLiteral("static.curve.lag.animated")},
                       kernelChannels,
                       "/trajectory/kernel_dynamics",
                       true,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));
    add(descriptors,
        makeDescriptor("h5:kernel_dynamics_psd",
                       "kernel_dynamics.psd",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "kernel_dynamics",
                       "Kernel power spectrum (PSD)",
                       SourceResidency::StartupLoaded,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::Spectrum,
                       none,
                       {},
                       {QStringLiteral("static.spectrum.power")},
                       kernelChannels,
                       "/trajectory/kernel_dynamics",
                       true,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));
    auto addKernelScalar = [&](const char* id, const char* conceptKey,
                               const char* label, const UnitSpec& units) {
        add(descriptors,
            makeDescriptor(id,
                           conceptKey,
                           SignalSourceKind::DenseH5Trajectory,
                           "TrajectoryH5",
                           "kernel_dynamics",
                           label,
                           SourceResidency::StartupLoaded,
                           SignalAxis::Atom,
                           SignalAxis::Atom,
                           SignalValueShape::PerClassBlock,
                           units,
                           perClassStripModes(),
                           perClassStaticModes(),
                           kernelChannels,
                           "/trajectory/kernel_dynamics",
                           true));
    };
    addKernelScalar("h5:kernel_dynamics_decay_time",
                    "kernel_dynamics.decay_time",
                    "Kernel decay time (per channel, ps)",
                    unit(UnitDimension::Time, "ps", "ps"));
    addKernelScalar("h5:kernel_dynamics_peak_freq",
                    "kernel_dynamics.peak_freq",
                    "Kernel peak frequency (per channel, 1/ps)",
                    unit(UnitDimension::Frequency, "1/ps", "1/ps"));
    addKernelScalar("h5:kernel_dynamics_spectral_centroid",
                    "kernel_dynamics.spectral_centroid",
                    "Kernel spectral centroid (per channel, 1/ps)",
                    unit(UnitDimension::Frequency, "1/ps", "1/ps"));

    // ── KernelCoherence (atom × matrix) ─────────────────────────────
    add(descriptors,
        makeDescriptor("h5:kernel_coherence",
                       "kernel_coherence.matrix",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "kernel_coherence",
                       "Kernel coherence matrix (per atom, 13×13 Pearson)",
                       SourceResidency::StartupLoaded,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::Matrix,
                       none,
                       {},
                       {QStringLiteral("static.chord.coupling"),
                        QStringLiteral("static.table")},
                       kernelChannels,
                       "/trajectory/kernel_coherence",
                       true,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));

    // ── Residue-axis dihedral autocorrelation (phi/psi only for v1) ─
    auto addDihedralScalar = [&](const char* id, const char* conceptKey,
                                  const char* label) {
        add(descriptors,
            makeDescriptor(id, conceptKey,
                           SignalSourceKind::DenseH5Trajectory,
                           "TrajectoryH5",
                           "dihedral_autocorrelation",
                           label,
                           SourceResidency::StartupLoaded,
                           SignalAxis::Residue,
                           SignalAxis::Residue,
                           SignalValueShape::Scalar,
                           unit(UnitDimension::Time, "ps", "ps"),
                           {},
                           {QStringLiteral("static.bar.sequence"),
                            QStringLiteral("static.table")},
                           scalarChannels(unit(UnitDimension::Time, "ps", "ps")),
                           "/trajectory/dihedral_autocorrelation",
                           true,
                           false,
                           SampleStatus::Valid,
                           GapReason::None));
    };
    addDihedralScalar("h5:dihedral_phi_corr_time", "dihedral.phi_corr_time",
                      "phi torsional decorrelation time (ps)");
    addDihedralScalar("h5:dihedral_psi_corr_time", "dihedral.psi_corr_time",
                      "psi torsional decorrelation time (ps)");

    auto addDihedralCurve = [&](const char* id, const char* conceptKey,
                                 const char* label) {
        add(descriptors,
            makeDescriptor(id, conceptKey,
                           SignalSourceKind::DenseH5Trajectory,
                           "TrajectoryH5",
                           "dihedral_autocorrelation",
                           label,
                           SourceResidency::StartupLoaded,
                           SignalAxis::Residue,
                           SignalAxis::Residue,
                           SignalValueShape::CurveOverLag,
                           none,
                           {},
                           {QStringLiteral("static.curve.lag.animated")},
                           scalarChannels(none),
                           "/trajectory/dihedral_autocorrelation",
                           true,
                           false,
                           SampleStatus::Valid,
                           GapReason::None));
    };
    addDihedralCurve("h5:dihedral_phi_acf", "dihedral.phi_acf", "phi torsional ACF");
    addDihedralCurve("h5:dihedral_psi_acf", "dihedral.psi_acf", "psi torsional ACF");

    // Chi[0..3] composite descriptors (L-2a, 2026-05-29). Option B per
    // user choice — one PerClassBlock scalar + one CurveOverLag curve,
    // each fanned across 4 chi channels. Matches the existing
    // kernel_dynamics composite shape; per-channel dispatch lives in
    // the controller's denseH5Plan branch + panel builders.
    QVector<ChannelDescriptor> chiChannels;
    for (int k = 0; k < 4; ++k) {
        const QString id = QStringLiteral("chi%1").arg(k);
        const QString label = QStringLiteral("Chi %1").arg(k);
        chiChannels.push_back(channel(id.toLatin1().constData(),
                                      label.toLatin1().constData(),
                                      SignalValueShape::Scalar, none));
    }
    add(descriptors,
        makeDescriptor("h5:dihedral_chi_corr_time",
                       "dihedral.chi_corr_time",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "dihedral_autocorrelation",
                       "chi[0..3] torsional decorrelation times (ps)",
                       SourceResidency::StartupLoaded,
                       SignalAxis::Residue,
                       SignalAxis::Residue,
                       SignalValueShape::PerClassBlock,
                       unit(UnitDimension::Time, "ps", "ps"),
                       {},
                       {QStringLiteral("static.bar.sequence"),
                        QStringLiteral("static.table")},
                       chiChannels,
                       "/trajectory/dihedral_autocorrelation",
                       true,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));
    add(descriptors,
        makeDescriptor("h5:dihedral_chi_acf",
                       "dihedral.chi_acf",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "dihedral_autocorrelation",
                       "chi[0..3] torsional ACF over lag",
                       SourceResidency::StartupLoaded,
                       SignalAxis::Residue,
                       SignalAxis::Residue,
                       SignalValueShape::CurveOverLag,
                       none,
                       {},
                       {QStringLiteral("static.curve.lag.animated")},
                       chiChannels,
                       "/trajectory/dihedral_autocorrelation",
                       true,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));

    // ── Bond-vector axis (Reorientational dynamics / Lipari-Szabo) ──
    // Nine descriptors all keyed to /trajectory/reorientational_dynamics:
    // - Five scalars (s2, tau_e, r1, r2, noe) via static.bar.sequence
    //   (SequenceBarPanel). Auto-compose into one panel when 2+ are
    //   active in the same dashboard panel (L-4 builder).
    // - Two TCF curves (body / lab frame) via static.curve.lag.animated
    //   (LagDecayPanel).
    // - Orientation tensor (Mat3 per vector) via static.tensor — the
    //   3-D ellipsoid glyph (L-3a math + revealTensor API; the trigger
    //   gesture is intentionally deferred per the planning conversation
    //   2026-05-29, so the glyph does not auto-fire on signal addition).
    // - Spectral density J(ω) at 5 KTB Larmor frequencies via
    //   static.fixed_freq (FixedFreqPanel, L-3b).
    auto addReorientScalar = [&](const char* id, const char* conceptKey,
                                  const char* label, const UnitSpec& units) {
        add(descriptors,
            makeDescriptor(id, conceptKey,
                           SignalSourceKind::DenseH5Trajectory,
                           "TrajectoryH5",
                           "reorientational_dynamics",
                           label,
                           SourceResidency::StartupLoaded,
                           SignalAxis::BondVector,
                           SignalAxis::BondVector,
                           SignalValueShape::Scalar,
                           units,
                           {},
                           {QStringLiteral("static.bar.sequence"),
                            QStringLiteral("static.table")},
                           scalarChannels(units),
                           "/trajectory/reorientational_dynamics",
                           true,
                           false,
                           SampleStatus::Valid,
                           GapReason::None));
    };
    addReorientScalar("h5:reorient_s2",  "reorient.s2",  "Reorientational S²", none);
    addReorientScalar("h5:reorient_tau_e","reorient.tau_e","Reorientational τ_e (ps)",
                      unit(UnitDimension::Time, "ps", "ps"));
    addReorientScalar("h5:reorient_r1",   "reorient.r1",  "15N R1 (NH only, s⁻¹)",
                      unit(UnitDimension::Frequency, "1/s", "1/s"));
    addReorientScalar("h5:reorient_r2",   "reorient.r2",  "15N R2 (NH only, s⁻¹)",
                      unit(UnitDimension::Frequency, "1/s", "1/s"));
    addReorientScalar("h5:reorient_noe",  "reorient.noe", "15N {1H} NOE (NH only)", none);

    auto addReorientCurve = [&](const char* id, const char* conceptKey, const char* label) {
        add(descriptors,
            makeDescriptor(id, conceptKey,
                           SignalSourceKind::DenseH5Trajectory,
                           "TrajectoryH5",
                           "reorientational_dynamics",
                           label,
                           SourceResidency::StartupLoaded,
                           SignalAxis::BondVector,
                           SignalAxis::BondVector,
                           SignalValueShape::CurveOverLag,
                           none,
                           {},
                           {QStringLiteral("static.curve.lag.animated")},
                           scalarChannels(none),
                           "/trajectory/reorientational_dynamics",
                           true,
                           false,
                           SampleStatus::Valid,
                           GapReason::None));
    };
    addReorientCurve("h5:reorient_acf_internal", "reorient.acf_internal",
                     "Reorientational TCF (body frame)");
    addReorientCurve("h5:reorient_acf_lab",      "reorient.acf_lab",
                     "Reorientational TCF (lab frame)");

    // L-3a (2026-05-29): per-vector Mat3 orientation tensor. The Mat3
    // payload feeds an ellipsoid glyph in the 3-D scene via
    // SceneRevealOverlay::revealTensor. static.table is kept as a
    // fallback so the user can inspect the raw 9 components even when
    // the glyph is disabled.
    add(descriptors,
        makeDescriptor("h5:reorient_orientation_tensor",
                       "reorient.orientation_tensor",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "reorientational_dynamics",
                       "Bond-frame orientation tensor ⟨u⊗u⟩ (Mat3 per vector)",
                       SourceResidency::StartupLoaded,
                       SignalAxis::BondVector,
                       SignalAxis::BondVector,
                       SignalValueShape::Mat3PerRow,
                       none,
                       {},
                       {QStringLiteral("static.tensor"),
                        QStringLiteral("static.table")},
                       scalarChannels(none),
                       "/trajectory/reorientational_dynamics",
                       true,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));

    // L-3b (2026-05-29): per-vector J(ω) sampled at the 5 KTB Larmor
    // frequencies (FixedFreqBlock shape). Rendered via FixedFreqPanel
    // (static.fixed_freq mode). NH-only — Cα-Hα and C=O rows carry
    // NaN per the producer; the panel filters those out at paint time.
    add(descriptors,
        makeDescriptor("h5:reorient_spectral_density",
                       "reorient.spectral_density",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "reorientational_dynamics",
                       "Spectral density J(ω) at 5 KTB Larmor combinations (NH only)",
                       SourceResidency::StartupLoaded,
                       SignalAxis::BondVector,
                       SignalAxis::BondVector,
                       SignalValueShape::FixedFreqBlock,
                       unit(UnitDimension::Time, "s", "s"),
                       {},
                       {QStringLiteral("static.fixed_freq"),
                        QStringLiteral("static.table")},
                       scalarChannels(unit(UnitDimension::Time, "s", "s")),
                       "/trajectory/reorientational_dynamics",
                       true,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));

    // ── Bond-vector axis (Lipari-Szabo / iRED) ──────────────────────
    // IRed S² is per-N-H-vector. Displayed as a per-residue sequence
    // bar (NMR-relaxation idiom) via static.bar.sequence on the
    // BondVector axis. The producer's vector identity (residue+kind)
    // is rebuilt into the QtIRedOrderParameters::identity table; the
    // sampler resolves a BondVectorAnchor(residue, kind=NH) to a row.
    add(descriptors,
        makeDescriptor("h5:ired_s2",
                       "ired.s2",
                       SignalSourceKind::DenseH5Trajectory,
                       "TrajectoryH5",
                       "lipari_szabo",
                       "iRED order parameter (N-H)",
                       SourceResidency::StartupLoaded,
                       SignalAxis::BondVector,
                       SignalAxis::BondVector,
                       SignalValueShape::Scalar,
                       none,
                       {},                                    // not a temporal strip
                       {QStringLiteral("static.bar.sequence"),
                        QStringLiteral("static.table")},
                       scalarChannels(none),
                       "/trajectory/ired_order_parameters",
                       true,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));
    add(descriptors, makeDescriptor("h5:dssp8_transition", "dssp_ss8.transition", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "dssp", "DSSP8 transition matrix", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::EventRecord, tag, eventStripModes(), eventStaticModes(), {}, "/trajectory/dssp8_transition", true));
    add(descriptors, makeDescriptor("h5:dihedral_bin_transition", "residue_dihedral.transition", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "planar_geometry", "Dihedral-bin transition matrix", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::EventRecord, tag, eventStripModes(), eventStaticModes(), {}, "/trajectory/dihedral_bin_transition", true));
    add(descriptors, makeDescriptor("h5:selections", "selections", SignalSourceKind::SelectionEvents, "TrajectoryH5", "selections", "Trajectory selections", SourceResidency::StartupLoaded, SignalAxis::Event, SignalAxis::Event, SignalValueShape::EventRecord, tag, eventStripModes(), eventStaticModes(), {}, "/trajectory/selections", true));
}

void addFrameNpy(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec none = unit(UnitDimension::Dimensionless, "", "");
    const UnitSpec length = unit(UnitDimension::Length, "A", "A");
    const UnitSpec angle = unit(UnitDimension::Angle, "degrees", "deg");
    const UnitSpec shielding = unit(UnitDimension::MagneticShielding, "ppm", "ppm");
    const UnitSpec ringShielding = unit(UnitDimension::MagneticShielding, "ppm_T_per_nA", "ppm_T_per_nA");
    const UnitSpec bfield = unit(UnitDimension::MagneticShielding, "T", "T");
    const UnitSpec efield = unit(UnitDimension::ElectricField, "V/A", "V/A");
    const UnitSpec efg = unit(UnitDimension::ElectricFieldGradient, "V/A^2", "V/A^2");
    const UnitSpec charge = unit(UnitDimension::Charge, "e", "e");
    const UnitSpec energy = unit(UnitDimension::Energy, "kJ/mol", "kJ/mol");
    const UnitSpec tag = unit(UnitDimension::Tag, "tag", "tag", 1.0, 0.0, false);

    auto npy = [&](const char* field,
                   const char* conceptKey,
                   const char* family,
                   const char* label,
                   SignalAxis axis,
                   SignalValueShape shape,
                   const UnitSpec& units,
                   const QStringList& temporalModes,
                   const QStringList& staticModes,
                   const QVector<ChannelDescriptor>& channels) {
        SignalDescriptor descriptor = makeDescriptor("npy:field",
                                                     conceptKey,
                                                     SignalSourceKind::FrameNpySnapshot,
                                                     "SDK_NPY",
                                                     family,
                                                     label,
                                                     SourceResidency::FrameLoaded,
                                                     axis,
                                                     axis,
                                                     shape,
                                                     units,
                                                     temporalModes,
                                                     staticModes,
                                                     channels,
                                                     field,
                                                     false,
                                                     true);
        descriptor.id = QStringLiteral("npy:%1").arg(QString::fromLatin1(field));
        add(descriptors, descriptor);
    };

    npy("pos", "positions", "identity", "Snapshot positions", SignalAxis::Atom, SignalValueShape::Vector3, length, vectorStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.table")}, vectorChannels(length));
    npy("element", "element", "identity", "Element identity", SignalAxis::Atom, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.topology"), QStringLiteral("static.category"), QStringLiteral("static.table")}, {channel("element", "Element", SignalValueShape::Category, tag)});
    npy("residue_index", "residue_index", "identity", "Residue index identity", SignalAxis::Atom, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.topology"), QStringLiteral("static.category"), QStringLiteral("static.table")}, {channel("residue_index", "Residue index", SignalValueShape::Category, tag)});
    npy("residue_type", "residue_type", "identity", "Residue type identity", SignalAxis::Atom, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.topology"), QStringLiteral("static.category"), QStringLiteral("static.table")}, {channel("residue_type", "Residue type", SignalValueShape::Category, tag)});
    npy("atoms_category_info", "atoms_category_info", "identity", "Atom category info", SignalAxis::Atom, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.topology"), QStringLiteral("static.category"), QStringLiteral("static.table")}, {});
    npy("ring_contributions", "ring_contributions", "identity", "Ring contributions", SignalAxis::RingContributionPair, SignalValueShape::TensorComponents, ringShielding, tensorStripModes(), {QStringLiteral("static.tensor"), QStringLiteral("static.per-class"), QStringLiteral("static.table")}, {});
    npy("ring_geometry", "ring_geometry", "identity", "Ring geometry", SignalAxis::AromaticRing, SignalValueShape::TensorComponents, length, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {});

    const struct TensorField {
        const char* field;
        const char* conceptKey;
        const char* family;
        const char* label;
        UnitSpec units;
    } tensorFields[] = {
        {"bs_shielding", "bs_shielding", "biot_savart", "Biot-Savart shielding", ringShielding},
        {"hm_shielding", "hm_shielding", "haigh_mallion", "Haigh-Mallion shielding", ringShielding},
        {"pq_shielding", "pq_shielding", "pi_quadrupole", "Pi quadrupole shielding", ringShielding},
        {"disp_shielding", "disp_shielding", "dispersion", "Dispersion shielding", ringShielding},
        {"ringchi_shielding", "ringchi_shielding", "ring_susceptibility", "Ring susceptibility shielding", shielding},
        {"mc_shielding", "mc_shielding", "mcconnell", "McConnell shielding", shielding},
        {"coulomb_shielding", "coulomb_shielding", "coulomb", "Coulomb shielding", shielding},
        {"hbond_shielding", "hbond_shielding", "hbond", "H-bond shielding", shielding},
        {"mopac_coulomb_shielding", "mopac_coulomb_shielding", "mopac_coulomb", "MOPAC Coulomb shielding", shielding},
        {"mopac_mc_shielding", "mopac_mc_shielding", "mopac_mcconnell", "MOPAC McConnell shielding", shielding},
        {"orca_total", "orca_total", "orca", "ORCA total shielding snapshot", shielding},
        {"orca_diamagnetic", "orca_diamagnetic", "orca", "ORCA diamagnetic shielding snapshot", shielding},
        {"orca_paramagnetic", "orca_paramagnetic", "orca", "ORCA paramagnetic shielding snapshot", shielding},
        {"tripeptide_bb_shielding", "tripeptide_bb_shielding", "tripeptide", "Tripeptide backbone shielding", shielding},
        {"tripeptide_neighbor_shielding", "tripeptide_neighbor_shielding", "tripeptide", "Tripeptide neighbor shielding", shielding},
        {"larsen_hbond_shielding", "larsen_hbond_shielding", "larsen_hbond", "Larsen H-bond shielding", shielding},
        {"larsen_hbond_1pHB_shielding", "larsen_hbond_1pHB_shielding", "larsen_hbond", "Larsen 1pHB shielding", shielding},
        {"larsen_hbond_2pHB_shielding", "larsen_hbond_2pHB_shielding", "larsen_hbond", "Larsen 2pHB shielding", shielding},
        {"larsen_hbond_1pHaB_shielding", "larsen_hbond_1pHaB_shielding", "larsen_hbond", "Larsen 1pHaB shielding", shielding},
        {"larsen_hbond_2pHaB_shielding", "larsen_hbond_2pHaB_shielding", "larsen_hbond", "Larsen 2pHaB shielding", shielding},
        {"larsen_hbond_diagnostic_CB_shielding", "larsen_hbond_diagnostic_CB_shielding", "larsen_hbond", "Larsen diagnostic CB shielding", shielding},
    };

    for (const TensorField& field : tensorFields)
        npy(field.field, field.conceptKey, field.family, field.label, SignalAxis::Atom, SignalValueShape::SphericalTensor, field.units, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(field.units));

    const struct PerClassField {
        const char* field;
        const char* conceptKey;
        const char* family;
        const char* label;
        UnitSpec units;
    } perClassFields[] = {
        {"bs_per_type_T0", "bs_per_type_T0", "biot_savart", "Biot-Savart per-type T0", ringShielding},
        {"bs_per_type_T2", "bs_per_type_T2", "biot_savart", "Biot-Savart per-type T2", ringShielding},
        {"hm_per_type_T0", "hm_per_type_T0", "haigh_mallion", "Haigh-Mallion per-type T0", ringShielding},
        {"hm_per_type_T2", "hm_per_type_T2", "haigh_mallion", "Haigh-Mallion per-type T2", ringShielding},
        {"pq_per_type_T0", "pq_per_type_T0", "pi_quadrupole", "Pi quadrupole per-type T0", ringShielding},
        {"pq_per_type_T2", "pq_per_type_T2", "pi_quadrupole", "Pi quadrupole per-type T2", ringShielding},
        {"disp_per_type_T0", "disp_per_type_T0", "dispersion", "Dispersion per-type T0", ringShielding},
        {"disp_per_type_T2", "disp_per_type_T2", "dispersion", "Dispersion per-type T2", ringShielding},
        {"mc_category_T2", "mc_category_T2", "mcconnell", "McConnell category T2", shielding},
        {"mopac_mc_category_T2", "mopac_mc_category_T2", "mopac_mcconnell", "MOPAC McConnell category T2", shielding},
    };
    for (const PerClassField& field : perClassFields)
        npy(field.field, field.conceptKey, field.family, field.label, SignalAxis::Atom, SignalValueShape::PerClassBlock, field.units, perClassStripModes(), perClassStaticModes(), {});

    npy("bs_total_B", "bs_total_B", "biot_savart", "Biot-Savart total B field", SignalAxis::Atom, SignalValueShape::Vector3, bfield, vectorStripModes(), vectorStaticModes(), vectorChannels(bfield));
    npy("bs_ring_counts", "bs_ring_counts", "biot_savart", "Biot-Savart ring counts", SignalAxis::Atom, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.per-class")}, perClassStaticModes(), countChannels());
    npy("coulomb_E", "coulomb_E", "coulomb", "Coulomb electric field", SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield));
    npy("coulomb_efg_backbone", "coulomb_efg_backbone", "coulomb", "Coulomb backbone EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("coulomb_efg_aromatic", "coulomb_efg_aromatic", "coulomb", "Coulomb aromatic EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("coulomb_scalars", "coulomb_scalars", "coulomb", "Coulomb scalar diagnostics", SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("hbond_scalars", "hbond_scalars", "hbond", "H-bond scalar diagnostics", SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("dssp_backbone", "dssp_backbone", "dssp", "DSSP backbone geometry", SignalAxis::Residue, SignalValueShape::PerClassBlock, angle, perClassStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {});
    npy("dssp_ss8", "dssp_ss8", "dssp", "DSSP SS8", SignalAxis::Residue, SignalValueShape::Category, tag, categoryStripModes(), categoryStaticModes(), {channel("ss8", "SS8", SignalValueShape::Category, tag)});
    npy("dssp_hbond_energy", "dssp_hbond_energy", "dssp", "DSSP H-bond energy", SignalAxis::Residue, SignalValueShape::PerClassBlock, energy, perClassStripModes(), perClassStaticModes(), {});
    npy("dssp_chi", "dssp_chi", "dssp", "DSSP chi angles", SignalAxis::Residue, SignalValueShape::PerClassBlock, angle, perClassStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {});
    npy("atom_sasa", "atom_sasa", "sasa", "Atom SASA", SignalAxis::Atom, SignalValueShape::Scalar, unit(UnitDimension::Length, "A^2", "A^2"), scalarStripModes(), scalarStaticModes(), scalarChannels(unit(UnitDimension::Length, "A^2", "A^2")));
    npy("sasa_normal", "sasa_normal", "sasa", "SASA normal", SignalAxis::Atom, SignalValueShape::Vector3, none, vectorStripModes(), vectorStaticModes(), vectorChannels(none));
    npy("water_efield", "water_efield", "water_field", "Water electric field", SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield));
    npy("water_efield_first", "water_efield_first", "water_field", "First-shell water electric field", SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield));
    npy("water_efg", "water_efg", "water_field", "Water EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("water_efg_first", "water_efg_first", "water_field", "First-shell water EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("water_shell_counts", "water_shell_counts", "water_field", "Water shell counts", SignalAxis::Atom, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.per-class")}, perClassStaticModes(), countChannels());
    npy("hydration_shell", "hydration_shell", "hydration", "Hydration shell", SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("water_polarization", "water_polarization", "water_polarization", "Water polarization", SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("eeq_charges", "eeq_charges", "eeq", "EEQ charges", SignalAxis::Atom, SignalValueShape::Scalar, charge, scalarStripModes(), scalarStaticModes(), scalarChannels(charge));
    npy("eeq_cn", "eeq_cn", "eeq", "EEQ coordination number", SignalAxis::Atom, SignalValueShape::Scalar, none, scalarStripModes(), scalarStaticModes(), scalarChannels(none));
    npy("gromacs_energy", "gromacs_energy", "gromacs", "Gromacs energy/runtime", SignalAxis::System, SignalValueShape::PerClassBlock, none, {QStringLiteral("strip.system"), QStringLiteral("strip.per-class")}, {QStringLiteral("static.system"), QStringLiteral("static.table")}, {});
    npy("bonded_energy", "bonded_energy", "bonded", "Bonded energy", SignalAxis::Atom, SignalValueShape::PerClassBlock, energy, perClassStripModes(), perClassStaticModes(), {});
    npy("mopac_charges", "mopac_charges", "mopac_core", "MOPAC charges", SignalAxis::Atom, SignalValueShape::Scalar, charge, scalarStripModes(), scalarStaticModes(), scalarChannels(charge));
    npy("mopac_scalars", "mopac_scalars", "mopac_core", "MOPAC scalar diagnostics", SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("mopac_bond_orders", "mopac_bond_orders", "mopac_core", "MOPAC bond orders", SignalAxis::Bond, SignalValueShape::Scalar, none, scalarStripModes(), {QStringLiteral("static.scalar"), QStringLiteral("static.topology"), QStringLiteral("static.table")}, scalarChannels(none));
    npy("mopac_global", "mopac_global", "mopac_core", "MOPAC global values", SignalAxis::System, SignalValueShape::PerClassBlock, none, {QStringLiteral("strip.system"), QStringLiteral("strip.per-class")}, {QStringLiteral("static.system"), QStringLiteral("static.table")}, {});
    npy("mopac_coulomb_E", "mopac_coulomb_E", "mopac_coulomb", "MOPAC Coulomb electric field", SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield));
    npy("mopac_coulomb_efg_backbone", "mopac_coulomb_efg_backbone", "mopac_coulomb", "MOPAC Coulomb backbone EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("mopac_coulomb_efg_aromatic", "mopac_coulomb_efg_aromatic", "mopac_coulomb", "MOPAC Coulomb aromatic EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("mopac_coulomb_scalars", "mopac_coulomb_scalars", "mopac_coulomb", "MOPAC Coulomb scalar diagnostics", SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("mopac_mc_scalars", "mopac_mc_scalars", "mopac_mcconnell", "MOPAC McConnell scalar diagnostics", SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("apbs_E", "apbs_E", "apbs", "APBS electric field", SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield));
    npy("apbs_efg", "apbs_efg", "apbs", "APBS EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("delta_shielding", "delta_shielding", "mutation_delta", "Mutation delta shielding", SignalAxis::MutationMatchPair, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding));
    npy("delta_scalars", "delta_scalars", "mutation_delta", "Mutation delta scalars", SignalAxis::MutationMatchPair, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("delta_apbs", "delta_apbs", "mutation_delta", "Mutation delta APBS", SignalAxis::MutationMatchPair, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("delta_ring_proximity", "delta_ring_proximity", "mutation_delta", "Mutation delta ring proximity", SignalAxis::MutationMatchPair, SignalValueShape::PerClassBlock, length, perClassStripModes(), perClassStaticModes(), {});
    npy("wt_shielding_diamagnetic", "wt_shielding_diamagnetic", "mutation_delta", "WT diamagnetic shielding", SignalAxis::MutationMatchPair, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding));
    npy("wt_shielding_paramagnetic", "wt_shielding_paramagnetic", "mutation_delta", "WT paramagnetic shielding", SignalAxis::MutationMatchPair, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding));
    npy("mut_shielding_diamagnetic", "mut_shielding_diamagnetic", "mutation_delta", "Mutant diamagnetic shielding", SignalAxis::MutationMatchPair, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding));
    npy("mut_shielding_paramagnetic", "mut_shielding_paramagnetic", "mutation_delta", "Mutant paramagnetic shielding", SignalAxis::MutationMatchPair, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding));
    npy("delta_shielding_diamagnetic", "delta_shielding_diamagnetic", "mutation_delta", "Delta diamagnetic shielding", SignalAxis::MutationMatchPair, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding));
    npy("delta_shielding_paramagnetic", "delta_shielding_paramagnetic", "mutation_delta", "Delta paramagnetic shielding", SignalAxis::MutationMatchPair, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding));
    npy("aimnet2_charges", "aimnet2_charges", "aimnet2", "AIMNet2 charges", SignalAxis::Atom, SignalValueShape::Scalar, charge, scalarStripModes(), scalarStaticModes(), scalarChannels(charge));
    npy("aimnet2_aim", "aimnet2_aim", "aimnet2", "AIMNet2 AIM descriptors", SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {});
    npy("aimnet2_efg", "aimnet2_efg", "aimnet2", "AIMNet2 EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("aimnet2_efg_aromatic", "aimnet2_efg_aromatic", "aimnet2", "AIMNet2 aromatic EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("aimnet2_efg_backbone", "aimnet2_efg_backbone", "aimnet2", "AIMNet2 backbone EFG", SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg));
    npy("aimnet2_charge_response_gradient", "aimnet2_charge_response_gradient", "aimnet2", "AIMNet2 charge-response gradient", SignalAxis::Atom, SignalValueShape::Vector3, charge, vectorStripModes(), vectorStaticModes(), vectorChannels(charge));
    npy("aimnet2_charge_response_gradient_scalar", "aimnet2_charge_response_gradient_scalar", "aimnet2", "AIMNet2 CRG scalar", SignalAxis::Atom, SignalValueShape::Scalar, charge, scalarStripModes(), scalarStaticModes(), scalarChannels(charge));
    npy("pyramidalization", "pyramidalization", "planar_geometry", "Pyramidalization", SignalAxis::Atom, SignalValueShape::Scalar, angle, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, angleChannels());
    npy("omega_actual", "omega_actual", "planar_geometry", "Omega actual", SignalAxis::Residue, SignalValueShape::Scalar, angle, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, angleChannels());
    npy("omega_deviation", "omega_deviation", "planar_geometry", "Omega deviation", SignalAxis::Residue, SignalValueShape::Scalar, angle, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, angleChannels());
    npy("aromatic_chi2", "aromatic_chi2", "planar_geometry", "Aromatic chi2", SignalAxis::AromaticRing, SignalValueShape::Scalar, angle, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, angleChannels());
    npy("pucker_Q", "pucker_Q", "planar_geometry", "Ring pucker Q", SignalAxis::SaturatedRing, SignalValueShape::Scalar, length, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, scalarChannels(length));
    npy("pucker_theta", "pucker_theta", "planar_geometry", "Ring pucker theta", SignalAxis::SaturatedRing, SignalValueShape::Scalar, angle, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, angleChannels());
    npy("omega_is_xpro", "omega_is_xpro", "planar_geometry", "Omega X-Proline flag", SignalAxis::Residue, SignalValueShape::Category, tag, categoryStripModes(), categoryStaticModes(), {channel("is_xpro", "Is X-Proline", SignalValueShape::Category, tag)});
    npy("tripeptide_bb_residual_vec", "tripeptide_bb_residual_vec", "tripeptide", "Tripeptide backbone residual vector", SignalAxis::Atom, SignalValueShape::Vector3, shielding, vectorStripModes(), vectorStaticModes(), vectorChannels(shielding));
    npy("tripeptide_bb_match_distance", "tripeptide_bb_match_distance", "tripeptide", "Tripeptide backbone match distance", SignalAxis::Atom, SignalValueShape::Scalar, length, scalarStripModes(), scalarStaticModes(), scalarChannels(length));
    npy("tripeptide_bb_method_tag", "tripeptide_bb_method_tag", "tripeptide", "Tripeptide backbone method tag", SignalAxis::Atom, SignalValueShape::Category, tag, categoryStripModes(), categoryStaticModes(), {channel("method", "Method", SignalValueShape::Category, tag)});
    npy("tripeptide_neighbor_residual_vec_prev", "tripeptide_neighbor_residual_vec_prev", "tripeptide", "Tripeptide previous-neighbor residual vector", SignalAxis::Atom, SignalValueShape::Vector3, shielding, vectorStripModes(), vectorStaticModes(), vectorChannels(shielding));
    npy("tripeptide_neighbor_residual_vec_next", "tripeptide_neighbor_residual_vec_next", "tripeptide", "Tripeptide next-neighbor residual vector", SignalAxis::Atom, SignalValueShape::Vector3, shielding, vectorStripModes(), vectorStaticModes(), vectorChannels(shielding));
    npy("larsen_hbond_water_term", "larsen_hbond_water_term", "larsen_hbond", "Larsen H-bond water term", SignalAxis::Atom, SignalValueShape::Scalar, shielding, scalarStripModes(), scalarStaticModes(), scalarChannels(shielding));
    npy("larsen_hbond_count", "larsen_hbond_count", "larsen_hbond", "Larsen H-bond count", SignalAxis::Atom, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.scalar")}, scalarStaticModes(), countChannels());
    npy("residues", "topology.residues", "topology", "Residue sidecar records", SignalAxis::Residue, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.topology"), QStringLiteral("static.table")}, {});
    npy("bonds", "topology.bonds", "topology", "Bond sidecar records", SignalAxis::Bond, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.topology"), QStringLiteral("static.table")}, {});
    npy("rings", "topology.rings", "topology", "Ring sidecar records", SignalAxis::Ring, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.topology"), QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {});
    npy("ring_membership", "topology.ring_membership", "topology", "Ring membership records", SignalAxis::RingMembership, SignalValueShape::Category, tag, categoryStripModes(), {QStringLiteral("static.topology"), QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {});
}

void addOrcaDft(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec shielding = unit(UnitDimension::MagneticShielding, "ppm", "ppm");
    add(descriptors,
        makeDescriptor("orca_dft:total",
                       "orca_total",
                       SignalSourceKind::OrcaDftFrame,
                       "ORCA",
                       "orca",
                       "ORCA total shielding",
                       SourceResidency::FrameLoaded,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::SphericalTensor,
                       shielding,
                       {QStringLiteral("strip.tensor.T0"), QStringLiteral("strip.spectrum")},
                       {QStringLiteral("static.tensor"), QStringLiteral("static.scalar"), QStringLiteral("static.table"), QStringLiteral("static.atomColor")},
                       sphericalTensorChannels(shielding),
                       "orca_total",
                       false,
                       true,
                       SampleStatus::Valid,
                       GapReason::None));
    add(descriptors, makeDescriptor("orca_dft:diamagnetic", "orca_diamagnetic", SignalSourceKind::OrcaDftFrame, "ORCA", "orca", "ORCA diamagnetic shielding", SourceResidency::FrameLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding), "orca_diamagnetic", false, true));
    add(descriptors, makeDescriptor("orca_dft:paramagnetic", "orca_paramagnetic", SignalSourceKind::OrcaDftFrame, "ORCA", "orca", "ORCA paramagnetic shielding", SourceResidency::FrameLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding), "orca_paramagnetic", false, true));
}

void addTopology(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec tag = unit(UnitDimension::Tag, "tag", "tag", 1.0, 0.0, false);
    const UnitSpec length = unit(UnitDimension::Length, "A", "A");
    add(descriptors, makeDescriptor("topology:atoms", "topology.atoms", SignalSourceKind::Topology, "TopologySidecar", "topology", "Atoms", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Category, tag, {}, {QStringLiteral("static.topology"), QStringLiteral("static.category"), QStringLiteral("static.table")}, {}, "atoms"));
    add(descriptors, makeDescriptor("topology:residues", "topology.residues", SignalSourceKind::Topology, "TopologySidecar", "topology", "Residues", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::Category, tag, {}, {QStringLiteral("static.topology"), QStringLiteral("static.category"), QStringLiteral("static.table")}, {}, "residues"));
    add(descriptors, makeDescriptor("topology:bonds", "topology.bonds", SignalSourceKind::Topology, "TopologySidecar", "topology", "Bonds", SourceResidency::StartupLoaded, SignalAxis::Bond, SignalAxis::Bond, SignalValueShape::Category, tag, {}, {QStringLiteral("static.topology"), QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {}, "bonds"));
    add(descriptors, makeDescriptor("topology:rings", "topology.rings", SignalSourceKind::Topology, "TopologySidecar", "topology", "Rings", SourceResidency::StartupLoaded, SignalAxis::Ring, SignalAxis::Ring, SignalValueShape::Category, tag, {}, {QStringLiteral("static.topology"), QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {}, "rings"));
    add(descriptors, makeDescriptor("topology:ring_membership", "topology.ring_membership", SignalSourceKind::Topology, "TopologySidecar", "topology", "Ring membership", SourceResidency::StartupLoaded, SignalAxis::RingMembership, SignalAxis::RingMembership, SignalValueShape::Category, tag, {}, {QStringLiteral("static.topology"), QStringLiteral("static.geometry"), QStringLiteral("static.table")}, {}, "ring_membership"));
    add(descriptors, makeDescriptor("topology:bond_length", "geometry.bond_length", SignalSourceKind::Topology, "DerivedTopology", "topology", "Bond length from topology positions", SourceResidency::Derived, SignalAxis::Bond, SignalAxis::Bond, SignalValueShape::Scalar, length, {QStringLiteral("strip.scalar"), QStringLiteral("strip.rollup")}, {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, scalarChannels(length), "bond_length"));
}

void addDerivedGeometry(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec length = unit(UnitDimension::Length, "A", "A");
    const UnitSpec degrees = unit(UnitDimension::Angle, "degrees", "deg");
    add(descriptors,
        makeDescriptor("geometry:distance",
                       "geometry.distance",
                       SignalSourceKind::DerivedGeometry,
                       "Derived",
                       "geometry",
                       "Distance",
                       SourceResidency::Derived,
                       SignalAxis::AtomTuple,
                       SignalAxis::AtomTuple,
                       SignalValueShape::Scalar,
                       length,
                       scalarStripModes(),
                       {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")},
                       scalarChannels(length),
                       "distance",
                       false,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));
    add(descriptors, makeDescriptor("geometry:angle", "geometry.angle", SignalSourceKind::DerivedGeometry, "Derived", "geometry", "Angle", SourceResidency::Derived, SignalAxis::AtomTuple, SignalAxis::AtomTuple, SignalValueShape::Scalar, degrees, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, angleChannels(), "angle", false, false, SampleStatus::Valid, GapReason::None));
    add(descriptors, makeDescriptor("geometry:dihedral", "geometry.dihedral", SignalSourceKind::DerivedGeometry, "Derived", "geometry", "Dihedral", SourceResidency::Derived, SignalAxis::AtomTuple, SignalAxis::AtomTuple, SignalValueShape::Scalar, degrees, scalarStripModes(), {QStringLiteral("static.geometry"), QStringLiteral("static.scalar"), QStringLiteral("static.table")}, angleChannels(), "dihedral", false, false, SampleStatus::Valid, GapReason::None));
    add(descriptors, makeDescriptor("geometry:atom_displacement", "geometry.atom_displacement", SignalSourceKind::DerivedGeometry, "Derived", "geometry", "Atom displacement", SourceResidency::Derived, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, length, vectorStripModes(), vectorStaticModes(), vectorChannels(length), "atom_displacement"));
}

void addSelectionEvents(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec tag = unit(UnitDimension::Tag, "tag", "tag", 1.0, 0.0, false);
    add(descriptors, makeDescriptor("events:selection_timeline", "selections", SignalSourceKind::SelectionEvents, "SelectionBag", "selections", "Selection timeline", SourceResidency::StartupLoaded, SignalAxis::Event, SignalAxis::Event, SignalValueShape::EventRecord, tag, eventStripModes(), eventStaticModes(), {}, "selection_timeline"));
    add(descriptors, makeDescriptor("events:selection_counts", "selection_counts", SignalSourceKind::SelectionEvents, "SelectionBag", "selections", "Selection counts", SourceResidency::StartupLoaded, SignalAxis::System, SignalAxis::System, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.system")}, {QStringLiteral("static.system"), QStringLiteral("static.table")}, countChannels(), "selection_counts"));
}

bool descriptorMatchesText(const SignalDescriptor& descriptor, const QString& text) {
    if (text.trimmed().isEmpty())
        return true;
    const QString haystack = QStringList{
                                 descriptor.id,
                                 descriptor.conceptKey,
                                 descriptor.importSet,
                                 descriptor.family,
                                 descriptor.label,
                                 descriptor.description,
                                 descriptor.storagePath,
                                 descriptor.tags.join(QStringLiteral(" ")),
                             }.join(QStringLiteral(" "));
    return haystack.contains(text, Qt::CaseInsensitive);
}

bool descriptorMatchesAxis(const SignalDescriptor& descriptor, SignalAxis axis) {
    return descriptor.nativeAxis == axis || descriptor.requiredAnchor == axis;
}

}  // namespace

TrajectorySignalCatalog::TrajectorySignalCatalog(QObject* parent)
    : QObject(parent)
    , descriptors_(BuildDescriptorCatalog()) {}

const SignalDescriptor* TrajectorySignalCatalog::findDescriptor(const QString& descriptorId) const {
    for (const SignalDescriptor& descriptor : descriptors_) {
        if (descriptor.id == descriptorId)
            return &descriptor;
    }
    return nullptr;
}

std::optional<SignalDescriptor> TrajectorySignalCatalog::descriptorById(const QString& descriptorId) const {
    const SignalDescriptor* descriptor = findDescriptor(descriptorId);
    if (!descriptor)
        return std::nullopt;
    return *descriptor;
}

QVector<SignalDescriptor> TrajectorySignalCatalog::filterDescriptors(const SignalDescriptorFilter& filter) const {
    QVector<SignalDescriptor> result;
    for (const SignalDescriptor& descriptor : descriptors_) {
        if (filter.sourceKind && descriptor.sourceKind != *filter.sourceKind)
            continue;
        if (filter.axis && !descriptorMatchesAxis(descriptor, *filter.axis))
            continue;
        if (filter.valueShape && descriptor.valueShape != *filter.valueShape)
            continue;
        if (!filter.displayModeId.isEmpty() && !SupportsDisplayMode(descriptor, filter.displayModeId))
            continue;
        if (!filter.includePending && descriptor.samplingStatus != SampleStatus::Valid)
            continue;
        if (!filter.includeTemporal && descriptor.temporal && !descriptor.staticDisplay)
            continue;
        if (!filter.includeStatic && descriptor.staticDisplay && !descriptor.temporal)
            continue;
        if (!descriptorMatchesText(descriptor, filter.text))
            continue;
        result.push_back(descriptor);
    }
    return result;
}

QVector<SignalDescriptor> TrajectorySignalCatalog::descriptorsForSource(SignalSourceKind sourceKind) const {
    SignalDescriptorFilter filter;
    filter.sourceKind = sourceKind;
    return filterDescriptors(filter);
}

QVector<SignalDescriptor> TrajectorySignalCatalog::descriptorsForAxis(SignalAxis axis) const {
    SignalDescriptorFilter filter;
    filter.axis = axis;
    return filterDescriptors(filter);
}

QVector<SignalDescriptor> TrajectorySignalCatalog::descriptorsForValueShape(SignalValueShape valueShape) const {
    SignalDescriptorFilter filter;
    filter.valueShape = valueShape;
    return filterDescriptors(filter);
}

QVector<SignalDescriptor> TrajectorySignalCatalog::descriptorsForDisplayMode(const QString& displayModeId) const {
    SignalDescriptorFilter filter;
    filter.displayModeId = displayModeId;
    return filterDescriptors(filter);
}

QVector<SignalDescriptor> TrajectorySignalCatalog::search(const QString& text) const {
    SignalDescriptorFilter filter;
    filter.text = text;
    return filterDescriptors(filter);
}

bool TrajectorySignalCatalog::supportsDisplayMode(const QString& descriptorId, const QString& displayModeId) const {
    const SignalDescriptor* descriptor = findDescriptor(descriptorId);
    return descriptor && SupportsDisplayMode(*descriptor, displayModeId);
}

bool TrajectorySignalCatalog::canBind(const DisplaySignalBinding& binding) const {
    const SignalDescriptor* descriptor = findDescriptor(binding.descriptorId);
    if (!descriptor)
        return false;
    if (descriptor->sourceKind != binding.sourceKind)
        return false;
    if (!binding.conceptKey.isEmpty() && descriptor->conceptKey != binding.conceptKey)
        return false;
    if (!SupportsDisplayMode(*descriptor, binding.displayModeId))
        return false;
    if (!binding.followsFocus && descriptor->requiredAnchor != SignalAxis::None
        && !AnchorMatchesAxis(binding.anchor, descriptor->requiredAnchor)) {
        return false;
    }
    return true;
}

bool TrajectorySignalCatalog::canSample(const DisplaySignalBinding& binding) const {
    const SignalDescriptor* descriptor = findDescriptor(binding.descriptorId);
    return descriptor && descriptor->temporal && binding.displayModeId.startsWith(QStringLiteral("strip."))
        && canBind(binding) && descriptor->samplingStatus == SampleStatus::Valid;
}

QVector<SignalDescriptor> TrajectorySignalCatalog::BuildDescriptorCatalog() {
    QVector<SignalDescriptor> descriptors;
    descriptors.reserve(160);
    addDenseH5(descriptors);
    addFrameNpy(descriptors);
    addOrcaDft(descriptors);
    addTopology(descriptors);
    addDerivedGeometry(descriptors);
    addSelectionEvents(descriptors);
    return descriptors;
}

}  // namespace h5reader::model
