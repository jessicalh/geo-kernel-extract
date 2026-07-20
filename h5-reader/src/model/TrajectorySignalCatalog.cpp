#include "TrajectorySignalCatalog.h"

#include "../io/FrameFieldPolicy.h"

#include <QSet>
#include <QStringList>

#include <array>
#include <optional>
#include <string_view>
#include <utility>

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

QString fromUtf8(std::string_view text) {
    return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size()));
}

std::optional<SignalAxis> signalAxisFor(io::NativeAxis axis) {
    switch (axis) {
    case io::NativeAxis::Atom:
        return SignalAxis::Atom;
    case io::NativeAxis::Residue:
        return SignalAxis::Residue;
    case io::NativeAxis::Bond:
        return SignalAxis::Bond;
    case io::NativeAxis::AromaticRing:
        return SignalAxis::AromaticRing;
    case io::NativeAxis::SaturatedRing:
        return SignalAxis::SaturatedRing;
    case io::NativeAxis::RingContributionPair:
        return SignalAxis::RingContributionPair;
    case io::NativeAxis::Ring:
        return SignalAxis::Ring;
    case io::NativeAxis::RingMembership:
        return SignalAxis::RingMembership;
    case io::NativeAxis::Protein:
        return SignalAxis::Protein;
    default:
        return std::nullopt;
    }
}

UnitSpec producerUnits(io::FieldKind kind) {
    const QString source = fromUtf8(io::FieldSpecFor(kind).units);
    UnitSpec result;
    result.sourceSymbol = source;
    result.displaySymbol = source;

    if (source == QStringLiteral("radians")) {
        result.dimension = UnitDimension::Angle;
        result.displaySymbol = QStringLiteral("deg");
        result.scaleToDisplay = 57.29577951308232;
    } else if (source == QStringLiteral("degrees")) {
        result.dimension = UnitDimension::Angle;
        result.displaySymbol = QStringLiteral("deg");
    } else if (source == QStringLiteral("ppm")
               || source == QStringLiteral("ppm_T_per_nA")) {
        result.dimension = UnitDimension::MagneticShielding;
    } else if (source == QStringLiteral("V/A")) {
        result.dimension = UnitDimension::ElectricField;
    } else if (source == QStringLiteral("V/A^2")) {
        result.dimension = UnitDimension::ElectricFieldGradient;
    } else if (source == QStringLiteral("e")) {
        result.dimension = UnitDimension::Charge;
    } else if (source == QStringLiteral("eV")
               || source == QStringLiteral("Hartree")
               || source == QStringLiteral("kcal/mol")
               || source == QStringLiteral("kJ/mol")) {
        result.dimension = UnitDimension::Energy;
    } else if (source == QStringLiteral("count")) {
        result.dimension = UnitDimension::Count;
    } else if (source == QStringLiteral("A") || source == QString::fromUtf8("\xC3\x85")) {
        result.dimension = UnitDimension::Length;
        result.displaySymbol = QStringLiteral("A");
    } else if (source == QStringLiteral("mask")
               || source.startsWith(QStringLiteral("enum:"))
               || source == QStringLiteral("index")) {
        result.dimension = UnitDimension::Tag;
        result.convertible = false;
    } else {
        result.dimension = UnitDimension::Dimensionless;
        result.convertible = false;
    }
    return result;
}

ChannelDescriptor channel(const char* id,
                          const char* label,
                          SignalValueShape shape,
                          const UnitSpec& sourceUnits,
                          const UnitSpec& displayUnits = UnitSpec{},
                          int firstSourceColumn = 0,
                          int sourceColumnCount = 0) {
    ChannelDescriptor descriptor;
    descriptor.id = QString::fromLatin1(id);
    descriptor.label = QString::fromLatin1(label);
    descriptor.valueShape = shape;
    descriptor.sourceUnits = sourceUnits;
    descriptor.defaultDisplayUnits = displayUnits.sourceSymbol.isEmpty() && displayUnits.displaySymbol.isEmpty()
                                         ? sourceUnits
                                         : displayUnits;
    descriptor.firstSourceColumn = firstSourceColumn;
    descriptor.sourceColumnCount = sourceColumnCount;
    return descriptor;
}

QVector<ChannelDescriptor> scalarChannels(const UnitSpec& units, int sourceColumn = 0) {
    return {channel("value", "Value", SignalValueShape::Scalar, units, {}, sourceColumn, 1)};
}

QVector<ChannelDescriptor> countChannels(int sourceColumn = 0) {
    const UnitSpec count = unit(UnitDimension::Count, "count", "count");
    return {channel("count", "Count", SignalValueShape::Count, count, {}, sourceColumn, 1)};
}

QVector<ChannelDescriptor> vectorChannels(const UnitSpec& units, int firstSourceColumn = 0) {
    return {
        channel("x", "X", SignalValueShape::Scalar, units, {}, firstSourceColumn, 3),
        channel("y", "Y", SignalValueShape::Scalar, units, {}, firstSourceColumn, 3),
        channel("z", "Z", SignalValueShape::Scalar, units, {}, firstSourceColumn, 3),
        channel("magnitude", "Magnitude", SignalValueShape::Scalar, units, {}, firstSourceColumn, 3),
    };
}

QVector<ChannelDescriptor> sphericalTensorChannels(const UnitSpec& units,
                                                    int firstSourceColumn = 0) {
    return {
        channel("T0", "T0", SignalValueShape::Scalar, units, {}, firstSourceColumn, 9),
        channel("T1", "T1", SignalValueShape::TensorComponents, units, {}, firstSourceColumn, 9),
        channel("T2", "T2", SignalValueShape::EfgT2, units, {}, firstSourceColumn, 9),
        channel("component", "Component", SignalValueShape::Scalar, units, {}, firstSourceColumn, 9),
    };
}

QVector<ChannelDescriptor> efgChannels(const UnitSpec& units, int firstSourceColumn = 0) {
    return {
        channel("T2", "T2", SignalValueShape::EfgT2, units, {}, firstSourceColumn, 5),
        channel("magnitude", "Magnitude", SignalValueShape::Scalar, units, {}, firstSourceColumn, 5),
        channel("component", "Component", SignalValueShape::Scalar, units, {}, firstSourceColumn, 5),
    };
}

QVector<ChannelDescriptor> angleChannels(int sourceColumn = 0) {
    const UnitSpec degrees = unit(UnitDimension::Angle, "degrees", "deg");
    return {channel("angle", "Angle", SignalValueShape::Scalar, degrees, {}, sourceColumn, 1)};
}

struct NamedChannel {
    const char* id;
    const char* label;
};

template <std::size_t N>
QVector<ChannelDescriptor> scalarBlockChannels(const std::array<NamedChannel, N>& names,
                                               const UnitSpec& units) {
    QVector<ChannelDescriptor> channels;
    channels.reserve(static_cast<qsizetype>(N));
    for (std::size_t i = 0; i < N; ++i) {
        channels.push_back(channel(names[i].id,
                                   names[i].label,
                                   SignalValueShape::Scalar,
                                   units,
                                   {},
                                   static_cast<int>(i),
                                   1));
    }
    return channels;
}

template <std::size_t N>
QVector<ChannelDescriptor> vectorBlockChannels(const std::array<NamedChannel, N>& names,
                                               const UnitSpec& units) {
    QVector<ChannelDescriptor> channels;
    channels.reserve(static_cast<qsizetype>(N));
    for (std::size_t i = 0; i < N; ++i) {
        channels.push_back(channel(names[i].id,
                                   names[i].label,
                                   SignalValueShape::Vector3,
                                   units,
                                   {},
                                   static_cast<int>(i * 3),
                                   3));
    }
    return channels;
}

template <std::size_t N>
QVector<ChannelDescriptor> t2BlockChannels(const std::array<NamedChannel, N>& names,
                                           const UnitSpec& units) {
    QVector<ChannelDescriptor> channels;
    channels.reserve(static_cast<qsizetype>(N));
    for (std::size_t i = 0; i < N; ++i) {
        channels.push_back(channel(names[i].id,
                                   names[i].label,
                                   SignalValueShape::EfgT2,
                                   units,
                                   {},
                                   static_cast<int>(i * 5),
                                   5));
    }
    return channels;
}

QStringList tensorStripModes() {
    return {
        QStringLiteral("strip.tensor.T0"),
        QStringLiteral("strip.tensor.T1"),
        QStringLiteral("strip.tensor.T2"),
        QStringLiteral("strip.tensor.component"),
    };
}

QStringList tensorStaticModes() {
    return {
        QStringLiteral("static.tensor"),
    };
}

QStringList efgStripModes() {
    return {
        QStringLiteral("strip.tensor.T2"),
        QStringLiteral("strip.tensor.component"),
    };
}

// Strip modes for a SPHERICAL-tensor field that is physically traceless +
// symmetric -- the Coulomb / MOPAC-Coulomb electric-field shielding
// contribution, whose isotropic (T0) and antisymmetric (T1) parts are
// identically zero (measured flat, span ~1e-16). Offer the rank-2 (T2) signal +
// component browse only, matching the EfgT2 H5 sibling; the on-disk snapshot
// stays a full spherical tensor, so the value SHAPE is unchanged.
QStringList tracelessTensorStripModes() {
    return {
        QStringLiteral("strip.tensor.T2"),
        QStringLiteral("strip.tensor.component"),
    };
}

QStringList efgStaticModes() {
    return {
        QStringLiteral("static.tensor"),
    };
}

QStringList scalarStripModes() {
    return {
        QStringLiteral("strip.scalar"),
    };
}

QStringList scalarStaticModes() {
    return {};
}

QStringList vectorStripModes() {
    return {
        QStringLiteral("strip.vector.component"),
        QStringLiteral("strip.vector.magnitude"),
    };
}

QStringList vectorStaticModes() {
    return {};
}

QStringList categoryStripModes() {
    return {
        QStringLiteral("strip.category"),
        QStringLiteral("strip.count"),
    };
}

QStringList categoryStaticModes() {
    return {};
}

QStringList perClassStripModes() {
    return {
        QStringLiteral("strip.per-class"),
        QStringLiteral("strip.scalar"),
    };
}

QStringList perClassStaticModes() {
    return {};
}

// rollupStripModes() removed: rollup-moment summaries are whole-trajectory
// statistics that render as a flat line in a temporal strip (see DisplayPolicy
// IsDashboardDisplayable). They are de-stripped pending a static mean +/- std
// readout; rollupStaticModes() stays as the slot that readout mode will fill.
QStringList rollupStaticModes() {
    return {};
}

QStringList eventStripModes() {
    return {
        QStringLiteral("strip.event"),
        QStringLiteral("strip.category"),
    };
}

QStringList eventStaticModes() {
    return {};
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
            QStringLiteral("/trajectory/mopac_coulomb_efg_time_series"),
            QStringLiteral("/trajectory/mopac_mc_shielding_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_1pHB_shielding_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_1pHaB_shielding_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_2pHB_shielding_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_2pHaB_shielding_time_series"),
            QStringLiteral("/trajectory/sasa_time_series"),
            QStringLiteral("/trajectory/aimnet2_charge_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_count_time_series"),
            QStringLiteral("/trajectory/larsen_hbond_water_term_time_series"),
            QStringLiteral("/trajectory/bonded_energy_time_series"),
            QStringLiteral("/trajectory/water_field_time_series"),
            QStringLiteral("/trajectory/hydration_shell_time_series"),
            QStringLiteral("/trajectory/hydration_geometry_time_series"),
            QStringLiteral("/trajectory/apbs_efield_time_series"),
            QStringLiteral("/trajectory/apbs_efg_time_series"),
            QStringLiteral("/trajectory/aimnet2_embedding_time_series"),
            QStringLiteral("/trajectory/aimnet2_charge_response_gradient_time_series"),
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
    case SignalSourceKind::ExperimentalShieldingMl:
        return true;
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
    const UnitSpec ringKernel = unit(UnitDimension::MagneticShielding, "ppm_T_per_nA", "ppm_T_per_nA");
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
	                       {},
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
        {"h5:bs_shielding_time_series", "bs_shielding", "biot_savart", "Biot-Savart unit-current kernel time series", "/trajectory/bs_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:hm_shielding_time_series", "hm_shielding", "haigh_mallion", "Haigh-Mallion unit-current kernel time series", "/trajectory/hm_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:mc_shielding_time_series", "mc_shielding", "mcconnell", "McConnell shielding time series", "/trajectory/mc_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:mopac_coulomb_efg_time_series", "mopac_coulomb_efg", "mopac_coulomb", "MOPAC-charge electric-field gradient time series", "/trajectory/mopac_coulomb_efg_time_series", SignalValueShape::EfgT2},
        {"h5:mopac_mc_shielding_time_series", "mopac_mc_shielding", "mopac_mcconnell", "MOPAC McConnell shielding time series", "/trajectory/mopac_mc_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:larsen_hbond_1pHB_shielding_time_series", "larsen_hbond_1pHB_shielding", "larsen_hbond", "Larsen 1pHB shielding time series", "/trajectory/larsen_hbond_1pHB_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:larsen_hbond_1pHaB_shielding_time_series", "larsen_hbond_1pHaB_shielding", "larsen_hbond", "Larsen 1pHaB shielding time series", "/trajectory/larsen_hbond_1pHaB_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:larsen_hbond_2pHB_shielding_time_series", "larsen_hbond_2pHB_shielding", "larsen_hbond", "Larsen 2pHB shielding time series", "/trajectory/larsen_hbond_2pHB_shielding_time_series", SignalValueShape::SphericalTensor},
        {"h5:larsen_hbond_2pHaB_shielding_time_series", "larsen_hbond_2pHaB_shielding", "larsen_hbond", "Larsen 2pHaB shielding time series", "/trajectory/larsen_hbond_2pHaB_shielding_time_series", SignalValueShape::SphericalTensor},
    };

    for (const TensorGroup& group : tensorGroups) {
        const bool efgOnly = group.shape == SignalValueShape::EfgT2;
        const QString conceptKey = QString::fromLatin1(group.conceptKey);
        const bool rawRingKernel = conceptKey == QStringLiteral("bs_shielding")
                                   || conceptKey == QStringLiteral("hm_shielding");
        const UnitSpec groupUnits =
            group.shape == SignalValueShape::EfgT2 ? efg
            : rawRingKernel                       ? ringKernel
                                                  : shielding;
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
                           groupUnits,
                           efgOnly ? efgStripModes() : tensorStripModes(),
                           efgOnly ? efgStaticModes() : tensorStaticModes(),
                           efgOnly ? efgChannels(groupUnits) : sphericalTensorChannels(groupUnits),
                           group.path,
                           true));
    }

    add(descriptors, makeDescriptor("h5:sasa_time_series", "atom_sasa", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "sasa", "SASA time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, unit(UnitDimension::Length, "A^2", "A^2"), scalarStripModes(), scalarStaticModes(), scalarChannels(unit(UnitDimension::Length, "A^2", "A^2")), "/trajectory/sasa_time_series", true));
    add(descriptors, makeDescriptor("h5:aimnet2_charge_time_series", "aimnet2_charges", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "aimnet2", "AIMNet2 charge time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, charge, scalarStripModes(), scalarStaticModes(), scalarChannels(charge), "/trajectory/aimnet2_charge_time_series", true));
    add(descriptors, makeDescriptor("h5:larsen_hbond_count_time_series", "larsen_hbond_count", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "larsen_hbond", "Larsen H-bond count time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.scalar")}, {}, countChannels(), "/trajectory/larsen_hbond_count_time_series", true));
    add(descriptors, makeDescriptor("h5:larsen_hbond_water_term_time_series", "larsen_hbond_water_term", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "larsen_hbond", "Larsen H-bond water term time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, shielding, scalarStripModes(), scalarStaticModes(), scalarChannels(shielding), "/trajectory/larsen_hbond_water_term_time_series", true));
    add(descriptors, makeDescriptor("h5:bonded_energy_time_series", "bonded_energy", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "bonded", "Bonded energy time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Scalar, energy, scalarStripModes(), scalarStaticModes(), scalarChannels(energy), "/trajectory/bonded_energy_time_series", true));

    add(descriptors, makeDescriptor("h5:water_field_efield_time_series", "water_efield", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "water_field", "Water electric field time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield), "/trajectory/water_field_time_series", true));
    add(descriptors, makeDescriptor("h5:water_field_efg_time_series", "water_efg", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "water_field", "Water EFG time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg), "/trajectory/water_field_time_series", true));
    add(descriptors, makeDescriptor("h5:water_shell_count_time_series", "water_shell_counts", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "water_field", "Water shell counts", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.scalar")}, {}, countChannels(), "/trajectory/water_field_time_series", true));
    add(descriptors, makeDescriptor("h5:hydration_shell_time_series", "hydration_shell", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "hydration", "Hydration shell time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {}, "/trajectory/hydration_shell_time_series", true));
    add(descriptors, makeDescriptor("h5:hydration_geometry_time_series", "hydration_geometry", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "hydration", "Hydration geometry time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::PerClassBlock, none, perClassStripModes(), perClassStaticModes(), {}, "/trajectory/hydration_geometry_time_series", true));
    add(descriptors, makeDescriptor("h5:apbs_efield_time_series", "apbs_E", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "apbs", "APBS electric field time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, efield, vectorStripModes(), vectorStaticModes(), vectorChannels(efield), "/trajectory/apbs_efield_time_series", true));
    add(descriptors, makeDescriptor("h5:apbs_efg_time_series", "apbs_efg", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "apbs", "APBS EFG time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::EfgT2, efg, efgStripModes(), efgStaticModes(), efgChannels(efg), "/trajectory/apbs_efg_time_series", true));
    add(descriptors, makeDescriptor("h5:aimnet2_embedding_time_series", "aimnet2_embedding", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "aimnet2", "AIMNet2 embedding time series", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Embedding, none, {}, {}, {}, "/trajectory/aimnet2_embedding_time_series", true));  // 256-d ML feature -> non-displayable (DisplayPolicy); offered no plottable mode
    add(descriptors, makeDescriptor("h5:aimnet2_charge_response_gradient_time_series", "aimnet2_charge_response_gradient", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "aimnet2", "AIMNet2 charge-response gradient", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, charge, vectorStripModes(), vectorStaticModes(), vectorChannels(charge), "/trajectory/aimnet2_charge_response_gradient_time_series", true));

    add(descriptors, makeDescriptor("h5:dihedral_time_series", "residue_dihedral", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "planar_geometry", "Residue dihedral time series", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::PerClassBlock, angle, perClassStripModes(), {QStringLiteral("static.newman")}, {channel("phi", "Phi", SignalValueShape::Scalar, angle), channel("psi", "Psi", SignalValueShape::Scalar, angle), channel("omega", "Omega", SignalValueShape::Scalar, angle), channel("chi", "Chi", SignalValueShape::Scalar, angle)}, "/trajectory/dihedral_time_series", true));
    add(descriptors, makeDescriptor("h5:dssp8_time_series", "dssp_ss8", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "dssp", "DSSP8 residue state time series", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::Category, tag, categoryStripModes(), {}, {channel("ss8", "SS8", SignalValueShape::Category, tag)}, "/trajectory/dssp8_time_series", true));
    add(descriptors, makeDescriptor("h5:j_coupling_time_series", "j_coupling", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "j_coupling", "J-coupling time series", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::PerClassBlock, unit(UnitDimension::Frequency, "Hz", "Hz"), perClassStripModes(), {}, {}, "/trajectory/j_coupling_time_series", true));
    add(descriptors, makeDescriptor("h5:ring_pucker_time_series", "ring_pucker", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "planar_geometry", "Ring pucker time series", SourceResidency::StartupLoaded, SignalAxis::Ring, SignalAxis::Ring, SignalValueShape::PerClassBlock, angle, perClassStripModes(), {}, {}, "/trajectory/ring_pucker_time_series", true));
    add(descriptors, makeDescriptor("h5:ring_neighbourhood_trajectory_stats", "ring_neighbourhood", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "ring_current", "Ring neighbourhood trajectory stats", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::PerClassBlock, length, perClassStripModes(), {}, {channel("distance", "Distance", SignalValueShape::Scalar, length), channel("rho", "Rho", SignalValueShape::Scalar, length), channel("z", "Z", SignalValueShape::Scalar, length), channel("in_plane_angle", "In-plane angle", SignalValueShape::Scalar, angle)}, "/trajectory/ring_neighbourhood_trajectory_stats", true));
    add(descriptors, makeDescriptor("h5:gromacs_energy_time_series", "gromacs_energy", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "gromacs", "Gromacs energy/runtime time series", SourceResidency::StartupLoaded, SignalAxis::System, SignalAxis::System, SignalValueShape::PerClassBlock, none, {QStringLiteral("strip.system"), QStringLiteral("strip.per-class"), QStringLiteral("strip.tensor.component")}, {}, {channel("energy", "Energy", SignalValueShape::Scalar, energy), channel("temperature", "Temperature", SignalValueShape::Scalar, unit(UnitDimension::Temperature, "K", "K")), channel("pressure", "Pressure", SignalValueShape::Scalar, unit(UnitDimension::Pressure, "bar", "bar")), channel("volume", "Volume", SignalValueShape::Scalar, unit(UnitDimension::Volume, "nm^3", "nm^3"))}, "/trajectory/gromacs_energy_time_series", true));
    add(descriptors, makeDescriptor("h5:rmsd_tracking", "rmsd_tracking", SignalSourceKind::DenseH5Trajectory, "TrajectoryH5", "rmsd", "RMSD tracking", SourceResidency::StartupLoaded, SignalAxis::System, SignalAxis::System, SignalValueShape::Scalar, length, {QStringLiteral("strip.system"), QStringLiteral("strip.scalar"), QStringLiteral("strip.event")}, {}, scalarChannels(length), "/trajectory/rmsd_tracking", true));

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
        {"h5:bs_welford", "bs_shielding.stats", "biot_savart", "Biot-Savart unit-current Welford rollup", "/trajectory/bs_welford", SignalAxis::Atom, ringKernel},
        {"h5:hm_welford", "hm_shielding.stats", "haigh_mallion", "Haigh-Mallion unit-current Welford rollup", "/trajectory/hm_welford", SignalAxis::Atom, ringKernel},
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
        {"h5:bs_t0_autocorrelation", "bs_shielding.autocorrelation", "biot_savart", "Biot-Savart unit-current T0 autocorrelation", "/trajectory/bs_t0_autocorrelation", SignalAxis::Atom, ringKernel},
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
                           {},   // de-stripped (DisplayPolicy): whole-traj summaries flat-line in a strip
                           rollupStaticModes(),
                           {},
                           group.path,
                           true));
    }

    // ── Per-atom × per-channel (KernelDynamics) ─────────────────────
    // The seven July kernel channels are stored by name in HDF5. Catalog declares
    // them as ChannelDescriptors so the strip-mode dispatch + the
    // panel renderers can index by channel id.
    static const struct { const char* id; const char* label; } kKernelChannels[] = {
        {"bs_T0", "BS T0"},     {"bs_absT2", "BS |T2|"},
        {"hm_T0", "HM T0"},     {"hm_absT2", "HM |T2|"},
        {"mc_T0", "MC T0"},     {"mc_absT2", "MC |T2|"},
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
                       "Kernel coherence matrix (per atom, 7x7 Pearson)",
                       SourceResidency::StartupLoaded,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::Matrix,
	                       none,
	                       {},
	                       {QStringLiteral("static.chord.coupling")},
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
	                           {QStringLiteral("static.bar.sequence")},
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
	                       {QStringLiteral("static.bar.sequence")},
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
    // - Orientation tensor (Mat3 per vector) via static.tensor -- shown as a
    //   focus-driven SCENE glyph (the shared TensorGlyphActor, the same ovaloid
    //   + arrows as the CSA glyph), not a dashboard panel; static.tensor stays
    //   tracked-but-hidden in the dashboard.
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
	                           {QStringLiteral("static.bar.sequence")},
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
    // payload feeds the tracked-hidden static.tensor mode; the 3-D glyph
    // is focus-driven by ReaderMainWindow rather than emitted as a dashboard
    // panel.
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
	                       {QStringLiteral("static.tensor")},
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
	                       {QStringLiteral("static.fixed_freq")},
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
	                       {QStringLiteral("static.bar.sequence")},
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
    auto npy = [&](io::FieldKind field,
                   const char* family,
                   const char* label,
                   SignalValueShape shape,
                   const QStringList& temporalModes,
                   const QStringList& staticModes,
                   const QVector<ChannelDescriptor>& channels) {
        const io::FieldSpec& spec = io::FieldSpecFor(field);
        const std::optional<SignalAxis> axis = signalAxisFor(spec.axis);
        Q_ASSERT(axis);
        Q_ASSERT(io::ShouldLoadFrameField(field));

        const UnitSpec units = producerUnits(field);
        SignalDescriptor descriptor = makeDescriptor("npy:field",
                                                     "field",
                                                     SignalSourceKind::FrameNpySnapshot,
                                                     "SDK_NPY",
                                                     family,
                                                     label,
                                                     SourceResidency::FrameLoaded,
                                                     axis.value_or(SignalAxis::None),
                                                     axis.value_or(SignalAxis::None),
                                                     shape,
                                                     units,
                                                     temporalModes,
                                                     staticModes,
                                                     channels,
                                                     "",
                                                     false,
                                                     true);
        const QString stem = fromUtf8(spec.stem);
        descriptor.id = QStringLiteral("npy:%1").arg(stem);
        descriptor.conceptKey = stem;
        descriptor.storagePath = stem;
        descriptor.description = fromUtf8(spec.description);
        descriptor.samplingStatus = SampleStatus::Valid;
        descriptor.samplingGapReason = GapReason::None;
        add(descriptors, descriptor);
    };

    const struct TensorField {
        io::FieldKind field;
        const char* family;
        const char* label;
    } tensorFields[] = {
        {io::FieldKind::BSShielding, "biot_savart", "Biot-Savart unit-current response"},
        {io::FieldKind::HMShielding, "haigh_mallion", "Haigh-Mallion geometry response"},
        {io::FieldKind::MOPACMcShielding, "mopac_mcconnell", "MOPAC-weighted McConnell response"},
        {io::FieldKind::LarsenHBondShielding, "larsen_hbond", "Larsen H-bond total"},
        {io::FieldKind::LarsenHBond1pHBShielding, "larsen_hbond", "Larsen primary amide-H contribution"},
        {io::FieldKind::LarsenHBond2pHBShielding, "larsen_hbond", "Larsen secondary amide-H contribution"},
        {io::FieldKind::LarsenHBond1pHaBShielding, "larsen_hbond", "Larsen primary alpha-H contribution"},
        {io::FieldKind::LarsenHBond2pHaBShielding, "larsen_hbond", "Larsen secondary alpha-H contribution"},
    };

    for (const TensorField& field : tensorFields) {
        const UnitSpec units = producerUnits(field.field);
        npy(field.field,
            field.family,
            field.label,
            SignalValueShape::SphericalTensor,
            tensorStripModes(),
            tensorStaticModes(),
            sphericalTensorChannels(units));
    }

    static constexpr std::array<NamedChannel, 8> kRingTypes{{
        {"phe_benzene", "PHE benzene"},
        {"tyr_phenol", "TYR phenol"},
        {"trp_benzene", "TRP benzene"},
        {"trp_pyrrole", "TRP pyrrole"},
        {"trp_perimeter", "TRP perimeter"},
        {"his_imidazole", "HIS imidazole"},
        {"hid_imidazole", "HID imidazole"},
        {"hie_imidazole", "HIE imidazole"},
    }};

    const struct RingTypeField {
        io::FieldKind field;
        const char* family;
        const char* label;
        SignalValueShape shape;
    } ringTypeFields[] = {
        {io::FieldKind::BSPerTypeT0, "biot_savart", "Biot-Savart response by ring type (T0)", SignalValueShape::PerClassBlock},
        {io::FieldKind::BSPerTypeT1, "biot_savart", "Biot-Savart response by ring type (T1)", SignalValueShape::Vector3},
        {io::FieldKind::BSPerTypeT2, "biot_savart", "Biot-Savart response by ring type (T2)", SignalValueShape::EfgT2},
        {io::FieldKind::HMPerTypeT0, "haigh_mallion", "Haigh-Mallion response by ring type (T0)", SignalValueShape::PerClassBlock},
        {io::FieldKind::HMPerTypeT1, "haigh_mallion", "Haigh-Mallion response by ring type (T1)", SignalValueShape::Vector3},
        {io::FieldKind::HMPerTypeT2, "haigh_mallion", "Haigh-Mallion response by ring type (T2)", SignalValueShape::EfgT2},
        {io::FieldKind::PQPerTypeT0, "pi_quadrupole", "Pi-quadrupole response by ring type", SignalValueShape::PerClassBlock},
        {io::FieldKind::RingChiPerTypeT0, "ring_susceptibility", "Ring-susceptibility geometry by ring type", SignalValueShape::PerClassBlock},
        {io::FieldKind::DispPerTypeT0, "dispersion", "Dispersion proximity by ring type", SignalValueShape::PerClassBlock},
    };

    for (const RingTypeField& field : ringTypeFields) {
        const UnitSpec units = producerUnits(field.field);
        if (field.shape == SignalValueShape::Vector3) {
            npy(field.field,
                field.family,
                field.label,
                field.shape,
                vectorStripModes(),
                vectorStaticModes(),
                vectorBlockChannels(kRingTypes, units));
        } else if (field.shape == SignalValueShape::EfgT2) {
            npy(field.field,
                field.family,
                field.label,
                field.shape,
                efgStripModes(),
                efgStaticModes(),
                t2BlockChannels(kRingTypes, units));
        } else {
            npy(field.field,
                field.family,
                field.label,
                field.shape,
                perClassStripModes(),
                perClassStaticModes(),
                scalarBlockChannels(kRingTypes, units));
        }
    }

    auto addScalar = [&](io::FieldKind field, const char* family, const char* label) {
        const UnitSpec units = producerUnits(field);
        npy(field,
            family,
            label,
            SignalValueShape::Scalar,
            scalarStripModes(),
            scalarStaticModes(),
            scalarChannels(units));
    };
    auto addCount = [&](io::FieldKind field, const char* family, const char* label) {
        npy(field,
            family,
            label,
            SignalValueShape::Count,
            {QStringLiteral("strip.count"), QStringLiteral("strip.scalar")},
            scalarStaticModes(),
            countChannels());
    };
    auto addVector = [&](io::FieldKind field, const char* family, const char* label) {
        const UnitSpec units = producerUnits(field);
        npy(field,
            family,
            label,
            SignalValueShape::Vector3,
            vectorStripModes(),
            vectorStaticModes(),
            vectorChannels(units));
    };
    auto addT2 = [&](io::FieldKind field, const char* family, const char* label) {
        const UnitSpec units = producerUnits(field);
        npy(field,
            family,
            label,
            SignalValueShape::EfgT2,
            efgStripModes(),
            efgStaticModes(),
            efgChannels(units));
    };
    auto addTensor = [&](io::FieldKind field, const char* family, const char* label) {
        const UnitSpec units = producerUnits(field);
        npy(field,
            family,
            label,
            SignalValueShape::SphericalTensor,
            tensorStripModes(),
            tensorStaticModes(),
            sphericalTensorChannels(units));
    };
    auto addTracelessTensor = [&](io::FieldKind field,
                                  const char* family,
                                  const char* label) {
        const UnitSpec units = producerUnits(field);
        npy(field,
            family,
            label,
            SignalValueShape::SphericalTensor,
            tracelessTensorStripModes(),
            tensorStaticModes(),
            sphericalTensorChannels(units));
    };
    auto addScalarBlock = [&](io::FieldKind field,
                              const char* family,
                              const char* label,
                              const QVector<ChannelDescriptor>& channels) {
        npy(field,
            family,
            label,
            SignalValueShape::PerClassBlock,
            perClassStripModes(),
            perClassStaticModes(),
            channels);
    };

    addVector(io::FieldKind::BSTotalB, "biot_savart", "Biot-Savart total magnetic field");
    static constexpr std::array<NamedChannel, 4> kRingCountShells{{
        {"within_3A", "Within 3 A"},
        {"within_5A", "Within 5 A"},
        {"within_8A", "Within 8 A"},
        {"within_12A", "Within 12 A"},
    }};
    addScalarBlock(io::FieldKind::BSRingCounts,
                   "biot_savart",
                   "Accepted ring counts",
                   scalarBlockChannels(kRingCountShells,
                                       producerUnits(io::FieldKind::BSRingCounts)));

    const struct ScalarField {
        io::FieldKind field;
        const char* family;
        const char* label;
    } scalarFields[] = {
        {io::FieldKind::TauNCAC, "local_geometry", "N-CA-C valence angle"},
        {io::FieldKind::AngleNCACB, "local_geometry", "N-CA-CB valence angle"},
        {io::FieldKind::AngleCBCAC, "local_geometry", "CB-CA-C valence angle"},
        {io::FieldKind::AngleCprevNCA, "local_geometry", "Previous C-N-CA valence angle"},
        {io::FieldKind::AngleCACNnext, "local_geometry", "CA-C-next N valence angle"},
        {io::FieldKind::CbDeviation, "local_geometry", "C-beta ideal-position deviation"},
        {io::FieldKind::AtomSASA, "sasa", "Atomic solvent-accessible area"},
        {io::FieldKind::AtomSASAFraction, "sasa", "Atomic solvent exposure fraction"},
        {io::FieldKind::FfPartialCharge, "force_field", "Prepared force-field charge"},
        {io::FieldKind::FfPbRadius, "force_field", "Poisson-Boltzmann radius"},
        {io::FieldKind::EEQCharges, "eeq", "EEQ charge"},
        {io::FieldKind::EEQCN, "eeq", "EEQ coordination number"},
        {io::FieldKind::EEQChiEff, "eeq", "EEQ effective electronegativity"},
        {io::FieldKind::APBSPhi, "apbs", "APBS reaction potential"},
        {io::FieldKind::AIMNet2Charges, "aimnet2", "AIMNet2 charge"},
        {io::FieldKind::AIMNet2ChargeResponseGradientScalar, "aimnet2", "AIMNet2 charge-response-gradient magnitude"},
        {io::FieldKind::AIMNet2EnergyMlp, "aimnet2", "AIMNet2 local learned energy"},
        {io::FieldKind::AIMNet2EnergyShiftedLocal, "aimnet2", "AIMNet2 shifted local energy"},
        {io::FieldKind::AIMNet2D3EDispAtom, "aimnet2", "AIMNet2 D3 atomic dispersion energy"},
        {io::FieldKind::AIMNet2D3CN, "aimnet2", "AIMNet2 D3 coordination number"},
        {io::FieldKind::Pyramidalization, "planar_geometry", "Out-of-plane displacement"},
        {io::FieldKind::OmegaActual, "planar_geometry", "Peptide omega"},
        {io::FieldKind::OmegaDeviation, "planar_geometry", "Peptide omega deviation from trans"},
        {io::FieldKind::AromaticChi2, "planar_geometry", "Aromatic chi2"},
        {io::FieldKind::PuckerQ, "planar_geometry", "Ring puckering amplitude"},
        {io::FieldKind::PuckerTheta, "planar_geometry", "Ring puckering phase"},
        {io::FieldKind::MOPACChargesFullPrecision, "mopac", "MOPAC Coulson charge"},
        {io::FieldKind::MOPACBondValenciesFullPrecision, "mopac", "MOPAC bond valency"},
        {io::FieldKind::MOPACAtomSPopulation, "mopac", "MOPAC atomic s population"},
        {io::FieldKind::MOPACAtomPPopulation, "mopac", "MOPAC atomic p population"},
        {io::FieldKind::MOPACAtomDPopulation, "mopac", "MOPAC atomic d population"},
        {io::FieldKind::MOPACLewisBondCount, "mopac", "MOPAC Lewis-bond count"},
        {io::FieldKind::MOPACMcCoSum, "mopac_mcconnell", "MOPAC McConnell C=O scalar sum"},
        {io::FieldKind::MOPACMcCNSum, "mopac_mcconnell", "MOPAC McConnell C-N scalar sum"},
        {io::FieldKind::MOPACMcSidechainSum, "mopac_mcconnell", "MOPAC McConnell side-chain scalar sum"},
        {io::FieldKind::MOPACMcAromaticSum, "mopac_mcconnell", "MOPAC McConnell aromatic scalar sum"},
        {io::FieldKind::MOPACMcNearestCoDist, "mopac_mcconnell", "Nearest MOPAC-weighted C=O distance"},
        {io::FieldKind::MOPACMcNearestCNDist, "mopac_mcconnell", "Nearest MOPAC-weighted C-N distance"},
        {io::FieldKind::LarsenHBondWaterTerm, "larsen_hbond", "Larsen water term"},
    };
    for (const ScalarField& field : scalarFields)
        addScalar(field.field, field.family, field.label);
    addCount(io::FieldKind::LarsenHBondCount, "larsen_hbond", "Larsen contributing-pair count");

    const struct VectorField {
        io::FieldKind field;
        const char* family;
        const char* label;
    } vectorFields[] = {
        {io::FieldKind::CbResidualVector, "local_geometry", "C-beta residual vector"},
        {io::FieldKind::SASANormal, "sasa", "Solvent-accessibility normal"},
        {io::FieldKind::HBondNearestDir, "hbond", "Nearest accepted H-bond direction"},
        {io::FieldKind::McNearestCoDir, "mcconnell", "Nearest C=O direction"},
        {io::FieldKind::CoulombE, "coulomb", "Force-field electric field"},
        {io::FieldKind::CoulombEBackbone, "coulomb", "Force-field backbone electric field"},
        {io::FieldKind::CoulombESidechain, "coulomb", "Force-field side-chain electric field"},
        {io::FieldKind::CoulombEAromatic, "coulomb", "Force-field aromatic electric field"},
        {io::FieldKind::EEQCoulombE, "eeq_coulomb", "EEQ electric field"},
        {io::FieldKind::MOPACCoulombE, "mopac_coulomb", "MOPAC-charge electric field"},
        {io::FieldKind::MOPACCoulombEBackbone, "mopac_coulomb", "MOPAC-charge backbone electric field"},
        {io::FieldKind::MOPACCoulombESidechain, "mopac_coulomb", "MOPAC-charge side-chain electric field"},
        {io::FieldKind::MOPACCoulombEAromatic, "mopac_coulomb", "MOPAC-charge aromatic electric field"},
        {io::FieldKind::APBSE, "apbs", "APBS reaction electric field"},
        {io::FieldKind::AIMNet2E, "aimnet2", "AIMNet2-charge electric field"},
        {io::FieldKind::AIMNet2EBackbone, "aimnet2", "AIMNet2-charge backbone electric field"},
        {io::FieldKind::AIMNet2ESidechain, "aimnet2", "AIMNet2-charge side-chain electric field"},
        {io::FieldKind::AIMNet2EAromatic, "aimnet2", "AIMNet2-charge aromatic electric field"},
        {io::FieldKind::AIMNet2ChargeResponseGradient, "aimnet2", "AIMNet2 charge-response gradient"},
        {io::FieldKind::WaterEfield, "water_field", "Water electric field"},
        {io::FieldKind::WaterEfieldFirst, "water_field", "First-shell water electric field"},
    };
    for (const VectorField& field : vectorFields)
        addVector(field.field, field.family, field.label);

    addTracelessTensor(io::FieldKind::CoulombEFG,
                       "coulomb",
                       "Force-field electric-field gradient (full tensor)");
    addTracelessTensor(io::FieldKind::EEQCoulombEFG,
                       "eeq_coulomb",
                       "EEQ electric-field gradient (full tensor)");
    addTracelessTensor(io::FieldKind::MOPACCoulombEFG,
                       "mopac_coulomb",
                       "MOPAC-charge electric-field gradient (full tensor)");

    const struct T2Field {
        io::FieldKind field;
        const char* family;
        const char* label;
    } t2Fields[] = {
        {io::FieldKind::CoulombEFGT2, "coulomb", "Force-field electric-field gradient"},
        {io::FieldKind::CoulombEFGBackbone, "coulomb", "Force-field backbone EFG"},
        {io::FieldKind::CoulombEFGSidechain, "coulomb", "Force-field side-chain EFG"},
        {io::FieldKind::CoulombEFGAromatic, "coulomb", "Force-field aromatic EFG"},
        {io::FieldKind::EEQCoulombEFGBackbone, "eeq_coulomb", "EEQ backbone EFG"},
        {io::FieldKind::EEQCoulombEFGSidechain, "eeq_coulomb", "EEQ side-chain EFG"},
        {io::FieldKind::EEQCoulombEFGAromatic, "eeq_coulomb", "EEQ aromatic EFG"},
        {io::FieldKind::MOPACCoulombEFGBackbone, "mopac_coulomb", "MOPAC-charge backbone EFG"},
        {io::FieldKind::MOPACCoulombEFGSidechain, "mopac_coulomb", "MOPAC-charge side-chain EFG"},
        {io::FieldKind::MOPACCoulombEFGAromatic, "mopac_coulomb", "MOPAC-charge aromatic EFG"},
        {io::FieldKind::APBSEFG, "apbs", "APBS reaction electric-field gradient"},
        {io::FieldKind::AIMNet2EFG, "aimnet2", "AIMNet2-charge electric-field gradient"},
        {io::FieldKind::AIMNet2EFGBackbone, "aimnet2", "AIMNet2-charge backbone EFG"},
        {io::FieldKind::AIMNet2EFGSidechain, "aimnet2", "AIMNet2-charge side-chain EFG"},
        {io::FieldKind::AIMNet2EFGAromatic, "aimnet2", "AIMNet2-charge aromatic EFG"},
        {io::FieldKind::WaterEFG, "water_field", "Water electric-field gradient"},
        {io::FieldKind::WaterEFGFirst, "water_field", "First-shell water electric-field gradient"},
    };
    for (const T2Field& field : t2Fields)
        addT2(field.field, field.family, field.label);

    const struct McTensorField {
        io::FieldKind field;
        const char* label;
    } mcTensorFields[] = {
        {io::FieldKind::McPeptideCoFixed, "Peptide C=O fixed anisotropy"},
        {io::FieldKind::McPeptideCNFixed, "Peptide C-N fixed anisotropy"},
        {io::FieldKind::McBackboneOtherFixed, "Other backbone fixed anisotropy"},
        {io::FieldKind::McSidechainCoFixed, "Side-chain C=O fixed anisotropy"},
        {io::FieldKind::McSidechainOtherFixed, "Other side-chain fixed anisotropy"},
        {io::FieldKind::McDisulfideFixed, "Disulfide fixed anisotropy"},
        {io::FieldKind::McAromaticFixed, "Aromatic fixed anisotropy"},
        {io::FieldKind::McBackboneXhFixed, "Backbone X-H fixed anisotropy"},
        {io::FieldKind::McSidechainXhFixed, "Side-chain X-H fixed anisotropy"},
        {io::FieldKind::McSHFixed, "S-H fixed anisotropy"},
        {io::FieldKind::McPeptideCoRhombic, "Peptide C=O rhombic contribution"},
        {io::FieldKind::McNearestCoT2, "Nearest C=O fixed anisotropy"},
        {io::FieldKind::McNearestCNT2, "Nearest C-N fixed anisotropy"},
    };
    for (const McTensorField& field : mcTensorFields)
        addTensor(field.field, "mcconnell", field.label);

    const struct MopacMcTensorField {
        io::FieldKind field;
        const char* label;
    } mopacMcTensorFields[] = {
        {io::FieldKind::MOPACMcNearestCoT2, "Nearest MOPAC-weighted C=O anisotropy"},
        {io::FieldKind::MOPACMcNearestCNT2, "Nearest MOPAC-weighted C-N anisotropy"},
        {io::FieldKind::MOPACMcBackboneTotal, "MOPAC-weighted backbone anisotropy"},
        {io::FieldKind::MOPACMcSidechainTotal, "MOPAC-weighted side-chain anisotropy"},
        {io::FieldKind::MOPACMcAromaticTotal, "MOPAC-weighted aromatic anisotropy"},
    };
    for (const MopacMcTensorField& field : mopacMcTensorFields)
        addTensor(field.field, "mopac_mcconnell", field.label);

    static constexpr std::array<NamedChannel, 2> kNearfieldCounts{{
        {"accepted", "Accepted"},
        {"rejected", "Rejected"},
    }};
    addScalarBlock(io::FieldKind::McNearfieldCounts,
                   "mcconnell",
                   "McConnell near-field source counts",
                   scalarBlockChannels(kNearfieldCounts,
                                       producerUnits(io::FieldKind::McNearfieldCounts)));

    static constexpr std::array<NamedChannel, 4> kFieldScalars{{
        {"magnitude", "Total magnitude"},
        {"parent_projection", "Parent-to-H projection"},
        {"backbone_projection", "Backbone projection"},
        {"aromatic_magnitude", "Aromatic magnitude"},
    }};
    for (const auto field : {io::FieldKind::CoulombScalars,
                             io::FieldKind::EEQCoulombScalars,
                             io::FieldKind::MOPACCoulombScalars}) {
        const char* family = field == io::FieldKind::CoulombScalars
                                 ? "coulomb"
                                 : (field == io::FieldKind::EEQCoulombScalars
                                        ? "eeq_coulomb"
                                        : "mopac_coulomb");
        const char* label = field == io::FieldKind::CoulombScalars
                                ? "Force-field electric-field summary"
                                : (field == io::FieldKind::EEQCoulombScalars
                                       ? "EEQ electric-field summary"
                                       : "MOPAC-charge electric-field summary");
        addScalarBlock(field,
                       family,
                       label,
                       scalarBlockChannels(kFieldScalars, producerUnits(field)));
    }

    static constexpr std::array<NamedChannel, 4> kHBondScalars{{
        {"nearest_distance", "Nearest H-to-target distance"},
        {"nearest_inverse_cube", "Nearest inverse-cube distance"},
        {"accepted_count", "Accepted source count"},
        {"angular_sum", "Angular inverse-cube sum"},
    }};
    const UnitSpec angstrom = unit(UnitDimension::Length, "A", "A");
    const UnitSpec inverseCube = unit(UnitDimension::Dimensionless, "A^-3", "A^-3", 1.0, 0.0, false);
    const UnitSpec count = unit(UnitDimension::Count, "count", "count");
    const QVector<ChannelDescriptor> hbondChannels{
        channel(kHBondScalars[0].id, kHBondScalars[0].label, SignalValueShape::Scalar, angstrom, {}, 0, 1),
        channel(kHBondScalars[1].id, kHBondScalars[1].label, SignalValueShape::Scalar, inverseCube, {}, 1, 1),
        channel(kHBondScalars[2].id, kHBondScalars[2].label, SignalValueShape::Count, count, {}, 2, 1),
        channel(kHBondScalars[3].id, kHBondScalars[3].label, SignalValueShape::Scalar, inverseCube, {}, 3, 1),
    };
    addScalarBlock(io::FieldKind::HBondScalars,
                   "hbond",
                   "Hydrogen-bond geometry summary",
                   hbondChannels);

    static constexpr std::array<NamedChannel, 8> kDsspClasses{{
        {"H", "Alpha helix"},
        {"G", "3-10 helix"},
        {"I", "Pi helix"},
        {"E", "Extended strand"},
        {"B", "Beta bridge"},
        {"T", "Turn"},
        {"S", "Bend"},
        {"C", "Coil"},
    }};
    addScalarBlock(io::FieldKind::DSSPSs8,
                   "dssp",
                   "DSSP secondary-structure class",
                   scalarBlockChannels(kDsspClasses, producerUnits(io::FieldKind::DSSPSs8)));
    addScalar(io::FieldKind::DSSPObserved, "dssp", "DSSP observed-state mask");
    addScalar(io::FieldKind::DSSPPpii, "dssp", "DSSP polyproline-II state");

    static constexpr std::array<NamedChannel, 4> kDsspHBondEnergies{{
        {"acceptor_0", "Acceptor 1"},
        {"acceptor_1", "Acceptor 2"},
        {"donor_0", "Donor 1"},
        {"donor_1", "Donor 2"},
    }};
    addScalarBlock(io::FieldKind::DSSPHBondEnergy,
                   "dssp",
                   "DSSP Coulombic H-bond assignment scores",
                   scalarBlockChannels(kDsspHBondEnergies,
                                       producerUnits(io::FieldKind::DSSPHBondEnergy)));

    static constexpr std::array<NamedChannel, 6> kTorsions{{
        {"phi", "Phi"},
        {"psi", "Psi"},
        {"chi1", "Chi1"},
        {"chi2", "Chi2"},
        {"chi3", "Chi3"},
        {"chi4", "Chi4"},
    }};
    addScalarBlock(io::FieldKind::DSSPTorsionAngle,
                   "dssp",
                   "Signed IUPAC torsion angles",
                   scalarBlockChannels(kTorsions,
                                       producerUnits(io::FieldKind::DSSPTorsionAngle)));

    static constexpr std::array<NamedChannel, 2> kWaterShellCounts{{
        {"first_shell", "First shell"},
        {"second_shell", "Second shell"},
    }};
    addScalarBlock(io::FieldKind::WaterShellCounts,
                   "water_field",
                   "Water-shell counts",
                   scalarBlockChannels(kWaterShellCounts,
                                       producerUnits(io::FieldKind::WaterShellCounts)));

    static constexpr std::array<NamedChannel, 2> kEeqHardness{{
        {"eta", "Eta"},
        {"self_hardness", "Self-hardness diagonal"},
    }};
    addScalarBlock(io::FieldKind::EEQHardness,
                   "eeq",
                   "EEQ hardness",
                   scalarBlockChannels(kEeqHardness,
                                       producerUnits(io::FieldKind::EEQHardness)));

    static constexpr std::array<NamedChannel, 3> kD3C6Stats{{
        {"sum", "C6 sum"},
        {"mean", "C6 mean"},
        {"maximum", "C6 maximum"},
    }};
    addScalarBlock(io::FieldKind::AIMNet2D3C6Stats,
                   "aimnet2",
                   "AIMNet2 D3 C6 summary",
                   scalarBlockChannels(kD3C6Stats,
                                       producerUnits(io::FieldKind::AIMNet2D3C6Stats)));
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
	                       {QStringLiteral("strip.tensor.T0")},
	                       {QStringLiteral("static.tensor")},
	                       sphericalTensorChannels(shielding),
                       "orca_total",
                       false,
                       true,
                       SampleStatus::Valid,
                       GapReason::None));
    add(descriptors, makeDescriptor("orca_dft:diamagnetic", "orca_diamagnetic", SignalSourceKind::OrcaDftFrame, "ORCA", "orca", "ORCA diamagnetic shielding", SourceResidency::FrameLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding), "orca_diamagnetic", false, true));
    add(descriptors, makeDescriptor("orca_dft:paramagnetic", "orca_paramagnetic", SignalSourceKind::OrcaDftFrame, "ORCA", "orca", "ORCA paramagnetic shielding", SourceResidency::FrameLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::SphericalTensor, shielding, tensorStripModes(), tensorStaticModes(), sphericalTensorChannels(shielding), "orca_paramagnetic", false, true));
}

void addExperimentalShieldingMl(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec shielding = unit(UnitDimension::MagneticShielding, "ppm", "ppm");
    add(descriptors,
        makeDescriptor("ml:experimental_shielding_iso",
                       "experimental_shielding_ml.iso",
                       SignalSourceKind::ExperimentalShieldingMl,
                       "Experimental Shielding ML",
                       "experimental_shielding_ml",
                       "Experimental ML isotropic shielding",
                       SourceResidency::Derived,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::Scalar,
                       shielding,
                       scalarStripModes(),
                       scalarStaticModes(),
                       scalarChannels(shielding),
                       "experimental_shielding_ml:iso",
                       false,
                       true));
    add(descriptors,
        makeDescriptor("ml:experimental_shielding_t2",
                       "experimental_shielding_ml.t2",
                       SignalSourceKind::ExperimentalShieldingMl,
                       "Experimental Shielding ML",
                       "experimental_shielding_ml",
                       "Experimental ML shielding tensor",
                       SourceResidency::Derived,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::EfgT2,
                       shielding,
                       efgStripModes(),
                       efgStaticModes(),
                       efgChannels(shielding),
                       "experimental_shielding_ml:t2",
                       false,
                       true));
    add(descriptors,
        makeDescriptor("ml:experimental_shielding_t2_norm",
                       "experimental_shielding_ml.t2_norm",
                       SignalSourceKind::ExperimentalShieldingMl,
                       "Experimental Shielding ML",
                       "experimental_shielding_ml",
                       "Experimental ML shielding |T2|",
                       SourceResidency::Derived,
                       SignalAxis::Atom,
                       SignalAxis::Atom,
                       SignalValueShape::Scalar,
                       shielding,
                       scalarStripModes(),
                       scalarStaticModes(),
                       scalarChannels(shielding),
                       "experimental_shielding_ml:t2_norm",
                       false,
                       true));
}

void addTopology(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec tag = unit(UnitDimension::Tag, "tag", "tag", 1.0, 0.0, false);
    const UnitSpec length = unit(UnitDimension::Length, "A", "A");
    add(descriptors, makeDescriptor("topology:atoms", "topology.atoms", SignalSourceKind::Topology, "TopologySidecar", "topology", "Atoms", SourceResidency::StartupLoaded, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Category, tag, {}, {}, {}, "atoms"));
    add(descriptors, makeDescriptor("topology:residues", "topology.residues", SignalSourceKind::Topology, "TopologySidecar", "topology", "Residues", SourceResidency::StartupLoaded, SignalAxis::Residue, SignalAxis::Residue, SignalValueShape::Category, tag, {}, {}, {}, "residues"));
    add(descriptors, makeDescriptor("topology:bonds", "topology.bonds", SignalSourceKind::Topology, "TopologySidecar", "topology", "Bonds", SourceResidency::StartupLoaded, SignalAxis::Bond, SignalAxis::Bond, SignalValueShape::Category, tag, {}, {}, {}, "bonds"));
    add(descriptors, makeDescriptor("topology:rings", "topology.rings", SignalSourceKind::Topology, "TopologySidecar", "topology", "Rings", SourceResidency::StartupLoaded, SignalAxis::Ring, SignalAxis::Ring, SignalValueShape::Category, tag, {}, {}, {}, "rings"));
    add(descriptors, makeDescriptor("topology:ring_membership", "topology.ring_membership", SignalSourceKind::Topology, "TopologySidecar", "topology", "Ring membership", SourceResidency::StartupLoaded, SignalAxis::RingMembership, SignalAxis::RingMembership, SignalValueShape::Category, tag, {}, {}, {}, "ring_membership"));
    add(descriptors, makeDescriptor("topology:bond_length", "geometry.bond_length", SignalSourceKind::Topology, "DerivedTopology", "topology", "Bond length from topology positions", SourceResidency::Derived, SignalAxis::Bond, SignalAxis::Bond, SignalValueShape::Scalar, length, {QStringLiteral("strip.scalar"), QStringLiteral("strip.rollup")}, {}, scalarChannels(length), "bond_length"));
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
	                       {},
	                       scalarChannels(length),
                       "distance",
                       false,
                       false,
                       SampleStatus::Valid,
                       GapReason::None));
    add(descriptors, makeDescriptor("geometry:angle", "geometry.angle", SignalSourceKind::DerivedGeometry, "Derived", "geometry", "Angle", SourceResidency::Derived, SignalAxis::AtomTuple, SignalAxis::AtomTuple, SignalValueShape::Scalar, degrees, scalarStripModes(), {}, angleChannels(), "angle", false, false, SampleStatus::Valid, GapReason::None));
    add(descriptors, makeDescriptor("geometry:dihedral", "geometry.dihedral", SignalSourceKind::DerivedGeometry, "Derived", "geometry", "Dihedral", SourceResidency::Derived, SignalAxis::AtomTuple, SignalAxis::AtomTuple, SignalValueShape::Scalar, degrees, scalarStripModes(), {}, angleChannels(), "dihedral", false, false, SampleStatus::Valid, GapReason::None));
    add(descriptors, makeDescriptor("geometry:atom_displacement", "geometry.atom_displacement", SignalSourceKind::DerivedGeometry, "Derived", "geometry", "Atom displacement", SourceResidency::Derived, SignalAxis::Atom, SignalAxis::Atom, SignalValueShape::Vector3, length, vectorStripModes(), vectorStaticModes(), vectorChannels(length), "atom_displacement"));
}

void addSelectionEvents(QVector<SignalDescriptor>& descriptors) {
    const UnitSpec tag = unit(UnitDimension::Tag, "tag", "tag", 1.0, 0.0, false);
    add(descriptors, makeDescriptor("events:selection_timeline", "selections", SignalSourceKind::SelectionEvents, "SelectionBag", "selections", "Selection timeline", SourceResidency::StartupLoaded, SignalAxis::Event, SignalAxis::Event, SignalValueShape::EventRecord, tag, eventStripModes(), eventStaticModes(), {}, "selection_timeline"));
    add(descriptors, makeDescriptor("events:selection_counts", "selection_counts", SignalSourceKind::SelectionEvents, "SelectionBag", "selections", "Selection counts", SourceResidency::StartupLoaded, SignalAxis::System, SignalAxis::System, SignalValueShape::Count, unit(UnitDimension::Count, "count", "count"), {QStringLiteral("strip.count"), QStringLiteral("strip.system")}, {}, countChannels(), "selection_counts"));
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
    , allDescriptors_(BuildDescriptorCatalog())
    , descriptors_(allDescriptors_) {}

void TrajectorySignalCatalog::setFieldAvailability(
    std::shared_ptr<const TrajectoryFieldAvailability> availability) {
    availability_ = std::move(availability);
    rebuildVisibleDescriptors();
}

void TrajectorySignalCatalog::rebuildVisibleDescriptors() {
    descriptors_.clear();
    descriptors_.reserve(allDescriptors_.size());
    for (const SignalDescriptor& descriptor : allDescriptors_) {
        if (!availability_ || availability_->allowsDescriptor(descriptor))
            descriptors_.push_back(descriptor);
    }
}

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
        && canBind(binding) && descriptor->samplingStatus == SampleStatus::Valid
        && (!availability_ || availability_->canSampleDescriptor(*descriptor));
}

QVector<SignalDescriptor> TrajectorySignalCatalog::BuildDescriptorCatalog() {
    QVector<SignalDescriptor> descriptors;
    descriptors.reserve(160);
    addDenseH5(descriptors);
    addFrameNpy(descriptors);
    addOrcaDft(descriptors);
    addExperimentalShieldingMl(descriptors);
    addTopology(descriptors);
    addDerivedGeometry(descriptors);
    addSelectionEvents(descriptors);
    return descriptors;
}

}  // namespace h5reader::model
