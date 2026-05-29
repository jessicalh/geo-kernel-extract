#include "DashboardSignal.h"

#include <QStringList>

#include <type_traits>

namespace h5reader::model {

namespace {

template <class... Ts>
struct Overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
Overloaded(Ts...) -> Overloaded<Ts...>;

QString indexLabel(const char* prefix, std::size_t index) {
    return QStringLiteral("%1 %2").arg(QString::fromLatin1(prefix)).arg(static_cast<qulonglong>(index));
}

}  // namespace

QString ToString(SignalSourceKind kind) {
    switch (kind) {
    case SignalSourceKind::DenseH5Trajectory:
        return QStringLiteral("Dense H5 trajectory");
    case SignalSourceKind::FrameNpySnapshot:
        return QStringLiteral("Frame NPY snapshot");
    case SignalSourceKind::OrcaDftFrame:
        return QStringLiteral("ORCA DFT frame");
    case SignalSourceKind::Topology:
        return QStringLiteral("Topology");
    case SignalSourceKind::DerivedGeometry:
        return QStringLiteral("Derived geometry");
    case SignalSourceKind::SelectionEvents:
        return QStringLiteral("Selections/events");
    }
    return QStringLiteral("Unknown source");
}

QString ToString(SourceResidency residency) {
    switch (residency) {
    case SourceResidency::StartupLoaded:
        return QStringLiteral("startup-loaded");
    case SourceResidency::FrameLoaded:
        return QStringLiteral("frame-loaded");
    case SourceResidency::Derived:
        return QStringLiteral("derived");
    }
    return QStringLiteral("unknown");
}

QString ToString(SignalAxis axis) {
    switch (axis) {
    case SignalAxis::None:
        return QStringLiteral("none");
    case SignalAxis::Atom:
        return QStringLiteral("atom");
    case SignalAxis::Residue:
        return QStringLiteral("residue");
    case SignalAxis::AtomTuple:
        return QStringLiteral("atom tuple");
    case SignalAxis::Bond:
        return QStringLiteral("bond");
    case SignalAxis::BondVector:
        return QStringLiteral("bond vector");
    case SignalAxis::Ring:
        return QStringLiteral("ring");
    case SignalAxis::AromaticRing:
        return QStringLiteral("aromatic ring");
    case SignalAxis::SaturatedRing:
        return QStringLiteral("saturated ring");
    case SignalAxis::RingContributionPair:
        return QStringLiteral("ring contribution pair");
    case SignalAxis::RingMembership:
        return QStringLiteral("ring membership");
    case SignalAxis::MutationMatchPair:
        return QStringLiteral("mutation match pair");
    case SignalAxis::Protein:
        return QStringLiteral("protein");
    case SignalAxis::System:
        return QStringLiteral("system");
    case SignalAxis::Event:
        return QStringLiteral("event");
    }
    return QStringLiteral("unknown");
}

QString ToString(SignalValueShape shape) {
    switch (shape) {
    case SignalValueShape::Scalar:
        return QStringLiteral("scalar");
    case SignalValueShape::Count:
        return QStringLiteral("count");
    case SignalValueShape::Category:
        return QStringLiteral("category");
    case SignalValueShape::Vector3:
        return QStringLiteral("vector3");
    case SignalValueShape::EfgT2:
        return QStringLiteral("EFG/T2 tensor");
    case SignalValueShape::SphericalTensor:
        return QStringLiteral("spherical tensor");
    case SignalValueShape::TensorComponents:
        return QStringLiteral("tensor components");
    case SignalValueShape::PerClassBlock:
        return QStringLiteral("per-class block");
    case SignalValueShape::Embedding:
        return QStringLiteral("embedding");
    case SignalValueShape::RollupMoments:
        return QStringLiteral("rollup moments");
    case SignalValueShape::EventRecord:
        return QStringLiteral("event record");
    case SignalValueShape::Matrix:
        return QStringLiteral("matrix");
    case SignalValueShape::Spectrum:
        return QStringLiteral("spectrum");
    case SignalValueShape::CurveOverLag:
        return QStringLiteral("curve over lag");
    case SignalValueShape::Mat3PerRow:
        return QStringLiteral("Mat3 per row");
    case SignalValueShape::FixedFreqBlock:
        return QStringLiteral("fixed-frequency block");
    }
    return QStringLiteral("unknown");
}

QString ToString(SampleStatus status) {
    switch (status) {
    case SampleStatus::Valid:
        return QStringLiteral("valid");
    case SampleStatus::Gap:
        return QStringLiteral("gap");
    case SampleStatus::NotAvailable:
        return QStringLiteral("not available");
    case SampleStatus::Invalid:
        return QStringLiteral("invalid");
    }
    return QStringLiteral("unknown");
}

QString ToString(GapReason reason) {
    switch (reason) {
    case GapReason::None:
        return QStringLiteral("none");
    case GapReason::SourceAbsent:
        return QStringLiteral("source absent");
    case GapReason::FrameSourceAbsent:
        return QStringLiteral("frame source absent");
    case GapReason::SourceMaskOff:
        return QStringLiteral("source mask off");
    case GapReason::AnchorUnavailable:
        return QStringLiteral("anchor unavailable");
    case GapReason::NotApplicable:
        return QStringLiteral("not applicable");
    case GapReason::NaNSentinel:
        return QStringLiteral("NaN sentinel");
    case GapReason::MalformedSource:
        return QStringLiteral("malformed source");
    case GapReason::Pending:
        return QStringLiteral("pending");
    }
    return QStringLiteral("unknown");
}

QString ToString(UnitDimension dimension) {
    switch (dimension) {
    case UnitDimension::Dimensionless:
        return QStringLiteral("dimensionless");
    case UnitDimension::Length:
        return QStringLiteral("length");
    case UnitDimension::Angle:
        return QStringLiteral("angle");
    case UnitDimension::MagneticShielding:
        return QStringLiteral("magnetic shielding");
    case UnitDimension::ElectricField:
        return QStringLiteral("electric field");
    case UnitDimension::ElectricFieldGradient:
        return QStringLiteral("electric field gradient");
    case UnitDimension::Charge:
        return QStringLiteral("charge");
    case UnitDimension::Energy:
        return QStringLiteral("energy");
    case UnitDimension::Temperature:
        return QStringLiteral("temperature");
    case UnitDimension::Pressure:
        return QStringLiteral("pressure");
    case UnitDimension::Volume:
        return QStringLiteral("volume");
    case UnitDimension::Time:
        return QStringLiteral("time");
    case UnitDimension::Frequency:
        return QStringLiteral("frequency");
    case UnitDimension::Power:
        return QStringLiteral("power");
    case UnitDimension::Count:
        return QStringLiteral("count");
    case UnitDimension::Tag:
        return QStringLiteral("tag");
    }
    return QStringLiteral("unknown");
}

SignalAxis AxisForAnchor(const SignalAnchor& anchor) {
    return std::visit(
        Overloaded{
            [](const NoneAnchor&) { return SignalAxis::None; },
            [](const AtomAnchor&) { return SignalAxis::Atom; },
            [](const ResidueAnchor&) { return SignalAxis::Residue; },
            [](const AtomTupleAnchor&) { return SignalAxis::AtomTuple; },
            [](const BondAnchor&) { return SignalAxis::Bond; },
            [](const BondVectorAnchor&) { return SignalAxis::BondVector; },
            [](const RingAnchor&) { return SignalAxis::Ring; },
            [](const AromaticRingAnchor&) { return SignalAxis::AromaticRing; },
            [](const SaturatedRingAnchor&) { return SignalAxis::SaturatedRing; },
            [](const RingContributionPairAnchor&) { return SignalAxis::RingContributionPair; },
            [](const RingMembershipAnchor&) { return SignalAxis::RingMembership; },
            [](const MutationMatchPairAnchor&) { return SignalAxis::MutationMatchPair; },
            [](const ProteinAnchor&) { return SignalAxis::Protein; },
            [](const SystemAnchor&) { return SignalAxis::System; },
            [](const EventAnchor&) { return SignalAxis::Event; },
        },
        anchor);
}

bool AxisCanSatisfy(SignalAxis selectedAxis, SignalAxis requiredAxis) {
    if (requiredAxis == SignalAxis::None)
        return true;
    if (selectedAxis == requiredAxis)
        return true;
    if (requiredAxis == SignalAxis::Ring
        && (selectedAxis == SignalAxis::AromaticRing || selectedAxis == SignalAxis::SaturatedRing))
        return true;
    // Residue-grouping ergonomic: a BondVector descriptor accepts a parent
    // Residue anchor (picking a residue surfaces its bond vectors as a
    // sub-list, analogous to how Ring axis accepts Aromatic/Saturated).
    if (requiredAxis == SignalAxis::BondVector && selectedAxis == SignalAxis::Residue)
        return true;
    return false;
}

bool AnchorMatchesAxis(const SignalAnchor& anchor, SignalAxis axis) {
    return AxisCanSatisfy(AxisForAnchor(anchor), axis);
}

QString AnchorLabel(const SignalAnchor& anchor) {
    return std::visit(
        Overloaded{
            [](const NoneAnchor&) { return QStringLiteral("No anchor"); },
            [](const AtomAnchor& a) { return indexLabel("Atom", a.atom); },
            [](const ResidueAnchor& a) { return indexLabel("Residue", a.residue); },
            [](const AtomTupleAnchor& a) {
                QStringList parts;
                parts.reserve(static_cast<int>(a.atoms.size()));
                for (const std::size_t atom : a.atoms)
                    parts.push_back(QString::number(static_cast<qulonglong>(atom)));
                return QStringLiteral("Atoms [%1]").arg(parts.join(QStringLiteral(", ")));
            },
            [](const BondAnchor& a) { return indexLabel("Bond", a.bond); },
            [](const BondVectorAnchor& a) {
                const char* kindName = nullptr;
                switch (a.kind) {
                case 1: kindName = "N-H"; break;
                case 2: kindName = "Cα-Hα"; break;
                case 3: kindName = "C=O"; break;
                default: kindName = "vector"; break;
                }
                return QStringLiteral("Residue %1 %2")
                    .arg(static_cast<qulonglong>(a.residue))
                    .arg(QString::fromUtf8(kindName));
            },
            [](const RingAnchor& a) { return indexLabel("Ring", a.ring); },
            [](const AromaticRingAnchor& a) { return indexLabel("Aromatic ring", a.ring); },
            [](const SaturatedRingAnchor& a) { return indexLabel("Saturated ring", a.ring); },
            [](const RingContributionPairAnchor& a) { return indexLabel("Ring contribution pair", a.pair); },
            [](const RingMembershipAnchor& a) { return indexLabel("Ring membership", a.membership); },
            [](const MutationMatchPairAnchor& a) { return indexLabel("Mutation match pair", a.pair); },
            [](const ProteinAnchor&) { return QStringLiteral("Protein"); },
            [](const SystemAnchor&) { return QStringLiteral("System"); },
            [](const EventAnchor&) { return QStringLiteral("Event"); },
        },
        anchor);
}

QStringList AllDisplayModes(const SignalDescriptor& descriptor) {
    QStringList modes;
    modes.reserve(descriptor.temporalModes.size() + descriptor.staticModes.size());
    for (const QString& mode : descriptor.temporalModes) {
        if (!modes.contains(mode))
            modes.push_back(mode);
    }
    for (const QString& mode : descriptor.staticModes) {
        if (!modes.contains(mode))
            modes.push_back(mode);
    }
    return modes;
}

bool SupportsDisplayMode(const SignalDescriptor& descriptor, const QString& displayModeId) {
    if (displayModeId.isEmpty())
        return true;
    return descriptor.temporalModes.contains(displayModeId) || descriptor.staticModes.contains(displayModeId);
}

}  // namespace h5reader::model
