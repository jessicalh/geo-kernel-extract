#pragma once

#include <QString>
#include <QStringList>
#include <QUuid>
#include <QVariant>
#include <QVector>

#include <cstddef>
#include <cstdint>
#include <variant>
#include <vector>

namespace h5reader::model {

enum class SignalSourceKind : std::uint8_t {
    DenseH5Trajectory = 0,
    FrameNpySnapshot,
    OrcaDftFrame,
    Topology,
    DerivedGeometry,
    SelectionEvents,
};

enum class SourceResidency : std::uint8_t {
    StartupLoaded = 0,
    FrameLoaded,
    Derived,
};

enum class SignalAxis : std::uint8_t {
    None = 0,
    Atom,
    Residue,
    AtomTuple,
    Bond,
    BondVector,
    Ring,
    AromaticRing,
    SaturatedRing,
    RingContributionPair,
    RingMembership,
    MutationMatchPair,
    Protein,
    System,
    Event,
};

enum class SignalValueShape : std::uint8_t {
    Scalar = 0,
    Count,
    Category,
    Vector3,
    EfgT2,
    SphericalTensor,
    TensorComponents,
    PerClassBlock,
    Embedding,
    RollupMoments,
    EventRecord,
    // Per-row N×N (KernelCoherence). Layout (row, channel, channel).
    Matrix,
    // Per-row 1D curve over a frequency grid (Spectrum). Distinct from
    // CurveOverLag so the catalog/REST inventory and the display-mode
    // dispatch can tell power-spectrum panels from lag-decay panels.
    Spectrum,
    // Per-row 1D curve over a lag grid (KernelDynamics ACF,
    // ReorientationalDynamics body/lab TCFs, DihedralAutocorrelation).
    CurveOverLag,
    // Per-row 3×3 dense tensor (ReorientationalDynamics body-frame
    // <u⊗u>). TensorComponents is reserved for the rank-2 EFG/CSA family
    // already in use; this shape is for a bare Mat3 attached to a
    // bond-vector row.
    Mat3PerRow,
    // Per-row K-vector at K externally-fixed frequencies
    // (ReorientationalDynamics J(ω) sampled at the 5 KTB Larmor
    // combinations). The K frequencies are descriptor metadata, not a
    // per-row axis.
    FixedFreqBlock,
};

enum class SampleStatus : std::uint8_t {
    Valid = 0,
    Gap,
    NotAvailable,
    Invalid,
};

enum class GapReason : std::uint8_t {
    None = 0,
    SourceAbsent,
    FrameSourceAbsent,
    SourceMaskOff,
    AnchorUnavailable,
    NotApplicable,
    NaNSentinel,
    MalformedSource,
    Pending,
};

enum class UnitDimension : std::uint8_t {
    Dimensionless = 0,
    Length,
    Angle,
    MagneticShielding,
    ElectricField,
    ElectricFieldGradient,
    Charge,
    Energy,
    Temperature,
    Pressure,
    Volume,
    Time,
    Frequency,
    Power,
    Count,
    Tag,
};

struct UnitSpec {
    UnitDimension dimension = UnitDimension::Dimensionless;
    QString sourceSymbol;
    QString displaySymbol;
    double scaleToDisplay = 1.0;
    double offsetToDisplay = 0.0;
    bool convertible = true;
};

struct ChannelDescriptor {
    QString id;
    QString label;
    SignalValueShape valueShape = SignalValueShape::Scalar;
    UnitSpec sourceUnits;
    UnitSpec defaultDisplayUnits;
};

struct SignalDescriptor {
    QString id;
    QString conceptKey;
    SignalSourceKind sourceKind = SignalSourceKind::DenseH5Trajectory;
    QString importSet;
    QString family;
    QString label;
    QString description;
    QString storagePath;
    UnitSpec sourceUnits;
    UnitSpec defaultDisplayUnits;

    SourceResidency residency = SourceResidency::StartupLoaded;
    SignalAxis nativeAxis = SignalAxis::None;
    SignalAxis requiredAnchor = SignalAxis::None;
    SignalValueShape valueShape = SignalValueShape::Scalar;

    bool temporal = false;
    bool staticDisplay = false;
    bool sourceAttachedMask = false;
    bool frameLocalPayload = false;
    bool finiteScalarRequired = true;

    SampleStatus samplingStatus = SampleStatus::Gap;
    GapReason samplingGapReason = GapReason::Pending;

    QStringList temporalModes;
    QStringList staticModes;
    QVector<ChannelDescriptor> channels;
    QStringList tags;
};

struct NoneAnchor {
    friend bool operator==(const NoneAnchor&, const NoneAnchor&) { return true; }
};
struct AtomAnchor {
    std::size_t atom = 0;
    friend bool operator==(const AtomAnchor& a, const AtomAnchor& b) { return a.atom == b.atom; }
};
struct ResidueAnchor {
    std::size_t residue = 0;
    friend bool operator==(const ResidueAnchor& a, const ResidueAnchor& b) { return a.residue == b.residue; }
};
struct AtomTupleAnchor {
    std::vector<std::size_t> atoms;
    friend bool operator==(const AtomTupleAnchor& a, const AtomTupleAnchor& b) { return a.atoms == b.atoms; }
};
struct BondAnchor {
    std::size_t bond = 0;
    friend bool operator==(const BondAnchor& a, const BondAnchor& b) { return a.bond == b.bond; }
};
// Semantic bond-vector identity (residue, kind). The producer-side
// IRedOrderParameter emits only N-H vectors; ReorientationalDynamics
// emits N-H / Cα-Hα / C=O. Each TR has its own M-row identity table,
// so an opaque row-index would not be invariant across TRs. The (residue,
// kind) tuple IS invariant: the per-TR sampler resolves it to its
// table's row at lookup time. Kind values match the producer's
// vector_kind enum: 1=NH, 2=CaHa, 3=CO; 0 = any/unspecified.
struct BondVectorAnchor {
    std::size_t residue = 0;
    std::uint8_t kind = 0;
    friend bool operator==(const BondVectorAnchor& a, const BondVectorAnchor& b) {
        return a.residue == b.residue && a.kind == b.kind;
    }
};
struct RingAnchor {
    std::size_t ring = 0;
    friend bool operator==(const RingAnchor& a, const RingAnchor& b) { return a.ring == b.ring; }
};
struct AromaticRingAnchor {
    std::size_t ring = 0;
    friend bool operator==(const AromaticRingAnchor& a, const AromaticRingAnchor& b) { return a.ring == b.ring; }
};
struct SaturatedRingAnchor {
    std::size_t ring = 0;
    friend bool operator==(const SaturatedRingAnchor& a, const SaturatedRingAnchor& b) { return a.ring == b.ring; }
};
struct RingContributionPairAnchor {
    std::size_t pair = 0;
    friend bool operator==(const RingContributionPairAnchor& a, const RingContributionPairAnchor& b)
    {
        return a.pair == b.pair;
    }
};
struct RingMembershipAnchor {
    std::size_t membership = 0;
    friend bool operator==(const RingMembershipAnchor& a, const RingMembershipAnchor& b)
    {
        return a.membership == b.membership;
    }
};
struct MutationMatchPairAnchor {
    std::size_t pair = 0;
    friend bool operator==(const MutationMatchPairAnchor& a, const MutationMatchPairAnchor& b)
    {
        return a.pair == b.pair;
    }
};
struct ProteinAnchor {
    friend bool operator==(const ProteinAnchor&, const ProteinAnchor&) { return true; }
};
struct SystemAnchor {
    friend bool operator==(const SystemAnchor&, const SystemAnchor&) { return true; }
};
struct EventAnchor {
    friend bool operator==(const EventAnchor&, const EventAnchor&) { return true; }
};

using SignalAnchor = std::variant<NoneAnchor,
                                  AtomAnchor,
                                  ResidueAnchor,
                                  AtomTupleAnchor,
                                  BondAnchor,
                                  BondVectorAnchor,
                                  RingAnchor,
                                  AromaticRingAnchor,
                                  SaturatedRingAnchor,
                                  RingContributionPairAnchor,
                                  RingMembershipAnchor,
                                  MutationMatchPairAnchor,
                                  ProteinAnchor,
                                  SystemAnchor,
                                  EventAnchor>;

struct DisplaySignalBinding {
    SignalSourceKind sourceKind = SignalSourceKind::DenseH5Trajectory;
    QString descriptorId;
    QString conceptKey;
    QString reducerId;
    QString displayModeId;
    SignalAnchor anchor = NoneAnchor{};
    bool followsFocus = false;
};

struct DashboardSignal {
    QUuid id;
    DisplaySignalBinding binding;
    QString label;
    QStringList displayModeIds;
    bool enabled = true;

    SignalAxis nativeAxis = SignalAxis::None;
    SignalAxis requiredAnchor = SignalAxis::None;
    SignalValueShape valueShape = SignalValueShape::Scalar;
};

struct SignalSample {
    SampleStatus status = SampleStatus::Gap;
    GapReason gapReason = GapReason::None;
    SignalValueShape shape = SignalValueShape::Scalar;
    QVariant value;
    QString diagnostic;
};

QString ToString(SignalSourceKind kind);
QString ToString(SourceResidency residency);
QString ToString(SignalAxis axis);
QString ToString(SignalValueShape shape);
QString ToString(SampleStatus status);
QString ToString(GapReason reason);
QString ToString(UnitDimension dimension);

SignalAxis AxisForAnchor(const SignalAnchor& anchor);
// True when an anchor whose axis is `selectedAxis` can satisfy a
// descriptor that requires `requiredAxis`. Encodes the widening rules
// (Aromatic/Saturated → Ring, Residue → BondVector). Used by both the
// controller-side anchor variant check (AnchorMatchesAxis below) and
// the picker dialog's filter (SignalDisplayDialog), so the two stay
// in step.
bool AxisCanSatisfy(SignalAxis selectedAxis, SignalAxis requiredAxis);
bool AnchorMatchesAxis(const SignalAnchor& anchor, SignalAxis axis);
QString AnchorLabel(const SignalAnchor& anchor);
QStringList AllDisplayModes(const SignalDescriptor& descriptor);
bool SupportsDisplayMode(const SignalDescriptor& descriptor, const QString& displayModeId);

}  // namespace h5reader::model

Q_DECLARE_METATYPE(h5reader::model::SignalSourceKind)
Q_DECLARE_METATYPE(h5reader::model::SourceResidency)
Q_DECLARE_METATYPE(h5reader::model::SignalAxis)
Q_DECLARE_METATYPE(h5reader::model::SignalValueShape)
Q_DECLARE_METATYPE(h5reader::model::SampleStatus)
Q_DECLARE_METATYPE(h5reader::model::GapReason)
Q_DECLARE_METATYPE(h5reader::model::UnitDimension)
