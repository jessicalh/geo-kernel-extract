#pragma once

#include "DashboardSignal.h"
#include "DisplayModeCapability.h"

#include <QString>
#include <QStringList>
#include <QVector>

#include <cstdint>

namespace h5reader::model {

class TrajectoryFieldAvailability;

enum class VisualizationType : std::uint8_t {
    TemporalStrip,
    TensorGlyph,
    AtomColor,
    SequenceBar,
    LagCurve,
    ChordCoupling,
    FixedFrequency,
    PowerSpectrum,
};

enum class DisplaySurface : std::uint8_t {
    Strip,
    Panel,
    Scene,
};

enum class StripComponent : std::uint8_t {
    Auto,
    VectorX,
    VectorY,
    VectorZ,
    VectorMagnitude,
    TensorT0,
    TensorT1,
    TensorT2,
    TensorComponent,
};

struct VisualizationContext {
    const TrajectoryFieldAvailability* availability = nullptr;
    bool hasTrajectory = false;
    bool hasDftStore = false;
    bool hasSceneOverlay = false;
    bool tensorGlyphGestureEnabled = false;
};

struct StripSeriesRequest {
    QString descriptorId;
    ChannelDescriptor channel;
    StripComponent component = StripComponent::Auto;
};

class VisualizationDefinition {
public:
    virtual ~VisualizationDefinition() = default;

    virtual VisualizationType type() const = 0;
    virtual QString label() const = 0;
    virtual DisplaySurface surface() const = 0;

    virtual bool supports(const SignalDescriptor& descriptor) const = 0;
    virtual bool isAvailable(const VisualizationContext& ctx,
                             const SignalDescriptor& descriptor) const = 0;
    virtual DisplayModeCapability capability() const = 0;
    virtual QStringList legacyModeIds() const = 0;
};

QString ToString(VisualizationType type);
QString ToString(DisplaySurface surface);

bool VisualizationDescriptorOffersMode(const SignalDescriptor& descriptor,
                                       const QString& modeId);
bool VisualizationDescriptorDataAvailable(const VisualizationContext& ctx,
                                          const SignalDescriptor& descriptor);

}  // namespace h5reader::model
