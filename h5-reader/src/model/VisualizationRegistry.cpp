#include "VisualizationRegistry.h"

#include "ChordCouplingVisualization.h"
#include "FixedFrequencyVisualization.h"
#include "LagCurveVisualization.h"
#include "NewmanVisualization.h"
#include "PowerSpectrumVisualization.h"
#include "SequenceBarVisualization.h"
#include "StripVisualization.h"
#include "TensorGlyphVisualization.h"
#include "TrajectorySignalCatalog.h"

#include <QLatin1StringView>
#include <QSet>

#include <array>

namespace h5reader::model {

namespace {

constexpr std::array<detail::DisplayModeCapabilityRow, 7> kStaticModeCapabilities{{
    {"static.bar.sequence", {true, true, true}},
    {"static.spectrum.power", {false, true, true}},
    {"static.curve.lag.animated", {true, true, true}},
    {"static.chord.coupling", {true, true, true}},
    {"static.fixed_freq", {true, true, true}},
    {"static.tensor", {false, false, true}},
    {"static.newman", {true, true, true}},
}};

bool isStripMode(const QString& modeId) {
    return modeId.startsWith(QStringLiteral("strip."));
}

}  // namespace

const VisualizationRegistry& VisualizationRegistry::instance() {
    static const VisualizationRegistry registry;
    return registry;
}

VisualizationRegistry::VisualizationRegistry() {
    defs_.push_back(std::make_unique<StripVisualization>());
    defs_.push_back(std::make_unique<SequenceBarVisualization>());
    defs_.push_back(std::make_unique<LagCurveVisualization>());
    defs_.push_back(std::make_unique<ChordCouplingVisualization>());
    defs_.push_back(std::make_unique<FixedFrequencyVisualization>());
    defs_.push_back(std::make_unique<PowerSpectrumVisualization>());
    defs_.push_back(std::make_unique<TensorGlyphVisualization>());
    defs_.push_back(std::make_unique<NewmanVisualization>());
}

QVector<const VisualizationDefinition*> VisualizationRegistry::definitions() const {
    QVector<const VisualizationDefinition*> out;
    out.reserve(static_cast<int>(defs_.size()));
    for (const auto& def : defs_)
        out.push_back(def.get());
    return out;
}

QVector<const VisualizationDefinition*>
VisualizationRegistry::supporting(const SignalDescriptor& descriptor) const {
    QVector<const VisualizationDefinition*> out;
    for (const auto& def : defs_) {
        if (def->supports(descriptor))
            out.push_back(def.get());
    }
    return out;
}

QVector<const VisualizationDefinition*>
VisualizationRegistry::visibleOfferable(const VisualizationContext& ctx,
                                        const SignalDescriptor& descriptor) const {
    QVector<const VisualizationDefinition*> out;
    for (const auto& def : defs_) {
        if (def->supports(descriptor)
            && def->isAvailable(ctx, descriptor)
            && def->capability().hasVisibleSurface) {
            out.push_back(def.get());
        }
    }
    return out;
}

QVector<const VisualizationDefinition*>
VisualizationRegistry::trackedButHidden(const SignalDescriptor& descriptor) const {
    QVector<const VisualizationDefinition*> out;
    for (const auto& def : defs_) {
        const DisplayModeCapability capability = def->capability();
        if (def->supports(descriptor)
            && capability.emitsPanelRef
            && !capability.hasVisibleSurface) {
            out.push_back(def.get());
        }
    }
    return out;
}

DisplayModeCapability VisualizationRegistry::capabilityForMode(const QString& modeId) const {
    if (isStripMode(modeId))
        return DisplayModeCapability{true, false, false};

    for (const detail::DisplayModeCapabilityRow& row : kStaticModeCapabilities) {
        if (modeId == QLatin1StringView(row.mode))
            return row.capability;
    }
    return {};
}

const VisualizationDefinition* VisualizationRegistry::definitionForMode(const QString& modeId) const {
    if (isStripMode(modeId))
        return nullptr;

    for (const auto& def : defs_) {
        if (def->legacyModeIds().contains(modeId))
            return def.get();
    }
    return nullptr;
}

QStringList VisualizationRegistry::unresolvedStaticModes(const TrajectorySignalCatalog& catalog) const {
    QSet<QString> unresolved;
    for (const SignalDescriptor& descriptor : catalog.allDescriptorList()) {
        for (const QString& mode : AllDisplayModes(descriptor)) {
            if (isStripMode(mode))
                continue;
            if (!definitionForMode(mode))
                unresolved.insert(mode);
        }
    }

    QStringList out;
    out.reserve(unresolved.size());
    for (const QString& mode : unresolved)
        out.push_back(mode);
    out.sort(Qt::CaseInsensitive);
    return out;
}

}  // namespace h5reader::model
