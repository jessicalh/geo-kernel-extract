#pragma once

#include <QLatin1StringView>
#include <QString>

#include <array>

namespace h5reader::model {

struct DisplayModeCapability {
    bool hasVisibleSurface = false;  // dialog: offer as a pickable display mode
    bool buildsPanelWidget = false;  // controller: build an AbstractStripPanel via setOwnedPanels
    bool emitsPanelRef = false;      // panel-model: track with a "panel" sentinel display ref
};

namespace detail {

struct DisplayModeCapabilityRow {
    const char* mode = "";
    DisplayModeCapability capability;
};

}  // namespace detail

inline DisplayModeCapability DisplayModeCapabilityFor(const QString& mode) {
    if (mode.startsWith(QStringLiteral("strip.")))
        return DisplayModeCapability{true, false, false};

    static constexpr std::array<detail::DisplayModeCapabilityRow, 6> kDisplayModeCapabilities{{
        { "static.bar.sequence", { true, true, true } },
        // PowerSpectrumPanel is renderable and tracked by a panel ref, but the
        // picker-visible gate intentionally stays closed for a later decision.
        { "static.spectrum.power", { false, true, true } },
        { "static.curve.lag.animated", { true, true, true } },
        { "static.chord.coupling", { true, true, true } },
        { "static.fixed_freq", { true, true, true } },
        // static.tensor is a scene-overlay mode, not a strip widget panel, but
        // it still needs a "panel" ref so ownership survives the deferred
        // trigger gesture.
        { "static.tensor", { false, false, true } },
    }};

    for (const detail::DisplayModeCapabilityRow& row : kDisplayModeCapabilities) {
        if (mode == QLatin1StringView(row.mode))
            return row.capability;
    }
    return {};
}

}  // namespace h5reader::model
