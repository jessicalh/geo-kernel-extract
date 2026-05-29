// SignalDictionary -- legacy reveal vocabulary for dashboard strip requests.
//
// Active signal identity lives in DashboardSignal/DisplaySignalBinding. This
// header remains only as a small adapter for the current strip renderer and VTK
// reveal overlay while those signals are fanned through Qt widgets.

#pragma once

#include "DashboardSignal.h"

#include <QString>

#include <cstdint>

namespace h5reader::model {

enum class SignalAnchorKind : uint8_t {
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

struct SignalBinding {
    QString descriptorId;
    QString conceptKey;
    SignalAnchor anchor = NoneAnchor{};
    bool followsFocus = false;

    friend bool operator==(const SignalBinding& a, const SignalBinding& b)
    {
        return a.descriptorId == b.descriptorId
               && a.conceptKey == b.conceptKey
               && a.anchor == b.anchor
               && a.followsFocus == b.followsFocus;
    }
};

}  // namespace h5reader::model
