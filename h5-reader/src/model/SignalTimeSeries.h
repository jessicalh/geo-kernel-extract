#pragma once

#include "DashboardSignal.h"
#include "StripChartChannel.h"

#include <QUuid>
#include <QString>

#include <cstddef>
#include <limits>
#include <vector>

namespace h5reader::model {

struct SignalChannelKey {
    QUuid signalId;
    QString descriptorId;
    QString displayModeId;
    QString channelId;

    QString stableId() const;
};

struct FrameSignalSample {
    SampleStatus status = SampleStatus::Gap;
    GapReason gapReason = GapReason::Pending;
    double value = std::numeric_limits<double>::quiet_NaN();

    static FrameSignalSample Valid(double x);
    static FrameSignalSample Gap(GapReason reason);
};

struct SignalBuffer {
    SignalChannelKey key;
    ChannelBuffer channel;
    std::vector<SampleStatus> statuses;
    std::vector<GapReason> gapReasons;

    long long lastFrame() const { return channel.lastFrame(); }
    void clear();
    void append(FrameSignalSample sample);
};

}  // namespace h5reader::model
