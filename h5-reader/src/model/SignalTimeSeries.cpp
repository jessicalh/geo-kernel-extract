#include "SignalTimeSeries.h"

#include <optional>

namespace h5reader::model {

QString SignalChannelKey::stableId() const {
    return QStringLiteral("%1|%2|%3|%4")
        .arg(signalId.toString(QUuid::WithoutBraces), descriptorId, displayModeId, channelId);
}

FrameSignalSample FrameSignalSample::Valid(double x) {
    FrameSignalSample sample;
    sample.status = SampleStatus::Valid;
    sample.gapReason = GapReason::None;
    sample.value = x;
    return sample;
}

FrameSignalSample FrameSignalSample::Gap(GapReason reason) {
    FrameSignalSample sample;
    sample.status = SampleStatus::Gap;
    sample.gapReason = reason;
    return sample;
}

void SignalBuffer::clear() {
    channel.clear();
    statuses.clear();
    gapReasons.clear();
}

void SignalBuffer::append(FrameSignalSample sample) {
    statuses.push_back(sample.status);
    gapReasons.push_back(sample.gapReason);
    if (sample.status == SampleStatus::Valid)
        channel.append(sample.value);
    else
        channel.append(std::nullopt);
}

}  // namespace h5reader::model
