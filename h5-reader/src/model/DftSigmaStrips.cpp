#include "DftSigmaStrips.h"

#include "Conformation.h"
#include "SpectralAnalysis.h"

#include <algorithm>
#include <cmath>

namespace h5reader::model {

namespace {
struct ValidSegment {
    std::size_t start = 0;
    std::size_t end = 0;
    std::size_t count = 0;
};

SignalKey DftSigmaKey() {
    return SignalKey{SignalFamily::DftShielding, QStringLiteral("sigma"), QStringLiteral("total.T0")};
}

bool IsDftSigmaAtomBinding(const SignalBinding& binding) {
    return binding.key == DftSigmaKey()
        && binding.anchorKind == SignalAnchorKind::Atom
        && binding.atom.has_value();
}

double DtPicosecondsForValidFrames(const Conformation* conformation, const std::vector<std::size_t>& frames) {
    if (!conformation || frames.size() < 2)
        return 0.0;
    const double dt = conformation->timePicoseconds(frames[1]) - conformation->timePicoseconds(frames[0]);
    if (dt <= 0.0 || !std::isfinite(dt))
        return 0.0;
    for (std::size_t i = 1; i < frames.size(); ++i) {
        if (frames[i] != frames[i - 1] + 1)
            return 0.0;
        const double step = conformation->timePicoseconds(frames[i]) - conformation->timePicoseconds(frames[i - 1]);
        if (std::abs(step - dt) > std::max(1e-9, std::abs(dt) * 1e-6))
            return 0.0;
    }
    return dt;
}

ValidSegment LongestContiguousValidSegment(const ChannelBuffer& buffer) {
    ValidSegment best;
    ValidSegment current;
    bool inRun = false;

    for (std::size_t frame = 0; frame < buffer.values.size(); ++frame) {
        const bool valid = frame < buffer.valid.size() && buffer.valid[frame];
        if (!valid) {
            if (inRun && current.count > best.count)
                best = current;
            inRun = false;
            current = {};
            continue;
        }
        if (!inRun) {
            current.start = frame;
            current.count = 0;
            inRun = true;
        }
        current.end = frame;
        ++current.count;
    }

    if (inRun && current.count > best.count)
        best = current;
    return best;
}

StripRenderData RenderBufferAsTimeData(const ChannelBuffer& buffer, const QString& readout) {
    StripRenderData data;
    data.label = buffer.label;
    data.xUnit = QStringLiteral("frame");
    data.yUnit = buffer.unit;
    data.readout = readout;
    data.domainKind = SignalDomainKind::FrameTime;
    data.points.reserve(static_cast<int>(buffer.size()));
    for (std::size_t frame = 0; frame < buffer.values.size(); ++frame) {
        if (frame >= buffer.valid.size() || !buffer.valid[frame])
            continue;
        data.points.push_back(QPointF(static_cast<double>(frame), buffer.values[frame]));
    }
    return data;
}
}  // namespace

DftSigmaAtomTimeStrip::DftSigmaAtomTimeStrip(const Conformation* conformation, DftShieldingStore* store)
    : conformation_(conformation),
      store_(store)
{
    buffer_.id = QStringLiteral("dft.sigma.total.T0.time");
    buffer_.label = QStringLiteral("DFT sigma total");
    buffer_.unit = QStringLiteral("ppm");
}

StripSpec DftSigmaAtomTimeStrip::spec() const {
    return StripSpec{
        QStringLiteral("dft.sigma.total.T0.time"),
        QStringLiteral("DFT sigma total"),
        DftSigmaKey(),
        SignalAnchorKind::Atom,
        SignalDomainKind::FrameTime,
    };
}

bool DftSigmaAtomTimeStrip::canBind(const SignalBinding& binding) const {
    return conformation_ && store_ && IsDftSigmaAtomBinding(binding);
}

void DftSigmaAtomTimeStrip::bind(const SignalBinding& binding) {
    bound_ = canBind(binding);
    binding_ = binding;
    buffer_.clear();
}

std::optional<double> DftSigmaAtomTimeStrip::sample(std::size_t frame) const {
    if (!bound_ || !conformation_ || !store_ || !binding_.atom)
        return std::nullopt;
    if (frame >= conformation_->frameCount())
        return std::nullopt;
    const std::size_t original = conformation_->originalFrameIndex(frame);
    store_->requestFrame(original);
    return store_->sample(original, *binding_.atom, DftPart::Total, DftScalar::IsotropicT0);
}

void DftSigmaAtomTimeStrip::extendToFrame(std::size_t frame) {
    if (!bound_)
        return;
    const std::size_t last = conformation_ ? std::min(frame, conformation_->frameCount() - 1) : frame;
    const long long from = buffer_.lastFrame() + 1;
    for (long long f = from; f < static_cast<long long>(last); ++f)
        buffer_.append(std::nullopt);
    if (from <= static_cast<long long>(last))
        buffer_.append(sample(last));
}

StripRenderData DftSigmaAtomTimeStrip::renderData() const {
    QString readout = QStringLiteral("—");
    for (std::size_t i = buffer_.values.size(); i > 0; --i) {
        const std::size_t idx = i - 1;
        if (idx < buffer_.valid.size() && buffer_.valid[idx]) {
            readout = QStringLiteral("f%1  %2 ppm").arg(idx + 1).arg(buffer_.values[idx], 0, 'f', 2);
            break;
        }
    }
    return RenderBufferAsTimeData(buffer_, readout);
}

DftSigmaAtomFftStrip::DftSigmaAtomFftStrip(const Conformation* conformation, DftShieldingStore* store)
    : conformation_(conformation),
      store_(store)
{
    buffer_.id = QStringLiteral("dft.sigma.total.T0.fft.input");
    buffer_.label = QStringLiteral("DFT sigma total");
    buffer_.unit = QStringLiteral("ppm");
}

StripSpec DftSigmaAtomFftStrip::spec() const {
    return StripSpec{
        QStringLiteral("dft.sigma.total.T0.fft"),
        QStringLiteral("DFT sigma total FFT"),
        DftSigmaKey(),
        SignalAnchorKind::Atom,
        SignalDomainKind::Frequency,
    };
}

bool DftSigmaAtomFftStrip::canBind(const SignalBinding& binding) const {
    return conformation_ && store_ && IsDftSigmaAtomBinding(binding);
}

void DftSigmaAtomFftStrip::bind(const SignalBinding& binding) {
    bound_ = canBind(binding);
    binding_ = binding;
    buffer_.clear();
    spectrum_.clear();
    readout_ = QStringLiteral("—");
}

std::optional<double> DftSigmaAtomFftStrip::sample(std::size_t frame) const {
    if (!bound_ || !conformation_ || !store_ || !binding_.atom)
        return std::nullopt;
    if (frame >= conformation_->frameCount())
        return std::nullopt;
    const std::size_t original = conformation_->originalFrameIndex(frame);
    store_->requestFrame(original);
    return store_->sample(original, *binding_.atom, DftPart::Total, DftScalar::IsotropicT0);
}

void DftSigmaAtomFftStrip::extendToFrame(std::size_t frame) {
    if (!bound_)
        return;
    const std::size_t last = conformation_ ? std::min(frame, conformation_->frameCount() - 1) : frame;
    const long long from = buffer_.lastFrame() + 1;
    bool changed = false;
    for (long long f = from; f < static_cast<long long>(last); ++f) {
        buffer_.append(std::nullopt);
        changed = true;
    }
    if (from <= static_cast<long long>(last)) {
        buffer_.append(sample(last));
        changed = true;
    }
    if (changed)
        rebuildSpectrum();
}

void DftSigmaAtomFftStrip::rebuildSpectrum() {
    const ValidSegment segment = LongestContiguousValidSegment(buffer_);
    std::vector<double> values;
    std::vector<std::size_t> frames;
    values.reserve(segment.count);
    frames.reserve(segment.count);
    if (segment.count >= 2) {
        for (std::size_t frame = segment.start; frame <= segment.end; ++frame) {
            if (frame >= buffer_.valid.size() || !buffer_.valid[frame])
                break;
            values.push_back(buffer_.values[frame]);
            frames.push_back(frame);
        }
    }

    if (values.size() < 2) {
        spectrum_.clear();
        readout_ = segment.count == 1
                       ? QStringLiteral("single valid frame")
                       : QStringLiteral("—");
        return;
    }

    if (values.size() != segment.count) {
        spectrum_.clear();
        readout_ = QStringLiteral("gapped samples");
        return;
    }

    const double dtPs = DtPicosecondsForValidFrames(conformation_, frames);
    if (dtPs <= 0.0) {
        spectrum_.clear();
        readout_ = QStringLiteral("non-uniform samples");
        return;
    }

    const PowerSpectrum spec = ComputePowerSpectrum(values, dtPs);
    spectrum_.clear();
    if (!spec.valid || spec.power.size() < 2) {
        readout_ = QStringLiteral("—");
        return;
    }

    spectrum_.reserve(static_cast<int>(spec.power.size()));
    for (std::size_t k = 1; k < spec.power.size(); ++k)
        spectrum_.push_back(QPointF(spec.frequencyPerNs[k], spec.power[k]));

    if (spec.dominantPeriodPs <= 0.0) {
        readout_ = QStringLiteral("f%1-f%2  no dominant period").arg(segment.start + 1).arg(segment.end + 1);
    } else {
        const QString period = spec.dominantPeriodPs >= 1000.0
                                   ? QStringLiteral("%1 ns").arg(spec.dominantPeriodPs / 1000.0, 0, 'f', 2)
                                   : QStringLiteral("%1 ps").arg(spec.dominantPeriodPs, 0, 'f', 0);
        readout_ = QStringLiteral("f%1-f%2  %3").arg(segment.start + 1).arg(segment.end + 1).arg(period);
    }
}

StripRenderData DftSigmaAtomFftStrip::renderData() const {
    StripRenderData data;
    data.label = QStringLiteral("DFT sigma total FFT");
    data.xUnit = QStringLiteral("1/ns");
    data.yUnit = QStringLiteral("power");
    data.readout = readout_;
    data.domainKind = SignalDomainKind::Frequency;
    data.points = spectrum_;
    return data;
}

}  // namespace h5reader::model
