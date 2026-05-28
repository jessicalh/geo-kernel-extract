#include "GeometryStrips.h"

#include "AtomSelection.h"
#include "Conformation.h"
#include "ConformationGeometry.h"
#include "QtProtein.h"
#include "SpectralAnalysis.h"

#include <QStringList>

#include <algorithm>
#include <cmath>
#include <vector>

namespace h5reader::model {

namespace {
struct ValidSegment {
    std::size_t start = 0;
    std::size_t end = 0;
    std::size_t count = 0;
};

GeometryKind KindForName(const QString& name) {
    if (name == QStringLiteral("distance"))
        return GeometryKind::Distance;
    if (name == QStringLiteral("angle"))
        return GeometryKind::Angle;
    if (name == QStringLiteral("dihedral"))
        return GeometryKind::Dihedral;
    return GeometryKind::None;
}

QString NameForKind(GeometryKind kind) {
    switch (kind) {
    case GeometryKind::Distance:
        return QStringLiteral("distance");
    case GeometryKind::Angle:
        return QStringLiteral("angle");
    case GeometryKind::Dihedral:
        return QStringLiteral("dihedral");
    case GeometryKind::None:
        break;
    }
    return QString();
}

QString DisplayNameForKind(GeometryKind kind) {
    return QString::fromLatin1(NameForGeometryKind(kind));
}

SignalKey KeyForKind(GeometryKind kind) {
    return SignalKey{SignalFamily::Geometry, NameForKind(kind), QString()};
}

GeometryKind KindForTupleSize(std::size_t count) {
    switch (count) {
    case 2:
        return GeometryKind::Distance;
    case 3:
        return GeometryKind::Angle;
    case 4:
        return GeometryKind::Dihedral;
    default:
        return GeometryKind::None;
    }
}

QString UnitForKind(GeometryKind kind) {
    return kind == GeometryKind::Distance ? QString::fromUtf8("Å") : QString::fromUtf8("°");
}

int DecimalsForKind(GeometryKind kind) {
    return kind == GeometryKind::Distance ? 3 : 1;
}

bool IsGeometryTupleBinding(const SignalBinding& binding) {
    if (binding.key.family != SignalFamily::Geometry
        || binding.anchorKind != SignalAnchorKind::AtomTuple)
        return false;

    const GeometryKind tupleKind = KindForTupleSize(binding.atomTuple.size());
    if (tupleKind == GeometryKind::None)
        return false;

    return KindForName(binding.key.name) == tupleKind;
}

QString ResidueDisplayLabel(const QtProtein* protein, std::size_t residueIdx) {
    if (!protein || residueIdx >= protein->residueCount())
        return QStringLiteral("residue %1").arg(residueIdx + 1);
    const auto& res = protein->residue(residueIdx);
    const QString res3 = protein->residueLabel(residueIdx,
                                               NamingConvention::Amber,
                                               NamingSource::Verbatim);
    const QString chain = res.address.chainId.isEmpty()
                              ? QString()
                              : QStringLiteral("%1:").arg(res.address.chainId);
    return QStringLiteral("%1%2%3").arg(chain, res3).arg(res.address.residueNumber);
}

QString AtomDisplayLabel(const QtProtein* protein, std::size_t atomIdx) {
    if (!protein || atomIdx >= protein->atomCount())
        return QStringLiteral("#%1").arg(atomIdx);
    const auto& atom = protein->atom(atomIdx);
    const QString atomName = protein->atomNames(atomIdx).amber;
    if (atom.residueIndex >= 0
        && static_cast<std::size_t>(atom.residueIndex) < protein->residueCount()) {
        return QStringLiteral("%1:%2").arg(
            ResidueDisplayLabel(protein, static_cast<std::size_t>(atom.residueIndex)),
            atomName);
    }
    return QStringLiteral("#%1:%2").arg(atomIdx).arg(atomName);
}

QString TupleDisplayLabel(const Conformation* conformation, const std::vector<std::size_t>& atoms) {
    const QtProtein* protein = conformation ? conformation->protein() : nullptr;
    QStringList labels;
    labels.reserve(static_cast<int>(atoms.size()));
    for (std::size_t atom : atoms)
        labels.push_back(AtomDisplayLabel(protein, atom));
    return labels.join(QStringLiteral(" - "));
}

QString LabelForBinding(const Conformation* conformation, const SignalBinding& binding) {
    const GeometryKind kind = KindForTupleSize(binding.atomTuple.size());
    const QString tuple = TupleDisplayLabel(conformation, binding.atomTuple);
    if (tuple.isEmpty())
        return DisplayNameForKind(kind);
    return QStringLiteral("%1: %2").arg(DisplayNameForKind(kind), tuple);
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

void UnwrapDegrees(std::vector<double>& values) {
    for (std::size_t k = 1; k < values.size(); ++k) {
        double d = values[k] - values[k - 1];
        while (d > 180.0) {
            values[k] -= 360.0;
            d -= 360.0;
        }
        while (d < -180.0) {
            values[k] += 360.0;
            d += 360.0;
        }
    }
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

GeometryTupleTimeStrip::GeometryTupleTimeStrip(const Conformation* conformation)
    : conformation_(conformation)
{
    buffer_.id = QStringLiteral("geometry.time");
    buffer_.label = QStringLiteral("Geometry");
    buffer_.unit = QString::fromUtf8("°");
}

StripSpec GeometryTupleTimeStrip::spec() const {
    const GeometryKind kind = bound_ ? KindForTupleSize(binding_.atomTuple.size()) : GeometryKind::None;
    return StripSpec{
        QStringLiteral("geometry.tuple.time"),
        bound_ ? DisplayNameForKind(kind) : QStringLiteral("Geometry"),
        bound_ ? KeyForKind(kind) : SignalKey{SignalFamily::Geometry, QString(), QString()},
        SignalAnchorKind::AtomTuple,
        SignalDomainKind::FrameTime,
    };
}

bool GeometryTupleTimeStrip::canBind(const SignalBinding& binding) const {
    return conformation_ && IsGeometryTupleBinding(binding);
}

void GeometryTupleTimeStrip::bind(const SignalBinding& binding) {
    bound_ = canBind(binding);
    binding_ = binding;
    configureBuffer();
    buffer_.clear();
    emit bindingChanged();
}

void GeometryTupleTimeStrip::configureBuffer() {
    const GeometryKind kind = KindForTupleSize(binding_.atomTuple.size());
    buffer_.id = QStringLiteral("geometry.%1.time").arg(NameForKind(kind));
    buffer_.label = bound_ ? LabelForBinding(conformation_, binding_) : QStringLiteral("Geometry");
    buffer_.unit = UnitForKind(kind);
}

std::optional<double> GeometryTupleTimeStrip::sample(std::size_t frame) const {
    if (!bound_ || !conformation_ || frame >= conformation_->frameCount())
        return std::nullopt;
    const GeometryMeasurement measurement = Measure(*conformation_, frame, binding_.atomTuple);
    const GeometryKind kind = KindForTupleSize(binding_.atomTuple.size());
    if (!measurement.valid || measurement.kind != kind)
        return std::nullopt;
    return measurement.value;
}

void GeometryTupleTimeStrip::extendToFrame(std::size_t frame) {
    if (!bound_ || !conformation_ || conformation_->frameCount() == 0)
        return;
    const std::size_t last = std::min(frame, conformation_->frameCount() - 1);
    bool changed = false;
    for (long long f = buffer_.lastFrame() + 1; f <= static_cast<long long>(last); ++f) {
        buffer_.append(sample(static_cast<std::size_t>(f)));
        changed = true;
    }
    if (changed)
        emit bufferChanged(last);
}

StripRenderData GeometryTupleTimeStrip::renderData() const {
    const GeometryKind kind = KindForTupleSize(binding_.atomTuple.size());
    QString readout = QStringLiteral("—");
    for (std::size_t i = buffer_.values.size(); i > 0; --i) {
        const std::size_t idx = i - 1;
        if (idx < buffer_.valid.size() && buffer_.valid[idx]) {
            readout = QStringLiteral("f%1  %2 %3")
                          .arg(idx + 1)
                          .arg(buffer_.values[idx], 0, 'f', DecimalsForKind(kind))
                          .arg(buffer_.unit);
            break;
        }
    }
    return RenderBufferAsTimeData(buffer_, readout);
}

GeometryTupleFftStrip::GeometryTupleFftStrip(const Conformation* conformation)
    : conformation_(conformation)
{
    buffer_.id = QStringLiteral("geometry.fft.input");
    buffer_.label = QStringLiteral("Geometry");
    buffer_.unit = QString::fromUtf8("°");
}

StripSpec GeometryTupleFftStrip::spec() const {
    const GeometryKind kind = bound_ ? KindForTupleSize(binding_.atomTuple.size()) : GeometryKind::None;
    return StripSpec{
        QStringLiteral("geometry.tuple.fft"),
        bound_ ? QStringLiteral("%1 FFT").arg(DisplayNameForKind(kind)) : QStringLiteral("Geometry FFT"),
        bound_ ? KeyForKind(kind) : SignalKey{SignalFamily::Geometry, QString(), QString()},
        SignalAnchorKind::AtomTuple,
        SignalDomainKind::Frequency,
    };
}

bool GeometryTupleFftStrip::canBind(const SignalBinding& binding) const {
    return conformation_ && IsGeometryTupleBinding(binding);
}

void GeometryTupleFftStrip::bind(const SignalBinding& binding) {
    bound_ = canBind(binding);
    binding_ = binding;
    configureBuffer();
    buffer_.clear();
    spectrum_.clear();
    readout_ = QStringLiteral("—");
    emit bindingChanged();
}

void GeometryTupleFftStrip::configureBuffer() {
    const GeometryKind kind = KindForTupleSize(binding_.atomTuple.size());
    buffer_.id = QStringLiteral("geometry.%1.fft.input").arg(NameForKind(kind));
    buffer_.label = bound_ ? LabelForBinding(conformation_, binding_) : QStringLiteral("Geometry");
    buffer_.unit = UnitForKind(kind);
}

std::optional<double> GeometryTupleFftStrip::sample(std::size_t frame) const {
    if (!bound_ || !conformation_ || frame >= conformation_->frameCount())
        return std::nullopt;
    const GeometryMeasurement measurement = Measure(*conformation_, frame, binding_.atomTuple);
    const GeometryKind kind = KindForTupleSize(binding_.atomTuple.size());
    if (!measurement.valid || measurement.kind != kind)
        return std::nullopt;
    return measurement.value;
}

void GeometryTupleFftStrip::extendToFrame(std::size_t frame) {
    if (!bound_ || !conformation_ || conformation_->frameCount() == 0)
        return;
    const std::size_t last = std::min(frame, conformation_->frameCount() - 1);
    bool changed = false;
    for (long long f = buffer_.lastFrame() + 1; f <= static_cast<long long>(last); ++f) {
        buffer_.append(sample(static_cast<std::size_t>(f)));
        changed = true;
    }
    if (changed)
        rebuildSpectrum();
    if (changed)
        emit bufferChanged(last);
}

void GeometryTupleFftStrip::rebuildSpectrum() {
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

    const GeometryKind kind = KindForTupleSize(binding_.atomTuple.size());
    if (kind == GeometryKind::Angle || kind == GeometryKind::Dihedral)
        UnwrapDegrees(values);

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

StripRenderData GeometryTupleFftStrip::renderData() const {
    StripRenderData data;
    data.label = QStringLiteral("Power spectrum");
    data.xUnit = QStringLiteral("1/ns");
    data.yUnit = QStringLiteral("power");
    data.readout = readout_;
    data.domainKind = SignalDomainKind::Frequency;
    data.points = spectrum_;
    return data;
}

}  // namespace h5reader::model
