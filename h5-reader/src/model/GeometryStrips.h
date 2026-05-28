// GeometryStrips -- concrete dashboard strips for atom-tuple geometry.
//
// Distance, angle, and dihedral strips deliberately use the existing
// ConformationGeometry::Measure path. The time strip and FFT strip bind to the
// same Geometry/AtomTuple signal and sample independently, matching the DFT
// strip convention.

#pragma once

#include "StripCalculation.h"
#include "StripChartChannel.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class Conformation;

class GeometryTupleTimeStrip final : public StripCalculation {
public:
    explicit GeometryTupleTimeStrip(const Conformation* conformation);

    StripSpec spec() const override;
    bool canBind(const SignalBinding& binding) const override;
    void bind(const SignalBinding& binding) override;
    void extendToFrame(std::size_t frame) override;
    StripRenderData renderData() const override;

    const ChannelBuffer& buffer() const { return buffer_; }

private:
    std::optional<double> sample(std::size_t frame) const;
    void configureBuffer();

    const Conformation* conformation_ = nullptr;
    SignalBinding binding_;
    ChannelBuffer buffer_;
    bool bound_ = false;
};

class GeometryTupleFftStrip final : public StripCalculation {
public:
    explicit GeometryTupleFftStrip(const Conformation* conformation);

    StripSpec spec() const override;
    bool canBind(const SignalBinding& binding) const override;
    void bind(const SignalBinding& binding) override;
    void extendToFrame(std::size_t frame) override;
    StripRenderData renderData() const override;

private:
    std::optional<double> sample(std::size_t frame) const;
    void configureBuffer();
    void rebuildSpectrum();

    const Conformation* conformation_ = nullptr;
    SignalBinding binding_;
    ChannelBuffer buffer_;
    QVector<QPointF> spectrum_;
    QString readout_ = QStringLiteral("—");
    bool bound_ = false;
};

}  // namespace h5reader::model
