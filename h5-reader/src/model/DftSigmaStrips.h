// DftSigmaStrips -- concrete dashboard strips for atom DFT sigma.
//
// These are deliberately concrete. The time strip and FFT strip bind to the
// same SignalKey/atom anchor and both sample DftShieldingStore directly; the
// FFT strip does not depend on the time strip being present.

#pragma once

#include "AbstractStrip.h"
#include "DftShieldingStore.h"
#include "StripChartChannel.h"

#include <QPointer>

#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::model {

class Conformation;

class DftSigmaAtomTimeStrip final : public AbstractStrip {
public:
    DftSigmaAtomTimeStrip(const Conformation* conformation, DftShieldingStore* store);

    StripSpec spec() const override;
    bool canBind(const SignalBinding& binding) const override;
    void bind(const SignalBinding& binding) override;
    void extendToFrame(std::size_t frame) override;
    StripRenderData renderData() const override;

    const ChannelBuffer& buffer() const { return buffer_; }

private:
    std::optional<double> sample(std::size_t frame) const;

    const Conformation* conformation_ = nullptr;
    QPointer<DftShieldingStore> store_;
    SignalBinding binding_;
    ChannelBuffer buffer_;
    bool bound_ = false;
};

class DftSigmaAtomFftStrip final : public AbstractStrip {
public:
    DftSigmaAtomFftStrip(const Conformation* conformation, DftShieldingStore* store);

    StripSpec spec() const override;
    bool canBind(const SignalBinding& binding) const override;
    void bind(const SignalBinding& binding) override;
    void extendToFrame(std::size_t frame) override;
    StripRenderData renderData() const override;

private:
    std::optional<double> sample(std::size_t frame) const;
    void rebuildSpectrum();

    const Conformation* conformation_ = nullptr;
    QPointer<DftShieldingStore> store_;
    SignalBinding binding_;
    ChannelBuffer buffer_;
    QVector<QPointF> spectrum_;
    QString readout_ = QStringLiteral("—");
    bool bound_ = false;
};

}  // namespace h5reader::model
