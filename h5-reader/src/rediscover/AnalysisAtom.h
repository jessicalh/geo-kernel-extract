// AnalysisAtom -- atom branch of the analysis-object model.

#pragma once

#include "AnalysisElement.h"

#include <QJsonObject>

#include <cstddef>
#include <memory>

namespace h5reader::rediscover {

class AnalysisAtom final : public AnalysisElement {
public:
    AnalysisAtom(const AnalysisObjectContext& context,
                 std::size_t atomIndex,
                 PerAtomSubstrateConfig config);

    ~AnalysisAtom() override;

    std::size_t atomIndex() const { return model_index_; }

    void Calculate(std::size_t step) override;
    QJsonObject Truth() const override;

    std::size_t sigmaFolds() const;
    std::size_t relationshipCount() const;
    std::size_t mappedBondCount() const;
    std::size_t mismatchEventCount() const;
    std::size_t boostCouplingCount() const;
    std::size_t boostSerialCount() const;
    bool oxygenGatePassed() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace h5reader::rediscover
