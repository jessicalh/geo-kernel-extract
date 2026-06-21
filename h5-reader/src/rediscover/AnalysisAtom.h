// AnalysisAtom -- atom branch of the analysis-object model.

#pragma once

#include "AnalysisElement.h"

#include <QJsonObject>
#include <QTextStream>

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
    // Compact per-context accumulator counts (replaces the deleted boost grid):
    // total characterized responses, and contexts the atom occupied.
    std::size_t accumulatorResponseCount() const;
    std::size_t accumulatorContextCount() const;
    bool oxygenGatePassed() const;
    std::size_t writeBoundedSigmaRows(QTextStream& out,
                                      const QString& datasetId,
                                      const QString& proteinId) const;
    std::size_t writeClassicalSourceTermRows(QTextStream& out,
                                             const QString& datasetId,
                                             const QString& proteinId) const;

    static void WriteBoundedSigmaHeader(QTextStream& out);
    static void WriteClassicalSourceTermHeader(QTextStream& out);
    static bool AssertMolCompOrder(QString* errOut = nullptr);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace h5reader::rediscover
