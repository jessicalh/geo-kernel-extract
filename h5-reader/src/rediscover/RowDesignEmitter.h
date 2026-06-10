#pragma once

#include "AnalysisBody.h"
#include "ConditioningSpec.h"
#include "RowDesignSink.h"

#include <cstddef>
#include <vector>

namespace h5reader::rediscover {

struct RowDesignStats {
    std::size_t rows = 0;
    std::size_t atoms = 0;
    std::size_t dftRows = 0;
    std::size_t dftPresent = 0;
    std::size_t phiPsiPresent = 0;
    std::size_t phiPsiEligible = 0;
    std::size_t phiPsiFiniteEligible = 0;
    std::size_t dsspPresent = 0;
    std::size_t embeddingPresent = 0;
    std::vector<std::size_t> populatedCounts;
};

bool EnsureRowDesignRingArrays(RunData& run, QString* err_out = nullptr);

RowDesignStats RunRowDesignEmit(const Body& body,
                                RowDesignSink& sink,
                                const ConditioningSpec& spec);

}  // namespace h5reader::rediscover
