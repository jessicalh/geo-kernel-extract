// RediscoveryExtraction — the thin extraction interface. NOT a framework: no
// scheduler, no dependency graph, no phases. The driver runs a
// vector<unique_ptr<RediscoveryExtraction>> in a plain loop (DESIGN.md). This
// is the QtRing-hierarchy idiom kept thin: each concrete extraction answers
// for its own stratum + schema + walk.

#pragma once

#include "AnalysisBody.h"
#include "RecordSink.h"

#include <QString>

namespace h5reader::rediscover {

class RediscoveryExtraction {
public:
    virtual ~RediscoveryExtraction() = default;

    // "ring_current" | "mcconnell" — names the output files + the stratum.
    virtual QString name() const = 0;

    // The column schema for this extraction's two row kinds.
    virtual FeatureSchema schema() const = 0;

    // Walk (target atom in this stratum, DFT frame row) over the run, building
    // one NeighborhoodRecord per case and handing it to the sink (both row
    // kinds). Returns the number of cases (records) emitted.
    virtual std::size_t extract(const Body& body, RecordSink& sink) const = 0;
};

}  // namespace h5reader::rediscover
