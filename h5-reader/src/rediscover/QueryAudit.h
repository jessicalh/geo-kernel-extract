// QueryAudit -- thin Parquet consumer proving the foundation query membrane.
//
// The audit table is a reader-side query result with its own schema. It gathers
// producer truth through RunQuery/Catalog::value(FieldKind), joins today's
// per_atom_substrate sidecars as consumed outputs, and never registers a reader
// result as a producer field.

#pragma once

#include "AnalysisBody.h"
#include "SubstrateParity.h"

#include <QString>

#include <cstddef>

namespace h5reader::rediscover {

struct QueryAuditOptions {
    QString oldOutputDir;   // Defaults to outputDir when empty.
    QString parquetPath;    // Defaults to outputDir/query_audit.parquet when empty.
    SubstrateParityOptions parityOptions;
};

struct QueryAuditStats {
    std::size_t rows = 0;
    std::size_t producer_mismatch_rows = 0;
    std::size_t old_output_mismatch_rows = 0;
    SubstrateParityStats old_output_parity;
};

bool WriteQueryAuditParquet(const Body& body,
                            const QString& outputDir,
                            const QueryAuditOptions& options,
                            QueryAuditStats* stats_out = nullptr,
                            QString* err_out = nullptr);

}  // namespace h5reader::rediscover
