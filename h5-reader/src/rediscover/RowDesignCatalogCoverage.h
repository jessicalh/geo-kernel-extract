#pragma once

#include "RowDesign.h"
#include "RunData.h"
#include "ScopedProducerCatalog.h"

#include "../io/QtFieldCatalog.gen.h"

#include <QJsonArray>
#include <QString>
#include <QStringList>

#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::rediscover {

struct RowDesignStats;

struct CatalogRowColumn {
    io::FieldKind kind;
    int component = 0;
    RowColumnSpec spec;
};

struct CatalogCoverageArtifacts {
    QJsonArray fields;
    QStringList sidecarFiles;
};

const std::vector<CatalogRowColumn>& CatalogRowColumns();
bool IsCatalogRowColumn(io::FieldKind kind);

bool EnsureCatalogRowColumnArrays(RunData& run, QString* err_out = nullptr);

std::optional<double> CatalogRowValue(const RunData& run,
                                      io::FieldKind kind,
                                      int component,
                                      std::size_t atom,
                                      std::size_t frame,
                                      std::size_t original);

bool WriteCatalogCoverageArtifacts(const QString& outDir,
                                   const std::vector<RunData>& runs,
                                   const std::vector<RowColumnSpec>& schema,
                                   const RowDesignStats& stats,
                                   CatalogCoverageArtifacts* artifacts,
                                   QString* err_out = nullptr);

}  // namespace h5reader::rediscover
