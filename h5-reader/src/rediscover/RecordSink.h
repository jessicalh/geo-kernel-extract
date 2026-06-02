// FeatureSchema / RecordSink — the CSV output surface. Two row kinds per case
// (lead decision), so two files per extraction:
//
//   <case>_sources.csv     — un-summed per-(atom, frame, source) rows, for a
//                            sum-pooling / equivariant fitter. One row per
//                            source; identity + frame + local frame + DFT
//                            target columns repeated per source.
//   <case>_aggregated.csv  — one row per (atom, frame): Σ_sources
//                            (3cos²θ−1)/r³ + per-ring-type / per-bond-category
//                            summed features, against the same DFT target. The
//                            well-posed input for a scalar / ridge fitter.
//
// Both files carry the SAME identity + DFT-target columns (the substrate
// contract). The schema names every column; the header row is the schema. The
// files are written atomically via QSaveFile.
//
// Plain data classes, no QObject. The sink streams rows (no full-substrate
// buffering) so the 1P9J run (~450 frames × hundreds of atoms) stays bounded.

#pragma once

#include "RediscoverTypes.h"

#include <QSaveFile>
#include <QStringList>
#include <QString>
#include <QTextStream>

#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::rediscover {

// A named, unit-tagged column. The schema is the column list; it is written as
// the CSV header and is the documentation of the row contract.
struct FeatureColumn {
    QString name;
    QString unit;  // "" for dimensionless / identity columns
};

enum class RelationshipKind { SourceSum, PerAtomFeature };
enum class SourceSchemaKind { None, Ring, Bond, Charge };

struct FeatureSchema {
    QString caseName;                    // "ring_current" | "mcconnell"
    RelationshipKind relationshipKind = RelationshipKind::SourceSum;
    SourceSchemaKind sourceSchemaKind = SourceSchemaKind::None;
    std::vector<FeatureColumn> sourceColumns;      // per-source row schema
    std::vector<FeatureColumn> aggregatedColumns;  // aggregated row schema
    std::size_t maxSourceSlots = 0;      // padding width for the aggregated row
    bool includeBareKernel = true;       // false when no producer cross-check exists
    bool includeMcLitKernel = false;     // McConnell Δχ-scaled local PCS tensor
};

// CSV sink: opens the two files, writes the header row, then accepts rows.
// Commit() flushes + atomically renames both. Returns false (and logs at the
// seam) on any I/O failure.
class RecordSink {
public:
    // outDir + caseName decide the two file paths
    // (<outDir>/<caseName>_sources.csv and _aggregated.csv).
    RecordSink(const QString& outDir, const FeatureSchema& schema);
    ~RecordSink();

    bool Ok() const { return ok_; }

    // Append the per-source rows for one neighbourhood record (one row per
    // source; emits nothing if the record has no sources).
    void WriteSourceRows(const NeighborhoodRecord& rec);

    // Append the single aggregated row for one neighbourhood record.
    // sumAll = Σ over all sources; sumValid = Σ over the producer-valid set
    // (self/bonded rings excluded); nSourcesValid = count of that set; perTypeSums
    // = per-ring-type / per-bond-category sums on the valid set; cutoff_A = the
    // recorded source-discovery cutoff (NaN when sources come from the H5).
    void WriteAggregatedRow(const NeighborhoodRecord& rec, double sumAll, double sumValid,
                            int nSourcesValid, const std::vector<double>& perTypeSums,
                            double cutoff_A);

    // Charge-dipole aggregated row: one row per (atom, frame), with
    // mu = Σ q_i (r_i - r_atom) in the target local frame.
    void WriteChargeDipoleAggregatedRow(const NeighborhoodRecord& rec, const Vec3& mu,
                                        double cutoff_A, const QString& chargeSource,
                                        bool excludeResidue);

    // Flush + atomic-commit both files. After Commit(), the sink is spent.
    bool Commit();

    std::size_t sourceRowsWritten() const { return sourceRows_; }
    std::size_t aggregatedRowsWritten() const { return aggRows_; }
    const QStringList& sidecarFiles() const { return sidecarFiles_; }

private:
    FeatureSchema schema_;
    QString sourcesPath_;
    QString aggregatedPath_;
    QStringList sidecarFiles_;
    std::unique_ptr<QSaveFile> sourcesFile_;
    std::unique_ptr<QSaveFile> aggregatedFile_;
    std::unique_ptr<QTextStream> sourcesOut_;
    std::unique_ptr<QTextStream> aggregatedOut_;
    bool ok_ = false;
    std::size_t sourceRows_ = 0;
    std::size_t aggRows_ = 0;
    std::vector<double> sourceTargetT2_;
    std::vector<double> sourceTargetLocalT2_;
    std::vector<double> sourceBareKernelT2_;
    std::vector<double> aggregatedTargetT2_;
    std::vector<double> aggregatedTargetLocalT2_;
    std::vector<double> aggregatedBareKernelT2_;
};

}  // namespace h5reader::rediscover
