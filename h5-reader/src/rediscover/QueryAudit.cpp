#include "QueryAudit.h"

#include "RunQuery.h"
#include "SphericalBasis.h"

#include "../io/QtNpyReader.h"
#include "../model/QtProtein.h"

#ifdef signals
#undef signals
#endif
#ifdef slots
#undef slots
#endif

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>

#include <QDir>
#include <QFileInfo>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <unordered_map>
#include <vector>

namespace h5reader::rediscover {
namespace {

constexpr double kAbsTol = 1.0e-9;
constexpr double kRelTol = 1.0e-9;

QString ss8Name(SecondaryStructure8 ss8) {
    switch (ss8) {
    case SecondaryStructure8::H: return QStringLiteral("H");
    case SecondaryStructure8::G: return QStringLiteral("G");
    case SecondaryStructure8::I: return QStringLiteral("I");
    case SecondaryStructure8::E: return QStringLiteral("E");
    case SecondaryStructure8::B: return QStringLiteral("B");
    case SecondaryStructure8::T: return QStringLiteral("T");
    case SecondaryStructure8::S: return QStringLiteral("S");
    case SecondaryStructure8::C: return QStringLiteral("C");
    case SecondaryStructure8::Unknown: return QStringLiteral("unknown");
    }
    return QStringLiteral("unknown");
}

QString ss3Name(SecondaryStructure3 ss3) {
    switch (ss3) {
    case SecondaryStructure3::Helix: return QStringLiteral("helix");
    case SecondaryStructure3::Sheet: return QStringLiteral("sheet");
    case SecondaryStructure3::Coil: return QStringLiteral("coil");
    case SecondaryStructure3::Unknown: return QStringLiteral("unknown");
    }
    return QStringLiteral("unknown");
}

QString signPolicyName(BackboneAngleSignPolicy policy) {
    switch (policy) {
    case BackboneAngleSignPolicy::Unresolved: return QStringLiteral("unresolved");
    case BackboneAngleSignPolicy::NpyIsIupac: return QStringLiteral("npy_is_iupac");
    case BackboneAngleSignPolicy::NpyIsNegatedIupac: return QStringLiteral("npy_is_negated_iupac");
    }
    return QStringLiteral("unresolved");
}

QString dihedralName(DihedralKind kind) {
    switch (kind) {
    case DihedralKind::Phi: return QStringLiteral("phi");
    case DihedralKind::Psi: return QStringLiteral("psi");
    case DihedralKind::Omega: return QStringLiteral("omega");
    case DihedralKind::Chi1: return QStringLiteral("chi1");
    case DihedralKind::Chi2: return QStringLiteral("chi2");
    case DihedralKind::Chi3: return QStringLiteral("chi3");
    case DihedralKind::Chi4: return QStringLiteral("chi4");
    case DihedralKind::Count: return QStringLiteral("count");
    }
    return QStringLiteral("unknown");
}

bool setError(QString* err_out, const QString& err) {
    if (err_out) *err_out = err;
    return false;
}

bool appendStatus(const arrow::Status& status, QString* err_out) {
    if (status.ok()) return true;
    return setError(err_out, QString::fromStdString(status.ToString()));
}

bool appendString(arrow::StringBuilder& builder, const QString& value, QString* err_out) {
    return appendStatus(builder.Append(value.toStdString()), err_out);
}

bool appendDouble(arrow::DoubleBuilder& builder, std::optional<double> value, QString* err_out) {
    if (!value) return appendStatus(builder.AppendNull(), err_out);
    return appendStatus(builder.Append(*value), err_out);
}

bool finishArray(arrow::ArrayBuilder& builder,
                 std::vector<std::shared_ptr<arrow::Array>>* arrays,
                 QString* err_out) {
    std::shared_ptr<arrow::Array> array;
    if (!appendStatus(builder.Finish(&array), err_out)) return false;
    arrays->push_back(std::move(array));
    return true;
}

bool near(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    const double diff = std::abs(a - b);
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return diff <= kAbsTol || diff <= kRelTol * scale;
}

double delta(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return 0.0;
    if (!std::isfinite(a) || !std::isfinite(b))
        return std::numeric_limits<double>::infinity();
    return std::abs(a - b);
}

std::optional<double> npyAt(const io::QtNpyReader::WidenedArray& array,
                            std::size_t row,
                            std::size_t col) {
    if (!array.ok || row >= array.rows || col >= array.cols) return std::nullopt;
    return array.data[row * array.cols + col];
}

bool loadArray(const QString& path,
               std::size_t rows,
               std::size_t cols,
               io::QtNpyReader::WidenedArray* out,
               QString* err_out) {
    if (!QFileInfo::exists(path))
        return setError(err_out, QStringLiteral("required old output sidecar missing: %1").arg(path));
    *out = io::QtNpyReader::ReadArrayWidened(path);
    if (!out->ok) return setError(err_out, out->error);
    if (out->rows != rows || out->cols != cols) {
        return setError(err_out,
                        QStringLiteral("%1 shape mismatch: got (%2,%3), expected (%4,%5)")
                            .arg(path)
                            .arg(static_cast<qulonglong>(out->rows))
                            .arg(static_cast<qulonglong>(out->cols))
                            .arg(static_cast<qulonglong>(rows))
                            .arg(static_cast<qulonglong>(cols)));
    }
    return true;
}

std::optional<model::Mat3> matFromValues(const std::vector<double>& values) {
    if (values.size() != 9) return std::nullopt;
    model::Mat3 out = model::Mat3::Zero();
    for (int i = 0; i < 9; ++i) {
        if (!std::isfinite(values[static_cast<std::size_t>(i)])) return std::nullopt;
        out(i / 3, i % 3) = values[static_cast<std::size_t>(i)];
    }
    return out;
}

std::array<std::optional<double>, 9> directProducerRaw(const Body& body,
                                                       std::size_t atom,
                                                       std::size_t frame) {
    std::array<std::optional<double>, 9> out;
    for (int c = 0; c < 9; ++c) {
        QString reason;
        out[static_cast<std::size_t>(c)] =
            body.catalog.value(body, io::FieldKind::OrcaTotal, atom, frame, c, &reason);
    }
    return out;
}

struct QueryAuditBuilders {
    arrow::Int64Builder row_id;
    arrow::Int64Builder atom_index;
    arrow::Int64Builder h5_row;
    arrow::Int64Builder frame_slot;
    arrow::Int64Builder original_frame_index;
    arrow::StringBuilder selector_labels;
    arrow::StringBuilder iupac_atom_name;
    arrow::StringBuilder chemical_category;
    arrow::StringBuilder ss8;
    arrow::StringBuilder ss3;
    arrow::BooleanBuilder ss_observed;
    arrow::StringBuilder backbone_angle_sign_policy;
    std::array<arrow::DoubleBuilder, static_cast<std::size_t>(DihedralKind::Count)> dihedral_rad;
    std::array<arrow::Int32Builder, static_cast<std::size_t>(DihedralKind::Count)> dihedral_bin;
    std::array<arrow::DoubleBuilder, 9> query_raw;
    std::array<arrow::DoubleBuilder, 9> producer_raw;
    arrow::BooleanBuilder query_matches_producer;
    arrow::DoubleBuilder query_producer_max_abs_delta;
    arrow::DoubleBuilder query_target_T0;
    std::array<arrow::DoubleBuilder, 3> query_target_T1;
    std::array<arrow::DoubleBuilder, 5> query_target_T2;
    arrow::Int64Builder old_row_id;
    arrow::DoubleBuilder old_target_T0;
    std::array<arrow::DoubleBuilder, 3> old_target_T1;
    std::array<arrow::DoubleBuilder, 5> old_target_T2;
    arrow::BooleanBuilder old_matches_query;
    arrow::DoubleBuilder old_target_max_abs_delta;
    arrow::BooleanBuilder old_output_parity_ok;
    arrow::Int32Builder old_output_parity_error_count;
};

std::shared_ptr<arrow::Schema> auditSchema() {
    std::vector<std::shared_ptr<arrow::Field>> fields = {
        arrow::field("row_id", arrow::int64(), false),
        arrow::field("atom_index", arrow::int64(), false),
        arrow::field("h5_row", arrow::int64(), false),
        arrow::field("frame_slot", arrow::int64(), false),
        arrow::field("original_frame_index", arrow::int64(), false),
        arrow::field("selector_labels", arrow::utf8(), true),
        arrow::field("iupac_atom_name", arrow::utf8(), true),
        arrow::field("chemical_category", arrow::utf8(), true),
        arrow::field("ss8", arrow::utf8(), false),
        arrow::field("ss3", arrow::utf8(), false),
        arrow::field("ss_observed", arrow::boolean(), false),
        arrow::field("backbone_angle_sign_policy", arrow::utf8(), false),
    };
    for (int k = 0; k < static_cast<int>(DihedralKind::Count); ++k) {
        const QString name = dihedralName(static_cast<DihedralKind>(k));
        fields.push_back(arrow::field(QStringLiteral("%1_rad").arg(name).toStdString(),
                                      arrow::float64(), true));
        fields.push_back(arrow::field(QStringLiteral("%1_fixed_bin").arg(name).toStdString(),
                                      arrow::int32(), false));
    }
    for (int i = 0; i < 9; ++i) {
        fields.push_back(arrow::field(QStringLiteral("query_orca_total_m%1%2").arg(i / 3).arg(i % 3).toStdString(),
                                      arrow::float64(), true));
    }
    for (int i = 0; i < 9; ++i) {
        fields.push_back(arrow::field(QStringLiteral("producer_orca_total_m%1%2").arg(i / 3).arg(i % 3).toStdString(),
                                      arrow::float64(), true));
    }
    fields.push_back(arrow::field("query_matches_producer", arrow::boolean(), false));
    fields.push_back(arrow::field("query_producer_max_abs_delta", arrow::float64(), false));
    fields.push_back(arrow::field("query_target_T0", arrow::float64(), true));
    for (int i = 0; i < 3; ++i)
        fields.push_back(arrow::field(QStringLiteral("query_target_T1_%1").arg(i).toStdString(),
                                      arrow::float64(), true));
    for (int i = 0; i < 5; ++i)
        fields.push_back(arrow::field(QStringLiteral("query_target_T2_%1").arg(i).toStdString(),
                                      arrow::float64(), true));
    fields.push_back(arrow::field("old_per_atom_row_id", arrow::int64(), false));
    fields.push_back(arrow::field("old_per_atom_target_T0", arrow::float64(), true));
    for (int i = 0; i < 3; ++i)
        fields.push_back(arrow::field(QStringLiteral("old_per_atom_target_T1_%1").arg(i).toStdString(),
                                      arrow::float64(), true));
    for (int i = 0; i < 5; ++i)
        fields.push_back(arrow::field(QStringLiteral("old_per_atom_target_T2_%1").arg(i).toStdString(),
                                      arrow::float64(), true));
    fields.push_back(arrow::field("old_per_atom_matches_query", arrow::boolean(), false));
    fields.push_back(arrow::field("old_per_atom_target_max_abs_delta", arrow::float64(), false));
    fields.push_back(arrow::field("old_output_parity_ok", arrow::boolean(), false));
    fields.push_back(arrow::field("old_output_parity_error_count", arrow::int32(), false));
    return arrow::schema(std::move(fields));
}

bool appendAuditRow(QueryAuditBuilders* b,
                    const Body& body,
                    const QueryRow& row,
                    std::size_t rowId,
                    std::size_t frameSlot,
                    std::size_t oldRowId,
                    const io::QtNpyReader::WidenedArray& oldT0,
                    const io::QtNpyReader::WidenedArray& oldT1,
                    const io::QtNpyReader::WidenedArray& oldT2,
                    bool oldParityOk,
                    int oldParityErrorCount,
                    QueryAuditStats* stats,
                    QString* err_out) {
    if (!appendStatus(b->row_id.Append(static_cast<int64_t>(rowId)), err_out)) return false;
    if (!appendStatus(b->atom_index.Append(static_cast<int64_t>(row.atom)), err_out)) return false;
    if (!appendStatus(b->h5_row.Append(static_cast<int64_t>(row.frame)), err_out)) return false;
    if (!appendStatus(b->frame_slot.Append(static_cast<int64_t>(frameSlot)), err_out)) return false;
    if (!appendStatus(b->original_frame_index.Append(static_cast<int64_t>(row.original_index)), err_out)) return false;
    if (!appendString(b->selector_labels, row.labels.join(QLatin1Char('|')), err_out)) return false;
    if (!appendString(b->iupac_atom_name, body.idx.iupacNames.nameForAtom(row.atom), err_out)) return false;
    if (!appendString(b->chemical_category, body.idx.chemicalCategories.categoryForAtom(row.atom), err_out)) return false;

    const SecondaryStructureState ss = body.idx.secondaryStructure.state(row.atom, row.frame);
    if (!appendString(b->ss8, ss8Name(ss.ss8), err_out)) return false;
    if (!appendString(b->ss3, ss3Name(ss.ss3), err_out)) return false;
    if (!appendStatus(b->ss_observed.Append(ss.observed), err_out)) return false;
    if (!appendString(b->backbone_angle_sign_policy,
                      signPolicyName(body.idx.dihedrals.backboneSignPolicy()), err_out)) return false;

    for (int k = 0; k < static_cast<int>(DihedralKind::Count); ++k) {
        const DihedralState d = body.idx.dihedrals.state(static_cast<DihedralKind>(k), row.atom, row.frame);
        if (!appendDouble(b->dihedral_rad[static_cast<std::size_t>(k)],
                          d.present ? std::optional<double>(d.radians) : std::nullopt,
                          err_out)) return false;
        if (!appendStatus(b->dihedral_bin[static_cast<std::size_t>(k)].Append(d.present ? d.fixed_bin : -1),
                          err_out)) return false;
    }

    const GatheredField* queryRaw = row.fields.empty() ? nullptr : &row.fields.front();
    const bool queryRawPresent = queryRaw && queryRaw->present && queryRaw->values.size() == 9;
    for (int i = 0; i < 9; ++i) {
        std::optional<double> v;
        if (queryRawPresent) v = queryRaw->values[static_cast<std::size_t>(i)];
        if (!appendDouble(b->query_raw[static_cast<std::size_t>(i)], v, err_out)) return false;
    }

    const std::array<std::optional<double>, 9> directRaw =
        directProducerRaw(body, row.atom, row.frame);
    for (int i = 0; i < 9; ++i) {
        if (!appendDouble(b->producer_raw[static_cast<std::size_t>(i)],
                          directRaw[static_cast<std::size_t>(i)], err_out)) return false;
    }

    bool queryProducerMatches = queryRawPresent;
    double queryProducerMaxDelta = 0.0;
    for (int i = 0; i < 9; ++i) {
        if (!queryRawPresent || !directRaw[static_cast<std::size_t>(i)]) {
            queryProducerMatches = false;
            queryProducerMaxDelta = std::numeric_limits<double>::infinity();
            continue;
        }
        const double q = queryRaw->values[static_cast<std::size_t>(i)];
        const double p = *directRaw[static_cast<std::size_t>(i)];
        queryProducerMaxDelta = std::max(queryProducerMaxDelta, delta(q, p));
        if (!near(q, p)) queryProducerMatches = false;
    }
    if (!queryProducerMatches) ++stats->producer_mismatch_rows;
    if (!appendStatus(b->query_matches_producer.Append(queryProducerMatches), err_out)) return false;
    if (!appendStatus(b->query_producer_max_abs_delta.Append(queryProducerMaxDelta), err_out)) return false;

    std::optional<model::SphericalTensor> queryTarget;
    if (queryRawPresent) {
        if (const std::optional<model::Mat3> m = matFromValues(queryRaw->values))
            queryTarget = DecomposeLibrary(*m);
    }
    if (!appendDouble(b->query_target_T0,
                      queryTarget ? std::optional<double>(queryTarget->T0) : std::nullopt,
                      err_out)) return false;
    for (int i = 0; i < 3; ++i) {
        if (!appendDouble(b->query_target_T1[static_cast<std::size_t>(i)],
                          queryTarget ? std::optional<double>(queryTarget->T1[static_cast<std::size_t>(i)])
                                      : std::nullopt,
                          err_out)) return false;
    }
    for (int i = 0; i < 5; ++i) {
        if (!appendDouble(b->query_target_T2[static_cast<std::size_t>(i)],
                          queryTarget ? std::optional<double>(queryTarget->T2[static_cast<std::size_t>(i)])
                                      : std::nullopt,
                          err_out)) return false;
    }

    if (!appendStatus(b->old_row_id.Append(static_cast<int64_t>(oldRowId)), err_out)) return false;
    const std::optional<double> old0 = npyAt(oldT0, oldRowId, 0);
    if (!appendDouble(b->old_target_T0, old0, err_out)) return false;
    std::array<std::optional<double>, 3> old1;
    for (int i = 0; i < 3; ++i) {
        old1[static_cast<std::size_t>(i)] = npyAt(oldT1, oldRowId, static_cast<std::size_t>(i));
        if (!appendDouble(b->old_target_T1[static_cast<std::size_t>(i)],
                          old1[static_cast<std::size_t>(i)], err_out)) return false;
    }
    std::array<std::optional<double>, 5> old2;
    for (int i = 0; i < 5; ++i) {
        old2[static_cast<std::size_t>(i)] = npyAt(oldT2, oldRowId, static_cast<std::size_t>(i));
        if (!appendDouble(b->old_target_T2[static_cast<std::size_t>(i)],
                          old2[static_cast<std::size_t>(i)], err_out)) return false;
    }

    bool oldMatches = queryTarget.has_value() && old0.has_value();
    double oldMaxDelta = 0.0;
    if (queryTarget && old0) {
        oldMaxDelta = std::max(oldMaxDelta, delta(*old0, queryTarget->T0));
        if (!near(*old0, queryTarget->T0)) oldMatches = false;
    } else {
        oldMaxDelta = std::numeric_limits<double>::infinity();
    }
    for (int i = 0; i < 3; ++i) {
        if (!queryTarget || !old1[static_cast<std::size_t>(i)]) {
            oldMatches = false;
            oldMaxDelta = std::numeric_limits<double>::infinity();
            continue;
        }
        const double old = *old1[static_cast<std::size_t>(i)];
        const double query = queryTarget->T1[static_cast<std::size_t>(i)];
        oldMaxDelta = std::max(oldMaxDelta, delta(old, query));
        if (!near(old, query)) oldMatches = false;
    }
    for (int i = 0; i < 5; ++i) {
        if (!queryTarget || !old2[static_cast<std::size_t>(i)]) {
            oldMatches = false;
            oldMaxDelta = std::numeric_limits<double>::infinity();
            continue;
        }
        const double old = *old2[static_cast<std::size_t>(i)];
        const double query = queryTarget->T2[static_cast<std::size_t>(i)];
        oldMaxDelta = std::max(oldMaxDelta, delta(old, query));
        if (!near(old, query)) oldMatches = false;
    }
    if (!oldMatches) ++stats->old_output_mismatch_rows;
    if (!appendStatus(b->old_matches_query.Append(oldMatches), err_out)) return false;
    if (!appendStatus(b->old_target_max_abs_delta.Append(oldMaxDelta), err_out)) return false;
    if (!appendStatus(b->old_output_parity_ok.Append(oldParityOk), err_out)) return false;
    if (!appendStatus(b->old_output_parity_error_count.Append(oldParityErrorCount), err_out)) return false;
    return true;
}

bool writeTable(QueryAuditBuilders* b,
                std::size_t rows,
                const QString& path,
                QString* err_out) {
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    arrays.reserve(auditSchema()->num_fields());
    if (!finishArray(b->row_id, &arrays, err_out)) return false;
    if (!finishArray(b->atom_index, &arrays, err_out)) return false;
    if (!finishArray(b->h5_row, &arrays, err_out)) return false;
    if (!finishArray(b->frame_slot, &arrays, err_out)) return false;
    if (!finishArray(b->original_frame_index, &arrays, err_out)) return false;
    if (!finishArray(b->selector_labels, &arrays, err_out)) return false;
    if (!finishArray(b->iupac_atom_name, &arrays, err_out)) return false;
    if (!finishArray(b->chemical_category, &arrays, err_out)) return false;
    if (!finishArray(b->ss8, &arrays, err_out)) return false;
    if (!finishArray(b->ss3, &arrays, err_out)) return false;
    if (!finishArray(b->ss_observed, &arrays, err_out)) return false;
    if (!finishArray(b->backbone_angle_sign_policy, &arrays, err_out)) return false;
    for (int k = 0; k < static_cast<int>(DihedralKind::Count); ++k) {
        if (!finishArray(b->dihedral_rad[static_cast<std::size_t>(k)], &arrays, err_out)) return false;
        if (!finishArray(b->dihedral_bin[static_cast<std::size_t>(k)], &arrays, err_out)) return false;
    }
    for (auto& builder : b->query_raw)
        if (!finishArray(builder, &arrays, err_out)) return false;
    for (auto& builder : b->producer_raw)
        if (!finishArray(builder, &arrays, err_out)) return false;
    if (!finishArray(b->query_matches_producer, &arrays, err_out)) return false;
    if (!finishArray(b->query_producer_max_abs_delta, &arrays, err_out)) return false;
    if (!finishArray(b->query_target_T0, &arrays, err_out)) return false;
    for (auto& builder : b->query_target_T1)
        if (!finishArray(builder, &arrays, err_out)) return false;
    for (auto& builder : b->query_target_T2)
        if (!finishArray(builder, &arrays, err_out)) return false;
    if (!finishArray(b->old_row_id, &arrays, err_out)) return false;
    if (!finishArray(b->old_target_T0, &arrays, err_out)) return false;
    for (auto& builder : b->old_target_T1)
        if (!finishArray(builder, &arrays, err_out)) return false;
    for (auto& builder : b->old_target_T2)
        if (!finishArray(builder, &arrays, err_out)) return false;
    if (!finishArray(b->old_matches_query, &arrays, err_out)) return false;
    if (!finishArray(b->old_target_max_abs_delta, &arrays, err_out)) return false;
    if (!finishArray(b->old_output_parity_ok, &arrays, err_out)) return false;
    if (!finishArray(b->old_output_parity_error_count, &arrays, err_out)) return false;

    const std::shared_ptr<arrow::Schema> schema = auditSchema();
    if (arrays.size() != static_cast<std::size_t>(schema->num_fields())) {
        return setError(err_out,
                        QStringLiteral("query audit schema mismatch: %1 arrays for %2 fields")
                            .arg(static_cast<qulonglong>(arrays.size()))
                            .arg(schema->num_fields()));
    }
    const std::shared_ptr<arrow::Table> table = arrow::Table::Make(schema, std::move(arrays), static_cast<int64_t>(rows));
    const QFileInfo fi(path);
    QDir().mkpath(fi.dir().absolutePath());
    const arrow::Result<std::shared_ptr<arrow::io::FileOutputStream>> outFile =
        arrow::io::FileOutputStream::Open(path.toStdString());
    if (!outFile.ok()) return appendStatus(outFile.status(), err_out);
    if (!appendStatus(parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), *outFile, 65536),
                      err_out)) return false;
    return appendStatus((*outFile)->Close(), err_out);
}

}  // namespace

bool WriteQueryAuditParquet(const Body& body,
                            const QString& outputDir,
                            const QueryAuditOptions& options,
                            QueryAuditStats* stats_out,
                            QString* err_out) {
    QueryAuditStats stats;
    if (!body.run.protein)
        return setError(err_out, QStringLiteral("query audit requires a loaded protein"));

    const std::size_t atomCount = body.run.protein->atomCount();
    const std::vector<std::size_t> dftRows = body.run.frameMap.dftRows();
    const std::size_t oldRows = atomCount * dftRows.size();
    std::unordered_map<std::size_t, std::size_t> frameSlot;
    frameSlot.reserve(dftRows.size());
    for (std::size_t i = 0; i < dftRows.size(); ++i)
        frameSlot[dftRows[i]] = i;

    const QString oldDir = options.oldOutputDir.isEmpty() ? outputDir : options.oldOutputDir;
    io::QtNpyReader::WidenedArray oldT0;
    io::QtNpyReader::WidenedArray oldT1;
    io::QtNpyReader::WidenedArray oldT2;
    if (!loadArray(QDir(oldDir).filePath(QStringLiteral("per_atom_substrate_target_T0.npy")),
                   oldRows, 1, &oldT0, err_out)) return false;
    if (!loadArray(QDir(oldDir).filePath(QStringLiteral("per_atom_substrate_target_T1.npy")),
                   oldRows, 3, &oldT1, err_out)) return false;
    if (!loadArray(QDir(oldDir).filePath(QStringLiteral("per_atom_substrate_target_T2.npy")),
                   oldRows, 5, &oldT2, err_out)) return false;

    stats.old_output_parity =
        AuditPerAtomSubstrateAgainstProducer(body, oldDir, options.parityOptions);
    const bool oldParityOk = stats.old_output_parity.ok();
    const int oldParityErrorCount = stats.old_output_parity.errors.size();

    Query query;
    query.domain = TraversalDomain::DftRows;
    query.where = {Selector::TwoPhase(
        QStringLiteral("all_atoms_orca_total"),
        [](const Body& b, std::size_t atom) {
            return b.run.protein && atom < b.run.protein->atomCount();
        },
        [](const Body& b, std::size_t atom, std::size_t frame) {
            QString reason;
            return b.catalog.present(b, io::FieldKind::OrcaTotal, atom, frame, 0, &reason);
        },
        [](const Body&, std::size_t, std::size_t) { return QStringLiteral("target_present"); })};
    query.gather = {FieldRef::Producer(io::FieldKind::OrcaTotal,
                                       QStringLiteral("query_orca_total"))};

    QueryAuditBuilders builders;
    std::size_t queryRowId = 0;
    QString appendErr;
    RunQuery(body, query, [&](const QueryRow& row) {
        if (!appendErr.isEmpty()) return;
        auto it = frameSlot.find(row.frame);
        if (it == frameSlot.end()) {
            appendErr = QStringLiteral("query audit row h5_row=%1 is not in dftRows")
                            .arg(static_cast<qulonglong>(row.frame));
            return;
        }
        const std::size_t slot = it->second;
        const std::size_t oldRowId = slot * atomCount + row.atom;
        if (!appendAuditRow(&builders, body, row, queryRowId, slot, oldRowId,
                            oldT0, oldT1, oldT2, oldParityOk, oldParityErrorCount,
                            &stats, &appendErr))
            return;
        ++queryRowId;
    });
    if (!appendErr.isEmpty()) return setError(err_out, appendErr);
    stats.rows = queryRowId;

    const QString parquetPath = options.parquetPath.isEmpty()
        ? QDir(outputDir).filePath(QStringLiteral("query_audit.parquet"))
        : options.parquetPath;
    if (!writeTable(&builders, stats.rows, parquetPath, err_out)) return false;

    if (stats_out) *stats_out = stats;
    if (!oldParityOk) {
        const QString first = stats.old_output_parity.errors.isEmpty()
            ? QString()
            : stats.old_output_parity.errors.front();
        return setError(err_out,
                        QStringLiteral("old per_atom_substrate parity failed after writing %1: %2 errors; first: %3")
                            .arg(parquetPath)
                            .arg(oldParityErrorCount)
                            .arg(first));
    }
    if (stats.producer_mismatch_rows != 0 || stats.old_output_mismatch_rows != 0) {
        return setError(err_out,
                        QStringLiteral("query audit mismatches after writing %1: producer_rows=%2 old_output_rows=%3")
                            .arg(parquetPath)
                            .arg(static_cast<qulonglong>(stats.producer_mismatch_rows))
                            .arg(static_cast<qulonglong>(stats.old_output_mismatch_rows)));
    }
    return true;
}

}  // namespace h5reader::rediscover
