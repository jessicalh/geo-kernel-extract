#include "DeltaRunData.h"

#include "ScopedProducerCatalog.h"

#include <algorithm>

namespace h5reader::rediscover {
namespace {

QString stem(io::FieldKind kind) {
    const io::FieldSpec& spec = io::FieldSpecFor(kind);
    return QString::fromUtf8(spec.stem.data(), static_cast<qsizetype>(spec.stem.size()));
}

const StaticNpyArray* deltaArray(const RunData& run, io::FieldKind kind) {
    const io::FieldSpec& spec = io::FieldSpecFor(kind);
    if (spec.group != io::FieldGroup::Delta)
        return nullptr;
    return run.producerArray(kind);
}

bool requireArray(const RunData& run,
                  io::FieldKind kind,
                  int requiredCols,
                  const StaticNpyArray** out,
                  DeltaRunData* audit,
                  QString* err_out) {
    const StaticNpyArray* a = deltaArray(run, kind);
    const QString name = stem(kind);
    audit->field_names << name;
    audit->array_names << (a ? a->path : QStringLiteral("<absent>"));
    if (!a) {
        const QString reason = QStringLiteral("%1.npy absent").arg(name);
        audit->absent_reasons << reason;
        if (err_out) *err_out = reason;
        return false;
    }
    if (requiredCols > 0 && a->cols != static_cast<std::size_t>(requiredCols)) {
        const QString reason = QStringLiteral("%1 has %2 cols; expected %3")
                                   .arg(name)
                                   .arg(a->cols)
                                   .arg(requiredCols);
        audit->absent_reasons << reason;
        if (err_out) *err_out = reason;
        return false;
    }
    if (a->rows == 0) {
        const QString reason = QStringLiteral("%1 has zero rows").arg(name);
        audit->absent_reasons << reason;
        if (err_out) *err_out = reason;
        return false;
    }
    *out = a;
    return true;
}

void optionalArray(const RunData& run,
                   io::FieldKind kind,
                   const StaticNpyArray** out,
                   DeltaRunData* audit) {
    const StaticNpyArray* a = deltaArray(run, kind);
    audit->field_names << stem(kind);
    audit->array_names << (a ? a->path : QStringLiteral("<absent>"));
    *out = a;
}

}  // namespace

QStringList DeltaFieldNames() {
    QStringList out;
    for (const io::FieldSpec& spec : io::kFieldCatalog) {
        if (spec.group == io::FieldGroup::Delta)
            out << QString::fromUtf8(spec.stem.data(), static_cast<qsizetype>(spec.stem.size()));
    }
    return out;
}

bool IsDeltaCatalogField(io::FieldKind kind) {
    return io::FieldSpecFor(kind).group == io::FieldGroup::Delta;
}

std::optional<DeltaRunData> LoadDeltaRunData(const RunData& run, QString* err_out) {
    DeltaRunData out;
    if (!requireArray(run, io::FieldKind::WtShieldingDiamagnetic, 9, &out.wt_dia, &out, err_out))
        return std::nullopt;
    if (!requireArray(run, io::FieldKind::WtShieldingParamagnetic, 9, &out.wt_para, &out, err_out))
        return std::nullopt;
    if (!requireArray(run, io::FieldKind::MutShieldingDiamagnetic, 9, &out.ala_dia, &out, err_out))
        return std::nullopt;
    if (!requireArray(run, io::FieldKind::MutShieldingParamagnetic, 9, &out.ala_para, &out, err_out))
        return std::nullopt;
    if (!requireArray(run, io::FieldKind::DeltaShielding, 9, &out.delta_total, &out, err_out))
        return std::nullopt;

    optionalArray(run, io::FieldKind::DeltaShieldingDiamagnetic, &out.delta_dia, &out);
    optionalArray(run, io::FieldKind::DeltaShieldingParamagnetic, &out.delta_para, &out);
    optionalArray(run, io::FieldKind::DeltaAPBS, &out.delta_apbs, &out);
    optionalArray(run, io::FieldKind::DeltaRingProximity, &out.delta_ring_proximity, &out);
    optionalArray(run, io::FieldKind::DeltaScalars, &out.delta_scalars, &out);

    const std::size_t rows = out.wt_dia->rows;
    const StaticNpyArray* required[] = {
        out.wt_dia, out.wt_para, out.ala_dia, out.ala_para, out.delta_total,
    };
    for (const StaticNpyArray* a : required) {
        if (a->rows != rows) {
            if (err_out) {
                *err_out = QStringLiteral("Delta MutationMatchPair axis mismatch: %1 rows vs %2")
                               .arg(a->rows)
                               .arg(rows);
            }
            return std::nullopt;
        }
    }

    out.sample_count = rows;
    out.matched_axis_count = rows;
    out.wt_n = rows;
    out.ala_n = rows;

    // Guard the scoped-catalog invariant: Delta remains excluded from the broad
    // producer catalogue and is only visible through this narrow loader.
    for (const io::FieldSpec& spec : io::kFieldCatalog) {
        if (spec.group == io::FieldGroup::Delta && IsScopedProducerField(spec)) {
            if (err_out) {
                *err_out = QStringLiteral("ScopedProducerCatalog exposes Delta field %1")
                               .arg(stem(spec.kind));
            }
            return std::nullopt;
        }
    }
    return out;
}

}  // namespace h5reader::rediscover
