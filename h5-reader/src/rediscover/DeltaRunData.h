#pragma once

#include "RunData.h"

#include "../io/QtFieldCatalog.gen.h"

#include <QString>
#include <QStringList>

#include <cstddef>
#include <optional>

namespace h5reader::rediscover {

struct DeltaRunData {
    const StaticNpyArray* wt_dia = nullptr;
    const StaticNpyArray* wt_para = nullptr;
    const StaticNpyArray* ala_dia = nullptr;
    const StaticNpyArray* ala_para = nullptr;
    const StaticNpyArray* delta_total = nullptr;
    const StaticNpyArray* delta_dia = nullptr;
    const StaticNpyArray* delta_para = nullptr;
    const StaticNpyArray* delta_apbs = nullptr;
    const StaticNpyArray* delta_ring_proximity = nullptr;
    const StaticNpyArray* delta_scalars = nullptr;

    QStringList field_names;
    QStringList array_names;
    QStringList absent_reasons;
    std::size_t sample_count = 0;
    std::size_t matched_axis_count = 0;
    std::size_t wt_n = 0;
    std::size_t ala_n = 0;
};

QStringList DeltaFieldNames();
bool IsDeltaCatalogField(io::FieldKind kind);

std::optional<DeltaRunData> LoadDeltaRunData(const RunData& run, QString* err_out = nullptr);

}  // namespace h5reader::rediscover
