#pragma once

#include "../io/QtFieldCatalog.gen.h"

#include <QString>

#include <cstddef>
#include <vector>

namespace h5reader::rediscover {

bool IsExcludedProducerGroup(io::FieldGroup group);
bool IsScopedProducerField(const io::FieldSpec& spec);
bool IsStructuredProducerField(const io::FieldSpec& spec);
std::size_t NominalComponentCount(const io::FieldSpec& spec);

QString FieldStem(const io::FieldSpec& spec);
QString FieldGroupName(io::FieldGroup group);
QString NativeAxisName(io::NativeAxis axis);

const std::vector<const io::FieldSpec*>& ScopedProducerCatalog();

}  // namespace h5reader::rediscover
