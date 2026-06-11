#pragma once

#include "../io/QtFieldCatalog.gen.h"

#include <vector>

namespace h5reader::rediscover {

bool IsScopedProducerField(const io::FieldSpec& spec);
const std::vector<const io::FieldSpec*>& ScopedProducerCatalog();

}  // namespace h5reader::rediscover
