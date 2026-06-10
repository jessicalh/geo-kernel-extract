#pragma once

#include "RunData.h"

#include <QString>

#include <optional>

namespace h5reader::rediscover {

class StaticRunData {
public:
    static std::optional<RunData> Load(const QString& poseDir, QString* err_out);
};

}  // namespace h5reader::rediscover
