#pragma once

#include <QString>

namespace h5reader::rediscover {

bool ValidateCanonical720Root(const QString& root, QString* err_out = nullptr);
bool ValidateCanonical720Pose(const QString& poseDir, QString* err_out = nullptr);
bool ValidateCanonical1p9jRoot(const QString& rootOrLgsPath, QString* err_out = nullptr);

}  // namespace h5reader::rediscover
