#include "TensorConventionGuard.h"

namespace h5reader::rediscover {

bool AssertPasShapeConventionEnv(QString* errOut) {
    const QString requested =
        qEnvironmentVariable("H5READER_PAS_SHAPE_CONVENTION_OVERRIDE");
    if (requested.isEmpty()
        || requested == QStringLiteral("haeberlen_distance_from_isotropic_v1")) {
        return true;
    }
    if (errOut) {
        *errOut = QStringLiteral("unsupported PAS shape convention '%1'; expected haeberlen_distance_from_isotropic_v1")
                      .arg(requested);
    }
    return false;
}

}  // namespace h5reader::rediscover
