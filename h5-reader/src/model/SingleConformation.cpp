// SingleConformation — implementation.

#include "SingleConformation.h"

#include "../io/QtFieldCatalog.gen.h"
#include "QtConformationSnapshot.h"

namespace h5reader::model {

SingleConformation::SingleConformation(const QtProtein* protein,
                                       std::shared_ptr<const QtConformationSnapshot> pose)
    : Conformation(protein), pose_(std::move(pose)) {
    setObjectName(QStringLiteral("SingleConformation"));
}

SingleConformation::~SingleConformation() = default;

Vec3 SingleConformation::atomPosition(std::size_t /*frame == 0*/, std::size_t atomIdx) const {
    if (!pose_)
        return Vec3::Zero();
    const auto& col = pose_->column(io::FieldKind::Pos);
    if (!col.present || col.cols < 3 || atomIdx >= static_cast<std::size_t>(col.rows))
        return Vec3::Zero();
    const double* r = col.row(atomIdx);
    return Vec3(r[0], r[1], r[2]);
}

std::shared_ptr<const QtConformationSnapshot> SingleConformation::loadSnapshot(std::size_t /*frame == 0*/) {
    return pose_;  // the one pose, already resident; the base caches it once
}

}  // namespace h5reader::model
