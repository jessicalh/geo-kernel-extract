// QtConformation implementation — thin wrapper around QtTrajectoryH5.

#include "QtConformation.h"

#include "QtProtein.h"

namespace h5reader::model {

QtConformation::QtConformation(const QtProtein* protein, std::unique_ptr<h5reader::io::QtTrajectoryH5> h5)
    : protein_(protein)
    , h5_(std::move(h5)) {}

double QtConformation::startTimePicoseconds() const {
    if (!h5_)
        return 0.0;
    const auto& ft = h5_->frameTimes();
    return ft.empty() ? 0.0 : ft.front();
}

double QtConformation::endTimePicoseconds() const {
    if (!h5_)
        return 0.0;
    const auto& ft = h5_->frameTimes();
    return ft.empty() ? 0.0 : ft.back();
}

std::size_t QtConformation::ringCount() const {
    return protein_ ? protein_->ringCount() : 0;
}

}  // namespace h5reader::model
