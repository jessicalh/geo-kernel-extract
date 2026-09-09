// A one-pose conformation backed by a frame NPY snapshot rather than
// trajectory H5 data.

#pragma once

#include "Conformation.h"

#include <cstddef>
#include <memory>

namespace h5reader::model {

class QtProtein;
class QtConformationSnapshot;

class SingleConformation final : public Conformation {
    Q_OBJECT

public:
    SingleConformation(const QtProtein* protein, std::shared_ptr<const QtConformationSnapshot> pose);
    ~SingleConformation() override;

    std::size_t frameCount() const override { return 1; }
    double timePicoseconds(std::size_t) const override { return 0.0; }
    Vec3 atomPosition(std::size_t frame, std::size_t atomIdx) const override;

protected:
    std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t frame) override;

private:
    std::shared_ptr<const QtConformationSnapshot> pose_;
};

}  // namespace h5reader::model
