#pragma once

#include "../model/Types.h"

namespace h5reader::physics {

struct McConnellSourceGeometry {
    model::BondCategory category = model::BondCategory::Unknown;
    model::Vec3 disp_local = model::Vec3::Zero();
    model::Vec3 bond_axis_local = model::Vec3::Zero();
    double r = 0.0;
};

bool McConnellLiteratureCategory(model::BondCategory category);
double McConnellDeltaChiQ(model::BondCategory category);
double McConnellMolarPrefactor();
double McConnellForwardContributionPpm(const model::SphericalTensor& literatureScaledKernel);
model::SphericalTensor McConnellSourceLiteratureKernelLocal(const McConnellSourceGeometry& source,
                                                            bool* present = nullptr);

}  // namespace h5reader::physics
