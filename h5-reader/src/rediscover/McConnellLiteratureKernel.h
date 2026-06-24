// McConnellLiteratureKernel — rediscover-side literature Δχ scaling for the
// McConnell bond-anisotropy tensor. This intentionally lives in rediscover, not
// the frozen library calculator, because the base library emits only unscaled
// Angstrom^-3 geometry.

#pragma once

#include "RediscoverTypes.h"

#include "../model/QtBond.h"

namespace h5reader::rediscover {

// q = Delta chi / (10^-6 cm^3 mol^-1). Provisional single-family
// Williamson-Asakura set, recorded in MCCONNELL_DCHI_LITERATURE.md.
bool McConnellLiteratureCategory(model::BondCategory category);
double McConnellDeltaChiQ(model::BondCategory category);
double McConnellMolarPrefactor();
double McConnellForwardContributionPpm(const model::SphericalTensor& literatureScaledKernel);

// Source-level shielding-signed, literature-scaled McConnell PCS tensor in the
// target atom's local frame. T0 is the signed isotropic ppm scalar; T2 is ppm in
// the rediscover/library component order.
model::SphericalTensor McConnellSourceLiteratureKernelLocal(const SourceSlot& source,
                                                            bool* present = nullptr);

}  // namespace h5reader::rediscover
