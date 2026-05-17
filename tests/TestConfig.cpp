#include "TestConfig.h"

#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"

#include <typeindex>

namespace nmr {
namespace test {

nmr::RunConfiguration BuildTestConfig(
        TestProfile profile,
        const std::string& name,
        std::size_t stride) {
    nmr::RunConfiguration config;
    config.SetName(name);
    auto& opts = config.MutablePerFrameRunOptions();

    // KernelOnly defaults: skip the expensive externals. The test layers
    // its source-calc requirement on top of this.
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;

    // Substrate every kernel test needs.
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::EnrichmentResult));

    if (profile == TestProfile::KernelWithDssp) {
        opts.skip_dssp = false;
    }

    config.SetStride(stride);
    return config;
}

}  // namespace test
}  // namespace nmr
