#pragma once
//
// TestConfig: named profile builder for test RunConfigurations.
//
// Test-suite realignment 2026-05-17: per-test BuildBaseConfig helpers with
// bespoke skip-flag combinations across ~20 test files were the surface area
// that today's adversarial review attention got consumed by. This module
// collapses that surface to two named profiles. See
// spec/plan/test-suite-realignment-deferred-2026-05-17.md for the full
// decision record.
//
// Two profiles cover today's TR tests:
//
//   TestProfile::KernelOnly
//     Geometry + SpatialIndex + EnrichmentResult attached. Other heavy
//     calcs (MOPAC, Coulomb, APBS, DSSP, AIMNet2) skipped unless the
//     test explicitly opts in. The test adds its own source-calc
//     RequireConformationResult on top.
//
//   TestProfile::KernelWithDssp
//     Like KernelOnly but DSSP is enabled. For tests of TRs whose source
//     calculator transitively depends on DsspResult (currently:
//     HBondResult-sourced TRs).
//
// Test code shape:
//
//     auto config = nmr::test::BuildTestConfig(
//         nmr::test::TestProfile::KernelOnly, "MyTestName", /*stride=*/300);
//     config.RequireConformationResult(typeid(SourceResult));
//     config.AddTrajectoryResultFactory([](const TrajectoryProtein& tp) {
//         return MyTrajectoryResult::Create(tp);
//     });
//
// Production-shape tests (when one lands) use
// `RunConfiguration::PerFrameExtractionSet()` directly — that's the
// already-named production factory; this module is for kernel-only test
// surface where running the full pipeline isn't justified.
//

#include "RunConfiguration.h"

#include <cstddef>
#include <string>

namespace nmr {
namespace test {

enum class TestProfile {
    // Minimum substrate for kernel TR tests. Geometry + SpatialIndex +
    // EnrichmentResult attached; MOPAC, Coulomb, APBS, DSSP, AIMNet2
    // skipped. Test adds its specific source on top.
    KernelOnly,

    // KernelOnly + DSSP enabled. For HBond-sourced TR tests (HBondResult
    // requires DsspResult).
    KernelWithDssp,
};

// Build a RunConfiguration matching the named profile. The test layers
// its specific source-calc requirement and TR factory on top of the
// returned config.
//
// `name` is for logs (sets RunConfiguration::SetName).
// `stride` defaults to 1; tests that want frame-0-only use a large stride
// (99999), tests that want a handful of frames use ~300.
nmr::RunConfiguration BuildTestConfig(
    TestProfile profile,
    const std::string& name,
    std::size_t stride = 1);

}  // namespace test
}  // namespace nmr
