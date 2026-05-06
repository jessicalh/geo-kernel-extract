#pragma once
//
// tests/test_naming_applicator_access.h
//
// Test-only access shim for NamingApplicator's private constructor.
// Codex Finding D3 (2026-05-06): the prior round's
// `NamingApplicator::CustomRules` constructor was public — declared
// "test-only" by comment but reachable from any production caller.
// This shim restricts test-only construction to friends of
// `nmr::test::NamingApplicatorTestAccess`, which is itself defined
// only in the tests/ tree.
//
// Production code MUST call `GlobalNamingApplicator()`; never include
// this header outside `tests/` and never link this header into
// `libnmr_shielding.a`.

#include "NamingRegistry.h"

namespace nmr {
namespace test {

class NamingApplicatorTestAccess {
public:
    /// Build a NamingApplicator with the supplied rule list installed
    /// in place of the production rules. The canonicality oracle and
    /// the Resolve() body are unchanged.
    static NamingApplicator MakeWithRules(std::vector<NamingRule> rules) {
        NamingApplicator::CustomRules custom;
        custom.rules = std::move(rules);
        return NamingApplicator{std::move(custom)};
    }
};

}  // namespace test
}  // namespace nmr
