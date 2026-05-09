#pragma once
//
// AIMNet2PolarisabilityResult: per-atom charge-polarisation gradient via
// autograd through the AIMNet2 TorchScript model.
//
// Computes coord.grad in a single backward pass on
//
//     L = sum_j q_j^2     (over non-sentinel atoms)
//
// and stores the per-atom 3-vector dL/d(r_i) plus its L2 norm. The
// naive choice L = sum(q) gives a near-zero gradient because AIMNet2
// enforces total-charge conservation (sum_j q_j is fixed by
// construction; d(const)/d(coord) = 0). The L2-of-charges scalar is
// the cheapest single-pass alternative whose gradient is physically
// meaningful — a charge-weighted per-atom polarisability — and is
// proportional to d(q_i)/d(r_i) when q_i dominates the local
// contribution.
//
// The exact diagonal d(q_i)/d(r_i) per atom would require N backward
// passes (one per atom). The single-pass L2-of-charges proxy is the
// first-test choice; future variants can swap the scalar objective
// without changing this Result's class shape.
//
// Per Amendment 2026-05-08(b) (spec/PLANNED_CALCULATORS_2026-04-22.md)
// + the .jpt requires_grad pre-flight check that PASSED on
// 2026-05-09 (tests/data/illustrative_peptides/aimnet2_requires_grad_check.py).
//
// Lifecycle: ConformationResult subclass, gated on a JobSpec test
// flag (--aimnet2-polarisability). Depends on AIMNet2Result for
// attach ordering, but does NOT share state with it — runs its own
// forward pass without NoGradGuard, then a single backward.
//
// Cost: roughly equal to AIMNet2Result itself (~5-6 s on a 4000-atom
// protein), since this Result re-runs the full forward pass with
// gradient tracking enabled. State-sharing with AIMNet2Result would
// save the second forward pass at the cost of coupling the two
// Compute paths; deferred to a later optimisation pass.
//
// CUDA mandatory (inherits from AIMNet2Model). No CPU fallback.
//

#include "ConformationResult.h"
#include "Types.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;
struct AIMNet2Model;  // defined in AIMNet2Result.h

class AIMNet2PolarisabilityResult : public ConformationResult {
public:
    std::string Name() const override { return "AIMNet2PolarisabilityResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: run the autograd backward pass and return the populated
    // Result. The caller attaches it to the conformation. Returns
    // nullptr on any failure (zero atoms, unknown elements, missing
    // 'charges' tensor in model output, undefined coord.grad).
    static std::unique_ptr<AIMNet2PolarisabilityResult> Compute(
        ProteinConformation& conf,
        AIMNet2Model& model);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
