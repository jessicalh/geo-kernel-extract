#pragma once
//
// AIMNet2ChargeResponseGradientResult: per-atom charge-response gradient
// d(Σ q²)/d(r) via autograd through the AIMNet2 TorchScript model.
// NOT a Buckingham polarisability α = ∂μ/∂E.
//
// Returned quantity: the per-atom 3-vector dL/d(r_i) and its L2 norm,
// NOT the exact charge-response diagonal d(q_i)/d(r_i).
//
// Computes coord.grad in a single backward pass on
//
//     L = sum_j q_j^2     (over non-sentinel atoms)
//
// The naive choice L = sum(q) gives a near-zero gradient because AIMNet2
// enforces total-charge conservation (sum_j q_j is fixed by
// construction; d(const)/d(coord) = 0). The L2-of-charges scalar is
// the cheapest single-pass objective whose gradient is non-trivial,
// and is proportional to d(q_i)/d(r_i) when q_i dominates the local
// contribution.
//
// The exact diagonal d(q_i)/d(r_i) per atom would require N backward
// passes (one per atom); the single-pass L2-of-charges proxy is the
// deliberate design choice. Future variants can swap the scalar
// objective without changing this Result's class shape.
//
// Prerequisite: the .jpt model must support requires_grad on coords.
//
// Lifecycle: ConformationResult subclass. Always-on for both the
// non-trajectory and trajectory pipelines; every run triggers an
// autograd backward (per frame in trajectory mode). Depends on
// AIMNet2Result for attach ordering, but does NOT share state with
// it — runs its own grad-tracking forward (no NoGradGuard), then a
// single backward. Cost is roughly equal to AIMNet2Result itself,
// since this Result re-runs the full forward pass with gradient
// tracking enabled.
//
// CUDA mandatory through the shared AIMNet2Model. No CPU fallback.
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

class AIMNet2ChargeResponseGradientResult : public ConformationResult {
public:
    std::string Name() const override { return "AIMNet2ChargeResponseGradientResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: run the autograd backward pass and return the populated
    // Result. The caller attaches it to the conformation. Returns
    // nullptr on checked failures (zero atoms, unsupported elements,
    // missing 'charges' tensor in model output, undefined or non-finite
    // coord.grad).
    static std::unique_ptr<AIMNet2ChargeResponseGradientResult> Compute(
        ProteinConformation& conf,
        AIMNet2Model& model);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr
