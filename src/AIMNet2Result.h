#pragma once
//
// AIMNet2Result: AIMNet2 neural network charge calculator via libtorch.
//
// Produces per-atom Hirshfeld charges, raw aim embedding, its committed
// element-conditioned projection, output-head energy diagnostics, D3
// diagnostics, and Coulomb E/EFG tensors from AIMNet2 charges. The field
// tensors are decomposed by source (backbone, sidechain, aromatic) using
// the same sign convention and kernel as CoulombResult.
//
// Married to AIMNet2 — no abstract interface, no factory pattern.
// The .jpt model is loaded once and shared across all conformations.
//
// AIMNet2Result does not compute charge sensitivity. The separate
// AIMNet2ChargeResponseGradientResult runs its own grad-tracking
// forward/backward pass and stores the charge-response gradient fields.
//
// CUDA mandatory. No CPU fallback.
//
// FAILURE POLICY: if AIMNet2 is requested and a checked model-load or
// Compute guard fails (CUDA unavailable, model corrupt, aim embedding
// missing, unsupported elements), the result is not attached and
// OperationRunner treats that as a hard failure. Silent degradation is
// not acceptable.
//

#include "ConformationResult.h"
#include "Types.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

#include <torch/script.h>

namespace nmr {

class ProteinConformation;
class Protein;

// AIMNET2_AIM_DIMS is defined in ConformationAtom.h (the canonical location,
// because the array extent depends on it and ConformationAtom.h must not
// include torch headers).

// Loaded AIMNet2 model — created once, shared across all conformations.
// Holds the TorchScript module on GPU and the model's short-range cutoff.
struct AIMNet2Model {
    torch::jit::Module module;
    double cutoff = 5.0;        // short-range, read from model attribute
    double cutoff_lr = 15.0;    // long-range DSF Coulomb, from TOML
    int max_nb = 128;           // max short-range neighbours per atom
    int max_nb_lr = 4096;       // max long-range neighbours per atom
    torch::Device device = torch::kCUDA;

    // Load model from .jpt file. Reads cutoff from model attribute.
    // cutoff_lr, max_nb, max_nb_lr from CalculatorConfig (TOML).
    // Returns nullptr on failure (CUDA unavailable, file not found, etc.)
    // with a clear error logged.
    static std::unique_ptr<AIMNet2Model> Load(const std::string& jpt_path);
};

class AIMNet2Result : public ConformationResult {
public:
    std::string Name() const override { return "AIMNet2Result"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute every AIMNet2-owned ConformationAtom field once.
    // WriteFeatures and trajectory reducers are read-back only.
    // model is the shared loaded model (created once at startup).
    // Returns nullptr for checked failures — never silently degrades.
    static std::unique_ptr<AIMNet2Result> Compute(
        ProteinConformation& conf,
        AIMNet2Model& model,
        int net_charge);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    double EnergyLocalSum() const { return energy_local_sum_; }
    double EnergyLongRangeCoulomb() const { return energy_lrcoulomb_; }
    double EnergyDftD3() const { return energy_dftd3_; }
    double EnergyTotal() const { return energy_total_; }
    double ConditionedNetCharge() const { return conditioned_net_charge_; }
    double NeutralConditioningFlag() const { return neutral_conditioning_flag_; }

    // Build the padded symmetric neighbour matrix for AIMNet2.
    // Returns (N+1, max_nb) int32 tensor, sentinel = N.
    //
    // Public so AIMNet2ChargeResponseGradientResult can reuse the convention
    // when building its own input dict (chained via Dependencies()
    // for ordering, but does not share state with this Result).
    static torch::Tensor BuildNeighbourMatrix(
        const ProteinConformation& conf,
        double cutoff_sq, int max_nb);

private:
    const ProteinConformation* conf_ = nullptr;

    // Protein-level one-molecule diagnostics. Per-atom products live on
    // ConformationAtom under the calculator's one-writer discipline.
    double energy_local_sum_ = 0.0;
    double energy_lrcoulomb_ = 0.0;
    double energy_dftd3_ = 0.0;
    double energy_total_ = 0.0;
    double conditioned_net_charge_ = 0.0;
    double neutral_conditioning_flag_ = 0.0;

    // Compute Coulomb E/EFG from AIMNet2 charges, decomposed by source.
    // Same sign and dipolar kernel as CoulombResult. Returns false after
    // logging if a non-finite model product would poison the result.
    static bool ComputeCoulombEFG(
        ProteinConformation& conf,
        double cutoff);
};

}  // namespace nmr
