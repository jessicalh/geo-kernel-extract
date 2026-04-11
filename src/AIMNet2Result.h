#pragma once
//
// AIMNet2Result: AIMNet2 neural network charge calculator via libtorch.
//
// Produces per-atom Hirshfeld charges, 256-dim aim embedding, and
// Coulomb EFG tensor from AIMNet2 charges. Decomposes EFG by source
// (backbone, aromatic) using the same kernel as CoulombResult.
//
// Married to AIMNet2 — no abstract interface, no factory pattern.
// The .jpt model is loaded once and shared across all conformations.
//
// charge_sensitivity: computed on first Compute call, stored on Atom
// (Protein level, permanent). 10 bulk 0.1A perturbations, per-atom
// charge variance. An intrinsic property of the atom's chemical
// identity — 112x variation across carbon atoms.
//
// CUDA mandatory. No CPU fallback.
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

// Loaded AIMNet2 model — created once, shared across all conformations.
// Holds the TorchScript module on GPU and the model's short-range cutoff.
struct AIMNet2Model {
    torch::jit::Module module;
    double cutoff = 5.0;        // short-range, read from model attribute
    double cutoff_lr = 15.0;    // long-range DSF Coulomb, from TOML
    int max_nb = 128;           // max short-range neighbours per atom
    int max_nb_lr = 512;        // max long-range neighbours per atom
    torch::Device device = torch::kCUDA;

    // Load model from .jpt file. Reads cutoff from model attribute.
    // cutoff_lr, max_nb, max_nb_lr from CalculatorConfig (TOML).
    static std::unique_ptr<AIMNet2Model> Load(const std::string& jpt_path);
};

class AIMNet2Result : public ConformationResult {
public:
    std::string Name() const override { return "AIMNet2Result"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute AIMNet2 charges, aim embedding, and Coulomb EFG.
    // model is the shared loaded model (created once at startup).
    static std::unique_ptr<AIMNet2Result> Compute(
        ProteinConformation& conf,
        AIMNet2Model& model);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;

    // Build the padded half-neighbour matrix for AIMNet2.
    // Returns (N+1, max_nb) int32 tensor, sentinel = N.
    static torch::Tensor BuildNeighbourMatrix(
        const ProteinConformation& conf,
        double cutoff_sq, int max_nb);

    // Compute Coulomb EFG from AIMNet2 charges, decomposed by source.
    // Same dipolar kernel as CoulombResult.
    static void ComputeCoulombEFG(
        ProteinConformation& conf,
        double cutoff);

    // Compute charge sensitivity on first run (lazy, stored on Atom).
    static void ComputeChargeSensitivity(
        ProteinConformation& conf,
        AIMNet2Model& model);
};

}  // namespace nmr
