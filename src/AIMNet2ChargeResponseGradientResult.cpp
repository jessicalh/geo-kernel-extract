#include "AIMNet2PolarisabilityResult.h"

#include "AIMNet2Result.h"
#include "ConformationAtom.h"
#include "ConformationResult.h"
#include "EnrichmentResult.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "SpatialIndexResult.h"

#include <torch/script.h>
#include <torch/cuda.h>

#include <cmath>

namespace nmr {


std::vector<std::type_index> AIMNet2PolarisabilityResult::Dependencies() const {
    return {
        std::type_index(typeid(AIMNet2Result)),
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(EnrichmentResult))
    };
}


std::unique_ptr<AIMNet2PolarisabilityResult>
AIMNet2PolarisabilityResult::Compute(
        ProteinConformation& conf,
        AIMNet2Model& model) {

    OperationLog::Scope scope("AIMNet2PolarisabilityResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();

    if (N == 0) {
        OperationLog::Error("AIMNet2PolarisabilityResult::Compute",
            "Zero atoms — cannot run polarisability backward.");
        return nullptr;
    }

    // Element guard — same set as AIMNet2Result. AIMNet2 has no
    // embedding for Z=0 (Element::Unknown).
    for (size_t i = 0; i < N; ++i) {
        Element e = protein.AtomAt(i).element;
        if (e != Element::H && e != Element::C && e != Element::N &&
            e != Element::O && e != Element::S) {
            OperationLog::Error("AIMNet2PolarisabilityResult::Compute",
                "Atom " + std::to_string(i) + " (" +
                protein.AtomAt(i).pdb_atom_name + " in residue " +
                std::to_string(protein.AtomAt(i).residue_index) +
                ") has Element::Unknown — AIMNet2 has no embedding.");
            return nullptr;
        }
    }

    auto result_ptr = std::make_unique<AIMNet2PolarisabilityResult>();
    result_ptr->conf_ = &conf;

    // GeometryChoice: one summary record naming the scalar objective.
    GeometryChoiceBuilder choices(conf);
    choices.Record(CalculatorId::AIMNet2, 0, "polarisability_backward",
        [&](GeometryChoice& gc) {
            AddNumber(gc, "objective_kind", 1.0, "L2_of_charges");
            AddNumber(gc, "atoms", static_cast<double>(N), "");
        });

    // ----------------------------------------------------------------------
    // Build input tensors. Duplicates AIMNet2Result's tensor build with
    // gradient enabled on `coord` (the leaf we differentiate against).
    // The neighbour matrices come from AIMNet2Result's shared static
    // helper so the input convention matches exactly.
    // ----------------------------------------------------------------------
    const int64_t N1 = static_cast<int64_t>(N + 1);

    auto coord_cpu = torch::zeros({N1, 3}, torch::kFloat32);
    auto coord_acc = coord_cpu.accessor<float, 2>();
    for (size_t i = 0; i < N; ++i) {
        Vec3 p = conf.PositionAt(i);
        coord_acc[i][0] = static_cast<float>(p.x());
        coord_acc[i][1] = static_cast<float>(p.y());
        coord_acc[i][2] = static_cast<float>(p.z());
    }
    // Row N stays zero (sentinel padding row).

    auto numbers_cpu = torch::zeros({N1}, torch::kInt64);
    auto num_acc = numbers_cpu.accessor<int64_t, 1>();
    for (size_t i = 0; i < N; ++i) {
        switch (protein.AtomAt(i).element) {
            case Element::H:  num_acc[i] = 1;  break;
            case Element::C:  num_acc[i] = 6;  break;
            case Element::N:  num_acc[i] = 7;  break;
            case Element::O:  num_acc[i] = 8;  break;
            case Element::S:  num_acc[i] = 16; break;
            default: break;  // unreachable: guarded above
        }
    }

    auto charge_cpu  = torch::zeros({1},  torch::kFloat32);
    auto mol_idx_cpu = torch::zeros({N1}, torch::kInt64);

    const double cutoff_sq = model.cutoff * model.cutoff;
    auto nbmat_cpu = AIMNet2Result::BuildNeighbourMatrix(
        conf, cutoff_sq, model.max_nb);

    const double cutoff_lr_sq = model.cutoff_lr * model.cutoff_lr;
    auto nbmat_lr_cpu = AIMNet2Result::BuildNeighbourMatrix(
        conf, cutoff_lr_sq, model.max_nb_lr);

    auto cutoff_lr_tensor = torch::tensor(
        {static_cast<float>(model.cutoff_lr)}, torch::kFloat32);

    // Move coord to GPU and enable grad. .detach() guarantees a leaf
    // tensor on the device before requesting gradient tracking; setting
    // requires_grad on a non-leaf raises in libtorch.
    auto coord_gpu = coord_cpu.to(model.device).detach().requires_grad_(true);

    c10::Dict<std::string, torch::Tensor> input_dict;
    input_dict.insert("coord",     coord_gpu);
    input_dict.insert("numbers",   numbers_cpu.to(model.device));
    input_dict.insert("charge",    charge_cpu.to(model.device));
    input_dict.insert("mol_idx",   mol_idx_cpu.to(model.device));
    input_dict.insert("nbmat",     nbmat_cpu.to(model.device));
    input_dict.insert("nbmat_lr",  nbmat_lr_cpu.to(model.device));
    input_dict.insert("cutoff_lr", cutoff_lr_tensor.to(model.device));

    OperationLog::Info(LogCalcOther,
        "AIMNet2PolarisabilityResult::Compute",
        "forward pass with grad-tracking enabled (no NoGradGuard)");

    // Forward pass. Gradient tape records the graph from coord_gpu
    // through the AIMNet2 forward to the charges output.
    auto output = model.module.forward({input_dict});
    auto output_dict = output.toGenericDict();

    if (!output_dict.contains("charges")) {
        OperationLog::Error("AIMNet2PolarisabilityResult::Compute",
            "Model output missing 'charges' tensor.");
        return nullptr;
    }
    auto charges_gpu = output_dict.at("charges").toTensor();

    // Scalar objective: L = sum_j q_j^2 over non-sentinel atoms.
    // sum(q) gives ~0 gradient under AIMNet2's charge-conservation
    // projection; the L2 of charges is the cheapest single-pass scalar
    // whose gradient
    //
    //     dL/d(r_i) = sum_j 2 q_j d(q_j)/d(r_i)
    //
    // is non-trivial. The exact d(q_i)/d(r_i) diagonal would require
    // N backward passes; this proxy is the first-test choice. See
    // header comment for the design rationale.
    auto charges_n = charges_gpu.slice(0, 0, static_cast<int64_t>(N));
    auto loss = (charges_n * charges_n).sum();

    loss.backward();

    if (!coord_gpu.grad().defined()) {
        OperationLog::Error("AIMNet2PolarisabilityResult::Compute",
            "coord.grad undefined after backward — autograd path is "
            "broken in the loaded .jpt model. Re-export with grad-"
            "tracking enabled.");
        return nullptr;
    }

    auto grad_cpu = coord_gpu.grad().to(torch::kCPU, torch::kFloat64);
    auto grad_acc = grad_cpu.accessor<double, 2>();

    double max_norm = 0.0;
    double sum_norm = 0.0;
    for (size_t i = 0; i < N; ++i) {
        Vec3 v(grad_acc[i][0], grad_acc[i][1], grad_acc[i][2]);
        const double s = v.norm();
        auto& ca = conf.MutableAtomAt(i);
        ca.aimnet2_polarisability_vector = v;
        ca.aimnet2_polarisability_scalar = s;
        if (s > max_norm) max_norm = s;
        sum_norm += s;
    }
    const double mean_norm = sum_norm / static_cast<double>(N);

    OperationLog::Info(LogCalcOther,
        "AIMNet2PolarisabilityResult::Compute",
        std::to_string(N) + " atoms; |dL/dr| max=" +
        std::to_string(max_norm) + " mean=" + std::to_string(mean_norm));

    return result_ptr;
}


int AIMNet2PolarisabilityResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    // (N, 3) float64 — per-atom polarisability gradient vector
    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const Vec3& v = conf.AtomAt(i).aimnet2_polarisability_vector;
            data[i * 3 + 0] = v.x();
            data[i * 3 + 1] = v.y();
            data[i * 3 + 2] = v.z();
        }
        NpyWriter::WriteFloat64(
            output_dir + "/aimnet2_polarisability.npy",
            data.data(), N, 3);
        written++;
    }

    // (N,) float64 — per-atom polarisability scalar (L2 norm of vector)
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i) {
            data[i] = conf.AtomAt(i).aimnet2_polarisability_scalar;
        }
        NpyWriter::WriteFloat64(
            output_dir + "/aimnet2_polarisability_scalar.npy",
            data.data(), N);
        written++;
    }

    return written;
}


}  // namespace nmr
