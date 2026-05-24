#include "AIMNet2ChargeResponseGradientResult.h"

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


std::vector<std::type_index> AIMNet2ChargeResponseGradientResult::Dependencies() const {
    return {
        std::type_index(typeid(AIMNet2Result)),
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(EnrichmentResult))
    };
}


std::unique_ptr<AIMNet2ChargeResponseGradientResult>
AIMNet2ChargeResponseGradientResult::Compute(
        ProteinConformation& conf,
        AIMNet2Model& model) {

    OperationLog::Scope scope("AIMNet2ChargeResponseGradientResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();

    if (N == 0) {
        OperationLog::Error("AIMNet2ChargeResponseGradientResult::Compute",
            "Zero atoms — cannot run charge-response backward.");
        return nullptr;
    }

    // Element guard — same set as AIMNet2Result. AIMNet2 has no
    // embedding for Z=0 (Element::Unknown).
    for (size_t i = 0; i < N; ++i) {
        Element e = protein.AtomAt(i).element;
        if (e != Element::H && e != Element::C && e != Element::N &&
            e != Element::O && e != Element::S) {
            OperationLog::Error("AIMNet2ChargeResponseGradientResult::Compute",
                "Atom " + std::to_string(i) + " (" +
                protein.AtomAt(i).pdb_atom_name + " in residue " +
                std::to_string(protein.AtomAt(i).residue_index) +
                ") has Element::Unknown — AIMNet2 has no embedding.");
            return nullptr;
        }
    }

    auto result_ptr = std::make_unique<AIMNet2ChargeResponseGradientResult>();
    result_ptr->conf_ = &conf;

    // GeometryChoice: one summary record naming the scalar objective.
    GeometryChoiceBuilder choices(conf);
    choices.Record(CalculatorId::AIMNet2, 0, "charge_response_gradient_backward",
        [&](GeometryChoice& gc) {
            // objective code: 1.0 = L2_of_charges (label string carries the name)
            AddNumber(gc, "objective_kind", 1.0, "L2_of_charges");
            AddNumber(gc, "atoms", static_cast<double>(N), "");
        });

    // AIMNet2 inputs. Build duplicates AIMNet2Result's tensor build with
    // gradient enabled on `coord` (the leaf we differentiate against);
    // neighbour matrices come from AIMNet2Result's shared static helper so
    // the input convention matches exactly.
    const int64_t padded_atom_count =
        static_cast<int64_t>(N + 1);  // +1 sentinel padding row (stays zero)

    auto coord_cpu = torch::zeros({padded_atom_count, 3}, torch::kFloat32);  // positions in Å
    auto coord_acc = coord_cpu.accessor<float, 2>();
    for (size_t i = 0; i < N; ++i) {
        Vec3 p = conf.PositionAt(i);
        coord_acc[i][0] = static_cast<float>(p.x());
        coord_acc[i][1] = static_cast<float>(p.y());
        coord_acc[i][2] = static_cast<float>(p.z());
    }
    // Row N stays zero (sentinel padding row).

    auto atomic_numbers_cpu = torch::zeros({padded_atom_count}, torch::kInt64);
    auto atomic_num_acc = atomic_numbers_cpu.accessor<int64_t, 1>();
    for (size_t i = 0; i < N; ++i) {
        switch (protein.AtomAt(i).element) {
            case Element::H:  atomic_num_acc[i] = 1;  break;
            case Element::C:  atomic_num_acc[i] = 6;  break;
            case Element::N:  atomic_num_acc[i] = 7;  break;
            case Element::O:  atomic_num_acc[i] = 8;  break;
            case Element::S:  atomic_num_acc[i] = 16; break;
            default: break;  // unreachable: guarded above
        }
    }

    // total system charge = 0 (neutral; required model input)
    auto total_charge_cpu = torch::zeros({1}, torch::kFloat32);
    // single molecule: all atoms in batch 0
    auto mol_idx_cpu = torch::zeros({padded_atom_count}, torch::kInt64);

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
    input_dict.insert("numbers",   atomic_numbers_cpu.to(model.device));
    input_dict.insert("charge",    total_charge_cpu.to(model.device));
    input_dict.insert("mol_idx",   mol_idx_cpu.to(model.device));
    input_dict.insert("nbmat",     nbmat_cpu.to(model.device));
    input_dict.insert("nbmat_lr",  nbmat_lr_cpu.to(model.device));
    input_dict.insert("cutoff_lr", cutoff_lr_tensor.to(model.device));

    OperationLog::Info(LogCalcOther,
        "AIMNet2ChargeResponseGradientResult::Compute",
        "forward pass with grad-tracking enabled (no NoGradGuard)");

    // Gradient tape records the graph from coord_gpu through the
    // AIMNet2 forward to the charges output.
    auto output = model.module.forward({input_dict});
    auto output_dict = output.toGenericDict();

    if (!output_dict.contains("charges")) {
        OperationLog::Error("AIMNet2ChargeResponseGradientResult::Compute",
            "Model output missing 'charges' tensor.");
        return nullptr;
    }
    auto charges_gpu = output_dict.at("charges").toTensor();

    // Scalar objective: L = sum_j q_j^2 over non-sentinel atoms, with
    //
    //     dL/d(r_i) = sum_j 2 q_j d(q_j)/d(r_i)
    //
    // sum(q) is conservation-degenerate (~0 gradient), so the L2 of
    // charges is the cheapest single-pass objective. See header for
    // the full rationale.
    auto real_atom_charges = charges_gpu.slice(0, 0, static_cast<int64_t>(N));
    auto charge_l2_objective = (real_atom_charges * real_atom_charges).sum();

    charge_l2_objective.backward();

    if (!coord_gpu.grad().defined()) {
        OperationLog::Error("AIMNet2ChargeResponseGradientResult::Compute",
            "coord.grad undefined after backward — autograd path is "
            "broken in the loaded .jpt model. Re-export with grad-"
            "tracking enabled.");
        return nullptr;
    }

    auto grad_cpu = coord_gpu.grad().to(torch::kCPU, torch::kFloat64);
    auto grad_acc = grad_cpu.accessor<double, 2>();

    // Finite-gradient guard: reject the entire frame on any non-finite
    // component. Defined-but-NaN/Inf gradients from a degenerate AIMNet2
    // backward would otherwise attach to the source result, set TS mask=1,
    // and poison the Welford running statistics (mean/M2/std/min/max).
    for (size_t i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (!std::isfinite(grad_acc[i][k])) {
                OperationLog::Error(
                    "AIMNet2ChargeResponseGradientResult::Compute",
                    "non-finite gradient at atom " + std::to_string(i) +
                    " component " + std::to_string(k) +
                    "; rejecting frame so TS mask=0 and Welford skips.");
                return nullptr;
            }
        }
    }

    // store per-atom gradient + accumulate norm stats for the log
    double max_norm = 0.0;
    double sum_norm = 0.0;
    for (size_t i = 0; i < N; ++i) {
        Vec3 grad_vec(grad_acc[i][0], grad_acc[i][1], grad_acc[i][2]);
        const double grad_norm = grad_vec.norm();
        auto& ca = conf.MutableAtomAt(i);
        ca.aimnet2_charge_response_gradient_vector = grad_vec;
        ca.aimnet2_charge_response_gradient_scalar = grad_norm;
        if (grad_norm > max_norm) max_norm = grad_norm;
        sum_norm += grad_norm;
    }
    const double mean_norm = sum_norm / static_cast<double>(N);

    OperationLog::Info(LogCalcOther,
        "AIMNet2ChargeResponseGradientResult::Compute",
        std::to_string(N) + " atoms; |dL/dr| max=" +
        std::to_string(max_norm) + " mean=" + std::to_string(mean_norm));

    return result_ptr;
}


int AIMNet2ChargeResponseGradientResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    // (N, 3) float64 — per-atom charge-response gradient dL/dr [e^2/Å]
    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const Vec3& v = conf.AtomAt(i).aimnet2_charge_response_gradient_vector;
            data[i * 3 + 0] = v.x();
            data[i * 3 + 1] = v.y();
            data[i * 3 + 2] = v.z();
        }
        NpyWriter::WriteFloat64(
            output_dir + "/aimnet2_charge_response_gradient.npy",
            data.data(), N, 3);
        written++;
    }

    // (N,) float64 — L2 norm of the gradient vector [e^2/Å]
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i) {
            data[i] = conf.AtomAt(i).aimnet2_charge_response_gradient_scalar;
        }
        NpyWriter::WriteFloat64(
            output_dir + "/aimnet2_charge_response_gradient_scalar.npy",
            data.data(), N);
        written++;
    }

    return written;
}


}  // namespace nmr
