#include "AIMNet2Result.h"
#include "CalculatorConfig.h"
#include "ConformationResult.h"
#include "EnrichmentResult.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "SpatialIndexResult.h"
#include "generated/AIMNet2AimProjection.h"

#include <torch/script.h>
#include <torch/cuda.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>

namespace nmr {

namespace {

// Exact constants used by scripts/extract_aimnet_extra_20260704.py.
constexpr double kAimnet2BohrPerAngstrom = 1.8897261258369282;
constexpr double kAimnet2HartreeToEv = 13.605693012183622;

c10::Dict<std::string, torch::Tensor> CloneTensorDict(
        const c10::impl::GenericDict& input) {
    c10::Dict<std::string, torch::Tensor> out;
    for (const auto& item : input) {
        if (item.value().isTensor()) {
            out.insert(item.key().toStringRef(),
                       item.value().toTensor().clone());
        }
    }
    return out;
}

c10::Dict<std::string, torch::Tensor> TensorDictFromIValue(
        const c10::IValue& value) {
    c10::Dict<std::string, torch::Tensor> out;
    const auto generic = value.toGenericDict();
    for (const auto& item : generic) {
        if (item.value().isTensor()) {
            out.insert(item.key().toStringRef(), item.value().toTensor());
        }
    }
    return out;
}

bool RunOutputHead(torch::jit::Module& outputs,
                   const char* head_name,
                   c10::Dict<std::string, torch::Tensor>& dict) {
    try {
        auto head = outputs.attr(head_name).toModule();
        dict = TensorDictFromIValue(head.forward({dict}));
        return true;
    } catch (const c10::Error& e) {
        OperationLog::Error(
            "AIMNet2Result::Compute",
            std::string("AIMNet2 output head failed: ") + head_name +
            " — " + e.what());
    } catch (const std::exception& e) {
        OperationLog::Error(
            "AIMNet2Result::Compute",
            std::string("AIMNet2 output head failed: ") + head_name +
            " — " + e.what());
    }
    return false;
}

double TensorScalar(const torch::Tensor& tensor) {
    return tensor.reshape({-1}).select(0, 0)
        .to(torch::kCPU, torch::kFloat64).item<double>();
}

torch::Tensor Dftd3AttrAsTensor(torch::jit::Module& dftd3,
                                const char* name,
                                const torch::Device& device) {
    const c10::IValue value = dftd3.attr(name);
    if (value.isTensor()) return value.toTensor().to(device);
    if (value.isDouble()) {
        return torch::tensor(
            value.toDouble(),
            torch::TensorOptions().dtype(torch::kFloat32).device(device));
    }
    if (value.isInt()) {
        return torch::tensor(
            static_cast<double>(value.toInt()),
            torch::TensorOptions().dtype(torch::kFloat32).device(device));
    }
    throw std::runtime_error(std::string("D3 attribute is not numeric: ") +
                             name);
}

std::size_t ProjectionElementSlot(Element element) {
    switch (element) {
        case Element::H:
            return static_cast<std::size_t>(
                Aimnet2AimProjectionElementSlot::H);
        case Element::C:
            return static_cast<std::size_t>(
                Aimnet2AimProjectionElementSlot::C);
        case Element::N:
            return static_cast<std::size_t>(
                Aimnet2AimProjectionElementSlot::N);
        case Element::O:
            return static_cast<std::size_t>(
                Aimnet2AimProjectionElementSlot::O);
        case Element::S:
            return static_cast<std::size_t>(
                Aimnet2AimProjectionElementSlot::S);
        case Element::Unknown:
            break;
    }
    throw std::runtime_error("unsupported element in AIMNet2 projection");
}

bool TensorPrefixIsFinite(const torch::Tensor& tensor,
                          std::int64_t rows) {
    return torch::isfinite(tensor.narrow(0, 0, rows)).all()
        .to(torch::kCPU).item<bool>();
}

}  // namespace

// ============================================================================
// EFG packing note
// ============================================================================

// PackST_AN2 (9-component SphericalTensor → flat double[9]) removed
// 2026-05-18 with the EFG schema rev: AIMNet2 EFG is T2-only 5-component
// post-schema-rev. WriteFeatures now packs T2[5] directly inline via
// the write_efg_t2 lambda. Removing the helper avoids unused-function
// warnings and signals the architectural change at the boundary.


// ============================================================================
// AIMNet2Model: load once, share across conformations
// ============================================================================

std::unique_ptr<AIMNet2Model> AIMNet2Model::Load(const std::string& jpt_path) {
    OperationLog::Scope scope("AIMNet2Model::Load", jpt_path);

    // Validate CUDA before attempting anything
    if (!torch::cuda::is_available()) {
        OperationLog::Error("AIMNet2Model::Load",
            "CUDA is not available. AIMNet2 requires a CUDA GPU.");
        return nullptr;
    }

    auto m = std::make_unique<AIMNet2Model>();

    try {
        m->module = torch::jit::load(jpt_path, torch::kCUDA);
        m->module.eval();
    } catch (const c10::Error& e) {
        OperationLog::Error("AIMNet2Model::Load",
            "Failed to load TorchScript model: " + std::string(e.what()));
        return nullptr;
    } catch (const std::exception& e) {
        OperationLog::Error("AIMNet2Model::Load",
            "Failed to load model: " + std::string(e.what()));
        return nullptr;
    }

    // Read cutoff from model attribute (confirmed: returns 5.0 for wb97m)
    try {
        m->cutoff = m->module.attr("cutoff").toDouble();
    } catch (...) {
        OperationLog::Error("AIMNet2Model::Load",
            "Model does not expose 'cutoff' attribute. Wrong .jpt file?");
        return nullptr;
    }

    // Long-range and neighbour limits from TOML
    m->cutoff_lr = CalculatorConfig::Get("aimnet2_cutoff_lr");
    m->max_nb    = static_cast<int>(CalculatorConfig::Get("aimnet2_max_nb"));
    m->max_nb_lr = static_cast<int>(CalculatorConfig::Get("aimnet2_max_nb_lr"));

    OperationLog::Info(LogCalcOther, "AIMNet2Model::Load",
        "cutoff=" + std::to_string(m->cutoff) +
        " cutoff_lr=" + std::to_string(m->cutoff_lr) +
        " max_nb=" + std::to_string(m->max_nb) +
        " max_nb_lr=" + std::to_string(m->max_nb_lr));

    return m;
}


// ============================================================================
// Neighbour matrix construction
// ============================================================================

torch::Tensor AIMNet2Result::BuildNeighbourMatrix(
        const ProteinConformation& conf,
        double cutoff_sq, int max_nb) {

    const size_t N = conf.AtomCount();

    // (N+1, max_nb) int32, initialised to sentinel N (matches Python numba output)
    auto nbmat = torch::full({static_cast<int64_t>(N + 1), static_cast<int64_t>(max_nb)},
                             static_cast<int32_t>(N), torch::kInt32);
    auto acc = nbmat.accessor<int32_t, 2>();

    // Per-atom neighbour count
    std::vector<int> count(N + 1, 0);

    const auto& spatial = conf.Result<SpatialIndexResult>();

    for (size_t i = 0; i < N; ++i) {
        Vec3 pos_i = conf.PositionAt(i);

        // nanoflann radius search returns atom indices within radius
        auto neighbours = spatial.AtomsWithinRadius(pos_i,
                                                     std::sqrt(cutoff_sq));
        for (size_t j : neighbours) {
            if (j <= i) continue;  // half-list: only j > i

            // Add j to row i and i to row j
            if (count[i] < max_nb) {
                acc[static_cast<int64_t>(i)][count[i]++] = static_cast<int32_t>(j);
            }
            if (count[j] < max_nb) {
                acc[static_cast<int64_t>(j)][count[j]++] = static_cast<int32_t>(i);
            }
        }
    }

    return nbmat;
}


// ============================================================================
// AIMNet2Result::Dependencies
// ============================================================================

std::vector<std::type_index> AIMNet2Result::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(EnrichmentResult))
    };
}


// ============================================================================
// AIMNet2Result::Compute
// ============================================================================

std::unique_ptr<AIMNet2Result> AIMNet2Result::Compute(
        ProteinConformation& conf,
        AIMNet2Model& model,
        int net_charge) {

    OperationLog::Scope scope("AIMNet2Result::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();

    // Guard: zero atoms (degenerate input)
    if (N == 0) {
        OperationLog::Error("AIMNet2Result::Compute",
            "Zero atoms — cannot run AIMNet2 on empty protein");
        return nullptr;
    }

    // Guard: check all elements are known BEFORE building tensors
    for (size_t i = 0; i < N; ++i) {
        Element e = protein.AtomAt(i).element;
        if (e != Element::H && e != Element::C && e != Element::N &&
            e != Element::O && e != Element::S) {
            OperationLog::Error("AIMNet2Result::Compute",
                "Atom " + std::to_string(i) + " (" +
                protein.AtomAt(i).pdb_atom_name + " in residue " +
                std::to_string(protein.AtomAt(i).residue_index) +
                ") has Element::Unknown. AIMNet2 has no embedding for Z=0. "
                "Fix the topology or exclude this atom.");
            return nullptr;
        }
    }

    auto result_ptr = std::make_unique<AIMNet2Result>();
    result_ptr->conf_ = &conf;

    // ------------------------------------------------------------------
    // 1. Record GeometryChoices
    // ------------------------------------------------------------------
    GeometryChoiceBuilder choices(conf);
    choices.Record(CalculatorId::AIMNet2, 0, "short_range_cutoff",
        [&](GeometryChoice& gc) {
            AddNumber(gc, "cutoff", model.cutoff, "A");
            AddNumber(gc, "source", 0.0, "jpt_model_attribute");
        });
    choices.Record(CalculatorId::AIMNet2, 0, "long_range_cutoff_lr",
        [&](GeometryChoice& gc) {
            AddNumber(gc, "cutoff_lr", model.cutoff_lr, "A");
        });
    choices.Record(CalculatorId::AIMNet2, 0, "neighbour_limits",
        [&](GeometryChoice& gc) {
            AddNumber(gc, "max_nb", static_cast<double>(model.max_nb), "");
            AddNumber(gc, "max_nb_lr", static_cast<double>(model.max_nb_lr), "");
        });

    // ------------------------------------------------------------------
    // 2. Build input tensors
    // ------------------------------------------------------------------
    torch::NoGradGuard no_grad;

    // All atom-level tensors are padded to (N+1) to match nbmat's sentinel row.
    // The model's prepare_data expects coord[nbmat] to be valid, and nbmat's
    // last row indexes position N (the sentinel). Python's pad_input does this.
    const int64_t N1 = static_cast<int64_t>(N + 1);

    // coord: (N+1, 3) float32, last row = zeros (sentinel padding)
    auto coord_cpu = torch::zeros({N1, 3}, torch::kFloat32);
    auto coord_acc = coord_cpu.accessor<float, 2>();
    for (size_t i = 0; i < N; ++i) {
        Vec3 p = conf.PositionAt(i);
        coord_acc[i][0] = static_cast<float>(p.x());
        coord_acc[i][1] = static_cast<float>(p.y());
        coord_acc[i][2] = static_cast<float>(p.z());
    }
    // Row N stays zero (sentinel)

    // numbers: (N+1,) int64 — atomic numbers, last = 0 (sentinel)
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

    // charge: (1,) float32 — the real system net charge (or 0 under the
    // neutral-conditioning knob; see OperationRunner). The model is trained to
    // condition on total charge, so this is its intended input, not a hack.
    auto charge_cpu = torch::full({1}, static_cast<float>(net_charge), torch::kFloat32);

    // mol_idx: (N+1,) int64 — all zeros (single molecule), sentinel = 0
    auto mol_idx_cpu = torch::zeros({N1}, torch::kInt64);

    // nbmat: (N+1, max_nb) int32 — short-range symmetric neighbour list
    double cutoff_sq = model.cutoff * model.cutoff;
    auto nbmat_cpu = BuildNeighbourMatrix(conf, cutoff_sq, model.max_nb);

    // nbmat_lr: (N+1, max_nb_lr) int32 — long-range symmetric neighbour list
    double cutoff_lr_sq = model.cutoff_lr * model.cutoff_lr;
    auto nbmat_lr_cpu = BuildNeighbourMatrix(conf, cutoff_lr_sq, model.max_nb_lr);

    // cutoff_lr: (1,) float32
    auto cutoff_lr_tensor = torch::tensor({static_cast<float>(model.cutoff_lr)},
                                           torch::kFloat32);

    // ------------------------------------------------------------------
    // 3. Build input dictionary and forward pass
    // ------------------------------------------------------------------
    auto to_gpu = [&](torch::Tensor t) { return t.to(model.device); };

    c10::Dict<std::string, torch::Tensor> input_dict;
    input_dict.insert("coord",     to_gpu(coord_cpu));
    input_dict.insert("numbers",   to_gpu(numbers_cpu));
    input_dict.insert("charge",    to_gpu(charge_cpu));
    input_dict.insert("mol_idx",   to_gpu(mol_idx_cpu));
    input_dict.insert("nbmat",     to_gpu(nbmat_cpu));
    input_dict.insert("nbmat_lr",  to_gpu(nbmat_lr_cpu));
    input_dict.insert("cutoff_lr", to_gpu(cutoff_lr_tensor));

    OperationLog::Info(LogCalcOther, "AIMNet2Result",
        "coord=" + std::to_string(coord_cpu.size(0)) + "x" + std::to_string(coord_cpu.size(1)) +
        " numbers=" + std::to_string(numbers_cpu.size(0)) +
        " mol_idx=" + std::to_string(mol_idx_cpu.size(0)) +
        " nbmat=" + std::to_string(nbmat_cpu.size(0)) + "x" + std::to_string(nbmat_cpu.size(1)) +
        " nbmat_lr=" + std::to_string(nbmat_lr_cpu.size(0)) + "x" + std::to_string(nbmat_lr_cpu.size(1)));

    // Forward pass
    auto output = model.module.forward({input_dict});

    const auto output_dict = output.toGenericDict();

    // ------------------------------------------------------------------
    // 4. Replay the model's output heads in their trained order.
    //
    // The full forward output is cloned for the mutable head chain. D3's
    // private diagnostic method intentionally receives the untouched full
    // output dictionary, matching the recovery script. Every per-atom value
    // is widened/stored here; emission is a pure read-back.
    // ------------------------------------------------------------------
    try {
        auto outputs = model.module.attr("outputs").toModule();
        auto dftd3 = outputs.attr("dftd3").toModule();
        auto energy_dict = CloneTensorDict(output_dict);

        if (!RunOutputHead(outputs, "energy_mlp", energy_dict)) {
            return nullptr;
        }
        if (!energy_dict.contains("energy")) {
            OperationLog::Error("AIMNet2Result::Compute",
                "energy_mlp output did not contain 'energy'");
            return nullptr;
        }
        const auto energy_mlp = energy_dict.at("energy").clone();

        if (!RunOutputHead(outputs, "atomic_shift", energy_dict)) {
            return nullptr;
        }
        if (!energy_dict.contains("energy")) {
            OperationLog::Error("AIMNet2Result::Compute",
                "atomic_shift output did not contain 'energy'");
            return nullptr;
        }
        const auto energy_shifted_local = energy_dict.at("energy").clone();

        if (!RunOutputHead(outputs, "atomic_sum", energy_dict)) {
            return nullptr;
        }
        const double local_sum = TensorScalar(energy_dict.at("energy"));

        if (!RunOutputHead(outputs, "lrcoulomb", energy_dict)) {
            return nullptr;
        }
        const double energy_after_lrcoulomb =
            TensorScalar(energy_dict.at("energy"));

        if (!RunOutputHead(outputs, "dftd3", energy_dict)) {
            return nullptr;
        }
        const double final_energy = TensorScalar(energy_dict.at("energy"));

        if (energy_mlp.dim() != 1 || energy_shifted_local.dim() != 1 ||
            energy_mlp.size(0) < N1 ||
            energy_shifted_local.size(0) < N1 ||
            !TensorPrefixIsFinite(energy_mlp, static_cast<std::int64_t>(N)) ||
            !TensorPrefixIsFinite(energy_shifted_local,
                                  static_cast<std::int64_t>(N)) ||
            !std::isfinite(local_sum) ||
            !std::isfinite(energy_after_lrcoulomb) ||
            !std::isfinite(final_energy)) {
            OperationLog::Error("AIMNet2Result::Compute",
                "non-finite or malformed AIMNet2 energy-head output; "
                "result not attached");
            return nullptr;
        }

        result_ptr->energy_local_sum_ = local_sum;
        result_ptr->energy_lrcoulomb_ = energy_after_lrcoulomb - local_sum;
        result_ptr->energy_dftd3_ = final_energy - energy_after_lrcoulomb;
        result_ptr->energy_total_ = final_energy;
        result_ptr->conditioned_net_charge_ = static_cast<double>(net_charge);
        result_ptr->neutral_conditioning_flag_ =
            CalculatorConfig::Get("charge_conditioning_neutral") != 0.0
                ? 1.0 : 0.0;

        const auto energy_mlp_cpu = energy_mlp
            .narrow(0, 0, static_cast<std::int64_t>(N))
            .to(torch::kCPU, torch::kFloat64).contiguous();
        const auto energy_shifted_cpu = energy_shifted_local
            .narrow(0, 0, static_cast<std::int64_t>(N))
            .to(torch::kCPU, torch::kFloat64).contiguous();
        const auto energy_mlp_acc = energy_mlp_cpu.accessor<double, 1>();
        const auto energy_shifted_acc =
            energy_shifted_cpu.accessor<double, 1>();
        for (std::size_t i = 0; i < N; ++i) {
            auto& atom = conf.MutableAtomAt(i);
            atom.aimnet2_energy_mlp =
                energy_mlp_acc[static_cast<std::int64_t>(i)];
            atom.aimnet2_energy_shifted_local =
                energy_shifted_acc[static_cast<std::int64_t>(i)];
        }

        // D3 diagnostics. _calc_c6ij consumes the untouched full model
        // output, exactly as the recovery script does.
        const auto c6_value = dftd3.get_method("_calc_c6ij")({output_dict});
        const auto c6_tuple = c6_value.toTuple();
        if (!c6_tuple || c6_tuple->elements().size() != 2) {
            OperationLog::Error("AIMNet2Result::Compute",
                "DFTD3._calc_c6ij returned an unexpected value");
            return nullptr;
        }
        const auto c6ij = c6_tuple->elements()[0].toTensor();
        const auto d_bohr = c6_tuple->elements()[1].toTensor();
        const auto device = c6ij.device();
        const auto dtype = c6ij.scalar_type();

        const auto numbers = output_dict.at("numbers").toTensor()
            .to(device, torch::kLong);
        const auto nb_lr = output_dict.at("nbmat_lr").toTensor()
            .to(device, torch::kLong);
        const auto nb_flat = nb_lr.reshape({-1});
        const auto mask_lr = torch::logical_not(
            output_dict.at("nb_pad_mask_lr").toTensor()
                .to(device, torch::kBool));

        const auto d_bohr_cn = output_dict.at("d_ij_lr").toTensor()
            .to(device, dtype) * kAimnet2BohrPerAngstrom;
        const auto rcov = Dftd3AttrAsTensor(dftd3, "rcov", device)
            .to(dtype);
        const auto cnmax = Dftd3AttrAsTensor(dftd3, "cnmax", device)
            .to(dtype);
        const auto r4r2 = Dftd3AttrAsTensor(dftd3, "r4r2", device)
            .to(dtype);
        const auto k1 = Dftd3AttrAsTensor(dftd3, "k1", device)
            .to(dtype);
        const auto a1 = Dftd3AttrAsTensor(dftd3, "a1", device)
            .to(dtype);
        const auto a2 = Dftd3AttrAsTensor(dftd3, "a2", device)
            .to(dtype);
        const auto s6 = Dftd3AttrAsTensor(dftd3, "s6", device)
            .to(dtype);
        const auto s8 = Dftd3AttrAsTensor(dftd3, "s8", device)
            .to(dtype);

        const auto rcov_i = rcov.index_select(0, numbers);
        const auto rcov_j = rcov_i.index_select(0, nb_flat)
            .reshape(nb_lr.sizes());
        const auto cn_raw = torch::sum(
            1.0 / (torch::exp(
                ((rcov_i.unsqueeze(1) + rcov_j) / d_bohr_cn - 1.0) * k1)
                + 1.0),
            1);
        const auto cn = torch::minimum(
            cn_raw, cnmax.index_select(0, numbers));

        const auto r4r2_i = r4r2.index_select(0, numbers);
        const auto r4r2_j = r4r2_i.index_select(0, nb_flat)
            .reshape(nb_lr.sizes());
        const auto rrij = 3.0 * r4r2_i.unsqueeze(1) * r4r2_j;
        const auto r0ij = a1 * torch::sqrt(rrij) + a2;
        const auto eij = c6ij * (
            s6 / (torch::pow(d_bohr, 6) + torch::pow(r0ij, 6)) +
            s8 * rrij /
                (torch::pow(d_bohr, 8) + torch::pow(r0ij, 8)));
        const auto e_disp_atom =
            -kAimnet2HartreeToEv * torch::sum(eij, 1);

        const auto valid_counts = torch::sum(mask_lr.to(dtype), 1);
        const auto c6_zeroed = torch::where(
            mask_lr, c6ij, torch::zeros_like(c6ij));
        const auto c6_sum = torch::sum(c6_zeroed, 1);
        const auto c6_mean = torch::where(
            valid_counts > 0.0,
            c6_sum / torch::where(valid_counts > 0.0,
                                  valid_counts,
                                  torch::ones_like(valid_counts)),
            torch::zeros_like(c6_sum));
        const auto c6_masked = torch::where(
            mask_lr, c6ij,
            torch::full_like(c6ij,
                -std::numeric_limits<double>::infinity()));
        auto c6_max = std::get<0>(torch::max(c6_masked, 1));
        c6_max = torch::where(valid_counts > 0.0,
                              c6_max,
                              torch::zeros_like(c6_sum));
        const auto c6_stats = torch::stack({c6_sum, c6_mean, c6_max}, 1);

        const double d3_sum = TensorScalar(e_disp_atom.sum());
        const double d3_increment = result_ptr->energy_dftd3_;
        const double d3_tolerance = std::max(
            1.0e-5,
            1.0e-6 * std::max(1.0, std::abs(d3_increment)));
        if (!std::isfinite(d3_sum) ||
            std::abs(d3_sum - d3_increment) > d3_tolerance) {
            OperationLog::Error("AIMNet2Result::Compute",
                "D3 atom-energy sum mismatch: sum=" +
                std::to_string(d3_sum) + " increment=" +
                std::to_string(d3_increment) + " tolerance=" +
                std::to_string(d3_tolerance));
            return nullptr;
        }

        if (!TensorPrefixIsFinite(e_disp_atom,
                                  static_cast<std::int64_t>(N)) ||
            !TensorPrefixIsFinite(cn, static_cast<std::int64_t>(N)) ||
            !TensorPrefixIsFinite(c6_stats,
                                  static_cast<std::int64_t>(N))) {
            OperationLog::Error("AIMNet2Result::Compute",
                "non-finite AIMNet2 D3 diagnostic; result not attached");
            return nullptr;
        }

        const auto e_disp_cpu = e_disp_atom
            .narrow(0, 0, static_cast<std::int64_t>(N))
            .to(torch::kCPU, torch::kFloat64).contiguous();
        const auto cn_cpu = cn.narrow(0, 0, static_cast<std::int64_t>(N))
            .to(torch::kCPU, torch::kFloat64).contiguous();
        const auto c6_cpu = c6_stats
            .narrow(0, 0, static_cast<std::int64_t>(N))
            .to(torch::kCPU, torch::kFloat64).contiguous();
        const auto e_disp_acc = e_disp_cpu.accessor<double, 1>();
        const auto cn_acc = cn_cpu.accessor<double, 1>();
        const auto c6_acc = c6_cpu.accessor<double, 2>();
        for (std::size_t i = 0; i < N; ++i) {
            const auto row = static_cast<std::int64_t>(i);
            auto& atom = conf.MutableAtomAt(i);
            atom.aimnet2_d3_e_disp_atom = e_disp_acc[row];
            atom.aimnet2_d3_cn = cn_acc[row];
            atom.aimnet2_d3_c6_stats = {
                c6_acc[row][0], c6_acc[row][1], c6_acc[row][2]
            };
        }
    } catch (const c10::Error& e) {
        OperationLog::Error("AIMNet2Result::Compute",
            "AIMNet2 output-head extraction failed: " +
            std::string(e.what()));
        return nullptr;
    } catch (const std::exception& e) {
        OperationLog::Error("AIMNet2Result::Compute",
            "AIMNet2 output-head extraction failed: " +
            std::string(e.what()));
        return nullptr;
    }

    // ------------------------------------------------------------------
    // 5. Extract charges from output
    // ------------------------------------------------------------------
    auto charges_gpu = output_dict.at("charges").toTensor();
    if (charges_gpu.dim() != 1 || charges_gpu.size(0) < N1 ||
        !TensorPrefixIsFinite(charges_gpu, static_cast<std::int64_t>(N))) {
        OperationLog::Error("AIMNet2Result::Compute",
            "non-finite or malformed AIMNet2 charge output; result not attached");
        return nullptr;
    }
    auto charges_cpu_tensor = charges_gpu
        .to(torch::kCPU, torch::kFloat64).contiguous();
    auto charges_acc = charges_cpu_tensor.accessor<double, 1>();

    // Store charges on ConformationAtom (first N elements, skip sentinel)
    for (size_t i = 0; i < N; ++i) {
        conf.MutableAtomAt(i).aimnet2_charge = charges_acc[i];
    }

    // ------------------------------------------------------------------
    // 6. Extract aim embedding and store its committed projection.
    // ------------------------------------------------------------------
    if (!output_dict.contains("aim")) {
        OperationLog::Error("AIMNet2Result::Compute",
            "Model output does not contain 'aim' embedding. "
            "This .jpt model may be an older version that does not expose "
            "the aim tensor. AIMNet2 features will be incomplete — aborting.");
        return nullptr;
    }

    {
        auto aim_gpu = output_dict.at("aim").toTensor();
        auto aim_cpu = aim_gpu.to(torch::kCPU, torch::kFloat32);

        if (aim_cpu.dim() != 2 || aim_cpu.size(0) < N1) {
            OperationLog::Error("AIMNet2Result::Compute",
                "aim embedding has malformed shape; expected (N+1, 256)");
            return nullptr;
        }
        const int64_t model_dims = aim_cpu.size(1);
        if (model_dims != static_cast<int64_t>(AIMNET2_AIM_DIMS)) {
            OperationLog::Error("AIMNet2Result::Compute",
                "aim embedding has " + std::to_string(model_dims) +
                " dims, expected " + std::to_string(AIMNET2_AIM_DIMS) +
                ". Model architecture mismatch.");
            return nullptr;
        }
        auto aim_acc = aim_cpu.accessor<float, 2>();

        for (size_t i = 0; i < N; ++i) {
            auto& atom = conf.MutableAtomAt(i);
            for (size_t d = 0; d < AIMNET2_AIM_DIMS; ++d) {
                const float value = aim_acc[i][d];
                if (!std::isfinite(value)) {
                    OperationLog::Error("AIMNet2Result::Compute",
                        "non-finite aim embedding at atom " +
                        std::to_string(i) + ", dim " + std::to_string(d) +
                        "; result not attached");
                    return nullptr;
                }
                atom.aimnet2_aim[d] = value;
            }

            // The projection is born here exactly once and becomes an
            // AIMNet2Result-owned ConformationAtom field. Neither
            // WriteFeatures nor trajectory statistics re-derive it.
            const std::size_t slot = ProjectionElementSlot(
                protein.AtomAt(i).element);
            for (std::size_t k = 0;
                 k < kAimnet2AimProjectionDims; ++k) {
                double projected = 0.0;
                for (std::size_t d = 0; d < AIMNET2_AIM_DIMS; ++d) {
                    projected += static_cast<double>(
                        kAimnet2AimProjectionBasis[slot][k][d]) *
                        static_cast<double>(atom.aimnet2_aim[d]);
                }
                if (!std::isfinite(projected)) {
                    OperationLog::Error("AIMNet2Result::Compute",
                        "non-finite aim projection at atom " +
                        std::to_string(i) + ", component " +
                        std::to_string(k) + "; result not attached");
                    return nullptr;
                }
                atom.aimnet2_aim_projection[k] =
                    static_cast<float>(projected);
            }
        }
    }

    // ------------------------------------------------------------------
    // 7. Charge-response gradient is a separate Result
    // ------------------------------------------------------------------
    // AIMNet2ChargeResponseGradientResult performs its own grad-tracking
    // forward/backward pass after AIMNet2Result attaches.

    // ------------------------------------------------------------------
    // 8. Coulomb E/EFG from AIMNet2 charges
    // ------------------------------------------------------------------
    if (!ComputeCoulombEFG(
            conf, CalculatorConfig::Get("aimnet2_coulomb_efg_cutoff"))) {
        return nullptr;
    }

    OperationLog::Info(LogCalcOther, "AIMNet2Result::Compute",
        std::to_string(N) + " atoms, charges range [" +
        std::to_string(charges_cpu_tensor
            .narrow(0, 0, static_cast<std::int64_t>(N))
            .min().item<double>()) + ", " +
        std::to_string(charges_cpu_tensor
            .narrow(0, 0, static_cast<std::int64_t>(N))
            .max().item<double>()) +
        "], energy/D3 diagnostics and aim projection extracted");

    return result_ptr;
}


// ============================================================================
// Coulomb E/EFG from AIMNet2 charges
// ============================================================================

bool AIMNet2Result::ComputeCoulombEFG(
        ProteinConformation& conf,
        double cutoff) {

    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();
    const auto& spatial = conf.Result<SpatialIndexResult>();

    const double charge_floor = CalculatorConfig::Get("coulomb_charge_noise_floor");
    const double singularity_guard = CalculatorConfig::Get("singularity_guard_distance");

    // Record GeometryChoices for the EFG computation
    GeometryChoiceBuilder choices(conf);
    choices.Record(CalculatorId::AIMNet2, 0, "aimnet2_coulomb_efg",
        [&](GeometryChoice& gc) {
            AddNumber(gc, "cutoff", cutoff, "A");
            AddNumber(gc, "charge_floor", charge_floor, "");
            AddNumber(gc, "singularity_guard", singularity_guard, "A");
        });

    // Pre-build atom classification (same as CoulombResult)
    std::vector<bool> is_backbone(N, false);
    std::vector<bool> is_aromatic(N, false);

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        auto mark = [&](size_t idx) {
            if (idx != Residue::NONE && idx < N) is_backbone[idx] = true;
        };
        mark(res.N); mark(res.CA); mark(res.C); mark(res.O);
        mark(res.H); mark(res.HA); mark(res.CB);
    }
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (size_t ai : protein.RingAt(ri).atom_indices) {
            if (ai < N) is_aromatic[ai] = true;
        }
    }

    // Coulomb E/EFG sum with AIMNet2 charges.
    // Sign convention: r = target - source and E = k q r / |r|^3.
    for (size_t i = 0; i < N; ++i) {
        Vec3 pos_i = conf.PositionAt(i);

        Vec3 E_total = Vec3::Zero();
        Vec3 E_backbone = Vec3::Zero();
        Vec3 E_sidechain = Vec3::Zero();
        Vec3 E_aromatic = Vec3::Zero();
        Mat3 EFG_total = Mat3::Zero();
        Mat3 EFG_backbone = Mat3::Zero();
        Mat3 EFG_sidechain = Mat3::Zero();
        Mat3 EFG_aromatic = Mat3::Zero();

        auto neighbours = spatial.AtomsWithinRadius(pos_i, cutoff);
        for (size_t j : neighbours) {
            if (j == i) continue;

            double q_j = conf.AtomAt(j).aimnet2_charge;
            if (std::abs(q_j) < charge_floor) continue;

            Vec3 r = pos_i - conf.PositionAt(j);
            double r_mag = r.norm();
            if (r_mag < singularity_guard) continue;

            double r3 = r_mag * r_mag * r_mag;
            double r5 = r3 * r_mag * r_mag;

            const Vec3 E_j = q_j * r / r3;
            // V_ab = q_j * (3 r_a r_b / r^5 - delta_ab / r^3)
            const Mat3 V_j = q_j * (3.0 * r * r.transpose() / r5
                                    - Mat3::Identity() / r3);

            E_total += E_j;
            EFG_total += V_j;

            if (is_aromatic[j]) {
                E_aromatic += E_j;
                EFG_aromatic += V_j;
            } else if (is_backbone[j]) {
                E_backbone += E_j;
                EFG_backbone += V_j;
            } else {
                E_sidechain += E_j;
                EFG_sidechain += V_j;
            }
        }

        // Apply Coulomb constant and traceless projection
        E_total       *= COULOMB_KE;
        E_backbone    *= COULOMB_KE;
        E_sidechain   *= COULOMB_KE;
        E_aromatic    *= COULOMB_KE;
        EFG_total    *= COULOMB_KE;
        EFG_backbone *= COULOMB_KE;
        EFG_sidechain *= COULOMB_KE;
        EFG_aromatic *= COULOMB_KE;

        EFG_total    -= (EFG_total.trace() / 3.0) * Mat3::Identity();
        EFG_backbone -= (EFG_backbone.trace() / 3.0) * Mat3::Identity();
        EFG_sidechain -= (EFG_sidechain.trace() / 3.0) * Mat3::Identity();
        EFG_aromatic -= (EFG_aromatic.trace() / 3.0) * Mat3::Identity();

        // AIMNet2 is a hard-failure calculator. A non-finite model charge
        // or coordinate must make the result absent (and therefore masked by
        // trajectory consumers), never be silently rewritten to zero. This
        // check deliberately does not sanitize the pre-existing EFG fields.
        if (!E_total.allFinite() || !E_backbone.allFinite() ||
            !E_sidechain.allFinite() || !E_aromatic.allFinite() ||
            !EFG_total.allFinite() || !EFG_backbone.allFinite() ||
            !EFG_sidechain.allFinite() || !EFG_aromatic.allFinite()) {
            OperationLog::Error("AIMNet2Result::ComputeCoulombEFG",
                "non-finite Coulomb E/EFG at target atom " +
                std::to_string(i) + "; AIMNet2Result will not attach");
            return false;
        }

        auto& ca = conf.MutableAtomAt(i);
        ca.aimnet2_E_total = E_total;
        ca.aimnet2_E_backbone = E_backbone;
        ca.aimnet2_E_sidechain = E_sidechain;
        ca.aimnet2_E_aromatic = E_aromatic;
        ca.aimnet2_EFG_total = EFG_total;
        ca.aimnet2_EFG_total_spherical = SphericalTensor::Decompose(EFG_total);
        ca.aimnet2_EFG_backbone = EFG_backbone;
        ca.aimnet2_EFG_backbone_spherical = SphericalTensor::Decompose(EFG_backbone);
        ca.aimnet2_EFG_sidechain = EFG_sidechain;
        ca.aimnet2_EFG_sidechain_spherical =
            SphericalTensor::Decompose(EFG_sidechain);
        ca.aimnet2_EFG_aromatic = EFG_aromatic;
        ca.aimnet2_EFG_aromatic_spherical = SphericalTensor::Decompose(EFG_aromatic);

        // Total shielding contribution = full 9-component SphericalTensor
        ca.aimnet2_shielding_contribution = ca.aimnet2_EFG_total_spherical;
    }

    return true;
}


// ============================================================================
// Charge-response gradients live in AIMNet2ChargeResponseGradientResult
// ============================================================================


// ============================================================================
// AIMNet2Result::WriteFeatures
// ============================================================================

int AIMNet2Result::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const size_t N = conf.AtomCount();
    int files_written = 0;
    auto record_write = [&](bool success, const std::string& filename) {
        if (success) {
            ++files_written;
            return;
        }
        OperationLog::Error(
            "AIMNet2Result::WriteFeatures",
            "failed to write " + output_dir + "/" + filename);
    };

    // aimnet2_charges.npy — (N,) float64
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).aimnet2_charge;
        record_write(
            NpyWriter::WriteFloat64(output_dir + "/aimnet2_charges.npy",
                                    data.data(), N),
            "aimnet2_charges.npy");
    }

    // aimnet2_aim.npy — (N, AIMNET2_AIM_DIMS) float32 (native torch precision)
    {
        std::vector<float> data(N * AIMNET2_AIM_DIMS);
        for (size_t i = 0; i < N; ++i)
            for (size_t d = 0; d < AIMNET2_AIM_DIMS; ++d)
                data[i * AIMNET2_AIM_DIMS + d] = conf.AtomAt(i).aimnet2_aim[d];
        record_write(
            NpyWriter::WriteFloat32(output_dir + "/aimnet2_aim.npy",
                                    data.data(), N, AIMNET2_AIM_DIMS),
            "aimnet2_aim.npy");
    }

    // aimnet2_efg.npy — (N, 5) float64, T2 only. AIMNet2 EFG from same
    // Coulomb-style outer-product physics as CoulombResult → T0+T1
    // structural zeros after the explicit traceless projection before
    // storage.
    auto write_efg_t2 = [&](const std::string& name,
                             SphericalTensor ConformationAtom::* member) {
        std::vector<double> data(N * 5);
        for (size_t i = 0; i < N; ++i) {
            (conf.AtomAt(i).*member).PackT2(&data[i*5]);
        }
        record_write(
            NpyWriter::WriteFloat64(output_dir + "/" + name,
                                    data.data(), N, 5),
            name);
    };
    write_efg_t2("aimnet2_efg.npy",          &ConformationAtom::aimnet2_EFG_total_spherical);
    write_efg_t2("aimnet2_efg_aromatic.npy", &ConformationAtom::aimnet2_EFG_aromatic_spherical);
    write_efg_t2("aimnet2_efg_backbone.npy", &ConformationAtom::aimnet2_EFG_backbone_spherical);
    write_efg_t2("aimnet2_efg_sidechain.npy", &ConformationAtom::aimnet2_EFG_sidechain_spherical);

    auto write_vec3 = [&](const std::string& name,
                          Vec3 ConformationAtom::* member) {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const Vec3& value = conf.AtomAt(i).*member;
            data[i*3 + 0] = value.x();
            data[i*3 + 1] = value.y();
            data[i*3 + 2] = value.z();
        }
        record_write(
            NpyWriter::WriteFloat64(
                output_dir + "/" + name, data.data(), N, 3),
            name);
    };
    write_vec3("aimnet2_E.npy", &ConformationAtom::aimnet2_E_total);
    write_vec3("aimnet2_E_backbone.npy", &ConformationAtom::aimnet2_E_backbone);
    write_vec3("aimnet2_E_sidechain.npy", &ConformationAtom::aimnet2_E_sidechain);
    write_vec3("aimnet2_E_aromatic.npy", &ConformationAtom::aimnet2_E_aromatic);

    auto write_scalar = [&](const std::string& name,
                            double ConformationAtom::* member) {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i) {
            data[i] = conf.AtomAt(i).*member;
        }
        record_write(
            NpyWriter::WriteFloat64(output_dir + "/" + name,
                                    data.data(), N),
            name);
    };
    write_scalar("aimnet2_energy_mlp.npy",
                 &ConformationAtom::aimnet2_energy_mlp);
    write_scalar("aimnet2_energy_shifted_local.npy",
                 &ConformationAtom::aimnet2_energy_shifted_local);

    {
        const std::array<double, 6> data = {
            energy_local_sum_,
            energy_lrcoulomb_,
            energy_dftd3_,
            energy_total_,
            conditioned_net_charge_,
            neutral_conditioning_flag_,
        };
        record_write(
            NpyWriter::WriteFloat64(
                output_dir + "/aimnet2_energy_terms.npy",
                data.data(), 1, data.size()),
            "aimnet2_energy_terms.npy");
    }

    write_scalar("aimnet2_d3_e_disp_atom.npy",
                 &ConformationAtom::aimnet2_d3_e_disp_atom);
    write_scalar("aimnet2_d3_cn.npy",
                 &ConformationAtom::aimnet2_d3_cn);

    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const auto& stats = conf.AtomAt(i).aimnet2_d3_c6_stats;
            data[i*3 + 0] = stats[0];
            data[i*3 + 1] = stats[1];
            data[i*3 + 2] = stats[2];
        }
        record_write(
            NpyWriter::WriteFloat64(
                output_dir + "/aimnet2_d3_c6_stats.npy",
                data.data(), N, 3),
            "aimnet2_d3_c6_stats.npy");
    }

    // Pure read-back of the projection born in Compute. In particular,
    // this writer never reads aimnet2_aim or the committed basis table.
    {
        std::vector<float> data(N * kAimnet2AimProjectionDims);
        for (size_t i = 0; i < N; ++i) {
            const auto& projection = conf.AtomAt(i).aimnet2_aim_projection;
            for (size_t k = 0; k < kAimnet2AimProjectionDims; ++k) {
                data[i*kAimnet2AimProjectionDims + k] = projection[k];
            }
        }
        record_write(
            NpyWriter::WriteFloat32(
                output_dir + "/aimnet2_aim_projection.npy",
                data.data(), N, kAimnet2AimProjectionDims),
            "aimnet2_aim_projection.npy");
    }

    // Charge-response-gradient NPYs are written by
    // AIMNet2ChargeResponseGradientResult.

    return files_written;
}


}  // namespace nmr
