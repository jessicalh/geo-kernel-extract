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

#include <torch/script.h>
#include <torch/cuda.h>
#include <cmath>
#include <random>

namespace nmr {

// ============================================================================
// Helper: pack SphericalTensor into 9-element array
// [T0, T1x, T1y, T1z, T2_0, T2_1, T2_2, T2_3, T2_4]
// ============================================================================

static void PackST_AN2(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    out[1] = st.T1[0]; out[2] = st.T1[1]; out[3] = st.T1[2];
    out[4] = st.T2[0]; out[5] = st.T2[1]; out[6] = st.T2[2];
    out[7] = st.T2[3]; out[8] = st.T2[4];
}


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
        AIMNet2Model& model) {

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

    // charge: (1,) float32 — net molecular charge
    auto charge_cpu = torch::zeros({1}, torch::kFloat32);

    // mol_idx: (N+1,) int64 — all zeros (single molecule), sentinel = 0
    auto mol_idx_cpu = torch::zeros({N1}, torch::kInt64);

    // nbmat: (N+1, max_nb) int32 — short-range half-neighbour list
    double cutoff_sq = model.cutoff * model.cutoff;
    auto nbmat_cpu = BuildNeighbourMatrix(conf, cutoff_sq, model.max_nb);

    // nbmat_lr: (N+1, max_nb_lr) int32 — long-range half-neighbour list
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

    // ------------------------------------------------------------------
    // 4. Extract charges from output
    // ------------------------------------------------------------------
    auto output_dict = output.toGenericDict();
    auto charges_gpu = output_dict.at("charges").toTensor();
    auto charges_cpu_tensor = charges_gpu.to(torch::kCPU, torch::kFloat64);
    auto charges_acc = charges_cpu_tensor.accessor<double, 1>();

    // Store charges on ConformationAtom (first N elements, skip sentinel)
    for (size_t i = 0; i < N; ++i) {
        conf.MutableAtomAt(i).aimnet2_charge = charges_acc[i];
    }

    // ------------------------------------------------------------------
    // 5. Extract aim embedding — HARD FAIL if missing
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
        auto aim_cpu = aim_gpu.to(torch::kCPU, torch::kFloat64);
        auto aim_acc = aim_cpu.accessor<double, 2>();

        int64_t model_dims = aim_cpu.size(1);
        if (model_dims != static_cast<int64_t>(AIMNET2_AIM_DIMS)) {
            OperationLog::Error("AIMNet2Result::Compute",
                "aim embedding has " + std::to_string(model_dims) +
                " dims, expected " + std::to_string(AIMNET2_AIM_DIMS) +
                ". Model architecture mismatch.");
            return nullptr;
        }

        for (size_t i = 0; i < N; ++i) {
            for (size_t d = 0; d < AIMNET2_AIM_DIMS; ++d) {
                conf.MutableAtomAt(i).aimnet2_aim[d] = aim_acc[i][d];
            }
        }
    }

    // ------------------------------------------------------------------
    // 6. Charge sensitivity: lazy init on Atom (topology, permanent)
    // ------------------------------------------------------------------
    if (protein.AtomAt(0).aimnet2_charge_sensitivity < 0.0) {
        ComputeChargeSensitivity(conf, model);
    }

    // ------------------------------------------------------------------
    // 7. Coulomb EFG from AIMNet2 charges
    // ------------------------------------------------------------------
    ComputeCoulombEFG(conf, CalculatorConfig::Get("aimnet2_coulomb_efg_cutoff"));

    OperationLog::Info(LogCalcOther, "AIMNet2Result::Compute",
        std::to_string(N) + " atoms, charges range [" +
        std::to_string(charges_cpu_tensor.min().item<double>()) + ", " +
        std::to_string(charges_cpu_tensor.max().item<double>()) + "], aim embedding extracted");

    return result_ptr;
}


// ============================================================================
// Coulomb EFG from AIMNet2 charges
// ============================================================================

void AIMNet2Result::ComputeCoulombEFG(
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

    // Coulomb EFG sum with AIMNet2 charges
    for (size_t i = 0; i < N; ++i) {
        Vec3 pos_i = conf.PositionAt(i);

        Mat3 EFG_total = Mat3::Zero();
        Mat3 EFG_backbone = Mat3::Zero();
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

            // V_ab = q_j * (3 r_a r_b / r^5 - delta_ab / r^3)
            Mat3 V_j = q_j * (3.0 * r * r.transpose() / r5
                              - Mat3::Identity() / r3);

            EFG_total += V_j;

            if (is_aromatic[j]) {
                EFG_aromatic += V_j;
            } else if (is_backbone[j]) {
                EFG_backbone += V_j;
            }
        }

        // Apply Coulomb constant and traceless projection
        EFG_total    *= COULOMB_KE;
        EFG_backbone *= COULOMB_KE;
        EFG_aromatic *= COULOMB_KE;

        EFG_total    -= (EFG_total.trace() / 3.0) * Mat3::Identity();
        EFG_backbone -= (EFG_backbone.trace() / 3.0) * Mat3::Identity();
        EFG_aromatic -= (EFG_aromatic.trace() / 3.0) * Mat3::Identity();

        auto& ca = conf.MutableAtomAt(i);
        ca.aimnet2_EFG_total = EFG_total;
        ca.aimnet2_EFG_total_spherical = SphericalTensor::Decompose(EFG_total);
        ca.aimnet2_EFG_backbone = EFG_backbone;
        ca.aimnet2_EFG_backbone_spherical = SphericalTensor::Decompose(EFG_backbone);
        ca.aimnet2_EFG_aromatic = EFG_aromatic;
        ca.aimnet2_EFG_aromatic_spherical = SphericalTensor::Decompose(EFG_aromatic);

        // Total shielding contribution = full 9-component SphericalTensor
        ca.aimnet2_shielding_contribution = ca.aimnet2_EFG_total_spherical;
    }
}


// ============================================================================
// Charge sensitivity: bulk perturbations, per-atom variance → Atom
// ============================================================================

void AIMNet2Result::ComputeChargeSensitivity(
        ProteinConformation& conf,
        AIMNet2Model& model) {

    OperationLog::Scope scope("AIMNet2Result::ChargeSensitivity",
        "atoms=" + std::to_string(conf.AtomCount()));

    // const_cast is sanctioned here: charge_sensitivity is an intrinsic
    // property of the atom's chemical identity, stored once on Atom
    // (Protein level). This is the one place a calculator writes to Protein.
    Protein& protein = const_cast<Protein&>(conf.ProteinRef());
    const size_t N = conf.AtomCount();
    const int n_perturb = static_cast<int>(
        CalculatorConfig::Get("aimnet2_sensitivity_n_perturbations"));
    const double displacement = CalculatorConfig::Get("aimnet2_sensitivity_displacement");

    // Accumulate per-atom charge mean and M2 (Welford online)
    std::vector<double> mean(N, 0.0);
    std::vector<double> M2(N, 0.0);

    std::mt19937 rng(42);  // deterministic seed for reproducibility
    std::normal_distribution<double> dist(0.0, displacement);

    // Build static inputs that don't change across perturbations.
    // All atom-level tensors padded to N+1 (same as Compute).
    const int64_t N1 = static_cast<int64_t>(N + 1);
    auto numbers_cpu = torch::zeros({N1}, torch::kInt64);
    auto num_acc = numbers_cpu.accessor<int64_t, 1>();
    for (size_t i = 0; i < N; ++i) {
        switch (protein.AtomAt(i).element) {
            case Element::H:  num_acc[i] = 1;  break;
            case Element::C:  num_acc[i] = 6;  break;
            case Element::N:  num_acc[i] = 7;  break;
            case Element::O:  num_acc[i] = 8;  break;
            case Element::S:  num_acc[i] = 16; break;
            default: break;  // unreachable: Compute guards for Unknown
        }
    }
    auto numbers_gpu = numbers_cpu.to(model.device);
    auto charge_gpu = torch::zeros({1}, torch::dtype(torch::kFloat32).device(model.device));
    auto mol_idx_gpu = torch::zeros({N1},
                                     torch::dtype(torch::kInt64).device(model.device));
    auto cutoff_lr_gpu = torch::tensor({static_cast<float>(model.cutoff_lr)},
                                        torch::dtype(torch::kFloat32).device(model.device));

    // Build neighbour matrices ONCE (0.1A perturbations do not change
    // the neighbour set at 5A cutoff — validated empirically).
    auto nbmat_gpu = BuildNeighbourMatrix(conf,
        model.cutoff * model.cutoff, model.max_nb).to(model.device);
    auto nbmat_lr_gpu = BuildNeighbourMatrix(conf,
        model.cutoff_lr * model.cutoff_lr, model.max_nb_lr).to(model.device);

    torch::NoGradGuard no_grad;

    for (int p = 0; p < n_perturb; ++p) {
        // Perturb all coordinates simultaneously (padded to N+1)
        auto coord_cpu = torch::zeros({N1, 3}, torch::kFloat32);
        auto coord_acc = coord_cpu.accessor<float, 2>();
        for (size_t i = 0; i < N; ++i) {
            Vec3 pos = conf.PositionAt(i);
            coord_acc[i][0] = static_cast<float>(pos.x() + dist(rng));
            coord_acc[i][1] = static_cast<float>(pos.y() + dist(rng));
            coord_acc[i][2] = static_cast<float>(pos.z() + dist(rng));
        }

        c10::Dict<std::string, torch::Tensor> input_dict;
        input_dict.insert("coord",     coord_cpu.to(model.device));
        input_dict.insert("numbers",   numbers_gpu);
        input_dict.insert("charge",    charge_gpu);
        input_dict.insert("mol_idx",   mol_idx_gpu);
        input_dict.insert("nbmat",     nbmat_gpu);
        input_dict.insert("nbmat_lr",  nbmat_lr_gpu);
        input_dict.insert("cutoff_lr", cutoff_lr_gpu);

        auto output = model.module.forward({input_dict});
        auto output_dict = output.toGenericDict();
        auto charges = output_dict.at("charges").toTensor().to(torch::kCPU, torch::kFloat64);
        auto q_acc = charges.accessor<double, 1>();

        // Welford online update
        for (size_t i = 0; i < N; ++i) {
            double q = q_acc[i];
            double delta = q - mean[i];
            mean[i] += delta / (p + 1);
            double delta2 = q - mean[i];
            M2[i] += delta * delta2;
        }
    }

    // Store variance on Atom (permanent)
    for (size_t i = 0; i < N; ++i) {
        protein.MutableAtomAt(i).aimnet2_charge_sensitivity =
            (n_perturb > 1) ? M2[i] / (n_perturb - 1) : 0.0;
    }

    OperationLog::Info(LogCalcOther, "AIMNet2Result::ChargeSensitivity",
        std::to_string(n_perturb) + " perturbations, " +
        std::to_string(N) + " atoms");
}


// ============================================================================
// AIMNet2Result::WriteFeatures
// ============================================================================

int AIMNet2Result::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const size_t N = conf.AtomCount();
    const Protein& protein = conf.ProteinRef();
    int files_written = 0;

    // aimnet2_charges.npy — (N,) float64
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).aimnet2_charge;
        NpyWriter::WriteFloat64(output_dir + "/aimnet2_charges.npy",
                                data.data(), N);
        files_written++;
    }

    // aimnet2_aim.npy — (N, AIMNET2_AIM_DIMS) float64
    {
        std::vector<double> data(N * AIMNET2_AIM_DIMS);
        for (size_t i = 0; i < N; ++i)
            for (size_t d = 0; d < AIMNET2_AIM_DIMS; ++d)
                data[i * AIMNET2_AIM_DIMS + d] = conf.AtomAt(i).aimnet2_aim[d];
        NpyWriter::WriteFloat64(output_dir + "/aimnet2_aim.npy",
                                data.data(), N, AIMNET2_AIM_DIMS);
        files_written++;
    }

    // aimnet2_efg.npy — (N, 9) float64, full SphericalTensor
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i)
            PackST_AN2(conf.AtomAt(i).aimnet2_EFG_total_spherical, &data[i * 9]);
        NpyWriter::WriteFloat64(output_dir + "/aimnet2_efg.npy",
                                data.data(), N, 9);
        files_written++;
    }

    // aimnet2_efg_aromatic.npy — (N, 9) float64
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i)
            PackST_AN2(conf.AtomAt(i).aimnet2_EFG_aromatic_spherical, &data[i * 9]);
        NpyWriter::WriteFloat64(output_dir + "/aimnet2_efg_aromatic.npy",
                                data.data(), N, 9);
        files_written++;
    }

    // aimnet2_efg_backbone.npy — (N, 9) float64
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i)
            PackST_AN2(conf.AtomAt(i).aimnet2_EFG_backbone_spherical, &data[i * 9]);
        NpyWriter::WriteFloat64(output_dir + "/aimnet2_efg_backbone.npy",
                                data.data(), N, 9);
        files_written++;
    }

    // aimnet2_charge_sensitivity.npy — (N,) float64, from Atom (topology)
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = protein.AtomAt(i).aimnet2_charge_sensitivity;
        NpyWriter::WriteFloat64(output_dir + "/aimnet2_charge_sensitivity.npy",
                                data.data(), N);
        files_written++;
    }

    return files_written;
}


}  // namespace nmr
