#include "MopacMcConnellResult.h"

#include "McConnellResult.h"
#include "MopacResult.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "ProteinConformation.h"

#include <limits>

namespace nmr {

std::vector<std::type_index> MopacMcConnellResult::Dependencies() const {
    return {
        std::type_index(typeid(MopacResult)),
        std::type_index(typeid(McConnellResult))
    };
}


std::unique_ptr<MopacMcConnellResult> MopacMcConnellResult::Compute(
        ProteinConformation& conf) {
    OperationLog::Scope scope("MopacMcConnellResult::Compute",
        "compatibility wrapper; atoms=" + std::to_string(conf.AtomCount()) +
        " bonds=" + std::to_string(conf.ProteinRef().BondCount()));

    if (!conf.HasResult<McConnellResult>()) {
        OperationLog::Warn(
            "MopacMcConnellResult::Compute",
            "McConnellResult is not attached; BO-channel legacy fields may "
            "be unset. AttachResult dependency checks will reject this result "
            "in production.");
    }

    auto result_ptr = std::make_unique<MopacMcConnellResult>();
    result_ptr->conf_ = &conf;
    OperationLog::Info(LogCalcMcConnell, "MopacMcConnellResult::Compute",
        "no separate tensor computation; uses McConnellResult BO channel "
        "(Wiberg-weighted D(r)Qhat)");
    return result_ptr;
}


int MopacMcConnellResult::WriteFeatures(const ProteinConformation& conf,
                                         const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    constexpr size_t kTensorCols = 9;
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    std::vector<double> co_sum(N);
    std::vector<double> cn_sum(N);
    std::vector<double> sidechain_sum(N);
    std::vector<double> aromatic_sum(N);
    std::vector<double> co_nearest(N);
    std::vector<double> nearest_co_dist(N);
    std::vector<double> nearest_cn_dist(N);
    std::vector<double> nearest_co_T2(N * kTensorCols, kNaN);
    std::vector<double> nearest_cn_T2(N * kTensorCols, kNaN);
    std::vector<double> backbone_total(N * kTensorCols);
    std::vector<double> sidechain_total(N * kTensorCols);
    std::vector<double> aromatic_total(N * kTensorCols);
    std::vector<double> shielding(N * kTensorCols);

    for (size_t i = 0; i < N; ++i) {
        const auto& atom = conf.AtomAt(i);
        co_sum[i] = atom.mopac_mc_co_sum;
        cn_sum[i] = atom.mopac_mc_cn_sum;
        sidechain_sum[i] = atom.mopac_mc_sidechain_sum;
        aromatic_sum[i] = atom.mopac_mc_aromatic_sum;
        co_nearest[i] = atom.mopac_mc_co_nearest;
        nearest_co_dist[i] = atom.mopac_mc_nearest_CO_dist;
        nearest_cn_dist[i] = atom.mopac_mc_nearest_CN_dist;
        if (atom.mopac_mc_nearest_CO_dist < NO_DATA_SENTINEL) {
            atom.mopac_mc_T2_CO_nearest.PackFull9(
                &nearest_co_T2[i * kTensorCols]);
        }
        if (atom.mopac_mc_nearest_CN_dist < NO_DATA_SENTINEL) {
            atom.mopac_mc_T2_CN_nearest.PackFull9(
                &nearest_cn_T2[i * kTensorCols]);
        }
        atom.mopac_mc_T2_backbone_total.PackFull9(
            &backbone_total[i * kTensorCols]);
        atom.mopac_mc_T2_sidechain_total.PackFull9(
            &sidechain_total[i * kTensorCols]);
        atom.mopac_mc_T2_aromatic_total.PackFull9(
            &aromatic_total[i * kTensorCols]);
        atom.mopac_mc_shielding_contribution.PackFull9(
            &shielding[i * kTensorCols]);
    }

    int files_written = 0;
    auto record_write = [&](bool success, const std::string& filename) {
        if (success) {
            ++files_written;
            return;
        }
        OperationLog::Error(
            "MopacMcConnellResult::WriteFeatures",
            "failed to write " + output_dir + "/" + filename);
    };
    auto write_scalar = [&](const std::string& stem,
                            const std::vector<double>& values) {
        record_write(
            NpyWriter::WriteFloat64(
                output_dir + "/" + stem + ".npy", values.data(), N),
            stem + ".npy");
    };
    auto write_tensor = [&](const std::string& stem,
                            const std::vector<double>& values) {
        record_write(
            NpyWriter::WriteFloat64(
                output_dir + "/" + stem + ".npy", values.data(),
                N, kTensorCols),
            stem + ".npy");
    };

    write_scalar("mopac_mc_co_sum", co_sum);
    write_scalar("mopac_mc_cn_sum", cn_sum);
    write_scalar("mopac_mc_sidechain_sum", sidechain_sum);
    write_scalar("mopac_mc_aromatic_sum", aromatic_sum);
    write_scalar("mopac_mc_co_nearest", co_nearest);
    write_scalar("mopac_mc_nearest_co_dist", nearest_co_dist);
    write_scalar("mopac_mc_nearest_cn_dist", nearest_cn_dist);
    write_tensor("mopac_mc_nearest_co_T2", nearest_co_T2);
    write_tensor("mopac_mc_nearest_cn_T2", nearest_cn_T2);
    write_tensor("mopac_mc_backbone_total", backbone_total);
    write_tensor("mopac_mc_sidechain_total", sidechain_total);
    write_tensor("mopac_mc_aromatic_total", aromatic_total);
    write_tensor("mopac_mc_shielding", shielding);
    return files_written;
}

}  // namespace nmr
