#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "Protein.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <filesystem>
#include <cstdint>

namespace fs = std::filesystem;

namespace nmr {

int ConformationResult::WriteAllFeatures(const ProteinConformation& conf,
                                          const std::string& output_dir) {
    OperationLog::Scope scope("ConformationResult::WriteAllFeatures",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " results=" + std::to_string(conf.AllResults().size()) +
        " dir=" + output_dir);

    fs::create_directories(output_dir);

    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();
    int total = 0;

    // Identity arrays — these belong to no single result.
    {
        std::vector<double> pos(N * 3);
        std::vector<int32_t> elem(N), res_idx(N), res_type(N);
        for (size_t i = 0; i < N; ++i) {
            Vec3 p = conf.PositionAt(i);
            pos[i*3+0] = p.x(); pos[i*3+1] = p.y(); pos[i*3+2] = p.z();
            const Atom& a = protein.AtomAt(i);
            elem[i] = AtomicNumberForElement(a.element);
            res_idx[i] = static_cast<int32_t>(a.residue_index);
            res_type[i] = static_cast<int32_t>(protein.ResidueAt(a.residue_index).type);
        }
        NpyWriter::WriteFloat64(output_dir + "/pos.npy", pos.data(), N, 3);
        NpyWriter::WriteInt32(output_dir + "/element.npy", elem.data(), N);
        NpyWriter::WriteInt32(output_dir + "/residue_index.npy", res_idx.data(), N);
        NpyWriter::WriteInt32(output_dir + "/residue_type.npy", res_type.data(), N);
        total += 4;
    }

    // Walk the conformation's accumulated results. Each one writes its own.
    for (const auto& [tid, result] : conf.AllResults()) {
        int n = result->WriteFeatures(conf, output_dir);
        if (n > 0) {
            OperationLog::Info(LogCalcOther, "WriteAllFeatures",
                result->Name() + " wrote " + std::to_string(n) + " arrays");
        }
        total += n;
    }

    OperationLog::Info(LogCalcOther, "WriteAllFeatures",
        "total: " + std::to_string(total) + " arrays");
    return total;
}

}  // namespace nmr
