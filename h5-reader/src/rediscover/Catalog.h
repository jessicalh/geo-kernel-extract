// Catalog — one typed entry per resident array and a uniform value() surface.

#pragma once

#include "../model/Types.h"

#include <QString>

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::rediscover {

class RunData;
struct Body;

using model::Mat3;
using model::SphericalTensor;
using model::Vec3;

enum class ArrayId : int {
    Positions = 0,
    KernelBs,
    KernelMc,
    RingNeighbourhood,
    ApbsEfg,
    ApbsEfield,
    Aimnet2Charge,
    Aimnet2ChargeRespScalar,
    Aimnet2ChargeRespVector,
    Aimnet2Embedding,
    Ff14sbCharge,
    MopacCharge,
    DftTotalRaw,
    DftDiaRaw,
    DftParaRaw,
};

enum class ArrayRank : int { Scalar, Vec3, T2_5, Tensor9, Embedding256, RingNbhd4 };
enum class ArrayDType : int { F64, F32, I32 };
enum class ArrayResidence : int { DenseH5, StaticTopol, SparseDftByOriginal, Absent };

struct AxisSpec {
    bool atom = false;
    bool frame = false;
    bool slot = false;
    bool comp = false;
    int comp_count = 0;
};

struct ArraySpec {
    ArrayId id;
    QString name;
    ArrayRank rank = ArrayRank::Scalar;
    AxisSpec axes;
    ArrayDType dtype = ArrayDType::F64;
    ArrayResidence residence = ArrayResidence::DenseH5;
    QString unit;
    bool available = false;
};

class Catalog {
public:
    explicit Catalog(const RunData& run);

    const ArraySpec& spec(ArrayId id) const;
    bool has(ArrayId id) const;
    bool present(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) const;

    double value(const Body& body, ArrayId id, std::size_t atom, std::size_t frame,
                 int slot = -1, int comp = -1) const;
    Vec3 valueVec3(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) const;
    std::array<double, 5> valueT2(const Body& body, ArrayId id, std::size_t atom,
                                  std::size_t frame) const;
    SphericalTensor valueTensor(const Body& body, ArrayId id, std::size_t atom,
                                std::size_t frame) const;
    const float* valueEmbedding(const Body& body, ArrayId id, std::size_t atom,
                                std::size_t frame, std::size_t& n_dims_out) const;

private:
    std::vector<ArraySpec> specs_;
};

}  // namespace h5reader::rediscover
