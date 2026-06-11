// Catalog — one typed entry per resident array and a uniform value() surface.

#pragma once

#include "../io/QtFieldCatalog.gen.h"
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
    MopacCharge,                 // producer mopac_charges.npy (static or flattened per-frame)
    MopacChargeWelfordMean,      // per-atom MOPAC charge Welford MEAN (static, no T)
    MopacCoulombShielding,       // T2: MOPAC-Coulomb-EFG-DERIVED shielding (NOT raw EFG)
    MopacCoulombEfield,
    MopacMcShielding,            // T2: MOPAC-charge McConnell bond-anisotropy shielding
    MopacVsFf14sbReconciliation, // scalar QC: cos(MOPAC-EFG-T2, FF14SB-EFG-T2)
    HbondScalars,
    HbondNearestDirection,
    HbondFlags,
    HbondCount,
    HmShielding,
    BSPerTypeT0,
    BSPerTypeT1,
    BSPerTypeT2,
    HMPerTypeT0,
    HMPerTypeT1,
    HMPerTypeT2,
    LarsenHBondShielding,
    WaterEfield,
    WaterNFirst,
    WaterNSecond,
    HydrationHalfShellAsymmetry,
    HydrationDipoleCos,
    Sasa,
    SasaNormal,
    EeqChargeMean,
    EeqCoordinationNumber,
    Aimnet2Efg,
    DsspBackbone,
    DsspChi,
    DsspSs8,
    DsspHBondEnergy,
    PuckerQ,
    PuckerTheta,
    OmegaActual,
    AromaticChi2,
    Pyramidalization,
    DftTotalRaw,
    DftDiaRaw,
    DftParaRaw,
    McPeptideCoFixed,
    McPeptideCoBo,
    McPeptideCoRhombic,
    McPeptideCnFixed,
    McPeptideCnBo,
    McBackboneOtherFixed,
    McBackboneOtherBo,
    McSidechainCoFixed,
    McSidechainCoBo,
    McSidechainOtherFixed,
    McSidechainOtherBo,
    McDisulfideFixed,
    McDisulfideBo,
    McAromaticZeroedFixed,
    McAromaticZeroedBo,
    McNearestCoT2,
    McNearestCnT2,
};

enum class ArrayRank : int {
    Scalar,
    Vec3,
    T2_5,
    Tensor9,
    PerTypeT0_8,
    PerTypeT1_24,
    PerTypeT2_40,
    Embedding256,
    RingNbhd4
};
enum class ArrayDType : int { F64, F32, I32 };
enum class ArrayResidence : int { DenseH5, StaticTopol, StaticNpy, SparseDftByOriginal, Absent };

enum class FieldProvider : int {
    StaticNpy,
    TrajectoryNpy,
    DenseH5,
    SparseDftByOriginal,
    TypedTopology,
    DatasetAbsent,
    Unsupported
};

QString FieldProviderName(FieldProvider provider);

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

std::optional<io::FieldKind> ProducerFieldFor(ArrayId id);
std::optional<ArrayId> ArrayIdForProducerField(io::FieldKind kind);
const io::FieldSpec* ProducerFieldSpecFor(ArrayId id);

class Catalog {
public:
    explicit Catalog(const RunData& run);

    const ArraySpec& spec(ArrayId id) const;
    bool has(ArrayId id) const;
    bool present(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) const;

    FieldProvider provider(const RunData& run, io::FieldKind kind,
                           QString* reason_out = nullptr) const;
    bool has(const RunData& run, io::FieldKind kind) const;
    bool present(const Body& body, io::FieldKind kind, std::size_t nativeRow,
                 std::size_t frame, int component = -1,
                 QString* reason_out = nullptr) const;
    std::optional<double> value(const Body& body, io::FieldKind kind,
                                std::size_t nativeRow, std::size_t frame,
                                int component, QString* reason_out = nullptr) const;
    std::size_t nativeRowCount(const Body& body, io::FieldKind kind,
                               std::size_t frame = 0) const;
    std::size_t componentCount(const Body& body, io::FieldKind kind,
                               std::size_t frame = 0) const;
    std::size_t frameCount(const Body& body, io::FieldKind kind) const;

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
