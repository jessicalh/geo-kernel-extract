#include "Catalog.h"

#include "AnalysisBody.h"
#include "RunData.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/Conformation.h"
#include "../model/DftShielding.h"
#include "../model/QtAimnet2Group.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtSpecialBuffers.h"
#include "../model/QtTimeSeriesBuffers.h"

#include <stdexcept>

namespace h5reader::rediscover {

namespace {

int ord(ArrayId id) { return static_cast<int>(id); }

AxisSpec axes(bool atom, bool frame, bool slot, int comp) {
    AxisSpec a;
    a.atom = atom;
    a.frame = frame;
    a.slot = slot;
    a.comp = comp > 0;
    a.comp_count = comp;
    return a;
}

void add(std::vector<ArraySpec>& specs, ArrayId id, const QString& name, ArrayRank rank,
         AxisSpec ax, ArrayResidence residence, const QString& unit, bool available,
         ArrayDType dtype = ArrayDType::F64) {
    const int idx = ord(id);
    if (idx >= static_cast<int>(specs.size())) specs.resize(static_cast<std::size_t>(idx + 1));
    specs[static_cast<std::size_t>(idx)] = {id, name, rank, ax, dtype, residence, unit, available};
}

double tensorComponent(const SphericalTensor& st, int comp) {
    if (comp == 0) return st.T0;
    if (comp >= 1 && comp <= 3) return st.T1[static_cast<std::size_t>(comp - 1)];
    if (comp >= 4 && comp <= 8) return st.T2[static_cast<std::size_t>(comp - 4)];
    return 0.0;
}

double matComponent(const Mat3& m, int comp) {
    if (comp < 0 || comp >= 9) return 0.0;
    return m(comp / 3, comp % 3);
}

const model::DftAtomShielding* dftAt(const Body& body, std::size_t atom, std::size_t frame) {
    if (frame >= body.run.frameMap.frameCount()) return nullptr;
    return body.run.dft.AtomShielding(atom, body.run.frameMap.originalIndex(frame));
}

}  // namespace

Catalog::Catalog(const RunData& run) {
    const io::QtTrajectoryH5* h5 = run.h5();
    const bool ff14 = run.protein && !run.protein->atoms().empty()
                      && run.protein->atoms().front().hasPartialCharge;
    add(specs_, ArrayId::Positions, QStringLiteral("positions"), ArrayRank::Vec3,
        axes(true, true, false, 3), ArrayResidence::DenseH5, QStringLiteral("Angstrom"),
        run.conformation != nullptr);
    add(specs_, ArrayId::KernelBs, QStringLiteral("kernel_bs"), ArrayRank::Tensor9,
        axes(true, true, false, 9), ArrayResidence::DenseH5, QStringLiteral("ppm"),
        h5 && h5->bsShielding());
    add(specs_, ArrayId::KernelMc, QStringLiteral("kernel_mc"), ArrayRank::Tensor9,
        axes(true, true, false, 9), ArrayResidence::DenseH5, QStringLiteral("ppm"),
        h5 && h5->mcShielding());
    add(specs_, ArrayId::RingNeighbourhood, QStringLiteral("ring_neighbourhood"), ArrayRank::RingNbhd4,
        axes(true, true, true, 4), ArrayResidence::DenseH5,
        QStringLiteral("Angstrom,Angstrom,Angstrom,radians"), h5 && h5->ringNeighbourhood());
    add(specs_, ArrayId::ApbsEfg, QStringLiteral("apbs_efg"), ArrayRank::T2_5,
        axes(true, true, false, 5), ArrayResidence::DenseH5, QStringLiteral("V/Angstrom^2"),
        h5 && h5->apbsEfg());
    add(specs_, ArrayId::ApbsEfield, QStringLiteral("apbs_efield"), ArrayRank::Vec3,
        axes(true, true, false, 3), ArrayResidence::DenseH5, QStringLiteral("V/Angstrom"),
        h5 && h5->apbsEfield());
    add(specs_, ArrayId::Aimnet2Charge, QStringLiteral("aimnet2_charge"), ArrayRank::Scalar,
        axes(true, true, false, 0), ArrayResidence::DenseH5, QStringLiteral("e"),
        h5 && h5->aimnet2Charge());
    add(specs_, ArrayId::Aimnet2ChargeRespScalar, QStringLiteral("aimnet2_charge_response_gradient_scalar"),
        ArrayRank::Scalar, axes(true, true, false, 0), ArrayResidence::DenseH5, QString(),
        h5 && h5->aimnet2ChargeResponseGradient());
    add(specs_, ArrayId::Aimnet2ChargeRespVector, QStringLiteral("aimnet2_charge_response_gradient_vector"),
        ArrayRank::Vec3, axes(true, true, false, 3), ArrayResidence::DenseH5, QString(),
        h5 && h5->aimnet2ChargeResponseGradient());
    add(specs_, ArrayId::Aimnet2Embedding, QStringLiteral("aimnet2_embedding"), ArrayRank::Embedding256,
        axes(true, true, false, 256), ArrayResidence::DenseH5, QString(),
        h5 && h5->aimnet2Embedding(), ArrayDType::F32);
    add(specs_, ArrayId::Ff14sbCharge, QStringLiteral("ff14sb_charge"), ArrayRank::Scalar,
        axes(true, false, false, 0), ArrayResidence::StaticTopol, QStringLiteral("e"), ff14);
    add(specs_, ArrayId::MopacCharge, QStringLiteral("mopac_charge"), ArrayRank::Scalar,
        axes(true, true, false, 0), ArrayResidence::Absent, QStringLiteral("e"), false);
    add(specs_, ArrayId::DftTotalRaw, QStringLiteral("dft_total_raw"), ArrayRank::Tensor9,
        axes(true, true, false, 9), ArrayResidence::SparseDftByOriginal, QStringLiteral("ppm"),
        run.dft.frameCount() > 0);
    add(specs_, ArrayId::DftDiaRaw, QStringLiteral("dft_dia_raw"), ArrayRank::Tensor9,
        axes(true, true, false, 9), ArrayResidence::SparseDftByOriginal, QStringLiteral("ppm"),
        run.dft.frameCount() > 0);
    add(specs_, ArrayId::DftParaRaw, QStringLiteral("dft_para_raw"), ArrayRank::Tensor9,
        axes(true, true, false, 9), ArrayResidence::SparseDftByOriginal, QStringLiteral("ppm"),
        run.dft.frameCount() > 0);
}

const ArraySpec& Catalog::spec(ArrayId id) const {
    const int idx = ord(id);
    if (idx < 0 || idx >= static_cast<int>(specs_.size())) throw std::out_of_range("Catalog::spec");
    return specs_[static_cast<std::size_t>(idx)];
}

bool Catalog::has(ArrayId id) const { return spec(id).available; }

bool Catalog::present(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) const {
    if (!has(id)) return false;
    switch (id) {
    case ArrayId::DftTotalRaw:
    case ArrayId::DftDiaRaw:
    case ArrayId::DftParaRaw:
        return dftAt(body, atom, frame) != nullptr;
    case ArrayId::Ff14sbCharge:
        return body.run.protein && atom < body.run.protein->atomCount()
               && body.run.protein->atom(atom).hasPartialCharge;
    case ArrayId::Aimnet2Charge:
        return body.run.h5() && body.run.h5()->aimnet2Charge()
               && atom < body.run.h5()->aimnet2Charge()->n_atoms
               && frame < body.run.h5()->aimnet2Charge()->n_frames
               && body.run.h5()->aimnet2Charge()->sourceAttachedAt(frame);
    case ArrayId::Aimnet2ChargeRespScalar:
    case ArrayId::Aimnet2ChargeRespVector:
        return body.run.h5() && body.run.h5()->aimnet2ChargeResponseGradient()
               && atom < body.run.h5()->aimnet2ChargeResponseGradient()->n_atoms
               && frame < body.run.h5()->aimnet2ChargeResponseGradient()->n_frames
               && body.run.h5()->aimnet2ChargeResponseGradient()->meta.sourceAttachedAt(frame);
    case ArrayId::Aimnet2Embedding:
        return body.run.h5() && body.run.h5()->aimnet2Embedding()
               && atom < body.run.h5()->aimnet2Embedding()->n_atoms
               && frame < body.run.h5()->aimnet2Embedding()->n_frames
               && body.run.h5()->aimnet2Embedding()->meta.sourceAttachedAt(frame);
    case ArrayId::MopacCharge:
        return false;
    default:
        return true;
    }
}

double Catalog::value(const Body& body, ArrayId id, std::size_t atom, std::size_t frame,
                      int slot, int comp) const {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    switch (id) {
    case ArrayId::Positions: {
        const Vec3 v = valueVec3(body, id, atom, frame);
        return comp >= 0 && comp < 3 ? v[comp] : 0.0;
    }
    case ArrayId::KernelBs:
        return h5 && h5->bsShielding() ? tensorComponent(h5->bsShielding()->at(atom, frame), comp) : 0.0;
    case ArrayId::KernelMc:
        return h5 && h5->mcShielding() ? tensorComponent(h5->mcShielding()->at(atom, frame), comp) : 0.0;
    case ArrayId::RingNeighbourhood:
        if (h5 && h5->ringNeighbourhood() && slot >= 0 && comp >= 0 && comp < 4)
            return h5->ringNeighbourhood()->at(atom, frame, static_cast<std::size_t>(slot))[static_cast<std::size_t>(comp)];
        return 0.0;
    case ArrayId::ApbsEfg:
        return comp >= 0 && comp < 5 ? valueT2(body, id, atom, frame)[static_cast<std::size_t>(comp)] : 0.0;
    case ArrayId::ApbsEfield:
    case ArrayId::Aimnet2ChargeRespVector: {
        const Vec3 v = valueVec3(body, id, atom, frame);
        return comp >= 0 && comp < 3 ? v[comp] : 0.0;
    }
    case ArrayId::Aimnet2Charge:
        return h5 && h5->aimnet2Charge() ? h5->aimnet2Charge()->at(atom, frame) : 0.0;
    case ArrayId::Aimnet2ChargeRespScalar:
        return h5 && h5->aimnet2ChargeResponseGradient()
                   ? h5->aimnet2ChargeResponseGradient()->scalarAt(atom, frame)
                   : 0.0;
    case ArrayId::Ff14sbCharge:
        return body.run.protein && atom < body.run.protein->atomCount()
                   ? body.run.protein->atom(atom).partialCharge
                   : 0.0;
    case ArrayId::DftTotalRaw:
        return dftAt(body, atom, frame) ? matComponent(dftAt(body, atom, frame)->total_raw, comp) : 0.0;
    case ArrayId::DftDiaRaw:
        return dftAt(body, atom, frame) ? matComponent(dftAt(body, atom, frame)->dia_raw, comp) : 0.0;
    case ArrayId::DftParaRaw:
        return dftAt(body, atom, frame) ? matComponent(dftAt(body, atom, frame)->para_raw, comp) : 0.0;
    case ArrayId::Aimnet2Embedding:
    case ArrayId::MopacCharge:
        return 0.0;
    }
    return 0.0;
}

Vec3 Catalog::valueVec3(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) const {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    switch (id) {
    case ArrayId::Positions:
        return body.run.conformation ? body.run.conformation->atomPosition(frame, atom) : Vec3::Zero();
    case ArrayId::ApbsEfield:
        return h5 && h5->apbsEfield() ? h5->apbsEfield()->at(atom, frame) : Vec3::Zero();
    case ArrayId::Aimnet2ChargeRespVector:
        return h5 && h5->aimnet2ChargeResponseGradient()
                   ? h5->aimnet2ChargeResponseGradient()->vecAt(atom, frame)
                   : Vec3::Zero();
    default:
        return Vec3::Zero();
    }
}

std::array<double, 5> Catalog::valueT2(const Body& body, ArrayId id, std::size_t atom,
                                       std::size_t frame) const {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (id == ArrayId::ApbsEfg && h5 && h5->apbsEfg()) return h5->apbsEfg()->at(atom, frame);
    return {};
}

SphericalTensor Catalog::valueTensor(const Body& body, ArrayId id, std::size_t atom,
                                     std::size_t frame) const {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (id == ArrayId::KernelBs && h5 && h5->bsShielding()) return h5->bsShielding()->at(atom, frame);
    if (id == ArrayId::KernelMc && h5 && h5->mcShielding()) return h5->mcShielding()->at(atom, frame);
    return {};
}

const float* Catalog::valueEmbedding(const Body& body, ArrayId id, std::size_t atom,
                                     std::size_t frame, std::size_t& n_dims_out) const {
    n_dims_out = 0;
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (id != ArrayId::Aimnet2Embedding || !h5 || !h5->aimnet2Embedding()) return nullptr;
    n_dims_out = h5->aimnet2Embedding()->n_dims;
    return h5->aimnet2Embedding()->dataAt(atom, frame);
}

}  // namespace h5reader::rediscover
