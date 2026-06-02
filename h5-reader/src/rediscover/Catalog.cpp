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

#include <cmath>
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
    // MOPAC per-frame charge is NOT a trajectory TR (only the Welford rollup is
    // on the dense H5). The Welford MEAN is the only MOPAC charge available here,
    // so it is a STATIC (atom-axis, no T) charge source — say so honestly.
    add(specs_, ArrayId::MopacChargeWelfordMean, QStringLiteral("mopac_charge_welford_mean"),
        ArrayRank::Scalar, axes(true, false, false, 0), ArrayResidence::StaticTopol,
        QStringLiteral("e"), h5 && h5->mopacChargeWelford());
    // MOPAC-Coulomb-EFG-DERIVED shielding (the moderate Stage-1 field/EFG leg) —
    // a contracted shielding T2, NOT the raw MOPAC Coulomb EFG tensor (that is a
    // per-atom NPY only, not on this trajectory substrate).
    add(specs_, ArrayId::MopacCoulombShielding, QStringLiteral("mopac_coulomb_shielding"),
        ArrayRank::T2_5, axes(true, true, false, 5), ArrayResidence::DenseH5,
        QStringLiteral("ppm"), h5 && h5->mopacCoulombShielding());
    // MOPAC-charge McConnell bond-anisotropy shielding (full shielding tensor; we
    // read its T2 leg, consistent with the e3nn T2 substrate).
    add(specs_, ArrayId::MopacMcShielding, QStringLiteral("mopac_mc_shielding"),
        ArrayRank::T2_5, axes(true, true, false, 5), ArrayResidence::DenseH5,
        QStringLiteral("ppm"), h5 && h5->mopacMcShielding());
    // The FullFat charge-source reconciliation probe output: cosine similarity of
    // the MOPAC-EFG-T2 vs FF14SB-EFG-T2 vectors. Diagnostic/QC, not a feature.
    add(specs_, ArrayId::MopacVsFf14sbReconciliation,
        QStringLiteral("mopac_vs_ff14sb_reconciliation"), ArrayRank::Scalar,
        axes(true, true, false, 0), ArrayResidence::DenseH5, QString(),
        h5 && h5->mopacVsFf14sbReconciliation());
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
    case ArrayId::MopacChargeWelfordMean: {
        const auto* welford = body.run.h5() ? body.run.h5()->mopacChargeWelford() : nullptr;
        return welford
               && atom < welford->n_atoms
               && atom < welford->value.size()
               && atom < welford->n_frames_per_atom.size()
               && welford->n_frames_per_atom[atom] > 0
               && std::isfinite(welford->value[atom].mean);
    }
    case ArrayId::MopacCoulombShielding:
        return body.run.h5() && body.run.h5()->mopacCoulombShielding()
               && atom < body.run.h5()->mopacCoulombShielding()->n_atoms
               && frame < body.run.h5()->mopacCoulombShielding()->n_frames
               && body.run.h5()->mopacCoulombShielding()->sourceAttachedAt(frame);
    case ArrayId::MopacMcShielding:
        return body.run.h5() && body.run.h5()->mopacMcShielding()
               && atom < body.run.h5()->mopacMcShielding()->n_atoms
               && frame < body.run.h5()->mopacMcShielding()->n_frames
               && body.run.h5()->mopacMcShielding()->sourceAttachedAt(frame);
    case ArrayId::MopacVsFf14sbReconciliation:
        return body.run.h5() && body.run.h5()->mopacVsFf14sbReconciliation()
               && atom < body.run.h5()->mopacVsFf14sbReconciliation()->n_atoms
               && frame < body.run.h5()->mopacVsFf14sbReconciliation()->n_frames
               && body.run.h5()->mopacVsFf14sbReconciliation()->sourceAttachedAt(frame);
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
    case ArrayId::MopacChargeWelfordMean:
        return present(body, id, atom, frame) && h5 && h5->mopacChargeWelford()
                   ? h5->mopacChargeWelford()->value[atom].mean
                   : 0.0;
    case ArrayId::MopacVsFf14sbReconciliation:
        return h5 && h5->mopacVsFf14sbReconciliation()
                   ? h5->mopacVsFf14sbReconciliation()->at(atom, frame)
                   : 0.0;
    case ArrayId::MopacCoulombShielding:
    case ArrayId::MopacMcShielding:
        return comp >= 0 && comp < 5 ? valueT2(body, id, atom, frame)[static_cast<std::size_t>(comp)]
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
    // MOPAC-Coulomb-EFG-DERIVED shielding is a T2-only TR (read directly).
    if (id == ArrayId::MopacCoulombShielding && h5 && h5->mopacCoulombShielding())
        return h5->mopacCoulombShielding()->at(atom, frame);
    // MOPAC-McConnell shielding is a full shielding tensor; project its T2 leg.
    if (id == ArrayId::MopacMcShielding && h5 && h5->mopacMcShielding())
        return h5->mopacMcShielding()->at(atom, frame).T2;
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
