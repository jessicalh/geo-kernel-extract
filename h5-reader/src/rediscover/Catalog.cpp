#include "Catalog.h"

#include "AnalysisBody.h"
#include "RunData.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/Conformation.h"
#include "../model/DftShielding.h"
#include "../model/QtAimnet2Group.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResultBlocks.h"
#include "../model/QtSpecialBuffers.h"
#include "../model/QtTimeSeriesBuffers.h"

#include <cmath>
#include <stdexcept>

namespace h5reader::rediscover {

std::optional<io::FieldKind> ProducerFieldFor(ArrayId id) {
    using io::FieldKind;
    switch (id) {
    case ArrayId::Positions: return FieldKind::Pos;
    case ArrayId::KernelBs: return FieldKind::BSShielding;
    case ArrayId::KernelMc: return FieldKind::McShielding;
    case ArrayId::ApbsEfg: return FieldKind::APBSEFG;
    case ArrayId::ApbsEfield: return FieldKind::APBSE;
    case ArrayId::Aimnet2Charge: return FieldKind::AIMNet2Charges;
    case ArrayId::Aimnet2ChargeRespScalar: return FieldKind::AIMNet2ChargeResponseGradientScalar;
    case ArrayId::Aimnet2ChargeRespVector: return FieldKind::AIMNet2ChargeResponseGradient;
    case ArrayId::Aimnet2Embedding: return FieldKind::AIMNet2Aim;
    case ArrayId::Ff14sbCharge: return FieldKind::FfPartialCharge;
    case ArrayId::MopacCharge: return FieldKind::MOPACCharges;
    case ArrayId::MopacCoulombShielding: return FieldKind::MOPACCoulombEFG;
    case ArrayId::MopacCoulombEfield: return FieldKind::MOPACCoulombE;
    case ArrayId::MopacMcShielding: return FieldKind::MOPACMcShielding;
    case ArrayId::HbondScalars: return FieldKind::HBondScalars;
    case ArrayId::HbondNearestDirection: return FieldKind::HBondNearestDir;
    case ArrayId::HbondFlags: return FieldKind::HBondFlags;
    case ArrayId::HbondCount: return FieldKind::LarsenHBondCount;
    case ArrayId::HmShielding: return FieldKind::HMShielding;
    case ArrayId::BSPerTypeT0: return FieldKind::BSPerTypeT0;
    case ArrayId::BSPerTypeT1: return FieldKind::BSPerTypeT1;
    case ArrayId::BSPerTypeT2: return FieldKind::BSPerTypeT2;
    case ArrayId::HMPerTypeT0: return FieldKind::HMPerTypeT0;
    case ArrayId::HMPerTypeT1: return FieldKind::HMPerTypeT1;
    case ArrayId::HMPerTypeT2: return FieldKind::HMPerTypeT2;
    case ArrayId::LarsenHBondShielding: return FieldKind::LarsenHBondShielding;
    case ArrayId::Sasa: return FieldKind::AtomSASA;
    case ArrayId::SasaNormal: return FieldKind::SASANormal;
    case ArrayId::EeqChargeMean: return FieldKind::EEQCharges;
    case ArrayId::EeqCoordinationNumber: return FieldKind::EEQCN;
    case ArrayId::Aimnet2Efg: return FieldKind::AIMNet2EFG;
    case ArrayId::DsspBackbone: return FieldKind::DSSPBackbone;
    case ArrayId::DsspChi: return FieldKind::DSSPChi;
    case ArrayId::DsspSs8: return FieldKind::DSSPSs8;
    case ArrayId::DsspHBondEnergy: return FieldKind::DSSPHBondEnergy;
    case ArrayId::PuckerQ: return FieldKind::PuckerQ;
    case ArrayId::PuckerTheta: return FieldKind::PuckerTheta;
    case ArrayId::OmegaActual: return FieldKind::OmegaActual;
    case ArrayId::AromaticChi2: return FieldKind::AromaticChi2;
    case ArrayId::Pyramidalization: return FieldKind::Pyramidalization;
    case ArrayId::DftTotalRaw: return FieldKind::OrcaTotal;
    case ArrayId::DftDiaRaw: return FieldKind::OrcaDiamagnetic;
    case ArrayId::DftParaRaw: return FieldKind::OrcaParamagnetic;
    case ArrayId::McPeptideCoFixed: return FieldKind::McPeptideCoFixed;
    case ArrayId::McPeptideCoBo: return FieldKind::McPeptideCoBo;
    case ArrayId::McPeptideCoRhombic: return FieldKind::McPeptideCoRhombic;
    case ArrayId::McPeptideCnFixed: return FieldKind::McPeptideCNFixed;
    case ArrayId::McPeptideCnBo: return FieldKind::McPeptideCNBo;
    case ArrayId::McBackboneOtherFixed: return FieldKind::McBackboneOtherFixed;
    case ArrayId::McBackboneOtherBo: return FieldKind::McBackboneOtherBo;
    case ArrayId::McSidechainCoFixed: return FieldKind::McSidechainCoFixed;
    case ArrayId::McSidechainCoBo: return FieldKind::McSidechainCoBo;
    case ArrayId::McSidechainOtherFixed: return FieldKind::McSidechainOtherFixed;
    case ArrayId::McSidechainOtherBo: return FieldKind::McSidechainOtherBo;
    case ArrayId::McDisulfideFixed: return FieldKind::McDisulfideFixed;
    case ArrayId::McDisulfideBo: return FieldKind::McDisulfideBo;
    case ArrayId::McAromaticZeroedFixed: return FieldKind::McAromaticZeroedFixed;
    case ArrayId::McAromaticZeroedBo: return FieldKind::McAromaticZeroedBo;
    case ArrayId::McNearestCoT2: return FieldKind::McNearestCoT2;
    case ArrayId::McNearestCnT2: return FieldKind::McNearestCNT2;
    default: return std::nullopt;
    }
}

std::optional<ArrayId> ArrayIdForProducerField(io::FieldKind kind) {
    switch (kind) {
    case io::FieldKind::Pos: return ArrayId::Positions;
    case io::FieldKind::BSShielding: return ArrayId::KernelBs;
    case io::FieldKind::McShielding: return ArrayId::KernelMc;
    case io::FieldKind::APBSEFG: return ArrayId::ApbsEfg;
    case io::FieldKind::APBSE: return ArrayId::ApbsEfield;
    case io::FieldKind::AIMNet2Charges: return ArrayId::Aimnet2Charge;
    case io::FieldKind::AIMNet2ChargeResponseGradientScalar: return ArrayId::Aimnet2ChargeRespScalar;
    case io::FieldKind::AIMNet2ChargeResponseGradient: return ArrayId::Aimnet2ChargeRespVector;
    case io::FieldKind::AIMNet2Aim: return ArrayId::Aimnet2Embedding;
    case io::FieldKind::FfPartialCharge: return ArrayId::Ff14sbCharge;
    case io::FieldKind::MOPACCharges: return ArrayId::MopacCharge;
    case io::FieldKind::MOPACCoulombEFG: return ArrayId::MopacCoulombShielding;
    case io::FieldKind::MOPACCoulombE: return ArrayId::MopacCoulombEfield;
    case io::FieldKind::MOPACMcShielding: return ArrayId::MopacMcShielding;
    case io::FieldKind::HBondScalars: return ArrayId::HbondScalars;
    case io::FieldKind::HBondNearestDir: return ArrayId::HbondNearestDirection;
    case io::FieldKind::HBondFlags: return ArrayId::HbondFlags;
    case io::FieldKind::LarsenHBondCount: return ArrayId::HbondCount;
    case io::FieldKind::HMShielding: return ArrayId::HmShielding;
    case io::FieldKind::BSPerTypeT0: return ArrayId::BSPerTypeT0;
    case io::FieldKind::BSPerTypeT1: return ArrayId::BSPerTypeT1;
    case io::FieldKind::BSPerTypeT2: return ArrayId::BSPerTypeT2;
    case io::FieldKind::HMPerTypeT0: return ArrayId::HMPerTypeT0;
    case io::FieldKind::HMPerTypeT1: return ArrayId::HMPerTypeT1;
    case io::FieldKind::HMPerTypeT2: return ArrayId::HMPerTypeT2;
    case io::FieldKind::LarsenHBondShielding: return ArrayId::LarsenHBondShielding;
    case io::FieldKind::AtomSASA: return ArrayId::Sasa;
    case io::FieldKind::SASANormal: return ArrayId::SasaNormal;
    case io::FieldKind::EEQCharges: return ArrayId::EeqChargeMean;
    case io::FieldKind::EEQCN: return ArrayId::EeqCoordinationNumber;
    case io::FieldKind::AIMNet2EFG: return ArrayId::Aimnet2Efg;
    case io::FieldKind::DSSPBackbone: return ArrayId::DsspBackbone;
    case io::FieldKind::DSSPChi: return ArrayId::DsspChi;
    case io::FieldKind::DSSPSs8: return ArrayId::DsspSs8;
    case io::FieldKind::DSSPHBondEnergy: return ArrayId::DsspHBondEnergy;
    case io::FieldKind::PuckerQ: return ArrayId::PuckerQ;
    case io::FieldKind::PuckerTheta: return ArrayId::PuckerTheta;
    case io::FieldKind::OmegaActual: return ArrayId::OmegaActual;
    case io::FieldKind::AromaticChi2: return ArrayId::AromaticChi2;
    case io::FieldKind::Pyramidalization: return ArrayId::Pyramidalization;
    case io::FieldKind::OrcaTotal: return ArrayId::DftTotalRaw;
    case io::FieldKind::OrcaDiamagnetic: return ArrayId::DftDiaRaw;
    case io::FieldKind::OrcaParamagnetic: return ArrayId::DftParaRaw;
    case io::FieldKind::McPeptideCoFixed: return ArrayId::McPeptideCoFixed;
    case io::FieldKind::McPeptideCoBo: return ArrayId::McPeptideCoBo;
    case io::FieldKind::McPeptideCoRhombic: return ArrayId::McPeptideCoRhombic;
    case io::FieldKind::McPeptideCNFixed: return ArrayId::McPeptideCnFixed;
    case io::FieldKind::McPeptideCNBo: return ArrayId::McPeptideCnBo;
    case io::FieldKind::McBackboneOtherFixed: return ArrayId::McBackboneOtherFixed;
    case io::FieldKind::McBackboneOtherBo: return ArrayId::McBackboneOtherBo;
    case io::FieldKind::McSidechainCoFixed: return ArrayId::McSidechainCoFixed;
    case io::FieldKind::McSidechainCoBo: return ArrayId::McSidechainCoBo;
    case io::FieldKind::McSidechainOtherFixed: return ArrayId::McSidechainOtherFixed;
    case io::FieldKind::McSidechainOtherBo: return ArrayId::McSidechainOtherBo;
    case io::FieldKind::McDisulfideFixed: return ArrayId::McDisulfideFixed;
    case io::FieldKind::McDisulfideBo: return ArrayId::McDisulfideBo;
    case io::FieldKind::McAromaticZeroedFixed: return ArrayId::McAromaticZeroedFixed;
    case io::FieldKind::McAromaticZeroedBo: return ArrayId::McAromaticZeroedBo;
    case io::FieldKind::McNearestCoT2: return ArrayId::McNearestCoT2;
    case io::FieldKind::McNearestCNT2: return ArrayId::McNearestCnT2;
    default: return std::nullopt;
    }
}

const io::FieldSpec* ProducerFieldSpecFor(ArrayId id) {
    const std::optional<io::FieldKind> kind = ProducerFieldFor(id);
    return kind ? &io::FieldSpecFor(*kind) : nullptr;
}

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

QString qsv(std::string_view s) {
    return QString::fromUtf8(s.data(), static_cast<qsizetype>(s.size()));
}

ArrayRank rankForProducerField(const io::FieldSpec& spec) {
    if (spec.cols == 256) return ArrayRank::Embedding256;
    if (spec.cols == 40) return ArrayRank::PerTypeT2_40;
    if (spec.cols == 24) return ArrayRank::PerTypeT1_24;
    if (spec.cols == 8) return ArrayRank::PerTypeT0_8;
    if (spec.cols == 9) return ArrayRank::Tensor9;
    if (spec.cols == 5) return ArrayRank::T2_5;
    if (spec.cols == 3) return ArrayRank::Vec3;
    return ArrayRank::Scalar;
}

AxisSpec axesForProducerField(const io::FieldSpec& spec) {
    AxisSpec a;
    a.atom = spec.axis == io::NativeAxis::Atom;
    a.frame = spec.axis != io::NativeAxis::Protein;
    a.comp_count = spec.cols > 0 ? spec.cols : 0;
    a.comp = a.comp_count > 1;
    return a;
}

void addProducer(std::vector<ArraySpec>& specs,
                 ArrayId id,
                 ArrayResidence residence,
                 bool available,
                 ArrayDType dtype = ArrayDType::F64) {
    const io::FieldSpec* field = ProducerFieldSpecFor(id);
    if (!field) return;
    add(specs, id, qsv(field->stem), rankForProducerField(*field),
        axesForProducerField(*field), residence, qsv(field->units), available, dtype);
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

bool shieldingPresent(const model::QtShieldingTimeSeries* ts, std::size_t atom, std::size_t frame) {
    return ts && atom < ts->n_atoms && frame < ts->n_frames && ts->sourceAttachedAt(frame);
}

bool scalarPresent(const model::QtScalarTimeSeries* ts, std::size_t atom, std::size_t frame) {
    return ts && atom < ts->n_atoms && frame < ts->n_frames && ts->sourceAttachedAt(frame);
}

bool scalarWelfordMeanPresent(const model::QtScalarWelford* w, std::size_t atom) {
    return w && atom < w->n_atoms && atom < w->value.size()
           && atom < w->n_frames_per_atom.size()
           && w->n_frames_per_atom[atom] > 0
           && std::isfinite(w->value[atom].mean);
}

const StaticNpyArray* staticAt(const RunData& run, ArrayId id) {
    const std::optional<io::FieldKind> kind = ProducerFieldFor(id);
    return kind ? run.producerArray(*kind) : nullptr;
}

std::size_t staticRow(const StaticNpyArray* a, std::size_t atom, std::size_t frame) {
    return a ? a->rowFor(atom, frame) : 0;
}

bool staticArrayPresent(const StaticNpyArray* a, std::size_t atom, std::size_t frame = 0,
                        int colsNeeded = 1) {
    if (!a || a->cols < static_cast<std::size_t>(colsNeeded)) return false;
    if (a->frameVarying
        && (atom >= a->atomsPerFrame || frame >= a->frameCount)) return false;
    return staticRow(a, atom, frame) < a->rows;
}

double staticValue(const StaticNpyArray* a, std::size_t atom, std::size_t frame = 0, int comp = 0) {
    if (!staticArrayPresent(a, atom, frame, comp + 1)) return 0.0;
    return a->value(staticRow(a, atom, frame), static_cast<std::size_t>(comp));
}

Vec3 staticVec3(const StaticNpyArray* a, std::size_t atom, std::size_t frame = 0) {
    if (!staticArrayPresent(a, atom, frame, 3)) return Vec3::Zero();
    const std::size_t r = staticRow(a, atom, frame);
    return Vec3(a->value(r, 0), a->value(r, 1), a->value(r, 2));
}

std::array<double, 5> staticT2(const StaticNpyArray* a, std::size_t atom, std::size_t frame = 0) {
    std::array<double, 5> out = {};
    if (!a || !staticArrayPresent(a, atom, frame)) return out;
    const std::size_t r = staticRow(a, atom, frame);
    if (a->cols >= 9) {
        for (std::size_t i = 0; i < 5; ++i) out[i] = a->value(r, 4 + i);
    } else if (a->cols >= 5) {
        for (std::size_t i = 0; i < 5; ++i) out[i] = a->value(r, i);
    }
    return out;
}

SphericalTensor staticTensor(const StaticNpyArray* a, std::size_t atom, std::size_t frame = 0) {
    if (!staticArrayPresent(a, atom, frame, 9)) return {};
    return model::UnpackSphericalTensor(a->rowData(staticRow(a, atom, frame)));
}

}  // namespace

Catalog::Catalog(const RunData& run) {
    const io::QtTrajectoryH5* h5 = run.h5();
    auto hasStatic = [&](ArrayId id) { return staticAt(run, id) != nullptr; };
    auto residence = [&](ArrayId id, ArrayResidence h5Residence = ArrayResidence::DenseH5) {
        return hasStatic(id) ? ArrayResidence::StaticNpy : h5Residence;
    };
    const bool ff14 = run.protein && !run.protein->atoms().empty()
                      && (run.protein->atoms().front().hasPartialCharge || hasStatic(ArrayId::Ff14sbCharge));
    addProducer(specs_, ArrayId::Positions, residence(ArrayId::Positions), run.conformation != nullptr);
    addProducer(specs_, ArrayId::KernelBs, residence(ArrayId::KernelBs),
                (h5 && h5->bsShielding()) || hasStatic(ArrayId::KernelBs));
    addProducer(specs_, ArrayId::KernelMc, residence(ArrayId::KernelMc),
                (h5 && h5->mcShielding()) || hasStatic(ArrayId::KernelMc));
    add(specs_, ArrayId::RingNeighbourhood, QStringLiteral("ring_neighbourhood"), ArrayRank::RingNbhd4,
        axes(true, true, true, 4), ArrayResidence::DenseH5,
        QStringLiteral("Angstrom,Angstrom,Angstrom,radians"), h5 && h5->ringNeighbourhood());
    addProducer(specs_, ArrayId::ApbsEfg, residence(ArrayId::ApbsEfg),
                (h5 && h5->apbsEfg()) || hasStatic(ArrayId::ApbsEfg));
    addProducer(specs_, ArrayId::ApbsEfield, residence(ArrayId::ApbsEfield),
                (h5 && h5->apbsEfield()) || hasStatic(ArrayId::ApbsEfield));
    addProducer(specs_, ArrayId::Aimnet2Charge, residence(ArrayId::Aimnet2Charge),
                (h5 && h5->aimnet2Charge()) || hasStatic(ArrayId::Aimnet2Charge));
    addProducer(specs_, ArrayId::Aimnet2ChargeRespScalar,
                residence(ArrayId::Aimnet2ChargeRespScalar),
                (h5 && h5->aimnet2ChargeResponseGradient()) || hasStatic(ArrayId::Aimnet2ChargeRespScalar));
    addProducer(specs_, ArrayId::Aimnet2ChargeRespVector,
                residence(ArrayId::Aimnet2ChargeRespVector),
                (h5 && h5->aimnet2ChargeResponseGradient()) || hasStatic(ArrayId::Aimnet2ChargeRespVector));
    addProducer(specs_, ArrayId::Aimnet2Embedding, residence(ArrayId::Aimnet2Embedding),
                (h5 && h5->aimnet2Embedding()) || hasStatic(ArrayId::Aimnet2Embedding), ArrayDType::F32);
    addProducer(specs_, ArrayId::Ff14sbCharge,
                hasStatic(ArrayId::Ff14sbCharge) ? ArrayResidence::StaticNpy : ArrayResidence::StaticTopol,
                ff14);
    addProducer(specs_, ArrayId::MopacCharge, residence(ArrayId::MopacCharge, ArrayResidence::Absent),
                hasStatic(ArrayId::MopacCharge));
    // MOPAC charge may be present as the producer's raw mopac_charges.npy
    // (registered above) or as a dense-H5 Welford rollup. Keep the rollup
    // advertised separately for old consumers.
    add(specs_, ArrayId::MopacChargeWelfordMean, QStringLiteral("mopac_charge_welford_mean"),
        ArrayRank::Scalar, axes(true, false, false, 0),
        residence(ArrayId::MopacChargeWelfordMean, ArrayResidence::StaticTopol),
        QStringLiteral("e"), (h5 && h5->mopacChargeWelford()) || hasStatic(ArrayId::MopacChargeWelfordMean));
    // MOPAC-Coulomb-EFG-DERIVED shielding (the moderate Stage-1 field/EFG leg) —
    // a contracted shielding T2, NOT the raw MOPAC Coulomb EFG tensor (that is a
    // per-atom NPY only, not on this trajectory substrate).
    addProducer(specs_, ArrayId::MopacCoulombShielding,
                residence(ArrayId::MopacCoulombShielding),
                (h5 && h5->mopacCoulombShielding()) || hasStatic(ArrayId::MopacCoulombShielding));
    addProducer(specs_, ArrayId::MopacCoulombEfield, residence(ArrayId::MopacCoulombEfield),
                hasStatic(ArrayId::MopacCoulombEfield));
    // MOPAC-charge McConnell bond-anisotropy shielding (full shielding tensor; we
    // read its T2 leg, consistent with the e3nn T2 substrate).
    addProducer(specs_, ArrayId::MopacMcShielding, residence(ArrayId::MopacMcShielding),
                (h5 && h5->mopacMcShielding()) || hasStatic(ArrayId::MopacMcShielding));
    // The FullFat charge-source reconciliation probe output: cosine similarity of
    // the MOPAC-EFG-T2 vs FF14SB-EFG-T2 vectors. Diagnostic/QC, not a feature.
    add(specs_, ArrayId::MopacVsFf14sbReconciliation,
        QStringLiteral("mopac_vs_ff14sb_reconciliation"), ArrayRank::Scalar,
        axes(true, true, false, 0), ArrayResidence::DenseH5, QString(),
        h5 && h5->mopacVsFf14sbReconciliation());
    addProducer(specs_, ArrayId::HbondScalars, residence(ArrayId::HbondScalars),
                hasStatic(ArrayId::HbondScalars));
    addProducer(specs_, ArrayId::HbondNearestDirection,
                residence(ArrayId::HbondNearestDirection),
                hasStatic(ArrayId::HbondNearestDirection));
    addProducer(specs_, ArrayId::HbondFlags, residence(ArrayId::HbondFlags),
                hasStatic(ArrayId::HbondFlags));
    addProducer(specs_, ArrayId::HbondCount, residence(ArrayId::HbondCount),
                (h5 && h5->larsenHBondCount()) || hasStatic(ArrayId::HbondCount));
    addProducer(specs_, ArrayId::HmShielding, residence(ArrayId::HmShielding),
                (h5 && h5->hmShielding()) || hasStatic(ArrayId::HmShielding));
    addProducer(specs_, ArrayId::BSPerTypeT0, residence(ArrayId::BSPerTypeT0),
                hasStatic(ArrayId::BSPerTypeT0));
    addProducer(specs_, ArrayId::BSPerTypeT1, residence(ArrayId::BSPerTypeT1),
                hasStatic(ArrayId::BSPerTypeT1));
    addProducer(specs_, ArrayId::BSPerTypeT2, residence(ArrayId::BSPerTypeT2),
                hasStatic(ArrayId::BSPerTypeT2));
    addProducer(specs_, ArrayId::HMPerTypeT0, residence(ArrayId::HMPerTypeT0),
                hasStatic(ArrayId::HMPerTypeT0));
    addProducer(specs_, ArrayId::HMPerTypeT1, residence(ArrayId::HMPerTypeT1),
                hasStatic(ArrayId::HMPerTypeT1));
    addProducer(specs_, ArrayId::HMPerTypeT2, residence(ArrayId::HMPerTypeT2),
                hasStatic(ArrayId::HMPerTypeT2));
    addProducer(specs_, ArrayId::LarsenHBondShielding, residence(ArrayId::LarsenHBondShielding),
                hasStatic(ArrayId::LarsenHBondShielding));
    add(specs_, ArrayId::WaterEfield, QStringLiteral("water_efield"),
        ArrayRank::Vec3, axes(true, true, false, 3), ArrayResidence::DenseH5,
        QStringLiteral("V/Angstrom"),
        h5 && h5->waterFieldTimeSeries() && h5->waterFieldTimeSeries()->hasEfield());
    add(specs_, ArrayId::WaterNFirst, QStringLiteral("water_n_first"),
        ArrayRank::Scalar, axes(true, true, false, 0), ArrayResidence::DenseH5,
        QStringLiteral("count"),
        h5 && h5->waterFieldTimeSeries() && h5->waterFieldTimeSeries()->hasNFirst());
    add(specs_, ArrayId::WaterNSecond, QStringLiteral("water_n_second"),
        ArrayRank::Scalar, axes(true, true, false, 0), ArrayResidence::DenseH5,
        QStringLiteral("count"),
        h5 && h5->waterFieldTimeSeries() && h5->waterFieldTimeSeries()->hasNSecond());
    add(specs_, ArrayId::HydrationHalfShellAsymmetry,
        QStringLiteral("hydration_half_shell_asymmetry"), ArrayRank::Scalar,
        axes(true, true, false, 0), ArrayResidence::DenseH5, QStringLiteral("fraction"),
        h5 && h5->hydrationShellTimeSeries()
            && h5->hydrationShellTimeSeries()->hasHalfShellAsymmetry());
    add(specs_, ArrayId::HydrationDipoleCos, QStringLiteral("hydration_dipole_cos"),
        ArrayRank::Scalar, axes(true, true, false, 0), ArrayResidence::DenseH5,
        QStringLiteral("cos_angle"),
        h5 && h5->hydrationShellTimeSeries()
            && h5->hydrationShellTimeSeries()->hasMeanWaterDipoleCos());
    addProducer(specs_, ArrayId::Sasa, residence(ArrayId::Sasa),
                (h5 && h5->sasa()) || hasStatic(ArrayId::Sasa));
    addProducer(specs_, ArrayId::SasaNormal, residence(ArrayId::SasaNormal),
                (h5 && h5->hydrationGeometryTimeSeries()
                 && h5->hydrationGeometryTimeSeries()->hasSurfaceNormal())
                    || hasStatic(ArrayId::SasaNormal));
    addProducer(specs_, ArrayId::EeqChargeMean,
                residence(ArrayId::EeqChargeMean, ArrayResidence::StaticTopol),
                (h5 && h5->eeqWelford()) || hasStatic(ArrayId::EeqChargeMean));
    addProducer(specs_, ArrayId::EeqCoordinationNumber,
                residence(ArrayId::EeqCoordinationNumber, ArrayResidence::Absent),
                hasStatic(ArrayId::EeqCoordinationNumber));
    addProducer(specs_, ArrayId::Aimnet2Efg, residence(ArrayId::Aimnet2Efg),
                hasStatic(ArrayId::Aimnet2Efg));
    addProducer(specs_, ArrayId::DsspBackbone, residence(ArrayId::DsspBackbone),
                hasStatic(ArrayId::DsspBackbone));
    addProducer(specs_, ArrayId::DsspChi, residence(ArrayId::DsspChi),
                hasStatic(ArrayId::DsspChi));
    addProducer(specs_, ArrayId::DsspSs8, residence(ArrayId::DsspSs8),
                hasStatic(ArrayId::DsspSs8));
    addProducer(specs_, ArrayId::DsspHBondEnergy, residence(ArrayId::DsspHBondEnergy),
                hasStatic(ArrayId::DsspHBondEnergy));
    addProducer(specs_, ArrayId::PuckerQ, residence(ArrayId::PuckerQ),
                hasStatic(ArrayId::PuckerQ));
    addProducer(specs_, ArrayId::PuckerTheta, residence(ArrayId::PuckerTheta),
                hasStatic(ArrayId::PuckerTheta));
    addProducer(specs_, ArrayId::OmegaActual, residence(ArrayId::OmegaActual),
                hasStatic(ArrayId::OmegaActual));
    addProducer(specs_, ArrayId::AromaticChi2, residence(ArrayId::AromaticChi2),
                hasStatic(ArrayId::AromaticChi2));
    addProducer(specs_, ArrayId::Pyramidalization,
                residence(ArrayId::Pyramidalization, ArrayResidence::Absent),
                hasStatic(ArrayId::Pyramidalization));
    addProducer(specs_, ArrayId::DftTotalRaw, ArrayResidence::SparseDftByOriginal,
                run.dft.frameCount() > 0);
    addProducer(specs_, ArrayId::DftDiaRaw, ArrayResidence::SparseDftByOriginal,
                run.dft.frameCount() > 0);
    addProducer(specs_, ArrayId::DftParaRaw, ArrayResidence::SparseDftByOriginal,
                run.dft.frameCount() > 0);
    auto addMcTensor = [&](ArrayId id) {
        addProducer(specs_, id, residence(id, ArrayResidence::Absent), hasStatic(id));
    };
    addMcTensor(ArrayId::McPeptideCoFixed);
    addMcTensor(ArrayId::McPeptideCoBo);
    addMcTensor(ArrayId::McPeptideCoRhombic);
    addMcTensor(ArrayId::McPeptideCnFixed);
    addMcTensor(ArrayId::McPeptideCnBo);
    addMcTensor(ArrayId::McBackboneOtherFixed);
    addMcTensor(ArrayId::McBackboneOtherBo);
    addMcTensor(ArrayId::McSidechainCoFixed);
    addMcTensor(ArrayId::McSidechainCoBo);
    addMcTensor(ArrayId::McSidechainOtherFixed);
    addMcTensor(ArrayId::McSidechainOtherBo);
    addMcTensor(ArrayId::McDisulfideFixed);
    addMcTensor(ArrayId::McDisulfideBo);
    addMcTensor(ArrayId::McAromaticZeroedFixed);
    addMcTensor(ArrayId::McAromaticZeroedBo);
    addMcTensor(ArrayId::McNearestCoT2);
    addMcTensor(ArrayId::McNearestCnT2);
}

const ArraySpec& Catalog::spec(ArrayId id) const {
    const int idx = ord(id);
    if (idx < 0 || idx >= static_cast<int>(specs_.size())) throw std::out_of_range("Catalog::spec");
    return specs_[static_cast<std::size_t>(idx)];
}

bool Catalog::has(ArrayId id) const { return spec(id).available; }

bool Catalog::present(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) const {
    if (!has(id)) return false;
    if (const StaticNpyArray* a = staticAt(body.run, id)) {
        if (id == ArrayId::Aimnet2Embedding) return atom < a->rows && a->cols == 256 && !a->floatValues.empty();
        return staticArrayPresent(a, atom, frame);
    }
    switch (id) {
    case ArrayId::DftTotalRaw:
    case ArrayId::DftDiaRaw:
    case ArrayId::DftParaRaw:
        return dftAt(body, atom, frame) != nullptr;
    case ArrayId::Ff14sbCharge:
        if (const StaticNpyArray* a = staticAt(body.run, id)) return staticArrayPresent(a, atom);
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
        return scalarWelfordMeanPresent(welford, atom);
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
    case ArrayId::HbondCount:
        return body.run.h5() && scalarPresent(body.run.h5()->larsenHBondCount(), atom, frame);
    case ArrayId::HmShielding:
        return body.run.h5() && shieldingPresent(body.run.h5()->hmShielding(), atom, frame);
    case ArrayId::WaterEfield:
    case ArrayId::WaterNFirst:
    case ArrayId::WaterNSecond: {
        const auto* ts = body.run.h5() ? body.run.h5()->waterFieldTimeSeries() : nullptr;
        if (!ts || atom >= ts->n_atoms || frame >= ts->n_frames || !ts->sourceAttachedAt(frame))
            return false;
        if (id == ArrayId::WaterEfield) return ts->hasEfield();
        if (id == ArrayId::WaterNFirst) return ts->hasNFirst();
        return ts->hasNSecond();
    }
    case ArrayId::HydrationHalfShellAsymmetry:
    case ArrayId::HydrationDipoleCos: {
        const auto* ts = body.run.h5() ? body.run.h5()->hydrationShellTimeSeries() : nullptr;
        if (!ts || atom >= ts->n_atoms || frame >= ts->n_frames || !ts->sourceAttachedAt(frame))
            return false;
        return id == ArrayId::HydrationHalfShellAsymmetry
                   ? ts->hasHalfShellAsymmetry()
                   : ts->hasMeanWaterDipoleCos();
    }
    case ArrayId::Sasa:
        return body.run.h5() && scalarPresent(body.run.h5()->sasa(), atom, frame);
    case ArrayId::SasaNormal: {
        const auto* ts = body.run.h5() ? body.run.h5()->hydrationGeometryTimeSeries() : nullptr;
        return ts && atom < ts->n_atoms && frame < ts->n_frames && ts->sourceAttachedAt(frame)
               && ts->hasSurfaceNormal();
    }
    case ArrayId::EeqChargeMean:
        return scalarWelfordMeanPresent(body.run.h5() ? body.run.h5()->eeqWelford() : nullptr, atom);
    case ArrayId::HbondScalars:
    case ArrayId::HbondNearestDirection:
    case ArrayId::HbondFlags:
    case ArrayId::EeqCoordinationNumber:
        return false;
    case ArrayId::MopacCoulombEfield:
    case ArrayId::Aimnet2Efg:
    case ArrayId::DsspBackbone:
    case ArrayId::DsspChi:
    case ArrayId::DsspSs8:
    case ArrayId::DsspHBondEnergy:
    case ArrayId::PuckerQ:
    case ArrayId::PuckerTheta:
    case ArrayId::OmegaActual:
    case ArrayId::AromaticChi2:
    case ArrayId::Pyramidalization:
    case ArrayId::LarsenHBondShielding:
        return false;
    default:
        return true;
    }
}

double Catalog::value(const Body& body, ArrayId id, std::size_t atom, std::size_t frame,
                      int slot, int comp) const {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (const StaticNpyArray* a = staticAt(body.run, id)) {
        if (id == ArrayId::Positions) {
            const Vec3 v = valueVec3(body, id, atom, frame);
            return comp >= 0 && comp < 3 ? v[comp] : 0.0;
        }
        if (id == ArrayId::Aimnet2Embedding) return 0.0;
        if (comp >= 0) return staticValue(a, atom, frame, comp);
        return staticValue(a, atom, frame, 0);
    }
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
    case ArrayId::Aimnet2ChargeRespVector:
    case ArrayId::HbondNearestDirection:
    case ArrayId::WaterEfield:
    case ArrayId::SasaNormal: {
        const Vec3 v = valueVec3(body, id, atom, frame);
        return comp >= 0 && comp < 3 ? v[comp] : 0.0;
    }
    case ArrayId::Aimnet2Charge:
        return h5 && h5->aimnet2Charge() ? h5->aimnet2Charge()->at(atom, frame) : 0.0;
    case ArrayId::Aimnet2ChargeRespScalar:
        return h5 && h5->aimnet2ChargeResponseGradient()
                   ? h5->aimnet2ChargeResponseGradient()->scalarAt(atom, frame)
                   : 0.0;
    case ArrayId::HbondScalars:
    case ArrayId::HbondFlags:
        return 0.0;
    case ArrayId::Ff14sbCharge:
        if (const StaticNpyArray* a = staticAt(body.run, id)) return staticValue(a, atom, frame);
        return body.run.protein && atom < body.run.protein->atomCount()
                   ? body.run.protein->atom(atom).partialCharge
                   : 0.0;
    case ArrayId::MopacChargeWelfordMean:
        return present(body, id, atom, frame) && h5 && h5->mopacChargeWelford()
                   ? h5->mopacChargeWelford()->value[atom].mean
                   : 0.0;
    case ArrayId::EeqChargeMean:
        return present(body, id, atom, frame) && h5 && h5->eeqWelford()
                   ? h5->eeqWelford()->value[atom].mean
                   : 0.0;
    case ArrayId::MopacVsFf14sbReconciliation:
        return h5 && h5->mopacVsFf14sbReconciliation()
                   ? h5->mopacVsFf14sbReconciliation()->at(atom, frame)
                   : 0.0;
    case ArrayId::MopacCoulombShielding:
    case ArrayId::Aimnet2Efg:
    case ArrayId::MopacMcShielding:
    case ArrayId::HmShielding:
    case ArrayId::LarsenHBondShielding:
        return comp >= 0 && comp < 5 ? valueT2(body, id, atom, frame)[static_cast<std::size_t>(comp)]
                                     : 0.0;
    case ArrayId::HbondCount:
        return h5 && h5->larsenHBondCount() ? h5->larsenHBondCount()->at(atom, frame) : 0.0;
    case ArrayId::WaterNFirst:
        return h5 && h5->waterFieldTimeSeries()
                   ? static_cast<double>(h5->waterFieldTimeSeries()->nFirstAt(atom, frame))
                   : 0.0;
    case ArrayId::WaterNSecond:
        return h5 && h5->waterFieldTimeSeries()
                   ? static_cast<double>(h5->waterFieldTimeSeries()->nSecondAt(atom, frame))
                   : 0.0;
    case ArrayId::HydrationHalfShellAsymmetry:
        return h5 && h5->hydrationShellTimeSeries()
                   ? h5->hydrationShellTimeSeries()->halfShellAsymmetryAt(atom, frame)
                   : 0.0;
    case ArrayId::HydrationDipoleCos:
        return h5 && h5->hydrationShellTimeSeries()
                   ? h5->hydrationShellTimeSeries()->meanWaterDipoleCosAt(atom, frame)
                   : 0.0;
    case ArrayId::Sasa:
        return h5 && h5->sasa() ? h5->sasa()->at(atom, frame) : 0.0;
    case ArrayId::DftTotalRaw:
        return dftAt(body, atom, frame) ? matComponent(dftAt(body, atom, frame)->total_raw, comp) : 0.0;
    case ArrayId::DftDiaRaw:
        return dftAt(body, atom, frame) ? matComponent(dftAt(body, atom, frame)->dia_raw, comp) : 0.0;
    case ArrayId::DftParaRaw:
        return dftAt(body, atom, frame) ? matComponent(dftAt(body, atom, frame)->para_raw, comp) : 0.0;
    case ArrayId::Aimnet2Embedding:
    case ArrayId::MopacCharge:
    case ArrayId::MopacCoulombEfield:
    case ArrayId::BSPerTypeT0:
    case ArrayId::BSPerTypeT1:
    case ArrayId::BSPerTypeT2:
    case ArrayId::HMPerTypeT0:
    case ArrayId::HMPerTypeT1:
    case ArrayId::HMPerTypeT2:
    case ArrayId::EeqCoordinationNumber:
    case ArrayId::DsspBackbone:
    case ArrayId::DsspChi:
    case ArrayId::DsspSs8:
    case ArrayId::DsspHBondEnergy:
    case ArrayId::PuckerQ:
    case ArrayId::PuckerTheta:
    case ArrayId::OmegaActual:
    case ArrayId::AromaticChi2:
    case ArrayId::Pyramidalization:
    case ArrayId::McPeptideCoFixed:
    case ArrayId::McPeptideCoBo:
    case ArrayId::McPeptideCoRhombic:
    case ArrayId::McPeptideCnFixed:
    case ArrayId::McPeptideCnBo:
    case ArrayId::McBackboneOtherFixed:
    case ArrayId::McBackboneOtherBo:
    case ArrayId::McSidechainCoFixed:
    case ArrayId::McSidechainCoBo:
    case ArrayId::McSidechainOtherFixed:
    case ArrayId::McSidechainOtherBo:
    case ArrayId::McDisulfideFixed:
    case ArrayId::McDisulfideBo:
    case ArrayId::McAromaticZeroedFixed:
    case ArrayId::McAromaticZeroedBo:
    case ArrayId::McNearestCoT2:
    case ArrayId::McNearestCnT2:
        return 0.0;
    }
    return 0.0;
}

Vec3 Catalog::valueVec3(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) const {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (const StaticNpyArray* a = staticAt(body.run, id)) return staticVec3(a, atom, frame);
    switch (id) {
    case ArrayId::Positions:
        return body.run.conformation ? body.run.conformation->atomPosition(frame, atom) : Vec3::Zero();
    case ArrayId::ApbsEfield:
        return h5 && h5->apbsEfield() ? h5->apbsEfield()->at(atom, frame) : Vec3::Zero();
    case ArrayId::Aimnet2ChargeRespVector:
        return h5 && h5->aimnet2ChargeResponseGradient()
                   ? h5->aimnet2ChargeResponseGradient()->vecAt(atom, frame)
                   : Vec3::Zero();
    case ArrayId::WaterEfield:
        return h5 && h5->waterFieldTimeSeries()
                   ? h5->waterFieldTimeSeries()->efieldAt(atom, frame)
                   : Vec3::Zero();
    case ArrayId::SasaNormal:
        return h5 && h5->hydrationGeometryTimeSeries()
                   ? h5->hydrationGeometryTimeSeries()->surfaceNormalAt(atom, frame)
                   : Vec3::Zero();
    default:
        return Vec3::Zero();
    }
}

std::array<double, 5> Catalog::valueT2(const Body& body, ArrayId id, std::size_t atom,
                                       std::size_t frame) const {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (const StaticNpyArray* a = staticAt(body.run, id)) return staticT2(a, atom, frame);
    if (id == ArrayId::ApbsEfg && h5 && h5->apbsEfg()) return h5->apbsEfg()->at(atom, frame);
    // MOPAC-Coulomb-EFG-DERIVED shielding is a T2-only TR (read directly).
    if (id == ArrayId::MopacCoulombShielding && h5 && h5->mopacCoulombShielding())
        return h5->mopacCoulombShielding()->at(atom, frame);
    // MOPAC-McConnell shielding is a full shielding tensor; project its T2 leg.
    if (id == ArrayId::MopacMcShielding && h5 && h5->mopacMcShielding())
        return h5->mopacMcShielding()->at(atom, frame).T2;
    if (id == ArrayId::HmShielding && h5 && h5->hmShielding())
        return h5->hmShielding()->at(atom, frame).T2;
    return {};
}

SphericalTensor Catalog::valueTensor(const Body& body, ArrayId id, std::size_t atom,
                                     std::size_t frame) const {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (const StaticNpyArray* a = staticAt(body.run, id)) return staticTensor(a, atom, frame);
    if (id == ArrayId::KernelBs && h5 && h5->bsShielding()) return h5->bsShielding()->at(atom, frame);
    if (id == ArrayId::KernelMc && h5 && h5->mcShielding()) return h5->mcShielding()->at(atom, frame);
    if (id == ArrayId::MopacMcShielding && h5 && h5->mopacMcShielding())
        return h5->mopacMcShielding()->at(atom, frame);
    if (id == ArrayId::HmShielding && h5 && h5->hmShielding()) return h5->hmShielding()->at(atom, frame);
    return {};
}

const float* Catalog::valueEmbedding(const Body& body, ArrayId id, std::size_t atom,
                                     std::size_t frame, std::size_t& n_dims_out) const {
    n_dims_out = 0;
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (const StaticNpyArray* a = staticAt(body.run, id)) {
        if (id != ArrayId::Aimnet2Embedding || atom >= a->rows || a->cols == 0
            || a->floatValues.empty())
            return nullptr;
        n_dims_out = a->cols;
        return a->floatRowData(atom);
    }
    if (id != ArrayId::Aimnet2Embedding || !h5 || !h5->aimnet2Embedding()) return nullptr;
    n_dims_out = h5->aimnet2Embedding()->n_dims;
    return h5->aimnet2Embedding()->dataAt(atom, frame);
}

}  // namespace h5reader::rediscover
