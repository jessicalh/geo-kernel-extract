#include "Catalog.h"

#include "AnalysisBody.h"
#include "RunData.h"
#include "ScopedProducerCatalog.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/Conformation.h"
#include "../model/DftShielding.h"
#include "../model/QtAimnet2Group.h"
#include "../model/QtAtom.h"
#include "../model/QtPerResidueBuffers.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtResultBlocks.h"
#include "../model/QtSpecialBuffers.h"
#include "../model/QtTimeSeriesBuffers.h"
#include "../model/QtTopology.h"

#include <cmath>
#include <limits>
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

QString FieldProviderName(FieldProvider provider) {
    switch (provider) {
    case FieldProvider::StaticProducerArray: return QStringLiteral("static_producer_array");
    case FieldProvider::DenseH5TimeSeries: return QStringLiteral("dense_h5_time_series");
    case FieldProvider::SparseDftByOriginal: return QStringLiteral("sparse_dft_by_original");
    case FieldProvider::TypedTopology: return QStringLiteral("typed_topology");
    case FieldProvider::DatasetAbsent: return QStringLiteral("dataset_absent");
    }
    return QStringLiteral("unknown");
}

QString ArrayResidenceName(ArrayResidence residence) {
    switch (residence) {
    case ArrayResidence::DenseH5: return QStringLiteral("dense_h5");
    case ArrayResidence::StaticTopol: return QStringLiteral("static_topology");
    case ArrayResidence::StaticNpy: return QStringLiteral("static_npy");
    case ArrayResidence::SparseDftByOriginal: return QStringLiteral("sparse_dft_by_original");
    case ArrayResidence::Absent: return QStringLiteral("absent");
    }
    return QStringLiteral("unknown");
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

int fieldOrd(io::FieldKind kind) { return static_cast<int>(kind); }

std::size_t nativeRowForStatic(const StaticNpyArray& a, std::size_t nativeRow, std::size_t frame) {
    return a.frameVarying ? frame * a.atomsPerFrame + nativeRow : nativeRow;
}

bool isTypedTopologyField(io::FieldKind kind) {
    switch (kind) {
    case io::FieldKind::Element:
    case io::FieldKind::ResidueIndex:
    case io::FieldKind::ResidueType:
    case io::FieldKind::FfPartialCharge:
    case io::FieldKind::AtomsCategoryInfo:
    case io::FieldKind::Residues:
    case io::FieldKind::Bonds:
    case io::FieldKind::Rings:
    case io::FieldKind::RingMembership:
        return true;
    default:
        return false;
    }
}

bool isDftField(io::FieldKind kind) {
    return kind == io::FieldKind::OrcaTotal
           || kind == io::FieldKind::OrcaDiamagnetic
           || kind == io::FieldKind::OrcaParamagnetic;
}

std::size_t topologyRows(const RunData& run, io::NativeAxis axis, io::FieldKind kind);

std::optional<std::size_t> fixedNativeAxisRows(const RunData& run,
                                               io::NativeAxis axis,
                                               io::FieldKind kind) {
    switch (axis) {
    case io::NativeAxis::Atom:
    case io::NativeAxis::Residue:
    case io::NativeAxis::Bond:
    case io::NativeAxis::Ring:
    case io::NativeAxis::RingMembership:
    case io::NativeAxis::AromaticRing:
    case io::NativeAxis::SaturatedRing:
    case io::NativeAxis::Protein:
        return topologyRows(run, axis, kind);
    default:
        return std::nullopt;
    }
}

bool topologyFieldAvailable(const RunData& run, const io::FieldSpec& field, std::size_t rows) {
    if (!run.protein) return false;
    if (field.kind == io::FieldKind::FfPartialCharge) {
        if (rows == 0) return false;
        for (std::size_t atom = 0; atom < rows; ++atom) {
            if (!run.protein->atom(atom).hasPartialCharge) return false;
        }
        return true;
    }
    if (IsStructuredProducerField(field)) return true;
    return rows > 0;
}

bool denseH5ProviderAvailable(const RunData& run, io::FieldKind kind) {
    const io::QtTrajectoryH5* h5 = run.h5();
    if (!h5) return false;
    switch (kind) {
    case io::FieldKind::Pos: return h5->positions() != nullptr;
    case io::FieldKind::BSShielding: return h5->bsShielding() != nullptr;
    case io::FieldKind::HMShielding: return h5->hmShielding() != nullptr;
    case io::FieldKind::APBSE: return h5->apbsEfield() != nullptr;
    case io::FieldKind::APBSEFG: return h5->apbsEfg() != nullptr;
    case io::FieldKind::AIMNet2Charges: return h5->aimnet2Charge() != nullptr;
    case io::FieldKind::AIMNet2Aim: return h5->aimnet2Embedding() != nullptr;
    case io::FieldKind::AIMNet2ChargeResponseGradient:
    case io::FieldKind::AIMNet2ChargeResponseGradientScalar:
        return h5->aimnet2ChargeResponseGradient() != nullptr;
    case io::FieldKind::AtomSASA: return h5->sasa() != nullptr;
    case io::FieldKind::SASANormal:
        return h5->hydrationGeometryTimeSeries()
               && h5->hydrationGeometryTimeSeries()->hasSurfaceNormal();
    case io::FieldKind::WaterEfield:
    case io::FieldKind::WaterEfieldFirst:
    case io::FieldKind::WaterEFG:
    case io::FieldKind::WaterEFGFirst:
    case io::FieldKind::WaterShellCounts:
        return h5->waterFieldTimeSeries() != nullptr;
    case io::FieldKind::HydrationShell:
        return h5->hydrationShellTimeSeries() != nullptr;
    case io::FieldKind::WaterPolarization:
        return h5->hydrationGeometryTimeSeries() != nullptr;
    case io::FieldKind::BondedEnergy:
        return h5->bondedEnergy() != nullptr;
    case io::FieldKind::GromacsEnergy:
        return h5->gromacsEnergy() != nullptr;
    case io::FieldKind::DSSPSs8:
    case io::FieldKind::DSSPHBondEnergy:
        return h5->dssp8() != nullptr;
    case io::FieldKind::DSSPChi:
        return h5->dihedrals() != nullptr;
    case io::FieldKind::OmegaActual:
    case io::FieldKind::OmegaDeviation:
    case io::FieldKind::OmegaIsXpro:
        return h5->dihedrals() != nullptr;
    case io::FieldKind::AromaticChi2:
    case io::FieldKind::PuckerQ:
    case io::FieldKind::PuckerTheta:
        return h5->ringPucker() != nullptr;
    default:
        return false;
    }
}

std::size_t staticNativeRows(const StaticNpyArray* a) {
    if (!a) return 0;
    return a->frameVarying ? a->atomsPerFrame : a->rows;
}

std::size_t topologyRows(const RunData& run, io::NativeAxis axis, io::FieldKind kind) {
    if (!run.protein) return 0;
    const model::QtTopology& topo = run.protein->topology();
    switch (kind) {
    case io::FieldKind::AtomsCategoryInfo:
        return run.protein->atomCount();
    default:
        break;
    }
    switch (axis) {
    case io::NativeAxis::Atom: return run.protein->atomCount();
    case io::NativeAxis::Residue: return run.protein->residueCount();
    case io::NativeAxis::Bond: return topo.bondCount();
    case io::NativeAxis::Ring: return topo.ringCount();
    case io::NativeAxis::RingMembership: return topo.ringMembershipCount();
    case io::NativeAxis::AromaticRing: return topo.aromaticRingCount();
    case io::NativeAxis::SaturatedRing: return topo.saturatedRingCount();
    case io::NativeAxis::Protein: return 1;
    default: return 0;
    }
}

std::size_t denseNativeRows(const RunData& run, io::FieldKind kind) {
    const io::QtTrajectoryH5* h5 = run.h5();
    if (!h5) return 0;
    switch (kind) {
    case io::FieldKind::Pos: return h5->positions() ? h5->positions()->n_atoms : 0;
    case io::FieldKind::BSShielding: return h5->bsShielding() ? h5->bsShielding()->n_atoms : 0;
    case io::FieldKind::HMShielding: return h5->hmShielding() ? h5->hmShielding()->n_atoms : 0;
    case io::FieldKind::APBSE: return h5->apbsEfield() ? h5->apbsEfield()->n_atoms : 0;
    case io::FieldKind::APBSEFG: return h5->apbsEfg() ? h5->apbsEfg()->n_atoms : 0;
    case io::FieldKind::AIMNet2Charges: return h5->aimnet2Charge() ? h5->aimnet2Charge()->n_atoms : 0;
    case io::FieldKind::AIMNet2Aim: return h5->aimnet2Embedding() ? h5->aimnet2Embedding()->n_atoms : 0;
    case io::FieldKind::AIMNet2ChargeResponseGradient:
    case io::FieldKind::AIMNet2ChargeResponseGradientScalar:
        return h5->aimnet2ChargeResponseGradient() ? h5->aimnet2ChargeResponseGradient()->n_atoms : 0;
    case io::FieldKind::AtomSASA: return h5->sasa() ? h5->sasa()->n_atoms : 0;
    case io::FieldKind::SASANormal:
        return h5->hydrationGeometryTimeSeries() ? h5->hydrationGeometryTimeSeries()->n_atoms : 0;
    case io::FieldKind::WaterEfield:
    case io::FieldKind::WaterEfieldFirst:
    case io::FieldKind::WaterEFG:
    case io::FieldKind::WaterEFGFirst:
    case io::FieldKind::WaterShellCounts:
        return h5->waterFieldTimeSeries() ? h5->waterFieldTimeSeries()->n_atoms : 0;
    case io::FieldKind::HydrationShell:
        return h5->hydrationShellTimeSeries() ? h5->hydrationShellTimeSeries()->n_atoms : 0;
    case io::FieldKind::WaterPolarization:
        return h5->hydrationGeometryTimeSeries() ? h5->hydrationGeometryTimeSeries()->n_atoms : 0;
    case io::FieldKind::BondedEnergy:
        return h5->bondedEnergy() ? h5->bondedEnergy()->n_atoms : 0;
    case io::FieldKind::GromacsEnergy:
        return h5->gromacsEnergy() ? 1 : 0;
    case io::FieldKind::DSSPSs8:
    case io::FieldKind::DSSPHBondEnergy:
        return h5->dssp8() ? h5->dssp8()->n_atoms : 0;
    case io::FieldKind::DSSPChi:
        return run.protein ? run.protein->atomCount() : 0;
    case io::FieldKind::OmegaActual:
    case io::FieldKind::OmegaDeviation:
    case io::FieldKind::OmegaIsXpro:
        return h5->dihedrals() ? h5->dihedrals()->n_residues : 0;
    case io::FieldKind::AromaticChi2:
        return h5->ringPucker() ? h5->ringPucker()->n_aromatic_rings : 0;
    case io::FieldKind::PuckerQ:
    case io::FieldKind::PuckerTheta:
        return h5->ringPucker() ? h5->ringPucker()->n_saturated_rings : 0;
    default:
        return 0;
    }
}

std::size_t denseFrameCount(const RunData& run, io::FieldKind kind) {
    const io::QtTrajectoryH5* h5 = run.h5();
    if (!h5) return 0;
    switch (kind) {
    case io::FieldKind::Pos: return h5->positions() ? h5->positions()->n_frames : 0;
    case io::FieldKind::BSShielding: return h5->bsShielding() ? h5->bsShielding()->n_frames : 0;
    case io::FieldKind::HMShielding: return h5->hmShielding() ? h5->hmShielding()->n_frames : 0;
    case io::FieldKind::APBSE: return h5->apbsEfield() ? h5->apbsEfield()->n_frames : 0;
    case io::FieldKind::APBSEFG: return h5->apbsEfg() ? h5->apbsEfg()->n_frames : 0;
    case io::FieldKind::AIMNet2Charges: return h5->aimnet2Charge() ? h5->aimnet2Charge()->n_frames : 0;
    case io::FieldKind::AIMNet2Aim: return h5->aimnet2Embedding() ? h5->aimnet2Embedding()->n_frames : 0;
    case io::FieldKind::AIMNet2ChargeResponseGradient:
    case io::FieldKind::AIMNet2ChargeResponseGradientScalar:
        return h5->aimnet2ChargeResponseGradient() ? h5->aimnet2ChargeResponseGradient()->n_frames : 0;
    case io::FieldKind::AtomSASA: return h5->sasa() ? h5->sasa()->n_frames : 0;
    case io::FieldKind::SASANormal:
        return h5->hydrationGeometryTimeSeries() ? h5->hydrationGeometryTimeSeries()->n_frames : 0;
    case io::FieldKind::WaterEfield:
    case io::FieldKind::WaterEfieldFirst:
    case io::FieldKind::WaterEFG:
    case io::FieldKind::WaterEFGFirst:
    case io::FieldKind::WaterShellCounts:
        return h5->waterFieldTimeSeries() ? h5->waterFieldTimeSeries()->n_frames : 0;
    case io::FieldKind::HydrationShell:
        return h5->hydrationShellTimeSeries() ? h5->hydrationShellTimeSeries()->n_frames : 0;
    case io::FieldKind::WaterPolarization:
        return h5->hydrationGeometryTimeSeries() ? h5->hydrationGeometryTimeSeries()->n_frames : 0;
    case io::FieldKind::BondedEnergy:
        return h5->bondedEnergy() ? h5->bondedEnergy()->n_frames : 0;
    case io::FieldKind::GromacsEnergy:
        return h5->gromacsEnergy() ? h5->gromacsEnergy()->n_frames : 0;
    case io::FieldKind::DSSPSs8:
    case io::FieldKind::DSSPHBondEnergy:
        return h5->dssp8() ? h5->dssp8()->n_frames : 0;
    case io::FieldKind::DSSPChi:
        return h5->dihedrals() ? h5->dihedrals()->n_frames : 0;
    case io::FieldKind::OmegaIsXpro:
        return 1;
    case io::FieldKind::OmegaActual:
    case io::FieldKind::OmegaDeviation:
        return h5->dihedrals() ? h5->dihedrals()->n_frames : 0;
    case io::FieldKind::AromaticChi2:
    case io::FieldKind::PuckerQ:
    case io::FieldKind::PuckerTheta:
        return h5->ringPucker() ? h5->ringPucker()->n_frames : 0;
    default:
        return 0;
    }
}

std::size_t denseComponentCount(const RunData& run, const io::FieldSpec& spec) {
    if (spec.kind == io::FieldKind::AIMNet2Aim) {
        const auto* emb = run.h5() ? run.h5()->aimnet2Embedding() : nullptr;
        return emb ? emb->n_dims : 0;
    }
    if (spec.cols > 0) return static_cast<std::size_t>(spec.cols);
    return 1;
}

FieldAccessSpec makeFieldAccessSpec(const RunData& run, const io::FieldSpec& field) {
    FieldAccessSpec out;
    out.kind = field.kind;
    out.stem = FieldStem(field);
    out.axis = field.axis;
    out.structured = IsStructuredProducerField(field);

    if (!IsScopedProducerField(field)) {
        out.absence_reason = QStringLiteral("out-of-scope");
        return out;
    }

    if (const StaticNpyArray* a = run.producerArray(field.kind)) {
        out.provider = FieldProvider::StaticProducerArray;
        out.residence = ArrayResidence::StaticNpy;
        out.native_rows = staticNativeRows(a);
        out.frames = a->frameVarying ? a->frameCount : 1;
        out.components = a->cols;
        const std::optional<std::size_t> expectedRows =
            fixedNativeAxisRows(run, field.axis, field.kind);
        const bool rowShapeOk = !expectedRows || out.native_rows > 0 || *expectedRows == 0;
        out.available = rowShapeOk && out.frames > 0 && out.components > 0;
        if (!out.available)
            out.absence_reason = rowShapeOk ? QStringLiteral("malformed-shape")
                                            : QStringLiteral("native-row-count-mismatch");
        return out;
    }

    const QString producerIssue = run.producerArrayIssue(field.kind);
    if (!producerIssue.isEmpty()) {
        out.provider = FieldProvider::DatasetAbsent;
        out.residence = ArrayResidence::Absent;
        out.components = field.cols > 0 ? static_cast<std::size_t>(field.cols) : 0;
        out.absence_reason = QStringLiteral("malformed-shape: %1").arg(producerIssue);
        return out;
    }

    if (isTypedTopologyField(field.kind)) {
        out.provider = FieldProvider::TypedTopology;
        out.residence = ArrayResidence::StaticTopol;
        out.native_rows = topologyRows(run, field.axis, field.kind);
        out.frames = 1;
        out.components = out.structured ? 0 : 1;
        out.available = topologyFieldAvailable(run, field, out.native_rows)
                        && (out.structured || out.components > 0);
        if (!out.available) out.absence_reason = QStringLiteral("topology-mismatch");
        return out;
    }

    if (isDftField(field.kind)) {
        out.provider = FieldProvider::SparseDftByOriginal;
        out.residence = ArrayResidence::SparseDftByOriginal;
        out.native_rows = run.protein ? run.protein->atomCount() : 0;
        out.frames = run.frameMap.dftRows().size();
        out.components = 9;
        out.available = out.native_rows > 0 && out.frames > 0;
        if (!out.available) out.absence_reason = QStringLiteral("not-produced-in-dataset");
        return out;
    }

    if (denseH5ProviderAvailable(run, field.kind)) {
        out.provider = FieldProvider::DenseH5TimeSeries;
        out.residence = ArrayResidence::DenseH5;
        out.native_rows = denseNativeRows(run, field.kind);
        out.frames = denseFrameCount(run, field.kind);
        out.components = denseComponentCount(run, field);
        const std::optional<std::size_t> expectedRows =
            fixedNativeAxisRows(run, field.axis, field.kind);
        const bool rowShapeOk = !expectedRows || *expectedRows == out.native_rows;
        out.available = rowShapeOk && out.frames > 0 && out.components > 0;
        if (!out.available)
            out.absence_reason = rowShapeOk ? QStringLiteral("malformed-shape")
                                            : QStringLiteral("native-row-count-mismatch");
        return out;
    }

    out.provider = FieldProvider::DatasetAbsent;
    out.residence = ArrayResidence::Absent;
    out.components = field.cols > 0 ? static_cast<std::size_t>(field.cols) : 0;
    out.absence_reason = run.poseKind() == PoseKind::Trajectory
                             ? QStringLiteral("not-produced-in-dataset")
                             : QStringLiteral("producer-array-absent");
    return out;
}

const std::vector<double>* gromacsChannel(const model::QtSystemEnergyTimeSeries& ts,
                                          std::size_t component,
                                          std::size_t* offset = nullptr) {
    if (offset) *offset = 0;
    switch (component) {
    case 0: return &ts.kinetic;
    case 1: return &ts.potential;
    case 2: return &ts.total_energy;
    case 3: return &ts.enthalpy;
    case 4: return &ts.temperature;
    case 5: return &ts.T_protein;
    case 6: return &ts.T_non_protein;
    case 7: return &ts.pressure;
    case 8: return &ts.density;
    case 9: return &ts.volume;
    case 10: return &ts.box_x;
    case 11: return &ts.box_y;
    case 12: return &ts.box_z;
    case 13: return &ts.bond;
    case 14: return &ts.angle;
    case 15: return &ts.proper_dih;
    case 16: return &ts.improper_dih;
    case 17: return &ts.urey_bradley;
    case 18: return &ts.cmap_dih;
    case 19: return &ts.coulomb_sr;
    case 20: return &ts.coulomb_recip;
    case 21: return &ts.coulomb_14;
    case 22: return &ts.lj_sr;
    case 23: return &ts.lj_14;
    case 24: return &ts.disper_corr;
    default:
        if (component >= 25 && component < 34) {
            if (offset) *offset = component - 25;
            return &ts.pressure_tensor;
        }
        if (component >= 34 && component < 43) {
            if (offset) *offset = component - 34;
            return &ts.virial;
        }
        return nullptr;
    }
}

bool gromacsComponentPresent(const model::QtSystemEnergyTimeSeries* ts,
                             std::size_t frame,
                             std::size_t component) {
    if (!ts || frame >= ts->n_frames || !ts->sourceAttachedAt(frame)) return false;
    std::size_t offset = 0;
    const std::vector<double>* ch = gromacsChannel(*ts, component, &offset);
    if (!ch) return false;
    if (component < 25) return ch->size() == ts->n_frames;
    return ch->size() == ts->n_frames * 9 && offset < 9;
}

double gromacsComponentValue(const model::QtSystemEnergyTimeSeries& ts,
                             std::size_t frame,
                             std::size_t component) {
    std::size_t offset = 0;
    const std::vector<double>* ch = gromacsChannel(ts, component, &offset);
    if (!ch) return 0.0;
    return component < 25 ? (*ch)[frame] : (*ch)[frame * 9 + offset];
}

const std::vector<double>* hydrationShellChannel(const model::QtHydrationShellTimeSeries& ts,
                                                 std::size_t component) {
    switch (component) {
    case 0: return &ts.half_shell_asymmetry;
    case 1: return &ts.mean_water_dipole_cos;
    case 2: return &ts.nearest_ion_charge;
    case 3: return &ts.nearest_ion_distance;
    default: return nullptr;
    }
}

const std::vector<double>* bondedEnergyChannel(const model::QtBondedEnergyTimeSeries& ts,
                                               std::size_t component) {
    switch (component) {
    case 0: return &ts.bond;
    case 1: return &ts.angle;
    case 2: return &ts.proper_dih;
    case 3: return &ts.improper_dih;
    case 4: return &ts.urey_bradley;
    case 5: return &ts.cmap_dih;
    case 6: return &ts.total;
    default: return nullptr;
    }
}

std::optional<std::size_t> residueForAtom(const RunData& run, std::size_t atom) {
    if (!run.protein || atom >= run.protein->atomCount()) return std::nullopt;
    const int32_t residue = run.protein->atom(atom).residueIndex;
    if (residue < 0) return std::nullopt;
    const std::size_t row = static_cast<std::size_t>(residue);
    if (row >= run.protein->residueCount()) return std::nullopt;
    return row;
}

std::optional<std::size_t> dsspResidueForAtom(const RunData& run,
                                              const model::QtDssp8TimeSeries* ts,
                                              std::size_t atom) {
    if (!ts || atom >= ts->n_atoms) return std::nullopt;
    if (ts->residue_index_per_atom.size() == ts->n_atoms) {
        const int32_t residue = ts->residue_index_per_atom[atom];
        if (residue < 0) return std::nullopt;
        const std::size_t row = static_cast<std::size_t>(residue);
        return row < ts->n_residues ? std::optional<std::size_t>(row) : std::nullopt;
    }
    const std::optional<std::size_t> row = residueForAtom(run, atom);
    return row && *row < ts->n_residues ? row : std::nullopt;
}

bool hydrationGeometryHasWaterPolarization(const model::QtHydrationGeometryTimeSeries& ts) {
    const std::size_t nt = ts.n_atoms * ts.n_frames;
    return ts.dipole_vector.size() == nt * 3
           && ts.surface_normal.size() == nt * 3
           && ts.half_shell_asymmetry.size() == nt
           && ts.dipole_alignment.size() == nt
           && ts.dipole_coherence.size() == nt
           && ts.first_shell_count.size() == nt;
}

double hydrationGeometryWaterPolarizationValue(const model::QtHydrationGeometryTimeSeries& ts,
                                               std::size_t atom,
                                               std::size_t frame,
                                               std::size_t component) {
    const std::size_t scalar = atom * ts.n_frames + frame;
    const std::size_t vec = scalar * 3;
    switch (component) {
    case 0:
    case 1:
    case 2:
        return ts.dipole_vector[vec + component];
    case 3:
    case 4:
    case 5:
        return ts.surface_normal[vec + component - 3];
    case 6:
        return ts.half_shell_asymmetry[scalar];
    case 7:
        return ts.dipole_alignment[scalar];
    case 8:
        return ts.dipole_coherence[scalar];
    case 9:
        return static_cast<double>(ts.first_shell_count[scalar]);
    default:
        return 0.0;
    }
}

double dsspSs8Value(model::DsspCode code, std::size_t component) {
    if (component >= 8 || code == model::DsspCode::Unknown) return 0.0;
    return static_cast<std::size_t>(static_cast<uint8_t>(code)) == component ? 1.0 : 0.0;
}

double dsspHBondEnergyValue(const model::QtDssp8TimeSeries& ts,
                            std::size_t residue,
                            std::size_t frame,
                            std::size_t component) {
    const std::size_t base = (residue * ts.n_frames + frame) * 2;
    if (component < 2) return ts.hbond_acceptor_energy[base + component];
    return ts.hbond_donor_energy[base + component - 2];
}

double dsspChiValue(const model::QtDihedralTimeSeries& ts,
                    std::size_t residue,
                    std::size_t frame,
                    std::size_t component) {
    const std::size_t chi = component / 3;
    const std::size_t slot = component % 3;
    const bool exists = ts.chi_exists.size() == ts.n_residues * 4
                        && ts.chi_exists[residue * 4 + chi] != 0;
    if (slot == 2) return exists ? 1.0 : 0.0;
    if (!exists) return 0.0;
    const double angle = ts.chiAt(residue, frame, static_cast<int>(chi));
    return slot == 0 ? std::cos(angle) : std::sin(angle);
}

bool denseFieldPresent(const RunData& run,
                       io::FieldKind kind,
                       std::size_t nativeRow,
                       std::size_t frame,
                       std::size_t component,
                       QString* reason) {
    const io::QtTrajectoryH5* h5 = run.h5();
    auto absent = [&](const QString& r) {
        if (reason) *reason = r;
        return false;
    };
    if (!h5) return absent(QStringLiteral("unsupported-in-residence"));
    switch (kind) {
    case io::FieldKind::Pos: {
        const auto* ts = h5->positions();
        if (!ts) return absent(QStringLiteral("not-produced-in-dataset"));
        if (nativeRow >= ts->n_atoms) return absent(QStringLiteral("native-row-out-of-range"));
        if (frame >= ts->n_frames) return absent(QStringLiteral("frame-out-of-range"));
        return component < 3 || absent(QStringLiteral("component-out-of-range"));
    }
    case io::FieldKind::BSShielding:
        return shieldingPresent(h5->bsShielding(), nativeRow, frame)
               && (component < 9 || absent(QStringLiteral("component-out-of-range")));
    case io::FieldKind::HMShielding:
        return shieldingPresent(h5->hmShielding(), nativeRow, frame)
               && (component < 9 || absent(QStringLiteral("component-out-of-range")));
    case io::FieldKind::APBSE:
        return h5->apbsEfield() && nativeRow < h5->apbsEfield()->n_atoms
               && frame < h5->apbsEfield()->n_frames
               && h5->apbsEfield()->sourceAttachedAt(frame)
               && (component < 3 || absent(QStringLiteral("component-out-of-range")));
    case io::FieldKind::APBSEFG:
        return h5->apbsEfg() && nativeRow < h5->apbsEfg()->n_atoms
               && frame < h5->apbsEfg()->n_frames
               && h5->apbsEfg()->sourceAttachedAt(frame)
               && (component < 5 || absent(QStringLiteral("component-out-of-range")));
    case io::FieldKind::AIMNet2Charges:
        return h5->aimnet2Charge() && nativeRow < h5->aimnet2Charge()->n_atoms
               && frame < h5->aimnet2Charge()->n_frames
               && h5->aimnet2Charge()->sourceAttachedAt(frame)
               && component == 0;
    case io::FieldKind::AIMNet2Aim:
        return h5->aimnet2Embedding() && nativeRow < h5->aimnet2Embedding()->n_atoms
               && frame < h5->aimnet2Embedding()->n_frames
               && component < h5->aimnet2Embedding()->n_dims
               && h5->aimnet2Embedding()->meta.sourceAttachedAt(frame);
    case io::FieldKind::AIMNet2ChargeResponseGradientScalar:
        return h5->aimnet2ChargeResponseGradient()
               && nativeRow < h5->aimnet2ChargeResponseGradient()->n_atoms
               && frame < h5->aimnet2ChargeResponseGradient()->n_frames
               && h5->aimnet2ChargeResponseGradient()->meta.sourceAttachedAt(frame)
               && component == 0;
    case io::FieldKind::AIMNet2ChargeResponseGradient:
        return h5->aimnet2ChargeResponseGradient()
               && nativeRow < h5->aimnet2ChargeResponseGradient()->n_atoms
               && frame < h5->aimnet2ChargeResponseGradient()->n_frames
               && h5->aimnet2ChargeResponseGradient()->meta.sourceAttachedAt(frame)
               && (component < 3 || absent(QStringLiteral("component-out-of-range")));
    case io::FieldKind::AtomSASA:
        return scalarPresent(h5->sasa(), nativeRow, frame) && component == 0;
    case io::FieldKind::SASANormal: {
        const auto* ts = h5->hydrationGeometryTimeSeries();
        return ts && nativeRow < ts->n_atoms && frame < ts->n_frames
               && ts->sourceAttachedAt(frame) && ts->hasSurfaceNormal()
               && (component < 3 || absent(QStringLiteral("component-out-of-range")));
    }
    case io::FieldKind::WaterEfield:
    case io::FieldKind::WaterEfieldFirst:
    case io::FieldKind::WaterEFG:
    case io::FieldKind::WaterEFGFirst:
    case io::FieldKind::WaterShellCounts: {
        const auto* ts = h5->waterFieldTimeSeries();
        if (!ts || nativeRow >= ts->n_atoms || frame >= ts->n_frames || !ts->sourceAttachedAt(frame))
            return absent(QStringLiteral("not-present"));
        if (kind == io::FieldKind::WaterEfield)
            return ts->efield.size() == ts->n_atoms * ts->n_frames * 3
                   && (component < 3 || absent(QStringLiteral("component-out-of-range")));
        if (kind == io::FieldKind::WaterEfieldFirst)
            return ts->efield_first.size() == ts->n_atoms * ts->n_frames * 3
                   && (component < 3 || absent(QStringLiteral("component-out-of-range")));
        if (kind == io::FieldKind::WaterEFG)
            return ts->efg.size() == ts->n_atoms * ts->n_frames * 5
                   && (component < 5 || absent(QStringLiteral("component-out-of-range")));
        if (kind == io::FieldKind::WaterEFGFirst)
            return ts->efg_first.size() == ts->n_atoms * ts->n_frames * 5
                   && (component < 5 || absent(QStringLiteral("component-out-of-range")));
        return ts->n_first.size() == ts->n_atoms * ts->n_frames
               && ts->n_second.size() == ts->n_atoms * ts->n_frames
               && (component < 2 || absent(QStringLiteral("component-out-of-range")));
    }
    case io::FieldKind::HydrationShell: {
        const auto* ts = h5->hydrationShellTimeSeries();
        const std::vector<double>* ch = ts ? hydrationShellChannel(*ts, component) : nullptr;
        return ts && ch && nativeRow < ts->n_atoms && frame < ts->n_frames
               && ts->sourceAttachedAt(frame) && ch->size() == ts->n_atoms * ts->n_frames;
    }
    case io::FieldKind::WaterPolarization: {
        const auto* ts = h5->hydrationGeometryTimeSeries();
        return ts && nativeRow < ts->n_atoms && frame < ts->n_frames
               && ts->sourceAttachedAt(frame)
               && hydrationGeometryHasWaterPolarization(*ts)
               && (component < 10 || absent(QStringLiteral("component-out-of-range")));
    }
    case io::FieldKind::BondedEnergy: {
        const auto* ts = h5->bondedEnergy();
        const std::vector<double>* ch = ts ? bondedEnergyChannel(*ts, component) : nullptr;
        return ts && ch && nativeRow < ts->n_atoms && frame < ts->n_frames
               && ts->sourceAttachedAt(frame) && ch->size() == ts->n_atoms * ts->n_frames;
    }
    case io::FieldKind::GromacsEnergy:
        if (nativeRow != 0) return absent(QStringLiteral("native-row-out-of-range"));
        return gromacsComponentPresent(h5->gromacsEnergy(), frame, component);
    case io::FieldKind::DSSPSs8: {
        const auto* ts = h5->dssp8();
        const std::optional<std::size_t> residue = dsspResidueForAtom(run, ts, nativeRow);
        if (!ts || !residue || frame >= ts->n_frames || !ts->sourceAttachedAt(frame))
            return absent(QStringLiteral("not-present"));
        if (ts->ss8_code.size() != ts->n_residues * ts->n_frames)
            return absent(QStringLiteral("malformed-shape"));
        return component < 8 || absent(QStringLiteral("component-out-of-range"));
    }
    case io::FieldKind::DSSPHBondEnergy: {
        const auto* ts = h5->dssp8();
        const std::optional<std::size_t> residue = dsspResidueForAtom(run, ts, nativeRow);
        if (!ts || !residue || frame >= ts->n_frames || !ts->sourceAttachedAt(frame))
            return absent(QStringLiteral("not-present"));
        const std::size_t need = ts->n_residues * ts->n_frames * 2;
        if (ts->hbond_acceptor_energy.size() != need || ts->hbond_donor_energy.size() != need)
            return absent(QStringLiteral("malformed-shape"));
        return component < 4 || absent(QStringLiteral("component-out-of-range"));
    }
    case io::FieldKind::DSSPChi: {
        const auto* ts = h5->dihedrals();
        const std::optional<std::size_t> residue = residueForAtom(run, nativeRow);
        if (!ts || !residue || *residue >= ts->n_residues
            || frame >= ts->n_frames || !ts->sourceAttachedAt(frame))
            return absent(QStringLiteral("not-present"));
        if (ts->chi.size() != ts->n_residues * ts->n_frames * 4
            || ts->chi_exists.size() != ts->n_residues * 4)
            return absent(QStringLiteral("malformed-shape"));
        return component < 12 || absent(QStringLiteral("component-out-of-range"));
    }
    case io::FieldKind::OmegaActual:
    case io::FieldKind::OmegaDeviation: {
        const auto* ts = h5->dihedrals();
        if (!ts || nativeRow >= ts->n_residues || frame >= ts->n_frames || !ts->sourceAttachedAt(frame))
            return absent(QStringLiteral("not-present"));
        if (kind == io::FieldKind::OmegaDeviation
            && ts->omega_deviation.size() != ts->n_residues * ts->n_frames)
            return absent(QStringLiteral("malformed-shape"));
        return component == 0;
    }
    case io::FieldKind::OmegaIsXpro: {
        const auto* ts = h5->dihedrals();
        return ts && nativeRow < ts->n_residues && frame == 0
               && nativeRow < ts->omega_is_xpro.size() && component == 0;
    }
    case io::FieldKind::AromaticChi2: {
        const auto* ts = h5->ringPucker();
        return ts && nativeRow < ts->n_aromatic_rings && frame < ts->n_frames
               && (ts->source_attached.empty()
                   || (frame < ts->source_attached.size() && ts->source_attached[frame] != 0))
               && component == 0;
    }
    case io::FieldKind::PuckerQ:
    case io::FieldKind::PuckerTheta: {
        const auto* ts = h5->ringPucker();
        if (!ts || nativeRow >= ts->n_saturated_rings || frame >= ts->n_frames)
            return absent(QStringLiteral("not-present"));
        if (!ts->source_attached.empty()
            && (frame >= ts->source_attached.size() || ts->source_attached[frame] == 0))
            return absent(QStringLiteral("frame-gap"));
        if (kind == io::FieldKind::PuckerTheta
            && ts->pucker_theta.size() != ts->n_saturated_rings * ts->n_frames)
            return absent(QStringLiteral("malformed-shape"));
        return component == 0;
    }
    default:
        return absent(QStringLiteral("not-produced-in-dataset"));
    }
}

double denseFieldValue(const RunData& run,
                       io::FieldKind kind,
                       std::size_t nativeRow,
                       std::size_t frame,
                       std::size_t component) {
    const io::QtTrajectoryH5* h5 = run.h5();
    if (!h5) return 0.0;
    switch (kind) {
    case io::FieldKind::Pos: {
        const Vec3 v = h5->positions()->at(nativeRow, frame);
        return v[static_cast<int>(component)];
    }
    case io::FieldKind::BSShielding:
        return tensorComponent(h5->bsShielding()->at(nativeRow, frame), static_cast<int>(component));
    case io::FieldKind::HMShielding:
        return tensorComponent(h5->hmShielding()->at(nativeRow, frame), static_cast<int>(component));
    case io::FieldKind::APBSE: {
        const Vec3 v = h5->apbsEfield()->at(nativeRow, frame);
        return v[static_cast<int>(component)];
    }
    case io::FieldKind::APBSEFG:
        return h5->apbsEfg()->at(nativeRow, frame)[component];
    case io::FieldKind::AIMNet2Charges:
        return h5->aimnet2Charge()->at(nativeRow, frame);
    case io::FieldKind::AIMNet2Aim:
        return static_cast<double>(h5->aimnet2Embedding()->dataAt(nativeRow, frame)[component]);
    case io::FieldKind::AIMNet2ChargeResponseGradientScalar:
        return h5->aimnet2ChargeResponseGradient()->scalarAt(nativeRow, frame);
    case io::FieldKind::AIMNet2ChargeResponseGradient: {
        const Vec3 v = h5->aimnet2ChargeResponseGradient()->vecAt(nativeRow, frame);
        return v[static_cast<int>(component)];
    }
    case io::FieldKind::AtomSASA:
        return h5->sasa()->at(nativeRow, frame);
    case io::FieldKind::SASANormal: {
        const Vec3 v = h5->hydrationGeometryTimeSeries()->surfaceNormalAt(nativeRow, frame);
        return v[static_cast<int>(component)];
    }
    case io::FieldKind::WaterEfield: {
        const auto* ts = h5->waterFieldTimeSeries();
        return ts->efield[(nativeRow * ts->n_frames + frame) * 3 + component];
    }
    case io::FieldKind::WaterEfieldFirst: {
        const auto* ts = h5->waterFieldTimeSeries();
        return ts->efield_first[(nativeRow * ts->n_frames + frame) * 3 + component];
    }
    case io::FieldKind::WaterEFG: {
        const auto* ts = h5->waterFieldTimeSeries();
        return ts->efg[(nativeRow * ts->n_frames + frame) * 5 + component];
    }
    case io::FieldKind::WaterEFGFirst: {
        const auto* ts = h5->waterFieldTimeSeries();
        return ts->efg_first[(nativeRow * ts->n_frames + frame) * 5 + component];
    }
    case io::FieldKind::WaterShellCounts: {
        const auto* ts = h5->waterFieldTimeSeries();
        const std::size_t idx = nativeRow * ts->n_frames + frame;
        return component == 0 ? static_cast<double>(ts->n_first[idx])
                              : static_cast<double>(ts->n_second[idx]);
    }
    case io::FieldKind::HydrationShell: {
        const auto* ts = h5->hydrationShellTimeSeries();
        return (*hydrationShellChannel(*ts, component))[nativeRow * ts->n_frames + frame];
    }
    case io::FieldKind::WaterPolarization:
        return hydrationGeometryWaterPolarizationValue(*h5->hydrationGeometryTimeSeries(),
                                                       nativeRow,
                                                       frame,
                                                       component);
    case io::FieldKind::BondedEnergy: {
        const auto* ts = h5->bondedEnergy();
        return (*bondedEnergyChannel(*ts, component))[nativeRow * ts->n_frames + frame];
    }
    case io::FieldKind::GromacsEnergy:
        return gromacsComponentValue(*h5->gromacsEnergy(), frame, component);
    case io::FieldKind::DSSPSs8: {
        const auto* ts = h5->dssp8();
        const std::optional<std::size_t> residue = dsspResidueForAtom(run, ts, nativeRow);
        return dsspSs8Value(ts->codeAt(*residue, frame), component);
    }
    case io::FieldKind::DSSPHBondEnergy: {
        const auto* ts = h5->dssp8();
        const std::optional<std::size_t> residue = dsspResidueForAtom(run, ts, nativeRow);
        return dsspHBondEnergyValue(*ts, *residue, frame, component);
    }
    case io::FieldKind::DSSPChi: {
        const auto* ts = h5->dihedrals();
        const std::optional<std::size_t> residue = residueForAtom(run, nativeRow);
        return dsspChiValue(*ts, *residue, frame, component);
    }
    case io::FieldKind::OmegaActual:
        return h5->dihedrals()->omegaAt(nativeRow, frame);
    case io::FieldKind::OmegaDeviation:
        return h5->dihedrals()->omega_deviation[nativeRow * h5->dihedrals()->n_frames + frame];
    case io::FieldKind::OmegaIsXpro:
        return static_cast<double>(h5->dihedrals()->omega_is_xpro[nativeRow]);
    case io::FieldKind::AromaticChi2:
        return h5->ringPucker()->aromaticChi2At(nativeRow, frame);
    case io::FieldKind::PuckerQ:
        return h5->ringPucker()->puckerQAt(nativeRow, frame);
    case io::FieldKind::PuckerTheta:
        return h5->ringPucker()->pucker_theta[nativeRow * h5->ringPucker()->n_frames + frame];
    default:
        return 0.0;
    }
}

std::optional<double> topologyValue(const RunData& run, io::FieldKind kind, std::size_t row) {
    if (!run.protein) return std::nullopt;
    switch (kind) {
    case io::FieldKind::Element:
        if (row >= run.protein->atomCount()) return std::nullopt;
        return static_cast<double>(run.protein->atom(row).AtomicNumber());
    case io::FieldKind::ResidueIndex:
        if (row >= run.protein->atomCount()) return std::nullopt;
        return static_cast<double>(run.protein->atom(row).residueIndex);
    case io::FieldKind::ResidueType: {
        if (row >= run.protein->atomCount()) return std::nullopt;
        const int32_t residue = run.protein->atom(row).residueIndex;
        if (residue < 0 || static_cast<std::size_t>(residue) >= run.protein->residueCount())
            return std::nullopt;
        return static_cast<double>(static_cast<int>(run.protein->residue(static_cast<std::size_t>(residue)).aminoAcid));
    }
    case io::FieldKind::FfPartialCharge:
        if (row >= run.protein->atomCount() || !run.protein->atom(row).hasPartialCharge)
            return std::nullopt;
        return run.protein->atom(row).partialCharge;
    default:
        return std::nullopt;
    }
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

    field_specs_.resize(static_cast<std::size_t>(io::FieldKind::Count));
    for (const io::FieldSpec* field : ScopedProducerCatalog()) {
        field_specs_[static_cast<std::size_t>(field->kind)] = makeFieldAccessSpec(run, *field);
    }
}

const ArraySpec& Catalog::spec(ArrayId id) const {
    const int idx = ord(id);
    if (idx < 0 || idx >= static_cast<int>(specs_.size())) throw std::out_of_range("Catalog::spec");
    return specs_[static_cast<std::size_t>(idx)];
}

bool Catalog::has(ArrayId id) const { return spec(id).available; }

const FieldAccessSpec& Catalog::fieldSpec(io::FieldKind kind) const {
    const int idx = fieldOrd(kind);
    if (idx < 0 || idx >= static_cast<int>(field_specs_.size()))
        throw std::out_of_range("Catalog::fieldSpec");
    return field_specs_[static_cast<std::size_t>(idx)];
}

bool Catalog::has(io::FieldKind kind) const { return fieldSpec(kind).available; }

FieldPresence Catalog::present(const Body& body,
                               io::FieldKind kind,
                               std::size_t native_row,
                               std::size_t frame,
                               std::size_t component) const {
    const FieldAccessSpec& fs = fieldSpec(kind);
    if (!fs.available) {
        return {false, fs.absence_reason.isEmpty() ? QStringLiteral("not-produced-in-dataset")
                                                   : fs.absence_reason};
    }
    if (fs.structured) {
        if (native_row >= fs.native_rows) return {false, QStringLiteral("native-row-out-of-range")};
        if (frame >= fs.frames) return {false, QStringLiteral("frame-out-of-range")};
        return {true, QString()};
    }
    if (component >= fs.components) return {false, QStringLiteral("component-out-of-range")};
    if (native_row >= fs.native_rows) return {false, QStringLiteral("native-row-out-of-range")};
    if (frame >= fs.frames && fs.provider != FieldProvider::SparseDftByOriginal)
        return {false, QStringLiteral("frame-out-of-range")};

    switch (fs.provider) {
    case FieldProvider::StaticProducerArray: {
        const StaticNpyArray* a = body.run.producerArray(kind);
        if (!a) return {false, QStringLiteral("producer-array-absent")};
        if (component >= a->cols) return {false, QStringLiteral("component-out-of-range")};
        if (a->frameVarying && frame >= a->frameCount) return {false, QStringLiteral("frame-out-of-range")};
        if (a->frameVarying && native_row >= a->atomsPerFrame)
            return {false, QStringLiteral("native-row-out-of-range")};
        const std::size_t row = nativeRowForStatic(*a, native_row, frame);
        if (row >= a->rows) return {false, QStringLiteral("native-row-out-of-range")};
        return {true, QString()};
    }
    case FieldProvider::TypedTopology:
        if (!topologyValue(body.run, kind, native_row)) return {false, QStringLiteral("topology-mismatch")};
        return {true, QString()};
    case FieldProvider::SparseDftByOriginal:
        if (!body.run.protein || native_row >= body.run.protein->atomCount())
            return {false, QStringLiteral("native-row-out-of-range")};
        if (frame >= body.run.frameMap.frameCount()) return {false, QStringLiteral("frame-out-of-range")};
        if (!dftAt(body, native_row, frame)) return {false, QStringLiteral("frame-gap")};
        return component < 9 ? FieldPresence{true, QString()}
                             : FieldPresence{false, QStringLiteral("component-out-of-range")};
    case FieldProvider::DenseH5TimeSeries: {
        QString reason;
        const bool ok = denseFieldPresent(body.run, kind, native_row, frame, component, &reason);
        return {ok, ok ? QString()
                       : (reason.isEmpty() ? QStringLiteral("not-present") : reason)};
    }
    case FieldProvider::DatasetAbsent:
        return {false, fs.absence_reason.isEmpty() ? QStringLiteral("not-produced-in-dataset")
                                                   : fs.absence_reason};
    }
    return {false, QStringLiteral("unsupported-in-residence")};
}

std::optional<double> Catalog::value(const Body& body,
                                     io::FieldKind kind,
                                     std::size_t native_row,
                                     std::size_t frame,
                                     std::size_t component) const {
    const FieldPresence p = present(body, kind, native_row, frame, component);
    if (!p.present) return std::nullopt;
    const FieldAccessSpec& fs = fieldSpec(kind);
    switch (fs.provider) {
    case FieldProvider::StaticProducerArray: {
        const StaticNpyArray* a = body.run.producerArray(kind);
        if (!a) return std::nullopt;
        const std::size_t row = nativeRowForStatic(*a, native_row, frame);
        return a->value(row, component);
    }
    case FieldProvider::TypedTopology:
        return topologyValue(body.run, kind, native_row);
    case FieldProvider::SparseDftByOriginal: {
        const model::DftAtomShielding* dft = dftAt(body, native_row, frame);
        if (!dft) return std::nullopt;
        if (kind == io::FieldKind::OrcaTotal) return matComponent(dft->total_raw, static_cast<int>(component));
        if (kind == io::FieldKind::OrcaDiamagnetic) return matComponent(dft->dia_raw, static_cast<int>(component));
        if (kind == io::FieldKind::OrcaParamagnetic) return matComponent(dft->para_raw, static_cast<int>(component));
        return std::nullopt;
    }
    case FieldProvider::DenseH5TimeSeries:
        return denseFieldValue(body.run, kind, native_row, frame, component);
    case FieldProvider::DatasetAbsent:
        return std::nullopt;
    }
    return std::nullopt;
}

bool Catalog::present(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) const {
    if (!has(id)) return false;
    if (const StaticNpyArray* a = staticAt(body.run, id)) {
        if (id == ArrayId::Aimnet2Embedding) return atom < a->rows && a->cols == 256 && !a->floatValues.empty();
        return staticArrayPresent(a, atom, frame);
    }
    switch (id) {
    case ArrayId::Positions:
        return body.run.protein && body.run.conformation && atom < body.run.protein->atomCount()
               && frame < body.run.conformation->frameCount();
    case ArrayId::KernelBs:
        return shieldingPresent(body.run.h5() ? body.run.h5()->bsShielding() : nullptr, atom, frame);
    case ArrayId::KernelMc:
        return shieldingPresent(body.run.h5() ? body.run.h5()->mcShielding() : nullptr, atom, frame);
    case ArrayId::ApbsEfg:
        return body.run.h5() && body.run.h5()->apbsEfg()
               && atom < body.run.h5()->apbsEfg()->n_atoms
               && frame < body.run.h5()->apbsEfg()->n_frames
               && body.run.h5()->apbsEfg()->sourceAttachedAt(frame);
    case ArrayId::ApbsEfield:
        return body.run.h5() && body.run.h5()->apbsEfield()
               && atom < body.run.h5()->apbsEfield()->n_atoms
               && frame < body.run.h5()->apbsEfield()->n_frames
               && body.run.h5()->apbsEfield()->sourceAttachedAt(frame);
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
        return false;
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
