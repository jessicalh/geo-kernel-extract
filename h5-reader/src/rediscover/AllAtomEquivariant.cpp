#include "AllAtomEquivariant.h"

#include "Catalog.h"
#include "ExtractionSupport.h"
#include "Relationship.h"
#include "SphericalBasis.h"
#include "SpatialIndexSet.h"
#include "Verbs.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtRing.h"
#include "../model/QtTimeSeriesBuffers.h"
#include "../model/QtTopology.h"

#include <QLoggingCategory>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cAllAtomEquivariant, "h5reader.rediscover.all_atom_equivariant")

const double kNaN = std::numeric_limits<double>::quiet_NaN();

bool validResidue(const model::QtProtein& p, int32_t r) {
    return r >= 0 && static_cast<std::size_t>(r) < p.residueCount();
}

bool finiteVec3(const Vec3& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

Vec3 unitOrZero(const Vec3& v) {
    const double n = v.norm();
    if (n > 1e-12 && std::isfinite(n)) return v / n;
    return Vec3::Zero();
}

double invR3(double r) {
    return r > 1e-12 ? 1.0 / (r * r * r) : kNaN;
}

double dipolarFrom(double cosTheta, double r) {
    return r > 1e-12 && std::isfinite(cosTheta)
               ? (3.0 * cosTheta * cosTheta - 1.0) / (r * r * r)
               : kNaN;
}

double t2Magnitude(const std::array<double, 5>& t2) {
    double s = 0.0;
    for (double v : t2) s += v * v;
    return std::sqrt(s);
}

QString bondCategoryName(model::BondCategory c) {
    switch (c) {
    case model::BondCategory::PeptideCO: return QStringLiteral("PeptideCO");
    case model::BondCategory::PeptideCN: return QStringLiteral("PeptideCN");
    case model::BondCategory::BackboneOther: return QStringLiteral("BackboneOther");
    case model::BondCategory::SidechainCO: return QStringLiteral("SidechainCO");
    case model::BondCategory::Aromatic: return QStringLiteral("Aromatic");
    case model::BondCategory::Disulfide: return QStringLiteral("Disulfide");
    case model::BondCategory::SidechainOther: return QStringLiteral("SidechainOther");
    case model::BondCategory::Unknown: return QStringLiteral("Unknown");
    }
    return QStringLiteral("Unknown");
}

QString ringTypeName(model::RingTypeIndex t) {
    switch (t) {
    case model::RingTypeIndex::PheBenzene: return QStringLiteral("PheBenzene");
    case model::RingTypeIndex::TyrPhenol: return QStringLiteral("TyrPhenol");
    case model::RingTypeIndex::TrpBenzene: return QStringLiteral("TrpBenzene");
    case model::RingTypeIndex::TrpPyrrole: return QStringLiteral("TrpPyrrole");
    case model::RingTypeIndex::TrpPerimeter: return QStringLiteral("TrpPerimeter");
    case model::RingTypeIndex::HisImidazole: return QStringLiteral("HisImidazole");
    case model::RingTypeIndex::HidImidazole: return QStringLiteral("HidImidazole");
    case model::RingTypeIndex::HieImidazole: return QStringLiteral("HieImidazole");
    case model::RingTypeIndex::ProPyrrolidine: return QStringLiteral("ProPyrrolidine");
    }
    return QStringLiteral("UnknownRing");
}

std::size_t heavyParent(const model::QtProtein& p, std::size_t atomIdx) {
    const model::QtAtom& a = p.atom(atomIdx);
    return a.parentAtomIndex >= 0 ? static_cast<std::size_t>(a.parentAtomIndex) : atomIdx;
}

bool ringContainsAtom(const model::QtRing& ring, int32_t atomIdx) {
    return std::find(ring.atomIndices.begin(), ring.atomIndices.end(), atomIdx)
           != ring.atomIndices.end();
}

void fillTargetIdentity(const Body& body, std::size_t atom, std::size_t row,
                        AllAtomEquivariantTargetRecord& out) {
    const model::QtProtein& p = *body.run.protein;
    const model::QtAtom& a = p.atom(atom);
    out.atom_index = a.atomIndex;
    out.residue_index = a.residueIndex;
    out.element = static_cast<int>(a.element);
    out.atom_name = p.atomLabel(atom, model::NamingConvention::Bmrb);
    if (validResidue(p, a.residueIndex)) {
        const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
        out.residue_number = r.address.residueNumber;
        out.amino_acid = static_cast<int>(r.aminoAcid);
    }
    out.h5_row = static_cast<int32_t>(row);
    out.original_index = static_cast<int32_t>(body.run.frameMap.originalIndex(row));
    out.time_ps = body.run.trajectory()->timePicoseconds(row);
}

void fillSourceTargetIdentity(const AllAtomEquivariantTargetRecord& target,
                              AllAtomEquivariantSourceRecord& out) {
    out.target_atom_index = target.atom_index;
    out.target_residue_index = target.residue_index;
    out.h5_row = target.h5_row;
    out.original_index = target.original_index;
    out.time_ps = target.time_ps;
    out.r = kNaN;
    out.inv_r3 = kNaN;
    out.cos_theta = kNaN;
    out.dipolar = kNaN;
    out.source_value = kNaN;
    out.source_value_2 = kNaN;
    out.value_vec_mag = kNaN;
    out.source_q_e = kNaN;
    out.q_over_r3 = kNaN;
}

void fillSourceAtomIdentity(const model::QtProtein& p, std::size_t srcAtom,
                            AllAtomEquivariantSourceRecord& out) {
    const model::QtAtom& a = p.atom(srcAtom);
    out.source_atom_index = a.atomIndex;
    out.source_residue_index = a.residueIndex;
    out.source_element = static_cast<int>(a.element);
    out.source_atom_name = p.atomLabel(srcAtom, model::NamingConvention::Bmrb);
    if (validResidue(p, a.residueIndex)) {
        const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
        out.source_residue_number = r.address.residueNumber;
    }
}

bool apbsEfieldPresent(const Body& body, std::size_t atom, std::size_t row) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtVec3TimeSeries* efield = h5 ? h5->apbsEfield() : nullptr;
    return efield && body.catalog.present(body, ArrayId::ApbsEfield, atom, row)
           && efield->sourceAttachedAt(row);
}

bool apbsEfgPresent(const Body& body, std::size_t atom, std::size_t row) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtT2TimeSeries* efg = h5 ? h5->apbsEfg() : nullptr;
    return efg && body.catalog.present(body, ArrayId::ApbsEfg, atom, row)
           && efg->sourceAttachedAt(row);
}

bool aimnetChargePresent(const Body& body, std::size_t atom, std::size_t row) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtScalarTimeSeries* charge = h5 ? h5->aimnet2Charge() : nullptr;
    return charge && body.catalog.present(body, ArrayId::Aimnet2Charge, atom, row)
           && charge->sourceAttachedAt(row);
}

bool aimnetCrgPresent(const Body& body, std::size_t atom, std::size_t row) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtAimnet2ChargeResponseGradientTimeSeries* crg =
        h5 ? h5->aimnet2ChargeResponseGradient() : nullptr;
    return crg && body.catalog.present(body, ArrayId::Aimnet2ChargeRespScalar, atom, row)
           && body.catalog.present(body, ArrayId::Aimnet2ChargeRespVector, atom, row)
           && crg->meta.sourceAttachedAt(row);
}

bool aimnetEmbeddingPresent(const Body& body, std::size_t atom, std::size_t row,
                            std::size_t& dims, const float*& ptr) {
    dims = 0;
    ptr = body.catalog.valueEmbedding(body, ArrayId::Aimnet2Embedding, atom, row, dims);
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtEmbeddingTimeSeries* emb = h5 ? h5->aimnet2Embedding() : nullptr;
    return ptr && emb && body.catalog.present(body, ArrayId::Aimnet2Embedding, atom, row)
           && emb->meta.sourceAttachedAt(row);
}

AllAtomEquivariantSourceRecord makeRingSource(const Body& body,
                                              const AllAtomEquivariantTargetRecord& target,
                                              std::size_t targetAtom,
                                              std::size_t ringIdx,
                                              std::size_t row) {
    AllAtomEquivariantSourceRecord out;
    fillSourceTargetIdentity(target, out);
    const model::QtProtein& p = *body.run.protein;
    const model::QtTopology& topo = p.topology();
    const model::QtRing& ring = topo.ringAt(ringIdx);
    const model::RingGeometry& g = verbs::ringGeom(body, ringIdx, row);
    const Vec3 targetPos = verbs::pos(body, targetAtom, row);
    const Vec3 disp = g.center - targetPos;
    const double r = disp.norm();
    const Vec3 nrm = unitOrZero(g.normal);
    const double cosT = r > 1e-12 ? disp.dot(nrm) / r : kNaN;

    out.mechanism = QStringLiteral("ring");
    out.source_kind = QStringLiteral("ring_center");
    out.category = ringTypeName(ring.TypeIndex());
    out.category_ord = ring.TypeIndexAsInt();
    out.source_id = static_cast<int32_t>(ringIdx);
    out.disp = disp;
    out.r = r;
    out.inv_r3 = invR3(r);
    out.cos_theta = cosT;
    out.dipolar = dipolarFrom(cosT, r);
    out.orientation_a = nrm;
    out.source_value = ring.LiteratureIntensity();
    out.source_units = QStringLiteral("nA/T");
    out.ring_index = static_cast<int32_t>(ringIdx);
    out.ring_type_index = ring.TypeIndexAsInt();
    out.ring_aromaticity = static_cast<int>(ring.Aromaticity());
    out.ring_size = ring.RingSizeValue();
    out.ring_nitrogen = ring.NitrogenCount();
    out.ring_fused = ring.IsFused();
    out.ring_intensity = ring.LiteratureIntensity();
    out.ring_jb_offset = ring.JohnsonBoveyLobeOffset();
    const int32_t targetHeavy = static_cast<int32_t>(heavyParent(p, targetAtom));
    out.source_is_self_or_bonded = ringContainsAtom(ring, targetHeavy);
    return out;
}

AllAtomEquivariantSourceRecord makeBondSource(const Body& body,
                                              const AllAtomEquivariantTargetRecord& target,
                                              std::size_t targetAtom,
                                              std::size_t bondTopoIdx,
                                              std::size_t row,
                                              double nearFieldRatio) {
    AllAtomEquivariantSourceRecord out;
    fillSourceTargetIdentity(target, out);
    const model::QtProtein& p = *body.run.protein;
    const model::QtBond& b = p.topology().bondAt(bondTopoIdx);
    const Vec3 posA = verbs::pos(body, static_cast<std::size_t>(b.atomIndexA), row);
    const Vec3 posB = verbs::pos(body, static_cast<std::size_t>(b.atomIndexB), row);
    const Vec3 axis = posB - posA;
    const double axisNorm = axis.norm();
    const Vec3 axisU = unitOrZero(axis);
    const Vec3 targetPos = verbs::pos(body, targetAtom, row);
    const Vec3 disp = 0.5 * (posA + posB) - targetPos;
    const double r = disp.norm();
    const double cosT = r > 1e-12 ? disp.dot(axisU) / r : kNaN;

    out.mechanism = QStringLiteral("bond");
    out.source_kind = QStringLiteral("bond_midpoint");
    out.category = bondCategoryName(b.category);
    out.category_ord = static_cast<int>(b.category);
    out.source_id = b.bondIndex;
    out.disp = disp;
    out.r = r;
    out.inv_r3 = invR3(r);
    out.cos_theta = cosT;
    out.dipolar = dipolarFrom(cosT, r);
    out.orientation_a = axisU;
    out.source_value = axisNorm;
    out.source_units = QStringLiteral("Angstrom");
    out.bond_index = b.bondIndex;
    out.bond_category = static_cast<int>(b.category);
    out.bond_order = static_cast<int>(b.order);
    out.bond_atom_a = b.atomIndexA;
    out.bond_atom_b = b.atomIndexB;
    out.bond_elem_a = static_cast<int>(p.atom(static_cast<std::size_t>(b.atomIndexA)).element);
    out.bond_elem_b = static_cast<int>(p.atom(static_cast<std::size_t>(b.atomIndexB)).element);
    out.bond_length = axisNorm;
    const bool endpointSelf = static_cast<int32_t>(targetAtom) == b.atomIndexA
                              || static_cast<int32_t>(targetAtom) == b.atomIndexB;
    const bool nearField = r <= axisNorm * nearFieldRatio;
    out.source_is_self_or_bonded = endpointSelf || nearField;
    return out;
}

AllAtomEquivariantSourceRecord makeChargeSource(const Body& body,
                                                const AllAtomEquivariantTargetRecord& target,
                                                std::size_t targetAtom,
                                                std::size_t sourceAtom,
                                                std::size_t row,
                                                ArrayId chargeArray,
                                                const QString& chargeName) {
    AllAtomEquivariantSourceRecord out;
    fillSourceTargetIdentity(target, out);
    const model::QtProtein& p = *body.run.protein;
    fillSourceAtomIdentity(p, sourceAtom, out);
    const Vec3 disp = verbs::pos(body, sourceAtom, row) - verbs::pos(body, targetAtom, row);
    const double r = disp.norm();
    const double inv = invR3(r);
    const double q = body.catalog.value(body, chargeArray, sourceAtom, row);

    out.mechanism = QStringLiteral("charge");
    out.source_kind = QStringLiteral("%1_charge_site").arg(chargeName);
    out.category = chargeName;
    out.category_ord = chargeName == QStringLiteral("ff14sb") ? 0 : 1;
    out.source_id = static_cast<int32_t>(sourceAtom);
    out.disp = disp;
    out.r = r;
    out.inv_r3 = inv;
    out.orientation_a = unitOrZero(disp);
    out.source_value = q;
    out.source_value_2 = std::isfinite(inv) ? q * inv : kNaN;
    out.source_units = QStringLiteral("e");
    out.charge_source = chargeName;
    out.source_q_e = q;
    out.q_over_r3 = out.source_value_2;
    return out;
}

AllAtomEquivariantSourceRecord makeVectorFeatureSource(
    const AllAtomEquivariantTargetRecord& target, const QString& mechanism,
    const QString& sourceKind, const QString& category, const Vec3& v, const QString& units) {
    AllAtomEquivariantSourceRecord out;
    fillSourceTargetIdentity(target, out);
    out.mechanism = mechanism;
    out.source_kind = sourceKind;
    out.category = category;
    out.source_id = target.atom_index;
    out.source_atom_index = target.atom_index;
    out.source_residue_index = target.residue_index;
    out.source_residue_number = target.residue_number;
    out.source_element = target.element;
    out.source_atom_name = target.atom_name;
    out.value_vec = v;
    out.value_vec_mag = v.norm();
    out.orientation_a = unitOrZero(v);
    out.source_value = out.value_vec_mag;
    out.source_units = units;
    return out;
}

AllAtomEquivariantSourceRecord makeTensorFeatureSource(
    const AllAtomEquivariantTargetRecord& target, const QString& mechanism,
    const QString& sourceKind, const QString& category, const std::array<double, 5>& t2,
    const QString& units) {
    AllAtomEquivariantSourceRecord out;
    fillSourceTargetIdentity(target, out);
    out.mechanism = mechanism;
    out.source_kind = sourceKind;
    out.category = category;
    out.source_id = target.atom_index;
    out.source_atom_index = target.atom_index;
    out.source_residue_index = target.residue_index;
    out.source_residue_number = target.residue_number;
    out.source_element = target.element;
    out.source_atom_name = target.atom_name;
    out.tensor_T2 = t2;
    out.tensor_present = true;
    out.source_value = t2Magnitude(t2);
    out.source_units = units;
    return out;
}

}  // namespace

AllAtomEquivariantStats RunAllAtomEquivariantEmit(const Body& body,
                                                  AllAtomEquivariantSink& sink,
                                                  const AllAtomEquivariantConfig& config) {
    AllAtomEquivariantStats stats;
    if (!body.run.protein || !body.run.trajectory()) return stats;

    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtEmbeddingTimeSeries* emb = h5 ? h5->aimnet2Embedding() : nullptr;
    if (emb && emb->n_dims != sink.embeddingDims()) {
        throw std::runtime_error(QStringLiteral("aimnet2_embedding dimension %1 != expected %2")
                                     .arg(static_cast<qulonglong>(emb->n_dims))
                                     .arg(static_cast<qulonglong>(sink.embeddingDims()))
                                     .toStdString());
    }

    const model::QtProtein& p = *body.run.protein;
    stats.atom_count = p.atomCount();
    stats.dft_rows = body.run.frameMap.dftRows().size();

    qCInfo(cAllAtomEquivariant).noquote()
        << "all_atom_equivariant | atoms=" << stats.atom_count
        << "| dft rows=" << stats.dft_rows
        << "| ring_cut=" << config.ring_cutoff_A
        << "| bond_cut=" << config.bond_cutoff_A
        << "| charge_cut=" << config.charge_cutoff_A
        << "| frame=molecular_lab_h5_orca_aligned"
        << "| local_frame=none"
        << "| ff14sb_charge=" << body.catalog.has(ArrayId::Ff14sbCharge)
        << "| aimnet2_charge=" << body.catalog.has(ArrayId::Aimnet2Charge)
        << "| apbs_efield=" << body.catalog.has(ArrayId::ApbsEfield)
        << "| apbs_efg=" << body.catalog.has(ArrayId::ApbsEfg)
        << "| aimnet2_embedding_dims=" << (emb ? emb->n_dims : 0);

    LocalFrame noLocalFrame;
    for (std::size_t row : body.run.frameMap.dftRows()) {
        const std::size_t orig = body.run.frameMap.originalIndex(row);
        for (std::size_t atom = 0; atom < p.atomCount(); ++atom) {
            AllAtomEquivariantTargetRecord target;
            fillTargetIdentity(body, atom, row, target);
            target.target = BuildTarget(body.run, atom, orig, noLocalFrame);

            target.apbs_efield_present = apbsEfieldPresent(body, atom, row);
            if (target.apbs_efield_present)
                target.apbs_efield_lab = body.catalog.valueVec3(body, ArrayId::ApbsEfield, atom, row);
            target.apbs_efg_present = apbsEfgPresent(body, atom, row);
            if (target.apbs_efg_present)
                target.apbs_efg_T2 = body.catalog.valueT2(body, ArrayId::ApbsEfg, atom, row);
            target.aimnet2_charge_present = aimnetChargePresent(body, atom, row);
            if (target.aimnet2_charge_present)
                target.aimnet2_charge = body.catalog.value(body, ArrayId::Aimnet2Charge, atom, row);
            target.aimnet2_crg_present = aimnetCrgPresent(body, atom, row);
            if (target.aimnet2_crg_present) {
                target.aimnet2_crg_scalar =
                    body.catalog.value(body, ArrayId::Aimnet2ChargeRespScalar, atom, row);
                target.aimnet2_crg_lab =
                    body.catalog.valueVec3(body, ArrayId::Aimnet2ChargeRespVector, atom, row);
            }
            const float* embeddingPtr = nullptr;
            std::size_t embeddingDims = 0;
            target.aimnet2_embedding_present =
                aimnetEmbeddingPresent(body, atom, row, embeddingDims, embeddingPtr);
            target.aimnet2_embedding = embeddingPtr;
            target.aimnet2_embedding_dims = embeddingDims;
            if (target.aimnet2_embedding_present && embeddingDims == sink.embeddingDims())
                ++stats.aimnet2_embedding_present;

            const int64_t rowId = sink.WriteTarget(target);
            ++stats.target_rows;
            if (target.target.present) ++stats.dft_present;

            for (const SourceRef& ref : verbs::near(body, CloudKind::RingCenters, atom, row,
                                                    config.ring_cutoff_A)) {
                if (ref.entity_index < 0) continue;
                const std::size_t ringIdx = static_cast<std::size_t>(ref.entity_index);
                if (ringIdx >= p.topology().ringCount()) continue;
                AllAtomEquivariantSourceRecord src =
                    makeRingSource(body, target, atom, ringIdx, row);
                sink.WriteSource(rowId, src);
                ++stats.source_rows;
                ++stats.ring_rows;
                if (src.ring_type_index >= 0
                    && src.ring_type_index < static_cast<int>(stats.ring_type_rows.size()))
                    ++stats.ring_type_rows[static_cast<std::size_t>(src.ring_type_index)];
            }

            for (const SourceRef& ref : verbs::near(body, CloudKind::AllBondMidpoints, atom, row,
                                                    config.bond_cutoff_A)) {
                if (ref.entity_index < 0) continue;
                const std::size_t bondIdx = static_cast<std::size_t>(ref.entity_index);
                if (bondIdx >= p.topology().bondCount()) continue;
                AllAtomEquivariantSourceRecord src =
                    makeBondSource(body, target, atom, bondIdx, row, config.mc_near_field_ratio);
                sink.WriteSource(rowId, src);
                ++stats.source_rows;
                ++stats.bond_rows;
                if (src.bond_category >= 0
                    && src.bond_category < static_cast<int>(stats.bond_category_rows.size()))
                    ++stats.bond_category_rows[static_cast<std::size_t>(src.bond_category)];
            }

            if (body.catalog.has(ArrayId::Ff14sbCharge)) {
                for (const SourceRef& ref : verbs::near(body, CloudKind::ChargeSites, atom, row,
                                                        config.charge_cutoff_A)) {
                    if (ref.entity_index < 0) continue;
                    const std::size_t srcAtom = static_cast<std::size_t>(ref.entity_index);
                    if (srcAtom >= p.atomCount() || srcAtom == atom) continue;
                    if (!body.catalog.present(body, ArrayId::Ff14sbCharge, srcAtom, row)) continue;
                    AllAtomEquivariantSourceRecord src =
                        makeChargeSource(body, target, atom, srcAtom, row, ArrayId::Ff14sbCharge,
                                         QStringLiteral("ff14sb"));
                    if (!(src.r > 1e-12) || !std::isfinite(src.source_q_e)) continue;
                    sink.WriteSource(rowId, src);
                    ++stats.source_rows;
                    ++stats.charge_ff14sb_rows;
                }
            }

            if (body.catalog.has(ArrayId::Aimnet2Charge)) {
                for (const SourceRef& ref : verbs::near(body, CloudKind::Atoms, atom, row,
                                                        config.charge_cutoff_A)) {
                    if (ref.entity_index < 0) continue;
                    const std::size_t srcAtom = static_cast<std::size_t>(ref.entity_index);
                    if (srcAtom >= p.atomCount() || srcAtom == atom) continue;
                    if (!body.catalog.present(body, ArrayId::Aimnet2Charge, srcAtom, row)) continue;
                    AllAtomEquivariantSourceRecord src =
                        makeChargeSource(body, target, atom, srcAtom, row, ArrayId::Aimnet2Charge,
                                         QStringLiteral("aimnet2"));
                    if (!(src.r > 1e-12) || !std::isfinite(src.source_q_e)) continue;
                    sink.WriteSource(rowId, src);
                    ++stats.source_rows;
                    ++stats.charge_aimnet2_rows;
                }
            }

            if (target.apbs_efield_present && finiteVec3(target.apbs_efield_lab)) {
                AllAtomEquivariantSourceRecord src = makeVectorFeatureSource(
                    target, QStringLiteral("field"), QStringLiteral("apbs_efield"),
                    QStringLiteral("buckingham_efield"), target.apbs_efield_lab,
                    QStringLiteral("V/Angstrom"));
                sink.WriteSource(rowId, src);
                ++stats.source_rows;
                ++stats.apbs_efield_rows;
            }

            if (target.apbs_efg_present) {
                AllAtomEquivariantSourceRecord src = makeTensorFeatureSource(
                    target, QStringLiteral("efg"), QStringLiteral("apbs_efg"),
                    QStringLiteral("apbs_efg"), target.apbs_efg_T2,
                    QStringLiteral("V/Angstrom^2"));
                sink.WriteSource(rowId, src);
                ++stats.source_rows;
                ++stats.apbs_efg_rows;
            }

            if ((target.aimnet2_charge_present && std::isfinite(target.aimnet2_charge))
                || (target.aimnet2_crg_present && finiteVec3(target.aimnet2_crg_lab))
                || target.aimnet2_embedding_present) {
                AllAtomEquivariantSourceRecord src = makeVectorFeatureSource(
                    target, QStringLiteral("aimnet2"), QStringLiteral("aimnet2_atom_feature"),
                    QStringLiteral("aimnet2"), target.aimnet2_crg_lab, QStringLiteral("e2/Angstrom"));
                src.source_value = target.aimnet2_charge_present ? target.aimnet2_charge : kNaN;
                src.source_value_2 =
                    target.aimnet2_crg_present ? target.aimnet2_crg_scalar : kNaN;
                src.charge_source = QStringLiteral("aimnet2");
                src.source_q_e = target.aimnet2_charge_present ? target.aimnet2_charge : kNaN;
                src.aimnet2_embedding_present = target.aimnet2_embedding_present;
                src.aimnet2_embedding_dims = target.aimnet2_embedding_dims;
                sink.WriteSource(rowId, src);
                ++stats.source_rows;
                ++stats.aimnet2_atom_rows;
            }
        }
    }

    qCInfo(cAllAtomEquivariant).noquote()
        << "all_atom_equivariant rows | targets=" << stats.target_rows
        << "| dft_present=" << stats.dft_present
        << "| sources=" << stats.source_rows
        << "| rings=" << stats.ring_rows
        << "| bonds=" << stats.bond_rows
        << "| charge_ff14sb=" << stats.charge_ff14sb_rows
        << "| charge_aimnet2=" << stats.charge_aimnet2_rows
        << "| apbs_E=" << stats.apbs_efield_rows
        << "| apbs_EFG=" << stats.apbs_efg_rows
        << "| aimnet2_atom=" << stats.aimnet2_atom_rows;
    return stats;
}

}  // namespace h5reader::rediscover
