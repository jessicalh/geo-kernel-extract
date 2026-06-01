#include "ChargeDipoleNeighborhood.h"

#include "AnalysisBody.h"
#include "Catalog.h"
#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "SpatialIndexSet.h"
#include "TemporalIndex.h"
#include "TypedAtomIndex.h"

#include "../model/Conformation.h"
#include "../model/QtProtein.h"

#include <QLoggingCategory>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cChargeDipole, "h5reader.rediscover.charge_dipole")

LocalFrame buildHN(const Body& body, std::size_t resIdx, std::size_t hIdx,
                   std::size_t frameRow) {
    const model::QtProtein& p = *body.run.protein;
    const model::Conformation& conf = *body.run.conformation;
    const model::QtResidue& r = p.residue(resIdx);
    const std::vector<int32_t> scope = body.idx.typedAtoms.residueScope(resIdx);

    TypedAtomSelector nSel;
    nSel.backboneRole = model::BackboneRole::Nitrogen;
    TypedAtomSelector caSel;
    caSel.backboneRole = model::BackboneRole::AlphaCarbon;

    QString err;
    const std::optional<int32_t> nIdx = body.idx.typedAtoms.selectUnique(scope, nSel, &err);
    if (!nIdx) {
        throw std::runtime_error(QStringLiteral("charge_dipole HN frame failed for residue %1: N lookup %2")
                                     .arg(resIdx)
                                     .arg(err)
                                     .toStdString());
    }
    const std::optional<int32_t> caIdx = body.idx.typedAtoms.selectUnique(scope, caSel, &err);
    if (!caIdx) {
        throw std::runtime_error(QStringLiteral("charge_dipole HN frame failed for residue %1: CA lookup %2")
                                     .arg(resIdx)
                                     .arg(err)
                                     .toStdString());
    }

    bool cPrevValid = false;
    Vec3 cPrev = Vec3::Zero();
    if (r.prevResidueIndex >= 0
        && static_cast<std::size_t>(r.prevResidueIndex) < p.residueCount()) {
        const std::vector<int32_t> prevScope =
            body.idx.typedAtoms.residueScope(static_cast<std::size_t>(r.prevResidueIndex));
        TypedAtomSelector cSel;
        cSel.backboneRole = model::BackboneRole::CarbonylCarbon;
        const std::optional<int32_t> cIdx = body.idx.typedAtoms.selectUnique(prevScope, cSel, &err);
        if (cIdx) {
            cPrev = conf.atomPosition(frameRow, static_cast<std::size_t>(*cIdx));
            cPrevValid = true;
        }
    }

    const Vec3 nPos = conf.atomPosition(frameRow, static_cast<std::size_t>(*nIdx));
    const Vec3 hPos = conf.atomPosition(frameRow, hIdx);
    const Vec3 caPos = conf.atomPosition(frameRow, static_cast<std::size_t>(*caIdx));
    return BuildHNFrame(nPos, hPos, caPos, cPrev, cPrevValid);
}

bool validResidueIndex(const model::QtProtein& p, int32_t resIdx) {
    return resIdx >= 0 && static_cast<std::size_t>(resIdx) < p.residueCount();
}

ArrayId chargeArrayForSource(const QString& charge_source) {
    if (charge_source == QStringLiteral("ff14sb")) return ArrayId::Ff14sbCharge;
    if (charge_source == QStringLiteral("aimnet2")) return ArrayId::Aimnet2Charge;
    if (charge_source == QStringLiteral("mopac")) return ArrayId::MopacCharge;
    return ArrayId::MopacCharge;
}

QString chargeSourceMissingMessage(const QString& charge_source) {
    if (charge_source == QStringLiteral("ff14sb"))
        return QStringLiteral("charge_dipole requires FF14SB charges, but topol.top charges were not loaded");
    if (charge_source == QStringLiteral("aimnet2"))
        return QStringLiteral("charge_dipole requires AIMNet2 Hirshfeld charges, but aimnet2_charge is absent");
    if (charge_source == QStringLiteral("mopac"))
        return QStringLiteral("charge_dipole requires charge_source=mopac, but per-frame MOPAC charges are absent");
    return QStringLiteral("charge_dipole received unknown charge_source=%1").arg(charge_source);
}

}  // namespace

FeatureSchema ChargeDipoleNeighborhood::schema() const {
    FeatureSchema s;
    s.caseName = name();
    s.relationshipKind = RelationshipKind::SourceSum;
    s.sourceSchemaKind = SourceSchemaKind::Charge;
    s.includeBareKernel = false;

    s.sourceColumns = IdentityColumns();
    const std::vector<FeatureColumn> srcGeom = {
        {QStringLiteral("source_kind"), {}},
        {QStringLiteral("charge_source"), {}},
        {QStringLiteral("source_atom_index"), {}},
        {QStringLiteral("source_residue_index"), {}},
        {QStringLiteral("source_residue_number"), {}},
        {QStringLiteral("source_amino_acid_ord"), {}},
        {QStringLiteral("source_element_ord"), {}},
        {QStringLiteral("source_atom_name"), {}},
        {QStringLiteral("source_q_e"), QStringLiteral("e")},
        {QStringLiteral("disp_local_x"), QStringLiteral("Angstrom")},
        {QStringLiteral("disp_local_y"), QStringLiteral("Angstrom")},
        {QStringLiteral("disp_local_z"), QStringLiteral("Angstrom")},
        {QStringLiteral("r"), QStringLiteral("Angstrom")},
    };
    for (const auto& c : srcGeom) s.sourceColumns.push_back(c);
    for (const auto& c : TargetColumns()) s.sourceColumns.push_back(c);

    s.aggregatedColumns = IdentityColumns();
    s.aggregatedColumns.push_back({QStringLiteral("n_sources"), {}});
    s.aggregatedColumns.push_back({QStringLiteral("charge_source"), {}});
    s.aggregatedColumns.push_back({QStringLiteral("exclude_residue"), {}});
    s.aggregatedColumns.push_back({QStringLiteral("cutoff_A"), QStringLiteral("Angstrom")});
    s.aggregatedColumns.push_back({QStringLiteral("mu_x"), QStringLiteral("e*Angstrom")});
    s.aggregatedColumns.push_back({QStringLiteral("mu_y"), QStringLiteral("e*Angstrom")});
    s.aggregatedColumns.push_back({QStringLiteral("mu_z"), QStringLiteral("e*Angstrom")});
    s.aggregatedColumns.push_back({QStringLiteral("mu_norm"), QStringLiteral("e*Angstrom")});
    for (const auto& c : TargetColumns()) s.aggregatedColumns.push_back(c);

    return s;
}

std::size_t ChargeDipoleNeighborhood::extract(const Body& body, RecordSink& sink) const {
    const ArrayId chargeArray = chargeArrayForSource(charge_source);
    if (!body.catalog.has(chargeArray)) {
        throw std::runtime_error(chargeSourceMissingMessage(charge_source).toStdString());
    }

    const RunData& run = body.run;
    const model::QtProtein& p = *run.protein;
    const model::Conformation& conf = *run.conformation;

    const auto stratum = [&]() {
        std::vector<std::size_t> atoms;
        for (std::size_t i = 0; i < p.atomCount(); ++i)
            if (p.atom(i).IsBackboneAmideHydrogen()) atoms.push_back(i);
        return atoms;
    }();

    const auto localFrame = [&](std::size_t atomIdx, std::size_t frameRow) {
        const model::QtAtom& a = p.atom(atomIdx);
        if (!validResidueIndex(p, a.residueIndex)) return LocalFrame{};
        return buildHN(body, static_cast<std::size_t>(a.residueIndex), atomIdx, frameRow);
    };

    const auto selectSources = [&](std::size_t frameRow, const Vec3& atomPos) {
        return body.idx.spatial.near(CloudKind::ChargeSites, frameRow, atomPos, cutoff_A);
    };

    const auto attachSource = [&](std::size_t targetAtom, std::size_t frameRow,
                                  const LocalFrame& frame, const Vec3& atomPos,
                                  const SourceRef& ref, SourceSlot& s) -> bool {
        if (ref.entity_index < 0) return false;
        const std::size_t srcAtom = static_cast<std::size_t>(ref.entity_index);
        if (srcAtom >= p.atomCount()) return false;
        if (!body.catalog.present(body, chargeArray, srcAtom, frameRow)) return false;

        const model::QtAtom& target = p.atom(targetAtom);
        const model::QtAtom& src = p.atom(srcAtom);
        if (exclude_residue && target.residueIndex >= 0 && src.residueIndex == target.residueIndex)
            return false;

        const double q = body.catalog.value(body, chargeArray, srcAtom, frameRow);
        if (!std::isfinite(q)) return false;

        const Vec3 dispLab = conf.atomPosition(frameRow, srcAtom) - atomPos;
        const double r = dispLab.norm();
        if (!(r > 1e-9)) return false;

        s.kind = SourceKind::Charge;
        s.source_charge_source = charge_source;
        s.source_q_e = q;
        s.source_atom_index = src.atomIndex;
        s.source_residue_index = src.residueIndex;
        s.source_element = static_cast<int>(src.element);
        s.source_atom_name = p.atomLabel(srcAtom, model::NamingConvention::Bmrb);
        if (validResidueIndex(p, src.residueIndex)) {
            const model::QtResidue& sr = p.residue(static_cast<std::size_t>(src.residueIndex));
            s.source_residue_number = sr.address.residueNumber;
            s.source_amino_acid = static_cast<int>(sr.aminoAcid);
        }
        s.disp_local = frame.is_valid ? frame.ToLocal(dispLab) : dispLab;
        s.r = r;
        return true;
    };

    const auto reduceMu = [](const std::vector<SourceSlot>& sources) {
        Vec3 mu = Vec3::Zero();
        for (const SourceSlot& s : sources) mu += s.source_q_e * s.disp_local;
        return mu;
    };

    qCInfo(cChargeDipole).noquote()
        << "charge_dipole stratum HN atoms=" << stratum.size()
        << "| charge_source=" << charge_source << "| cutoff_A=" << cutoff_A
        << "| exclude_residue=" << exclude_residue
        << "| dft rows=" << run.frameMap.dftRows().size();

    std::size_t cases = 0;
    for (std::size_t row : run.frameMap.dftRows()) {
        for (std::size_t atomIdx : stratum) {
            const FrameWindow instant = body.idx.temporal.range(atomIdx, row, 0, 0);
            if (instant.size() != 1 || !instant.contains(row)) {
                throw std::runtime_error("charge_dipole temporal index failed to produce the instantaneous frame");
            }
            const std::size_t frameRow = row;
            const std::size_t orig = run.frameMap.originalIndex(frameRow);
            const LocalFrame frame = localFrame(atomIdx, frameRow);

            NeighborhoodRecord rec;
            FillIdentity(rec, run, atomIdx, frameRow, name(), frame);
            rec.target = BuildTarget(run, atomIdx, orig, frame);

            const Vec3 atomPos = conf.atomPosition(frameRow, atomIdx);
            for (const SourceRef& ref : selectSources(frameRow, atomPos)) {
                SourceSlot s;
                if (attachSource(atomIdx, frameRow, frame, atomPos, ref, s)) {
                    rec.sources.push_back(std::move(s));
                }
            }

            const Vec3 mu = reduceMu(rec.sources);
            sink.WriteSourceRows(rec);
            sink.WriteChargeDipoleAggregatedRow(rec, mu, cutoff_A, charge_source, exclude_residue);
            ++cases;
        }
    }

    qCInfo(cChargeDipole).noquote() << "charge_dipole cases=" << cases;
    return cases;
}

}  // namespace h5reader::rediscover
