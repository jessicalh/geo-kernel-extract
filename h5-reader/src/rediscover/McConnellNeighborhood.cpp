#include "McConnellNeighborhood.h"

#include "AnalysisBody.h"
#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "McConnellLiteratureKernel.h"
#include "SpatialIndexSet.h"

#include "../model/QtBond.h"
#include "../model/QtProtein.h"
#include "../io/QtTrajectoryH5.h"

#include <QLoggingCategory>

#include <cmath>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cMc, "h5reader.rediscover.mc")

// The four anisotropic categories, in a fixed order for the per-category sum
// columns. Ordinals from model::BondCategory.
struct CatCol { model::BondCategory cat; const char* tag; };
const CatCol kCatCols[] = {
    {model::BondCategory::PeptideCO, "peptide_co"},
    {model::BondCategory::PeptideCN, "peptide_cn"},
    {model::BondCategory::SidechainCO, "sidechain_co"},
    {model::BondCategory::Aromatic, "aromatic"},
};
constexpr int kCatCount = 4;

int catColumn(int categoryOrd) {
    for (int i = 0; i < kCatCount; ++i)
        if (static_cast<int>(kCatCols[i].cat) == categoryOrd) return i;
    return -1;
}

// Build the HN amide-plane frame for the HN atom `hIdx` in residue `resIdx`.
LocalFrame buildHN(const model::QtProtein& p, const model::Conformation& conf,
                   std::size_t resIdx, std::size_t hIdx, std::size_t frame) {
    const model::QtResidue& r = p.residue(resIdx);
    if (r.N == model::QtResidue::NONE || r.CA == model::QtResidue::NONE)
        return {};  // invalid frame
    const Vec3 nPos = conf.atomPosition(frame, static_cast<std::size_t>(r.N));
    const Vec3 hPos = conf.atomPosition(frame, hIdx);
    const Vec3 caPos = conf.atomPosition(frame, static_cast<std::size_t>(r.CA));

    bool cPrevValid = false;
    Vec3 cPrev = Vec3::Zero();
    if (r.prevResidueIndex >= 0
        && static_cast<std::size_t>(r.prevResidueIndex) < p.residueCount()) {
        const model::QtResidue& prev = p.residue(static_cast<std::size_t>(r.prevResidueIndex));
        if (prev.C != model::QtResidue::NONE) {
            cPrev = conf.atomPosition(frame, static_cast<std::size_t>(prev.C));
            cPrevValid = true;
        }
    }
    return BuildHNFrame(nPos, hPos, caPos, cPrev, cPrevValid);
}

}  // namespace

FeatureSchema McConnellNeighborhood::schema() const {
    FeatureSchema s;
    s.caseName = name();
    s.relationshipKind = RelationshipKind::SourceSum;
    s.sourceSchemaKind = SourceSchemaKind::Bond;
    s.includeMcLitKernel = true;

    s.sourceColumns = IdentityColumns();
    const std::vector<FeatureColumn> srcGeom = {
        {QStringLiteral("source_kind"), {}},
        {QStringLiteral("disp_local_x"), QStringLiteral("Angstrom")},
        {QStringLiteral("disp_local_y"), QStringLiteral("Angstrom")},
        {QStringLiteral("disp_local_z"), QStringLiteral("Angstrom")},
        {QStringLiteral("r"), QStringLiteral("Angstrom")},
        {QStringLiteral("cos_theta_bond_axis"), {}},
        {QStringLiteral("dipolar_3cos2m1_over_r3"), QStringLiteral("Angstrom^-3")},
        {QStringLiteral("bond_category"), {}},
        {QStringLiteral("bond_order"), {}},
        {QStringLiteral("bond_elem_a"), {}},
        {QStringLiteral("bond_elem_b"), {}},
        {QStringLiteral("bond_index"), {}},
        {QStringLiteral("bond_atom_a"), {}},
        {QStringLiteral("bond_atom_b"), {}},
        {QStringLiteral("bond_axis_local_x"), {}},
        {QStringLiteral("bond_axis_local_y"), {}},
        {QStringLiteral("bond_axis_local_z"), {}},
        {QStringLiteral("mc_lit_kernel_present"), {}},
        {QStringLiteral("mc_lit_T0"), QStringLiteral("ppm")},
        {QStringLiteral("mc_lit_T2_local_0"), QStringLiteral("ppm")},
        {QStringLiteral("mc_lit_T2_local_1"), QStringLiteral("ppm")},
        {QStringLiteral("mc_lit_T2_local_2"), QStringLiteral("ppm")},
        {QStringLiteral("mc_lit_T2_local_3"), QStringLiteral("ppm")},
        {QStringLiteral("mc_lit_T2_local_4"), QStringLiteral("ppm")},
    };
    for (const auto& c : srcGeom) s.sourceColumns.push_back(c);
    for (const auto& c : BareKernelColumns()) s.sourceColumns.push_back(c);
    for (const auto& c : TargetColumns()) s.sourceColumns.push_back(c);

    // McConnell sources are discovered bonds (no self-ring concept), so
    // n_sources_valid == n_sources and sum_dipolar_producer_valid == _all;
    // cutoff_A records the (required, non-default) source-discovery cutoff.
    s.aggregatedColumns = IdentityColumns();
    s.aggregatedColumns.push_back({QStringLiteral("n_sources"), {}});
    s.aggregatedColumns.push_back({QStringLiteral("n_sources_valid"), {}});
    s.aggregatedColumns.push_back({QStringLiteral("sum_dipolar_all"), QStringLiteral("Angstrom^-3")});
    s.aggregatedColumns.push_back({QStringLiteral("sum_dipolar_producer_valid"), QStringLiteral("Angstrom^-3")});
    for (const auto& cc : kCatCols)
        s.aggregatedColumns.push_back(
            {QStringLiteral("sum_dipolar_%1").arg(QString::fromLatin1(cc.tag)),
             QStringLiteral("Angstrom^-3")});
    s.aggregatedColumns.push_back({QStringLiteral("cutoff_A"), QStringLiteral("Angstrom")});
    s.aggregatedColumns.push_back({QStringLiteral("mc_lit_kernel_present"), {}});
    s.aggregatedColumns.push_back({QStringLiteral("mc_lit_T0"), QStringLiteral("ppm")});
    for (int i = 0; i < 5; ++i)
        s.aggregatedColumns.push_back(
            {QStringLiteral("mc_lit_T2_local_%1").arg(i), QStringLiteral("ppm")});
    for (const auto& c : BareKernelColumns()) s.aggregatedColumns.push_back(c);
    for (const auto& c : TargetColumns()) s.aggregatedColumns.push_back(c);

    return s;
}

std::size_t McConnellNeighborhood::extract(const Body& body, RecordSink& sink) const {
    const RunData& run = body.run;
    const model::QtProtein& p = *run.protein;
    const model::Conformation& conf = *run.conformation;
    const io::QtTrajectoryH5* h5 = run.h5();
    const model::QtShieldingTimeSeries* mc = h5 ? h5->mcShielding() : nullptr;

    std::vector<std::size_t> stratum;
    for (std::size_t i = 0; i < p.atomCount(); ++i)
        if (p.atom(i).IsBackboneAmideHydrogen()) stratum.push_back(i);
    qCInfo(cMc).noquote() << "mcconnell stratum HN atoms=" << stratum.size()
                          << "| cutoff_A=" << cutoff_A << "| dft rows=" << run.frameMap.dftRows().size();

    std::size_t cases = 0;
    for (std::size_t row : run.frameMap.dftRows()) {
        const std::size_t orig = run.frameMap.originalIndex(row);

        for (std::size_t atomIdx : stratum) {
            const model::QtAtom& a = p.atom(atomIdx);
            if (a.residueIndex < 0) continue;
            const LocalFrame frame =
                buildHN(p, conf, static_cast<std::size_t>(a.residueIndex), atomIdx, row);

            NeighborhoodRecord rec;
            FillIdentity(rec, run, atomIdx, row, name(), frame);
            if (mc) {
                rec.bare_kernel = mc->at(atomIdx, row);
                rec.bare_kernel_present = true;
            }
            rec.target = BuildTarget(run, atomIdx, orig, frame);

            const Vec3 atomPos = conf.atomPosition(row, atomIdx);

            std::vector<double> perCat(kCatCount, 0.0);
            double sumDipolar = 0.0;
            int nFinite = 0;

            for (const SourceRef& ref : body.idx.spatial.near(CloudKind::BondMidpoints, row, atomPos, cutoff_A)) {
                if (ref.entity_index < 0) continue;
                const model::QtBond& b = p.topology().bondAt(static_cast<std::size_t>(ref.entity_index));
                if (b.atomIndexA < 0 || b.atomIndexB < 0) continue;
                const Vec3 posA = conf.atomPosition(row, static_cast<std::size_t>(b.atomIndexA));
                const Vec3 posB = conf.atomPosition(row, static_cast<std::size_t>(b.atomIndexB));
                const Vec3 midpoint = 0.5 * (posA + posB);
                const Vec3 axis = (posB - posA);
                const double axisNorm = axis.norm();
                if (!(axisNorm > 1e-9)) continue;
                const Vec3 axisU = axis / axisNorm;

                // r and θ are about the BOND MIDPOINT and the BOND AXIS.
                const Vec3 disp = midpoint - atomPos;  // target ← midpoint (lab)
                const double r = disp.norm();
                if (!(r > 1e-6)) continue;  // exclude a bond whose own H is the atom
                const double cosT = disp.dot(axisU) / r;
                const double dipolar = (3.0 * cosT * cosT - 1.0) / (r * r * r);

                SourceSlot s;
                s.kind = SourceKind::Bond;
                s.disp_local = frame.is_valid ? frame.ToLocal(disp) : disp;
                s.r = r;
                s.cos_theta = cosT;
                s.dipolar = dipolar;
                s.bond_category = static_cast<int>(b.category);
                s.bond_order = static_cast<int>(b.order);
                s.bond_elem_a = static_cast<int>(p.atom(static_cast<std::size_t>(b.atomIndexA)).element);
                s.bond_elem_b = static_cast<int>(p.atom(static_cast<std::size_t>(b.atomIndexB)).element);
                s.bond_index = b.bondIndex;
                s.bond_atom_a = b.atomIndexA;
                s.bond_atom_b = b.atomIndexB;
                // The unit bond axis in the local frame — what an equivariant fit
                // needs to reconstruct the McConnell tensor (the scalar cosθ alone
                // throws the axis orientation away).
                s.bond_axis_local = frame.is_valid ? frame.ToLocal(axisU) : axisU;
                if (frame.is_valid) {
                    s.bond_mc_lit_kernel =
                        McConnellSourceLiteratureKernelLocal(s, &s.bond_mc_lit_kernel_present);
                }

                rec.sources.push_back(s);
                if (std::isfinite(dipolar)) {
                    sumDipolar += dipolar;
                    ++nFinite;
                    const int cc = catColumn(static_cast<int>(b.category));
                    if (cc >= 0) perCat[static_cast<std::size_t>(cc)] += dipolar;
                }
            }

            sink.WriteSourceRows(rec);
            // No self/bonded concept for bonds: valid == all. cutoff_A recorded.
            sink.WriteAggregatedRow(rec, sumDipolar, sumDipolar, nFinite, perCat, cutoff_A);
            ++cases;
        }
    }
    qCInfo(cMc).noquote() << "mcconnell cases=" << cases;
    return cases;
}

}  // namespace h5reader::rediscover
