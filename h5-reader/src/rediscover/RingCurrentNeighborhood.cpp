#include "RingCurrentNeighborhood.h"

#include "AnalysisBody.h"
#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "RingCurrentKernel.h"
#include "TypedAtomIndex.h"

#include "../model/ConformationGeometry.h"
#include "../model/QtProtein.h"
#include "../io/QtTrajectoryH5.h"

#include <QLoggingCategory>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cRing, "h5reader.rediscover.ring")

// Per-aromatic-ring-type summed dipolar feature columns (the aggregated row's
// per-ring-type sums). 8 aromatic types (RingTypeIndex 0..7); the saturated
// PRO ring (8) carries no ring current and is excluded.
constexpr int kAromaticTypeCount = model::kAromaticRingTypeCount;  // 8

// Find an aromatic ring the atom (its heavy parent) belongs to. Returns the
// absolute rings.npy row, or -1. Used for the atom's own ring-normal frame.
int ownAromaticRing(const model::QtProtein& p, std::size_t atomIdx) {
    const model::QtAtom& a = p.atom(atomIdx);
    const std::size_t heavy =
        a.parentAtomIndex >= 0 ? static_cast<std::size_t>(a.parentAtomIndex) : atomIdx;
    const model::QtTopology& topo = p.topology();
    for (int memb : topo.ringMembershipsForAtom(heavy)) {
        const model::QtRingMembership& m = topo.ringMembershipAt(static_cast<std::size_t>(memb));
        if (m.ringId < 0) continue;
        const model::QtRing& ring = topo.ringAt(static_cast<std::size_t>(m.ringId));
        if (ring.IsAromatic()) return m.ringId;
    }
    return -1;
}

}  // namespace

FeatureSchema RingCurrentNeighborhood::schema() const {
    FeatureSchema s;
    s.caseName = name();
    s.relationshipKind = RelationshipKind::SourceSum;
    s.sourceSchemaKind = SourceSchemaKind::Ring;

    // ── per-source columns ──
    s.sourceColumns = IdentityColumns();
    const std::vector<FeatureColumn> srcGeom = {
        {QStringLiteral("source_kind"), {}},
        {QStringLiteral("disp_local_x"), QStringLiteral("Angstrom")},
        {QStringLiteral("disp_local_y"), QStringLiteral("Angstrom")},
        {QStringLiteral("disp_local_z"), QStringLiteral("Angstrom")},
        {QStringLiteral("r"), QStringLiteral("Angstrom")},
        {QStringLiteral("cos_theta"), {}},
        {QStringLiteral("dipolar_3cos2m1_over_r3"), QStringLiteral("Angstrom^-3")},
        {QStringLiteral("ring_z"), QStringLiteral("Angstrom")},
        {QStringLiteral("ring_rho"), QStringLiteral("Angstrom")},
        {QStringLiteral("ring_in_plane_angle"), QStringLiteral("radians")},
        {QStringLiteral("ring_index"), {}},
        {QStringLiteral("is_self_or_bonded"), {}},
        {QStringLiteral("ring_type_index"), {}},
        {QStringLiteral("ring_intensity"), QStringLiteral("nA/T")},
        {QStringLiteral("ring_nitrogen_count"), {}},
        {QStringLiteral("ring_jb_offset"), QStringLiteral("Angstrom")},
        {QStringLiteral("ring_aromaticity_ord"), {}},
        {QStringLiteral("ring_size"), {}},
        {QStringLiteral("ring_fused"), {}},
        {QStringLiteral("source_normal_local_x"), {}},
        {QStringLiteral("source_normal_local_y"), {}},
        {QStringLiteral("source_normal_local_z"), {}},
        {QStringLiteral("jb_kernel_present"), {}},
        {QStringLiteral("jb_unit_T0"), QStringLiteral("ppm_T_per_nA")},
        {QStringLiteral("jb_unit_T2_local_0"), QStringLiteral("ppm_T_per_nA")},
        {QStringLiteral("jb_unit_T2_local_1"), QStringLiteral("ppm_T_per_nA")},
        {QStringLiteral("jb_unit_T2_local_2"), QStringLiteral("ppm_T_per_nA")},
        {QStringLiteral("jb_unit_T2_local_3"), QStringLiteral("ppm_T_per_nA")},
        {QStringLiteral("jb_unit_T2_local_4"), QStringLiteral("ppm_T_per_nA")},
        {QStringLiteral("jb_T0"), QStringLiteral("ppm")},
        {QStringLiteral("jb_T2_local_0"), QStringLiteral("ppm")},
        {QStringLiteral("jb_T2_local_1"), QStringLiteral("ppm")},
        {QStringLiteral("jb_T2_local_2"), QStringLiteral("ppm")},
        {QStringLiteral("jb_T2_local_3"), QStringLiteral("ppm")},
        {QStringLiteral("jb_T2_local_4"), QStringLiteral("ppm")},
    };
    for (const auto& c : srcGeom) s.sourceColumns.push_back(c);
    for (const auto& c : BareKernelColumns()) s.sourceColumns.push_back(c);
    for (const auto& c : TargetColumns()) s.sourceColumns.push_back(c);

    // ── aggregated columns ──
    // sum_dipolar_all = Σ over all neighbourhood rings (incl. self/bonded, which
    // the producer excludes); sum_dipolar_producer_valid = Σ over the
    // through-space (producer-valid) set; per-ring-type sums on the valid set.
    s.aggregatedColumns = IdentityColumns();
    s.aggregatedColumns.push_back({QStringLiteral("n_sources"), {}});
    s.aggregatedColumns.push_back({QStringLiteral("n_sources_valid"), {}});
    s.aggregatedColumns.push_back({QStringLiteral("sum_dipolar_all"), QStringLiteral("Angstrom^-3")});
    s.aggregatedColumns.push_back({QStringLiteral("sum_dipolar_producer_valid"), QStringLiteral("Angstrom^-3")});
    for (int t = 0; t < kAromaticTypeCount; ++t)
        s.aggregatedColumns.push_back(
            {QStringLiteral("sum_dipolar_ringtype_%1").arg(t), QStringLiteral("Angstrom^-3")});
    s.aggregatedColumns.push_back({QStringLiteral("cutoff_A"), QStringLiteral("Angstrom")});
    for (const auto& c : BareKernelColumns()) s.aggregatedColumns.push_back(c);
    for (const auto& c : TargetColumns()) s.aggregatedColumns.push_back(c);

    return s;
}

std::size_t RingCurrentNeighborhood::extract(const Body& body, RecordSink& sink) const {
    const RunData& run = body.run;
    const model::QtProtein& p = *run.protein;
    const model::Conformation& conf = *run.conformation;
    const io::QtTrajectoryH5* h5 = run.h5();
    const model::QtRingNeighbourhoodTimeSeries* rn = h5 ? h5->ringNeighbourhood() : nullptr;
    const model::QtShieldingTimeSeries* bs = h5 ? h5->bsShielding() : nullptr;

    if (!rn) {
        qCWarning(cRing) << "no ring-neighbourhood TS in H5 — ring-current extraction empty";
        return 0;
    }

    // Stratum atoms (aromatic ring-facing H), computed once.
    std::vector<std::size_t> stratum;
    for (std::size_t i = 0; i < p.atomCount(); ++i)
        if (p.atom(i).IsAromaticRingHydrogen()) stratum.push_back(i);
    qCInfo(cRing).noquote() << "ring-current stratum atoms=" << stratum.size()
                            << "| dft rows=" << run.frameMap.dftRows().size();

    const model::QtTopology& topo = p.topology();

    // Per stratum atom, precompute the producer's self/bonded exclusion set: the
    // aromatic rings the H's parent heavy atom belongs to (self), and the union
    // of those rings' atoms (catches fused partners, which share a bridgehead).
    // A source ring is self/bonded if it IS an own ring or shares any atom.
    std::unordered_map<std::size_t, std::unordered_set<int32_t>> ownRingsByAtom;
    std::unordered_map<std::size_t, std::unordered_set<int32_t>> ownAtomsByAtom;
    for (std::size_t atomIdx : stratum) {
        const model::QtAtom& a = p.atom(atomIdx);
        const std::size_t heavy =
            a.parentAtomIndex >= 0 ? static_cast<std::size_t>(a.parentAtomIndex) : atomIdx;
        std::unordered_set<int32_t> rings, atoms;
        for (int memb : topo.ringMembershipsForAtom(heavy)) {
            const model::QtRingMembership& m = topo.ringMembershipAt(static_cast<std::size_t>(memb));
            if (m.ringId < 0) continue;
            const model::QtRing& ring = topo.ringAt(static_cast<std::size_t>(m.ringId));
            if (!ring.IsAromatic()) continue;
            rings.insert(m.ringId);
            for (int32_t ra : ring.atomIndices) atoms.insert(ra);
        }
        ownRingsByAtom[atomIdx] = std::move(rings);
        ownAtomsByAtom[atomIdx] = std::move(atoms);
    }

    constexpr double kNoCutoff = std::numeric_limits<double>::quiet_NaN();  // sources from H5
    std::size_t cases = 0;
    for (std::size_t row : run.frameMap.dftRows()) {
        const std::size_t orig = run.frameMap.originalIndex(row);
        for (std::size_t atomIdx : stratum) {
            // The atom's own ring-normal frame, anchored on the chemistry-typed
            // CG/CD2 (substrate_conventions), with the normal flipped to the
            // canonical ring-traversal direction so the azimuthal frame and the
            // T2 components are stable across frames/tautomers.
            LocalFrame frame;
            int anchorIdx = -1;
            const int ownRing = ownAromaticRing(p, atomIdx);
            if (ownRing >= 0) {
                model::RingGeometry g = body.idx.ringGeometry.at(static_cast<std::size_t>(ownRing), row);
                const model::QtRing& ring = topo.ringAt(static_cast<std::size_t>(ownRing));
                TypedAtomSelector sel;
                sel.element = model::Element::C;
                sel.locant = (ring.TypeIndex() == model::RingTypeIndex::TrpBenzene)
                                 ? model::Locant::Delta
                                 : model::Locant::Gamma;
                QString anchorErr;
                const std::optional<int32_t> anchor = body.idx.typedAtoms.selectUnique(
                    ring.atomIndices, sel, &anchorErr);
                if (!anchor) {
                    throw std::runtime_error(
                        QStringLiteral("ring_current typed frame anchor failed for atom %1 ring %2: %3")
                            .arg(atomIdx)
                            .arg(ownRing)
                            .arg(anchorErr)
                            .toStdString());
                }
                anchorIdx = *anchor;
                const Vec3 anchorPos = conf.atomPosition(row, static_cast<std::size_t>(anchorIdx));
                frame = BuildAromaticHFrame(g.center, g.normal, anchorPos);
            }

            NeighborhoodRecord rec;
            FillIdentity(rec, run, atomIdx, row, name(), frame);
            rec.frame_anchor_atom_index = anchorIdx;
            if (bs) {
                rec.bare_kernel = bs->at(atomIdx, row);
                rec.bare_kernel_present = true;
            }
            rec.target = BuildTarget(run, atomIdx, orig, frame);

            const Vec3 atomPos = conf.atomPosition(row, atomIdx);
            const std::unordered_set<int32_t>& ownRings = ownRingsByAtom[atomIdx];
            const std::unordered_set<int32_t>& ownAtoms = ownAtomsByAtom[atomIdx];

            // Two aggregated sums: ALL neighbourhood rings, and the producer-valid
            // (through-space) subset that drops self/bonded rings. Per-type sums
            // are on the valid subset.
            std::vector<double> perType(static_cast<std::size_t>(kAromaticTypeCount), 0.0);
            double sumAll = 0.0, sumValid = 0.0;
            int nValid = 0;

            // Walk the H5 ring-neighbourhood slots (frame-0 membership snapshot).
            for (std::size_t slot = 0; slot < rn->n_slots; ++slot) {
                const int32_t srcRingId = rn->ringIndexAt(atomIdx, slot);
                if (srcRingId < 0) continue;  // unused slot
                const std::array<double, 4> ch = rn->at(atomIdx, row, slot);
                const double distance = ch[0];
                const double rho = ch[1];
                const double z = ch[2];
                const double inPlane = ch[3];
                if (!std::isfinite(distance) || !(distance > 1e-6)) continue;

                const double cosT = z / distance;
                const double dipolar = (3.0 * cosT * cosT - 1.0) / (distance * distance * distance);

                const model::QtRing& sring = topo.ringAt(static_cast<std::size_t>(srcRingId));

                // Self/bonded: own ring, or any shared atom with an own ring (fused).
                bool selfOrBonded = ownRings.count(srcRingId) > 0;
                if (!selfOrBonded)
                    for (int32_t ra : sring.atomIndices)
                        if (ownAtoms.count(ra)) { selfOrBonded = true; break; }

                SourceSlot s;
                s.kind = SourceKind::Ring;
                s.r = distance;
                s.cos_theta = cosT;
                s.dipolar = dipolar;
                s.ring_z = z;
                s.ring_rho = rho;
                s.ring_in_plane_angle = inPlane;
                s.ring_index = srcRingId;
                s.is_self_or_bonded = selfOrBonded;

                // Lab displacement target→source-ring-center, in the atom's frame.
                const model::RingGeometry sg = body.idx.ringGeometry.at(static_cast<std::size_t>(srcRingId), row);
                const Vec3 dispLab = sg.center - atomPos;
                s.disp_local = frame.is_valid ? frame.ToLocal(dispLab) : dispLab;
                // The source ring's unit normal in the target's local frame — the
                // dipole axis for the l=2 ring-current tensor (sign-irrelevant).
                const Vec3 nrm = sg.normal.norm() > 1e-9 ? sg.normal.normalized() : sg.normal;
                s.source_normal_local = frame.is_valid ? frame.ToLocal(nrm) : nrm;

                // Ring physics from the typed QtRing virtual.
                s.ring_type_index = sring.TypeIndexAsInt();
                s.ring_intensity = sring.LiteratureIntensity();
                s.ring_nitrogen = sring.NitrogenCount();
                s.ring_jb_offset = sring.JohnsonBoveyLobeOffset();
                s.ring_aromaticity = static_cast<int>(sring.Aromaticity());
                s.ring_size = sring.RingSizeValue();
                s.ring_fused = sring.IsFused();
                s.ring_jb_unit_kernel =
                    JohnsonBoveySourceUnitKernelLocal(body, frame, atomPos,
                                                      static_cast<std::size_t>(srcRingId), sring, row);
                s.ring_jb_kernel = ScaleSphericalTensor(s.ring_jb_unit_kernel,
                                                        sring.LiteratureIntensity());
                s.ring_jb_kernel_present = true;

                rec.sources.push_back(s);
                if (std::isfinite(dipolar)) {
                    sumAll += dipolar;
                    if (!selfOrBonded) {
                        sumValid += dipolar;
                        ++nValid;
                        if (s.ring_type_index >= 0 && s.ring_type_index < kAromaticTypeCount)
                            perType[static_cast<std::size_t>(s.ring_type_index)] += dipolar;
                    }
                }
            }

            sink.WriteSourceRows(rec);
            sink.WriteAggregatedRow(rec, sumAll, sumValid, nValid, perType, kNoCutoff);
            ++cases;
        }
    }
    qCInfo(cRing).noquote() << "ring-current cases=" << cases;
    return cases;
}

}  // namespace h5reader::rediscover
