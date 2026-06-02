#include "ComposedRelationships.h"

#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "McConnellLiteratureKernel.h"
#include "McConnellNeighborhood.h"
#include "RingCurrentKernel.h"
#include "RingCurrentNeighborhood.h"
#include "Verbs.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtProtein.h"
#include "../model/QtRing.h"
#include "../model/QtTopology.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace h5reader::rediscover {

namespace {

constexpr int kAromaticTypeCount = model::kAromaticRingTypeCount;  // 8

// McConnell's four anisotropic categories, in the schema's fixed column order
// (must match McConnellNeighborhood's kCatCols ordering exactly).
const model::BondCategory kMcCats[] = {
    model::BondCategory::PeptideCO,
    model::BondCategory::PeptideCN,
    model::BondCategory::SidechainCO,
    model::BondCategory::Aromatic,
};
constexpr int kMcCatCount = 4;

int mcCatColumn(int categoryOrd) {
    for (int i = 0; i < kMcCatCount; ++i)
        if (static_cast<int>(kMcCats[i]) == categoryOrd) return i;
    return -1;
}

// ── ring_current closures ──────────────────────────────────────────────────

// LocalFrameFn: the aromatic-H ring-normal frame, anchored on the typed CG
// (PHE/TYR/HIS/TRP-pyrrole) or CD2 (TRP-benzene) — verbs::atomOf(ring,
// {C, Γ|Δ}); FAIL LOUD on a non-unique anchor (the IUPAC-trap guard). Mirrors
// RingCurrentNeighborhood::extract's per-atom frame block exactly.
FrameResult ringFrameFn(const Body& body, std::size_t atom, std::size_t frame) {
    FrameResult fr;
    const model::QtProtein& p = *body.run.protein;
    const model::QtTopology& topo = p.topology();
    const int ownRing = verbs::ownAromaticRing(body, atom);
    if (ownRing < 0) return fr;  // frame invalid; invariants (z/r/cosθ) unaffected

    const model::RingGeometry& g = verbs::ringGeom(body, static_cast<std::size_t>(ownRing), frame);
    const model::QtRing& ring = topo.ringAt(static_cast<std::size_t>(ownRing));
    TypedAtomSelector sel;
    sel.element = model::Element::C;
    sel.locant = (ring.TypeIndex() == model::RingTypeIndex::TrpBenzene) ? model::Locant::Delta
                                                                        : model::Locant::Gamma;
    QString anchorErr;
    const std::optional<int32_t> anchor = verbs::atomOf(body, ring.atomIndices, sel, &anchorErr);
    if (!anchor) {
        throw std::runtime_error(
            QStringLiteral("ring_current typed frame anchor failed for atom %1 ring %2: %3")
                .arg(atom)
                .arg(ownRing)
                .arg(anchorErr)
                .toStdString());
    }
    fr.anchor_atom_index = *anchor;
    const Vec3 anchorPos = verbs::pos(body, static_cast<std::size_t>(*anchor), frame);
    fr.frame = BuildAromaticHFrame(g.center, g.normal, anchorPos);
    return fr;
}

// Attacher (ring geometry + identity, from the H5 slot + RingGeometryCache +
// the QtRing virtuals). Mirrors the per-slot body of RingCurrentNeighborhood.
void ringAttacher(const Body& body, const AtomState& st, const FrameResult& fr,
                  const RawSource& raw, SourceSlot& s) {
    if (raw.kind != SourceKind::Ring) return;
    const verbs::RingSlot& rs = raw.ring;
    const model::QtTopology& topo = body.run.protein->topology();
    const model::QtRing& sring = topo.ringAt(static_cast<std::size_t>(rs.ring_index));

    const double cosT = rs.z / rs.distance;
    const double dipolar = (3.0 * cosT * cosT - 1.0) / (rs.distance * rs.distance * rs.distance);

    s.kind = SourceKind::Ring;
    s.r = rs.distance;
    s.cos_theta = cosT;
    s.dipolar = dipolar;
    s.ring_z = rs.z;
    s.ring_rho = rs.rho;
    s.ring_in_plane_angle = rs.in_plane_angle;
    s.ring_index = rs.ring_index;

    // Lab displacement target→source-ring-center, expressed in the local frame.
    const model::RingGeometry& sg =
        verbs::ringGeom(body, static_cast<std::size_t>(rs.ring_index), st.frame);
    const Vec3 dispLab = sg.center - st.pos;
    s.disp_local = fr.frame.is_valid ? fr.frame.ToLocal(dispLab) : dispLab;
    const Vec3 nrm = sg.normal.norm() > 1e-9 ? sg.normal.normalized() : sg.normal;
    s.source_normal_local = fr.frame.is_valid ? fr.frame.ToLocal(nrm) : nrm;

    // Ring physics from the typed QtRing virtuals (objects answer for themselves).
    s.ring_type_index = sring.TypeIndexAsInt();
    s.ring_intensity = sring.LiteratureIntensity();
    s.ring_nitrogen = sring.NitrogenCount();
    s.ring_jb_offset = sring.JohnsonBoveyLobeOffset();
    s.ring_aromaticity = static_cast<int>(sring.Aromaticity());
    s.ring_size = sring.RingSizeValue();
    s.ring_fused = sring.IsFused();
    s.ring_jb_unit_kernel =
        JohnsonBoveySourceUnitKernelLocal(body, fr.frame, st.pos,
                                          static_cast<std::size_t>(rs.ring_index), sring, st.frame);
    s.ring_jb_kernel = ScaleSphericalTensor(s.ring_jb_unit_kernel,
                                            sring.LiteratureIntensity());
    s.ring_jb_kernel_present = true;

    // is_self_or_bonded is set by the per-source classifier (it needs the atom's
    // own-ring / own-atom sets, prepped once per atom); the attacher leaves it
    // default false. The reducer then reads the flag off the slot.
}

// ClassifierPrep: build the atom's own-ring id set + own-ring atom set (the
// producer's self/bonded exclusion) once per atom. verbs::ringsOf /
// verbs::ownRingAtoms read the resident topology reverse index — no
// re-derivation of connectivity.
ClassifierContext ringClassifierPrep(const Body& body, std::size_t atom) {
    ClassifierContext ctx;
    ctx.ownRings = verbs::ringsOf(body, atom);
    ctx.ownAtoms = verbs::ownRingAtoms(body, atom);
    return ctx;
}

// Classifier (SourceClassifier): stamp is_self_or_bonded onto each ring source
// — own ring, or any atom shared with an own ring (fused via a bridgehead).
// Runs before the per-source row is written, so the flag lands on the row.
void ringClassifier(const Body& body, const ClassifierContext& ctx, SourceSlot& s) {
    if (s.kind != SourceKind::Ring || s.ring_index < 0) return;
    bool selfOrBonded = false;
    for (int32_t r : ctx.ownRings)
        if (r == s.ring_index) { selfOrBonded = true; break; }
    if (!selfOrBonded) {
        const model::QtRing& sring =
            body.run.protein->topology().ringAt(static_cast<std::size_t>(s.ring_index));
        for (int32_t ra : sring.atomIndices) {
            bool hit = false;
            for (int32_t oa : ctx.ownAtoms)
                if (oa == ra) { hit = true; break; }
            if (hit) { selfOrBonded = true; break; }
        }
    }
    s.is_self_or_bonded = selfOrBonded;
}

// Reducer: the dual sum (all vs producer-valid) + per-ring-type sums on the
// valid set, reading the classifier's is_self_or_bonded off each slot. This is
// the STATE-CARRYING fold; the classification is already on the slots.
AggregateResult ringReducer(const Body& /*body*/, std::size_t /*atom*/,
                            const std::vector<SourceSlot>& sources) {
    AggregateResult agg;
    agg.per_type.assign(static_cast<std::size_t>(kAromaticTypeCount), 0.0);
    agg.cutoff_A = std::numeric_limits<double>::quiet_NaN();  // sources from H5

    for (const SourceSlot& s : sources) {
        if (!std::isfinite(s.dipolar)) continue;
        agg.sum_all += s.dipolar;
        if (!s.is_self_or_bonded) {
            agg.sum_valid += s.dipolar;
            ++agg.n_valid;
            if (s.ring_type_index >= 0 && s.ring_type_index < kAromaticTypeCount)
                agg.per_type[static_cast<std::size_t>(s.ring_type_index)] += s.dipolar;
        }
    }
    return agg;
}

// BareKernelFn: the producer's bs kernel at (atom, frame).
void ringBareKernel(const Body& body, std::size_t atom, std::size_t frame,
                    NeighborhoodRecord& rec) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (h5 && h5->bsShielding()) {
        rec.bare_kernel = h5->bsShielding()->at(atom, frame);
        rec.bare_kernel_present = true;
    }
}

// ── mcconnell closures ─────────────────────────────────────────────────────

// LocalFrameFn: the HN amide-plane frame from the residue's N / H / CA / C_prev
// (N-terminus fallback to CA). Mirrors McConnellNeighborhood::buildHN exactly.
FrameResult mcFrameFn(const Body& body, std::size_t atom, std::size_t frame) {
    FrameResult fr;
    const model::QtProtein& p = *body.run.protein;
    const model::QtAtom& a = p.atom(atom);
    if (a.residueIndex < 0) return fr;
    const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
    if (r.N == model::QtResidue::NONE || r.CA == model::QtResidue::NONE) return fr;

    const Vec3 nPos = verbs::pos(body, static_cast<std::size_t>(r.N), frame);
    const Vec3 hPos = verbs::pos(body, atom, frame);
    const Vec3 caPos = verbs::pos(body, static_cast<std::size_t>(r.CA), frame);

    bool cPrevValid = false;
    Vec3 cPrev = Vec3::Zero();
    if (r.prevResidueIndex >= 0 && static_cast<std::size_t>(r.prevResidueIndex) < p.residueCount()) {
        const model::QtResidue& prev = p.residue(static_cast<std::size_t>(r.prevResidueIndex));
        if (prev.C != model::QtResidue::NONE) {
            cPrev = verbs::pos(body, static_cast<std::size_t>(prev.C), frame);
            cPrevValid = true;
        }
    }
    fr.frame = BuildHNFrame(nPos, hPos, caPos, cPrev, cPrevValid);
    fr.anchor_atom_index = -1;  // HN has no typed azimuth anchor
    return fr;
}

// Attacher (bond geometry + identity + axis). Mirrors the per-bond body of
// McConnellNeighborhood. A degenerate axis or r<=1e-6 marks the slot invalid
// (dipolar = NaN) so the reducer drops it — matching the cell's `continue`.
void mcAttacher(const Body& body, const AtomState& st, const FrameResult& fr,
                const RawSource& raw, SourceSlot& s) {
    if (raw.kind != SourceKind::Bond || raw.ref.entity_index < 0) {
        s.dipolar = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const model::QtProtein& p = *body.run.protein;
    const model::QtBond& b = p.topology().bondAt(static_cast<std::size_t>(raw.ref.entity_index));
    if (b.atomIndexA < 0 || b.atomIndexB < 0) {
        s.dipolar = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const Vec3 posA = verbs::pos(body, static_cast<std::size_t>(b.atomIndexA), st.frame);
    const Vec3 posB = verbs::pos(body, static_cast<std::size_t>(b.atomIndexB), st.frame);
    const Vec3 midpoint = 0.5 * (posA + posB);
    const Vec3 axis = posB - posA;
    const double axisNorm = axis.norm();
    if (!(axisNorm > 1e-9)) {
        s.dipolar = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const Vec3 axisU = axis / axisNorm;

    const Vec3 disp = midpoint - st.pos;  // target ← midpoint (lab)
    const double r = disp.norm();
    if (!(r > 1e-6)) {  // a bond whose own H is the atom
        s.dipolar = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const double cosT = disp.dot(axisU) / r;
    const double dipolar = (3.0 * cosT * cosT - 1.0) / (r * r * r);

    s.kind = SourceKind::Bond;
    s.disp_local = fr.frame.is_valid ? fr.frame.ToLocal(disp) : disp;
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
    s.bond_axis_local = fr.frame.is_valid ? fr.frame.ToLocal(axisU) : axisU;
    if (fr.frame.is_valid) {
        s.bond_mc_lit_kernel =
            McConnellSourceLiteratureKernelLocal(s, &s.bond_mc_lit_kernel_present);
    }
}

// Reducer: the plain per-category sum, no self/bonded concept (valid == all).
AggregateResult mcReducer(const Body& /*body*/, std::size_t /*atom*/,
                          const std::vector<SourceSlot>& sources, double cutoff_A) {
    AggregateResult agg;
    agg.per_type.assign(static_cast<std::size_t>(kMcCatCount), 0.0);
    agg.cutoff_A = cutoff_A;
    for (const SourceSlot& s : sources) {
        if (!std::isfinite(s.dipolar)) continue;
        agg.sum_all += s.dipolar;
        ++agg.n_valid;
        const int cc = mcCatColumn(s.bond_category);
        if (cc >= 0) agg.per_type[static_cast<std::size_t>(cc)] += s.dipolar;
    }
    agg.sum_valid = agg.sum_all;  // no self/bonded for bonds
    return agg;
}

void mcBareKernel(const Body& body, std::size_t atom, std::size_t frame, NeighborhoodRecord& rec) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    if (h5 && h5->mcShielding()) {
        rec.bare_kernel = h5->mcShielding()->at(atom, frame);
        rec.bare_kernel_present = true;
    }
}

}  // namespace

Relationship MakeRingCurrentRelationship() {
    Relationship rel;
    rel.name = QStringLiteral("ring_current");
    // Reuse the proven cell schema verbatim (byte-identical CSV headers) — the
    // schema is static data; the cell stays as the reference oracle.
    rel.schema = RingCurrentNeighborhood().schema();

    rel.stratum = atomsWhere([](const model::QtAtom& a) { return a.IsAromaticRingHydrogen(); });
    rel.frame_fn = &ringFrameFn;
    rel.selectors = {slotsBackend()};   // the frozen-membership ring backend
    rel.attachers = {&ringAttacher};
    rel.classifier_prep = &ringClassifierPrep;
    rel.classifier = &ringClassifier;
    rel.reducer = &ringReducer;
    rel.target_fn = [](const Body& body, std::size_t atom, std::size_t orig,
                       const LocalFrame& f) { return BuildTarget(body.run, atom, orig, f); };
    rel.bare_kernel_fn = &ringBareKernel;
    return rel;
}

Relationship MakeMcConnellRelationship(double cutoff_A) {
    Relationship rel;
    rel.name = QStringLiteral("mcconnell");
    {
        McConnellNeighborhood cell;
        cell.cutoff_A = cutoff_A;
        rel.schema = cell.schema();
    }

    // Stratum: backbone amide HN with a resolved residue. The procedural cell
    // `continue`s (emits no case) for an HN whose residueIndex < 0; residueIndex
    // is static, so folding it into the stratum reproduces that exactly (every
    // backbone HN has a residue in practice, but faithful == faithful).
    rel.stratum = atomsWhere([](const model::QtAtom& a) {
        return a.IsBackboneAmideHydrogen() && a.residueIndex >= 0;
    });
    rel.frame_fn = &mcFrameFn;
    rel.selectors = {nearBackend(CloudKind::BondMidpoints, cutoff_A)};
    rel.attachers = {&mcAttacher};
    // Drop the geometrically-rejected bonds (degenerate axis / coincident
    // midpoint / own-H): mcAttacher marked them with a NaN dipolar. This is the
    // cell's `continue` — those sources are never emitted as rows.
    rel.source_filter = [](const SourceSlot& s) { return std::isfinite(s.dipolar); };
    rel.reducer = [cutoff_A](const Body& body, std::size_t atom,
                             const std::vector<SourceSlot>& srcs) {
        return mcReducer(body, atom, srcs, cutoff_A);
    };
    rel.target_fn = [](const Body& body, std::size_t atom, std::size_t orig,
                       const LocalFrame& f) { return BuildTarget(body.run, atom, orig, f); };
    rel.bare_kernel_fn = &mcBareKernel;
    return rel;
}

}  // namespace h5reader::rediscover
