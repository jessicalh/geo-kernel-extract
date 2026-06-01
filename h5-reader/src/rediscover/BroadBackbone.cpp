#include "BroadBackbone.h"

#include "Catalog.h"
#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "RingGeometryCache.h"
#include "SpatialIndexSet.h"
#include "Verbs.h"

#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtRing.h"
#include "../model/QtTopology.h"

#include <QLoggingCategory>

#include <cmath>
#include <limits>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cBroad, "h5reader.rediscover.broad_backbone")

const double kNaN = std::numeric_limits<double>::quiet_NaN();

bool validResidue(const model::QtProtein& p, int32_t r) {
    return r >= 0 && static_cast<std::size_t>(r) < p.residueCount();
}

// ── Backbone-class-dispatching LocalFrameFn (the new conceptual piece) ───────
// Dispatches on the atom's typed BackboneRole (plus the GLY-HA Locant case) and
// builds the matching backbone frame from TYPED anchor atoms — the QtResidue
// N/CA/C/O/HA cache, built at load from BackboneRole+Locant (collision-safe,
// never positional; the IUPAC-revert trap). Reuses the resident topology; no
// re-parse, no name scan, no positional index. The previous-residue carbonyl C
// (for the N interior frame) comes from the typed prevResidueIndex link.
FrameResult backboneFrameFn(const Body& body, std::size_t atom, std::size_t frame) {
    FrameResult fr;
    const model::QtProtein& p = *body.run.protein;
    const model::QtAtom& a = p.atom(atom);
    if (!validResidue(p, a.residueIndex)) return fr;  // invalid; invariants unaffected
    const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));

    auto posOf = [&](int32_t ai) {
        return verbs::pos(body, static_cast<std::size_t>(ai), frame);
    };

    // N frame: z = N→CA, in-plane ref = N→C_prev (interior) or N→C_own (N-term).
    if (a.IsBackboneNitrogen()) {
        if (r.N == model::QtResidue::NONE || r.CA == model::QtResidue::NONE) return fr;
        const Vec3 nPos = posOf(r.N);
        const Vec3 caPos = posOf(r.CA);
        bool cPrevValid = false;
        Vec3 cRef = Vec3::Zero();
        if (validResidue(p, r.prevResidueIndex)) {
            const model::QtResidue& prev = p.residue(static_cast<std::size_t>(r.prevResidueIndex));
            if (prev.C != model::QtResidue::NONE) {
                cRef = posOf(prev.C);
                cPrevValid = true;
            }
        }
        if (!cPrevValid && r.C != model::QtResidue::NONE) cRef = posOf(r.C);  // N-term fallback
        fr.frame = BuildBackboneNFrame(nPos, caPos, cRef, cPrevValid);
        fr.anchor_atom_index = cPrevValid && validResidue(p, r.prevResidueIndex)
                                   ? p.residue(static_cast<std::size_t>(r.prevResidueIndex)).C
                                   : r.C;
        return fr;
    }

    // Cα frame: z = bisector(Cα→N, Cα→C), x = Cα→N.
    if (a.IsBackboneAlphaCarbon()) {
        if (r.CA == model::QtResidue::NONE || r.N == model::QtResidue::NONE
            || r.C == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCaFrame(posOf(r.CA), posOf(r.N), posOf(r.C));
        fr.anchor_atom_index = r.N;
        return fr;
    }

    // Carbonyl C frame: z = C→O, x in-plane ref C→CA.
    if (a.IsBackboneCarbonylCarbon()) {
        if (r.C == model::QtResidue::NONE || r.O == model::QtResidue::NONE
            || r.CA == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCarbonylCFrame(posOf(r.C), posOf(r.O), posOf(r.CA));
        fr.anchor_atom_index = r.CA;
        return fr;
    }

    // Carbonyl O frame: z = O→C, x in-plane ref C→CA.
    if (a.IsBackboneCarbonylOxygen()) {
        if (r.O == model::QtResidue::NONE || r.C == model::QtResidue::NONE
            || r.CA == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCarbonylOFrame(posOf(r.O), posOf(r.C), posOf(r.CA));
        fr.anchor_atom_index = r.CA;
        return fr;
    }

    // Backbone amide H (HN): reuse the existing HN amide-plane frame (untouched
    // builder — z = N→H, x in amide plane via C_prev). This is the SAME frame the
    // mcconnell oracle uses; broad-backbone shares it for the HN class.
    if (a.IsBackboneAmideHydrogen()) {
        if (r.N == model::QtResidue::NONE || r.CA == model::QtResidue::NONE) return fr;
        const Vec3 nPos = posOf(r.N);
        const Vec3 hPos = verbs::pos(body, atom, frame);
        const Vec3 caPos = posOf(r.CA);
        bool cPrevValid = false;
        Vec3 cPrev = Vec3::Zero();
        if (validResidue(p, r.prevResidueIndex)) {
            const model::QtResidue& prev = p.residue(static_cast<std::size_t>(r.prevResidueIndex));
            if (prev.C != model::QtResidue::NONE) {
                cPrev = posOf(prev.C);
                cPrevValid = true;
            }
        }
        fr.frame = BuildHNFrame(nPos, hPos, caPos, cPrev, cPrevValid);
        fr.anchor_atom_index = -1;
        return fr;
    }

    // HA (Cα-H): non-GLY via BackboneRole::AlphaHydrogen; GLY HA2/HA3 via the
    // IsAnyAlphaHydrogen predicate (Locant::Alpha + BackboneRole::None per
    // Markley). z = Cα→HA, x = Cα→N.
    if (a.IsAnyAlphaHydrogen()) {
        if (r.CA == model::QtResidue::NONE || r.N == model::QtResidue::NONE) return fr;
        fr.frame = BuildBackboneHaFrame(verbs::pos(body, atom, frame), posOf(r.CA), posOf(r.N));
        fr.anchor_atom_index = r.N;
        return fr;
    }

    return fr;  // not a recognised backbone class (shouldn't happen given the stratum)
}

// ── Ring attacher (GENERAL KD ring-centres backend, NOT the slots backend) ───
// Source = a ring discovered by nearBackend(RingCenters): raw.ref.entity_index
// is the absolute ring index. Geometry from the resident RingGeometryCache; ring
// physics from the typed QtRing virtuals. is_self_or_bonded flags the atom's own
// ring (for a backbone atom in PRO/HIS/etc. whose residue carries a ring).
void ringAttacher(const Body& body, const AtomState& st, const FrameResult& fr,
                  const RawSource& raw, SourceSlot& s) {
    if (raw.kind != SourceKind::Ring) return;
    if (raw.ref.entity_index < 0) {
        s.dipolar = kNaN;
        return;
    }
    const std::size_t ringIdx = static_cast<std::size_t>(raw.ref.entity_index);
    const model::QtTopology& topo = body.run.protein->topology();
    if (ringIdx >= topo.ringCount()) {
        s.dipolar = kNaN;
        return;
    }
    const model::RingGeometry& g = verbs::ringGeom(body, ringIdx, st.frame);
    const Vec3 dispLab = g.center - st.pos;
    const double r = dispLab.norm();
    if (!(r > 1e-6)) {
        s.dipolar = kNaN;
        return;
    }
    const Vec3 nrm = g.normal.norm() > 1e-9 ? g.normal.normalized() : g.normal;
    const double cosT = dispLab.dot(nrm) / r;
    const double dipolar = (3.0 * cosT * cosT - 1.0) / (r * r * r);

    s.kind = SourceKind::Ring;
    s.r = r;
    s.cos_theta = cosT;
    s.dipolar = dipolar;
    s.disp_local = fr.frame.is_valid ? fr.frame.ToLocal(dispLab) : dispLab;
    s.source_normal_local = fr.frame.is_valid ? fr.frame.ToLocal(nrm) : nrm;
    s.ring_index = static_cast<int32_t>(ringIdx);

    const model::QtRing& sring = topo.ringAt(ringIdx);
    s.ring_type_index = sring.TypeIndexAsInt();
    s.ring_intensity = sring.LiteratureIntensity();
    s.ring_nitrogen = sring.NitrogenCount();
    s.ring_jb_offset = sring.JohnsonBoveyLobeOffset();
    s.ring_aromaticity = static_cast<int>(sring.Aromaticity());
    s.ring_size = sring.RingSizeValue();
    s.ring_fused = sring.IsFused();

    // Self/bonded: the source ring belongs to the target atom's own residue
    // (a backbone atom whose residue carries a ring — PRO pyrrolidine etc.).
    const model::QtAtom& target = body.run.protein->atom(st.atom);
    s.is_self_or_bonded = (target.residueIndex >= 0
                           && sring.parentResidueIndex == target.residueIndex);
}

// ── Bond attacher (anisotropic bond-midpoints KD; mirrors mcAttacher) ────────
void bondAttacher(const Body& body, const AtomState& st, const FrameResult& fr,
                  const RawSource& raw, SourceSlot& s) {
    if (raw.kind != SourceKind::Bond) return;
    if (raw.ref.entity_index < 0) {
        s.dipolar = kNaN;
        return;
    }
    const model::QtProtein& p = *body.run.protein;
    const model::QtBond& b = p.topology().bondAt(static_cast<std::size_t>(raw.ref.entity_index));
    if (b.atomIndexA < 0 || b.atomIndexB < 0) {
        s.dipolar = kNaN;
        return;
    }
    const Vec3 posA = verbs::pos(body, static_cast<std::size_t>(b.atomIndexA), st.frame);
    const Vec3 posB = verbs::pos(body, static_cast<std::size_t>(b.atomIndexB), st.frame);
    const Vec3 axis = posB - posA;
    const double axisNorm = axis.norm();
    if (!(axisNorm > 1e-9)) {
        s.dipolar = kNaN;
        return;
    }
    const Vec3 axisU = axis / axisNorm;
    const Vec3 disp = 0.5 * (posA + posB) - st.pos;
    const double r = disp.norm();
    if (!(r > 1e-6)) {  // a bond whose own atom is the target
        s.dipolar = kNaN;
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
}

// ── Charge attacher (FF14SB charge-sites KD; mirrors ChargeDipoleNeighborhood) ─
// Source = a charged atom from nearBackend(ChargeSites): raw.ref.entity_index is
// the source atom index. Fills q + the displacement in the local frame. A
// rejected charge (charges absent, self atom, self residue when exclude_residue,
// non-finite q, or coincident) is marked with source_atom_index = -1 so the
// broad source_filter drops it (charge slots legitimately carry dipolar==0, so
// the NaN-dipolar reject the ring/bond attachers use does not apply here).
// Curried over charge_source + exclude_residue.
Attacher makeChargeAttacher(const QString& charge_source, bool exclude_residue) {
    return [charge_source, exclude_residue](const Body& body, const AtomState& st,
                                            const FrameResult& fr, const RawSource& raw,
                                            SourceSlot& s) {
        if (raw.kind != SourceKind::Charge) return;  // not ours; leave the slot
        s.kind = SourceKind::Charge;
        s.source_atom_index = -1;  // rejected until proven valid (the filter sentinel)
        if (!body.catalog.has(ArrayId::Ff14sbCharge) || raw.ref.entity_index < 0) return;
        const model::QtProtein& p = *body.run.protein;
        const std::size_t srcAtom = static_cast<std::size_t>(raw.ref.entity_index);
        if (srcAtom >= p.atomCount() || srcAtom == st.atom) return;  // never self
        if (!body.catalog.present(body, ArrayId::Ff14sbCharge, srcAtom, st.frame)) return;
        const model::QtAtom& target = p.atom(st.atom);
        const model::QtAtom& src = p.atom(srcAtom);
        if (exclude_residue && target.residueIndex >= 0
            && src.residueIndex == target.residueIndex)
            return;  // through-space: drop the target's own residue
        const double q = body.catalog.value(body, ArrayId::Ff14sbCharge, srcAtom, st.frame);
        if (!std::isfinite(q)) return;
        const Vec3 dispLab =
            body.run.conformation->atomPosition(st.frame, srcAtom) - st.pos;
        const double r = dispLab.norm();
        if (!(r > 1e-9)) return;

        s.source_charge_source = charge_source;
        s.source_q_e = q;
        s.source_atom_index = src.atomIndex;
        s.source_residue_index = src.residueIndex;
        s.source_element = static_cast<int>(src.element);
        if (validResidue(p, src.residueIndex)) {
            const model::QtResidue& sr = p.residue(static_cast<std::size_t>(src.residueIndex));
            s.source_residue_number = sr.address.residueNumber;
            s.source_amino_acid = static_cast<int>(sr.aminoAcid);
        }
        s.disp_local = fr.frame.is_valid ? fr.frame.ToLocal(dispLab) : dispLab;
        s.r = r;
    };
}

}  // namespace

BroadRelationship MakeBroadBackboneRelationship(double ring_cutoff_A, double bond_cutoff_A,
                                                double charge_cutoff_A,
                                                const QString& charge_source,
                                                bool exclude_residue) {
    BroadRelationship brel;
    brel.ring_cutoff_A = ring_cutoff_A;
    brel.bond_cutoff_A = bond_cutoff_A;
    brel.charge_cutoff_A = charge_cutoff_A;
    brel.charge_source = charge_source;

    Relationship& rel = brel.rel;
    rel.name = QStringLiteral("broad_backbone");

    // Stratum = EVERY backbone atom (N/CA/C/O/HN/HA, every residue), via the
    // typed BackboneRole predicate (chemistry, not a name scan). IsBackbone()
    // covers N/CA/C/O/HN and non-GLY HA; the GLY HA2/HA3 pair carries
    // Locant::Alpha + BackboneRole::None, so add IsAnyAlphaHydrogen explicitly.
    rel.stratum = atomsWhere([](const model::QtAtom& a) {
        return a.IsBackbone() || a.IsAnyAlphaHydrogen();
    });

    rel.frame_fn = &backboneFrameFn;

    // Heterogeneous selector LIST — the whole point of the engine. The GENERAL
    // KD backends (ring-centres / bond-midpoints / charge-sites), each curried
    // with its recorded cutoff. NOT the aromatic-H slots backend.
    rel.selectors = {
        nearBackend(CloudKind::RingCenters, ring_cutoff_A),
        nearBackend(CloudKind::BondMidpoints, bond_cutoff_A),
        nearBackend(CloudKind::ChargeSites, charge_cutoff_A),
    };

    // Per-mechanism attachers (a LIST — each branches on the typed source kind).
    // Ring + bond write geometry/identity; charge writes q + disp_local + the
    // FIELD inputs, curried over charge_source + exclude_residue. This is the
    // heterogeneous-attacher composition SURFACE_DESIGN's "more lambdas, no
    // special case" calls for.
    rel.attachers = {&ringAttacher, &bondAttacher,
                     makeChargeAttacher(charge_source, exclude_residue)};

    // Drop geometrically/chemically rejected sources before they become rows
    // (the cell's `continue`). Ring/bond rejects carry NaN dipolar; a rejected
    // charge carries source_atom_index < 0 (the attacher's sentinel — charge
    // slots legitimately have dipolar==0, so NaN can't be the charge reject).
    rel.source_filter = [](const SourceSlot& s) {
        if (s.kind == SourceKind::Charge) return s.source_atom_index >= 0;
        return std::isfinite(s.dipolar);
    };

    rel.target_fn = [](const Body& body, std::size_t atom, std::size_t orig,
                       const LocalFrame& f) { return BuildTarget(body.run, atom, orig, f); };
    // No producer bare-kernel cross-check for the broad case (heterogeneous; no
    // single producer kernel covers all backbone atoms × all mechanisms).
    rel.bare_kernel_fn = {};

    // The broad reducer: per-mechanism dipolar sums + the local Coulomb FIELD.
    brel.broad_reducer = [ring_cutoff_A, bond_cutoff_A, charge_cutoff_A, charge_source](
                             const Body& /*body*/, std::size_t /*atom*/, const FrameResult& /*fr*/,
                             const std::vector<SourceSlot>& sources) {
        return ReduceBroadBackboneSources(sources, ring_cutoff_A, bond_cutoff_A,
                                          charge_cutoff_A, charge_source);
    };

    return brel;
}

std::size_t RunBroadBackbone(const BroadRelationship& brel, const Body& body,
                             BroadBackboneSink& sink) {
    const RunData& run = body.run;
    const Relationship& rel = brel.rel;

    // The stratum closure (currying applied once).
    const std::vector<std::size_t> stratum = rel.stratum(body);
    qCInfo(cBroad).noquote() << "broad_backbone stratum atoms=" << stratum.size()
                             << "| ring_cut=" << brel.ring_cutoff_A
                             << "| bond_cut=" << brel.bond_cutoff_A
                             << "| charge_cut=" << brel.charge_cutoff_A
                             << "| charge_source=" << brel.charge_source
                             << "| dft rows=" << run.frameMap.dftRows().size();

    std::size_t cases = 0;
    // The SAME (atom, frame) traversal + curried-closure protocol as
    // RelationshipEngine::RunRelationship (stratum → frame_fn → selectors →
    // attachers → source_filter). It diverges ONLY in the reducer's output shape
    // (BroadAggregate, not the scalar-sum AggregateResult) and the sink
    // (BroadBackboneSink's two-kind target-once carrier, not RecordSink). That
    // divergence is the precise latent narrowness the broad case surfaced.
    for (std::size_t row : run.frameMap.dftRows()) {
        const std::size_t orig = run.frameMap.originalIndex(row);
        for (std::size_t atom : stratum) {
            const FrameResult fr = rel.frame_fn(body, atom, row);

            NeighborhoodRecord rec;
            FillIdentity(rec, run, atom, row, rel.name, fr.frame);
            rec.frame_anchor_atom_index = fr.anchor_atom_index;
            rec.target = rel.target_fn(body, atom, orig, fr.frame);

            const AtomState st = verbs::at(body, atom, row);

            // src = flatten(sel(...) for sel in selectors) — heterogeneous list.
            std::vector<RawSource> sources;
            for (const SourceSelector& sel : rel.selectors) {
                std::vector<RawSource> part = sel(body, atom, row);
                sources.insert(sources.end(), part.begin(), part.end());
            }

            // Inner map: every attacher on every source (each branches on the
            // typed source kind), then the source_filter. Identical to the engine.
            rec.sources.reserve(sources.size());
            for (const RawSource& raw : sources) {
                SourceSlot slot;
                for (const Attacher& attach : rel.attachers) attach(body, st, fr, raw, slot);
                if (rel.source_filter && !rel.source_filter(slot)) continue;
                rec.sources.push_back(slot);
            }

            // rec = broad_reducer(SourceSet): per-mechanism sums + the local field.
            const BroadAggregate agg = brel.broad_reducer(body, atom, fr, rec.sources);
            sink.Write(rec, agg);
            ++cases;
        }
    }
    qCInfo(cBroad).noquote() << "broad_backbone cases=" << cases;
    return cases;
}

}  // namespace h5reader::rediscover
