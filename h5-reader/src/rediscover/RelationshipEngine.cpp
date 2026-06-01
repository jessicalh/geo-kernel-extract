#include "RelationshipEngine.h"

#include "ExtractionSupport.h"
#include "Verbs.h"

#include "../model/QtProtein.h"

namespace h5reader::rediscover {

std::size_t RunRelationship(const Relationship& rel, const Body& body, RecordSink& sink) {
    const RunData& run = body.run;

    // The stratum closure: captured body arrives here (currying applied once).
    const std::vector<std::size_t> stratum = rel.stratum(body);

    std::size_t cases = 0;

    // Outer map: the (atom, frame) index space, DFT-row-outer × atom-inner (the
    // cache-friendly schedule; order-free by referential transparency, picked to
    // match the oracle's emission order so the gate's file compare is exact).
    for (std::size_t row : run.frameMap.dftRows()) {
        const std::size_t orig = run.frameMap.originalIndex(row);

        for (std::size_t atom : stratum) {
            // lf = frame_fn(body, atom, frame): the curried local-frame closure.
            const FrameResult fr = rel.frame_fn(body, atom, row);

            NeighborhoodRecord rec;
            FillIdentity(rec, run, atom, row, rel.name, fr.frame);
            rec.frame_anchor_atom_index = fr.anchor_atom_index;

            // The producer's bare-kernel cross-check (optional curried closure).
            if (rel.bare_kernel_fn) rel.bare_kernel_fn(body, atom, row, rec);

            // The DFT target, in the recorded local frame (curried target closure).
            rec.target = rel.target_fn(body, atom, orig, fr.frame);

            const AtomState st = verbs::at(body, atom, row);

            // src = flatten( sel(body, atom, frame) for sel in selectors ):
            // the heterogeneous source list, each selector a curried closure.
            std::vector<RawSource> sources;
            for (const SourceSelector& sel : rel.selectors) {
                std::vector<RawSource> part = sel(body, atom, row);
                sources.insert(sources.end(), part.begin(), part.end());
            }

            // The per-atom classifier context (built once per atom; the state the
            // per-source classifier closes over — the own-ring / own-atom sets for
            // ring_current, nothing for mcconnell). SURFACE_DESIGN's SourceClassifier.
            ClassifierContext cctx;
            if (rel.classifier_prep) cctx = rel.classifier_prep(body, atom);

            // Inner map: for each source, run every attacher (geometry ⊕ identity
            // ⊕ kernel — each writes its fields into the slot), then the classifier
            // (validity / key). The attachers are curried over the frame closure;
            // the source + atom-state arrive here. The classifier runs before the
            // slot is recorded, so is_self_or_bonded lands on the per-source row —
            // matching the procedural cell, which sets it inline before pushing.
            rec.sources.reserve(sources.size());
            for (const RawSource& raw : sources) {
                SourceSlot slot;
                for (const Attacher& attach : rel.attachers) attach(body, st, fr, raw, slot);
                // A source_filter reject is the cell's geometric `continue`: the
                // slot never becomes a row (no NaN source row emitted).
                if (rel.source_filter && !rel.source_filter(slot)) continue;
                if (rel.classifier) rel.classifier(body, cctx, slot);
                rec.sources.push_back(slot);
            }

            // rec = reducer(SourceSet): the state-carrying fold over the attached
            // sources → the aggregated-row parameters (sumAll / sumValid / nValid
            // / per-type / cutoff). Identical to the procedural aggregation.
            const AggregateResult agg = rel.reducer(body, atom, rec.sources);

            // sink.fold(rec): stream both row kinds (untouched RecordSink).
            sink.WriteSourceRows(rec);
            sink.WriteAggregatedRow(rec, agg.sum_all, agg.sum_valid, agg.n_valid, agg.per_type,
                                    agg.cutoff_A);
            ++cases;
        }
    }
    return cases;
}

}  // namespace h5reader::rediscover
