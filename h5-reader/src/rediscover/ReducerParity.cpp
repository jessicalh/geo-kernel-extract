#include "ReducerParity.h"

#include "ComposedRelationships.h"
#include "Verbs.h"

#include "../model/QtBond.h"

#include <algorithm>
#include <limits>
#include <unordered_set>

namespace h5reader::rediscover {
namespace {

const model::BondCategory kMcCats[] = {
    model::BondCategory::PeptideCO,
    model::BondCategory::PeptideCN,
    model::BondCategory::SidechainCO,
    model::BondCategory::Aromatic,
};
constexpr int kMcCatCount = 4;

int mcCatColumn(int categoryOrd) {
    for (int i = 0; i < kMcCatCount; ++i) {
        if (static_cast<int>(kMcCats[i]) == categoryOrd) return i;
    }
    return -1;
}

bool containsAtom(const std::unordered_set<std::size_t>& atoms, std::size_t atom) {
    return atoms.find(atom) != atoms.end();
}

ReducerAggregateSnapshot fromAggregateResult(const AggregateResult& agg, bool source_present) {
    ReducerAggregateSnapshot out;
    out.present = source_present;
    out.sum_all = agg.sum_all;
    out.sum_valid = agg.sum_valid;
    out.n_valid = agg.n_valid;
    out.per_type = agg.per_type;
    return out;
}

ReducerAggregateSnapshot composedAggregate(const Relationship& rel,
                                           const Body& body,
                                           std::size_t atom,
                                           std::size_t frame) {
    const FrameResult fr = rel.frame_fn(body, atom, frame);
    const AtomState st = verbs::at(body, atom, frame);
    std::vector<RawSource> rawSources;
    for (const SourceSelector& sel : rel.selectors) {
        std::vector<RawSource> part = sel(body, atom, frame);
        rawSources.insert(rawSources.end(), part.begin(), part.end());
    }

    ClassifierContext cctx;
    if (rel.classifier_prep) cctx = rel.classifier_prep(body, atom);

    std::vector<SourceSlot> sourceSlots;
    sourceSlots.reserve(rawSources.size());
    for (const RawSource& raw : rawSources) {
        SourceSlot slot;
        for (const Attacher& attach : rel.attachers) attach(body, st, fr, raw, slot);
        if (rel.source_filter && !rel.source_filter(slot)) continue;
        if (rel.classifier) rel.classifier(body, cctx, slot);
        sourceSlots.push_back(slot);
    }
    return fromAggregateResult(rel.reducer(body, atom, sourceSlots), !sourceSlots.empty());
}

ReducerAggregateSnapshot perAtomRingAggregate(const Body& body,
                                              std::size_t atom,
                                              std::size_t frame,
                                              const PerAtomSubstrateConfig& config,
                                              const LocalFrame& localFrame) {
    ReducerAggregateSnapshot out;
    out.per_type.assign(8, 0.0);
    const std::vector<PairContribution> pairs =
        PerAtomRowPairContributions(body, atom, frame, config, localFrame);
    for (const PairContribution& p : pairs) {
        if (p.mechanism != QStringLiteral("ring_jb") || !std::isfinite(p.dipolar))
            continue;
        out.present = true;
        out.sum_all += p.dipolar;
        if (!(p.pointer_flags & SelfOrBondedFlag)) {
            out.sum_valid += p.dipolar;
            ++out.n_valid;
            if (p.source_category_ord >= 0
                && static_cast<std::size_t>(p.source_category_ord) < out.per_type.size()) {
                out.per_type[static_cast<std::size_t>(p.source_category_ord)] += p.dipolar;
            }
        }
    }
    return out;
}

ReducerAggregateSnapshot perAtomMcAggregate(const Body& body,
                                            std::size_t atom,
                                            std::size_t frame,
                                            const PerAtomSubstrateConfig& config,
                                            const LocalFrame& localFrame) {
    ReducerAggregateSnapshot out;
    out.per_type.assign(4, 0.0);
    const std::vector<PairContribution> pairs =
        PerAtomRowPairContributions(body, atom, frame, config, localFrame);
    for (const PairContribution& p : pairs) {
        if (p.mechanism != QStringLiteral("mc_lit_valid") || !std::isfinite(p.dipolar))
            continue;
        out.present = true;
        out.sum_all += p.dipolar;
        if (!(p.pointer_flags & SelfOrBondedFlag)) {
            out.sum_valid += p.dipolar;
            ++out.n_valid;
            const int cc = mcCatColumn(p.source_category_ord);
            if (cc >= 0) {
                out.per_type[static_cast<std::size_t>(cc)] += p.dipolar;
            }
        }
    }
    return out;
}

}  // namespace

ReducerParityStats AuditReducerParity(const Body& body,
                                      const std::vector<std::size_t>& atoms,
                                      const std::vector<std::size_t>& frames,
                                      const PerAtomSubstrateConfig& config,
                                      const ReducerParityOptions& options) {
    ReducerParityStats stats;
    if (!body.run.protein) return stats;

    const Relationship ring = MakeRingCurrentRelationship();
    const Relationship mc = MakeMcConnellRelationship(config.bond_cutoff_A);
    const std::vector<std::size_t> ringStratum = ring.stratum(body);
    const std::vector<std::size_t> mcStratum = mc.stratum(body);
    const std::unordered_set<std::size_t> ringAtoms(ringStratum.begin(), ringStratum.end());
    const std::unordered_set<std::size_t> mcAtoms(mcStratum.begin(), mcStratum.end());

    for (std::size_t frame : frames) {
        for (std::size_t atom : atoms) {
            if (options.max_cases != 0 && stats.cases_checked >= options.max_cases)
                return stats;

            bool anyChecked = false;
            if (containsAtom(ringAtoms, atom)) {
                const ReducerAggregateSnapshot ringComposed =
                    composedAggregate(ring, body, atom, frame);
                const LocalFrame ringFrame = ring.frame_fn(body, atom, frame).frame;
                const ReducerAggregateSnapshot ringPerAtom =
                    perAtomRingAggregate(body, atom, frame, config, ringFrame);
                CompareReducerAggregateSnapshots(ringComposed, ringPerAtom,
                                                 QStringLiteral("ring atom=%1 frame=%2").arg(atom).arg(frame),
                                                 options, &stats.mismatches);
                ++stats.ring_cases_checked;
                anyChecked = true;
            }

            if (containsAtom(mcAtoms, atom)) {
                const ReducerAggregateSnapshot mcComposed =
                    composedAggregate(mc, body, atom, frame);
                const LocalFrame mcFrame = mc.frame_fn(body, atom, frame).frame;
                const ReducerAggregateSnapshot mcPerAtom =
                    perAtomMcAggregate(body, atom, frame, config, mcFrame);
                CompareReducerAggregateSnapshots(mcComposed, mcPerAtom,
                                                 QStringLiteral("mc atom=%1 frame=%2").arg(atom).arg(frame),
                                                 options, &stats.mismatches);
                ++stats.mc_cases_checked;
                anyChecked = true;
            }
            if (anyChecked) ++stats.cases_checked;
        }
    }
    return stats;
}

}  // namespace h5reader::rediscover
