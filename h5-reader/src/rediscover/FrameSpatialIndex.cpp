#include "FrameSpatialIndex.h"

#include "../model/Conformation.h"
#include "../model/QtProtein.h"

namespace h5reader::rediscover {

namespace {

// The McConnell anisotropic-bond categories (DESIGN.md): peptide C=O, peptide
// C–N, sidechain C=O, aromatic. Other categories (backbone-other, disulfide,
// sidechain-other) carry no through-space anisotropy term here.
bool isAnisotropicCategory(model::BondCategory c) {
    switch (c) {
    case model::BondCategory::PeptideCO:
    case model::BondCategory::PeptideCN:
    case model::BondCategory::SidechainCO:
    case model::BondCategory::Aromatic:
        return true;
    default:
        return false;
    }
}

}  // namespace

FrameSpatialIndex::FrameSpatialIndex(const model::QtProtein& protein,
                                     const model::Conformation& conf, std::size_t frame) {
    const model::QtTopology& topo = protein.topology();
    bonds_.reserve(topo.bondCount());
    for (std::size_t bi = 0; bi < topo.bondCount(); ++bi) {
        const model::QtBond& b = topo.bondAt(bi);
        if (!isAnisotropicCategory(b.category)) continue;
        if (b.atomIndexA < 0 || b.atomIndexB < 0) continue;

        AnisoBond ab;
        ab.bondIndex = b.bondIndex;
        ab.atomIndexA = b.atomIndexA;
        ab.atomIndexB = b.atomIndexB;
        ab.category = static_cast<int>(b.category);
        ab.order = static_cast<int>(b.order);
        ab.elemA = static_cast<int>(protein.atom(static_cast<std::size_t>(b.atomIndexA)).element);
        ab.elemB = static_cast<int>(protein.atom(static_cast<std::size_t>(b.atomIndexB)).element);
        ab.posA = conf.atomPosition(frame, static_cast<std::size_t>(b.atomIndexA));
        ab.posB = conf.atomPosition(frame, static_cast<std::size_t>(b.atomIndexB));
        ab.midpoint = 0.5 * (ab.posA + ab.posB);
        bonds_.push_back(ab);
    }

    tree_ = std::make_unique<KdTree>(3, *this);
    tree_->buildIndex();
}

std::vector<std::size_t> FrameSpatialIndex::Within(const Vec3& query, double cutoff) const {
    std::vector<std::size_t> out;
    if (!tree_ || bonds_.empty()) return out;
    const double query_pt[3] = {query.x(), query.y(), query.z()};
    // L2_Simple_Adaptor distances are SQUARED — pass the squared radius.
    const double r2 = cutoff * cutoff;
    std::vector<nanoflann::ResultItem<std::size_t, double>> hits;
    const std::size_t nHits = tree_->radiusSearch(query_pt, r2, hits);  // count == hits.size()
    (void)nHits;  // hits filled by reference; capture the [[nodiscard]] return
    out.reserve(hits.size());
    for (const auto& h : hits) out.push_back(h.first);
    return out;
}

}  // namespace h5reader::rediscover
