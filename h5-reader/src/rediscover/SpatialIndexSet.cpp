#include "SpatialIndexSet.h"

#include "RunData.h"

#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"
#include "../model/QtBond.h"
#include "../model/QtProtein.h"

#include <cmath>
#include <stdexcept>

namespace h5reader::rediscover {

namespace {

int ord(CloudKind kind) { return static_cast<int>(kind); }

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

CloudTree::CloudTree(CloudKind kind, std::vector<CloudPoint> points)
    : kind_(kind), points_(std::move(points)) {
    tree_ = std::make_unique<KdTree>(3, *this);
    tree_->buildIndex();
}

std::vector<SourceRef> CloudTree::near(const Vec3& query, double cutoff) const {
    std::vector<SourceRef> out;
    if (!tree_ || points_.empty() || !(cutoff >= 0.0)) return out;
    const double query_pt[3] = {query.x(), query.y(), query.z()};
    const double r2 = cutoff * cutoff;
    std::vector<nanoflann::ResultItem<std::size_t, double>> hits;
    const std::size_t nHits = tree_->radiusSearch(query_pt, r2, hits);
    (void)nHits;
    out.reserve(hits.size());
    for (const auto& h : hits) out.push_back(points_[h.first].ref);
    return out;
}

std::vector<SourceRef> CloudTree::range(const Vec3& query, double innerCutoff,
                                        double outerCutoff) const {
    std::vector<SourceRef> out;
    if (!(outerCutoff >= innerCutoff) || innerCutoff < 0.0) return out;
    const double inner2 = innerCutoff * innerCutoff;
    for (SourceRef ref : near(query, outerCutoff)) {
        const Vec3 d = points_[static_cast<std::size_t>(ref.cloud_index)].point - query;
        if (d.squaredNorm() >= inner2) out.push_back(ref);
    }
    return out;
}

std::vector<SourceRef> CloudTree::nearAtom(std::size_t atom, double cutoff) const {
    if (atom >= points_.size()) return {};
    return near(points_[atom].point, cutoff);
}

std::vector<SourceRef> CloudTree::rangeAtom(std::size_t atom, double innerCutoff,
                                            double outerCutoff) const {
    if (atom >= points_.size()) return {};
    return range(points_[atom].point, innerCutoff, outerCutoff);
}

SpatialIndexSet::SpatialIndexSet(const RunData& run) {
    if (!run.protein || !run.conformation) return;
    const model::QtProtein& p = *run.protein;
    const model::Conformation& conf = *run.conformation;
    const model::QtTopology& topo = p.topology();
    frameCount_ = conf.frameCount();
    for (auto& v : trees_) v.reserve(frameCount_);

    for (std::size_t frame = 0; frame < frameCount_; ++frame) {
        std::vector<CloudPoint> atoms;
        atoms.reserve(p.atomCount());
        for (std::size_t ai = 0; ai < p.atomCount(); ++ai) {
            atoms.push_back({conf.atomPosition(frame, ai),
                             {CloudKind::Atoms, static_cast<int32_t>(atoms.size()),
                              static_cast<int32_t>(ai)}});
        }

        std::vector<CloudPoint> bonds;
        bonds.reserve(topo.bondCount());
        std::vector<CloudPoint> allBonds;
        allBonds.reserve(topo.bondCount());
        for (std::size_t bi = 0; bi < topo.bondCount(); ++bi) {
            const model::QtBond& b = topo.bondAt(bi);
            if (b.atomIndexA < 0 || b.atomIndexB < 0) continue;
            const Vec3 a = conf.atomPosition(frame, static_cast<std::size_t>(b.atomIndexA));
            const Vec3 c = conf.atomPosition(frame, static_cast<std::size_t>(b.atomIndexB));
            const Vec3 midpoint = 0.5 * (a + c);
            allBonds.push_back({midpoint,
                                {CloudKind::AllBondMidpoints,
                                 static_cast<int32_t>(allBonds.size()),
                                 static_cast<int32_t>(bi)}});
            if (isAnisotropicCategory(b.category)) {
                bonds.push_back({midpoint,
                                 {CloudKind::BondMidpoints,
                                  static_cast<int32_t>(bonds.size()),
                                  static_cast<int32_t>(bi)}});
            }
        }

        std::vector<CloudPoint> rings;
        rings.reserve(p.ringCount());
        for (std::size_t ri = 0; ri < p.ringCount(); ++ri) {
            const model::RingGeometry g = model::RingGeometryAt(conf, ri, frame);
            rings.push_back({g.center,
                             {CloudKind::RingCenters, static_cast<int32_t>(rings.size()),
                              static_cast<int32_t>(ri)}});
        }

        std::vector<CloudPoint> charges;
        charges.reserve(p.atomCount());
        for (std::size_t ai = 0; ai < p.atomCount(); ++ai) {
            if (!p.atom(ai).hasPartialCharge) continue;
            charges.push_back({conf.atomPosition(frame, ai),
                               {CloudKind::ChargeSites, static_cast<int32_t>(charges.size()),
                                static_cast<int32_t>(ai)}});
        }

        trees_[ord(CloudKind::Atoms)].emplace_back(CloudKind::Atoms, std::move(atoms));
        trees_[ord(CloudKind::BondMidpoints)].emplace_back(CloudKind::BondMidpoints, std::move(bonds));
        trees_[ord(CloudKind::RingCenters)].emplace_back(CloudKind::RingCenters, std::move(rings));
        trees_[ord(CloudKind::ChargeSites)].emplace_back(CloudKind::ChargeSites, std::move(charges));
        trees_[ord(CloudKind::AllBondMidpoints)].emplace_back(CloudKind::AllBondMidpoints,
                                                              std::move(allBonds));
    }
}

const CloudTree& SpatialIndexSet::tree(CloudKind kind, std::size_t frame) const {
    const auto& v = trees_[ord(kind)];
    if (frame >= v.size()) throw std::out_of_range("SpatialIndexSet::tree frame out of range");
    return v[frame];
}

std::vector<SourceRef> SpatialIndexSet::near(CloudKind kind, std::size_t frame,
                                             const Vec3& query, double cutoff) const {
    return tree(kind, frame).near(query, cutoff);
}

std::vector<SourceRef> SpatialIndexSet::range(CloudKind kind, std::size_t frame,
                                              const Vec3& query, double innerCutoff,
                                              double outerCutoff) const {
    return tree(kind, frame).range(query, innerCutoff, outerCutoff);
}

}  // namespace h5reader::rediscover
