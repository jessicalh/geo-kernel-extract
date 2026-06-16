#include "AnalysisAtom.h"

#include "ExtractionSupport.h"
#include "SphericalBasis.h"
#include "Verbs.h"

#include "../io/QtTrajectoryH5.h"
#include "../io/QtFieldCatalog.gen.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtRing.h"
#include "../model/QtTopology.h"

#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/statistics/linear_regression.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>

#include <QByteArray>
#include <QStringList>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <unordered_map>
#include <utility>

namespace h5reader::rediscover {
namespace {

namespace bmstats = boost::math::statistics;

constexpr std::size_t kRollingWindow = 32;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kCoulombKe = 14.3996;

double norm3(const std::array<double, 3>& values) {
    double s = 0.0;
    for (double v : values) s += v * v;
    return std::sqrt(s);
}

double norm5(const std::array<double, 5>& values) {
    double s = 0.0;
    for (double v : values) s += v * v;
    return std::sqrt(s);
}

double tensorFrobenius(const Mat3& m) {
    return std::sqrt((m.array() * m.array()).sum());
}

bool finite(double v) {
    return std::isfinite(v);
}

std::string keyString(const QString& q) {
    const QByteArray bytes = q.toUtf8();
    return std::string(bytes.constData(), static_cast<std::size_t>(bytes.size()));
}

QString qjoin(const QStringList& parts) {
    return parts.join(QStringLiteral("|"));
}

bool validResidue(const model::QtProtein& protein, int32_t residue_index) {
    return residue_index >= 0
           && static_cast<std::size_t>(residue_index) < protein.residueCount();
}

QString residueLabel(const model::QtProtein& protein, int32_t residue_index) {
    if (!validResidue(protein, residue_index)) return QStringLiteral("residue=unknown");
    const model::QtResidue& residue = protein.residue(static_cast<std::size_t>(residue_index));
    return QStringLiteral("residue=%1:aa=%2")
        .arg(residue.address.residueNumber)
        .arg(static_cast<int>(residue.aminoAcid));
}

QString residueSetLabel(const model::QtProtein& protein, const std::set<int32_t>& residues) {
    if (residues.empty()) return QStringLiteral("residue=unknown");
    QStringList parts;
    for (int32_t residue_index : residues) parts << residueLabel(protein, residue_index);
    return parts.join(QStringLiteral("+"));
}

QString targetChemistryCategory(const Body& body, std::size_t atom, std::size_t row) {
    int element = -1;
    if (body.run.protein && atom < body.run.protein->atomCount())
        element = static_cast<int>(body.run.protein->atom(atom).element);

    int hybridisation = -1;
    if (body.catalog.present(body, io::FieldKind::EnrichmentHybridisation, atom, row, -1)) {
        const std::optional<double> v =
            body.catalog.value(body, io::FieldKind::EnrichmentHybridisation, atom, row, -1);
        if (v && std::isfinite(*v)) hybridisation = static_cast<int>(*v);
    }

    return qjoin({QStringLiteral("hyb=%1").arg(hybridisation),
                  QStringLiteral("element=%1").arg(element)});
}

bool ringContainsTargetHeavy(const Body& body, std::size_t target_atom, int32_t ring_id) {
    if (!body.run.protein || ring_id < 0) return false;
    const model::QtProtein& protein = *body.run.protein;
    if (static_cast<std::size_t>(ring_id) >= protein.topology().ringCount()) return false;
    const std::size_t heavy = verbs::heavyParent(body, target_atom);
    const model::QtRing& ring = protein.topology().ringAt(static_cast<std::size_t>(ring_id));
    return std::find(ring.atomIndices.begin(), ring.atomIndices.end(),
                     static_cast<int32_t>(heavy)) != ring.atomIndices.end();
}

bool bondHasTargetEndpoint(const Body& body, std::size_t target_atom, int32_t bond_id) {
    if (!body.run.protein || bond_id < 0) return false;
    const model::QtProtein& protein = *body.run.protein;
    if (static_cast<std::size_t>(bond_id) >= protein.topology().bondCount()) return false;
    const model::QtBond& bond = protein.topology().bondAt(static_cast<std::size_t>(bond_id));
    return bond.atomIndexA == static_cast<int32_t>(target_atom)
           || bond.atomIndexB == static_cast<int32_t>(target_atom);
}

QString relationshipScope(const Body& body, std::size_t target_atom, const PairContribution& pair) {
    const bool selfOrBonded = (pair.pointer_flags & SelfOrBondedFlag) != 0;
    if (!selfOrBonded) return QStringLiteral("scope=through_space");

    const bool self =
        (pair.source_atom_index == static_cast<int32_t>(target_atom))
        || (pair.source_kind == QStringLiteral("ring_center")
            && ringContainsTargetHeavy(body, target_atom, pair.source_id))
        || (pair.source_kind == QStringLiteral("bond_midpoint")
            && bondHasTargetEndpoint(body, target_atom, pair.source_id));
    if (self) return QStringLiteral("scope=self");
    if (pair.pointer_flags & NearFieldFlag) return QStringLiteral("scope=bonded_near_field");
    return QStringLiteral("scope=bonded");
}

QString contactedResidueLabel(const Body& body, const PairContribution& pair) {
    if (!body.run.protein) return QStringLiteral("residue=unknown");
    const model::QtProtein& protein = *body.run.protein;
    if (pair.source_atom_index >= 0
        && static_cast<std::size_t>(pair.source_atom_index) < protein.atomCount()) {
        return residueLabel(protein, protein.atom(static_cast<std::size_t>(pair.source_atom_index)).residueIndex);
    }

    std::set<int32_t> residues;
    if (pair.source_kind == QStringLiteral("bond_midpoint") && pair.source_id >= 0
        && static_cast<std::size_t>(pair.source_id) < protein.topology().bondCount()) {
        const model::QtBond& bond = protein.topology().bondAt(static_cast<std::size_t>(pair.source_id));
        if (bond.atomIndexA >= 0 && static_cast<std::size_t>(bond.atomIndexA) < protein.atomCount())
            residues.insert(protein.atom(static_cast<std::size_t>(bond.atomIndexA)).residueIndex);
        if (bond.atomIndexB >= 0 && static_cast<std::size_t>(bond.atomIndexB) < protein.atomCount())
            residues.insert(protein.atom(static_cast<std::size_t>(bond.atomIndexB)).residueIndex);
    } else if (pair.source_kind == QStringLiteral("ring_center") && pair.source_id >= 0
               && static_cast<std::size_t>(pair.source_id) < protein.topology().ringCount()) {
        const model::QtRing& ring = protein.topology().ringAt(static_cast<std::size_t>(pair.source_id));
        for (int32_t atom : ring.atomIndices) {
            if (atom >= 0 && static_cast<std::size_t>(atom) < protein.atomCount())
                residues.insert(protein.atom(static_cast<std::size_t>(atom)).residueIndex);
        }
    }
    return residueSetLabel(protein, residues);
}

std::size_t hashCombine(std::size_t seed, std::size_t value) {
    return seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
}

struct RelationshipKeyBase {
    std::string chemistry;
    std::string contacted_residue;
    std::string scope;
    std::string mechanism;
    std::string source_kind;
    int32_t source_id = -1;
    int source_category_ord = -1;
};

struct RelationshipFoldKey {
    RelationshipKeyBase base;
    std::string facet;

    bool operator==(const RelationshipFoldKey& other) const {
        return base.chemistry == other.base.chemistry
               && base.contacted_residue == other.base.contacted_residue
               && base.scope == other.base.scope
               && base.mechanism == other.base.mechanism
               && base.source_kind == other.base.source_kind
               && base.source_id == other.base.source_id
               && base.source_category_ord == other.base.source_category_ord
               && facet == other.facet;
    }
};

struct RelationshipFoldKeyHash {
    std::size_t operator()(const RelationshipFoldKey& key) const {
        std::size_t seed = 0;
        const std::hash<std::string> hs;
        const std::hash<int> hi;
        seed = hashCombine(seed, hs(key.base.chemistry));
        seed = hashCombine(seed, hs(key.base.contacted_residue));
        seed = hashCombine(seed, hs(key.base.scope));
        seed = hashCombine(seed, hs(key.base.mechanism));
        seed = hashCombine(seed, hs(key.base.source_kind));
        seed = hashCombine(seed, hi(key.base.source_id));
        seed = hashCombine(seed, hi(key.base.source_category_ord));
        seed = hashCombine(seed, hs(key.facet));
        return seed;
    }
};

struct SelfFoldKey {
    std::string chemistry;
    std::string property;

    bool operator==(const SelfFoldKey& other) const {
        return chemistry == other.chemistry && property == other.property;
    }
};

struct SelfFoldKeyHash {
    std::size_t operator()(const SelfFoldKey& key) const {
        std::size_t seed = 0;
        const std::hash<std::string> hs;
        seed = hashCombine(seed, hs(key.chemistry));
        seed = hashCombine(seed, hs(key.property));
        return seed;
    }
};

RelationshipKeyBase relationshipKeyBase(const Body& body,
                                        std::size_t atom,
                                        std::size_t row,
                                        const PairContribution& pair) {
    return RelationshipKeyBase{
        keyString(targetChemistryCategory(body, atom, row)),
        keyString(contactedResidueLabel(body, pair)),
        keyString(relationshipScope(body, atom, pair)),
        keyString(pair.mechanism),
        keyString(pair.source_kind),
        pair.source_id,
        pair.source_category_ord,
    };
}

RelationshipKeyBase larsenRelationshipKeyBase(const Body& body,
                                              std::size_t atom,
                                              std::size_t row,
                                              const QString& source_kind,
                                              int32_t source_id) {
    return RelationshipKeyBase{
        keyString(targetChemistryCategory(body, atom, row)),
        keyString(QStringLiteral("residue=self")),
        keyString(QStringLiteral("scope=self")),
        keyString(QStringLiteral("larsen_hbond")),
        keyString(source_kind),
        source_id,
        static_cast<int>(source_id),
    };
}

QString relationshipLabel(const RelationshipFoldKey& key) {
    return qjoin({QStringLiteral("relationship"),
                  QString::fromStdString(key.base.chemistry),
                  QString::fromStdString(key.base.contacted_residue),
                  QString::fromStdString(key.base.scope),
                  QStringLiteral("mechanism=%1").arg(QString::fromStdString(key.base.mechanism)),
                  QStringLiteral("source_kind=%1").arg(QString::fromStdString(key.base.source_kind)),
                  QStringLiteral("source_id=%1").arg(key.base.source_id),
                  QStringLiteral("source_category=%1").arg(key.base.source_category_ord),
                  QStringLiteral("facet=%1").arg(QString::fromStdString(key.facet))});
}

SelfFoldKey selfFoldKey(const Body& body, std::size_t atom, std::size_t row, const QString& property) {
    return SelfFoldKey{keyString(targetChemistryCategory(body, atom, row)), keyString(property)};
}

QString selfLabel(const SelfFoldKey& key) {
    return qjoin({QStringLiteral("self_property"),
                  QString::fromStdString(key.chemistry),
                  QStringLiteral("residue=self"),
                  QStringLiteral("scope=self"),
                  QStringLiteral("property=%1").arg(QString::fromStdString(key.property))});
}

struct SeriesSample {
    std::uint32_t frame_slot = 0;
    double value = kNaN;
};

struct BatchStats {
    std::size_t n = 0;
    double mean_value = kNaN;
    double variance_value = kNaN;
    double covariance_with_sigma = kNaN;
    double slope_sigma_on_value = kNaN;
    double rolling_slope_sigma_on_value = kNaN;
    double interpolated_midpoint_value = kNaN;
};

struct SeriesOrgan {
    void fold(std::uint32_t frame_slot, double value) {
        if (!finite(value)) return;
        samples.push_back(SeriesSample{frame_slot, value});
    }

    std::size_t count() const { return samples.size(); }

    BatchStats deriveBatch(const std::vector<double>& sigma_by_frame) const {
        BatchStats out;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> t;
        x.reserve(samples.size());
        y.reserve(samples.size());
        t.reserve(samples.size());
        for (const SeriesSample& sample : samples) {
            if (sample.frame_slot >= sigma_by_frame.size()) continue;
            const double sigma = sigma_by_frame[sample.frame_slot];
            if (!finite(sample.value) || !finite(sigma)) continue;
            x.push_back(sample.value);
            y.push_back(sigma);
            t.push_back(static_cast<double>(sample.frame_slot));
        }
        out.n = x.size();
        if (x.empty()) return out;
        out.mean_value = bmstats::mean(x);
        out.variance_value = x.size() > 1 ? bmstats::variance(x) : 0.0;
        if (x.size() > 1) {
            try {
                out.covariance_with_sigma = bmstats::covariance(x, y);
                out.slope_sigma_on_value = bmstats::simple_ordinary_least_squares(x, y).second;
            } catch (const std::domain_error&) {
            }
            const std::size_t first = x.size() > kRollingWindow ? x.size() - kRollingWindow : 0;
            if (x.size() - first > 1) {
                std::vector<double> rx(x.begin() + static_cast<std::ptrdiff_t>(first), x.end());
                std::vector<double> ry(y.begin() + static_cast<std::ptrdiff_t>(first), y.end());
                try {
                    out.rolling_slope_sigma_on_value =
                        bmstats::simple_ordinary_least_squares(rx, ry).second;
                } catch (const std::domain_error&) {
                }
            }
        }
        if (x.size() >= 4) {
            try {
                boost::math::interpolators::pchip<std::vector<double>> smooth(std::move(t),
                                                                              std::move(x));
                const double mid = 0.5 * (samples.front().frame_slot + samples.back().frame_slot);
                out.interpolated_midpoint_value = smooth(mid);
            } catch (const std::domain_error&) {
            }
        }
        return out;
    }

    std::vector<SeriesSample> samples;
};

template <typename Key, typename Hash>
struct FoldBook {
    template <typename Labeler>
    void fold(const Key& key,
              std::uint32_t frame_slot,
              double sigma_iso,
              double value,
              Labeler labeler) {
        if (!finite(sigma_iso) || !finite(value)) return;
        auto it = organs.find(key);
        if (it == organs.end()) it = organs.emplace(key, SeriesOrgan{}).first;
        it->second.fold(frame_slot, value);
        if (!sample) {
            sample = std::make_pair(labeler(key), value);
            sample_key = key;
        }
    }

    std::size_t organCount() const { return organs.size(); }

    std::size_t foldedObservations() const {
        std::size_t out = 0;
        for (const auto& item : organs) out += item.second.count();
        return out;
    }

    BatchStats sampleStats(const std::vector<double>& sigma_by_frame) const {
        if (!sample_key) return {};
        const auto it = organs.find(*sample_key);
        if (it != organs.end()) return it->second.deriveBatch(sigma_by_frame);
        return {};
    }

    std::unordered_map<Key, SeriesOrgan, Hash> organs;
    std::optional<std::pair<QString, double>> sample;
    std::optional<Key> sample_key;
};

using RelationshipBook = FoldBook<RelationshipFoldKey, RelationshipFoldKeyHash>;
using SelfBook = FoldBook<SelfFoldKey, SelfFoldKeyHash>;

void foldSelfKey(SelfBook& book,
                 const SelfFoldKey& key,
                 std::uint32_t frame_slot,
                 double sigma,
                 double value) {
    book.fold(key, frame_slot, sigma, value, selfLabel);
}

void foldSelfProperty(SelfBook& book,
                      const Body& body,
                      std::size_t atom,
                      std::size_t row,
                      std::uint32_t frame_slot,
                      const QString& property,
                      double sigma,
                      double value) {
    foldSelfKey(book, selfFoldKey(body, atom, row, property), frame_slot, sigma, value);
}

void foldT2Components(SelfBook& book,
                      const Body& body,
                      std::size_t atom,
                      std::size_t row,
                      std::uint32_t frame_slot,
                      const QString& base,
                      double sigma,
                      const std::array<double, 5>& t2) {
    foldSelfProperty(book, body, atom, row, frame_slot, base + QStringLiteral(".norm"), sigma,
                     norm5(t2));
    for (std::size_t i = 0; i < t2.size(); ++i)
        foldSelfProperty(book, body, atom, row, frame_slot,
                         base + QStringLiteral(".c%1").arg(static_cast<int>(i)),
                         sigma, t2[i]);
}

void foldT1Components(SelfBook& book,
                      const Body& body,
                      std::size_t atom,
                      std::size_t row,
                      std::uint32_t frame_slot,
                      const QString& base,
                      double sigma,
                      const std::array<double, 3>& t1) {
    foldSelfProperty(book, body, atom, row, frame_slot, base + QStringLiteral(".norm"), sigma,
                     norm3(t1));
    for (std::size_t i = 0; i < t1.size(); ++i)
        foldSelfProperty(book, body, atom, row, frame_slot,
                         base + QStringLiteral(".c%1").arg(static_cast<int>(i)),
                         sigma, t1[i]);
}

}  // namespace

AnalysisAtom::AnalysisAtom(const Body& body, std::size_t atom_index)
    : body_(body), atom_index_(atom_index) {}

class NeighbourhoodAccumulatingAnalysisAtom::Impl {
public:
    explicit Impl(PerAtomSubstrateConfig config) : config_(std::move(config)) {}

    void observe(const Body& body, std::size_t atom, std::size_t h5_row) {
        const std::uint32_t frame_slot = static_cast<std::uint32_t>(frame_events_);
        ++frame_events_;
        AtomFrameRecord& record = ensureFrameRecord(frame_slot, h5_row);
        if (h5_row >= body.run.frameMap.frameCount()) return;
        const std::size_t original = body.run.frameMap.originalIndex(h5_row);
        record.original_frame = original;
        const DftTarget target = BuildTarget(body.run, atom, original, LocalFrame{});
        if (!target.present || !std::isfinite(target.total_decomp.T0)) return;
        ++dft_present_;
        const double sigma_iso = target.total_decomp.T0;
        record.dft_present = true;
        record.total_decomp = target.total_decomp;
        record.dia_decomp = target.dia_decomp;
        record.para_decomp = target.para_decomp;
        record.total_frobenius = tensorFrobenius(target.total_raw);
        if (frame_slot >= sigma_by_frame_.size())
            sigma_by_frame_.resize(static_cast<std::size_t>(frame_slot) + 1, kNaN);
        sigma_by_frame_[frame_slot] = sigma_iso;

        foldSelfProperties(body, atom, h5_row, frame_slot, record, target, sigma_iso);
        foldRelationships(body, atom, h5_row, frame_slot, record, sigma_iso);
        foldLarsenHbondRelationships(body, atom, h5_row, frame_slot, record, sigma_iso);
    }

    AnalysisAtomTruth diagnostics() const {
        AnalysisAtomTruth out;
        out.frame_events = frame_events_;
        out.dft_present = dft_present_;
        out.relationship_folds = relationship_book_.foldedObservations();
        out.self_property_folds = self_book_.foldedObservations();
        out.relationship_organs = relationship_book_.organCount();
        out.self_property_organs = self_book_.organCount();
        out.sigma_folds = sigma_folds_;
        out.dihedral_folds = dihedral_folds_;
        out.mopac_scalar_folds = mopac_scalar_folds_;
        out.efg_folds = efg_folds_;
        if (relationship_book_.sample) {
            out.sample_label = relationship_book_.sample->first;
            out.sample_value = relationship_book_.sample->second;
            out.has_sample = true;
            const BatchStats stats = relationship_book_.sampleStats(sigma_by_frame_);
            (void)stats;
        } else if (self_book_.sample) {
            out.sample_label = self_book_.sample->first;
            out.sample_value = self_book_.sample->second;
            out.has_sample = true;
            const BatchStats stats = self_book_.sampleStats(sigma_by_frame_);
            (void)stats;
        }
        return out;
    }

private:
    struct DihedralFrameValue {
        std::string name;
        double radians = kNaN;
        double sin_value = kNaN;
        double cos_value = kNaN;
        int fixed_bin = -1;
    };

    struct AtomFrameRecord {
        std::size_t h5_row = 0;
        std::size_t original_frame = 0;
        bool dft_present = false;
        model::SphericalTensor total_decomp;
        model::SphericalTensor dia_decomp;
        model::SphericalTensor para_decomp;
        double total_frobenius = kNaN;
        std::vector<DihedralFrameValue> dihedrals;
        std::vector<int32_t> partner_source_ids;
    };

    AtomFrameRecord& ensureFrameRecord(std::uint32_t frame_slot, std::size_t h5_row) {
        if (frame_slot >= frame_records_.size())
            frame_records_.resize(static_cast<std::size_t>(frame_slot) + 1);
        AtomFrameRecord& record = frame_records_[frame_slot];
        record.h5_row = h5_row;
        if (frame_slot >= sigma_by_frame_.size())
            sigma_by_frame_.resize(static_cast<std::size_t>(frame_slot) + 1, kNaN);
        return record;
    }

    void addPartnerSource(AtomFrameRecord& record, int32_t source_id) {
        if (source_id < 0) return;
        if (std::find(record.partner_source_ids.begin(), record.partner_source_ids.end(), source_id)
            == record.partner_source_ids.end())
            record.partner_source_ids.push_back(source_id);
    }

    void foldRelationshipFacet(const RelationshipKeyBase& base,
                               std::uint32_t frame_slot,
                               double sigma,
                               const QString& facet,
                               double value) {
        RelationshipFoldKey key{base, keyString(facet)};
        relationship_book_.fold(key, frame_slot, sigma, value, relationshipLabel);
    }

    void foldRelationshipTensorFacets(const RelationshipKeyBase& base,
                                      std::uint32_t frame_slot,
                                      double sigma,
                                      double t0,
                                      const std::array<double, 5>& t2) {
        foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("kernel_T0"), t0);
        foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("kernel_T2_norm"), norm5(t2));
        for (std::size_t i = 0; i < t2.size(); ++i) {
            foldRelationshipFacet(base, frame_slot, sigma,
                                  QStringLiteral("kernel_T2_c%1").arg(static_cast<int>(i)),
                                  t2[i]);
        }
        foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("contribution"), norm5(t2));
    }

    void foldSelfFacet(const Body& body,
                       std::size_t atom,
                       std::size_t row,
                       std::uint32_t frame_slot,
                       const QString& property,
                       double sigma,
                       double value) {
        foldSelfProperty(self_book_, body, atom, row, frame_slot, property, sigma, value);
    }

    std::optional<model::SphericalTensor> mopacPerSourceEfg(const Body& body,
                                                           std::size_t source_atom,
                                                           std::size_t row,
                                                           const PairContribution& pair) const {
        if (!(pair.r > 1.0e-12)
            || !body.catalog.present(body, ArrayId::MopacChargeWelfordMean, source_atom, row))
            return std::nullopt;
        const double q = body.catalog.value(body, ArrayId::MopacChargeWelfordMean, source_atom, row);
        if (!std::isfinite(q)) return std::nullopt;
        const double r3 = pair.r * pair.r * pair.r;
        const double r5 = r3 * pair.r * pair.r;
        Mat3 efg = q * (3.0 * pair.disp * pair.disp.transpose() / r5 - Mat3::Identity() / r3);
        efg *= kCoulombKe;
        efg -= (efg.trace() / 3.0) * Mat3::Identity();
        const model::SphericalTensor st = DecomposeLibrary(efg);
        if (!std::isfinite(st.T0)) return std::nullopt;
        for (double v : st.T2)
            if (!std::isfinite(v)) return std::nullopt;
        return st;
    }

    void foldRelationships(const Body& body,
                           std::size_t atom,
                           std::size_t row,
                           std::uint32_t frame_slot,
                           AtomFrameRecord& record,
                           double sigma) {
        Vec3 ff14sb_field = Vec3::Zero();
        Vec3 mopac_field = Vec3::Zero();
        bool ff14sb_field_present = false;
        bool mopac_field_present = false;

        const std::vector<PairContribution> pairs =
            PerAtomRowPairContributions(body, atom, row, config_, LocalFrame{});
        for (const PairContribution& pair : pairs) {
            if (pair.pointer_flags & PresentFlag) addPartnerSource(record, pair.source_id);
            const RelationshipKeyBase base = relationshipKeyBase(body, atom, row, pair);
            foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("distance_A"), pair.r);
            foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("inv_r3"), pair.inv_r3);
            foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("cos_theta"), pair.cos_theta);
            foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("dipolar"), pair.dipolar);
            foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("kernel_T0"), pair.kernel_T0);
            foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("kernel_T2_norm"),
                                  norm5(pair.kernel_T2));
            for (std::size_t i = 0; i < pair.kernel_T2.size(); ++i) {
                foldRelationshipFacet(base, frame_slot, sigma,
                                      QStringLiteral("kernel_T2_c%1").arg(static_cast<int>(i)),
                                      pair.kernel_T2[i]);
            }
            foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("contribution"),
                                  pair.contribution);

            if (!body.run.protein || pair.source_atom_index < 0
                || static_cast<std::size_t>(pair.source_atom_index) >= body.run.protein->atomCount()
                || !(pair.r > 1.0e-12)) {
                continue;
            }
            const std::size_t source_atom = static_cast<std::size_t>(pair.source_atom_index);
            if (pair.mechanism == QStringLiteral("charge_q_over_r3")
                && body.catalog.present(body, ArrayId::Ff14sbCharge, source_atom, row)) {
                const double q = body.catalog.value(body, ArrayId::Ff14sbCharge, source_atom, row);
                if (std::isfinite(q)) {
                    ff14sb_field += (-q / (pair.r * pair.r * pair.r)) * pair.disp;
                    ff14sb_field_present = true;
                }
            } else if (pair.mechanism == QStringLiteral("field_mopac_coulomb")
                       && body.catalog.present(body, ArrayId::MopacChargeWelfordMean, source_atom, row)) {
                const double q = body.catalog.value(body, ArrayId::MopacChargeWelfordMean, source_atom, row);
                if (std::isfinite(q)) {
                    mopac_field += (-q / (pair.r * pair.r * pair.r)) * pair.disp;
                    mopac_field_present = true;
                }
                if (const std::optional<model::SphericalTensor> efg =
                        mopacPerSourceEfg(body, source_atom, row, pair)) {
                    RelationshipKeyBase efg_base = base;
                    efg_base.mechanism = keyString(QStringLiteral("mopac_charge_q_over_r3"));
                    foldRelationshipTensorFacets(efg_base, frame_slot, sigma, efg->T0, efg->T2);
                }
            }
        }

        if (ff14sb_field_present)
            foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("net_ff14sb_field_norm"), sigma,
                          ff14sb_field.norm());
        if (mopac_field_present)
            foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("net_mopac_field_norm"), sigma,
                          mopac_field.norm());
    }

    void foldSelfProperties(const Body& body,
                            std::size_t atom,
                            std::size_t row,
                            std::uint32_t frame_slot,
                            AtomFrameRecord& record,
                            const DftTarget& target,
                            double sigma) {
        foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("sigma.total_T0"), sigma,
                      target.total_decomp.T0);
        ++sigma_folds_;
        foldT1Components(self_book_, body, atom, row, frame_slot,
                         QStringLiteral("sigma.total_T1"), sigma, target.total_decomp.T1);
        sigma_folds_ += 4;
        foldT2Components(self_book_, body, atom, row, frame_slot,
                         QStringLiteral("sigma.total_T2"), sigma, target.total_decomp.T2);
        sigma_folds_ += 6;
        foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("sigma.total_raw_frobenius"), sigma,
                      tensorFrobenius(target.total_raw));
        foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("sigma.dia_T0"), sigma,
                      target.dia_decomp.T0);
        foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("sigma.para_T0"), sigma,
                      target.para_decomp.T0);
        sigma_folds_ += 3;

        foldDihedrals(body, atom, row, frame_slot, record, sigma);
        foldMopacScalars(body, atom, row, frame_slot, sigma);
        foldEfgFacets(body, atom, row, frame_slot, sigma);
    }

    void foldDihedrals(const Body& body,
                       std::size_t atom,
                       std::size_t row,
                       std::uint32_t frame_slot,
                       AtomFrameRecord& record,
                       double sigma) {
        static constexpr std::array<DihedralKind, 7> kinds = {
            DihedralKind::Phi,   DihedralKind::Psi,  DihedralKind::Omega,
            DihedralKind::Chi1,  DihedralKind::Chi2, DihedralKind::Chi3,
            DihedralKind::Chi4,
        };
        static constexpr std::array<const char*, 7> names = {
            "phi", "psi", "omega", "chi1", "chi2", "chi3", "chi4",
        };
        for (std::size_t i = 0; i < kinds.size(); ++i) {
            const DihedralState st = body.idx.dihedrals.state(kinds[i], atom, row);
            if (!st.present || !std::isfinite(st.radians)) continue;
            const QString prefix = QStringLiteral("dihedral.%1").arg(QString::fromLatin1(names[i]));
            const double sin_value = std::sin(st.radians);
            const double cos_value = std::cos(st.radians);
            record.dihedrals.push_back(DihedralFrameValue{
                names[i],
                st.radians,
                sin_value,
                cos_value,
                st.fixed_bin,
            });
            foldSelfFacet(body, atom, row, frame_slot, prefix + QStringLiteral(".radians"), sigma, st.radians);
            foldSelfFacet(body, atom, row, frame_slot, prefix + QStringLiteral(".sin"), sigma, sin_value);
            foldSelfFacet(body, atom, row, frame_slot, prefix + QStringLiteral(".cos"), sigma, cos_value);
            foldSelfFacet(body, atom, row, frame_slot, prefix + QStringLiteral(".fixed_bin"), sigma,
                          static_cast<double>(st.fixed_bin));
            dihedral_folds_ += 4;
        }
    }

    void foldMopacScalars(const Body& body,
                          std::size_t atom,
                          std::size_t row,
                          std::uint32_t frame_slot,
                          double sigma) {
        if (!body.catalog.present(body, io::FieldKind::MOPACScalars, atom, row, 0)) return;
        const std::optional<double> charge =
            body.catalog.value(body, io::FieldKind::MOPACScalars, atom, row, 0);
        const std::optional<double> s_pop =
            body.catalog.value(body, io::FieldKind::MOPACScalars, atom, row, 1);
        const std::optional<double> p_pop =
            body.catalog.value(body, io::FieldKind::MOPACScalars, atom, row, 2);
        const std::optional<double> valency =
            body.catalog.value(body, io::FieldKind::MOPACScalars, atom, row, 3);
        if (charge) {
            foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("mopac.charge"), sigma, *charge);
            ++mopac_scalar_folds_;
        }
        if (s_pop) {
            foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("mopac.s_population"), sigma, *s_pop);
            ++mopac_scalar_folds_;
        }
        if (p_pop) {
            foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("mopac.p_population"), sigma, *p_pop);
            ++mopac_scalar_folds_;
        }
        if (s_pop && p_pop && std::isfinite(*s_pop) && std::isfinite(*p_pop)
            && (*s_pop + *p_pop) > 1.0e-12) {
            foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("mopac.s_character_fraction"), sigma,
                          *s_pop / (*s_pop + *p_pop));
            ++mopac_scalar_folds_;
        }
        if (valency) {
            foldSelfFacet(body, atom, row, frame_slot, QStringLiteral("mopac.valency"), sigma, *valency);
            ++mopac_scalar_folds_;
        }
    }

    void foldEfgFacets(const Body& body,
                       std::size_t atom,
                       std::size_t row,
                       std::uint32_t frame_slot,
                       double sigma) {
        auto foldT2Array = [&](ArrayId id, const QString& prefix) {
            if (!body.catalog.present(body, id, atom, row)) return;
            const std::array<double, 5> t2 = body.catalog.valueT2(body, id, atom, row);
            foldT2Components(self_book_, body, atom, row, frame_slot, prefix, sigma, t2);
            efg_folds_ += 6;
        };
        foldT2Array(ArrayId::ApbsEfg, QStringLiteral("efg.apbs_t2"));
        foldT2Array(ArrayId::Aimnet2Efg, QStringLiteral("efg.aimnet2_t2"));
        foldT2Array(ArrayId::MopacCoulombShielding,
                    QStringLiteral("efg.mopac_coulomb_derived_shielding_t2"));
        foldT2Array(ArrayId::MopacMcShielding, QStringLiteral("efg.mopac_mc_shielding_t2"));
    }

    bool finiteTensor(const model::SphericalTensor& tensor) const {
        if (!std::isfinite(tensor.T0)) return false;
        for (double v : tensor.T1)
            if (!std::isfinite(v)) return false;
        for (double v : tensor.T2)
            if (!std::isfinite(v)) return false;
        return true;
    }

    bool nonzeroLarsenTensor(const model::SphericalTensor& tensor) const {
        return std::abs(tensor.T0) > 1.0e-12 || norm3(tensor.T1) > 1.0e-12
               || norm5(tensor.T2) > 1.0e-12;
    }

    void foldLarsenTensor(const Body& body,
                          std::size_t atom,
                          std::size_t row,
                          std::uint32_t frame_slot,
                          AtomFrameRecord& record,
                          double sigma,
                          const model::QtShieldingTimeSeries* ts,
                          const QString& source_kind,
                          int32_t source_id) {
        if (!ts || atom >= ts->n_atoms || row >= ts->n_frames || !ts->sourceAttachedAt(row))
            return;
        const model::SphericalTensor tensor = ts->at(atom, row);
        if (!finiteTensor(tensor) || !nonzeroLarsenTensor(tensor)) return;
        const RelationshipKeyBase base =
            larsenRelationshipKeyBase(body, atom, row, source_kind, source_id);
        foldRelationshipTensorFacets(base, frame_slot, sigma, tensor.T0, tensor.T2);
        addPartnerSource(record, source_id);
    }

    void foldLarsenHbondRelationships(const Body& body,
                                      std::size_t atom,
                                      std::size_t row,
                                      std::uint32_t frame_slot,
                                      AtomFrameRecord& record,
                                      double sigma) {
        const io::QtTrajectoryH5* h5 = body.run.h5();
        if (!h5) return;
        foldLarsenTensor(body, atom, row, frame_slot, record, sigma,
                         h5->larsenHBond1pHBShielding(),
                         QStringLiteral("larsen_hbond_1pHB"), 0);
        foldLarsenTensor(body, atom, row, frame_slot, record, sigma,
                         h5->larsenHBond2pHBShielding(),
                         QStringLiteral("larsen_hbond_2pHB"), 1);
        foldLarsenTensor(body, atom, row, frame_slot, record, sigma,
                         h5->larsenHBond1pHaBShielding(),
                         QStringLiteral("larsen_hbond_1pHaB"), 2);
        foldLarsenTensor(body, atom, row, frame_slot, record, sigma,
                         h5->larsenHBond2pHaBShielding(),
                         QStringLiteral("larsen_hbond_2pHaB"), 3);

        const model::QtScalarTimeSeries* count = h5->larsenHBondCount();
        if (count && atom < count->n_atoms && row < count->n_frames && count->sourceAttachedAt(row)) {
            const double value = count->at(atom, row);
            if (std::isfinite(value) && value > 0.0) {
                const RelationshipKeyBase base =
                    larsenRelationshipKeyBase(body, atom, row,
                                              QStringLiteral("larsen_hbond_count"), 4);
                foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("count"), value);
                addPartnerSource(record, 4);
            }
        }

        const model::QtScalarTimeSeries* water = h5->larsenHBondWaterTerm();
        if (water && atom < water->n_atoms && row < water->n_frames && water->sourceAttachedAt(row)) {
            const double value = water->at(atom, row);
            if (std::isfinite(value) && std::abs(value) > 1.0e-12) {
                const RelationshipKeyBase base =
                    larsenRelationshipKeyBase(body, atom, row,
                                              QStringLiteral("larsen_hbond_water_term"), 5);
                foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("kernel_T0"), value);
                foldRelationshipFacet(base, frame_slot, sigma, QStringLiteral("contribution"),
                                      std::abs(value));
                addPartnerSource(record, 5);
            }
        }
    }

    PerAtomSubstrateConfig config_;
    RelationshipBook relationship_book_;
    SelfBook self_book_;
    std::vector<double> sigma_by_frame_;
    std::vector<AtomFrameRecord> frame_records_;
    std::size_t frame_events_ = 0;
    std::size_t dft_present_ = 0;
    std::size_t sigma_folds_ = 0;
    std::size_t dihedral_folds_ = 0;
    std::size_t mopac_scalar_folds_ = 0;
    std::size_t efg_folds_ = 0;
};

NeighbourhoodAccumulatingAnalysisAtom::NeighbourhoodAccumulatingAnalysisAtom(
    const Body& body, std::size_t atom_index, PerAtomSubstrateConfig config)
    : AnalysisAtom(body, atom_index), impl_(std::make_unique<Impl>(std::move(config))) {}

void NeighbourhoodAccumulatingAnalysisAtom::observeFrame(std::size_t h5_row) {
    impl_->observe(body_, atom_index_, h5_row);
}

AnalysisAtomTruth NeighbourhoodAccumulatingAnalysisAtom::diagnostics() const {
    return impl_->diagnostics();
}

AnalysisAtomDiagnostics RunAnalysisAtomFirstPass(const Body& body,
                                                 const PerAtomSubstrateConfig& config) {
    AnalysisAtomDiagnostics diagnostics;
    if (!body.run.protein) return diagnostics;

    const std::size_t atom_count = body.run.protein->atomCount();
    const std::vector<std::size_t>& rows = body.run.frameMap.dftRows();
    diagnostics.atom_count = atom_count;
    diagnostics.dft_rows = rows.size();

    std::vector<std::unique_ptr<AnalysisAtom>> atoms;
    atoms.reserve(atom_count);
    for (std::size_t atom = 0; atom < atom_count; ++atom) {
        atoms.push_back(std::make_unique<NeighbourhoodAccumulatingAnalysisAtom>(body, atom, config));
    }

    for (std::size_t row : rows) {
        for (const std::unique_ptr<AnalysisAtom>& atom : atoms) atom->observeFrame(row);
    }

    for (const std::unique_ptr<AnalysisAtom>& atom : atoms) {
        const AnalysisAtomTruth truth = atom->diagnostics();
        diagnostics.frame_events += truth.frame_events;
        diagnostics.dft_present += truth.dft_present;
        diagnostics.relationship_folds += truth.relationship_folds;
        diagnostics.self_property_folds += truth.self_property_folds;
        diagnostics.relationship_organs += truth.relationship_organs;
        diagnostics.self_property_organs += truth.self_property_organs;
        diagnostics.sigma_folds += truth.sigma_folds;
        diagnostics.dihedral_folds += truth.dihedral_folds;
        diagnostics.mopac_scalar_folds += truth.mopac_scalar_folds;
        diagnostics.efg_folds += truth.efg_folds;
        if (!diagnostics.has_sample && truth.has_sample) {
            diagnostics.sample_label = truth.sample_label;
            diagnostics.sample_value = truth.sample_value;
            diagnostics.has_sample = true;
        }
    }
    return diagnostics;
}

}  // namespace h5reader::rediscover
