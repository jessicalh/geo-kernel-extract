#include "AnalysisAtom.h"

#include "ExtractionSupport.h"
#include "Verbs.h"

#include "../io/QtFieldCatalog.gen.h"
#include "../model/QtAtom.h"
#include "../model/QtBond.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtRing.h"
#include "../model/QtTopology.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>
#include <boost/accumulators/statistics/rolling_variance.hpp>
#include <boost/accumulators/statistics/rolling_window.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <QByteArray>
#include <QStringList>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <set>
#include <unordered_map>
#include <utility>

namespace h5reader::rediscover {
namespace {

namespace ba = boost::accumulators;

constexpr std::size_t kRollingWindow = 32;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

using ScalarAccumulator =
    ba::accumulator_set<double,
                        ba::stats<ba::tag::count,
                                  ba::tag::mean,
                                  ba::tag::variance,
                                  ba::tag::min,
                                  ba::tag::max,
                                  ba::tag::rolling_mean,
                                  ba::tag::rolling_variance,
                                  ba::tag::extended_p_square>>;

using CovarianceAccumulator =
    ba::accumulator_set<double,
                        ba::stats<ba::tag::covariance<double, ba::tag::covariate1>>>;

const std::vector<double>& coverageProbabilities() {
    static const std::vector<double> probs{0.05, 0.25, 0.50, 0.75, 0.95};
    return probs;
}

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

QString relationshipCategory(const Body& body,
                             std::size_t atom,
                             std::size_t row,
                             const PairContribution& pair,
                             const QString& facet) {
    return qjoin({QStringLiteral("relationship"),
                  targetChemistryCategory(body, atom, row),
                  contactedResidueLabel(body, pair),
                  relationshipScope(body, atom, pair),
                  QStringLiteral("mechanism=%1").arg(pair.mechanism),
                  QStringLiteral("source_kind=%1").arg(pair.source_kind),
                  QStringLiteral("source_id=%1").arg(pair.source_id),
                  QStringLiteral("source_category=%1").arg(pair.source_category_ord),
                  QStringLiteral("facet=%1").arg(facet)});
}

QString selfCategory(const Body& body, std::size_t atom, std::size_t row, const QString& property) {
    return qjoin({QStringLiteral("self_property"),
                  targetChemistryCategory(body, atom, row),
                  QStringLiteral("residue=self"),
                  QStringLiteral("scope=self"),
                  QStringLiteral("property=%1").arg(property)});
}

struct ScalarOrgan {
    ScalarOrgan()
        : value(ba::tag::rolling_window::window_size = kRollingWindow,
                ba::tag::extended_p_square::probabilities = coverageProbabilities()),
          sigma(ba::tag::rolling_window::window_size = kRollingWindow,
                ba::tag::extended_p_square::probabilities = coverageProbabilities()) {}

    void fold(double sigma_iso, double x) {
        if (!std::isfinite(sigma_iso) || !std::isfinite(x)) return;
        value(x);
        sigma(sigma_iso);
        covariance(sigma_iso, ba::covariate1 = x);
        ++n;

        const double dx = x - mean_x;
        mean_x += dx / static_cast<double>(n);
        const double dy = sigma_iso - mean_sigma;
        mean_sigma += dy / static_cast<double>(n);
        c_xy += dx * (sigma_iso - mean_sigma);
        m2_x += dx * (x - mean_x);

        rolling.emplace_back(sigma_iso, x);
        rolling_sum_sigma += sigma_iso;
        rolling_sum_x += x;
        rolling_sum_sigmax += sigma_iso * x;
        rolling_sum_x2 += x * x;
        if (rolling.size() > kRollingWindow) {
            const auto old = rolling.front();
            rolling.pop_front();
            rolling_sum_sigma -= old.first;
            rolling_sum_x -= old.second;
            rolling_sum_sigmax -= old.first * old.second;
            rolling_sum_x2 -= old.second * old.second;
        }
    }

    std::size_t count() const { return n; }

    double slope() const {
        return (n > 1 && m2_x > 0.0) ? c_xy / m2_x : kNaN;
    }

    double rollingSlope() const {
        const std::size_t rn = rolling.size();
        if (rn < 2) return kNaN;
        const double inv_n = 1.0 / static_cast<double>(rn);
        const double cov = rolling_sum_sigmax - rolling_sum_sigma * rolling_sum_x * inv_n;
        const double var_x = rolling_sum_x2 - rolling_sum_x * rolling_sum_x * inv_n;
        return var_x > 0.0 ? cov / var_x : kNaN;
    }

    ScalarAccumulator value;
    ScalarAccumulator sigma;
    CovarianceAccumulator covariance;
    std::size_t n = 0;
    double mean_x = 0.0;
    double mean_sigma = 0.0;
    double c_xy = 0.0;
    double m2_x = 0.0;
    std::deque<std::pair<double, double>> rolling;
    double rolling_sum_sigma = 0.0;
    double rolling_sum_x = 0.0;
    double rolling_sum_sigmax = 0.0;
    double rolling_sum_x2 = 0.0;
};

struct FoldBook {
    void fold(const QString& key, double sigma_iso, double value) {
        if (!finite(sigma_iso) || !finite(value)) return;
        const std::string skey = keyString(key);
        organs[skey].fold(sigma_iso, value);
        ++fold_count;
        if (!sample) sample = std::make_pair(key, value);
    }

    std::size_t organCount() const { return organs.size(); }

    std::size_t foldedObservations() const {
        std::size_t out = 0;
        for (const auto& item : organs) out += item.second.count();
        return out;
    }

    std::unordered_map<std::string, ScalarOrgan> organs;
    std::size_t fold_count = 0;
    std::optional<std::pair<QString, double>> sample;
};

void foldT2Components(FoldBook& book,
                      const QString& base,
                      double sigma,
                      const std::array<double, 5>& t2) {
    book.fold(base + QStringLiteral(".norm"), sigma, norm5(t2));
    for (std::size_t i = 0; i < t2.size(); ++i)
        book.fold(base + QStringLiteral(".c%1").arg(static_cast<int>(i)), sigma, t2[i]);
}

void foldT1Components(FoldBook& book,
                      const QString& base,
                      double sigma,
                      const std::array<double, 3>& t1) {
    book.fold(base + QStringLiteral(".norm"), sigma, norm3(t1));
    for (std::size_t i = 0; i < t1.size(); ++i)
        book.fold(base + QStringLiteral(".c%1").arg(static_cast<int>(i)), sigma, t1[i]);
}

}  // namespace

AnalysisAtom::AnalysisAtom(const Body& body, std::size_t atom_index)
    : body_(body), atom_index_(atom_index) {}

class NeighbourhoodAccumulatingAnalysisAtom::Impl {
public:
    explicit Impl(PerAtomSubstrateConfig config) : config_(std::move(config)) {}

    void observe(const Body& body, std::size_t atom, std::size_t h5_row) {
        ++frame_events_;
        if (h5_row >= body.run.frameMap.frameCount()) return;
        const std::size_t original = body.run.frameMap.originalIndex(h5_row);
        const DftTarget target = BuildTarget(body.run, atom, original, LocalFrame{});
        if (!target.present || !std::isfinite(target.total_decomp.T0)) return;
        ++dft_present_;
        const double sigma_iso = target.total_decomp.T0;

        foldSelfProperties(body, atom, h5_row, target, sigma_iso);
        foldRelationships(body, atom, h5_row, sigma_iso);
    }

    void contribute(AnalysisAtomDiagnostics& out) const {
        out.frame_events += frame_events_;
        out.dft_present += dft_present_;
        out.relationship_folds += relationship_book_.foldedObservations();
        out.self_property_folds += self_book_.foldedObservations();
        out.relationship_organs += relationship_book_.organCount();
        out.self_property_organs += self_book_.organCount();
        out.sigma_folds += sigma_folds_;
        out.dihedral_folds += dihedral_folds_;
        out.mopac_scalar_folds += mopac_scalar_folds_;
        out.efg_folds += efg_folds_;
        if (!out.has_sample) {
            if (relationship_book_.sample) {
                out.sample_label = relationship_book_.sample->first;
                out.sample_value = relationship_book_.sample->second;
                out.has_sample = true;
            } else if (self_book_.sample) {
                out.sample_label = self_book_.sample->first;
                out.sample_value = self_book_.sample->second;
                out.has_sample = true;
            }
        }
    }

private:
    void foldRelationshipFacet(const Body& body,
                               std::size_t atom,
                               std::size_t row,
                               const PairContribution& pair,
                               const QString& facet,
                               double sigma,
                               double value) {
        relationship_book_.fold(relationshipCategory(body, atom, row, pair, facet), sigma, value);
    }

    void foldSelfFacet(const Body& body,
                       std::size_t atom,
                       std::size_t row,
                       const QString& property,
                       double sigma,
                       double value) {
        self_book_.fold(selfCategory(body, atom, row, property), sigma, value);
    }

    void foldRelationships(const Body& body, std::size_t atom, std::size_t row, double sigma) {
        Vec3 ff14sb_field = Vec3::Zero();
        Vec3 mopac_field = Vec3::Zero();
        bool ff14sb_field_present = false;
        bool mopac_field_present = false;

        const std::vector<PairContribution> pairs =
            PerAtomRowPairContributions(body, atom, row, config_, LocalFrame{});
        for (const PairContribution& pair : pairs) {
            foldRelationshipFacet(body, atom, row, pair, QStringLiteral("distance_A"), sigma, pair.r);
            foldRelationshipFacet(body, atom, row, pair, QStringLiteral("inv_r3"), sigma, pair.inv_r3);
            foldRelationshipFacet(body, atom, row, pair, QStringLiteral("cos_theta"), sigma, pair.cos_theta);
            foldRelationshipFacet(body, atom, row, pair, QStringLiteral("dipolar"), sigma, pair.dipolar);
            foldRelationshipFacet(body, atom, row, pair, QStringLiteral("kernel_T0"), sigma, pair.kernel_T0);
            foldRelationshipFacet(body, atom, row, pair, QStringLiteral("kernel_T2_norm"), sigma,
                                  norm5(pair.kernel_T2));
            for (std::size_t i = 0; i < pair.kernel_T2.size(); ++i) {
                foldRelationshipFacet(body, atom, row, pair,
                                      QStringLiteral("kernel_T2_c%1").arg(static_cast<int>(i)),
                                      sigma, pair.kernel_T2[i]);
            }
            foldRelationshipFacet(body, atom, row, pair, QStringLiteral("contribution"), sigma,
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
            }
        }

        if (ff14sb_field_present)
            foldSelfFacet(body, atom, row, QStringLiteral("net_ff14sb_field_norm"), sigma,
                          ff14sb_field.norm());
        if (mopac_field_present)
            foldSelfFacet(body, atom, row, QStringLiteral("net_mopac_field_norm"), sigma,
                          mopac_field.norm());
    }

    void foldSelfProperties(const Body& body,
                            std::size_t atom,
                            std::size_t row,
                            const DftTarget& target,
                            double sigma) {
        foldSelfFacet(body, atom, row, QStringLiteral("sigma.total_T0"), sigma,
                      target.total_decomp.T0);
        ++sigma_folds_;
        foldT1Components(self_book_, selfCategory(body, atom, row, QStringLiteral("sigma.total_T1")),
                         sigma, target.total_decomp.T1);
        sigma_folds_ += 4;
        foldT2Components(self_book_, selfCategory(body, atom, row, QStringLiteral("sigma.total_T2")),
                         sigma, target.total_decomp.T2);
        sigma_folds_ += 6;
        foldSelfFacet(body, atom, row, QStringLiteral("sigma.total_raw_frobenius"), sigma,
                      tensorFrobenius(target.total_raw));
        foldSelfFacet(body, atom, row, QStringLiteral("sigma.dia_T0"), sigma,
                      target.dia_decomp.T0);
        foldSelfFacet(body, atom, row, QStringLiteral("sigma.para_T0"), sigma,
                      target.para_decomp.T0);
        sigma_folds_ += 3;

        foldDihedrals(body, atom, row, sigma);
        foldMopacScalars(body, atom, row, sigma);
        foldEfgFacets(body, atom, row, sigma);
    }

    void foldDihedrals(const Body& body, std::size_t atom, std::size_t row, double sigma) {
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
            foldSelfFacet(body, atom, row, prefix + QStringLiteral(".radians"), sigma, st.radians);
            foldSelfFacet(body, atom, row, prefix + QStringLiteral(".sin"), sigma, std::sin(st.radians));
            foldSelfFacet(body, atom, row, prefix + QStringLiteral(".cos"), sigma, std::cos(st.radians));
            foldSelfFacet(body, atom, row, prefix + QStringLiteral(".fixed_bin"), sigma,
                          static_cast<double>(st.fixed_bin));
            dihedral_folds_ += 4;
        }
    }

    void foldMopacScalars(const Body& body, std::size_t atom, std::size_t row, double sigma) {
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
            foldSelfFacet(body, atom, row, QStringLiteral("mopac.charge"), sigma, *charge);
            ++mopac_scalar_folds_;
        }
        if (s_pop) {
            foldSelfFacet(body, atom, row, QStringLiteral("mopac.s_population"), sigma, *s_pop);
            ++mopac_scalar_folds_;
        }
        if (p_pop) {
            foldSelfFacet(body, atom, row, QStringLiteral("mopac.p_population"), sigma, *p_pop);
            ++mopac_scalar_folds_;
        }
        if (s_pop && p_pop && std::isfinite(*s_pop) && std::isfinite(*p_pop)
            && (*s_pop + *p_pop) > 1.0e-12) {
            foldSelfFacet(body, atom, row, QStringLiteral("mopac.s_character_fraction"), sigma,
                          *s_pop / (*s_pop + *p_pop));
            ++mopac_scalar_folds_;
        }
        if (valency) {
            foldSelfFacet(body, atom, row, QStringLiteral("mopac.valency"), sigma, *valency);
            ++mopac_scalar_folds_;
        }
    }

    void foldEfgFacets(const Body& body, std::size_t atom, std::size_t row, double sigma) {
        auto foldT2Array = [&](ArrayId id, const QString& prefix) {
            if (!body.catalog.present(body, id, atom, row)) return;
            const std::array<double, 5> t2 = body.catalog.valueT2(body, id, atom, row);
            foldT2Components(self_book_, selfCategory(body, atom, row, prefix), sigma, t2);
            efg_folds_ += 6;
        };
        foldT2Array(ArrayId::ApbsEfg, QStringLiteral("efg.apbs_t2"));
        foldT2Array(ArrayId::Aimnet2Efg, QStringLiteral("efg.aimnet2_t2"));
        foldT2Array(ArrayId::MopacCoulombShielding,
                    QStringLiteral("efg.mopac_coulomb_derived_shielding_t2"));
        foldT2Array(ArrayId::MopacMcShielding, QStringLiteral("efg.mopac_mc_shielding_t2"));
    }

    PerAtomSubstrateConfig config_;
    FoldBook relationship_book_;
    FoldBook self_book_;
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

void NeighbourhoodAccumulatingAnalysisAtom::contributeDiagnostics(
    AnalysisAtomDiagnostics& out) const {
    impl_->contribute(out);
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

    for (const std::unique_ptr<AnalysisAtom>& atom : atoms) atom->contributeDiagnostics(diagnostics);
    return diagnostics;
}

}  // namespace h5reader::rediscover
