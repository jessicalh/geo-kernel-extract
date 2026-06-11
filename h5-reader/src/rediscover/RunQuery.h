// RunQuery -- foundation query layer over the resident producer spine.
//
// This is reader-side query machinery. It reads producer FieldKind inputs
// through Catalog::present/value and emits query rows to caller-owned sinks;
// it does not add reader outputs to the producer catalog.

#pragma once

#include "AnalysisBody.h"

#include "../io/QtFieldCatalog.gen.h"

#include <QString>

#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace h5reader::rediscover {

enum class TraversalDomain {
    AllFrames,
    DftRows,
    TargetPresentRows,
};

enum class GatherRank {
    Scalar,
    Vec3,
    T2_5,
    Tensor9,
    Dynamic,
};

struct FieldRef {
    io::FieldKind kind = io::FieldKind::Pos;
    QString name;
    GatherRank rank = GatherRank::Dynamic;
    int component = -1;  // -1 gathers all components implied by rank/field.

    static FieldRef Producer(io::FieldKind kind, QString name = {}, int component = -1);
};

struct GatheredField {
    FieldRef ref;
    bool present = false;
    QString absence_reason;
    std::vector<double> values;
};

GatherRank RankForField(const io::FieldSpec& spec);
std::size_t ExpectedComponentCount(GatherRank rank, std::size_t catalogComponentCount);
GatheredField GatherField(const Body& body, const FieldRef& ref,
                          std::size_t nativeRow, std::size_t frame);

class Selector {
public:
    using AtomPredicate = std::function<bool(const Body&, std::size_t atom)>;
    using FramePredicate = std::function<bool(const Body&, std::size_t atom, std::size_t frame)>;
    using LabelFn = std::function<QString(const Body&, std::size_t atom, std::size_t frame)>;

    Selector() = default;

    static Selector All(QString label = QStringLiteral("all"));
    static Selector Atom(QString name, AtomPredicate predicate, LabelFn label = {});
    static Selector Frame(QString name, FramePredicate predicate, LabelFn label = {});
    static Selector TwoPhase(QString name, AtomPredicate atomPredicate,
                             FramePredicate framePredicate, LabelFn label = {});

    static Selector And(std::vector<Selector> selectors, QString name = {});
    static Selector Or(std::vector<Selector> selectors, QString name = {});

    const QString& name() const { return name_; }
    bool hasAtomPrefilter() const { return static_cast<bool>(atomPredicate_); }
    bool hasFrameFilter() const { return static_cast<bool>(framePredicate_); }
    bool atomAccepts(const Body& body, std::size_t atom) const;
    bool frameAccepts(const Body& body, std::size_t atom, std::size_t frame) const;
    QString label(const Body& body, std::size_t atom, std::size_t frame) const;

private:
    QString name_;
    AtomPredicate atomPredicate_;
    FramePredicate framePredicate_;
    LabelFn label_;
};

std::vector<std::size_t> ApplyAtomPrefilters(const Body& body,
                                             const std::vector<std::size_t>& scope,
                                             const std::vector<Selector>& selectors);
bool AcceptsFrameSelectors(const Body& body,
                           std::size_t atom,
                           std::size_t frame,
                           const std::vector<Selector>& selectors);
QStringList SelectorLabels(const Body& body,
                           std::size_t atom,
                           std::size_t frame,
                           const std::vector<Selector>& selectors);

struct Query {
    TraversalDomain domain = TraversalDomain::DftRows;
    std::vector<Selector> where;
    std::vector<FieldRef> gather;
};

struct QueryRow {
    std::size_t atom = 0;
    std::size_t frame = 0;
    std::size_t original_index = 0;
    QStringList labels;
    std::vector<GatheredField> fields;
};

using QueryRowSink = std::function<void(const QueryRow&)>;

std::vector<std::size_t> FramesForTraversal(const Body& body, TraversalDomain domain);
std::size_t RunQuery(const Body& body, const Query& query, const QueryRowSink& sink);

Selector TargetPresentSelector(io::FieldKind target = io::FieldKind::OrcaTotal);
Selector IupacNameSelector(QString iupacName);
Selector ChemicalCategorySelector(QString category);
Selector SecondaryStructureSelector(SecondaryStructure3 ss3);
Selector SecondaryStructureSelector(SecondaryStructure8 ss8);
Selector DihedralBinSelector(DihedralKind kind, int fixedBin);
Selector DihedralRangeSelector(DihedralKind kind, double loRadians, double hiRadians);

}  // namespace h5reader::rediscover
