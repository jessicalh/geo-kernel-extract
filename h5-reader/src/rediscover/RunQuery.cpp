#include "RunQuery.h"

#include "../model/QtProtein.h"

#include <algorithm>

namespace h5reader::rediscover {
namespace {

QString qsv(std::string_view s) {
    return QString::fromUtf8(s.data(), static_cast<qsizetype>(s.size()));
}

std::size_t fieldComponentCount(const Body& body, io::FieldKind kind, std::size_t frame) {
    return body.catalog.componentCount(body, kind, frame);
}

}  // namespace

FieldRef FieldRef::Producer(io::FieldKind kind, QString name, int component) {
    const io::FieldSpec& spec = io::FieldSpecFor(kind);
    FieldRef ref;
    ref.kind = kind;
    ref.name = name.isEmpty() ? qsv(spec.stem) : std::move(name);
    ref.rank = RankForField(spec);
    ref.component = component;
    return ref;
}

GatherRank RankForField(const io::FieldSpec& spec) {
    if (spec.cols == 1 || spec.cols < 0) return GatherRank::Scalar;
    if (spec.cols == 3) return GatherRank::Vec3;
    if (spec.cols == 5) return GatherRank::T2_5;
    if (spec.cols == 9) return GatherRank::Tensor9;
    return GatherRank::Dynamic;
}

std::size_t ExpectedComponentCount(GatherRank rank, std::size_t catalogComponentCount) {
    switch (rank) {
    case GatherRank::Scalar: return 1;
    case GatherRank::Vec3: return 3;
    case GatherRank::T2_5: return 5;
    case GatherRank::Tensor9: return 9;
    case GatherRank::Dynamic: return catalogComponentCount == 0 ? 1 : catalogComponentCount;
    }
    return catalogComponentCount == 0 ? 1 : catalogComponentCount;
}

GatheredField GatherField(const Body& body, const FieldRef& ref,
                          std::size_t nativeRow, std::size_t frame) {
    GatheredField out;
    out.ref = ref;
    const std::size_t catalogComponents = fieldComponentCount(body, ref.kind, frame);
    const std::size_t expected = ExpectedComponentCount(ref.rank, catalogComponents);
    if (ref.component >= 0) {
        QString reason;
        const std::optional<double> v =
            body.catalog.value(body, ref.kind, nativeRow, frame, ref.component, &reason);
        if (!v) {
            out.absence_reason = reason;
            return out;
        }
        out.present = true;
        out.values.push_back(*v);
        return out;
    }

    out.values.reserve(expected);
    for (std::size_t c = 0; c < expected; ++c) {
        QString reason;
        const std::optional<double> v =
            body.catalog.value(body, ref.kind, nativeRow, frame, static_cast<int>(c), &reason);
        if (!v) {
            out.values.clear();
            out.absence_reason = reason;
            return out;
        }
        out.values.push_back(*v);
    }
    out.present = true;
    return out;
}

Selector Selector::All(QString label) {
    Selector s;
    s.name_ = std::move(label);
    s.label_ = [name = s.name_](const Body&, std::size_t, std::size_t) { return name; };
    return s;
}

Selector Selector::Atom(QString name, AtomPredicate predicate, LabelFn label) {
    Selector s;
    s.name_ = std::move(name);
    s.atomPredicate_ = std::move(predicate);
    s.label_ = std::move(label);
    return s;
}

Selector Selector::Frame(QString name, FramePredicate predicate, LabelFn label) {
    Selector s;
    s.name_ = std::move(name);
    s.framePredicate_ = std::move(predicate);
    s.label_ = std::move(label);
    return s;
}

Selector Selector::TwoPhase(QString name, AtomPredicate atomPredicate,
                            FramePredicate framePredicate, LabelFn label) {
    Selector s;
    s.name_ = std::move(name);
    s.atomPredicate_ = std::move(atomPredicate);
    s.framePredicate_ = std::move(framePredicate);
    s.label_ = std::move(label);
    return s;
}

Selector Selector::And(std::vector<Selector> selectors, QString name) {
    Selector s;
    s.name_ = name.isEmpty() ? QStringLiteral("and") : std::move(name);
    s.atomPredicate_ = [selectors](const Body& body, std::size_t atom) {
        for (const Selector& selector : selectors)
            if (!selector.atomAccepts(body, atom)) return false;
        return true;
    };
    s.framePredicate_ = [selectors](const Body& body, std::size_t atom, std::size_t frame) {
        for (const Selector& selector : selectors)
            if (!selector.frameAccepts(body, atom, frame)) return false;
        return true;
    };
    s.label_ = [selectors](const Body& body, std::size_t atom, std::size_t frame) {
        QStringList labels;
        for (const Selector& selector : selectors) {
            const QString label = selector.label(body, atom, frame);
            if (!label.isEmpty()) labels << label;
        }
        return labels.join(QLatin1Char('&'));
    };
    return s;
}

Selector Selector::Or(std::vector<Selector> selectors, QString name) {
    Selector s;
    s.name_ = name.isEmpty() ? QStringLiteral("or") : std::move(name);
    s.atomPredicate_ = [selectors](const Body& body, std::size_t atom) {
        for (const Selector& selector : selectors)
            if (selector.atomAccepts(body, atom)) return true;
        return false;
    };
    s.framePredicate_ = [selectors](const Body& body, std::size_t atom, std::size_t frame) {
        for (const Selector& selector : selectors)
            if (selector.frameAccepts(body, atom, frame)) return true;
        return false;
    };
    s.label_ = [selectors](const Body& body, std::size_t atom, std::size_t frame) {
        QStringList labels;
        for (const Selector& selector : selectors) {
            if (!selector.atomAccepts(body, atom) || !selector.frameAccepts(body, atom, frame))
                continue;
            const QString label = selector.label(body, atom, frame);
            if (!label.isEmpty()) labels << label;
        }
        return labels.join(QLatin1Char('|'));
    };
    return s;
}

bool Selector::atomAccepts(const Body& body, std::size_t atom) const {
    return atomPredicate_ ? atomPredicate_(body, atom) : true;
}

bool Selector::frameAccepts(const Body& body, std::size_t atom, std::size_t frame) const {
    return framePredicate_ ? framePredicate_(body, atom, frame) : true;
}

QString Selector::label(const Body& body, std::size_t atom, std::size_t frame) const {
    if (label_) return label_(body, atom, frame);
    return name_;
}

std::vector<std::size_t> ApplyAtomPrefilters(const Body& body,
                                             const std::vector<std::size_t>& scope,
                                             const std::vector<Selector>& selectors) {
    std::vector<std::size_t> out;
    out.reserve(scope.size());
    for (std::size_t atom : scope) {
        bool ok = true;
        for (const Selector& selector : selectors) {
            if (!selector.atomAccepts(body, atom)) {
                ok = false;
                break;
            }
        }
        if (ok) out.push_back(atom);
    }
    return out;
}

bool AcceptsFrameSelectors(const Body& body,
                           std::size_t atom,
                           std::size_t frame,
                           const std::vector<Selector>& selectors) {
    for (const Selector& selector : selectors) {
        if (!selector.frameAccepts(body, atom, frame)) return false;
    }
    return true;
}

QStringList SelectorLabels(const Body& body,
                           std::size_t atom,
                           std::size_t frame,
                           const std::vector<Selector>& selectors) {
    QStringList out;
    for (const Selector& selector : selectors) {
        const QString label = selector.label(body, atom, frame);
        if (!label.isEmpty()) out << label;
    }
    return out;
}

std::vector<std::size_t> FramesForTraversal(const Body& body, TraversalDomain domain) {
    switch (domain) {
    case TraversalDomain::AllFrames: {
        std::vector<std::size_t> frames(body.run.frameMap.frameCount());
        for (std::size_t i = 0; i < frames.size(); ++i) frames[i] = i;
        return frames;
    }
    case TraversalDomain::DftRows:
        return body.run.frameMap.dftRows();
    case TraversalDomain::TargetPresentRows: {
        std::vector<std::size_t> frames;
        for (std::size_t row = 0; row < body.run.frameMap.frameCount(); ++row) {
            const std::size_t atomCount = body.run.protein ? body.run.protein->atomCount() : 0;
            for (std::size_t atom = 0; atom < atomCount; ++atom) {
                QString reason;
                if (body.catalog.present(body, io::FieldKind::OrcaTotal, atom, row, 0, &reason)) {
                    frames.push_back(row);
                    break;
                }
            }
        }
        return frames;
    }
    }
    return {};
}

std::size_t RunQuery(const Body& body, const Query& query, const QueryRowSink& sink) {
    if (!body.run.protein) return 0;
    std::vector<std::size_t> allAtoms(body.run.protein->atomCount());
    for (std::size_t atom = 0; atom < allAtoms.size(); ++atom) allAtoms[atom] = atom;
    const std::vector<std::size_t> atoms = ApplyAtomPrefilters(body, allAtoms, query.where);
    const std::vector<std::size_t> frames = FramesForTraversal(body, query.domain);

    std::size_t rows = 0;
    for (std::size_t frame : frames) {
        const std::size_t original = body.run.frameMap.originalIndex(frame);
        for (std::size_t atom : atoms) {
            if (!AcceptsFrameSelectors(body, atom, frame, query.where)) continue;
            QueryRow row;
            row.atom = atom;
            row.frame = frame;
            row.original_index = original;
            row.labels = SelectorLabels(body, atom, frame, query.where);
            row.fields.reserve(query.gather.size());
            for (const FieldRef& ref : query.gather)
                row.fields.push_back(GatherField(body, ref, atom, frame));
            sink(row);
            ++rows;
        }
    }
    return rows;
}

Selector TargetPresentSelector(io::FieldKind target) {
    return Selector::Frame(
        QStringLiteral("target_present"),
        [target](const Body& body, std::size_t atom, std::size_t frame) {
            QString reason;
            return body.catalog.present(body, target, atom, frame, 0, &reason);
        },
        [](const Body&, std::size_t, std::size_t) { return QStringLiteral("target_present"); });
}

}  // namespace h5reader::rediscover
