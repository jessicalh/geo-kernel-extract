#include "TypedAtomIndex.h"

#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"

#include <QStringList>

namespace h5reader::rediscover {

TypedAtomIndex::TypedAtomIndex(const model::QtProtein& protein) : protein_(&protein) {}

std::vector<int32_t> TypedAtomIndex::residueScope(std::size_t residueIndex) const {
    if (!protein_ || residueIndex >= protein_->residueCount()) return {};
    return protein_->residue(residueIndex).atomIndices;
}

std::vector<int32_t> TypedAtomIndex::select(const std::vector<int32_t>& scope,
                                            const TypedAtomSelector& selector) const {
    std::vector<int32_t> out;
    if (!protein_) return out;
    for (int32_t ai : scope) {
        if (ai < 0 || static_cast<std::size_t>(ai) >= protein_->atomCount()) continue;
        const model::QtAtom& a = protein_->atom(static_cast<std::size_t>(ai));
        if (selector.element && a.element != *selector.element) continue;
        if (selector.locant && a.locant != *selector.locant) continue;
        if (selector.branch && !(a.branch == *selector.branch)) continue;
        if (selector.diIndex && a.diIndex != *selector.diIndex) continue;
        if (selector.backboneRole && a.backboneRole != *selector.backboneRole) continue;
        out.push_back(ai);
    }
    return out;
}

std::optional<int32_t> TypedAtomIndex::selectUnique(const std::vector<int32_t>& scope,
                                                    const TypedAtomSelector& selector,
                                                    QString* err_out) const {
    const std::vector<int32_t> matches = select(scope, selector);
    if (matches.size() == 1) return matches.front();
    if (err_out) {
        QStringList ids;
        for (int32_t m : matches) ids << QString::number(m);
        *err_out = QStringLiteral("typed atom selector expected one match, got %1 [%2]")
                       .arg(matches.size())
                       .arg(ids.join(QStringLiteral(",")));
    }
    return std::nullopt;
}

}  // namespace h5reader::rediscover
