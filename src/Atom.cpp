#include "Atom.h"

namespace nmr {

std::unique_ptr<Atom> Atom::Create(Element elem) {
    auto a = std::make_unique<Atom>();
    a->element = elem;
    return a;
}

std::unique_ptr<Atom> Atom::Create(const std::string& sym) {
    return Create(ElementFromSymbol(sym));
}

}  // namespace nmr
