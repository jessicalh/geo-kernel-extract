#include "ScopedProducerCatalog.h"

namespace h5reader::rediscover {
namespace {

bool excludedGroup(io::FieldGroup g) {
    switch (g) {
    case io::FieldGroup::LarsenHBond:
    case io::FieldGroup::Tripeptide:
    case io::FieldGroup::Delta:
    case io::FieldGroup::Rediscover:
    case io::FieldGroup::McConnellLegacy:
    case io::FieldGroup::MOPACMcConnellLegacy:
    case io::FieldGroup::MOPACCoulombLegacy:
    case io::FieldGroup::CoulombLegacy:
        return true;
    default:
        return false;
    }
}

}  // namespace

bool IsScopedProducerField(const io::FieldSpec& spec) {
    return !excludedGroup(spec.group);
}

const std::vector<const io::FieldSpec*>& ScopedProducerCatalog() {
    static const std::vector<const io::FieldSpec*> fields = [] {
        std::vector<const io::FieldSpec*> out;
        for (const io::FieldSpec& spec : io::kFieldCatalog) {
            if (IsScopedProducerField(spec)) out.push_back(&spec);
        }
        return out;
    }();
    return fields;
}

}  // namespace h5reader::rediscover
