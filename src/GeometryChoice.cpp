#include "GeometryChoice.h"
#include "ProteinConformation.h"

namespace nmr {

GeometryChoiceBuilder::GeometryChoiceBuilder(ProteinConformation& conf)
    : conf_(conf) {}

void GeometryChoiceBuilder::Record(CalculatorId calculator,
                                   size_t group_key,
                                   const char* label,
                                   std::function<void(GeometryChoice&)> populate) {
    GeometryChoice gc;
    gc.label_ = label;
    gc.calculator_ = calculator;
    gc.group_key_ = group_key;
    populate(gc);
    conf_.geometry_choices.push_back(std::move(gc));
}

}  // namespace nmr
