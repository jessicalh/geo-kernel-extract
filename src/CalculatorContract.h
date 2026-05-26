#pragma once
//
// CalculatorContract: a calculator's declared topology + charge-table
// contract. The two are separate type parameters so topology and charge
// models stay orthogonal.
//

#include "ProteinConformation.h"
#include "Protein.h"
#include "ForceFieldChargeTable.h"
#include <type_traits>

namespace nmr {

struct NoChargeTable {};

template<class TopologyT, class ChargeTableT = NoChargeTable>
struct CalculatorContract {
    static_assert(std::is_base_of_v<ProteinTopology, TopologyT>,
                  "calculator topology must derive from ProteinTopology");
    using Topology = TopologyT;
    using ChargeTable = ChargeTableT;
};

template<class CalculatorT>
const typename CalculatorT::Contract::Topology&
RequiredTopology(const ProteinConformation& conf) {
    return conf.ProteinRef()
        .template TopologyAs<typename CalculatorT::Contract::Topology>();
}

template<class CalculatorT>
const typename CalculatorT::Contract::ChargeTable&
RequiredChargeTable(const ProteinConformation& conf) {
    using ChargeTableT = typename CalculatorT::Contract::ChargeTable;
    static_assert(std::is_same_v<ChargeTableT, ForceFieldChargeTable>,
                  "only ForceFieldChargeTable is wired in this landing");
    return conf.ProteinRef().ForceFieldCharges();
}

}  // namespace nmr
