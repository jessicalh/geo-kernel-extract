#include "CalculationArea.h"
#include "Protein.h"

namespace nmr {


// ============================================================================
// RadialThreshold
// ============================================================================

RadialThreshold::RadialThreshold(const char* name, const char* description,
                                 DomainKind domain, double radius_angstroms,
                                 Sense sense)
    : name_(name)
    , description_(description)
    , domain_(domain)
    , radius_(radius_angstroms)
    , sense_(sense) {}

const char* RadialThreshold::Name() const { return name_; }
const char* RadialThreshold::Description() const { return description_; }
DomainKind RadialThreshold::Domain() const { return domain_; }
double RadialThreshold::Radius() const { return radius_; }
Sense RadialThreshold::ThresholdSense() const { return sense_; }


// ============================================================================
// SourceRelativeExclusion
// ============================================================================

SourceRelativeExclusion::SourceRelativeExclusion(const char* name,
                                                 const char* description,
                                                 double factor)
    : name_(name)
    , description_(description)
    , factor_(factor) {}

const char* SourceRelativeExclusion::Name() const { return name_; }
const char* SourceRelativeExclusion::Description() const { return description_; }
DomainKind SourceRelativeExclusion::Domain() const { return DomainKind::SourceRelative; }
double SourceRelativeExclusion::Factor() const { return factor_; }


// ============================================================================
// RingBondedExclusion
//
// Walks: Protein -> RingAt(ri) -> atom_indices -> AtomAt(vi) -> bond_indices
//        -> BondAt(bi) -> atom_index_a / atom_index_b
//
// For each ring, the exclusion set contains:
//   - Every vertex atom of the ring
//   - Every atom covalently bonded to any vertex
//
// Mirrors RingBondedExclusionFilter's constructor in
// KernelEvaluationFilter.cpp.
// ============================================================================

RingBondedExclusion::RingBondedExclusion(const char* name,
                                         const char* description,
                                         const Protein& protein)
    : name_(name)
    , description_(description) {
    size_t n_rings = protein.RingCount();
    ring_excluded_.resize(n_rings);

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        auto& excluded = ring_excluded_[ri];

        for (size_t vi : ring.atom_indices) {
            excluded.insert(vi);  // the vertex itself

            const auto& atom = protein.AtomAt(vi);
            for (size_t bi : atom.bond_indices) {
                const auto& bond = protein.BondAt(bi);
                excluded.insert(bond.atom_index_a);
                excluded.insert(bond.atom_index_b);
            }
        }
    }
}

const char* RingBondedExclusion::Name() const { return name_; }
const char* RingBondedExclusion::Description() const { return description_; }
DomainKind RingBondedExclusion::Domain() const { return DomainKind::Topological; }

static const std::set<size_t> kEmptyAtomSet;

const std::set<size_t>& RingBondedExclusion::ExcludedAtoms(size_t ring_index) const {
    if (ring_index >= ring_excluded_.size()) return kEmptyAtomSet;
    return ring_excluded_[ring_index];
}

size_t RingBondedExclusion::RingCount() const { return ring_excluded_.size(); }


// ============================================================================
// SelfSourceExclusion
// ============================================================================

SelfSourceExclusion::SelfSourceExclusion(const char* name,
                                         const char* description)
    : name_(name)
    , description_(description) {}

const char* SelfSourceExclusion::Name() const { return name_; }
const char* SelfSourceExclusion::Description() const { return description_; }
DomainKind SelfSourceExclusion::Domain() const { return DomainKind::Topological; }


// ============================================================================
// ShellBoundary
// ============================================================================

ShellBoundary::ShellBoundary(const char* name, const char* description,
                             double bin_edge_angstroms)
    : name_(name)
    , description_(description)
    , bin_edge_(bin_edge_angstroms) {}

const char* ShellBoundary::Name() const { return name_; }
const char* ShellBoundary::Description() const { return description_; }
DomainKind ShellBoundary::Domain() const { return DomainKind::Shell; }
double ShellBoundary::BinEdge() const { return bin_edge_; }


// ============================================================================
// SequenceGate
// ============================================================================

SequenceGate::SequenceGate(const char* name, const char* description,
                           int min_separation)
    : name_(name)
    , description_(description)
    , min_separation_(min_separation) {}

const char* SequenceGate::Name() const { return name_; }
const char* SequenceGate::Description() const { return description_; }
DomainKind SequenceGate::Domain() const { return DomainKind::Sequence; }
int SequenceGate::MinSeparation() const { return min_separation_; }


// ============================================================================
// SwitchingFunction
// ============================================================================

SwitchingFunction::SwitchingFunction(const char* name, const char* description,
                                     double onset_angstroms,
                                     double cutoff_angstroms)
    : name_(name)
    , description_(description)
    , onset_(onset_angstroms)
    , cutoff_(cutoff_angstroms) {}

const char* SwitchingFunction::Name() const { return name_; }
const char* SwitchingFunction::Description() const { return description_; }
DomainKind SwitchingFunction::Domain() const { return DomainKind::Switching; }
double SwitchingFunction::Onset() const { return onset_; }
double SwitchingFunction::Cutoff() const { return cutoff_; }


// ============================================================================
// DecayFunction
// ============================================================================

DecayFunction::DecayFunction(const char* name, const char* description,
                             double decay_length, const char* unit)
    : name_(name)
    , description_(description)
    , decay_length_(decay_length)
    , unit_(unit) {}

const char* DecayFunction::Name() const { return name_; }
const char* DecayFunction::Description() const { return description_; }
DomainKind DecayFunction::Domain() const { return DomainKind::Decay; }
double DecayFunction::DecayLength() const { return decay_length_; }
const char* DecayFunction::Unit() const { return unit_; }


// ============================================================================
// RingCurrent
// ============================================================================

RingCurrent::RingCurrent(const char* name, const char* description,
                         RingTypeIndex ring_type, double intensity_nanoamperes)
    : name_(name)
    , description_(description)
    , ring_type_(ring_type)
    , intensity_(intensity_nanoamperes) {}

const char* RingCurrent::Name() const { return name_; }
const char* RingCurrent::Description() const { return description_; }
DomainKind RingCurrent::Domain() const { return DomainKind::RingMagnitude; }
RingTypeIndex RingCurrent::RingType() const { return ring_type_; }
double RingCurrent::Intensity() const { return intensity_; }


// ============================================================================
// LobeOffset
// ============================================================================

LobeOffset::LobeOffset(const char* name, const char* description,
                       RingTypeIndex ring_type, double offset_angstroms)
    : name_(name)
    , description_(description)
    , ring_type_(ring_type)
    , offset_(offset_angstroms) {}

const char* LobeOffset::Name() const { return name_; }
const char* LobeOffset::Description() const { return description_; }
DomainKind LobeOffset::Domain() const { return DomainKind::RingGeometry; }
RingTypeIndex LobeOffset::RingType() const { return ring_type_; }
double LobeOffset::Offset() const { return offset_; }


// ============================================================================
// NumericalAccuracy
// ============================================================================

NumericalAccuracy::NumericalAccuracy(const char* name, const char* description,
                                     double threshold_angstroms)
    : name_(name)
    , description_(description)
    , threshold_(threshold_angstroms) {}

const char* NumericalAccuracy::Name() const { return name_; }
const char* NumericalAccuracy::Description() const { return description_; }
DomainKind NumericalAccuracy::Domain() const { return DomainKind::Numerical; }
double NumericalAccuracy::Threshold() const { return threshold_; }


// ============================================================================
// ValueClamp
// ============================================================================

ValueClamp::ValueClamp(const char* name, const char* description,
                       double limit, const char* unit)
    : name_(name)
    , description_(description)
    , limit_(limit)
    , unit_(unit) {}

const char* ValueClamp::Name() const { return name_; }
const char* ValueClamp::Description() const { return description_; }
DomainKind ValueClamp::Domain() const { return DomainKind::Numerical; }
double ValueClamp::Limit() const { return limit_; }
const char* ValueClamp::Unit() const { return unit_; }


// ============================================================================
// ValueGate
// ============================================================================

ValueGate::ValueGate(const char* name, const char* description,
                     double floor, const char* unit)
    : name_(name)
    , description_(description)
    , floor_(floor)
    , unit_(unit) {}

const char* ValueGate::Name() const { return name_; }
const char* ValueGate::Description() const { return description_; }
DomainKind ValueGate::Domain() const { return DomainKind::ValueThreshold; }
double ValueGate::Floor() const { return floor_; }
const char* ValueGate::Unit() const { return unit_; }


// ============================================================================
// SentinelValue
// ============================================================================

SentinelValue::SentinelValue(const char* name, const char* description,
                             double value)
    : name_(name)
    , description_(description)
    , value_(value) {}

const char* SentinelValue::Name() const { return name_; }
const char* SentinelValue::Description() const { return description_; }
DomainKind SentinelValue::Domain() const { return DomainKind::Sentinel; }
double SentinelValue::Value() const { return value_; }


// ============================================================================
// WholeDomain
// ============================================================================

WholeDomain::WholeDomain(const char* name, const char* description)
    : name_(name)
    , description_(description) {}

const char* WholeDomain::Name() const { return name_; }
const char* WholeDomain::Description() const { return description_; }
DomainKind WholeDomain::Domain() const { return DomainKind::WholeDomain; }



}  // namespace nmr
