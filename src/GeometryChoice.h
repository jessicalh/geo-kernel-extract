#pragma once
//
// GeometryChoice: a record of one geometric decision a calculator made —
// which model entities (ConformationAtom*/Ring*/Bond*) and named numbers
// were involved, and whether each was included or excluded. The conformation
// owns a flat vector<GeometryChoice>, populated during Compute() via
// GeometryChoiceBuilder.
//
// Entity roles:
//   Source  — the ring/bond/atom generating the field
//   Target  — the atom being shielded or evaluated
//   Context — referenced but not directly source or target (e.g. bond endpoints)
//
// Outcome:
//   Included     — entity passed the gate
//   Excluded     — entity failed the gate
//   Triggered    — a guard or clamp fired (singularity, sanity, refinement)
//   NotTriggered — a guard was checked but did not fire
//

#include "Types.h"
#include <string>
#include <vector>
#include <functional>

namespace nmr {

// Forward declarations — these are the live model objects.
class Ring;
class ConformationAtom;
struct Bond;
class ProteinConformation;


enum class EntityRole   { Source, Target, Context };
enum class EntityOutcome { Included, Excluded, Triggered, NotTriggered };


// GeometryEntity — one entry in the bag; exactly one of {atom, ring, bond}
// is non-null.
struct GeometryEntity {
    const ConformationAtom* atom = nullptr;
    const Ring*             ring = nullptr;
    const Bond*             bond = nullptr;
    size_t                  atom_index = SIZE_MAX;  // into protein atom list

    EntityRole    role    = EntityRole::Target;
    EntityOutcome outcome = EntityOutcome::Included;

    // Optional: which filter rejected this entity (empty if not rejected)
    std::string   filter_name;
};


struct NamedNumber {
    std::string name;
    double      value = 0.0;
    std::string unit;
};


class GeometryChoice {
    friend class GeometryChoiceBuilder;
    friend void AddAtom(GeometryChoice&, const ConformationAtom*, size_t,
                        EntityRole, EntityOutcome, const char*);
    friend void AddRing(GeometryChoice&, const Ring*,
                        EntityRole, EntityOutcome, const char*);
    friend void AddBond(GeometryChoice&, const Bond*,
                        EntityRole, EntityOutcome, const char*);
    friend void AddNumber(GeometryChoice&, const char*, double, const char*);
    friend void SetSampler(GeometryChoice&, std::function<SphericalTensor(Vec3)>);

public:
    const std::string&              Label()      const { return label_; }
    CalculatorId                    Calculator() const { return calculator_; }
    size_t                          GroupKey()   const { return group_key_; }
    const std::vector<GeometryEntity>& Entities() const { return entities_; }
    const std::vector<NamedNumber>&    Numbers()  const { return numbers_; }

    // Optional field sampler: a stored callback evaluating a SphericalTensor
    // at a 3D point. SampleAt returns a default SphericalTensor if unset.
    bool HasSampler() const { return sampler_ != nullptr; }
    SphericalTensor SampleAt(Vec3 point) const {
        return sampler_ ? sampler_(point) : SphericalTensor{};
    }

private:
    GeometryChoice() = default;

    std::string label_;
    CalculatorId calculator_ = CalculatorId::BiotSavart;
    size_t group_key_ = 0;
    std::vector<GeometryEntity> entities_;
    std::vector<NamedNumber>    numbers_;
    std::function<SphericalTensor(Vec3)> sampler_;
};


// GeometryChoiceBuilder: entities and numbers are added ONLY inside the
// populate lambda passed to Record().
class GeometryChoiceBuilder {
public:
    explicit GeometryChoiceBuilder(ProteinConformation& conf);

    // Record one geometric decision. The populate lambda fills the choice,
    // then it goes straight onto the conformation's list. No commit step.
    void Record(CalculatorId calculator,
                size_t group_key,
                const char* label,
                std::function<void(GeometryChoice&)> populate);

private:
    ProteinConformation& conf_;
};


// Adders for use inside the populate lambda.
inline void AddAtom(GeometryChoice& gc,
                    const ConformationAtom* atom,
                    size_t atom_index,
                    EntityRole role,
                    EntityOutcome outcome,
                    const char* filter = nullptr)
{
    GeometryEntity e;
    e.atom = atom;
    e.atom_index = atom_index;
    e.role = role;
    e.outcome = outcome;
    if (filter) e.filter_name = filter;
    gc.entities_.push_back(std::move(e));
}

inline void AddRing(GeometryChoice& gc,
                    const Ring* ring,
                    EntityRole role,
                    EntityOutcome outcome,
                    const char* filter = nullptr)
{
    GeometryEntity e;
    e.ring = ring;
    e.role = role;
    e.outcome = outcome;
    if (filter) e.filter_name = filter;
    gc.entities_.push_back(std::move(e));
}

inline void AddBond(GeometryChoice& gc,
                    const Bond* bond,
                    EntityRole role,
                    EntityOutcome outcome,
                    const char* filter = nullptr)
{
    GeometryEntity e;
    e.bond = bond;
    e.role = role;
    e.outcome = outcome;
    if (filter) e.filter_name = filter;
    gc.entities_.push_back(std::move(e));
}

inline void AddNumber(GeometryChoice& gc,
                      const char* name,
                      double value,
                      const char* unit)
{
    gc.numbers_.push_back({name, value, unit});
}

inline void SetSampler(GeometryChoice& gc,
                       std::function<SphericalTensor(Vec3)> sampler)
{
    gc.sampler_ = std::move(sampler);
}


}  // namespace nmr
