#pragma once
//
// GeometryChoice: runtime record of one geometric decision made by a calculator.
//
// A GeometryChoice is a bag of model objects — ConformationAtom*, Ring*, Bond*,
// and named numbers — recording what a calculator decided, which entities were
// involved, and whether each was included or excluded.
//
// The conformation owns a flat vector<GeometryChoice>. Calculators populate it
// during Compute() via the GeometryChoiceBuilder factory, which enforces that
// every choice is populated within a lambda (keeping the recording code visually
// separate from the physics).
//
// The UI walks the list, follows pointers back into the live model, and draws.
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


// ============================================================================
// Enums
// ============================================================================

enum class EntityRole   { Source, Target, Context };
enum class EntityOutcome { Included, Excluded, Triggered, NotTriggered };


// ============================================================================
// GeometryEntity — one entry in the bag.
//
// Exactly one of {atom, ring, bond} is non-null per entry.
// Named numbers use the NamedNumber struct instead.
// ============================================================================

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


// ============================================================================
// NamedNumber — a numeric value with a name and unit.
//
// "horizon", 15.0, "A"
// "intensity", -12.0, "nA"
// "distance", 7.3, "A"
// ============================================================================

struct NamedNumber {
    std::string name;
    double      value = 0.0;
    std::string unit;
};


// ============================================================================
// GeometryChoice — the bag itself.
// ============================================================================

class GeometryChoice {
    friend class GeometryChoiceBuilder;
    friend void AddAtom(GeometryChoice&, const ConformationAtom*, size_t,
                        EntityRole, EntityOutcome, const char*);
    friend void AddRing(GeometryChoice&, const Ring*,
                        EntityRole, EntityOutcome, const char*);
    friend void AddBond(GeometryChoice&, const Bond*,
                        EntityRole, EntityOutcome, const char*);
    friend void AddNumber(GeometryChoice&, const char*, double, const char*);

public:
    // Read-only access for the UI
    const std::string&              Label()      const { return label_; }
    CalculatorId                    Calculator() const { return calculator_; }
    size_t                          GroupKey()   const { return group_key_; }
    const std::vector<GeometryEntity>& Entities() const { return entities_; }
    const std::vector<NamedNumber>&    Numbers()  const { return numbers_; }

private:
    GeometryChoice() = default;

    std::string label_;
    CalculatorId calculator_ = CalculatorId::BiotSavart;
    size_t group_key_ = 0;
    std::vector<GeometryEntity> entities_;
    std::vector<NamedNumber>    numbers_;
};


// ============================================================================
// GeometryChoiceBuilder — factory that enforces population in a lambda.
//
// Usage in a calculator's Compute():
//
//   GeometryChoiceBuilder choices(conf);
//
//   choices.Record(CalculatorId::BiotSavart, group, "ring horizon",
//       [&](GeometryChoice& gc) {
//           AddRing(gc, ring, EntityRole::Source, EntityOutcome::Included);
//           AddAtom(gc, &conf.AtomAt(ai), ai, EntityRole::Target, EntityOutcome::Included);
//           AddNumber(gc, "horizon", 15.0, "A");
//           AddNumber(gc, "distance", dist, "A");
//       });
//
// The lambda body is the ONLY place entities and numbers are added.
// This keeps recording code visually offset from physics code.
// ============================================================================

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


// ============================================================================
// Convenience methods on GeometryChoice for use inside the populate lambda.
// ============================================================================

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


}  // namespace nmr
